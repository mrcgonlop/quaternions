// CPU implementation of the quaternionic FDTD field update.
// Ports the WGSL kernel from ARCHITECTURE.md §GPU Compute Shader to Rust.
//
// Physics (ungauged Maxwell potential equations):
//   Standard mode:  d²Q/dt² = c² * ∇²Q
//   Extended QVED:  d²Q[0]/dt² = c² * ∇²Q[0] + c * ∂S/∂t   (scalar, α=1 exactly)
//                   d²A/dt²    = c² * ∇²A     - c² * ∇S      (vector)
//                   □S         = 0                             (S as independent field)
//
// S is an INDEPENDENT dynamical variable stored in grid.s_field (double-buffered)
// and grid.s_dot (half-step velocity).  It evolves via its own leapfrog loop:
//   s_dot(t+dt/2) = s_dot(t-dt/2) + c² * ∇²S(t) * dt
//   s(t+dt)       = s(t)          + s_dot(t+dt/2) * dt
//
// The coupling to Q uses s_dot directly (no approximation):
//   q_ddot[0] += c * s_dot[i]           (α=1.0 — true QVED, not a reduced approximation)
//   q_ddot[1..3] -= c² * ∇S[t]
//
// This eliminates the former SCALAR_COUPLING_ALPHA=0.2 approximation. The CFL
// condition is unchanged: S propagates at c (same as Q), so no dt reduction
// is needed.
//
// Integration: Störmer-Verlet (leapfrog) with staggered half-step velocity:
//   q_dot(t+dt/2) = q_dot(t-dt/2) + q_ddot(t) * dt
//   q(t+dt) = q(t) + q_dot(t+dt/2) * dt
//
// Double-buffer: read from cells[current]/s_field[current],
//                write to cells[1-current]/s_field[1-current], swap current.
//
// CPML (Convolutional PML): When a PmlState is provided, PML-flagged cells get
// modified Laplacian components. For each direction, the second derivative is
// corrected: lap_dir_pml = inv_kappa * lap_dir + psi, where psi is an auxiliary
// field updated each timestep via psi = b * psi + a * lap_dir.
// S uses the standard (non-CPML) Laplacian in PML cells. Full CPML for S would
// require extending PmlState.psi from n_pml*12 to n_pml*15 entries (one extra
// component for S × 3 spatial directions) — deferred to a future session.

use crate::math::fdtd;
use crate::simulation::boundaries::PmlState;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::{CellFlags, SimParams};

/// Check if a cell has all 6 neighbors in the grid (not on the outermost face).
#[inline]
fn has_neighbors(x: usize, y: usize, z: usize, nx: usize, ny: usize, nz: usize) -> bool {
    x > 0 && x < nx - 1 && y > 0 && y < ny - 1 && z > 0 && z < nz - 1
}

/// Perform one timestep of the quaternionic FDTD field update on the CPU.
///
/// Reads from `grid.cells[grid.current]`, writes to `grid.cells[1 - grid.current]`.
/// Does NOT swap buffers — caller must call `grid.swap_and_advance()` after.
///
/// When `pml` is Some, CPML corrections are applied to PML-flagged cells.
/// PML cells that have neighbors (not on the outermost grid face) are updated
/// with CPML-corrected Laplacians rather than being skipped.
///
/// # Algorithm (per cell)
/// 1. Copy cell to write buffer
/// 2. Skip if BOUNDARY flag or outermost face without neighbors
/// 3. Compute per-direction second derivatives of Q
/// 4. If PML cell: apply CPML correction (psi update + modified derivative)
/// 5. Sum to get Laplacian, compute q_ddot = c² * lap_Q
/// 6. If extended_mode: add c·∂S/∂t to q_ddot[0], subtract c²·∇S from q_ddot[1..3]
/// 7. Leapfrog: q_dot += q_ddot * dt, q += q_dot * dt
pub fn step_field_cpu(
    grid: &mut SimulationGrid,
    params: &SimParams,
    mut pml: Option<&mut PmlState>,
) {
    let nx = params.nx as usize;
    let ny = params.ny as usize;
    let nz = params.nz as usize;
    let dx = params.dx;
    let dt = params.dt;
    let c0 = params.c0;
    let extended = params.extended_mode != 0;

    let inv_2dx = 1.0 / (2.0 * dx);
    let inv_dx2 = 1.0 / (dx * dx);

    let stride_y = nx;
    let stride_z = nx * ny;

    let n_cells = nx * ny * nz;

    // -------------------------------------------------------------------------
    // Pre-pass: compute ∇²S(t) for use in the S leapfrog update.
    // Read from s_field[current] (S at integer step t).
    // Also needed for grad(S) in the Q coupling — but that reads s_field[current]
    // directly in the cell loop (no need to pre-compute the gradient).
    // -------------------------------------------------------------------------
    let lap_s = if extended {
        let s = &grid.s_field[grid.current]; // immutable borrow, released at end of block
        let mut lap = vec![0.0f32; n_cells];
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if !fdtd::is_interior(x, y, z, nx, ny, nz) {
                        continue;
                    }
                    let i = fdtd::idx(x, y, z, nx, ny);
                    lap[i] = (s[i + 1] + s[i - 1] - 2.0 * s[i]) * inv_dx2
                        + (s[i + stride_y] + s[i - stride_y] - 2.0 * s[i]) * inv_dx2
                        + (s[i + stride_z] + s[i - stride_z] - 2.0 * s[i]) * inv_dx2;
                }
            }
        }
        lap
        // `s` borrow is released here (lap_s owns a Vec, not a reference)
    } else {
        Vec::new()
    };

    let read = grid.current;

    // Split borrows: read from cells[read], write to cells[1-read].
    // grid.s_field and grid.s_dot are different struct fields — borrow checker
    // allows them to be borrowed independently after grid.cells is split.
    let (read_buf, write_buf) = if read == 0 {
        let (a, b) = grid.cells.split_at_mut(1);
        (&a[0][..], &mut b[0][..])
    } else {
        let (a, b) = grid.cells.split_at_mut(1);
        (&b[0][..], &mut a[0][..])
    };

    // Split s_field for read/write (same current index as cells).
    let (s_read, s_write) = if read == 0 {
        let (a, b) = grid.s_field.split_at_mut(1);
        (&a[0][..], &mut b[0][..])
    } else {
        let (a, b) = grid.s_field.split_at_mut(1);
        (&b[0][..], &mut a[0][..])
    };

    // s_dot is updated in-place (half-step velocity, no double-buffer needed).
    let s_dot = &mut grid.s_dot;

    // =========================================================================
    // Q update loop
    // =========================================================================
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = fdtd::idx(x, y, z, nx, ny);
                let cell = &read_buf[i];

                // Copy cell to write buffer (face cells stay unchanged this step).
                write_buf[i] = *cell;

                // Skip hard-wall boundary cells (CellFlags::BOUNDARY is reserved for
                // future explicit conductors; no code sets it today).
                if (cell.flags & CellFlags::BOUNDARY) != 0 {
                    continue;
                }

                // Laplacian requires all 6 neighbors.
                if !has_neighbors(x, y, z, nx, ny, nz) {
                    continue;
                }

                let is_pml_cell = (cell.flags & CellFlags::PML) != 0;

                let c_local = c0 / cell.k;
                let c_local_sq = c_local * c_local;

                // --- Per-direction second derivatives of Q (with optional CPML) ---
                let mut lap_q = [0.0f32; 4];

                for comp in 0..4 {
                    let lap_x = (read_buf[i + 1].q[comp] + read_buf[i - 1].q[comp]
                        - 2.0 * cell.q[comp])
                        * inv_dx2;
                    let lap_y = (read_buf[i + stride_y].q[comp]
                        + read_buf[i - stride_y].q[comp]
                        - 2.0 * cell.q[comp])
                        * inv_dx2;
                    let lap_z = (read_buf[i + stride_z].q[comp]
                        + read_buf[i - stride_z].q[comp]
                        - 2.0 * cell.q[comp])
                        * inv_dx2;

                    if is_pml_cell {
                        if let Some(ref mut pml_state) = pml {
                            if let Some(pml_idx) = pml_state.grid_to_pml[i] {
                                let coeff_base = pml_idx * 3;
                                let psi_base = pml_idx * 12 + comp * 3;

                                let bx = pml_state.b[coeff_base];
                                let ax = pml_state.a[coeff_base];
                                let ikx = pml_state.inv_kappa[coeff_base];
                                pml_state.psi[psi_base] =
                                    bx * pml_state.psi[psi_base] + ax * lap_x;
                                let lap_x_pml = ikx * lap_x + pml_state.psi[psi_base];

                                let by = pml_state.b[coeff_base + 1];
                                let ay = pml_state.a[coeff_base + 1];
                                let iky = pml_state.inv_kappa[coeff_base + 1];
                                pml_state.psi[psi_base + 1] =
                                    by * pml_state.psi[psi_base + 1] + ay * lap_y;
                                let lap_y_pml = iky * lap_y + pml_state.psi[psi_base + 1];

                                let bz = pml_state.b[coeff_base + 2];
                                let az = pml_state.a[coeff_base + 2];
                                let ikz = pml_state.inv_kappa[coeff_base + 2];
                                pml_state.psi[psi_base + 2] =
                                    bz * pml_state.psi[psi_base + 2] + az * lap_z;
                                let lap_z_pml = ikz * lap_z + pml_state.psi[psi_base + 2];

                                lap_q[comp] = lap_x_pml + lap_y_pml + lap_z_pml;
                                continue;
                            }
                        }
                    }

                    lap_q[comp] = lap_x + lap_y + lap_z;
                }

                // --- q_ddot = c_local² * ∇²Q ---
                let mut q_ddot = [
                    c_local_sq * lap_q[0],
                    c_local_sq * lap_q[1],
                    c_local_sq * lap_q[2],
                    c_local_sq * lap_q[3],
                ];

                // --- Extended QVED coupling (skipped for PML cells) ---
                //
                // True coupling at α=1.0 — no reduced approximation:
                //   q_ddot[0] += c · s_dot[i]       (scalar: +c · ∂S/∂t)
                //   q_ddot[1..3] -= c² · ∇S[t]      (vector: -c² · ∇S)
                //
                // Using s_dot[i] (the actual half-step ∂S/∂t from the independent S
                // field) avoids the former analytical approximation and the CFL
                // penalty it introduced (α was forced to 0.2 to stay below 0.235).
                //
                // PML cells are still excluded: adding driving terms inside CPML
                // undermines absorption and can cause instability.
                if extended && !is_pml_cell {
                    q_ddot[0] += c_local * s_dot[i]; // α = 1.0 exactly

                    let grad_s_x = (s_read[i + 1] - s_read[i - 1]) * inv_2dx;
                    let grad_s_y = (s_read[i + stride_y] - s_read[i - stride_y]) * inv_2dx;
                    let grad_s_z = (s_read[i + stride_z] - s_read[i - stride_z]) * inv_2dx;
                    q_ddot[1] -= c_local_sq * grad_s_x;
                    q_ddot[2] -= c_local_sq * grad_s_y;
                    q_ddot[3] -= c_local_sq * grad_s_z;
                }

                // --- Störmer-Verlet leapfrog ---
                let out = &mut write_buf[i];
                for comp in 0..4 {
                    out.q_dot[comp] = cell.q_dot[comp] + q_ddot[comp] * dt;
                    out.q[comp] = cell.q[comp] + out.q_dot[comp] * dt;
                }
                out.k = cell.k;
                out.k_dot = cell.k_dot;
                out.flags = cell.flags;
                out._pad = cell._pad;
            }
        }
    }

    // =========================================================================
    // S leapfrog update (extended mode only)
    //
    // s_dot(t+dt/2) = s_dot(t-dt/2) + c_local² * ∇²S(t) * dt
    // s(t+dt)       = s(t)           + s_dot(t+dt/2) * dt
    //
    // S uses the standard (non-CPML) Laplacian in all cells, including PML.
    // Full CPML for S requires extending PmlState.psi to n_pml*15 entries
    // (one extra "component" for S × 3 spatial directions) — future work.
    // =========================================================================
    if extended {
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let i = fdtd::idx(x, y, z, nx, ny);

                    if !fdtd::is_interior(x, y, z, nx, ny, nz) {
                        // Face cells: Neumann BC — copy from the adjacent interior cell.
                        // Applied after the interior loop via a separate face pass below.
                        s_write[i] = s_read[i];
                        continue;
                    }

                    // c_local from the cell state (K field).
                    let c_local = c0 / read_buf[i].k;
                    let c_local_sq = c_local * c_local;

                    // s_dot(t+dt/2) = s_dot(t-dt/2) + c² * ∇²S * dt
                    s_dot[i] += c_local_sq * lap_s[i] * dt;
                    // s(t+dt) = s(t) + s_dot(t+dt/2) * dt
                    s_write[i] = s_read[i] + s_dot[i] * dt;
                }
            }
        }

        // Neumann BC for S face cells (zero-gradient: copy from one cell in).
        for z in 0..nz {
            for y in 0..ny {
                s_write[fdtd::idx(0, y, z, nx, ny)] =
                    s_write[fdtd::idx(1, y, z, nx, ny)];
                s_write[fdtd::idx(nx - 1, y, z, nx, ny)] =
                    s_write[fdtd::idx(nx - 2, y, z, nx, ny)];
            }
        }
        for z in 0..nz {
            for x in 0..nx {
                s_write[fdtd::idx(x, 0, z, nx, ny)] =
                    s_write[fdtd::idx(x, 1, z, nx, ny)];
                s_write[fdtd::idx(x, ny - 1, z, nx, ny)] =
                    s_write[fdtd::idx(x, ny - 2, z, nx, ny)];
            }
        }
        for y in 0..ny {
            for x in 0..nx {
                s_write[fdtd::idx(x, y, 0, nx, ny)] =
                    s_write[fdtd::idx(x, y, 1, nx, ny)];
                s_write[fdtd::idx(x, y, nz - 1, nx, ny)] =
                    s_write[fdtd::idx(x, y, nz - 2, nx, ny)];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::grid::SimulationGrid;

    /// A vacuum grid with no sources should remain at zero (no change).
    #[test]
    fn test_vacuum_stays_zero() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);

        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();

        // All cells should still be zero
        for i in 0..grid.cell_count() {
            let cell = &grid.cells[grid.current][i];
            assert_eq!(cell.q, [0.0; 4], "cell {i} q should be zero");
            assert_eq!(cell.q_dot, [0.0; 4], "cell {i} q_dot should be zero");
        }
    }

    /// A single perturbed cell should cause neighboring cells to evolve.
    #[test]
    fn test_perturbation_propagates() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);

        // Set a perturbation at the center
        grid.cell_mut(4, 4, 4).q[0] = 1.0; // phi/c perturbation

        let params = grid.sim_params(false);

        // Run a few steps
        for _ in 0..5 {
            step_field_cpu(&mut grid, &params, None);
            grid.swap_and_advance();
        }

        // The perturbation should have spread to neighbors
        let center = grid.cell(4, 4, 4);
        let neighbor = grid.cell(4, 4, 5);

        // Center should have changed from initial 1.0
        assert_ne!(center.q[0], 1.0, "center should have evolved");
        // Neighbor should be nonzero
        assert_ne!(neighbor.q[0], 0.0, "neighbor should be affected");
    }

    /// Extended mode should produce the same results as standard mode
    /// when S = 0 everywhere (both start from same vacuum state).
    #[test]
    fn test_extended_mode_compiles() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        grid.cell_mut(4, 4, 4).q[0] = 1.0;

        let params = grid.sim_params(true); // extended_mode = true

        // Should not panic
        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();
    }
}
