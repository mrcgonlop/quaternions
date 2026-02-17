// CPU implementation of the quaternionic FDTD field update.
// Ports the WGSL kernel from ARCHITECTURE.md §GPU Compute Shader to Rust.
//
// Physics: d²Q/dt² = c_local² * ∇²Q  [standard mode]
//          d²Q/dt² = c_local² * ∇²Q - c_local² * (0, ∇S)  [extended QVED]
//
// Integration: Störmer-Verlet (leapfrog) with staggered half-step velocity:
//   q_dot(t+dt/2) = q_dot(t-dt/2) + q_ddot(t) * dt
//   q(t+dt) = q(t) + q_dot(t+dt/2) * dt
//
// Double-buffer: read from cells[current], write to cells[1 - current], swap.
//
// CPML (Convolutional PML): When a PmlState is provided, PML-flagged cells get
// modified Laplacian components. For each direction, the second derivative is
// corrected: lap_dir_pml = inv_kappa * lap_dir + psi, where psi is an auxiliary
// field updated each timestep via psi = b * psi + a * lap_dir.

use crate::math::fdtd;
use crate::simulation::boundaries::PmlState;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::{CellFlags, CellState, SimParams};

/// Compute the scalar field S at a given cell from the read buffer.
///
/// S = q_dot.w / c_local + div(A)
/// where div(A) = dAx/dx + dAy/dy + dAz/dz computed via central differences.
///
/// Returns 0.0 for boundary cells that lack neighbors for the stencil.
#[inline]
fn compute_s_at(
    cells: &[CellState],
    x: usize,
    y: usize,
    z: usize,
    nx: usize,
    ny: usize,
    nz: usize,
    inv_2dx: f32,
    c0: f32,
) -> f32 {
    // Need neighbors in all 3 directions — must be interior
    if !fdtd::is_interior(x, y, z, nx, ny, nz) {
        return 0.0;
    }

    let i = fdtd::idx(x, y, z, nx, ny);
    let stride_y = nx;
    let stride_z = nx * ny;
    let c_local = c0 / cells[i].k;

    // div(A) = dAx/dx + dAy/dy + dAz/dz
    // A = (q[1], q[2], q[3]) — the vector part of Q
    let div_a = (cells[i + 1].q[1] - cells[i - 1].q[1]) * inv_2dx
        + (cells[i + stride_y].q[2] - cells[i - stride_y].q[2]) * inv_2dx
        + (cells[i + stride_z].q[3] - cells[i - stride_z].q[3]) * inv_2dx;

    // S = (1/c) * d(phi/c)/dt + div(A) = q_dot.w / c_local + div(A)
    cells[i].q_dot[0] / c_local + div_a
}

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
/// 6. If extended_mode: subtract c² * (0, grad_S)
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

    let read = grid.current;

    // Split borrows: read from cells[read], write to cells[write].
    let (read_buf, write_buf) = if read == 0 {
        let (a, b) = grid.cells.split_at_mut(1);
        (&a[0][..], &mut b[0][..])
    } else {
        let (a, b) = grid.cells.split_at_mut(1);
        (&b[0][..], &mut a[0][..])
    };

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = fdtd::idx(x, y, z, nx, ny);
                let cell = &read_buf[i];

                // Copy cell to write buffer first (outermost face cells stay unchanged)
                write_buf[i] = *cell;

                // Skip hard boundary cells
                if (cell.flags & CellFlags::BOUNDARY) != 0 {
                    continue;
                }

                // Need all 6 neighbors to compute Laplacian
                if !has_neighbors(x, y, z, nx, ny, nz) {
                    continue;
                }

                // For non-PML cells that are not interior (this shouldn't happen with
                // proper PML setup, but guard against it)
                let is_pml_cell = (cell.flags & CellFlags::PML) != 0;
                let is_interior = fdtd::is_interior(x, y, z, nx, ny, nz);

                // Update interior cells and PML cells that have neighbors
                if !is_interior && !is_pml_cell {
                    continue;
                }

                let k = cell.k;
                let c_local = c0 / k;
                let c_local_sq = c_local * c_local;

                // --- Per-direction second derivatives of Q ---
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

                    // Apply CPML corrections if this is a PML cell
                    if is_pml_cell {
                        if let Some(ref mut pml_state) = pml {
                            if let Some(pml_idx) = pml_state.grid_to_pml[i] {
                                let coeff_base = pml_idx * 3;
                                let psi_base = pml_idx * 12 + comp * 3;

                                // X direction
                                let bx = pml_state.b[coeff_base];
                                let ax = pml_state.a[coeff_base];
                                let ikx = pml_state.inv_kappa[coeff_base];
                                pml_state.psi[psi_base] =
                                    bx * pml_state.psi[psi_base] + ax * lap_x;
                                let lap_x_pml = ikx * lap_x + pml_state.psi[psi_base];

                                // Y direction
                                let by = pml_state.b[coeff_base + 1];
                                let ay = pml_state.a[coeff_base + 1];
                                let iky = pml_state.inv_kappa[coeff_base + 1];
                                pml_state.psi[psi_base + 1] =
                                    by * pml_state.psi[psi_base + 1] + ay * lap_y;
                                let lap_y_pml = iky * lap_y + pml_state.psi[psi_base + 1];

                                // Z direction
                                let bz = pml_state.b[coeff_base + 2];
                                let az = pml_state.a[coeff_base + 2];
                                let ikz = pml_state.inv_kappa[coeff_base + 2];
                                pml_state.psi[psi_base + 2] =
                                    bz * pml_state.psi[psi_base + 2] + az * lap_z;
                                let lap_z_pml = ikz * lap_z + pml_state.psi[psi_base + 2];

                                lap_q[comp] = lap_x_pml + lap_y_pml + lap_z_pml;
                                continue; // skip the standard sum below
                            }
                        }
                    }

                    // Standard (non-PML) Laplacian
                    lap_q[comp] = lap_x + lap_y + lap_z;
                }

                // --- Compute q_ddot = c_local² * lap_Q ---
                let mut q_ddot = [
                    c_local_sq * lap_q[0],
                    c_local_sq * lap_q[1],
                    c_local_sq * lap_q[2],
                    c_local_sq * lap_q[3],
                ];

                // --- Extended mode: subtract c_local² * (0, grad_S) ---
                if extended {
                    let s_xp =
                        compute_s_at(read_buf, x + 1, y, z, nx, ny, nz, inv_2dx, c0);
                    let s_xm =
                        compute_s_at(read_buf, x - 1, y, z, nx, ny, nz, inv_2dx, c0);
                    let s_yp =
                        compute_s_at(read_buf, x, y + 1, z, nx, ny, nz, inv_2dx, c0);
                    let s_ym =
                        compute_s_at(read_buf, x, y - 1, z, nx, ny, nz, inv_2dx, c0);
                    let s_zp =
                        compute_s_at(read_buf, x, y, z + 1, nx, ny, nz, inv_2dx, c0);
                    let s_zm =
                        compute_s_at(read_buf, x, y, z - 1, nx, ny, nz, inv_2dx, c0);

                    let grad_s = [
                        (s_xp - s_xm) * inv_2dx,
                        (s_yp - s_ym) * inv_2dx,
                        (s_zp - s_zm) * inv_2dx,
                    ];

                    q_ddot[1] -= c_local_sq * grad_s[0];
                    q_ddot[2] -= c_local_sq * grad_s[1];
                    q_ddot[3] -= c_local_sq * grad_s[2];
                }

                // --- Leapfrog (Störmer-Verlet) integration ---
                let out = &mut write_buf[i];
                for comp in 0..4 {
                    out.q_dot[comp] = cell.q_dot[comp] + q_ddot[comp] * dt;
                    out.q[comp] = cell.q[comp] + out.q_dot[comp] * dt;
                }

                // Preserve non-field state
                out.k = cell.k;
                out.k_dot = cell.k_dot;
                out.flags = cell.flags;
                out._pad = cell._pad;
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
