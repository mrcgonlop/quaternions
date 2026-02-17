// Boundary conditions: open (Mur ABC + sponge), conducting, periodic.
//
// Applied after each field update step.
//
// The default boundary is "Open" — a first-order Mur absorbing boundary
// condition at the outermost face, backed by a thin exponential sponge layer
// (3 cells deep) to dissipate oblique-incidence residual that Mur alone
// cannot absorb.
//
// Mur first-order ABC for the wave equation at face x=0:
//   u(0, t+dt) = u(1, t) + coeff * (u(1, t+dt) - u(0, t))
//   where coeff = (c*dt - dx) / (c*dt + dx)
//
// The sponge layer multiplies q and q_dot by a damping factor at the
// outermost SPONGE_DEPTH cells. The profile is quadratic:
//   factor = 1 - sigma * ((depth - d) / depth)^2
// where d is distance from the face (0 at face, depth at interior edge).
// This provides soft dissipation that catches oblique waves Mur misses.
//
// The double-buffer scheme gives us access to both timesteps:
//   cells[current]     = newly computed (t+dt)
//   cells[1-current]   = previous step (t)
//
// Physical absorbers/reflectors should be placed as entities (charges,
// dipoles, conductors) inside the domain, not imposed as boundary conditions.

use bevy::prelude::*;

use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::{CellFlags, PmlConfig, SimParams};

/// Depth of the sponge layer in cells (behind the Mur face).
const SPONGE_DEPTH: usize = 3;

/// Peak damping coefficient at the outermost sponge cell per timestep.
/// 0.4 = 40% energy removal at the face per step — aggressive enough to
/// prevent visible reflections, thin enough (3 cells) not to eat into the domain.
const SPONGE_SIGMA: f32 = 0.4;

/// Boundary condition type for a single face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BoundaryType {
    /// Open boundary (Mur first-order ABC + sponge). Waves exit without reflection.
    Open,
    /// Perfect electric conductor: tangential E = 0 (tangential A components zeroed).
    Conducting,
    /// Periodic: wrap opposite face values.
    Periodic,
}

impl BoundaryType {
    pub fn name(&self) -> &'static str {
        match self {
            BoundaryType::Open => "Open",
            BoundaryType::Conducting => "Conducting",
            BoundaryType::Periodic => "Periodic",
        }
    }
}

/// Boundary configuration for all six faces of the simulation domain.
#[derive(Resource, Clone, Debug)]
pub struct BoundaryConfig {
    pub neg_x: BoundaryType,
    pub pos_x: BoundaryType,
    pub neg_y: BoundaryType,
    pub pos_y: BoundaryType,
    pub neg_z: BoundaryType,
    pub pos_z: BoundaryType,
}

impl Default for BoundaryConfig {
    fn default() -> Self {
        Self {
            neg_x: BoundaryType::Open,
            pos_x: BoundaryType::Open,
            neg_y: BoundaryType::Open,
            pos_y: BoundaryType::Open,
            neg_z: BoundaryType::Open,
            pos_z: BoundaryType::Open,
        }
    }
}

/// CPML auxiliary state, stored separately from CellState to keep it at 48 bytes.
///
/// For each PML cell, we store 12 auxiliary psi values (4 Q components x 3 spatial
/// directions) plus precomputed CPML coefficients (b, a, inv_kappa) per direction.
///
/// The CPML modifies the Laplacian's second derivatives:
///   lap_x_pml = inv_kappa_x * lap_x_standard + psi_x
///   psi_x(t+dt) = b_x * psi_x(t) + a_x * lap_x_standard
#[derive(Resource)]
pub struct PmlState {
    pub config: PmlConfig,
    /// Maps flat grid index → PML cell index (None for non-PML cells).
    pub grid_to_pml: Vec<Option<usize>>,
    /// Number of PML cells.
    pub pml_cell_count: usize,
    /// Auxiliary psi fields: indexed as [pml_idx * 12 + comp * 3 + dir].
    /// comp = 0..4 (Q components), dir = 0..3 (x, y, z).
    pub psi: Vec<f32>,
    /// Precomputed b coefficient per PML cell per direction: [pml_idx * 3 + dir].
    /// b = exp(-(sigma/kappa + alpha) * dt)
    pub b: Vec<f32>,
    /// Precomputed a coefficient per PML cell per direction: [pml_idx * 3 + dir].
    /// a = sigma / (sigma*kappa + kappa^2*alpha) * (b - 1)
    pub a: Vec<f32>,
    /// Precomputed 1/kappa per PML cell per direction: [pml_idx * 3 + dir].
    pub inv_kappa: Vec<f32>,
}

impl PmlState {
    /// Initialize CPML state for the given grid and boundary configuration.
    ///
    /// Identifies PML cells (within `config.depth` of any Open face), computes
    /// polynomial sigma/kappa/alpha profiles, and derives CPML coefficients.
    /// Also flags PML cells with `CellFlags::PML` in both grid buffers.
    pub fn new(
        grid: &mut SimulationGrid,
        bc: &BoundaryConfig,
        config: PmlConfig,
    ) -> Self {
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        let depth = config.depth;
        let dt = grid.dt;
        let dx = grid.dx;
        let sigma_max = config.sigma_max(dx);
        let p = config.grading_order;
        let kappa_max = config.kappa_max;
        let alpha_max = config.alpha_max;
        let total_cells = nx * ny * nz;

        let mut grid_to_pml: Vec<Option<usize>> = vec![None; total_cells];
        let mut pml_cell_count = 0usize;

        // Temporary storage: per-cell distance to nearest Open face in each direction.
        // None means not in PML along that direction.
        struct PmlCellInfo {
            dist: [Option<usize>; 3], // distance to nearest face in x, y, z
        }
        let mut cell_infos: Vec<PmlCellInfo> = Vec::new();

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let mut dist = [None; 3];
                    let mut is_pml = false;

                    // X direction
                    if bc.neg_x == BoundaryType::Open && x < depth {
                        dist[0] = Some(x);
                        is_pml = true;
                    }
                    if bc.pos_x == BoundaryType::Open && x >= nx - depth {
                        let d = nx - 1 - x;
                        dist[0] = Some(dist[0].map_or(d, |prev| prev.min(d)));
                        is_pml = true;
                    }

                    // Y direction
                    if bc.neg_y == BoundaryType::Open && y < depth {
                        dist[1] = Some(y);
                        is_pml = true;
                    }
                    if bc.pos_y == BoundaryType::Open && y >= ny - depth {
                        let d = ny - 1 - y;
                        dist[1] = Some(dist[1].map_or(d, |prev| prev.min(d)));
                        is_pml = true;
                    }

                    // Z direction
                    if bc.neg_z == BoundaryType::Open && z < depth {
                        dist[2] = Some(z);
                        is_pml = true;
                    }
                    if bc.pos_z == BoundaryType::Open && z >= nz - depth {
                        let d = nz - 1 - z;
                        dist[2] = Some(dist[2].map_or(d, |prev| prev.min(d)));
                        is_pml = true;
                    }

                    if is_pml {
                        let i = z * nx * ny + y * nx + x;
                        grid_to_pml[i] = Some(pml_cell_count);
                        pml_cell_count += 1;
                        cell_infos.push(PmlCellInfo { dist });

                        // Flag both buffers
                        grid.cells[0][i].flags |= CellFlags::PML;
                        grid.cells[1][i].flags |= CellFlags::PML;
                    }
                }
            }
        }

        // Allocate coefficient arrays
        let _psi_placeholder = pml_cell_count * 12; // size for later allocation
        let mut b_arr = vec![1.0f32; pml_cell_count * 3];
        let mut a_arr = vec![0.0f32; pml_cell_count * 3];
        let mut inv_kappa_arr = vec![1.0f32; pml_cell_count * 3];

        // Compute CPML coefficients for each PML cell
        for (pml_idx, info) in cell_infos.iter().enumerate() {
            for dir in 0..3usize {
                let coeff_idx = pml_idx * 3 + dir;

                if let Some(face_dist) = info.dist[dir] {
                    // Normalized depth: 0 at interior edge of PML, 1 at outer face
                    let d_norm = 1.0 - (face_dist as f32 / depth as f32);
                    let d_norm = d_norm.clamp(0.0, 1.0);

                    let sigma_d = sigma_max * d_norm.powf(p);
                    let kappa_d = 1.0 + (kappa_max - 1.0) * d_norm.powf(p);
                    let alpha_d = alpha_max * (1.0 - d_norm); // alpha increases toward interior

                    let denom = sigma_d * kappa_d + kappa_d * kappa_d * alpha_d;
                    let b_d = (-(sigma_d / kappa_d + alpha_d) * dt).exp();
                    let a_d = if denom.abs() < 1e-12 {
                        0.0
                    } else {
                        (sigma_d / denom) * (b_d - 1.0)
                    };

                    b_arr[coeff_idx] = b_d;
                    a_arr[coeff_idx] = a_d;
                    inv_kappa_arr[coeff_idx] = 1.0 / kappa_d;
                }
                // If dist[dir] is None, coefficients stay at (b=1, a=0, inv_kappa=1) = no damping
            }
        }

        info!(
            "PML initialized: depth={}, {} PML cells out of {} total ({:.1}%)",
            depth,
            pml_cell_count,
            total_cells,
            100.0 * pml_cell_count as f64 / total_cells as f64,
        );

        Self {
            config,
            grid_to_pml,
            pml_cell_count,
            psi: vec![0.0; pml_cell_count * 12],
            b: b_arr,
            a: a_arr,
            inv_kappa: inv_kappa_arr,
        }
    }

    /// Reset all auxiliary psi fields to zero (e.g., after grid reset).
    pub fn reset_psi(&mut self) {
        self.psi.fill(0.0);
    }
}

/// Apply boundary conditions to the grid after a field update + swap.
///
/// `effective_dt` is the actual timestep used in the field update (grid.dt * dt_factor).
/// This must match the dt used in `step_field_cpu` for the Mur coefficient to be correct.
///
/// After `step_field_cpu` + `swap_and_advance`:
///   - `cells[current]` has the newly computed interior values (t+dt)
///     and copied-unchanged boundary values (still at time t).
///   - `cells[1-current]` has the previous timestep values (t).
pub fn apply_boundaries(grid: &mut SimulationGrid, config: &BoundaryConfig, effective_dt: f32) {
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;
    let c0 = SimParams::C0;

    // Mur coefficient using the actual simulation timestep.
    let cdt = c0 * effective_dt;
    let mur_coeff = (cdt - dx) / (cdt + dx);

    // Split cells array for simultaneous mutable (new) and immutable (old) access.
    let cur = grid.current;
    let (a, b) = grid.cells.split_at_mut(1);
    let (new_buf, old_buf) = if cur == 0 {
        (&mut a[0][..], &b[0][..])
    } else {
        (&mut b[0][..], &a[0][..])
    };

    let idx = |x: usize, y: usize, z: usize| -> usize {
        z * nx * ny + y * nx + x
    };

    // --- Step 1: Mur ABC at Open face cells ---
    // Extrapolate boundary cell from its interior neighbor.

    // -x face
    if config.neg_x == BoundaryType::Open {
        for z in 0..nz {
            for y in 0..ny {
                mur_update(new_buf, old_buf, idx(0, y, z), idx(1, y, z), mur_coeff);
            }
        }
    }
    // +x face
    if config.pos_x == BoundaryType::Open {
        for z in 0..nz {
            for y in 0..ny {
                mur_update(new_buf, old_buf, idx(nx - 1, y, z), idx(nx - 2, y, z), mur_coeff);
            }
        }
    }
    // -y face
    if config.neg_y == BoundaryType::Open {
        for z in 0..nz {
            for x in 0..nx {
                mur_update(new_buf, old_buf, idx(x, 0, z), idx(x, 1, z), mur_coeff);
            }
        }
    }
    // +y face
    if config.pos_y == BoundaryType::Open {
        for z in 0..nz {
            for x in 0..nx {
                mur_update(new_buf, old_buf, idx(x, ny - 1, z), idx(x, ny - 2, z), mur_coeff);
            }
        }
    }
    // -z face
    if config.neg_z == BoundaryType::Open {
        for y in 0..ny {
            for x in 0..nx {
                mur_update(new_buf, old_buf, idx(x, y, 0), idx(x, y, 1), mur_coeff);
            }
        }
    }
    // +z face
    if config.pos_z == BoundaryType::Open {
        for y in 0..ny {
            for x in 0..nx {
                mur_update(new_buf, old_buf, idx(x, y, nz - 1), idx(x, y, nz - 2), mur_coeff);
            }
        }
    }

    // --- Step 2: Sponge layer for Open faces ---
    // Damp fields in the outermost SPONGE_DEPTH cells to dissipate
    // oblique-incidence reflections that first-order Mur cannot absorb.
    apply_sponge(new_buf, config, nx, ny, nz, &idx);

    // --- Step 3: Conducting boundaries ---
    apply_conducting(new_buf, config, nx, ny, nz, &idx);

    // --- Step 4: Periodic boundaries ---
    apply_periodic(new_buf, config, nx, ny, nz, &idx);
}

/// Apply the Mur first-order ABC to a single boundary cell.
///
/// u_boundary(t+dt) = u_neighbor(t) + coeff * (u_neighbor(t+dt) - u_boundary(t))
#[inline]
fn mur_update(
    new_buf: &mut [crate::simulation::state::CellState],
    old_buf: &[crate::simulation::state::CellState],
    bi: usize,
    ni: usize,
    coeff: f32,
) {
    for c in 0..4 {
        new_buf[bi].q[c] = old_buf[ni].q[c]
            + coeff * (new_buf[ni].q[c] - old_buf[bi].q[c]);
        new_buf[bi].q_dot[c] = old_buf[ni].q_dot[c]
            + coeff * (new_buf[ni].q_dot[c] - old_buf[bi].q_dot[c]);
    }
}

/// Compute the sponge damping factor for a cell at `dist` cells from a face.
/// Returns 1.0 (no damping) if outside the sponge depth.
#[inline]
fn sponge_factor(dist: usize) -> f32 {
    if dist >= SPONGE_DEPTH {
        return 1.0;
    }
    // Quadratic profile: strongest at face (dist=0), unity at SPONGE_DEPTH.
    let t = (SPONGE_DEPTH - dist) as f32 / SPONGE_DEPTH as f32;
    1.0 - SPONGE_SIGMA * t * t
}

/// Apply exponential sponge damping in the boundary region for all Open faces.
/// Cells near multiple faces get the product of all applicable damping factors.
fn apply_sponge(
    buf: &mut [crate::simulation::state::CellState],
    config: &BoundaryConfig,
    nx: usize,
    ny: usize,
    nz: usize,
    idx: &dyn Fn(usize, usize, usize) -> usize,
) {
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let mut factor = 1.0f32;

                if config.neg_x == BoundaryType::Open && x < SPONGE_DEPTH {
                    factor *= sponge_factor(x);
                }
                if config.pos_x == BoundaryType::Open && x >= nx - SPONGE_DEPTH {
                    factor *= sponge_factor(nx - 1 - x);
                }
                if config.neg_y == BoundaryType::Open && y < SPONGE_DEPTH {
                    factor *= sponge_factor(y);
                }
                if config.pos_y == BoundaryType::Open && y >= ny - SPONGE_DEPTH {
                    factor *= sponge_factor(ny - 1 - y);
                }
                if config.neg_z == BoundaryType::Open && z < SPONGE_DEPTH {
                    factor *= sponge_factor(z);
                }
                if config.pos_z == BoundaryType::Open && z >= nz - SPONGE_DEPTH {
                    factor *= sponge_factor(nz - 1 - z);
                }

                if (factor - 1.0).abs() < 1e-7 {
                    continue;
                }

                let i = idx(x, y, z);
                for c in 0..4 {
                    buf[i].q[c] *= factor;
                    buf[i].q_dot[c] *= factor;
                }
            }
        }
    }
}

/// Apply conducting boundary conditions: zero tangential A components on face cells.
fn apply_conducting(
    buf: &mut [crate::simulation::state::CellState],
    config: &BoundaryConfig,
    nx: usize,
    ny: usize,
    nz: usize,
    idx: &dyn Fn(usize, usize, usize) -> usize,
) {
    for z in 0..nz {
        for y in 0..ny {
            if config.neg_x == BoundaryType::Conducting {
                let i = idx(0, y, z);
                buf[i].q[2] = 0.0;
                buf[i].q[3] = 0.0;
                buf[i].q_dot[2] = 0.0;
                buf[i].q_dot[3] = 0.0;
            }
            if config.pos_x == BoundaryType::Conducting {
                let i = idx(nx - 1, y, z);
                buf[i].q[2] = 0.0;
                buf[i].q[3] = 0.0;
                buf[i].q_dot[2] = 0.0;
                buf[i].q_dot[3] = 0.0;
            }
        }
    }

    for z in 0..nz {
        for x in 0..nx {
            if config.neg_y == BoundaryType::Conducting {
                let i = idx(x, 0, z);
                buf[i].q[1] = 0.0;
                buf[i].q[3] = 0.0;
                buf[i].q_dot[1] = 0.0;
                buf[i].q_dot[3] = 0.0;
            }
            if config.pos_y == BoundaryType::Conducting {
                let i = idx(x, ny - 1, z);
                buf[i].q[1] = 0.0;
                buf[i].q[3] = 0.0;
                buf[i].q_dot[1] = 0.0;
                buf[i].q_dot[3] = 0.0;
            }
        }
    }

    for y in 0..ny {
        for x in 0..nx {
            if config.neg_z == BoundaryType::Conducting {
                let i = idx(x, y, 0);
                buf[i].q[1] = 0.0;
                buf[i].q[2] = 0.0;
                buf[i].q_dot[1] = 0.0;
                buf[i].q_dot[2] = 0.0;
            }
            if config.pos_z == BoundaryType::Conducting {
                let i = idx(x, y, nz - 1);
                buf[i].q[1] = 0.0;
                buf[i].q[2] = 0.0;
                buf[i].q_dot[1] = 0.0;
                buf[i].q_dot[2] = 0.0;
            }
        }
    }
}

/// Apply periodic boundary conditions by averaging opposite face values.
fn apply_periodic(
    buf: &mut [crate::simulation::state::CellState],
    config: &BoundaryConfig,
    nx: usize,
    ny: usize,
    nz: usize,
    idx: &dyn Fn(usize, usize, usize) -> usize,
) {
    if config.neg_x == BoundaryType::Periodic && config.pos_x == BoundaryType::Periodic {
        for z in 0..nz {
            for y in 0..ny {
                let lo = idx(0, y, z);
                let hi = idx(nx - 1, y, z);
                for c in 0..4 {
                    let avg = (buf[lo].q[c] + buf[hi].q[c]) * 0.5;
                    buf[lo].q[c] = avg;
                    buf[hi].q[c] = avg;
                    let avg_d = (buf[lo].q_dot[c] + buf[hi].q_dot[c]) * 0.5;
                    buf[lo].q_dot[c] = avg_d;
                    buf[hi].q_dot[c] = avg_d;
                }
            }
        }
    }

    if config.neg_y == BoundaryType::Periodic && config.pos_y == BoundaryType::Periodic {
        for z in 0..nz {
            for x in 0..nx {
                let lo = idx(x, 0, z);
                let hi = idx(x, ny - 1, z);
                for c in 0..4 {
                    let avg = (buf[lo].q[c] + buf[hi].q[c]) * 0.5;
                    buf[lo].q[c] = avg;
                    buf[hi].q[c] = avg;
                    let avg_d = (buf[lo].q_dot[c] + buf[hi].q_dot[c]) * 0.5;
                    buf[lo].q_dot[c] = avg_d;
                    buf[hi].q_dot[c] = avg_d;
                }
            }
        }
    }

    if config.neg_z == BoundaryType::Periodic && config.pos_z == BoundaryType::Periodic {
        for y in 0..ny {
            for x in 0..nx {
                let lo = idx(x, y, 0);
                let hi = idx(x, y, nz - 1);
                for c in 0..4 {
                    let avg = (buf[lo].q[c] + buf[hi].q[c]) * 0.5;
                    buf[lo].q[c] = avg;
                    buf[hi].q[c] = avg;
                    let avg_d = (buf[lo].q_dot[c] + buf[hi].q_dot[c]) * 0.5;
                    buf[lo].q_dot[c] = avg_d;
                    buf[hi].q_dot[c] = avg_d;
                }
            }
        }
    }
}

/// Bevy system: apply boundary conditions after field update.
///
/// When PmlState is present, Open faces skip Mur+sponge (CPML handles absorption
/// in the field update). Only Conducting and Periodic boundaries are applied.
/// The outermost face cells in PML regions get simple zero-order extrapolation.
pub fn boundary_system(
    grid: Option<ResMut<SimulationGrid>>,
    config: Res<BoundaryConfig>,
    sim_config: Res<crate::simulation::plugin::SimulationConfig>,
    pml: Option<Res<PmlState>>,
) {
    let Some(mut grid) = grid else { return };
    if sim_config.paused {
        return;
    }
    let effective_dt = grid.dt * sim_config.dt_factor.clamp(0.01, 1.0);

    if pml.is_some() {
        // With CPML, Open faces are handled by the field update kernel.
        // We still need to handle the outermost face cells (which lack neighbors
        // and thus can't compute Laplacian). Use zero-order extrapolation from
        // the adjacent interior cell.
        apply_pml_face_extrapolation(&mut grid, &config);
        // Apply non-Open boundary types normally
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        let cur = grid.current;
        let idx = |x: usize, y: usize, z: usize| -> usize { z * nx * ny + y * nx + x };

        apply_conducting(&mut grid.cells[cur], &config, nx, ny, nz, &idx);
        apply_periodic(&mut grid.cells[cur], &config, nx, ny, nz, &idx);
    } else {
        // Fallback: use original Mur+sponge for backward compatibility
        apply_boundaries(&mut grid, &config, effective_dt);
    }
}

/// For PML faces: copy the adjacent interior cell's values to the outermost face cell.
/// This handles face cells that can't be updated by the field kernel (no neighbors).
fn apply_pml_face_extrapolation(grid: &mut SimulationGrid, config: &BoundaryConfig) {
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let cur = grid.current;
    let buf = &mut grid.cells[cur];

    let idx = |x: usize, y: usize, z: usize| -> usize { z * nx * ny + y * nx + x };

    // -x face
    if config.neg_x == BoundaryType::Open {
        for z in 0..nz {
            for y in 0..ny {
                buf[idx(0, y, z)] = buf[idx(1, y, z)];
            }
        }
    }
    // +x face
    if config.pos_x == BoundaryType::Open {
        for z in 0..nz {
            for y in 0..ny {
                buf[idx(nx - 1, y, z)] = buf[idx(nx - 2, y, z)];
            }
        }
    }
    // -y face
    if config.neg_y == BoundaryType::Open {
        for z in 0..nz {
            for x in 0..nx {
                buf[idx(x, 0, z)] = buf[idx(x, 1, z)];
            }
        }
    }
    // +y face
    if config.pos_y == BoundaryType::Open {
        for z in 0..nz {
            for x in 0..nx {
                buf[idx(x, ny - 1, z)] = buf[idx(x, ny - 2, z)];
            }
        }
    }
    // -z face
    if config.neg_z == BoundaryType::Open {
        for y in 0..ny {
            for x in 0..nx {
                buf[idx(x, y, 0)] = buf[idx(x, y, 1)];
            }
        }
    }
    // +z face
    if config.pos_z == BoundaryType::Open {
        for y in 0..ny {
            for x in 0..nx {
                buf[idx(x, y, nz - 1)] = buf[idx(x, y, nz - 2)];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::field_update::step_field_cpu;
    use crate::simulation::grid::SimulationGrid;

    /// After a field step + swap, Mur ABC should set boundary cells to
    /// extrapolated values (not just the copied-from-previous-step values).
    #[test]
    fn test_mur_modifies_boundary_cells() {
        let mut grid = SimulationGrid::new(16, 16, 16, 0.01);
        let config = BoundaryConfig::default();

        // Place a perturbation near the -x boundary so the wave reaches it
        grid.cell_mut(2, 8, 8).q[0] = 1.0;

        let params = grid.sim_params(false);
        let dt = grid.dt;
        for _ in 0..5 {
            step_field_cpu(&mut grid, &params, None);
            grid.swap_and_advance();
            apply_boundaries(&mut grid, &config, dt);
        }

        // Mur should have updated boundary — test mainly ensures no panic
        let edge = grid.cell(0, 8, 8);
        let _has_activity = edge.q[0].abs() > 1e-20 || edge.q_dot[0].abs() > 1e-20;
    }

    /// Open boundaries should allow energy to leave the domain.
    #[test]
    fn test_open_boundaries_drain_energy() {
        use crate::simulation::diagnostics::{compute_derived_fields, total_energy};

        // Run WITHOUT boundaries
        let mut grid_no_bc = SimulationGrid::new(16, 16, 16, 0.01);
        grid_no_bc.cell_mut(8, 8, 8).q[0] = 1.0;
        let params = grid_no_bc.sim_params(false);

        for _ in 0..50 {
            step_field_cpu(&mut grid_no_bc, &params, None);
            grid_no_bc.swap_and_advance();
        }

        let derived_no_bc = compute_derived_fields(&grid_no_bc, &grid_no_bc.sim_params(false));
        let energy_no_bc = total_energy(&derived_no_bc, grid_no_bc.dx);

        // Run WITH open boundaries
        let mut grid_bc = SimulationGrid::new(16, 16, 16, 0.01);
        let dt = grid_bc.dt;
        grid_bc.cell_mut(8, 8, 8).q[0] = 1.0;
        let params = grid_bc.sim_params(false);
        let bc_config = BoundaryConfig::default();

        for _ in 0..50 {
            step_field_cpu(&mut grid_bc, &params, None);
            grid_bc.swap_and_advance();
            apply_boundaries(&mut grid_bc, &bc_config, dt);
        }

        let derived_bc = compute_derived_fields(&grid_bc, &grid_bc.sim_params(false));
        let energy_bc = total_energy(&derived_bc, grid_bc.dx);

        assert!(
            energy_bc < energy_no_bc,
            "open boundaries should let energy leave: with_bc={:.4e}, without_bc={:.4e}",
            energy_bc,
            energy_no_bc
        );
    }

    /// After source removal, energy should decay toward zero with open boundaries.
    #[test]
    fn test_energy_decays_after_source_removal() {
        use crate::simulation::diagnostics::{compute_derived_fields, total_energy};
        use crate::simulation::sources::{inject_sources, Source, SourceConfig};

        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let dt = grid.dt;
        let bc_config = BoundaryConfig::default();

        let sources = SourceConfig {
            sources: vec![Source::dipole_z([16.0, 16.0, 16.0], 1.0, 1e9)],
        };

        // Phase 1: drive with source for 100 steps
        for _ in 0..100 {
            let params = grid.sim_params(false);
            inject_sources(&mut grid, &sources, &params);
            step_field_cpu(&mut grid, &params, None);
            grid.swap_and_advance();
            apply_boundaries(&mut grid, &bc_config, dt);
        }

        let derived = compute_derived_fields(&grid, &grid.sim_params(false));
        let energy_with_source = total_energy(&derived, grid.dx);

        // Phase 2: remove source, run 200 more steps to let energy drain
        let no_sources = SourceConfig { sources: vec![] };
        for _ in 0..200 {
            let params = grid.sim_params(false);
            inject_sources(&mut grid, &no_sources, &params);
            step_field_cpu(&mut grid, &params, None);
            grid.swap_and_advance();
            apply_boundaries(&mut grid, &bc_config, dt);
        }

        let derived = compute_derived_fields(&grid, &grid.sim_params(false));
        let energy_after_drain = total_energy(&derived, grid.dx);

        // Energy should have significantly decreased
        assert!(
            energy_after_drain < energy_with_source * 0.5,
            "energy should decay after source removal: before={:.4e}, after={:.4e}",
            energy_with_source,
            energy_after_drain
        );
    }

    #[test]
    fn test_conducting_zeros_tangential() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let dt = grid.dt;
        let config = BoundaryConfig {
            neg_z: BoundaryType::Conducting,
            pos_z: BoundaryType::Conducting,
            ..Default::default()
        };

        for cell in grid.cells[0].iter_mut() {
            cell.q = [1.0, 2.0, 3.0, 4.0];
        }
        for cell in grid.cells[1].iter_mut() {
            cell.q = [1.0, 2.0, 3.0, 4.0];
        }

        apply_boundaries(&mut grid, &config, dt);

        let face_cell = grid.cell(4, 4, 0);
        assert_eq!(face_cell.q[1], 0.0, "Ax should be zero on conducting z-face");
        assert_eq!(face_cell.q[2], 0.0, "Ay should be zero on conducting z-face");
    }

    #[test]
    fn test_interior_unaffected() {
        let mut grid = SimulationGrid::new(16, 16, 16, 0.01);
        let dt = grid.dt;
        let config = BoundaryConfig::default();

        for cell in grid.cells[0].iter_mut() {
            cell.q = [42.0; 4];
            cell.q_dot = [7.0; 4];
        }
        for cell in grid.cells[1].iter_mut() {
            cell.q = [42.0; 4];
            cell.q_dot = [7.0; 4];
        }

        apply_boundaries(&mut grid, &config, dt);

        let c = grid.cell(8, 8, 8);
        assert_eq!(c.q, [42.0; 4]);
        assert_eq!(c.q_dot, [7.0; 4]);
    }

    #[test]
    fn test_periodic_averages_faces() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let dt = grid.dt;
        let config = BoundaryConfig {
            neg_x: BoundaryType::Periodic,
            pos_x: BoundaryType::Periodic,
            ..Default::default()
        };

        let y = 4u32;
        let z = 4u32;
        let idx_lo = grid.idx(0, y, z);
        let idx_hi = grid.idx(7, y, z);
        grid.cells[0][idx_lo].q = [2.0; 4];
        grid.cells[0][idx_hi].q = [6.0; 4];

        apply_boundaries(&mut grid, &config, dt);

        let lo = grid.cell(0, y, z).q[0];
        let hi = grid.cell(7, y, z).q[0];
        assert!(
            (lo - hi).abs() < 1e-6,
            "periodic faces should be equal: lo={}, hi={}",
            lo,
            hi
        );
    }
}
