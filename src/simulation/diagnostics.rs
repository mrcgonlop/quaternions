// Derived field computation: E, B, S, energy density, Poynting vector.
// Computes observable electromagnetic quantities from the quaternionic potential Q.
//
// Physics:
//   E = -grad(phi) - dA/dt = -c * grad(Q.w) - Q_dot.vector()
//   B = curl(A) = curl(Q.vector())
//   S = Q_dot.w / c_local + div(A)   (Lorenz gauge term — zero in standard EM)
//   energy_density = 0.5 * epsilon_0 * |E|^2 + 0.5 / mu_0 * |B|^2 + 0.5 * epsilon_0 * S^2
//   poynting = (1/mu_0) * E x B

use bevy::prelude::*;

use crate::math::fdtd;
use crate::math::vector_field::{self, Vec3f};
use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::{CellFlags, DerivedFields, SimParams};

// Physical constants (SI)
const EPSILON_0: f32 = 8.854e-12; // F/m
const MU_0: f32 = 1.2566e-6; // H/m
const INV_MU_0: f32 = 1.0 / MU_0;

/// Summary diagnostics computed each frame for UI display.
#[derive(Resource, Default, Clone, Debug)]
pub struct DiagnosticsState {
    /// Derived fields for every cell (recomputed each frame).
    pub fields: Vec<DerivedFields>,

    /// Total electromagnetic energy summed over all cells (Joules).
    pub total_energy: f64,

    /// Maximum |E| across all cells.
    pub max_e: f32,

    /// Maximum |B| across all cells.
    pub max_b: f32,

    /// Maximum |S| across all cells.
    pub max_s: f32,
}

/// Compute derived electromagnetic fields from the current grid state.
///
/// For each cell, computes E, B, S, energy density, and Poynting vector
/// from the quaternionic potential Q and its time derivative Q_dot.
///
/// Boundary cells get zero-valued derived fields since central differences
/// require neighbors on all sides.
pub fn compute_derived_fields(grid: &SimulationGrid, params: &SimParams) -> Vec<DerivedFields> {
    let nx = params.nx as usize;
    let ny = params.ny as usize;
    let nz = params.nz as usize;
    let dx = params.dx;
    let c0 = params.c0;
    let n = nx * ny * nz;

    let inv_2dx = 1.0 / (2.0 * dx);

    let cells = grid.read_buf();
    let mut derived = vec![DerivedFields::default(); n];

    // Extract component arrays for use with vector_field helpers.
    // We need: phi (q[0] * c), Ax (q[1]), Ay (q[2]), Az (q[3])
    // For gradient/curl/divergence operations, we need flat scalar arrays.
    //
    // To avoid allocating separate arrays, we compute derivatives inline
    // using the same central-difference approach as vector_field.rs.

    let stride_y = nx;
    let stride_z = nx * ny;

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = fdtd::idx(x, y, z, nx, ny);

                // Boundary cells: leave as default (zeros)
                if !fdtd::is_interior(x, y, z, nx, ny, nz) {
                    derived[i].c_local = c0 / cells[i].k;
                    continue;
                }

                let cell = &cells[i];
                let c_local = c0 / cell.k;

                // --- E = -c * grad(Q.w) - Q_dot.vector() ---
                // grad(Q.w) via central differences on q[0]
                let grad_qw: Vec3f = [
                    (cells[i + 1].q[0] - cells[i - 1].q[0]) * inv_2dx,
                    (cells[i + stride_y].q[0] - cells[i - stride_y].q[0]) * inv_2dx,
                    (cells[i + stride_z].q[0] - cells[i - stride_z].q[0]) * inv_2dx,
                ];

                // E = -c * grad(phi/c) - dA/dt = -c * grad(Q.w) - Q_dot.vector()
                let e: Vec3f = [
                    -c_local * grad_qw[0] - cell.q_dot[1],
                    -c_local * grad_qw[1] - cell.q_dot[2],
                    -c_local * grad_qw[2] - cell.q_dot[3],
                ];

                // --- B = curl(A) = curl(Q.vector()) ---
                // curl computed inline from q[1], q[2], q[3] components
                let b: Vec3f = [
                    // dAz/dy - dAy/dz
                    (cells[i + stride_y].q[3] - cells[i - stride_y].q[3]) * inv_2dx
                        - (cells[i + stride_z].q[2] - cells[i - stride_z].q[2]) * inv_2dx,
                    // dAx/dz - dAz/dx
                    (cells[i + stride_z].q[1] - cells[i - stride_z].q[1]) * inv_2dx
                        - (cells[i + 1].q[3] - cells[i - 1].q[3]) * inv_2dx,
                    // dAy/dx - dAx/dy
                    (cells[i + 1].q[2] - cells[i - 1].q[2]) * inv_2dx
                        - (cells[i + stride_y].q[1] - cells[i - stride_y].q[1]) * inv_2dx,
                ];

                // --- S = Q_dot.w / c_local + div(A) ---
                let div_a = (cells[i + 1].q[1] - cells[i - 1].q[1]) * inv_2dx
                    + (cells[i + stride_y].q[2] - cells[i - stride_y].q[2]) * inv_2dx
                    + (cells[i + stride_z].q[3] - cells[i - stride_z].q[3]) * inv_2dx;
                let s = cell.q_dot[0] / c_local + div_a;

                // --- Energy density ---
                // u = 0.5 * epsilon_0 * |E|^2 + 0.5 / mu_0 * |B|^2 + 0.5 * epsilon_0 * S^2
                let e_sq = vector_field::magnitude_sq(e);
                let b_sq = vector_field::magnitude_sq(b);
                let energy = 0.5 * EPSILON_0 * e_sq + 0.5 * INV_MU_0 * b_sq + 0.5 * EPSILON_0 * s * s;

                // --- Poynting vector: (1/mu_0) * E x B ---
                let poynting = vector_field::scale(vector_field::cross(e, b), INV_MU_0);

                derived[i] = DerivedFields {
                    e,
                    b,
                    s,
                    energy_density: energy,
                    poynting,
                    c_local,
                };
            }
        }
    }

    derived
}

/// Compute total electromagnetic energy by summing energy_density * dx^3 over all cells.
///
/// When `cells` is provided, PML-flagged cells are excluded from the sum so that
/// the artificial PML absorption doesn't appear in energy diagnostics.
pub fn total_energy(derived: &[DerivedFields], dx: f32) -> f64 {
    let cell_volume = (dx as f64).powi(3);
    derived
        .iter()
        .map(|d| d.energy_density as f64 * cell_volume)
        .sum()
}

/// Compute total energy excluding PML cells.
pub fn total_energy_excluding_pml(
    derived: &[DerivedFields],
    grid: &SimulationGrid,
    dx: f32,
) -> f64 {
    let cell_volume = (dx as f64).powi(3);
    let cells = grid.read_buf();
    derived
        .iter()
        .zip(cells.iter())
        .filter(|(_, c)| (c.flags & CellFlags::PML) == 0)
        .map(|(d, _)| d.energy_density as f64 * cell_volume)
        .sum()
}

/// Find the maximum |S| across all cells.
pub fn max_s(derived: &[DerivedFields]) -> f32 {
    derived
        .iter()
        .map(|d| d.s.abs())
        .fold(0.0f32, f32::max)
}

/// Find the maximum |E| across all cells.
pub fn max_e(derived: &[DerivedFields]) -> f32 {
    derived
        .iter()
        .map(|d| vector_field::magnitude_sq(d.e))
        .fold(0.0f32, f32::max)
        .sqrt()
}

/// Find the maximum |B| across all cells.
pub fn max_b(derived: &[DerivedFields]) -> f32 {
    derived
        .iter()
        .map(|d| vector_field::magnitude_sq(d.b))
        .fold(0.0f32, f32::max)
        .sqrt()
}

/// Bevy system: compute derived fields and diagnostics each frame.
///
/// When PmlState is present, energy is computed excluding PML cells.
pub fn diagnostics_system(
    grid: Option<Res<SimulationGrid>>,
    config: Res<crate::simulation::plugin::SimulationConfig>,
    mut diag: ResMut<DiagnosticsState>,
    pml: Option<Res<crate::simulation::boundaries::PmlState>>,
) {
    let Some(grid) = grid else { return };

    let params = grid.sim_params(config.extended_mode);
    let fields = compute_derived_fields(&grid, &params);

    diag.total_energy = if pml.is_some() {
        total_energy_excluding_pml(&fields, &grid, grid.dx)
    } else {
        total_energy(&fields, grid.dx)
    };
    diag.max_e = max_e(&fields);
    diag.max_b = max_b(&fields);
    diag.max_s = max_s(&fields);
    diag.fields = fields;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::grid::SimulationGrid;
    /// Vacuum grid (all zeros) should produce zero derived fields everywhere.
    #[test]
    fn test_vacuum_derived_fields() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        assert_eq!(derived.len(), grid.cell_count());
        for d in &derived {
            assert_eq!(d.e, [0.0; 3]);
            assert_eq!(d.b, [0.0; 3]);
            assert_eq!(d.s, 0.0);
            assert_eq!(d.energy_density, 0.0);
            assert_eq!(d.poynting, [0.0; 3]);
        }
    }

    /// Total energy of vacuum should be zero.
    #[test]
    fn test_vacuum_energy_zero() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);
        let energy = total_energy(&derived, grid.dx);
        assert_eq!(energy, 0.0);
    }

    /// A nonzero phi gradient should produce a nonzero E field.
    #[test]
    fn test_phi_gradient_produces_e_field() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);

        // Set a linear phi ramp along x: phi/c = Q.w = x * dx
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let i = grid.idx(x as u32, y as u32, z as u32);
                    grid.cells[0][i].q[0] = x as f32 * grid.dx;
                }
            }
        }

        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        // Interior cell should have E_x = -c * d(Q.w)/dx = -c * 1.0 (slope is 1)
        let i = grid.idx(4, 4, 4);
        let c0 = SimParams::C0;
        // E_x = -c * grad(Q.w).x = -c * (Q.w[x+1] - Q.w[x-1]) / (2*dx)
        // With Q.w = x*dx: grad = dx / (2*dx) * 2 = 1.0  (central diff of linear = exact slope)
        // Actually: (5*dx - 3*dx) / (2*dx) = 1.0
        // So E_x = -c * 1.0
        assert!(
            (derived[i].e[0] + c0).abs() < c0 * 1e-4,
            "E_x should be -c0, got {}",
            derived[i].e[0]
        );
        // E_y and E_z should be zero (no phi variation in y,z)
        assert!(derived[i].e[1].abs() < 1e-3, "E_y should be ~0");
        assert!(derived[i].e[2].abs() < 1e-3, "E_z should be ~0");
    }

    /// A nonzero vector potential curl should produce a nonzero B field.
    #[test]
    fn test_curl_a_produces_b_field() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);

        // Set Az = y * dx (linear in y), so B_x = dAz/dy = 1.0
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let i = grid.idx(x as u32, y as u32, z as u32);
                    grid.cells[0][i].q[3] = y as f32 * grid.dx; // Az = q[3]
                }
            }
        }

        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        let i = grid.idx(4, 4, 4);
        // B_x = dAz/dy = 1.0 (exact for linear field with central diff)
        assert!(
            (derived[i].b[0] - 1.0).abs() < 1e-4,
            "B_x should be 1.0, got {}",
            derived[i].b[0]
        );
        assert!(derived[i].b[1].abs() < 1e-4, "B_y should be ~0");
        assert!(derived[i].b[2].abs() < 1e-4, "B_z should be ~0");
    }

    /// Nonzero fields should produce positive energy.
    #[test]
    fn test_nonzero_fields_have_energy() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);

        // Set a phi perturbation to create E field
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let i = grid.idx(x as u32, y as u32, z as u32);
                    grid.cells[0][i].q[0] = x as f32 * grid.dx;
                }
            }
        }

        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);
        let energy = total_energy(&derived, grid.dx);
        assert!(energy > 0.0, "energy should be positive for nonzero E field");
    }

    /// max_s should find the peak scalar field value.
    #[test]
    fn test_max_s_detection() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);

        // Inject q_dot.w at center to create nonzero S
        let idx = grid.idx(4, 4, 4);
        grid.cells[0][idx].q_dot[0] = 100.0;

        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);
        let s_max = max_s(&derived);
        assert!(s_max > 0.0, "max_s should detect nonzero S field");
    }
}
