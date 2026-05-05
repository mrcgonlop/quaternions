// Derived field computation: E, B, S, energy density, Poynting vector.
// Computes observable electromagnetic quantities from the quaternionic potential Q.
//
// Physics:
//   E = -grad(phi) - dA/dt = -c * grad(Q.w) - Q_dot.vector()
//   B = curl(A) = curl(Q.vector())
//   S = Q_dot.w / c_local + div(A)   (Lorenz gauge term — zero in standard EM)
//   energy_density = 0.5 * epsilon_0 * |E|^2 + 0.5 / mu_0 * |B|^2 + 0.5 * epsilon_0 * S^2
//   poynting = (1/mu_0) * E x B

use std::collections::VecDeque;

use bevy::prelude::*;

use crate::math::fdtd;
use crate::math::quaternion::Quat;
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

    /// Maximum K (vacuum refractive index) across all non-PML interior cells.
    /// K=1.0 is standard vacuum; values > 1 indicate polarized vacuum regions.
    pub max_k: f32,

    /// Mean K across all non-PML interior cells.
    pub mean_k: f32,

    /// Topological charge (Baryon number / Skyrmion winding number).
    /// Integer for true topological configurations; near-zero for trivial fields.
    pub topological_charge: f32,
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

                // --- S: independent dynamical field in extended mode; Lorenz gauge term otherwise ---
                let div_a = (cells[i + 1].q[1] - cells[i - 1].q[1]) * inv_2dx
                    + (cells[i + stride_y].q[2] - cells[i - stride_y].q[2]) * inv_2dx
                    + (cells[i + stride_z].q[3] - cells[i - stride_z].q[3]) * inv_2dx;
                let s = if params.extended_mode != 0 {
                    // True QVED: S is a free dynamical field evolved by □S = 0.
                    // Read it directly from the grid instead of the Lorenz gauge formula.
                    grid.s_read()[i]
                } else {
                    cell.q_dot[0] / c_local + div_a
                };

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
/// Includes all cells (including PML). Use `total_energy_excluding_pml` to exclude
/// PML cells so that artificial absorption does not appear in energy diagnostics.
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

/// Find the maximum K (vacuum refractive index) across all non-PML interior cells.
///
/// Returns 1.0 (vacuum default) if the grid has no interior cells.
pub fn max_k_field(grid: &SimulationGrid) -> f32 {
    let cells = grid.read_buf();
    cells
        .iter()
        .filter(|c| (c.flags & CellFlags::PML) == 0)
        .map(|c| c.k)
        .fold(1.0f32, f32::max)
}

/// Find the mean K across all non-PML interior cells.
///
/// Returns 1.0 if the grid has no non-PML cells.
pub fn mean_k_field(grid: &SimulationGrid) -> f32 {
    let cells = grid.read_buf();
    let (sum, count) = cells
        .iter()
        .filter(|c| (c.flags & CellFlags::PML) == 0)
        .fold((0.0f64, 0usize), |(s, n), c| (s + c.k as f64, n + 1));
    if count == 0 {
        1.0
    } else {
        (sum / count as f64) as f32
    }
}

/// Find the maximum |B| across all cells.
pub fn max_b(derived: &[DerivedFields]) -> f32 {
    derived
        .iter()
        .map(|d| vector_field::magnitude_sq(d.b))
        .fold(0.0f32, f32::max)
        .sqrt()
}

/// Compute the topological charge (Baryon / Skyrmion number) of the quaternionic field.
///
/// The charge density is:
///   ρ_topo = (1/24π²) ε^{ijk} Tr( L_i · L_j · L_k )
/// where L_i = (∂_i U) * U⁻¹ are left Maurer-Cartan forms and U = Q/|Q| ∈ S³.
///
/// The integral ∫ ρ_topo dV = n ∈ ℤ for topologically non-trivial configurations.
pub fn compute_topological_charge(grid: &SimulationGrid) -> f32 {
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;
    let inv_2dx = 1.0 / (2.0 * dx);
    let cells = grid.read_buf();

    let stride_y = nx;
    let stride_z = nx * ny;

    // The full ε^{ijk} Tr(L_i L_j L_k) reduces to 6 * scalar_part(L_x [L_y, L_z])
    // because: (a) cyclic trace collapses 6 terms to 3 identical pairs,
    //          (b) SU(2) Tr = 2 * quaternion scalar_part.
    // So (1/24π²) * 6 = 1/(4π²).
    let prefactor = 1.0 / (4.0 * std::f32::consts::PI * std::f32::consts::PI);
    let cell_volume = dx * dx * dx;

    let mut charge = 0.0f64;

    for z in 1..(nz - 1) {
        for y in 1..(ny - 1) {
            for x in 1..(nx - 1) {
                let i = fdtd::idx(x, y, z, nx, ny);

                // Skip PML cells
                if (cells[i].flags & CellFlags::PML) != 0 {
                    continue;
                }

                // Build quaternion from cell Q field
                let q_i = Quat::new(cells[i].q[0], cells[i].q[1], cells[i].q[2], cells[i].q[3]);
                let norm = q_i.norm();

                // Skip near-zero field (vacuum or degenerate)
                if norm < 1e-12 {
                    continue;
                }

                // U = Q / |Q| (unit quaternion on S³)
                let u = q_i * (1.0 / norm);
                let u_conj = u.conj();

                // Central differences of the unit quaternion field: ∂U/∂x, ∂U/∂y, ∂U/∂z
                // We need to normalize neighbors too for consistent S³ mapping
                let neighbors = [
                    (i + 1, i - 1),                    // x+1, x-1
                    (i + stride_y, i - stride_y),      // y+1, y-1
                    (i + stride_z, i - stride_z),      // z+1, z-1
                ];

                let mut du = [Quat::zero(); 3]; // ∂U/∂x, ∂U/∂y, ∂U/∂z
                let mut skip = false;

                for (axis, &(ip, im)) in neighbors.iter().enumerate() {
                    let qp = Quat::new(cells[ip].q[0], cells[ip].q[1], cells[ip].q[2], cells[ip].q[3]);
                    let qm = Quat::new(cells[im].q[0], cells[im].q[1], cells[im].q[2], cells[im].q[3]);

                    let np = qp.norm();
                    let nm = qm.norm();

                    if np < 1e-12 || nm < 1e-12 {
                        skip = true;
                        break;
                    }

                    let up = qp * (1.0 / np);
                    let um = qm * (1.0 / nm);

                    du[axis] = (up - um) * inv_2dx;
                }

                if skip {
                    continue;
                }

                // Left Maurer-Cartan forms: L_i = (∂_i U) * U⁻¹ = (∂_i U) * U*
                // (For unit quaternions, U⁻¹ = U*)
                let lx = du[0].hamilton(u_conj);
                let ly = du[1].hamilton(u_conj);
                let lz = du[2].hamilton(u_conj);

                // ρ_topo = (1/24π²) * scalar_part( L_x * (L_y * L_z - L_z * L_y) )
                // The commutator [L_y, L_z] = L_y*L_z - L_z*L_y
                let commutator = ly.hamilton(lz) - lz.hamilton(ly);
                let rho = prefactor * lx.hamilton(commutator).scalar();

                charge += rho as f64 * cell_volume as f64;
            }
        }
    }

    charge as f32
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
    diag.max_k = max_k_field(&grid);
    diag.mean_k = mean_k_field(&grid);

    // Topological charge is expensive — throttle to every 10 steps
    if grid.iteration % 10 == 0 {
        diag.topological_charge = compute_topological_charge(&grid);
    }

    diag.fields = fields;
}

// ---------------------------------------------------------------------------
// Probes: point-sampled time series for quantitative analysis
// ---------------------------------------------------------------------------

/// A scalar field quantity that a `Probe` can record.
///
/// This is a subset of the visualization `FieldQuantity` enum, restricted to
/// quantities that are meaningful as scalar time series at a single cell.
/// Defined locally to avoid a `simulation -> visualization` dependency.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ProbeField {
    SField,
    EMagnitude,
    BMagnitude,
    Phi,
    Ax,
    Ay,
    Az,
    EnergyDensity,
}

impl ProbeField {
    pub const ALL: &'static [ProbeField] = &[
        ProbeField::SField,
        ProbeField::EMagnitude,
        ProbeField::BMagnitude,
        ProbeField::Phi,
        ProbeField::Ax,
        ProbeField::Ay,
        ProbeField::Az,
        ProbeField::EnergyDensity,
    ];

    pub fn name(self) -> &'static str {
        match self {
            ProbeField::SField => "S",
            ProbeField::EMagnitude => "|E|",
            ProbeField::BMagnitude => "|B|",
            ProbeField::Phi => "Phi",
            ProbeField::Ax => "Ax",
            ProbeField::Ay => "Ay",
            ProbeField::Az => "Az",
            ProbeField::EnergyDensity => "u",
        }
    }
}

/// Sample a `ProbeField` at a flat cell index, returning a scalar value.
///
/// Uses the most recently computed `DerivedFields` for E, B, S, energy, and
/// reads Phi / A components directly from the grid.
pub fn sample_probe_field(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    idx: usize,
    field: ProbeField,
) -> f32 {
    match field {
        ProbeField::SField => diag.fields[idx].s,
        ProbeField::EMagnitude => vector_field::magnitude_sq(diag.fields[idx].e).sqrt(),
        ProbeField::BMagnitude => vector_field::magnitude_sq(diag.fields[idx].b).sqrt(),
        ProbeField::Phi => grid.read_buf()[idx].q[0],
        ProbeField::Ax => grid.read_buf()[idx].q[1],
        ProbeField::Ay => grid.read_buf()[idx].q[2],
        ProbeField::Az => grid.read_buf()[idx].q[3],
        ProbeField::EnergyDensity => diag.fields[idx].energy_density,
    }
}

/// A measurement probe: a fixed grid location plus a ring-buffer time series
/// of sampled field values. Populated each simulation step by `probe_system`.
#[derive(Clone, Debug)]
pub struct Probe {
    /// Human-readable label displayed in the UI.
    pub label: String,
    /// Grid-coordinate position (integer indices, clamped to grid on sample).
    pub position: [u32; 3],
    /// Which field quantity to record.
    pub field: ProbeField,
    /// Ring buffer of (time_seconds, value) samples, oldest first.
    pub history: VecDeque<(f32, f32)>,
}

impl Probe {
    pub fn new(label: impl Into<String>, position: [u32; 3], field: ProbeField) -> Self {
        Self {
            label: label.into(),
            position,
            field,
            history: VecDeque::new(),
        }
    }

    /// Drop all recorded samples.
    pub fn clear(&mut self) {
        self.history.clear();
    }

    /// Most recent sample value, or 0.0 if the probe has never been sampled.
    pub fn latest(&self) -> f32 {
        self.history.back().map(|(_, v)| *v).unwrap_or(0.0)
    }
}

/// Collection of probes. Resource inserted by `SimulationPlugin`.
#[derive(Resource, Clone, Debug)]
pub struct ProbeSet {
    pub probes: Vec<Probe>,
    /// Maximum number of samples retained per probe.
    pub max_history: usize,
}

impl Default for ProbeSet {
    fn default() -> Self {
        Self {
            probes: Vec::new(),
            max_history: 4096,
        }
    }
}

impl ProbeSet {
    pub fn clear(&mut self) {
        self.probes.clear();
    }

    pub fn clear_histories(&mut self) {
        for p in &mut self.probes {
            p.clear();
        }
    }

    pub fn push(&mut self, probe: Probe) {
        self.probes.push(probe);
    }
}

/// Discrete Fourier transform of a probe's time series.
///
/// Returns `Vec<(frequency_hz, magnitude)>` for frequencies in
/// `[0, f_nyquist]`. Uses a naive O(N²) DFT — probes are short (typically a
/// few thousand samples) so this runs in well under a millisecond and keeps
/// us off of an fft dependency.
///
/// If the probe has fewer than 2 samples or the sample interval is
/// degenerate, returns an empty vector.
pub fn probe_fft(probe: &Probe) -> Vec<(f32, f32)> {
    let n = probe.history.len();
    if n < 2 {
        return Vec::new();
    }

    let (t0, _) = probe.history[0];
    let (t_last, _) = probe.history[n - 1];
    let total = t_last - t0;
    if total <= 0.0 {
        return Vec::new();
    }
    let dt = total / (n - 1) as f32;
    if dt <= 0.0 {
        return Vec::new();
    }

    // Detrend by subtracting the mean — removes the DC spike that would
    // otherwise dominate the spectrum and hide real oscillation peaks.
    let mean: f32 = probe.history.iter().map(|(_, v)| *v).sum::<f32>() / n as f32;

    let half = n / 2;
    let mut out = Vec::with_capacity(half + 1);
    let inv_n = 1.0 / n as f32;

    for k in 0..=half {
        let mut re = 0.0f32;
        let mut im = 0.0f32;
        let omega = -2.0 * std::f32::consts::PI * k as f32 / n as f32;
        for (i, (_, v)) in probe.history.iter().enumerate() {
            let sample = *v - mean;
            let theta = omega * i as f32;
            re += sample * theta.cos();
            im += sample * theta.sin();
        }
        let magnitude = (re * re + im * im).sqrt() * inv_n;
        let freq = k as f32 / (n as f32 * dt);
        out.push((freq, magnitude));
    }

    out
}

/// Bevy system: sample every probe at the current simulation time and append
/// to its history ring buffer. Runs after `diagnostics_system` so the sampled
/// `DerivedFields` are fresh for this step.
pub fn probe_system(
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    mut probes: ResMut<ProbeSet>,
) {
    let Some(grid) = grid else { return };
    if diag.fields.is_empty() {
        return;
    }

    let max_history = probes.max_history;
    let t = grid.time as f32;

    for probe in &mut probes.probes {
        let x = probe.position[0].min(grid.nx.saturating_sub(1));
        let y = probe.position[1].min(grid.ny.saturating_sub(1));
        let z = probe.position[2].min(grid.nz.saturating_sub(1));
        let idx = grid.idx(x, y, z);
        let value = sample_probe_field(&grid, &diag, idx, probe.field);

        probe.history.push_back((t, value));
        while probe.history.len() > max_history {
            probe.history.pop_front();
        }
    }
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

    /// probe_fft should recover the frequency of a pure sinusoid within bin
    /// resolution after detrending.
    #[test]
    fn test_probe_fft_recovers_sinusoid() {
        use std::f32::consts::PI;

        let n = 512;
        let dt = 1.0e-3_f32; // 1 kHz sample rate
        let target_freq = 50.0_f32;

        let mut probe = Probe::new("sin", [0, 0, 0], ProbeField::SField);
        for i in 0..n {
            let t = i as f32 * dt;
            let v = (2.0 * PI * target_freq * t).sin();
            probe.history.push_back((t, v));
        }

        let spectrum = probe_fft(&probe);
        assert!(!spectrum.is_empty());

        // Peak should be at (or adjacent to) the target frequency bin.
        let df = spectrum[1].0 - spectrum[0].0;
        let (peak_freq, _) = spectrum
            .iter()
            .copied()
            .fold((0.0f32, 0.0f32), |acc, (f, m)| {
                if m > acc.1 { (f, m) } else { acc }
            });
        assert!(
            (peak_freq - target_freq).abs() <= df,
            "peak freq {} should be within one bin ({}) of {}",
            peak_freq,
            df,
            target_freq,
        );
    }

    /// probe_system should append one sample per call and cap history at max_history.
    #[test]
    fn test_probe_system_ring_buffer_cap() {
        let mut probe = Probe::new("p", [4, 4, 4], ProbeField::Phi);
        for i in 0..10 {
            probe.history.push_back((i as f32, i as f32));
        }
        let max_history = 6;
        while probe.history.len() > max_history {
            probe.history.pop_front();
        }
        assert_eq!(probe.history.len(), 6);
        assert_eq!(probe.history.front().unwrap().1, 4.0);
        assert_eq!(probe.history.back().unwrap().1, 9.0);
    }
}
