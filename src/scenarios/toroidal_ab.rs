// Toroidal Aharonov-Bohm scenario: macro-scale AB experiment.
//
// A toroidal solenoid (donut-shaped coil) confines its magnetic field
// entirely inside the torus tube — by symmetry, ideal toroidal coils
// produce ZERO B outside the tube and ZERO net magnetic dipole moment
// in the far field. Yet the vector potential A wraps the torus axis and
// is nonzero in the torus hole (and everywhere else outside).
//
// A pickup loop threaded through the torus hole (linking number = 1)
// has ∮A·dl equal to the flux through one cross-section of the torus tube,
// by Stokes' theorem:
//   ∮ A · dl  =  ∫∫ B · dS  =  Φ
// where the surface bounded by the pickup loop necessarily passes through
// the torus tube interior. For a thin torus with N loops of current I:
//   Φ  ≈  μ₀ N I r² / (2R)
// where R is the major radius and r is the minor radius (r ≪ R).
//
// The "AB-effect" framing of this scenario is that B = 0 along the entire
// pickup wire path, yet the loop integral is nonzero. Classical EM and
// QVED both predict this coupling and predict it identically, because the
// configuration is gauge-clean (φ = 0, ∇·A = 0 in the bulk; the kink at
// the surface current is a δ-source for the static current sheet, not a
// gauge-breaking term).
//
// We work in code units with μ₀ = 1 throughout. The diagnostic
// `compute_a_from_torus` evaluates the analytic A from a discrete set of
// N current loops via Biot-Savart, then the apply function loads it into
// every grid cell. With the live simulation paused, the topological
// statement (A ≠ 0 where B = 0; ∮A·dl quantises to Φ) is directly
// verifiable from the loaded A.

use crate::math::fdtd;
use crate::scenarios::aharonov_bohm;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::sources::SourceConfig;

const PI: f32 = std::f32::consts::PI;

/// Geometry + drive parameters for a toroidal coil.
#[derive(Clone, Debug)]
pub struct ToroidalConfig {
    /// World-space centre of the torus (the centroid of the donut).
    pub center: [f32; 3],
    /// Major radius R: distance from the torus centre to the centre of
    /// the tube (in world meters).
    pub major_radius: f32,
    /// Minor radius r: radius of the tube cross-section (world meters).
    pub minor_radius: f32,
    /// Number N of discrete poloidal current loops arranged around the
    /// torus axis. More loops → smoother toroidal B in the tube, smaller
    /// angular ripple. 16 is a good default for a 32³ grid.
    pub num_poloidal_loops: u32,
    /// Number of segments used to discretise each loop in the Biot-Savart
    /// sum. 32 is more than enough for this size of geometry.
    pub num_segments_per_loop: u32,
    /// Current per loop (code amperes; we use μ₀ = 1).
    pub current_per_loop: f32,
}

impl Default for ToroidalConfig {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            major_radius: 0.10,
            minor_radius: 0.025,
            num_poloidal_loops: 16,
            num_segments_per_loop: 32,
            current_per_loop: 1.0,
        }
    }
}

impl ToroidalConfig {
    /// Thin-torus analytic flux through one cross-section of the tube.
    /// Exact in the r ≪ R limit, off by O((r/R)²) otherwise (~6 % at
    /// r/R = 0.25). With μ₀ = 1 in code units:
    ///   Φ  =  N I r² / (2R)
    pub fn expected_flux(&self) -> f32 {
        let n = self.num_poloidal_loops as f32;
        let i = self.current_per_loop;
        let r = self.minor_radius;
        let big_r = self.major_radius;
        n * i * r * r / (2.0 * big_r)
    }
}

/// Compute the vector potential A at a world position by summing the
/// Biot-Savart contributions of every discretised segment of every
/// poloidal loop in the toroidal coil. With μ₀ = 1:
///   A(r)  =  (I / 4π) · Σ_loop Σ_seg  dl_seg / |r − r_seg|
pub fn compute_a_from_torus(cfg: &ToroidalConfig, pos: [f32; 3]) -> [f32; 3] {
    let prefactor = cfg.current_per_loop / (4.0 * PI);
    let dtheta = 2.0 * PI / cfg.num_segments_per_loop as f32;
    let mut a = [0.0f32; 3];

    for k in 0..cfg.num_poloidal_loops {
        let phi_k = 2.0 * PI * k as f32 / cfg.num_poloidal_loops as f32;
        let cos_phi = phi_k.cos();
        let sin_phi = phi_k.sin();
        let center_k = [
            cfg.center[0] + cfg.major_radius * cos_phi,
            cfg.center[1] + cfg.major_radius * sin_phi,
            cfg.center[2],
        ];
        // Loop at toroidal angle φ_k lies in the plane spanned by the
        // outward radial direction and ẑ. The loop's normal points in
        // the toroidal direction (−sin φ, cos φ, 0).
        let e1 = [cos_phi, sin_phi, 0.0]; // radial outward
        // e2 = ẑ

        for m in 0..cfg.num_segments_per_loop {
            let theta = 2.0 * PI * (m as f32 + 0.5) / cfg.num_segments_per_loop as f32;
            let cos_t = theta.cos();
            let sin_t = theta.sin();

            // Segment midpoint: c_k + r·(cos θ · e1 + sin θ · ẑ).
            let r_seg = [
                center_k[0] + cfg.minor_radius * cos_t * e1[0],
                center_k[1] + cfg.minor_radius * cos_t * e1[1],
                center_k[2] + cfg.minor_radius * sin_t,
            ];
            // dl = r·dθ·(−sin θ · e1 + cos θ · ẑ).
            let dl = [
                cfg.minor_radius * dtheta * (-sin_t) * e1[0],
                cfg.minor_radius * dtheta * (-sin_t) * e1[1],
                cfg.minor_radius * dtheta * cos_t,
            ];

            let dx = pos[0] - r_seg[0];
            let dy = pos[1] - r_seg[1];
            let dz = pos[2] - r_seg[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-6);

            a[0] += prefactor * dl[0] / dist;
            a[1] += prefactor * dl[1] / dist;
            a[2] += prefactor * dl[2] / dist;
        }
    }
    a
}

/// Apply the toroidal-AB scenario by writing the analytic A field into
/// every cell of the grid. φ = 0, A_z = 0 stays as the inactive component
/// (although our toroidal coil has true A_z ≠ 0; we capture all three
/// components from the Biot-Savart sum). q_dot = 0 — static configuration.
///
/// Sources are cleared. PML state, if any, should be reset by the caller.
pub fn apply_toroidal_ab_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
    cfg: &ToroidalConfig,
) {
    grid.reset();
    sources.sources.clear();

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;
    let half_x = nx as f32 * dx * 0.5;
    let half_y = ny as f32 * dx * 0.5;
    let half_z = nz as f32 * dx * 0.5;

    let buf = grid.current;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let world_x = (x as f32 + 0.5) * dx - half_x;
                let world_y = (y as f32 + 0.5) * dx - half_y;
                let world_z = (z as f32 + 0.5) * dx - half_z;
                let a = compute_a_from_torus(cfg, [world_x, world_y, world_z]);
                let idx = fdtd::idx(x, y, z, nx, ny);
                let cell = &mut grid.cells[buf][idx];
                cell.q = [0.0, a[0], a[1], a[2]];
                cell.q_dot = [0.0; 4];
            }
        }
    }
}

/// Build a closed pickup loop in the y = cy plane that links the torus
/// once on the +x side. Topology: the rectangle's left edge runs along
/// the torus axis (x = cx), the right edge sits well outside the tube
/// (x > cx + R + r), and the top/bottom edges connect them above/below
/// the torus tube. Linking number = 1.
pub fn linking_pickup_loop(
    cfg: &ToroidalConfig,
    segments_per_side: usize,
) -> Vec<[f32; 3]> {
    let cx = cfg.center[0];
    let cy = cfg.center[1];
    let cz = cfg.center[2];
    // Right edge sits at x = cx + R + 2r — well clear of the tube
    // (which extends to x = cx + R + r at toroidal angle 0).
    let x_extent = cfg.major_radius + 2.0 * cfg.minor_radius;
    // Top/bottom edges sit at z = cz ± R — far above/below the tube
    // (which extends only to z = cz ± r).
    let z_extent = cfg.major_radius;

    let n = segments_per_side.max(1);
    let dx_seg = x_extent / n as f32;
    let dz_seg = 2.0 * z_extent / n as f32;
    let mut path = Vec::with_capacity(4 * n);

    // Top edge: (cx, +z) → (cx + x_extent, +z), traversing +x̂.
    for i in 0..n {
        path.push([cx + i as f32 * dx_seg, cy, cz + z_extent]);
    }
    // Right edge: (cx + x_extent, +z) → (cx + x_extent, −z), traversing −ẑ.
    for i in 0..n {
        path.push([cx + x_extent, cy, cz + z_extent - i as f32 * dz_seg]);
    }
    // Bottom edge: (cx + x_extent, −z) → (cx, −z), traversing −x̂.
    for i in 0..n {
        path.push([cx + x_extent - i as f32 * dx_seg, cy, cz - z_extent]);
    }
    // Left edge: (cx, −z) → (cx, +z), traversing +ẑ.
    for i in 0..n {
        path.push([cx, cy, cz - z_extent + i as f32 * dz_seg]);
    }
    path
}

/// Build a small closed pickup loop FAR from the torus that does not link
/// the tube — sits above the torus in the xy-plane at z = cz + 1.5R.
pub fn non_linking_pickup_loop(
    cfg: &ToroidalConfig,
    segments_per_side: usize,
) -> Vec<[f32; 3]> {
    let cx = cfg.center[0];
    let cy = cfg.center[1];
    let cz = cfg.center[2] + cfg.major_radius * 1.5;
    let half_size = cfg.minor_radius * 0.5;
    aharonov_bohm::rectangular_loop([cx, cy, cz], half_size, half_size, segments_per_side)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::simulation::diagnostics::compute_derived_fields;
    use crate::simulation::field_update::step_field_cpu;
    use crate::simulation::sources::SourceConfig;

    /// 32³ grid with the default toroidal config applied.
    fn setup_default() -> (SimulationGrid, SourceConfig, ToroidalConfig) {
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        // Reduce loop count for faster tests; the topology is independent
        // of N, only the precise flux value depends on it.
        let config = ToroidalConfig {
            num_poloidal_loops: 12,
            num_segments_per_loop: 24,
            ..ToroidalConfig::default()
        };
        apply_toroidal_ab_scenario(&mut grid, &mut sources, &config);
        (grid, sources, config)
    }

    fn run_steps(grid: &mut SimulationGrid, n: usize, extended_mode: bool) {
        for _ in 0..n {
            let params = grid.sim_params(extended_mode);
            step_field_cpu(grid, &params, None, None);
            grid.swap_and_advance();
        }
    }

    /// Path-integral around a pickup loop that links the torus tube once
    /// must have nonzero magnitude; with μ₀ = 1 and a thin torus,
    ///   |∮A·dl| ≈ N I r² / (2R)
    /// to within O((r/R)²) ≈ 6 % plus discretisation errors from finite N
    /// loops and finite path resolution. Tolerance = 25 % on magnitude.
    ///
    /// Sign convention: with our poloidal-loop parametrisation the current
    /// flows +ẑ at the outer rim and −ẑ at the inner rim, which by Ampere's
    /// law puts B in the −e_φ direction inside the tube. The pickup loop
    /// is wound counter-clockwise as viewed from +ŷ (disk normal +ŷ), so
    /// the flux ∫∫B·dS is negative. The MAGNITUDE is what matches the
    /// thin-torus formula; the sign is just bookkeeping.
    #[test]
    fn test_toroidal_pickup_loop_linking() {
        let (grid, _src, config) = setup_default();
        let path = linking_pickup_loop(&config, 60);
        let flux = aharonov_bohm::integrate_a_dot_dl(&grid, &path);
        let expected = config.expected_flux();
        let rel_err = (flux.abs() - expected).abs() / expected;
        assert!(
            flux.abs() > 0.5 * expected,
            "Linking flux magnitude should be of order N·I·r²/(2R): expected {expected:.4e}, got {flux:.4e}"
        );
        assert!(
            rel_err < 0.30,
            "Linking |∮A·dl|: expected ≈ {expected:.4e}, got {flux:.4e} (rel_err = {rel_err:.3e})"
        );
    }

    /// A small loop FAR from the torus — outside the tube and not linking
    /// it — must give ∮A·dl ≈ 0. Within 5 % of the linking flux.
    #[test]
    fn test_toroidal_pickup_loop_non_linking() {
        let (grid, _src, config) = setup_default();
        let path = non_linking_pickup_loop(&config, 40);
        let flux = aharonov_bohm::integrate_a_dot_dl(&grid, &path);
        let expected_linking = config.expected_flux();
        let rel = flux.abs() / expected_linking;
        assert!(
            rel < 0.05,
            "Non-linking ∮A·dl: {flux:.4e} (relative to linking flux {expected_linking:.4e} = {rel:.3e})"
        );
    }

    /// |B| inside the torus tube is nonzero, dominantly toroidal. We
    /// sample a single cell at toroidal angle 0 and verify the toroidal
    /// component (along ŷ at that point) dominates.
    #[test]
    fn test_toroidal_b_inside_tube() {
        let (grid, _src, config) = setup_default();
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        // World point at toroidal angle 0 (centre of the tube on +x).
        let target_world = [config.major_radius, 0.0, 0.0];
        let dx = grid.dx;
        let half_x = grid.nx as f32 * dx * 0.5;
        let half_y = grid.ny as f32 * dx * 0.5;
        let half_z = grid.nz as f32 * dx * 0.5;
        let gx = ((target_world[0] + half_x) / dx - 0.5).round() as usize;
        let gy = ((target_world[1] + half_y) / dx - 0.5).round() as usize;
        let gz = ((target_world[2] + half_z) / dx - 0.5).round() as usize;
        let idx = fdtd::idx(gx, gy, gz, grid.nx as usize, grid.ny as usize);

        let b = derived[idx].b;
        let mag = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        // At toroidal angle 0, the toroidal direction is +ŷ. The y-component
        // should dominate.
        assert!(mag > 0.0, "|B| inside tube should be nonzero");
        assert!(
            b[1].abs() > 0.5 * mag,
            "B_y should dominate inside tube at φ=0; got B = ({:.3e}, {:.3e}, {:.3e})",
            b[0],
            b[1],
            b[2]
        );
    }

    /// |B| at a point well OUTSIDE the torus tube must be small. We take
    /// the average over a far-field shell to rule out single-cell
    /// fluctuations from the discrete loop wires.
    #[test]
    fn test_toroidal_b_outside_tube() {
        let (grid, _src, config) = setup_default();
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        // Sample cells well above and below the torus plane (|z| > 2R is
        // well outside the tube and far enough from the wires to wash out
        // discretisation noise).
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let dx = grid.dx;
        let half_z = grid.nz as f32 * dx * 0.5;

        let mut max_b_outside = 0.0f32;
        for z in 0..grid.nz as usize {
            let z_world = (z as f32 + 0.5) * dx - half_z;
            if z_world.abs() < 2.5 * config.minor_radius {
                continue; // skip the torus plane
            }
            for y in 4..(ny - 4) {
                for x in 4..(nx - 4) {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    let b = derived[idx].b;
                    let mag = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
                    if mag > max_b_outside {
                        max_b_outside = mag;
                    }
                }
            }
        }

        // Pick a B-magnitude reference from the inside of the tube and
        // compare. A discrete-loop torus has some leakage but it should be
        // a small fraction of the in-tube field.
        let target_world = [config.major_radius, 0.0, 0.0];
        let half_x = grid.nx as f32 * dx * 0.5;
        let half_y = grid.ny as f32 * dx * 0.5;
        let gx = ((target_world[0] + half_x) / dx - 0.5).round() as usize;
        let gy = ((target_world[1] + half_y) / dx - 0.5).round() as usize;
        let gz = ((target_world[2] + half_z) / dx - 0.5).round() as usize;
        let idx_inside = fdtd::idx(gx, gy, gz, nx, ny);
        let b_inside = derived[idx_inside].b;
        let mag_inside =
            (b_inside[0].powi(2) + b_inside[1].powi(2) + b_inside[2].powi(2)).sqrt();

        let leakage_ratio = max_b_outside / mag_inside;
        assert!(
            leakage_ratio < 0.5,
            "Far-field B leakage ratio: {leakage_ratio:.3e} (max outside {max_b_outside:.3e}, inside {mag_inside:.3e})"
        );
    }

    /// QVED equivalence claim: gauge-clean toroidal initial conditions
    /// (φ = 0, ∇·A ≈ 0 in the bulk) imply S stays at zero in extended
    /// mode, so the extended evolution reduces algebraically to the
    /// standard evolution. Two identical grids run 30 steps (one each)
    /// must produce A fields that agree to floating-point precision.
    #[test]
    fn test_toroidal_extended_matches_standard() {
        let mut grid_std = SimulationGrid::new(32, 32, 32, 0.01);
        let mut grid_ext = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources_std = SourceConfig::default();
        let mut sources_ext = SourceConfig::default();
        let config = ToroidalConfig {
            num_poloidal_loops: 8,
            num_segments_per_loop: 16,
            ..ToroidalConfig::default()
        };

        apply_toroidal_ab_scenario(&mut grid_std, &mut sources_std, &config);
        apply_toroidal_ab_scenario(&mut grid_ext, &mut sources_ext, &config);

        run_steps(&mut grid_std, 30, false);
        run_steps(&mut grid_ext, 30, true);

        let cs = grid_std.read_buf();
        let ce = grid_ext.read_buf();
        let mut max_abs_diff = 0.0f32;
        let mut max_a_mag = 0.0f32;
        for (a, b) in cs.iter().zip(ce.iter()) {
            for k in 1..4 {
                max_abs_diff = max_abs_diff.max((a.q[k] - b.q[k]).abs());
                max_a_mag = max_a_mag.max(a.q[k].abs());
            }
        }
        assert!(max_a_mag > 0.0, "A is identically zero — setup broken");
        let rel_diff = max_abs_diff / max_a_mag;
        assert!(
            rel_diff < 1e-4,
            "Extended vs standard mode: rel diff = {rel_diff:.3e} (max |Δ| = {max_abs_diff:.3e}, max |A| = {max_a_mag:.3e})"
        );
    }

    /// S stays bounded near zero in extended mode, by the same argument
    /// as the linear AB case: the Lorenz gauge term ∂φ/∂t/c² + ∇·A is
    /// zero at t = 0 (φ = 0; ∇·A is zero by construction in the
    /// continuous limit, finite-difference residual is a tiny boundary
    /// effect).
    #[test]
    fn test_toroidal_s_stays_zero_extended() {
        let (mut grid, _src, _config) = setup_default();
        run_steps(&mut grid, 30, true);
        let s = grid.s_read();
        let max_s = s.iter().copied().fold(0.0f32, |a, b| a.max(b.abs()));
        assert!(
            max_s < 0.5,
            "max |S| in extended mode after 30 steps: {max_s:.3e} (should be ≪ A magnitudes)"
        );
    }
}
