// Aharonov-Bohm scenario: solenoid with B confined inside, A extending outside.
//
// Physics: the AB effect is the canonical demonstration that the
// electromagnetic potentials carry physical information that the fields E, B
// alone do not. A long solenoid with a steady current produces:
//   * Inside (ρ < R):  B = B₀ ẑ, A_θ = B₀ ρ / 2     (uniform B, linear A)
//   * Outside (ρ ≥ R): B = 0,    A_θ = B₀ R² / (2ρ)  (zero B, but A ≠ 0)
//
// Stokes' theorem gives the topological hallmark:
//   ∮ A · dl = ∫∫ B · dS  =  π R² B₀  ≡  Φ
// for ANY closed path enclosing the solenoid, and 0 for any path that does
// not. The value of the loop integral depends only on the topology of the
// path relative to the flux tube, not on the path's shape — that is the
// Aharonov-Bohm effect, classically.
//
// Quantum mechanics adds the observation that the wavefunction phase shift
// `Δϕ = (q/ħ) ∮ A · dl` modulates an electron-interferometer pattern. We
// can't simulate the quantum interference here (this is a classical FDTD),
// but we CAN demonstrate the underlying classical gauge-potential structure
// — A ≠ 0 where B = 0, and the loop integral quantises to Φ.
//
// "Satisfied in both formulations":
//   The analytic A above is divergence-free everywhere except for the
//   surface current at ρ = R. With φ ≡ 0 and ∇·A = 0, the QVED scalar
//   field S = (1/c²)∂φ/∂t + ∇·A is identically zero at t = 0 with zero
//   time-derivative. Since □S = −ρ/ε₀ and there's no charge density
//   in this configuration, S stays zero forever, and the extended-mode
//   evolution reduces EXACTLY to the standard-mode evolution. The AB
//   loop integral is therefore preserved identically in both formulations.

use crate::math::fdtd;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::sources::SourceConfig;

/// Scenario configuration for an Aharonov-Bohm solenoid.
#[derive(Clone, Debug)]
pub struct AbConfig {
    /// Solenoid axis position in fractional grid-cell coordinates (xy plane).
    /// `cell + 0.5` is the cell centre, so `[15.5, 15.5]` puts the axis at
    /// world (0, 0) for a 32×32×32 grid.
    pub axis_xy_cells: [f32; 2],
    /// Solenoid radius in cell widths. Should be ≥ 3 to resolve the inside,
    /// ≤ ~nx/3 to leave room for an exterior region.
    pub radius_cells: f32,
    /// Magnetic field amplitude inside the solenoid (code units).
    pub b0: f32,
}

impl Default for AbConfig {
    fn default() -> Self {
        Self {
            axis_xy_cells: [15.5, 15.5], // grid centre for 32³
            radius_cells: 6.0,
            b0: 1.0,
        }
    }
}

impl AbConfig {
    /// Expected magnetic flux through the solenoid (= ∮A·dl for any
    /// encircling path) in world units.
    pub fn expected_flux(&self, grid: &SimulationGrid) -> f32 {
        let r_world = self.radius_cells * grid.dx;
        std::f32::consts::PI * r_world * r_world * self.b0
    }
}

/// Apply the Aharonov-Bohm solenoid scenario to a grid.
///
/// Initialises Q at every cell to the analytic solenoid potential:
///   φ = 0                 (no scalar potential)
///   A_x = -B₀·dy/2        for ρ < R   (cylindrical inside)
///   A_y = +B₀·dx/2        for ρ < R
///   A_x = -B₀R²·dy/(2ρ²)  for ρ ≥ R   (1/ρ outside)
///   A_y = +B₀R²·dx/(2ρ²)  for ρ ≥ R
///   A_z = 0
/// where (dx, dy) are world-space offsets from the solenoid axis.
///
/// Time-derivatives are zero — the configuration is static. ∇·A is zero by
/// construction except at the boundary ρ = R, so S stays at zero in
/// extended mode (see module header).
///
/// Existing source list is cleared; the AB configuration is driven only
/// by the initial A field, with no on-going injection.
pub fn apply_ab_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
    config: &AbConfig,
) {
    grid.reset();
    sources.sources.clear();

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;

    let cx = config.axis_xy_cells[0];
    let cy = config.axis_xy_cells[1];
    let r_world = config.radius_cells * dx;
    let r2_world = r_world * r_world;
    let b0 = config.b0;

    let buf = grid.current;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                // World-space offsets from the solenoid axis. Cell centres
                // sit at (cell + 0.5) in fractional grid coordinates.
                let off_x = ((x as f32 + 0.5) - cx) * dx;
                let off_y = ((y as f32 + 0.5) - cy) * dx;
                let rho2 = off_x * off_x + off_y * off_y;

                let (a_x, a_y) = if rho2 < r2_world {
                    // Inside: A = (-B₀·y/2, +B₀·x/2, 0). Curl gives B_z = B₀.
                    (-0.5 * b0 * off_y, 0.5 * b0 * off_x)
                } else {
                    // Outside: A = (-B₀R²·y/(2ρ²), +B₀R²·x/(2ρ²), 0). Curl
                    // gives B = 0; A is the harmonic 1/ρ continuation.
                    let factor = 0.5 * b0 * r2_world / rho2.max(1e-30);
                    (-factor * off_y, factor * off_x)
                };

                let idx = fdtd::idx(x, y, z, nx, ny);
                let cell = &mut grid.cells[buf][idx];
                cell.q = [0.0, a_x, a_y, 0.0];
                cell.q_dot = [0.0; 4];
            }
        }
    }
}

/// Sample the vector potential A trilinearly at a world-space point.
/// Returns `[0, 0, 0]` if the point is outside the safe interpolation
/// region (≤ 1 cell from any face).
pub fn sample_a_at_world(grid: &SimulationGrid, world_pos: [f32; 3]) -> [f32; 3] {
    let nx = grid.nx as i32;
    let ny = grid.ny as i32;
    let nz = grid.nz as i32;
    let dx = grid.dx;
    let half_x = nx as f32 * dx * 0.5;
    let half_y = ny as f32 * dx * 0.5;
    let half_z = nz as f32 * dx * 0.5;

    let gx = (world_pos[0] + half_x) / dx - 0.5;
    let gy = (world_pos[1] + half_y) / dx - 0.5;
    let gz = (world_pos[2] + half_z) / dx - 0.5;

    let i = gx.floor() as i32;
    let j = gy.floor() as i32;
    let k = gz.floor() as i32;

    if i < 0 || j < 0 || k < 0 || i + 1 >= nx || j + 1 >= ny || k + 1 >= nz {
        return [0.0; 3];
    }

    let fx = gx - i as f32;
    let fy = gy - j as f32;
    let fz = gz - k as f32;
    let nx_us = nx as usize;
    let ny_us = ny as usize;
    let cells = grid.read_buf();

    let mut a = [0.0f32; 3];
    for di in 0..2 {
        for dj in 0..2 {
            for dk in 0..2 {
                let ix = (i + di) as usize;
                let iy = (j + dj) as usize;
                let iz = (k + dk) as usize;
                let idx = fdtd::idx(ix, iy, iz, nx_us, ny_us);
                let wx = if di == 0 { 1.0 - fx } else { fx };
                let wy = if dj == 0 { 1.0 - fy } else { fy };
                let wz = if dk == 0 { 1.0 - fz } else { fz };
                let w = wx * wy * wz;
                a[0] += w * cells[idx].q[1];
                a[1] += w * cells[idx].q[2];
                a[2] += w * cells[idx].q[3];
            }
        }
    }
    a
}

/// Numerically integrate `∮ A · dl` along a closed polyline in world space.
///
/// `path` is a list of vertices; the loop is closed by adding a final
/// segment from `path.last()` back to `path[0]`. Each segment contributes
/// `A(midpoint) · (p_{i+1} − p_i)` (midpoint rule, 2nd-order accurate).
pub fn integrate_a_dot_dl(grid: &SimulationGrid, path: &[[f32; 3]]) -> f32 {
    if path.len() < 2 {
        return 0.0;
    }
    let mut sum = 0.0f32;
    for i in 0..path.len() {
        let p0 = path[i];
        let p1 = path[(i + 1) % path.len()];
        let mid = [
            0.5 * (p0[0] + p1[0]),
            0.5 * (p0[1] + p1[1]),
            0.5 * (p0[2] + p1[2]),
        ];
        let a = sample_a_at_world(grid, mid);
        let dl = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        sum += a[0] * dl[0] + a[1] * dl[1] + a[2] * dl[2];
    }
    sum
}

/// Build a closed rectangular loop of `4 × segments_per_side` vertices in
/// the xy-plane, centred on `centre = [cx, cy, cz]`, with half-widths
/// `half_x`, `half_y`. Wound counter-clockwise (positive flux through +ẑ).
pub fn rectangular_loop(
    centre: [f32; 3],
    half_x: f32,
    half_y: f32,
    segments_per_side: usize,
) -> Vec<[f32; 3]> {
    let n = segments_per_side.max(1);
    let cx = centre[0];
    let cy = centre[1];
    let cz = centre[2];
    let mut path = Vec::with_capacity(4 * n);

    let dx_seg = 2.0 * half_x / n as f32;
    let dy_seg = 2.0 * half_y / n as f32;

    // Right edge: (cx + half_x, cy − half_y → cy + half_y)
    for i in 0..n {
        path.push([cx + half_x, cy - half_y + i as f32 * dy_seg, cz]);
    }
    // Top edge: (cx + half_x → cx − half_x, cy + half_y)
    for i in 0..n {
        path.push([cx + half_x - i as f32 * dx_seg, cy + half_y, cz]);
    }
    // Left edge
    for i in 0..n {
        path.push([cx - half_x, cy + half_y - i as f32 * dy_seg, cz]);
    }
    // Bottom edge
    for i in 0..n {
        path.push([cx - half_x + i as f32 * dx_seg, cy - half_y, cz]);
    }
    path
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::simulation::diagnostics::compute_derived_fields;
    use crate::simulation::field_update::step_field_cpu;
    use crate::simulation::sources::SourceConfig;

    /// Convenience: 32³ grid with the default AB config applied.
    fn setup_default() -> (SimulationGrid, SourceConfig, AbConfig) {
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        let config = AbConfig::default();
        apply_ab_scenario(&mut grid, &mut sources, &config);
        (grid, sources, config)
    }

    /// Run `n` simulation steps with no sources, no PML, no vacuum config.
    /// `step_field_cpu` is the same kernel that runs in the live simulation.
    fn run_steps(grid: &mut SimulationGrid, n: usize, extended_mode: bool) {
        for _ in 0..n {
            let params = grid.sim_params(extended_mode);
            step_field_cpu(grid, &params, None, None);
            grid.swap_and_advance();
        }
    }

    /// B_z averaged over cells well inside the solenoid (ρ < R/2) should
    /// be ≈ B₀. Cells at ρ ≈ R suffer from the surface-current discontinuity
    /// in finite differences, so we sample only the deep interior.
    #[test]
    fn test_ab_b_field_inside() {
        let (grid, _src, config) = setup_default();
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        let cx = config.axis_xy_cells[0];
        let cy = config.axis_xy_cells[1];
        let r_inner = config.radius_cells * 0.5;

        let mut bz_sum = 0.0f32;
        let mut count = 0usize;
        for z in 4..28 {
            for y in 0..32 {
                for x in 0..32 {
                    let rx = (x as f32 + 0.5) - cx;
                    let ry = (y as f32 + 0.5) - cy;
                    if rx * rx + ry * ry < r_inner * r_inner {
                        let idx = fdtd::idx(x, y, z, 32, 32);
                        bz_sum += derived[idx].b[2];
                        count += 1;
                    }
                }
            }
        }
        assert!(count > 0, "no cells in the deep interior");
        let bz_avg = bz_sum / count as f32;
        let rel_err = (bz_avg - config.b0).abs() / config.b0;
        assert!(
            rel_err < 0.02,
            "B_z inside: expected {}, got {bz_avg} (rel_err = {rel_err:.3e}, n = {count})",
            config.b0
        );
    }

    /// Outside the solenoid (ρ > 1.5 R), |B| should be very small. We
    /// allow a few percent of B₀ as numerical residual from the boundary.
    #[test]
    fn test_ab_b_field_outside() {
        let (grid, _src, config) = setup_default();
        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);

        let cx = config.axis_xy_cells[0];
        let cy = config.axis_xy_cells[1];
        let r_outer = config.radius_cells * 1.5;

        let mut max_b = 0.0f32;
        for z in 4..28 {
            for y in 1..31 {
                for x in 1..31 {
                    let rx = (x as f32 + 0.5) - cx;
                    let ry = (y as f32 + 0.5) - cy;
                    if rx * rx + ry * ry > r_outer * r_outer {
                        let idx = fdtd::idx(x, y, z, 32, 32);
                        let b = derived[idx].b;
                        let mag = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
                        if mag > max_b {
                            max_b = mag;
                        }
                    }
                }
            }
        }
        assert!(
            max_b < 0.05 * config.b0,
            "max |B| outside: {max_b} (should be ≪ {})",
            config.b0
        );
    }

    /// |A| at a deep-exterior point must be NONZERO — that is the AB
    /// hallmark, A ≠ 0 where B = 0.
    #[test]
    fn test_ab_a_nonzero_outside() {
        let (grid, _src, config) = setup_default();

        let r_world = config.radius_cells * grid.dx;
        // World point at ρ ≈ 2R from axis (axis is at world (0, 0) for the
        // default 32³ config with cx=cy=15.5).
        let pos = [2.0 * r_world, 0.0, 0.0];
        let a = sample_a_at_world(&grid, pos);
        let mag = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();

        // Analytic value at ρ = 2R: A_θ = B₀R²/(2·2R) = B₀R/4.
        let expected = config.b0 * r_world / 4.0;
        let rel_err = (mag - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "|A| at ρ=2R: expected {expected}, got {mag} (rel_err = {rel_err:.3e})"
        );
    }

    /// ∮A·dl on a square loop encircling the solenoid should equal the
    /// flux Φ = πR²B₀ (Stokes' theorem). Tolerance allows for trilinear
    /// interpolation and finite-segment discretisation error.
    #[test]
    fn test_ab_path_integral_encircling() {
        let (grid, _src, config) = setup_default();

        // Square path well outside the solenoid radius (R = 6 cells = 0.06 m;
        // half-side 0.10 m = 10 cells). Centre is the solenoid axis at world
        // (0, 0) for cx = cy = 15.5.
        let path = rectangular_loop([0.0, 0.0, 0.0], 0.10, 0.10, 80);
        let flux = integrate_a_dot_dl(&grid, &path);
        let expected = config.expected_flux(&grid);
        let rel_err = (flux - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "∮A·dl encircling: expected {expected:.4e}, got {flux:.4e} (rel_err = {rel_err:.3e})"
        );
    }

    /// ∮A·dl on a square loop FAR from the solenoid (no enclosed flux)
    /// should vanish — within the same tolerance scaled by the expected
    /// encircling flux.
    #[test]
    fn test_ab_path_integral_non_encircling() {
        let (grid, _src, config) = setup_default();

        // Small loop far in +x, well clear of the flux tube.
        // Solenoid is centred at (0, 0); push the loop centre to (+0.10, 0)
        // with half-side 0.02 m — entirely outside the safe interpolation
        // region near the +x face? Let me check: 32³ at dx=0.01 spans
        // [-0.16, +0.16]. Centre at x=0.10 with half-side 0.02 → loop
        // occupies x ∈ [0.08, 0.12], well within the grid.
        let path = rectangular_loop([0.10, 0.0, 0.0], 0.02, 0.02, 40);
        let flux = integrate_a_dot_dl(&grid, &path);
        let expected_encircling = config.expected_flux(&grid);
        let rel_err = flux.abs() / expected_encircling;
        assert!(
            rel_err < 0.02,
            "∮A·dl non-encircling: {flux:.4e} (relative to encircling flux = {rel_err:.3e})"
        );
    }

    /// Short-time flux preservation: the analytic AB potential is
    /// divergence-free with ∇²A = 0 in the bulk, but it has a kink at
    /// ρ = R that produces a δ-function-like Laplacian under finite
    /// differences. That kink seeds an outgoing wave at speed c. As long
    /// as the wave hasn't propagated to the test loop, the flux is
    /// preserved; we verify this for ~5 steps before the wave reaches
    /// the loop. Long-time drift is real physics — a static AB
    /// configuration without a sustaining current is not a steady-state
    /// of the homogeneous wave equation.
    #[test]
    fn test_ab_short_time_flux_preservation() {
        let (mut grid, _src, config) = setup_default();
        let path = rectangular_loop([0.0, 0.0, 0.0], 0.10, 0.10, 80);
        let flux_initial = integrate_a_dot_dl(&grid, &path);

        run_steps(&mut grid, 5, false);

        let flux_final = integrate_a_dot_dl(&grid, &path);
        let expected = config.expected_flux(&grid);
        let rel_drift = (flux_final - flux_initial).abs() / expected;
        assert!(
            rel_drift < 0.05,
            "Flux drift after 5 steps: |Δ|/Φ = {rel_drift:.3e} (initial {flux_initial:.4e}, final {flux_final:.4e}, Φ = {expected:.4e})"
        );
    }

    /// Standard ↔ extended mode equivalence — the central AB-in-quaternion
    /// claim. With φ=0 and ∇·A=0 at t=0, S stays identically zero in
    /// extended mode (verified by `test_ab_s_field_stays_zero_extended`),
    /// so the extended-mode update collapses algebraically to the
    /// standard-mode update. Both runs of the same number of steps must
    /// yield bitwise-similar A fields, modulo floating-point noise from
    /// the slightly different code paths.
    ///
    /// This is what "AB satisfied in both formulations" means precisely:
    /// the QVED extension introduces no new physics beyond conventional EM
    /// for a gauge-clean configuration.
    #[test]
    fn test_ab_extended_matches_standard() {
        let mut grid_std = SimulationGrid::new(32, 32, 32, 0.01);
        let mut grid_ext = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources_std = SourceConfig::default();
        let mut sources_ext = SourceConfig::default();
        let config = AbConfig::default();

        apply_ab_scenario(&mut grid_std, &mut sources_std, &config);
        apply_ab_scenario(&mut grid_ext, &mut sources_ext, &config);

        run_steps(&mut grid_std, 50, false);
        run_steps(&mut grid_ext, 50, true);

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
        // Extremely tight: the two code paths execute the same algebraic
        // operations on A when S = 0, so drift should be at f32 epsilon
        // accumulated over 50 steps. Tolerance allows for that without
        // hiding any real disagreement.
        assert!(
            rel_diff < 1e-4,
            "Extended vs standard mode after 50 steps: rel diff = {rel_diff:.3e} (max |Δ| = {max_abs_diff:.3e}, max |A| = {max_a_mag:.3e})"
        );
    }

    /// Extended-mode S stays bounded near zero — gauge-clean initial
    /// conditions (∇·A = 0, φ = 0) give S ≡ 0 at t = 0 with zero
    /// time-derivative. With no sources, S obeys □S = 0 and stays zero.
    /// Numerical noise from the boundary discontinuity may seed a tiny
    /// S signal; the assertion is that it stays small.
    #[test]
    fn test_ab_s_field_stays_zero_extended() {
        let (mut grid, _src, _config) = setup_default();
        run_steps(&mut grid, 50, true);

        let s = grid.s_read();
        let max_s = s.iter().copied().fold(0.0f32, |a, b| a.max(b.abs()));
        // S should be orders of magnitude below B₀ = 1.0. A loose threshold
        // (any value < 0.1) catches the catastrophic case where extended
        // mode disrupts the AB structure; tighter values may flag the
        // expected boundary-ringing on a coarse 32³ grid.
        assert!(
            max_s < 0.1,
            "max |S| in extended mode after 50 steps: {max_s:.3e} (should be ≪ B₀ = 1)"
        );
    }
}
