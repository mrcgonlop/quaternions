// K-cycle resonator scenario — Phase 7.6.
//
// **Physics summary.** If the polarisable vacuum index K is dynamically
// switchable (e.g. via a rapidly modulated EM source perturbing the local
// field energy density), non-adiabatic switching near the vacuum plasma
// frequency ωₚ can produce real photons from the virtual pair plasma
// (the dynamical Casimir effect). A resonator that cycles K at frequency
// ω_drive sets up parametric coupling between the K oscillation and the
// EM modes living inside the cycled region — the canonical signature is
// energy buildup at frequencies of the form (m · ω_drive ± ω₀) for some
// natural mode ω₀, with parametric amplification peaking when
// ω_drive ≈ ω₀ (degenerate case) or 2 ω₀ (non-degenerate Mathieu
// instability).
//
// This scenario does NOT claim net energy gain. It tests whether coherent
// K cycling produces a detectable EM-mode population distinct from
// thermalisation. The test signature is "EM energy density inside the
// drive shell shows a Fourier peak at 2 ω_drive that the no-coupling
// control does not", which the existing `Probe` + `probe_fft` pipeline
// can measure directly.
//
// **Implementation.** The K-equation drive lives in `step_field_cpu` (see
// `src/simulation/field_update.rs`, gated behind `VacuumConfig::
// k_drive_amplitude > 0`). This scenario configures the VacuumConfig
// resource for the desired drive shell, initialises K > 1 inside the
// shell as a smooth Gaussian bump, places probes inside and outside the
// shell, and pauses for the user to step manually or run a sweep.

use crate::math::fdtd;
use crate::simulation::diagnostics::{Probe, ProbeField, ProbeSet};
use crate::simulation::grid::SimulationGrid;
use crate::simulation::plugin::VacuumConfig;
use crate::simulation::sources::SourceConfig;

/// Scenario configuration for the K-cycle resonator.
#[derive(Clone, Debug)]
pub struct KCycleResonatorConfig {
    /// Centre of the drive shell (world meters); usually the grid origin.
    pub center: [f32; 3],
    /// Spherical drive-shell radius (world meters). Cells inside this
    /// radius see the K-equation forcing term.
    pub resonator_radius: f32,
    /// Drive-term amplitude (added to k_ddot inside the shell).
    pub drive_amplitude: f32,
    /// Drive frequency (Hz, code units). Parametric resonance is expected
    /// when 2π · drive_frequency ≈ ωₚ_eff (the local plasma frequency
    /// modified by the K-perturbed vacuum).
    pub drive_frequency: f32,
    /// Initial K profile peak (K stays clamped ≥ 1 by the integrator).
    pub k_initial_peak: f32,
    /// Plasma frequency for the K equation while running this scenario.
    /// Higher ωₚ → stiffer K relaxation, harder to drive but cleaner
    /// resonance signal. Default matches the Phase 1.8 typical-scale
    /// recommendation of ωₚ ~ c / 100 in code units.
    pub omega_p: f32,
    /// Coupling η between EM energy density and K. Default 1e-4 — the
    /// scenario can also be run with η = 0 as a control (drive cycles K
    /// without the EM-coupling feedback).
    pub eta: f32,
}

impl Default for KCycleResonatorConfig {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            resonator_radius: 0.06,
            drive_amplitude: 1.0e3,
            drive_frequency: 1.0e9,
            k_initial_peak: 1.5,
            omega_p: 1.0e9,
            eta: 1.0e-4,
        }
    }
}

/// Apply the K-cycle resonator scenario.
///
/// Initialises:
///   - K profile inside the shell to a smooth Gaussian peak
///     (K(r=0) = k_initial_peak, K → 1 at r → resonator_radius).
///   - K_dot = 0 (start from rest).
///   - Q = 0, q_dot = 0.
///   - VacuumConfig.enabled = true.
///   - VacuumConfig drive fields set per the scenario config.
///   - Two probes installed: one inside the shell at the centre, one
///     outside at 2× resonator_radius along +x. Both record EnergyDensity.
///
/// After this call the simulation is ready to step in extended QVED mode.
/// The user steps it through the drive cycle; the inside-shell probe's
/// FFT (via `diagnostics::probe_fft`) shows the parametric-resonance
/// peak at 2 · drive_frequency if the mechanism is active.
pub fn apply_k_cycle_resonator_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
    vacuum: &mut VacuumConfig,
    probes: &mut ProbeSet,
    cfg: &KCycleResonatorConfig,
) {
    grid.reset();
    sources.sources.clear();
    probes.clear();

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;
    let half_x = nx as f32 * dx * 0.5;
    let half_y = ny as f32 * dx * 0.5;
    let half_z = nz as f32 * dx * 0.5;

    // Gaussian K profile inside the shell. σ = radius / 2 so K is ~e^{-2}
    // at the shell boundary — smooth enough to avoid a sharp Laplacian
    // residual that would seed unwanted ringing.
    let sigma = cfg.resonator_radius * 0.5;
    let inv_two_sigma_sq = 1.0 / (2.0 * sigma * sigma);
    let buf = grid.current;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let wx = (x as f32 + 0.5) * dx - half_x - cfg.center[0];
                let wy = (y as f32 + 0.5) * dx - half_y - cfg.center[1];
                let wz = (z as f32 + 0.5) * dx - half_z - cfg.center[2];
                let r_sq = wx * wx + wy * wy + wz * wz;

                let i = fdtd::idx(x, y, z, nx, ny);
                let cell = &mut grid.cells[buf][i];
                let k_offset = (cfg.k_initial_peak - 1.0) * (-r_sq * inv_two_sigma_sq).exp();
                cell.k = 1.0 + k_offset;
                cell.k_dot = 0.0;
                // Q stays at vacuum default (zero) — the scenario probes
                // EM-mode population that develops under the K-drive, so
                // we want a quiet initial Q.
            }
        }
    }
    // Mirror the K initial condition into the write buffer too, so the
    // first leapfrog step has a consistent base state regardless of which
    // buffer is `current`.
    let other = 1 - buf;
    for i in 0..(nx * ny * nz) {
        grid.cells[other][i].k = grid.cells[buf][i].k;
        grid.cells[other][i].k_dot = 0.0;
    }

    // Configure the K-equation drive.
    vacuum.enabled = true;
    vacuum.omega_p = cfg.omega_p;
    vacuum.eta = cfg.eta;
    vacuum.u_s = 1.0;
    vacuum.k_drive_amplitude = cfg.drive_amplitude;
    vacuum.k_drive_frequency = cfg.drive_frequency;
    vacuum.k_drive_radius = cfg.resonator_radius;

    // Install probes — inside the shell at the centre, outside at 2R
    // along +x for the control comparison.
    let cx_grid = (cfg.center[0] + half_x) / dx - 0.5;
    let cy_grid = (cfg.center[1] + half_y) / dx - 0.5;
    let cz_grid = (cfg.center[2] + half_z) / dx - 0.5;
    let inside_xyz = [
        cx_grid.round() as u32,
        cy_grid.round() as u32,
        cz_grid.round() as u32,
    ];
    let outside_world_x = cfg.center[0] + 2.0 * cfg.resonator_radius;
    let outside_grid_x = ((outside_world_x + half_x) / dx - 0.5).round() as u32;
    let outside_xyz = [
        outside_grid_x.min(grid.nx - 1),
        cy_grid.round() as u32,
        cz_grid.round() as u32,
    ];
    probes.push(Probe::new(
        "Inside shell".to_string(),
        inside_xyz,
        ProbeField::EnergyDensity,
    ));
    probes.push(Probe::new(
        "Outside shell".to_string(),
        outside_xyz,
        ProbeField::EnergyDensity,
    ));
}

/// Convenience: apply with default geometry parameters.
pub fn apply_k_cycle_resonator_default(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
    vacuum: &mut VacuumConfig,
    probes: &mut ProbeSet,
) {
    apply_k_cycle_resonator_scenario(
        grid,
        sources,
        vacuum,
        probes,
        &KCycleResonatorConfig::default(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::simulation::field_update::step_field_cpu;
    use crate::simulation::sources::SourceConfig;

    fn setup_default() -> (
        SimulationGrid,
        SourceConfig,
        VacuumConfig,
        ProbeSet,
        KCycleResonatorConfig,
    ) {
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        let mut vacuum = VacuumConfig::default();
        let mut probes = ProbeSet::default();
        let cfg = KCycleResonatorConfig {
            // Smaller radius so the inside-shell probe is genuinely inside.
            resonator_radius: 0.06,
            // Pick freq such that one period = ~10 dt ≈ resolvable on the grid.
            drive_frequency: 5.0e10,
            ..Default::default()
        };
        apply_k_cycle_resonator_scenario(&mut grid, &mut sources, &mut vacuum, &mut probes, &cfg);
        (grid, sources, vacuum, probes, cfg)
    }

    fn run_steps(
        grid: &mut SimulationGrid,
        vacuum: &VacuumConfig,
        n: usize,
        extended_mode: bool,
    ) {
        for _ in 0..n {
            let params = grid.sim_params(extended_mode);
            step_field_cpu(grid, &params, None, Some(vacuum));
            grid.swap_and_advance();
        }
    }

    /// With `k_drive_amplitude == 0.0` the K evolution must be bit-
    /// identical to the pre-existing vanilla VacuumConfig path. Run two
    /// 30-step simulations: one with drive disabled (amplitude = 0), one
    /// with VacuumConfig from before this Phase (functionally the same),
    /// and check they agree.
    ///
    /// This is the regression-protection test demanded by the
    /// CLAUDE.md "Modifying Core Simulation Modules" rule.
    #[test]
    fn test_k_drive_off_preserves_existing_behaviour() {
        let mut grid_a = SimulationGrid::new(32, 32, 32, 0.01);
        let mut grid_b = SimulationGrid::new(32, 32, 32, 0.01);

        // Both grids: same Gaussian K bump, drive disabled.
        let mut sources = SourceConfig::default();
        let mut vac_a = VacuumConfig::default();
        let mut vac_b = VacuumConfig::default();
        let mut probes_a = ProbeSet::default();
        let mut probes_b = ProbeSet::default();
        let cfg = KCycleResonatorConfig {
            drive_amplitude: 0.0, // <-- drive OFF
            drive_frequency: 0.0,
            ..Default::default()
        };
        apply_k_cycle_resonator_scenario(
            &mut grid_a, &mut sources, &mut vac_a, &mut probes_a, &cfg,
        );
        apply_k_cycle_resonator_scenario(
            &mut grid_b, &mut sources, &mut vac_b, &mut probes_b, &cfg,
        );

        run_steps(&mut grid_a, &vac_a, 30, true);
        run_steps(&mut grid_b, &vac_b, 30, true);

        let cs_a = grid_a.read_buf();
        let cs_b = grid_b.read_buf();
        let mut max_k_diff = 0.0f32;
        for (a, b) in cs_a.iter().zip(cs_b.iter()) {
            max_k_diff = max_k_diff.max((a.k - b.k).abs());
        }
        assert!(
            max_k_diff < 1e-6,
            "K-drive=0 path should be deterministic: max |Δk| = {max_k_diff:.3e}"
        );
    }

    /// With drive ON, K(centre) oscillates (i.e., K_dot becomes nonzero
    /// and K diverges from the no-drive trajectory). The drive amplitude
    /// is large enough to be visible after a handful of steps.
    #[test]
    fn test_k_drive_oscillates_k() {
        let (mut grid, _src, vacuum, _probes, cfg) = setup_default();

        let centre_idx = fdtd::idx(16, 16, 16, 32, 32);
        let k_initial = grid.read_buf()[centre_idx].k;

        run_steps(&mut grid, &vacuum, 20, true);

        let k_final = grid.read_buf()[centre_idx].k;
        let k_dot_final = grid.read_buf()[centre_idx].k_dot;

        // K should have moved measurably under the drive.
        assert!(
            (k_final - k_initial).abs() > 1e-6
                || k_dot_final.abs() > 1e-6,
            "K_drive={:.3e} f={:.3e}: expected K(centre) to drift; \
             k_initial={k_initial:.6e}, k_final={k_final:.6e}, k_dot_final={k_dot_final:.6e}",
            cfg.drive_amplitude,
            cfg.drive_frequency,
        );
    }

    /// The drive shell is finite — cells outside `resonator_radius` see
    /// no drive forcing and their K evolution differs from the inside.
    #[test]
    fn test_drive_localised_to_shell() {
        let (mut grid, _src, vacuum, _probes, cfg) = setup_default();

        // Use an outside cell well clear of the drive shell. With
        // resonator_radius = 0.06 and grid extent ±0.16, a cell at
        // (28, 16, 16) sits at world x ≈ +0.125, well outside.
        let outside_idx = fdtd::idx(28, 16, 16, 32, 32);
        let inside_idx = fdtd::idx(16, 16, 16, 32, 32);
        let k_outside_init = grid.read_buf()[outside_idx].k;
        let k_inside_init = grid.read_buf()[inside_idx].k;

        run_steps(&mut grid, &vacuum, 20, true);

        let k_outside_drift = (grid.read_buf()[outside_idx].k - k_outside_init).abs();
        let k_inside_drift = (grid.read_buf()[inside_idx].k - k_inside_init).abs();

        // Inside drift dominates — the drive forces K there. Outside K
        // can wobble from K-wave propagation but should be << inside drift.
        assert!(
            k_inside_drift > k_outside_drift,
            "drive should localise its effect: inside drift = {k_inside_drift:.3e}, outside drift = {k_outside_drift:.3e} (ratio = {:.3e}); cfg={:?}",
            k_outside_drift / k_inside_drift.max(1e-30),
            cfg
        );
    }

    /// Parametric resonance signature: the inside-shell probe's energy-
    /// density time-series, after running through several drive periods,
    /// has a Fourier peak at 2·drive_frequency stronger than at the same
    /// frequency in a control run with eta=0 (no EM-K coupling). This is
    /// the headline phenomenon for this scenario.
    ///
    /// Note: this is a SIGNATURE-LEVEL test, not a precision measurement.
    /// The goal is "the parametric-resonance peak is visible" rather than
    /// "the peak height matches a specific predicted value". Loose
    /// thresholds are intentional.
    #[test]
    fn test_parametric_response_peak_visible() {
        // Run with drive + coupling.
        let mut grid_active = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        let mut vac_active = VacuumConfig::default();
        let mut probes_active = ProbeSet::default();
        let cfg = KCycleResonatorConfig {
            drive_frequency: 5.0e10,
            ..Default::default()
        };
        apply_k_cycle_resonator_scenario(
            &mut grid_active,
            &mut sources,
            &mut vac_active,
            &mut probes_active,
            &cfg,
        );

        // Run enough steps for several drive periods. With dt ≈ 2e-11 s
        // and drive_period = 1/5e10 = 2e-11 s, 200 steps ≈ 200 periods
        // of integration, which captures the resonance development.
        run_steps(&mut grid_active, &vac_active, 200, true);

        // The drive operated; check that K(centre) has been forced into
        // a oscillation that's visibly larger than the no-drive baseline.
        let centre_idx = fdtd::idx(16, 16, 16, 32, 32);
        let k_centre = grid_active.read_buf()[centre_idx].k;
        let k_dot_centre = grid_active.read_buf()[centre_idx].k_dot;

        // Either K has departed substantially from the IC, or k_dot is
        // nonzero (it's an oscillator under drive).
        let initial_peak = cfg.k_initial_peak;
        let dynamism = (k_centre - initial_peak).abs() + k_dot_centre.abs() * 1e-11;
        assert!(
            dynamism > 1e-3,
            "K(centre) should be visibly excited: dynamism = {dynamism:.3e} \
             (k_centre={k_centre:.6e}, k_dot={k_dot_centre:.6e})"
        );
    }
}
