// VacuumK scenario: demonstrate polarizable vacuum K-field evolution.
//
// Places a dense oscillating charge cluster at the grid center with
// vacuum_config.enabled = true. The strong E field excites virtual e⁺e⁻
// pairs, raising K above 1 in the high-field region. When K > √2, the
// Weber force sign flips for like charges at sub-luminal velocities,
// enabling self-sustaining EVO (Exotic Vacuum Object) formation.
//
// Expected behavior:
//   - K rises above 1 in the cluster region (visible in KVacuum slice)
//   - Cluster remains stable longer than in K=1 (standard vacuum)
//   - In extended mode, the K pocket also modulates S propagation speed

use crate::simulation::grid::SimulationGrid;
use crate::simulation::plugin::VacuumConfig;
use crate::simulation::sources::{Source, SourceConfig};
use crate::simulation::state::SimParams;

/// Apply the VacuumK scenario to an existing (reset) grid.
///
/// Sets up a dense oscillating source at center to drive K field evolution,
/// and configures VacuumConfig with physically motivated parameters.
///
/// `omega_p_scale`: fraction of (c/dx) to use as ωₚ. At 0.01, K evolves
/// on a timescale ~100 dt steps — observable within a short run.
pub fn apply_vacuum_k_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
    vacuum_config: &mut VacuumConfig,
) {
    grid.reset();
    sources.sources.clear();

    let c0 = SimParams::C0;
    let dx = grid.dx;

    // Oscillating dipole cluster at center, high amplitude to drive significant K rise.
    // Frequency: wavelength ≈ 10 cells (higher than standard dipole for stronger E field).
    let center = grid.nx as f32 / 2.0;
    let freq = c0 / (10.0 * dx);
    let amplitude = 10.0; // 10× standard amplitude for strong local field

    sources.sources.push(Source::dipole_z(
        [center, center, center],
        amplitude,
        freq,
    ));

    // VacuumConfig: use ωₚ ≈ c/(100·dx) so K evolves visibly within ~100 steps.
    // eta = 1e-4 (QED Euler-Heisenberg coupling), u_s = 1.0 (code units).
    let omega_p = c0 / (100.0 * dx);
    vacuum_config.enabled = true;
    vacuum_config.omega_p = omega_p;
    vacuum_config.eta = 1e-4;
    vacuum_config.u_s = 1.0;
}
