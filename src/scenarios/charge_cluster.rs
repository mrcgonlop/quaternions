// Charge-cluster (EVO / Shoulders' "exotic vacuum object") stability search
// (Tier-3 stub).
//
// **Background.** Kenneth Shoulders reported in the 1980s–90s that high-
// voltage breakdown experiments produced micron-scale charge clusters
// containing ~10¹¹ electrons in dynamically stable configurations,
// surviving despite Coulomb repulsion that should disperse them in
// femtoseconds. Mainstream physics has no accepted explanation; the
// observations remain disputed.
//
// **What this scenario tests.** Whether the (Q, K) coupled system with
// strong vacuum polarisation in the cluster interior provides enough
// confining pressure to stabilise a like-charge cluster against Coulomb
// expansion. Per README §3.4, K > √2 inside the cluster reverses the
// effective Weber force between drifting like charges from repulsive to
// attractive — providing a candidate confinement mechanism.
//
// **Geometry.** N like-charged particles initialised in a small spherical
// volume with random thermal velocities. K initialised either at vacuum
// (K = 1, see if K self-organises into a stable shell) or with a pre-
// formed shell (K_inside > √2, see if it persists).
//
// **Diagnostic.**
//   * RMS cluster radius vs time — a stable cluster maintains finite
//     radius; an unstable one expands without bound.
//   * Total charge (should be conserved exactly).
//   * K(r) profile — does a self-consistent K shell form?
//   * Per-particle Weber force history — does the v²/c² · K² correction
//     shift to attractive at observed drift speeds?
//
// **Expected observables (if Shoulders' claims have any QVED basis).**
//   * Below a critical η (VacuumConfig coupling), cluster expands ballistically.
//   * Above critical η, K self-organises into a shell with K_max > √2,
//     and the cluster radius oscillates around a finite value.
//
// **Parameter ranges to explore.**
//   * N from 8 to 1024 particles (limited by N² Weber complexity).
//   * Initial cluster radius from 1 to 5 cells.
//   * Initial drift speed from 0 to 0.5 c.
//   * VacuumConfig.eta from 0 to 1.0 (search for the critical value).

use crate::simulation::particles::ParticleSystem;
use crate::simulation::plugin::VacuumConfig;
use crate::simulation::weber::ForceMode;

/// Initialisation parameters for a charge-cluster stability run.
#[derive(Clone, Debug)]
pub struct ChargeClusterConfig {
    /// World-space centre of the initial cluster.
    pub center: [f32; 3],
    /// Number of charged particles to seed.
    pub num_particles: u32,
    /// Initial cluster radius (world meters).
    pub initial_radius: f32,
    /// Charge per particle (code units; like-charged for repulsion test).
    pub charge_per_particle: f32,
    /// Mass per particle (code units; lighter → faster expansion / orbit).
    pub mass_per_particle: f32,
    /// Initial RMS thermal velocity per particle (m/s, code units).
    pub initial_thermal_speed: f32,
    /// Random seed for reproducible initial conditions.
    pub seed: u64,
}

impl Default for ChargeClusterConfig {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            num_particles: 64,
            initial_radius: 0.01,
            charge_per_particle: 1.0,
            mass_per_particle: 1.0,
            initial_thermal_speed: 1.0e6,
            // "C1A5_5EED" — pun-encoded reproducibility marker.
            seed: 0xC1A5_5EED_C1A5_5EED_u64,
        }
    }
}

/// Apply the charge-cluster scenario (NOT IMPLEMENTED — Tier-3 stub).
///
/// PSEUDOCODE:
///   1. Clear `ParticleSystem`.
///   2. Sample N random positions uniformly in a sphere of radius
///      `initial_radius` around `center`.
///   3. Sample N random velocities from an isotropic Maxwell distribution
///      with rms speed `initial_thermal_speed`.
///   4. Push N particles with the given charge and mass.
///   5. Set `*force_mode = ForceMode::Weber` (or `Both` to include any
///      grid-coupled Lorentz contribution from external fields).
///   6. Enable VacuumConfig with elevated eta so the cluster's field
///      energy can polarise vacuum strongly enough for K > √2.
///   7. Run for many cluster-traversal times; record RMS radius and K(r).
///
/// Implementation requires:
///   * A seeded random number generator (deterministic for reproducible
///     stability searches).
///   * A cluster-radius-vs-time diagnostic resource and UI plot.
///   * A K(r) radial-profile diagnostic.
pub fn apply_charge_cluster_scenario(
    _particles: &mut ParticleSystem,
    _force_mode: &mut ForceMode,
    _vacuum: &mut VacuumConfig,
    _cfg: &ChargeClusterConfig,
) {
    todo!(
        "Charge-cluster (EVO) stability scenario (Phase 7.4 stub) — needs a \
         seeded RNG, RMS-radius and K(r) diagnostics, and an elevated-η \
         VacuumConfig regime. See module docstring for the full pseudocode \
         and parameter ranges."
    );
}
