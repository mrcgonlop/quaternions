// Brown / Biefeld asymmetric capacitor in a polarisable vacuum (Tier-3 stub).
//
// **Background.** T. T. Brown reported in the 1920s–60s that asymmetric
// high-voltage capacitors (one large electrode, one small electrode,
// separated by a dielectric) experienced a small unidirectional thrust
// when charged to several kV. Mainstream explanation: ion wind from corona
// discharge at the small electrode plus electrostatic asymmetric pressure.
// Speculative QVED interpretation: the strong field asymmetry polarises
// the vacuum (K(x) > 1 in the high-field region), creating a K-gradient
// pressure on the device.
//
// **What this scenario tests.** Whether the (Q, K) coupled system, with
// VacuumConfig parameters scaled to make the polarisable-vacuum response
// observable, predicts a net force on the device. NOT a verification of
// Brown's claims — there is no reproducible experimental evidence at the
// magnitude claimed by Brown's ultra-high-voltage demonstrations. The
// simulation produces a phenomenological prediction that one can compare
// against experimental upper bounds (e.g., Tajmar's NASA-funded
// replications giving null results at ~1 % of the claimed thrust).
//
// **Geometry.** Asymmetric parallel-plate capacitor with axes along ẑ:
//   * Large electrode: a wide disk at z = 0 (radius R_large).
//   * Small electrode: a narrow disk at z = d (radius R_small ≪ R_large).
//   * Voltage V applied between them (set as initial conditions on φ at
//     the conductor cells, hold via `CellFlags::CONDUCTOR`).
//
// **Diagnostic.** Net axial force on the device = ∫(stress tensor) · ẑ
// over a closed surface enclosing both electrodes. With dynamic K (Phase
// 1.8 already gives this in the K field via VacuumConfig), the K-gradient
// contribution adds to the standard Maxwell stress.
//
// **Expected observables.**
//   * K(x) > 1 between the plates, peaking near the small electrode where
//     |E| is highest.
//   * Net force has a component aligned with the small-electrode axis.
//   * The force scales sub-linearly with V (saturating as K saturates).
//
// **Parameter ranges to explore.**
//   * V from 1 to 100 kV (in scaled code units).
//   * R_large / R_small from 2 to 20.
//   * d from 0.5 to 5 cm (relative to grid resolution).
//   * VacuumConfig.eta scaling factor — sweep to find when the K-gradient
//     contribution rises above the standard Maxwell-stress baseline.

use crate::simulation::grid::SimulationGrid;
use crate::simulation::plugin::VacuumConfig;
use crate::simulation::sources::SourceConfig;

/// Geometry parameters for the asymmetric capacitor.
#[derive(Clone, Debug)]
pub struct BrownConfig {
    /// World-space centre of the device (midplane between electrodes).
    pub center: [f32; 3],
    /// Plate separation along ẑ (world meters).
    pub plate_separation: f32,
    /// Large-electrode radius (world meters).
    pub large_radius: f32,
    /// Small-electrode radius (world meters).
    pub small_radius: f32,
    /// Applied voltage between plates (code units).
    pub voltage: f32,
}

impl Default for BrownConfig {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            plate_separation: 0.04,
            large_radius: 0.06,
            small_radius: 0.012,
            voltage: 1.0,
        }
    }
}

/// Apply the Brown capacitor scenario (NOT IMPLEMENTED — Tier-3 stub).
///
/// PSEUDOCODE:
///   1. Reset grid, clear sources.
///   2. Mark conductor cells: a disk of radius `large_radius` at z = -d/2,
///      a disk of radius `small_radius` at z = +d/2 (cell flags
///      `CellFlags::CONDUCTOR`).
///   3. Initialise φ = +V/2 on the small electrode, φ = −V/2 on the large
///      electrode, φ = 0 elsewhere (linear ramp inside the gap as a
///      first-iterate; let the field solver relax to the correct boundary
///      conditions over a few hundred steps).
///   4. Enable VacuumConfig (Phase 1.8 K dynamics) so the polarisable
///      vacuum responds to the field energy density.
///   5. Run until K relaxes to a quasi-steady configuration.
///   6. Compute the net axial force via Maxwell stress + K-gradient
///      pressure on a closed surface enclosing the device.
///
/// Implementation requires:
///   * A boundary-condition update inside `step_field_cpu` that holds φ
///     fixed on conductor cells (currently unimplemented).
///   * A force-integration diagnostic on a user-specified surface.
pub fn apply_brown_capacitor_scenario(
    _grid: &mut SimulationGrid,
    _sources: &mut SourceConfig,
    _vacuum: &mut VacuumConfig,
    _cfg: &BrownConfig,
) {
    todo!(
        "Brown capacitor scenario (Phase 7.4 stub) — needs conductor BC \
         enforcement in step_field_cpu and a Maxwell-stress force diagnostic. \
         See module docstring for the full pseudocode and parameter ranges."
    );
}
