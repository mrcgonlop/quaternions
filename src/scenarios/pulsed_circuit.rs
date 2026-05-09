// Pulsed-circuit scalar-mode excitation (Tier-3 stub).
//
// **Background.** Sharply interrupted currents — particularly on time
// scales below the local L/R relaxation time — should excite the QVED
// scalar mode S more efficiently than steady or slowly-varying drive.
// The classical signature is the ringing impedance spike on a circuit
// breaker; the QVED prediction is an additional radiation channel into
// S that conventional Maxwell theory can't carry.
//
// **What this scenario tests.** Whether a step-function current
// interruption produces a measurable scalar-wave pulse at probes placed
// outside the circuit, with amplitude and time-dependence distinct from
// the classical EM near-field that always accompanies di/dt events.
//
// **Geometry.** A simple straight-wire current path through the grid
// driven by a switchable source:
//   * Wire along x̂ centred on the grid, length ≪ grid extent.
//   * Current I(t) = I₀ for t < t_switch, then drops to 0 with rise
//     time τ ≪ L/c (where L is the wire length).
//   * Probes positioned away from the wire to record S, |E|, |B|.
//
// **Diagnostic.** Time-series of S, |E|, |B| at the probes. The classical
// EM signal (dipole-radiation-like burst from the current interruption)
// shows in |E| and |B|; the QVED-specific signal shows in S only when
// extended_mode is enabled.
//
// **Expected observables.**
//   * Standard EM: probes see a transverse-EM pulse (E ⊥ propagation,
//     B ⊥ propagation), no S signal.
//   * QVED: same transverse-EM pulse plus a longitudinal S pulse that
//     propagates at speed c with amplitude ∝ d²I/dt² at the source.
//   * Pulse delay = (probe distance) / c, identical for both transverse
//     and longitudinal channels.
//
// **Parameter ranges to explore.**
//   * Rise time τ from 0.1 to 10 cells / c.
//   * I₀ from 1e-3 to 1.0 in code units.
//   * Wire length from 4 to 16 cells.
//   * Probe distance from 4 to 16 cells.

use crate::simulation::grid::SimulationGrid;
use crate::simulation::sources::SourceConfig;

/// Geometry + drive parameters for the pulsed-circuit scenario.
#[derive(Clone, Debug)]
pub struct PulsedCircuitConfig {
    /// World-space centre of the circuit.
    pub center: [f32; 3],
    /// Wire length (world meters).
    pub wire_length: f32,
    /// Steady-state current amplitude before switch (code units).
    pub initial_current: f32,
    /// Time at which the current interrupts (seconds, sim time).
    pub switch_time: f32,
    /// Rise time of the interruption — small τ produces large dI/dt
    /// and a stronger scalar-mode excitation.
    pub rise_time: f32,
}

impl Default for PulsedCircuitConfig {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            wire_length: 0.08,
            initial_current: 1.0,
            switch_time: 1.0e-10,
            rise_time: 1.0e-12,
        }
    }
}

/// Apply the pulsed-circuit scenario (NOT IMPLEMENTED — Tier-3 stub).
///
/// PSEUDOCODE:
///   1. Reset grid, clear sources.
///   2. Add a `CurrentPulse` source along the wire: amplitude follows a
///      step function dropping at `switch_time` over `rise_time`.
///   3. Add probes at multiple distances along ŷ to record S, |E|, |B|.
///   4. Run for ~3× wire-length-over-c so the propagating pulse reaches
///      all probes.
///   5. Compare probe time-series in standard vs extended mode — the S
///      channel should be silent in standard, populated in extended.
///
/// Implementation requires:
///   * A SwitchedSource SourceType (or extension of CurrentPulse) that
///     models a step-function current.
///   * A multi-probe array configured automatically by the scenario.
pub fn apply_pulsed_circuit_scenario(
    _grid: &mut SimulationGrid,
    _sources: &mut SourceConfig,
    _cfg: &PulsedCircuitConfig,
) {
    todo!(
        "Pulsed-circuit scalar-mode scenario (Phase 7.4 stub) — needs a \
         switched-current source variant and an automatic multi-probe \
         install. See module docstring for the full pseudocode."
    );
}
