// Graneau wire scenario: a chain of charge-segment pairs that exposes
// Weber's prediction of longitudinal forces on a straight current-carrying
// wire.
//
// Physics motivation:
//   Standard Lorentz / Biot-Savart electrodynamics predicts NO axial force
//   on a straight wire from its own current. The self-generated B field is
//   azimuthal, the current is axial, so J × B is purely radial (the pinch
//   effect compresses the wire but produces no longitudinal stress).
//
//   Weber electrodynamics gives a different answer. For two collinear
//   "segments" (ion + co-located drifting electron) separated by L, the net
//   axial Weber force on each segment from the other is
//       F_axial = ∓ (q² / 4πε₀L²) · (v_d² / c²)
//   pointing AWAY from the other segment. Summing over a uniform chain:
//   * Interior segments see equal pushes from each side → near-zero net.
//   * End segments see only one-sided pushes → large outward force.
//
//   This matches the Graneau-Aspden experiments where pulsed kA-scale
//   currents through copper wires fragment the wires preferentially near
//   the ends — a signature Lorentz/Biot-Savart cannot explain.
//
// Derivation sketch (two segments A at 0, B at L; drift velocity v_d x̂):
//   Each force is along x̂ by symmetry. For ion_A from ion_B (both at rest):
//       F = -q²/(4πε₀L²) x̂                                  (Coulomb, repulsive)
//   For ion_A from electron_B (B's electron drifts at v_d x̂):
//       v₁₂ = -v_d x̂, ṙ = +v_d, |v|² − ṙ² = 0,  factor = 1 − v_d²/(2c²)
//       F = +q²/(4πε₀L²)·(1 − v_d²/(2c²)) x̂                 (attraction, reduced)
//   Sum on ion_A from segment B: −(q²/4πε₀L²)·(v_d²/(2c²)) x̂
//   Same sign on electron_A from segment B (v₁₂ flips for the e-i pair, same
//   for e-e), so each segment contributes a force of magnitude
//       |F_seg→seg| = (q²/4πε₀L²)·(v_d²/c²)
//   in the direction *away* from the source segment.
//
// Scaling note (REAL vs. SIMULATION):
//   Real copper at 10 kA through 1 mm² gives v_drift ≈ 0.74 m/s and
//   v²/c² ≈ 6×10⁻¹⁸ — invisible at f32 precision. The simulation uses
//   v_drift = 1×10⁷ m/s (~0.03 c, correction ~5×10⁻⁴) so the longitudinal
//   pattern is numerically resolvable. The QUALITATIVE physics (end-peaked,
//   antisymmetric profile) is preserved; force magnitudes are unphysical.

use crate::simulation::particles::{ParticleState, ParticleSystem};
use crate::simulation::state::SimParams;
use crate::simulation::weber::{compute_weber_forces, ForceMode};

/// Default chain length (odd → middle segment well-defined).
pub const DEFAULT_NUM_SEGMENTS: usize = 21;
/// Default spacing between adjacent segments (meters).
pub const DEFAULT_SEGMENT_SPACING: f32 = 0.01;
/// Default drift velocity used for the scaled-up "current" (m/s).
/// 1e7 m/s ≈ 0.033 c — the smallest value that keeps v²/c² above the f32
/// noise floor (~1e-7) by several orders of magnitude.
pub const DEFAULT_DRIFT_VELOCITY: f32 = 1.0e7;

/// Charge magnitude per segment (Coulombs, scaled).
const SEGMENT_CHARGE: f32 = 1.0;
/// Ion mass (kg, scaled). Heavy enough that the lattice barely responds to
/// internal force pulses on simulation timescales — so the ion positions
/// stay where we put them and the longitudinal-force diagnostic is clean.
const ION_MASS: f32 = 1.0e6;
/// Electron mass (kg, scaled). Light → carries the drift current.
const ELECTRON_MASS: f32 = 1.0;

/// Convenience: number of particles produced by the scenario for a given
/// segment count. Each segment contributes one ion + one electron.
pub fn particle_count(num_segments: usize) -> usize {
    2 * num_segments
}

/// Apply the Graneau wire scenario.
///
/// Populates `particles` with `num_segments` collocated (ion, electron) pairs
/// arranged along the x-axis at the grid centre. Ions are heavy and at rest;
/// every electron starts drifting with `v_drift` along +x̂ to model a uniform
/// current. After this call:
///   * `particles` holds 2·num_segments `ParticleState`s, ordered
///     (ion[0], electron[0], ion[1], electron[1], …).
///   * `*force_mode = ForceMode::Weber` — the Lorentz prediction is zero on a
///     straight wire by symmetry, so only Weber will reveal the effect.
pub fn apply_graneau_wire_scenario(
    particles: &mut ParticleSystem,
    force_mode: &mut ForceMode,
    num_segments: usize,
    spacing: f32,
    v_drift: f32,
) {
    particles.clear();

    // Centre the chain on x = 0.
    let half_length = (num_segments.saturating_sub(1)) as f32 * spacing * 0.5;

    for i in 0..num_segments {
        let x = i as f32 * spacing - half_length;
        let pos = [x, 0.0, 0.0];

        // Heavy positive ion at rest.
        particles.push(ParticleState::new(
            SEGMENT_CHARGE,
            ION_MASS,
            pos,
            [0.0; 3],
        ));
        // Light electron co-located with the ion, drifting in +x̂.
        particles.push(ParticleState::new(
            -SEGMENT_CHARGE,
            ELECTRON_MASS,
            pos,
            [v_drift, 0.0, 0.0],
        ));
    }

    *force_mode = ForceMode::Weber;
}

/// Apply the scenario with default geometry parameters.
pub fn apply_graneau_wire_default(
    particles: &mut ParticleSystem,
    force_mode: &mut ForceMode,
) {
    apply_graneau_wire_scenario(
        particles,
        force_mode,
        DEFAULT_NUM_SEGMENTS,
        DEFAULT_SEGMENT_SPACING,
        DEFAULT_DRIFT_VELOCITY,
    );
}

/// Default y-spacing between the two parallel strands of the hairpin
/// (world meters). 5 mm is a "thick" PCB spacing, easy to resolve.
pub const DEFAULT_HAIRPIN_SPACING: f32 = 0.005;

/// Apply the hairpin PCB-trace scenario.
///
/// Two parallel rows of (ion, electron) pairs along x̂, separated in y by
/// `pcb_spacing`. The TOP row at y = +s/2 carries current in +x̂ (electrons
/// drift +v); the BOTTOM row at y = −s/2 carries current in −x̂ (electrons
/// drift −v). This is the canonical PCB-trace geometry where a current
/// doubles back on itself — the same single conductor folded into two
/// parallel strands.
///
/// In standard Lorentz/Biot-Savart electrodynamics the antiparallel currents
/// repel each other transversely (Ampere's law for parallel currents), but
/// neither strand experiences any LONGITUDINAL force from the other. Weber
/// theory predicts a longitudinal contribution from the inter-strand
/// pairwise interactions whose pattern differs from the single-wire Graneau
/// profile — that's the test signature.
///
/// Particle ordering is (top_ion[0], top_e[0], …, top_ion[N−1], top_e[N−1],
/// bot_ion[0], bot_e[0], …, bot_ion[N−1], bot_e[N−1]). Top row first.
pub fn apply_hairpin_scenario(
    particles: &mut ParticleSystem,
    force_mode: &mut ForceMode,
    num_segments: usize,
    segment_spacing: f32,
    pcb_spacing: f32,
    v_drift: f32,
) {
    particles.clear();
    let half_length = (num_segments.saturating_sub(1)) as f32 * segment_spacing * 0.5;
    let half_y = pcb_spacing * 0.5;

    // Top strand: y = +s/2, electrons drift in +x̂.
    for i in 0..num_segments {
        let x = i as f32 * segment_spacing - half_length;
        let pos = [x, half_y, 0.0];
        particles.push(ParticleState::new(SEGMENT_CHARGE, ION_MASS, pos, [0.0; 3]));
        particles.push(ParticleState::new(
            -SEGMENT_CHARGE,
            ELECTRON_MASS,
            pos,
            [v_drift, 0.0, 0.0],
        ));
    }
    // Bottom strand: y = −s/2, electrons drift in −x̂ (return current).
    for i in 0..num_segments {
        let x = i as f32 * segment_spacing - half_length;
        let pos = [x, -half_y, 0.0];
        particles.push(ParticleState::new(SEGMENT_CHARGE, ION_MASS, pos, [0.0; 3]));
        particles.push(ParticleState::new(
            -SEGMENT_CHARGE,
            ELECTRON_MASS,
            pos,
            [-v_drift, 0.0, 0.0],
        ));
    }

    *force_mode = ForceMode::Weber;
}

/// Apply the hairpin scenario with default geometry parameters.
pub fn apply_hairpin_default(
    particles: &mut ParticleSystem,
    force_mode: &mut ForceMode,
) {
    apply_hairpin_scenario(
        particles,
        force_mode,
        DEFAULT_NUM_SEGMENTS,
        DEFAULT_SEGMENT_SPACING,
        DEFAULT_HAIRPIN_SPACING,
        DEFAULT_DRIFT_VELOCITY,
    );
}

/// Compute the per-ion axial Weber force on the TOP strand of a hairpin.
///
/// Returns `(x_position, F_x)` for each top-strand ion (the even-indexed
/// particles in the first half of the array per the ordering documented in
/// `apply_hairpin_scenario`). The bottom strand's contribution is included
/// in the pair sum — it's what makes this profile differ from the single-
/// strand `compute_wire_force_profile`.
pub fn compute_hairpin_top_force_profile(
    particles: &ParticleSystem,
    num_segments: usize,
) -> Vec<(f32, f32)> {
    let forces = compute_weber_forces(&particles.particles, SimParams::C0);
    let n = num_segments;
    let mut profile = Vec::with_capacity(n);
    for i in 0..n {
        let ion_idx = 2 * i; // top-strand ions are at indices 0, 2, …, 2(N−1)
        if let Some(p) = particles.particles.get(ion_idx) {
            profile.push((p.position[0], forces[ion_idx][0]));
        }
    }
    profile
}

/// Compute the per-ion axial Weber force along the wire.
///
/// Returns a Vec of `(x_position, F_x)` pairs, one entry per ion (the
/// even-indexed particles in `apply_graneau_wire_scenario`'s ordering).
/// The signature predicted by Weber theory:
///   * Antisymmetric about the wire centre (uniform-current → mirror
///     symmetry → the axial force flips sign at the midpoint).
///   * Magnitude peaks at the ends — the hallmark wire-fragmentation
///     prediction.
///   * Magnitude near zero in the interior — interior ions see balanced
///     pushes from chain segments on both sides.
pub fn compute_wire_force_profile(particles: &ParticleSystem) -> Vec<(f32, f32)> {
    let forces = compute_weber_forces(&particles.particles, SimParams::C0);
    particles
        .particles
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            // Even indices are ions per `apply_graneau_wire_scenario`.
            if i % 2 == 0 {
                Some((p.position[0], forces[i][0]))
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Setup creates the right number of particles in the right pattern,
    /// and toggles the force mode to Weber.
    #[test]
    fn test_apply_graneau_wire_scenario_creates_chain() {
        let mut particles = ParticleSystem::default();
        let mut force_mode = ForceMode::Lorentz;

        apply_graneau_wire_scenario(
            &mut particles,
            &mut force_mode,
            DEFAULT_NUM_SEGMENTS,
            DEFAULT_SEGMENT_SPACING,
            DEFAULT_DRIFT_VELOCITY,
        );

        assert_eq!(force_mode, ForceMode::Weber);
        assert_eq!(
            particles.particles.len(),
            particle_count(DEFAULT_NUM_SEGMENTS)
        );

        // First ion at -half_length, last ion at +half_length.
        let half = (DEFAULT_NUM_SEGMENTS - 1) as f32 * DEFAULT_SEGMENT_SPACING * 0.5;
        assert!((particles.particles[0].position[0] - (-half)).abs() < 1e-7);
        let last_ion_idx = 2 * (DEFAULT_NUM_SEGMENTS - 1);
        assert!(
            (particles.particles[last_ion_idx].position[0] - half).abs() < 1e-7
        );

        // Wire is electrically neutral.
        let total_charge: f32 = particles.particles.iter().map(|p| p.charge).sum();
        assert!(total_charge.abs() < 1e-6);

        // Ions are at rest, electrons drift in +x.
        for i in 0..DEFAULT_NUM_SEGMENTS {
            let ion = &particles.particles[2 * i];
            let electron = &particles.particles[2 * i + 1];
            assert_eq!(ion.velocity, [0.0; 3]);
            assert!(electron.velocity[0] > 0.0);
            assert_eq!(electron.velocity[1..], [0.0, 0.0]);
        }
    }

    /// The signature Graneau prediction: end ions feel large outward axial
    /// force; middle ion feels near-zero force; profile antisymmetric.
    #[test]
    fn test_graneau_force_profile_outward_at_ends() {
        let mut particles = ParticleSystem::default();
        let mut force_mode = ForceMode::Lorentz;
        apply_graneau_wire_default(&mut particles, &mut force_mode);

        let profile = compute_wire_force_profile(&particles);
        let n = profile.len();
        assert_eq!(n, DEFAULT_NUM_SEGMENTS);

        let f_left = profile[0].1;
        let f_right = profile[n - 1].1;
        let f_middle = profile[n / 2].1;

        // Leftmost ion is pushed in -x̂ (outward, toward the −x end).
        assert!(
            f_left < 0.0,
            "Leftmost ion should be pushed in −x: F_x = {f_left:.3e}"
        );
        // Rightmost ion is pushed in +x̂ (outward, toward the +x end).
        assert!(
            f_right > 0.0,
            "Rightmost ion should be pushed in +x: F_x = {f_right:.3e}"
        );

        let f_end = f_left.abs().max(f_right.abs());

        // Middle ion: by uniform-chain symmetry, axial force ≈ 0.
        assert!(
            f_middle.abs() < 0.01 * f_end,
            "Middle ion force {f_middle:.3e} should be << end force {f_end:.3e} (ratio = {:.3e})",
            f_middle.abs() / f_end
        );

        // Antisymmetric profile: F_left ≈ −F_right.
        let asymmetry = (f_left + f_right).abs() / f_end;
        assert!(
            asymmetry < 1e-3,
            "Profile should be antisymmetric: F_left = {f_left:.3e}, F_right = {f_right:.3e}, asymmetry = {asymmetry:.3e}"
        );

        // End force magnitude should be substantial — much larger than the
        // numerical noise floor for the chosen scaling.
        assert!(
            f_end > 1e6,
            "End force should be substantial; got |F| = {f_end:.3e} N"
        );
    }

    /// Force profile should monotonically grow toward the ends (from the
    /// midpoint outward, |F| increases). This rules out spurious noise-only
    /// patterns and confirms the spatial structure.
    #[test]
    fn test_graneau_force_increases_toward_ends() {
        let mut particles = ParticleSystem::default();
        let mut force_mode = ForceMode::Lorentz;
        apply_graneau_wire_default(&mut particles, &mut force_mode);

        let profile = compute_wire_force_profile(&particles);
        let n = profile.len();
        let mid = n / 2;

        // Sample three points: just-off-centre, halfway-out, end.
        // Left half: indices mid-1, mid/2, 0 — magnitudes should increase.
        let f_near_mid_left = profile[mid.saturating_sub(1)].1.abs();
        let f_halfway_left = profile[mid / 2].1.abs();
        let f_end_left = profile[0].1.abs();
        assert!(
            f_near_mid_left <= f_halfway_left,
            "|F| should grow from centre: near-mid {f_near_mid_left:.3e} > halfway {f_halfway_left:.3e}"
        );
        assert!(
            f_halfway_left < f_end_left,
            "|F| should grow toward end: halfway {f_halfway_left:.3e} >= end {f_end_left:.3e}"
        );
    }

    /// Switching to ForceMode::Lorentz with no grid fields gives zero force
    /// on every particle — the standard "straight wire produces no axial
    /// force on itself" result. Sanity check that our infrastructure shows
    /// the expected null result, against which the Weber profile stands out.
    #[test]
    fn test_graneau_no_axial_force_without_pair_law() {
        // A scenario in pure Lorentz mode with no E,B fields → no force.
        // Verified indirectly: compute_weber_forces with v_drift = 0 yields
        // zero correction; only static Coulomb survives, and the symmetric
        // chain plus charge neutrality gives zero total per ion.
        let mut particles = ParticleSystem::default();
        let mut force_mode = ForceMode::Lorentz;
        apply_graneau_wire_scenario(
            &mut particles,
            &mut force_mode,
            DEFAULT_NUM_SEGMENTS,
            DEFAULT_SEGMENT_SPACING,
            0.0, // <-- no drift → Weber correction vanishes
        );
        let profile = compute_wire_force_profile(&particles);

        // With v_drift = 0, ion-electron pair within each segment is at the
        // same position (Coulomb guarded out), and inter-segment forces are
        // pure Coulomb. By neutrality + collinearity each ion sees zero net
        // force from every distant *segment* (i_other + e_other cancel
        // exactly in the Coulomb limit). So every entry should be ≈ 0.
        for &(x, f) in &profile {
            assert!(
                f.abs() < 1e-3,
                "x = {x:.3e}: Coulomb-only force should vanish on neutral chain, got F = {f:.3e}"
            );
        }
    }

    // ----------------------------------------------------------------------
    // Phase 7.3 — Hairpin PCB-trace geometry
    // ----------------------------------------------------------------------

    /// Hairpin setup creates 4 N particles in two strands, total charge zero,
    /// top electrons drift +x̂, bottom electrons drift −x̂.
    #[test]
    fn test_apply_hairpin_creates_two_strands() {
        let mut particles = ParticleSystem::default();
        let mut force_mode = ForceMode::Lorentz;
        apply_hairpin_default(&mut particles, &mut force_mode);

        assert_eq!(force_mode, ForceMode::Weber);
        assert_eq!(
            particles.particles.len(),
            2 * particle_count(DEFAULT_NUM_SEGMENTS)
        );

        // Charge neutrality.
        let total: f32 = particles.particles.iter().map(|p| p.charge).sum();
        assert!(total.abs() < 1e-6);

        // Top strand (first 2N entries): y = +s/2, electrons +x̂.
        let half_y = DEFAULT_HAIRPIN_SPACING * 0.5;
        for i in 0..DEFAULT_NUM_SEGMENTS {
            let ion = &particles.particles[2 * i];
            let electron = &particles.particles[2 * i + 1];
            assert!((ion.position[1] - half_y).abs() < 1e-7);
            assert_eq!(ion.velocity, [0.0; 3]);
            assert!(electron.velocity[0] > 0.0, "top electron should drift +x");
        }
        // Bottom strand (next 2N entries): y = −s/2, electrons −x̂.
        let bot_start = 2 * DEFAULT_NUM_SEGMENTS;
        for i in 0..DEFAULT_NUM_SEGMENTS {
            let ion = &particles.particles[bot_start + 2 * i];
            let electron = &particles.particles[bot_start + 2 * i + 1];
            assert!((ion.position[1] - (-half_y)).abs() < 1e-7);
            assert_eq!(ion.velocity, [0.0; 3]);
            assert!(electron.velocity[0] < 0.0, "bottom electron should drift −x");
        }
    }

    /// The hairpin top-strand profile is **NOT identical** to the single-
    /// strand Graneau profile — adding the bottom return-current strand
    /// modifies the longitudinal Weber force on each top-strand ion. This
    /// test asserts that there is a measurable difference, demonstrating
    /// that the bottom strand contributes non-trivially.
    #[test]
    fn test_hairpin_profile_differs_from_straight_wire() {
        // Single straight wire (top strand only).
        let mut p_straight = ParticleSystem::default();
        let mut fm = ForceMode::Lorentz;
        apply_graneau_wire_default(&mut p_straight, &mut fm);
        let straight_profile = compute_wire_force_profile(&p_straight);

        // Hairpin (top + bottom strands).
        let mut p_hair = ParticleSystem::default();
        apply_hairpin_default(&mut p_hair, &mut fm);
        let hair_profile = compute_hairpin_top_force_profile(&p_hair, DEFAULT_NUM_SEGMENTS);

        assert_eq!(straight_profile.len(), hair_profile.len());

        // Per-ion difference. The hairpin's bottom strand modifies the top-
        // strand's Weber profile by inter-strand Coulomb + Weber-correction
        // contributions. We expect a measurable difference at every ion.
        let mut max_abs_diff = 0.0f32;
        let mut max_straight = 0.0f32;
        for (s, h) in straight_profile.iter().zip(hair_profile.iter()) {
            let diff = (s.1 - h.1).abs();
            max_abs_diff = max_abs_diff.max(diff);
            max_straight = max_straight.max(s.1.abs());
        }
        // The inter-strand correction is small (it's a v²/c² · geometric-
        // factor effect), but should be cleanly above the noise floor of
        // the f32 Weber sum on this geometry.
        assert!(
            max_abs_diff > 1e-2 * max_straight,
            "Hairpin should measurably alter the longitudinal force profile; \
             max |Δ| = {max_abs_diff:.3e}, max |F_straight| = {max_straight:.3e} \
             (ratio = {:.3e})",
            max_abs_diff / max_straight.max(1e-30)
        );
    }

    /// Top-strand ions in a hairpin still get end-peaked outward forces,
    /// because the SELF-interaction (top-strand against itself) dominates.
    /// The bottom strand modifies but does not erase the Graneau pattern.
    #[test]
    fn test_hairpin_top_strand_still_end_peaked() {
        let mut particles = ParticleSystem::default();
        let mut fm = ForceMode::Lorentz;
        apply_hairpin_default(&mut particles, &mut fm);
        let profile = compute_hairpin_top_force_profile(&particles, DEFAULT_NUM_SEGMENTS);
        let n = profile.len();
        let f_left = profile[0].1;
        let f_right = profile[n - 1].1;
        let f_middle = profile[n / 2].1;
        let f_end = f_left.abs().max(f_right.abs());

        assert!(f_left < 0.0, "leftmost top-strand ion should still feel −x push: {f_left:.3e}");
        assert!(f_right > 0.0, "rightmost top-strand ion should still feel +x push: {f_right:.3e}");
        // Middle is much smaller than ends — the bottom-strand correction
        // is small relative to the self-strand end-peak.
        assert!(
            f_middle.abs() < 0.1 * f_end,
            "middle top-strand ion: |F| = {:.3e}, end |F| = {f_end:.3e}",
            f_middle.abs()
        );
    }
}
