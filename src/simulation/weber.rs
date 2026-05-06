// Weber electrodynamics: direct particle-particle force law extending Coulomb.
//
// Force on particle 1 due to particle 2:
//     F₁ = (q₁q₂ / 4πε₀ r²) r̂ [1 − ṙ²/(2c²) + r·r̈/c²]
// where:
//     r₁₂ = r₁ − r₂   (displacement of 1 from 2; r̂ points from 2 toward 1)
//     r   = |r₁₂|
//     v₁₂ = v₁ − v₂
//     a₁₂ = a₁ − a₂   (uses prev-step acceleration in this explicit scheme)
//     ṙ   = r̂ · v₁₂                         (relative radial velocity)
//     r̈   = r̂ · a₁₂ + (|v₁₂|² − ṙ²) / r    (relative radial acceleration)
//
// Three terms inside the brackets:
//   1.  [1]                    standard Coulomb interaction.
//   2.  [−ṙ²/(2c²)]            velocity correction — always reduces |F| for
//                              charges in relative motion.
//   3.  [+r·r̈/c²]              acceleration correction — produces the
//                              longitudinal forces along current direction
//                              that Lorentz / Biot-Savart cannot account for.
//
// For static charges (v = a = 0) the correction vanishes and the law reduces
// to Coulomb. For drifting charges in a wire the velocity correction is tiny
// (v_drift²/c² ≈ 10⁻¹⁸ for copper at 1 mA), but the collective effect over
// many pairs in a chain produces the longitudinal force measured in Graneau-
// type wire-fragmentation experiments (Phase 3.3).
//
// Explicit acceleration scheme:
//   a₁₂ involves the present-step acceleration, which depends on the very
//   force we're computing. To break the implicit dependence we read the
//   PREVIOUS step's acceleration from each particle's `prev_acceleration`
//   field. The Lorentz Boris pusher and `step_particles_weber` both keep this
//   field updated with each particle's actual per-step acceleration. For
//   v ≪ c the explicit scheme is first-order accurate and stable; if higher
//   accuracy becomes important, switch to a predictor-corrector. See
//   TODO.md §Known Issues.

use bevy::prelude::*;

use crate::math::vector_field;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::particles::{self, ParticleState};

/// Vacuum permittivity (F/m). Matches the value used in `diagnostics.rs`.
const EPSILON_0: f32 = 8.854e-12;

/// Selects which force law drives the particle integrator each substep.
///
/// Stored as a Bevy resource so the UI can toggle it at runtime; the
/// `simulation_step_system` dispatches on the value before calling the
/// matching pusher.
#[derive(Resource, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ForceMode {
    /// Standard EM: Lorentz force F = q(E + v × B), Boris pusher, fields
    /// sampled trilinearly from the grid. (Phase 3.1.)
    #[default]
    Lorentz,
    /// Weber pair force only: no E/B coupling, no field sampling. Pure
    /// action-at-a-distance between charged particles.
    Weber,
    /// Apply BOTH forces additively (Boris pusher + Weber kick). Useful for
    /// validating that Weber reduces to Lorentz in the appropriate limit, or
    /// for scenarios where both grid-coupled radiation and Weber longitudinal
    /// forces matter.
    Both,
}

impl ForceMode {
    pub const ALL: &'static [ForceMode] =
        &[ForceMode::Lorentz, ForceMode::Weber, ForceMode::Both];

    pub fn name(self) -> &'static str {
        match self {
            ForceMode::Lorentz => "Lorentz",
            ForceMode::Weber => "Weber",
            ForceMode::Both => "Both",
        }
    }
}

/// Weber force on particle 1 due to particle 2.
///
/// Inputs are the relative quantities r₁₂ = r₁ − r₂, v₁₂ = v₁ − v₂, a₁₂ =
/// a₁ − a₂. Returns the force vector along ±r̂.
///
/// Returns zero if particles are coincident (|r₁₂| < ε) — caller is
/// responsible for not putting two particles at the same position.
pub fn weber_force(
    q1: f32,
    q2: f32,
    r12: [f32; 3],
    v12: [f32; 3],
    a12: [f32; 3],
    c0: f32,
) -> [f32; 3] {
    let r_sq = vector_field::magnitude_sq(r12);
    if r_sq < 1e-30 {
        return [0.0; 3];
    }
    let r = r_sq.sqrt();
    let inv_r = 1.0 / r;
    let r_hat = vector_field::scale(r12, inv_r);

    let r_dot = vector_field::dot(r_hat, v12);
    let v_sq = vector_field::magnitude_sq(v12);
    let r_ddot = vector_field::dot(r_hat, a12) + (v_sq - r_dot * r_dot) * inv_r;

    let c_sq = c0 * c0;
    let coulomb_mag = q1 * q2 / (4.0 * std::f32::consts::PI * EPSILON_0 * r_sq);
    let weber_factor = 1.0 - r_dot * r_dot / (2.0 * c_sq) + r * r_ddot / c_sq;

    vector_field::scale(r_hat, coulomb_mag * weber_factor)
}

/// Total Weber force on every particle from every other particle.
///
/// O(N²) in pair count but uses Newton's third law (F₁₂ = −F₂₁) to compute
/// each pair only once: the inner loop runs over j > i and updates both
/// `forces[i]` and `forces[j]`.
pub fn compute_weber_forces(particles: &[ParticleState], c0: f32) -> Vec<[f32; 3]> {
    let n = particles.len();
    let mut forces = vec![[0.0f32; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let r_ij =
                vector_field::sub(particles[i].position, particles[j].position);
            let v_ij =
                vector_field::sub(particles[i].velocity, particles[j].velocity);
            let a_ij = vector_field::sub(
                particles[i].prev_acceleration,
                particles[j].prev_acceleration,
            );
            let f = weber_force(
                particles[i].charge,
                particles[j].charge,
                r_ij,
                v_ij,
                a_ij,
                c0,
            );
            forces[i] = vector_field::add(forces[i], f);
            forces[j] = vector_field::sub(forces[j], f);
        }
    }
    forces
}

/// Advance every particle by one step using ONLY the Weber force.
///
/// Kick-drift (semi-implicit Euler) integrator:
///     v_new = v + (F/m)·dt
///     x_new = x + v_new·dt
/// `prev_acceleration` is set to F/m so the next step's Weber force has it
/// available for the r̈ term.
pub fn step_particles_weber(
    particles: &mut [ParticleState],
    c0: f32,
    dt: f32,
) {
    let forces = compute_weber_forces(particles, c0);
    for (p, f) in particles.iter_mut().zip(forces.iter()) {
        let inv_m = 1.0 / p.mass;
        let a = [f[0] * inv_m, f[1] * inv_m, f[2] * inv_m];
        p.velocity[0] += a[0] * dt;
        p.velocity[1] += a[1] * dt;
        p.velocity[2] += a[2] * dt;
        p.position[0] += p.velocity[0] * dt;
        p.position[1] += p.velocity[1] * dt;
        p.position[2] += p.velocity[2] * dt;
        p.prev_acceleration = a;
    }
}

/// Advance every particle by one step applying BOTH Lorentz and Weber forces.
///
/// Procedure:
///   1. Compute Weber pair forces using positions at the start of the step.
///   2. Run the Lorentz Boris pusher: this updates velocity AND position and
///      sets `prev_acceleration` to the Lorentz acceleration over the step.
///   3. Apply the Weber kick: `v += (F_weber/m)·dt`. We do NOT add another
///      drift term here — the Boris pusher already advanced position by the
///      full dt. For small Weber corrections (v ≪ c) the position-update
///      lag is bounded by the relative magnitude of Weber/Lorentz and is
///      negligible.
///   4. Add Weber's contribution into `prev_acceleration` so next step's
///      Weber force computation uses the total acceleration from this step.
pub fn step_particles_both(
    particles: &mut [ParticleState],
    grid: &SimulationGrid,
    c0: f32,
    dt: f32,
) {
    let weber_forces = compute_weber_forces(particles, c0);
    particles::step_particles(particles, grid, c0, dt);
    for (p, f) in particles.iter_mut().zip(weber_forces.iter()) {
        let inv_m = 1.0 / p.mass;
        p.velocity[0] += f[0] * inv_m * dt;
        p.velocity[1] += f[1] * inv_m * dt;
        p.velocity[2] += f[2] * inv_m * dt;
        p.prev_acceleration[0] += f[0] * inv_m;
        p.prev_acceleration[1] += f[1] * inv_m;
        p.prev_acceleration[2] += f[2] * inv_m;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::state::SimParams;

    const C0: f32 = SimParams::C0;

    /// Static charges (v = a = 0): Weber force collapses to Coulomb's law
    /// exactly — the correction factor [1 − ṙ²/(2c²) + r·r̈/c²] = 1.
    #[test]
    fn test_weber_static_limit_equals_coulomb() {
        let q1: f32 = 1.0;
        let q2: f32 = 1.0;
        let r12 = [1.0_f32, 0.0, 0.0];
        let r: f32 = 1.0;

        let f = weber_force(q1, q2, r12, [0.0; 3], [0.0; 3], C0);

        let coulomb = q1 * q2 / (4.0 * std::f32::consts::PI * EPSILON_0 * r * r);
        let expected = [coulomb, 0.0, 0.0];

        for axis in 0..3 {
            let rel_err = (f[axis] - expected[axis]).abs() / coulomb.abs();
            assert!(
                rel_err < 1e-5,
                "axis {axis}: expected {}, got {}",
                expected[axis],
                f[axis]
            );
        }
    }

    /// Two like charges moving directly apart: |F_Weber| < |F_Coulomb|.
    /// Specifically, with v₁₂ = (2v, 0, 0) along +x̂, the correction factor
    /// becomes 1 − 2v²/c². (a₁₂ = 0 → r̈ from the (|v|² − ṙ²)/r term
    /// vanishes since v₁₂ ∥ r̂.)
    ///
    /// Uses unit charges to keep magnitudes in the f32 normal range —
    /// elementary-charge scales (q ≈ 1.6e-19) push squared force magnitudes
    /// into f32 subnormals where the 2v²/c² ≈ 2e-5 relative correction is
    /// numerically lost.
    #[test]
    fn test_weber_moving_apart_less_than_coulomb() {
        let q1: f32 = 1.0;
        let q2: f32 = 1.0;
        let r12 = [1.0_f32, 0.0, 0.0];
        let v: f32 = 1.0e6;
        let v12 = [2.0 * v, 0.0, 0.0];

        let f_weber = weber_force(q1, q2, r12, v12, [0.0; 3], C0);
        let f_coulomb = weber_force(q1, q2, r12, [0.0; 3], [0.0; 3], C0);

        // Both forces are along +x̂; compare components directly to avoid
        // sqrt(magnitude²) underflow even at unit-charge scales.
        assert!(
            f_weber[0] < f_coulomb[0],
            "F_Weber.x = {}, F_Coulomb.x = {} — Weber should be smaller when charges separate",
            f_weber[0],
            f_coulomb[0]
        );

        let expected_factor = 1.0 - 2.0 * v * v / (C0 * C0);
        let actual_factor = f_weber[0] / f_coulomb[0];
        let rel_err = ((actual_factor - expected_factor) / expected_factor).abs();
        assert!(
            rel_err < 1e-3,
            "Velocity correction factor: expected {expected_factor}, got {actual_factor} (rel_err = {rel_err:.3e})"
        );
    }

    /// Force is collinear with r̂ regardless of v₁₂, a₁₂ orientation.
    /// Cross product F × r₁₂ should vanish to numerical precision.
    ///
    /// Unit charges, meter-scale separation — keeps |F| in normal f32 range
    /// so the cross-product / magnitude check has stable precision.
    #[test]
    fn test_weber_force_along_r_hat() {
        let r12 = [3.0_f32, 4.0, 0.0];
        let v12 = [1.0e5_f32, 2.0e5, 3.0e5];
        let a12 = [1.0e10_f32, -2.0e10, 1.0e10];
        let q1: f32 = 1.0;
        let q2: f32 = -1.0;

        let f = weber_force(q1, q2, r12, v12, a12, C0);

        let cross = vector_field::cross(f, r12);
        let f_mag = vector_field::magnitude(f);
        let r_mag = vector_field::magnitude(r12);
        let cross_mag = vector_field::magnitude(cross);
        // sin(angle between F and r̂). Tolerance is loose because the cross
        // product computes a small difference of large products in f32.
        let sin_theta = cross_mag / (f_mag * r_mag);
        assert!(
            sin_theta < 1e-5,
            "F not aligned with r̂: sin(θ) = {sin_theta:.3e}"
        );
    }

    /// Newton's third law: F on 1 due to 2 = −F on 2 due to 1.
    /// Achieved automatically because the force formula is symmetric under
    /// the simultaneous swaps q₁↔q₂, r₁₂→−r₁₂, v₁₂→−v₁₂, a₁₂→−a₁₂.
    #[test]
    fn test_weber_newton_third_law() {
        let r12 = [3.0_f32, 4.0, -1.0];
        let v12 = [1.0e5_f32, -2.0e5, 3.0e5];
        let a12 = [1.0e10_f32, 2.0e10, -1.0e10];
        let q1: f32 = 1.0;
        let q2: f32 = 1.5;

        let f12 = weber_force(q1, q2, r12, v12, a12, C0);
        let r21 = [-r12[0], -r12[1], -r12[2]];
        let v21 = [-v12[0], -v12[1], -v12[2]];
        let a21 = [-a12[0], -a12[1], -a12[2]];
        let f21 = weber_force(q2, q1, r21, v21, a21, C0);

        for axis in 0..3 {
            let sum = f12[axis] + f21[axis];
            let scale = f12[axis].abs().max(f21[axis].abs()).max(1e-30);
            assert!(
                (sum / scale).abs() < 1e-5,
                "axis {axis}: F₁₂ = {}, F₂₁ = {}, sum = {sum}",
                f12[axis],
                f21[axis]
            );
        }
    }

    /// Two like charges at rest, no fields: under Weber-only integration
    /// they should accelerate apart along the line connecting them.
    #[test]
    fn test_step_particles_weber_two_charges_repel() {
        let q = 1.6e-19;
        let m = 9.1e-31;
        let r0: f32 = 1.0e-9;

        let mut particles = vec![
            ParticleState::new(q, m, [-r0 * 0.5, 0.0, 0.0], [0.0; 3]),
            ParticleState::new(q, m, [r0 * 0.5, 0.0, 0.0], [0.0; 3]),
        ];

        let dt = 1.0e-18;
        for _ in 0..100 {
            step_particles_weber(&mut particles, C0, dt);
        }

        // Particle 0 must move in -x; particle 1 in +x; both at nontrivial speed.
        assert!(
            particles[0].position[0] < -r0 * 0.5,
            "Particle 0 didn't move in -x: pos = {:?}",
            particles[0].position
        );
        assert!(
            particles[1].position[0] > r0 * 0.5,
            "Particle 1 didn't move in +x: pos = {:?}",
            particles[1].position
        );
        assert!(particles[0].velocity[0] < 0.0);
        assert!(particles[1].velocity[0] > 0.0);
        // Force is purely along x — no transverse motion.
        assert!(
            particles[0].velocity[1].abs()
                < 1.0e-3 * particles[0].velocity[0].abs(),
            "Unexpected y-velocity: {}",
            particles[0].velocity[1]
        );
        assert!(
            particles[0].velocity[2].abs()
                < 1.0e-3 * particles[0].velocity[0].abs(),
            "Unexpected z-velocity: {}",
            particles[0].velocity[2]
        );
        // Newton's third law: total momentum should remain zero.
        let p_total_x = particles[0].velocity[0] + particles[1].velocity[0];
        let p_scale = particles[0].velocity[0].abs().max(1.0);
        assert!(
            (p_total_x / p_scale).abs() < 1.0e-4,
            "Total momentum should vanish: got {p_total_x}"
        );
    }
}
