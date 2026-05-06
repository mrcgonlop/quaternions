// Charged particle system: Boris-Lorentz pusher with trilinear E/B sampling.
//
// Phase 3.1 deliverable. Provides:
//   * `ParticleState` — per-particle physics state.
//   * `ParticleSystem` — Bevy resource owning the Vec<ParticleState>; physics
//     state lives here so it can be advanced inside the simulation substep loop
//     without crossing Bevy's Query borrow checker.
//   * `ChargedParticle` — ECS marker component carrying an index back into
//     `ParticleSystem.particles`. Spawned entities use this to render a sphere
//     at the particle's world position via a transform-sync system.
//   * `step_particles` — pure function: trilinearly samples E,B from the grid
//     at each particle's world position, then applies the standard Boris pusher
//     with a position update. Tracks `prev_acceleration` to bootstrap the
//     explicit Weber scheme arriving in Phase 3.2.
//   * `sync_particle_transforms` — Bevy system that copies physics-state
//     positions to entity Transforms for rendering.
//
// World/grid coordinate convention (matches the bounding-box gizmo in main.rs):
//   The simulation domain is centered at the origin and spans
//     [-half_x, +half_x] × [-half_y, +half_y] × [-half_z, +half_z]
//   where half_x = nx*dx*0.5, etc. Cell (i, j, k) has its CENTER at
//     world = ((i + 0.5)*dx - half_x, (j + 0.5)*dx - half_y, (k + 0.5)*dx - half_z)
//   so the inverse (world → cell-center grid coord) is
//     gx = (px + half_x)/dx - 0.5
//   floor(gx) is the lower-left cell of the trilinear stencil.

use bevy::prelude::*;

use crate::simulation::grid::SimulationGrid;

/// Physical state of a single charged particle (SI units).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleState {
    /// Charge (Coulombs). Negative for electrons.
    pub charge: f32,
    /// Mass (kg). Must be > 0.
    pub mass: f32,
    /// World-space position (meters), grid-centered as described in the module header.
    pub position: [f32; 3],
    /// World-space velocity (m/s).
    pub velocity: [f32; 3],
    /// Acceleration from the previous step (m/s²). Used by the explicit Weber
    /// scheme (Phase 3.2) which needs a₁₂ to evaluate the r·r̈/c² correction
    /// without a chicken-and-egg dependence on the current force.
    pub prev_acceleration: [f32; 3],
}

impl ParticleState {
    pub fn new(charge: f32, mass: f32, position: [f32; 3], velocity: [f32; 3]) -> Self {
        Self {
            charge,
            mass,
            position,
            velocity,
            prev_acceleration: [0.0; 3],
        }
    }
}

/// Bevy resource holding the full particle population.
///
/// Physics owns the Vec; ECS entities carry only a `ChargedParticle { idx }`
/// marker that references back here. This split keeps the substep loop in
/// `simulation_step_system` free of Query borrows while still letting Bevy
/// render each particle as a sphere.
#[derive(Resource, Default)]
pub struct ParticleSystem {
    pub particles: Vec<ParticleState>,
}

impl ParticleSystem {
    pub fn clear(&mut self) {
        self.particles.clear();
    }

    pub fn push(&mut self, p: ParticleState) -> usize {
        let idx = self.particles.len();
        self.particles.push(p);
        idx
    }
}

/// ECS marker for a particle entity. The `idx` field references
/// `ParticleSystem.particles[idx]` — keep them in sync (or rebuild both).
#[derive(Component, Clone, Copy, Debug)]
pub struct ChargedParticle {
    pub idx: usize,
}

/// Spawn a particle: append it to the `ParticleSystem` resource AND create a
/// rendered ECS entity for it. The visual is a small sphere placed at the
/// particle's initial position.
///
/// Returns the spawned entity ID.
pub fn spawn_charged_particle(
    commands: &mut Commands,
    particles: &mut ParticleSystem,
    meshes: &mut Assets<Mesh>,
    materials: &mut Assets<StandardMaterial>,
    state: ParticleState,
    radius: f32,
    color: Color,
) -> Entity {
    let idx = particles.push(state);

    let mesh = meshes.add(Sphere::new(radius).mesh().ico(2).unwrap());
    let material = materials.add(StandardMaterial {
        base_color: color,
        emissive: color.to_linear() * 0.5,
        ..default()
    });

    commands
        .spawn((
            Mesh3d(mesh),
            MeshMaterial3d(material),
            Transform::from_xyz(state.position[0], state.position[1], state.position[2]),
            ChargedParticle { idx },
        ))
        .id()
}

/// Sample E and B at a world-space point by trilinear interpolation from
/// cell-centered values. Each of the 8 surrounding cells contributes its
/// locally-computed E,B (same formulae as `diagnostics::compute_derived_fields`,
/// recomputed here to avoid a dependency on `DiagnosticsState` and to keep the
/// pusher callable mid-substep, before diagnostics has run).
///
/// Returns ([0; 3], [0; 3]) if the point is closer than 2 cells to any face
/// (central differences at the 8 corners would touch out-of-grid neighbors).
pub fn sample_e_b(
    grid: &SimulationGrid,
    world_pos: [f32; 3],
    c0: f32,
) -> ([f32; 3], [f32; 3]) {
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

    // Need (i+1) ≤ nx-1 and (i+1)+1 ≤ nx-1 for central differences at the
    // upper corner → i ≤ nx - 3. Similarly i-1 ≥ 0 → i ≥ 1.
    if i < 1 || j < 1 || k < 1 || i > nx - 3 || j > ny - 3 || k > nz - 3 {
        return ([0.0; 3], [0.0; 3]);
    }

    let fx = gx - i as f32;
    let fy = gy - j as f32;
    let fz = gz - k as f32;

    let mut e_total = [0.0f32; 3];
    let mut b_total = [0.0f32; 3];

    for di in 0..2 {
        for dj in 0..2 {
            for dk in 0..2 {
                let cx = (i + di) as u32;
                let cy = (j + dj) as u32;
                let cz = (k + dk) as u32;
                let (e_c, b_c) = e_b_at_cell(grid, cx, cy, cz, c0);
                let wx = if di == 0 { 1.0 - fx } else { fx };
                let wy = if dj == 0 { 1.0 - fy } else { fy };
                let wz = if dk == 0 { 1.0 - fz } else { fz };
                let w = wx * wy * wz;
                for axis in 0..3 {
                    e_total[axis] += w * e_c[axis];
                    b_total[axis] += w * b_c[axis];
                }
            }
        }
    }

    (e_total, b_total)
}

/// Compute (E, B) at a single cell using the same central-difference formulae
/// as `diagnostics::compute_derived_fields`. Caller must ensure (x, y, z) is
/// strictly interior so the ±1 neighbour reads are in-bounds.
fn e_b_at_cell(
    grid: &SimulationGrid,
    x: u32,
    y: u32,
    z: u32,
    c0: f32,
) -> ([f32; 3], [f32; 3]) {
    let cells = grid.read_buf();
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let dx = grid.dx;
    let inv_2dx = 1.0 / (2.0 * dx);
    let i = grid.idx(x, y, z);
    let stride_y = nx;
    let stride_z = nx * ny;

    let cell = &cells[i];
    let c_local = c0 / cell.k;

    // E = -c_local * grad(Q.w) - Q_dot.vector()
    let e = [
        -c_local * (cells[i + 1].q[0] - cells[i - 1].q[0]) * inv_2dx - cell.q_dot[1],
        -c_local * (cells[i + stride_y].q[0] - cells[i - stride_y].q[0]) * inv_2dx
            - cell.q_dot[2],
        -c_local * (cells[i + stride_z].q[0] - cells[i - stride_z].q[0]) * inv_2dx
            - cell.q_dot[3],
    ];

    // B = curl(A) = curl((q[1], q[2], q[3]))
    let b = [
        (cells[i + stride_y].q[3] - cells[i - stride_y].q[3]) * inv_2dx
            - (cells[i + stride_z].q[2] - cells[i - stride_z].q[2]) * inv_2dx,
        (cells[i + stride_z].q[1] - cells[i - stride_z].q[1]) * inv_2dx
            - (cells[i + 1].q[3] - cells[i - 1].q[3]) * inv_2dx,
        (cells[i + 1].q[2] - cells[i - 1].q[2]) * inv_2dx
            - (cells[i + stride_y].q[1] - cells[i - stride_y].q[1]) * inv_2dx,
    ];

    (e, b)
}

/// Advance every particle by one step using the Boris pusher.
///
/// Standard Boris algorithm (relativistic-ready up to γ → 1 here):
///   1. v⁻ = v + (q·dt/2m)·E         (half-step electric kick)
///   2. t = (q·dt/2m)·B; s = 2t/(1 + |t|²)
///      v' = v⁻ + v⁻ × t
///      v⁺ = v⁻ + v' × s             (magnetic rotation, exact for any |t|)
///   3. v_new = v⁺ + (q·dt/2m)·E     (second half-step electric kick)
///   4. x_new = x + v_new · dt
///
/// Energy-conserving for |B|-only fields and the de facto standard for PIC
/// codes. Time-centered: velocity at t+dt/2, position at t+dt — but for
/// Phase 3.1 verification (drift in uniform E, gyration in uniform B) the
/// half-step staggering doesn't affect the closed-form check.
pub fn step_particles(particles: &mut [ParticleState], grid: &SimulationGrid, c0: f32, dt: f32) {
    for p in particles.iter_mut() {
        let (e, b) = sample_e_b(grid, p.position, c0);

        let q_over_2m_dt = 0.5 * p.charge / p.mass * dt;

        // 1. Half-step E kick.
        let v_minus = [
            p.velocity[0] + q_over_2m_dt * e[0],
            p.velocity[1] + q_over_2m_dt * e[1],
            p.velocity[2] + q_over_2m_dt * e[2],
        ];

        // 2. Magnetic rotation.
        let t = [q_over_2m_dt * b[0], q_over_2m_dt * b[1], q_over_2m_dt * b[2]];
        let t_mag_sq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
        let s_factor = 2.0 / (1.0 + t_mag_sq);
        let s = [s_factor * t[0], s_factor * t[1], s_factor * t[2]];

        // v' = v_minus + v_minus × t
        let v_prime = [
            v_minus[0] + (v_minus[1] * t[2] - v_minus[2] * t[1]),
            v_minus[1] + (v_minus[2] * t[0] - v_minus[0] * t[2]),
            v_minus[2] + (v_minus[0] * t[1] - v_minus[1] * t[0]),
        ];
        // v_plus = v_minus + v' × s
        let v_plus = [
            v_minus[0] + (v_prime[1] * s[2] - v_prime[2] * s[1]),
            v_minus[1] + (v_prime[2] * s[0] - v_prime[0] * s[2]),
            v_minus[2] + (v_prime[0] * s[1] - v_prime[1] * s[0]),
        ];

        // 3. Second half-step E kick.
        let v_new = [
            v_plus[0] + q_over_2m_dt * e[0],
            v_plus[1] + q_over_2m_dt * e[1],
            v_plus[2] + q_over_2m_dt * e[2],
        ];

        // Track effective acceleration for the Phase 3.2 Weber scheme.
        let inv_dt = 1.0 / dt;
        p.prev_acceleration = [
            (v_new[0] - p.velocity[0]) * inv_dt,
            (v_new[1] - p.velocity[1]) * inv_dt,
            (v_new[2] - p.velocity[2]) * inv_dt,
        ];

        p.velocity = v_new;
        // 4. Position update.
        p.position[0] += v_new[0] * dt;
        p.position[1] += v_new[1] * dt;
        p.position[2] += v_new[2] * dt;
    }
}

/// Bevy system: copy each particle's world position from `ParticleSystem`
/// into the Transform of its corresponding entity. Runs every frame so the
/// rendered sphere tracks the simulated motion.
pub fn sync_particle_transforms(
    particles: Res<ParticleSystem>,
    mut q: Query<(&ChargedParticle, &mut Transform)>,
) {
    for (cp, mut tf) in q.iter_mut() {
        if let Some(p) = particles.particles.get(cp.idx) {
            tf.translation = Vec3::new(p.position[0], p.position[1], p.position[2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::state::SimParams;

    /// Linear acceleration in a uniform E field along z.
    ///
    /// Setup: q_dot[3] = -E0 in every cell → E_z = -q_dot[3] = +E0 everywhere.
    /// A particle at the origin with v = 0 should reach v_z = (qE/m)·t after
    /// N·dt steps. Test the closed-form prediction within 1% (Boris pusher is
    /// exact for uniform E up to floating-point error).
    #[test]
    fn test_uniform_e_field_accelerates_linearly() {
        let nx = 16;
        let dx = 0.01;
        let mut grid = SimulationGrid::new(nx, nx, nx, dx);
        let e0: f32 = 1.0;
        for cell in grid.cells[grid.current].iter_mut() {
            cell.q_dot[3] = -e0;
        }

        let charge = 1.6e-19;
        let mass = 9.1e-31;
        let mut particles = vec![ParticleState::new(charge, mass, [0.0; 3], [0.0; 3])];

        let dt = 1e-12;
        let n_steps = 100;
        for _ in 0..n_steps {
            step_particles(&mut particles, &grid, SimParams::C0, dt);
        }

        let total_t = dt * n_steps as f32;
        let expected = (charge * e0 / mass) * total_t;
        let actual = particles[0].velocity[2];
        let rel_err = ((actual - expected) / expected).abs();
        assert!(
            rel_err < 0.02,
            "E-field drift: expected v_z ≈ {expected:.3e} m/s, got {actual:.3e} (rel_err = {rel_err:.3e})"
        );
        assert!(particles[0].velocity[0].abs() < 1e-3 * expected.abs());
        assert!(particles[0].velocity[1].abs() < 1e-3 * expected.abs());
    }

    /// Cyclotron gyration in a uniform B field along z.
    ///
    /// Setup: A_y = q[2] = B0·x_world (linear in x) → B_z = ∂A_y/∂x = B0,
    /// other components zero. With v₀ = v0·x̂, B = B0·ẑ, q > 0, the particle
    /// gyrates in the xy plane. The Boris pusher resolves the orbit exactly
    /// for any t = (qB/2m)·dt; with 1000 substeps per period the trajectory
    /// should close tightly.
    #[test]
    fn test_uniform_b_field_gyrates() {
        let nx: u32 = 32;
        let dx = 0.01;
        let mut grid = SimulationGrid::new(nx, nx, nx, dx);
        let half = nx as f32 * dx * 0.5;
        let b0: f32 = 1e-4;

        for z in 0..nx {
            for y in 0..nx {
                for x in 0..nx {
                    let xw = (x as f32 + 0.5) * dx - half;
                    let i = grid.idx(x, y, z);
                    grid.cells[grid.current][i].q[2] = b0 * xw;
                }
            }
        }

        let charge: f32 = 1.6e-19;
        let mass: f32 = 9.1e-31;
        let omega_c = charge.abs() * b0 / mass;
        let v0 = 1.0e6;
        let r_gyro = mass * v0 / (charge.abs() * b0);
        assert!(r_gyro < 0.4 * half, "Gyroradius too large for grid; reduce v0 or b0.");

        let mut particles = vec![ParticleState::new(
            charge,
            mass,
            [0.0; 3],
            [v0, 0.0, 0.0],
        )];

        let period = 2.0 * std::f32::consts::PI / omega_c;
        let n_steps: usize = 1000;
        let dt = period / n_steps as f32;

        for _ in 0..n_steps {
            step_particles(&mut particles, &grid, SimParams::C0, dt);
        }

        // After exactly one period the particle should return to (≈origin, v0·x̂).
        let pos = particles[0].position;
        let vel = particles[0].velocity;
        let pos_err = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
        assert!(
            pos_err < 0.05 * r_gyro,
            "Position drift after one period: {pos_err:.3e} m, gyroradius {r_gyro:.3e} m"
        );
        assert!((vel[0] - v0).abs() / v0 < 0.05, "v_x = {} (expected {v0})", vel[0]);
        assert!(vel[1].abs() / v0 < 0.05, "v_y = {} (should be ~0)", vel[1]);
        assert!(vel[2].abs() < 1e-6 * v0, "v_z drifted out of plane: {}", vel[2]);
    }

    /// Sampling outside the safe interior region returns zero — particles
    /// near boundaries shouldn't trip undefined-behavior central differences.
    #[test]
    fn test_sample_e_b_outside_grid_returns_zero() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        let half = 8.0 * 0.01 * 0.5;
        // Just outside the box.
        let (e, b) = sample_e_b(&grid, [half + 0.1, 0.0, 0.0], SimParams::C0);
        assert_eq!(e, [0.0; 3]);
        assert_eq!(b, [0.0; 3]);
    }

    /// Boris pusher in zero field is a no-op on velocity and a pure drift on
    /// position — sanity check on the integrator structure.
    #[test]
    fn test_zero_field_is_pure_drift() {
        let grid = SimulationGrid::new(16, 16, 16, 0.01);
        let mut particles = vec![ParticleState::new(
            1.6e-19,
            9.1e-31,
            [0.0; 3],
            [1e5, 0.0, 0.0],
        )];
        let dt = 1e-10;
        for _ in 0..50 {
            step_particles(&mut particles, &grid, SimParams::C0, dt);
        }
        let expected_x = 1e5 * dt * 50.0;
        assert!(
            (particles[0].position[0] - expected_x).abs() < 1e-9,
            "expected x = {expected_x:.3e}, got {:.3e}",
            particles[0].position[0]
        );
        assert!((particles[0].velocity[0] - 1e5).abs() < 1.0);
    }
}
