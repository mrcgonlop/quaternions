// Bevy plugin: registers simulation resources and systems.

use bevy::prelude::*;

use crate::simulation::boundaries::{self, BoundaryConfig, PmlState};
use crate::simulation::diagnostics::{self, DiagnosticsState, ProbeSet};
use crate::simulation::field_update::step_field_cpu;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::particles::{self, ParticleSystem};
use crate::simulation::sources::{self, SourceConfig};
use crate::simulation::state::PmlConfig;
use crate::simulation::weber::{self, ForceMode};

/// Polarizable vacuum (K field) configuration.
///
/// Controls the dynamic K field that models QED vacuum polarization.
/// When enabled, the vacuum refractive index K(x,t) evolves as a damped
/// scalar wave driven by the local electromagnetic energy density:
///   ∂²K/∂t² = c_local² ∇²K − ωₚ²(K − 1) + η · u_field / u_s
/// where ωₚ is the vacuum plasma frequency, η is the QED coupling,
/// and u_s is the characteristic energy normalization.
///
/// The K-cycle resonator scenario (Phase 7.6) adds an external sinusoidal
/// drive term to the K-update inside a configurable spherical region —
/// when `k_drive_amplitude > 0` and a cell sits within `k_drive_radius` of
/// the grid centre, `k_drive_amplitude · sin(2π · k_drive_frequency · t)`
/// is added to k_dot each step. Default zero amplitude leaves the existing
/// K dynamics untouched.
#[derive(Resource, Clone, Debug)]
pub struct VacuumConfig {
    /// Master switch: set false to skip K update entirely (K stays constant).
    pub enabled: bool,
    /// Vacuum plasma frequency (rad/s in simulation units).
    /// ωₚ = mₑc²/ℏ ≈ 7.8×10²⁰ rad/s in SI; use ~c/100 for visible simulation dynamics.
    pub omega_p: f32,
    /// QED coupling constant: η ≈ 2α/(45π) ≈ 1e-4 (Euler-Heisenberg).
    pub eta: f32,
    /// Characteristic energy density normalization u_s = ε₀(mₑc²/e)².
    /// Prevents numerical overflow; default 1.0 means u_field is in code units.
    pub u_s: f32,
    /// Phase 7.6: amplitude of an external sinusoidal forcing on the K
    /// equation, applied inside the spherical drive region. Default 0.0
    /// (no drive). When > 0 the K-cycle resonator scenario uses this to
    /// model a rapidly switched EM source modulating the local refractive
    /// index. Drive frequency is `k_drive_frequency`.
    pub k_drive_amplitude: f32,
    /// Phase 7.6: drive frequency in code-units rad/(2π·s) (i.e., Hz). The
    /// drive enters the K leapfrog as
    ///   k_dot += k_drive_amplitude · sin(2π · k_drive_frequency · t)
    /// for cells inside the drive region.
    pub k_drive_frequency: f32,
    /// Phase 7.6: spherical drive region radius in world meters, centred
    /// on the grid origin. When 0, no drive is applied even if the
    /// amplitude is nonzero.
    pub k_drive_radius: f32,
}

impl Default for VacuumConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            omega_p: 0.0,
            eta: 1e-4,
            u_s: 1.0,
            k_drive_amplitude: 0.0,
            k_drive_frequency: 0.0,
            k_drive_radius: 0.0,
        }
    }
}

/// Configuration resource for the simulation.
#[derive(Resource)]
pub struct SimulationConfig {
    /// Grid dimensions (cells per axis).
    pub grid_size: u32,

    /// Cell spacing in meters.
    pub dx: f32,

    /// Whether the simulation is paused.
    pub paused: bool,

    /// Number of simulation steps to run per frame.
    pub steps_per_frame: u32,

    /// Toggle between standard EM and QVED extended mode.
    pub extended_mode: bool,

    /// Timestep scale factor (0.01–1.0). Multiplied with the CFL-derived dt.
    /// Lower values give finer temporal resolution at the cost of slower progression.
    /// 1.0 = full CFL timestep (maximum stable), 0.1 = 10× finer resolution.
    pub dt_factor: f32,

    /// When true, run exactly one step this frame even if paused, then clear the flag.
    /// Set by the Step button in the UI.
    pub step_requested: bool,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            grid_size: 32,
            dx: 0.01, // 1 cm cells
            paused: true,
            steps_per_frame: 1,
            extended_mode: false,
            dt_factor: 1.0,
            step_requested: false,
        }
    }
}

/// System set for all simulation systems. Visualization should run after this.
#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct SimulationSet;

/// Plugin that initializes the simulation grid and registers simulation systems.
pub struct SimulationPlugin;

impl Plugin for SimulationPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SimulationConfig>()
            .init_resource::<SourceConfig>()
            .init_resource::<DiagnosticsState>()
            .init_resource::<BoundaryConfig>()
            .init_resource::<VacuumConfig>()
            .init_resource::<ProbeSet>()
            .init_resource::<ParticleSystem>()
            .init_resource::<ForceMode>()
            .add_systems(Startup, init_grid)
            .add_systems(
                Update,
                (
                    // Source injection is now inside simulation_step_system's substep loop
                    // so that sources are sampled at the correct time each substep.
                    simulation_step_system,
                    boundaries::boundary_system,
                    diagnostics::diagnostics_system,
                    diagnostics::probe_system,
                    particles::sync_particle_transforms,
                )
                    .chain()
                    .in_set(SimulationSet),
            );
    }
}

/// Startup system: create the SimulationGrid and PmlState resources.
fn init_grid(
    mut commands: Commands,
    config: Res<SimulationConfig>,
    bc_config: Res<BoundaryConfig>,
) {
    let mut grid = SimulationGrid::new(
        config.grid_size,
        config.grid_size,
        config.grid_size,
        config.dx,
    );
    info!(
        "Simulation grid initialized: {}³ cells, dx={:.4} m, dt={:.4e} s",
        config.grid_size, config.dx, grid.dt
    );

    // Initialize CPML state for Open boundary faces
    let pml_config = PmlConfig::default();
    let pml_state = PmlState::new(&mut grid, &bc_config, pml_config);
    commands.insert_resource(pml_state);

    commands.insert_resource(grid);
}

/// Per-frame simulation step system.
/// Runs `steps_per_frame` iterations of the CPU field update when not paused.
/// Uses `dt_factor` to scale the timestep below the CFL maximum.
///
/// Source injection is performed inside the substep loop so that sources
/// are sampled at the correct simulation time for each substep, rather than
/// only once per frame (which would leave most substeps undriven at high
/// steps_per_frame values).
fn simulation_step_system(
    grid: Option<ResMut<SimulationGrid>>,
    mut config: ResMut<SimulationConfig>,
    source_config: Res<SourceConfig>,
    vacuum_config: Res<VacuumConfig>,
    pml: Option<ResMut<crate::simulation::boundaries::PmlState>>,
    mut particles: ResMut<ParticleSystem>,
    force_mode: Res<ForceMode>,
) {
    let Some(mut grid) = grid else { return };

    let single_step = config.step_requested;
    if config.paused && !single_step {
        return;
    }

    let effective_dt = grid.dt * config.dt_factor.clamp(0.01, 1.0);
    let step_count = if single_step { 1 } else { config.steps_per_frame };

    let mut pml = pml;
    for _ in 0..step_count {
        let params = grid.sim_params_with_dt(config.extended_mode, effective_dt);
        // Push particles using the field at time t BEFORE advancing it. After
        // step_field_cpu + swap, the particles and the field will both be at
        // t+dt, keeping their integration step-aligned.
        if !particles.particles.is_empty() {
            match *force_mode {
                ForceMode::Lorentz => particles::step_particles(
                    &mut particles.particles,
                    &grid,
                    params.c0,
                    effective_dt,
                ),
                ForceMode::Weber => weber::step_particles_weber(
                    &mut particles.particles,
                    params.c0,
                    effective_dt,
                ),
                ForceMode::Both => weber::step_particles_both(
                    &mut particles.particles,
                    &grid,
                    params.c0,
                    effective_dt,
                ),
            }
        }
        sources::inject_sources(&mut grid, &source_config, &params);
        step_field_cpu(&mut grid, &params, pml.as_deref_mut(), Some(&*vacuum_config));
        grid.swap_and_advance_with_dt(effective_dt);
    }

    if single_step {
        config.step_requested = false;
    }
}
