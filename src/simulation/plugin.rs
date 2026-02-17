// Bevy plugin: registers simulation resources and systems.

use bevy::prelude::*;

use crate::simulation::boundaries::{self, BoundaryConfig, PmlState};
use crate::simulation::diagnostics::{self, DiagnosticsState};
use crate::simulation::field_update::step_field_cpu;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::sources::{self, SourceConfig};
use crate::simulation::state::PmlConfig;

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
            .add_systems(Startup, init_grid)
            .add_systems(
                Update,
                (
                    sources::source_injection_system,
                    simulation_step_system,
                    boundaries::boundary_system,
                    diagnostics::diagnostics_system,
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
fn simulation_step_system(
    grid: Option<ResMut<SimulationGrid>>,
    config: Res<SimulationConfig>,
    pml: Option<ResMut<crate::simulation::boundaries::PmlState>>,
) {
    let Some(mut grid) = grid else { return };
    if config.paused {
        return;
    }

    let effective_dt = grid.dt * config.dt_factor.clamp(0.01, 1.0);
    let params = grid.sim_params_with_dt(config.extended_mode, effective_dt);

    // Convert Option<ResMut<PmlState>> to Option<&mut PmlState> for the field update
    let mut pml = pml;
    for _ in 0..config.steps_per_frame {
        step_field_cpu(&mut grid, &params, pml.as_deref_mut());
        grid.swap_and_advance_with_dt(effective_dt);
    }
}
