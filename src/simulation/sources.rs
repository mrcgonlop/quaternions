// Source injection: currents, charges, waveforms.
// Sources drive the simulation by injecting current density J into q_dot
// or setting phi at source cell locations.
//
// Runs BEFORE the field update system each frame.

use bevy::prelude::*;

use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::SimParams;

/// Type of electromagnetic source.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceType {
    /// Static point charge: sets phi at source cell.
    PointCharge,
    /// Oscillating electric dipole: sinusoidal current along a given axis.
    OscillatingDipole {
        /// Dipole axis: 0=x, 1=y, 2=z
        axis: usize,
    },
    /// Pulsed current: Gaussian envelope current injection.
    CurrentPulse {
        /// Dipole axis: 0=x, 1=y, 2=z
        axis: usize,
        /// Width of the Gaussian envelope in seconds.
        sigma: f32,
        /// Center time of the Gaussian peak in seconds.
        t_center: f32,
    },
}

/// A single electromagnetic source.
#[derive(Clone, Debug)]
pub struct Source {
    /// Grid position (fractional coordinates allowed; nearest cell used).
    pub position: [f32; 3],

    /// Source type determines how injection works.
    pub source_type: SourceType,

    /// Amplitude (units depend on source type):
    ///   PointCharge: charge in Coulombs (sets phi = q / (4πε₀r) at source cell)
    ///   OscillatingDipole: peak current density J0 (A/m²)
    ///   CurrentPulse: peak current density J0 (A/m²)
    pub amplitude: f32,

    /// Oscillation frequency in Hz (used by OscillatingDipole).
    pub frequency: f32,

    /// Phase offset in radians.
    pub phase: f32,

    /// Whether this source is currently active.
    pub active: bool,
}

impl Source {
    /// Create a z-directed oscillating dipole at the given grid position.
    pub fn dipole_z(position: [f32; 3], amplitude: f32, frequency: f32) -> Self {
        Self {
            position,
            source_type: SourceType::OscillatingDipole { axis: 2 },
            amplitude,
            frequency,
            phase: 0.0,
            active: true,
        }
    }

    /// Create a static point charge at the given grid position.
    pub fn point_charge(position: [f32; 3], amplitude: f32) -> Self {
        Self {
            position,
            source_type: SourceType::PointCharge,
            amplitude,
            frequency: 0.0,
            phase: 0.0,
            active: true,
        }
    }
}

/// Bevy resource holding all active sources.
#[derive(Resource, Default)]
pub struct SourceConfig {
    pub sources: Vec<Source>,
}

/// Inject all active sources into the simulation grid.
///
/// For each source:
/// - OscillatingDipole: injects J = J0 * sin(2πft + phase) into q_dot along the dipole axis.
///   q_dot[axis+1] += (J0/ε₀) * sin(ωt + phase) * dt
///   The +1 offset maps axis 0,1,2 → q_dot components 1,2,3 (vector part of Q).
///
/// - PointCharge: sets q[0] = phi/c at the source cell.
///   phi = amplitude / (4πε₀ * dx) → q[0] = phi / c0
///   (Approximate: sets the potential at the source cell directly.)
///
/// - CurrentPulse: Gaussian-envelope current injection.
///   J(t) = J0 * exp(-(t - t_center)² / (2σ²)) along axis.
pub fn inject_sources(grid: &mut SimulationGrid, sources: &SourceConfig, params: &SimParams) {
    let time = grid.time as f32;
    let dt = params.dt;
    let c0 = params.c0;
    let dx = grid.dx;
    // ε₀ = 1 / (μ₀ * c²), μ₀ = 4π×10⁻⁷
    let epsilon_0: f32 = 8.854e-12;

    for source in &sources.sources {
        if !source.active {
            continue;
        }

        // Convert position to nearest grid cell
        let gx = source.position[0].round() as u32;
        let gy = source.position[1].round() as u32;
        let gz = source.position[2].round() as u32;

        // Bounds check
        if gx >= grid.nx || gy >= grid.ny || gz >= grid.nz {
            continue;
        }

        let cell = grid.cell_mut(gx, gy, gz);

        match &source.source_type {
            SourceType::OscillatingDipole { axis } => {
                // J(t) = J0 * sin(2πft + phase)
                // The source drives ∂A/∂t at the source location.
                // For the vector potential: q_dot[axis+1] += (J0/ε₀) * sin(ωt) * dt
                // The axis+1 maps x,y,z axes → q components 1,2,3.
                let omega = 2.0 * std::f32::consts::PI * source.frequency;
                let j = source.amplitude * (omega * time + source.phase).sin();

                // Inject as acceleration of the vector potential component
                cell.q_dot[axis + 1] += (j / epsilon_0) * dt;
            }
            SourceType::PointCharge => {
                // Static charge: set phi/c at the source cell.
                // phi ≈ q / (4πε₀ * dx) (monopole approximation at the cell)
                let phi = source.amplitude / (4.0 * std::f32::consts::PI * epsilon_0 * dx);
                cell.q[0] = phi / c0;
            }
            SourceType::CurrentPulse { axis, sigma, t_center } => {
                // Gaussian-envelope current pulse
                // J(t) = J0 * exp(-(t - t_center)² / (2σ²))
                let t_diff = time - t_center;
                let j = source.amplitude * (-t_diff * t_diff / (2.0 * sigma * sigma)).exp();

                cell.q_dot[axis + 1] += (j / epsilon_0) * dt;
            }
        }
    }
}

/// Bevy system that injects sources before the field update.
pub fn source_injection_system(
    mut grid: Option<ResMut<SimulationGrid>>,
    sources: Option<Res<SourceConfig>>,
    config: Res<crate::simulation::plugin::SimulationConfig>,
) {
    let Some(ref mut grid) = grid else { return };
    let Some(ref sources) = sources else { return };
    if config.paused || sources.sources.is_empty() {
        return;
    }

    let effective_dt = grid.dt * config.dt_factor.clamp(0.01, 1.0);
    let params = grid.sim_params_with_dt(config.extended_mode, effective_dt);
    inject_sources(grid, sources, &params);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::grid::SimulationGrid;

    #[test]
    fn test_dipole_source_injects() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);

        let sources = SourceConfig {
            sources: vec![Source::dipole_z([4.0, 4.0, 4.0], 1.0, 1e9)],
        };

        // Before injection, q_dot should be zero
        assert_eq!(grid.cell(4, 4, 4).q_dot[3], 0.0);

        inject_sources(&mut grid, &sources, &params);

        // After injection at t=0 with phase=0, sin(0)=0 so no change at t=0.
        // Let's advance time and try again.
        grid.time = 0.25e-9; // quarter period for 1 GHz
        inject_sources(&mut grid, &sources, &params);

        // Now q_dot[3] (Az component) should be nonzero
        assert_ne!(grid.cell(4, 4, 4).q_dot[3], 0.0, "dipole should inject into q_dot");
    }

    #[test]
    fn test_point_charge_sets_phi() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);

        let sources = SourceConfig {
            sources: vec![Source::point_charge([4.0, 4.0, 4.0], 1e-10)],
        };

        inject_sources(&mut grid, &sources, &params);

        // q[0] = phi/c should be nonzero
        assert_ne!(grid.cell(4, 4, 4).q[0], 0.0, "point charge should set phi");
    }

    #[test]
    fn test_out_of_bounds_source_ignored() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);

        let sources = SourceConfig {
            sources: vec![Source::dipole_z([100.0, 100.0, 100.0], 1.0, 1e9)],
        };

        // Should not panic
        inject_sources(&mut grid, &sources, &params);
    }

    #[test]
    fn test_inactive_source_skipped() {
        let mut grid = SimulationGrid::new(8, 8, 8, 0.01);
        let params = grid.sim_params(false);

        let mut source = Source::dipole_z([4.0, 4.0, 4.0], 1.0, 1e9);
        source.active = false;

        let sources = SourceConfig {
            sources: vec![source],
        };

        grid.time = 0.25e-9;
        inject_sources(&mut grid, &sources, &params);

        assert_eq!(grid.cell(4, 4, 4).q_dot[3], 0.0, "inactive source should not inject");
    }
}
