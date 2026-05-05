// Bifilar coil scenario: scalar wave transmitter/receiver.
//
// A bifilar coil is wound so that adjacent wires carry current in opposite
// directions. The B fields from adjacent wires cancel (to first order), but
// the A (vector potential) fields add. This creates a source with strong A
// but weak B — ideal for exciting scalar/longitudinal modes because the
// transverse (B-based) radiation is suppressed while the potential-based
// effects survive.
//
// Two sub-scenarios:
//   1. Single bifilar coil: demonstrates B cancellation + S field excitation
//   2. Bifilar pair: transmitter + receiver separated by distance, measure
//      S field coupling in extended mode vs standard mode
//
// NOTE: At 64³ / 0.5m domain, dx ≈ 8mm is too coarse for realistic wire
// spacing (~1mm). We use an idealized geometry: two parallel current lines
// separated by 1 cell (the minimum resolvable spacing). This captures the
// essential physics — B cancellation at distances >> wire spacing — while
// remaining grid-resolvable.

use crate::simulation::diagnostics::{Probe, ProbeField, ProbeSet};
use crate::simulation::grid::SimulationGrid;
use crate::simulation::sources::{Source, SourceConfig, SourceType};
use crate::simulation::state::SimParams;

/// Create a bifilar source set: two opposing-current dipoles + scalar drive.
///
/// The scalar drive represents the net ∂(∇·A)/∂t produced by the opposing
/// currents. In a real bifilar geometry, the B-cancellation channels energy
/// into the scalar/longitudinal mode; we model this explicitly with a
/// ScalarDrive source at the midpoint.
fn bifilar_sources(
    center_x: f32,
    center_y: f32,
    center_z: f32,
    amplitude: f32,
    frequency: f32,
    driven: bool,
) -> Vec<Source> {
    let amp = if driven { amplitude } else { 0.0 };
    let s_amp = if driven { amplitude * 0.5 } else { 0.0 };

    vec![
        // Wire 1: z-directed current
        Source {
            position: [center_x - 0.5, center_y, center_z],
            source_type: SourceType::OscillatingDipole { axis: 2 },
            amplitude: amp,
            frequency,
            phase: 0.0,
            active: true,
        },
        // Wire 2: opposite z-directed current
        Source {
            position: [center_x + 0.5, center_y, center_z],
            source_type: SourceType::OscillatingDipole { axis: 2 },
            amplitude: -amp,
            frequency,
            phase: 0.0,
            active: true,
        },
        // Scalar drive at midpoint — models the ∂(∇·A)/∂t excitation
        Source {
            position: [center_x, center_y, center_z],
            source_type: SourceType::ScalarDrive,
            amplitude: s_amp,
            frequency,
            phase: 0.0,
            active: true,
        },
    ]
}

/// Apply the single bifilar coil scenario to an existing grid.
///
/// Places two parallel z-directed current lines at the grid center,
/// separated by 1 cell in the x-direction, carrying current in opposite
/// directions. The opposing currents produce:
///   - B fields that largely cancel in the far field
///   - A fields that add (both wires contribute A in the same direction
///     because A is proportional to current direction, and the geometric
///     arrangement means the vector potentials reinforce)
///   - Strong ∂A/∂t oscillation → drives S = (1/c²)∂φ/∂t + ∇·A
///
/// In extended mode, the S field propagates as a scalar longitudinal wave.
/// In standard mode, S is constrained to zero (Lorenz gauge) and the
/// B-cancellation means very little transverse radiation escapes.
pub fn apply_bifilar_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
) {
    grid.reset();
    sources.sources.clear();

    let c0 = SimParams::C0;
    let dx = grid.dx;

    // Frequency: wavelength ≈ 16 cells — well-resolved on grid
    let wavelength_cells = 16.0;
    let frequency = c0 / (wavelength_cells * dx);

    // Amplitude: moderately strong to produce visible S field
    let amplitude = 5.0;

    let cx = grid.nx as f32 / 2.0;
    let cy = grid.ny as f32 / 2.0;
    let cz = grid.nz as f32 / 2.0;

    sources.sources.extend(bifilar_sources(cx, cy, cz, amplitude, frequency, true));
}

/// Apply the bifilar pair scenario: transmitter + receiver.
///
/// Places a bifilar transmitter at 1/4 of the grid and a bifilar receiver
/// at 3/4. Both are z-directed current pairs separated by 1 cell in x.
///
/// The transmitter is driven; the receiver is passive (amplitude = 0).
/// In extended mode, the S field excited by the transmitter propagates
/// to the receiver location. Measurement probes (added in task 2.4) will
/// show S-field signal at the receiver in extended mode but not in standard.
///
/// For now, this sets up the geometry. Quantitative measurement requires
/// the probe system from task 2.4.
pub fn apply_bifilar_pair_scenario(
    grid: &mut SimulationGrid,
    sources: &mut SourceConfig,
) {
    grid.reset();
    sources.sources.clear();

    let c0 = SimParams::C0;
    let dx = grid.dx;

    // Frequency: wavelength ≈ 16 cells
    let wavelength_cells = 16.0;
    let frequency = c0 / (wavelength_cells * dx);
    let amplitude = 5.0;

    let cy = grid.ny as f32 / 2.0;
    let cz = grid.nz as f32 / 2.0;

    // Transmitter at x = nx/4 (driven)
    let tx_x = grid.nx as f32 / 4.0;
    sources.sources.extend(bifilar_sources(tx_x, cy, cz, amplitude, frequency, true));

    // Receiver at x = 3*nx/4 (passive — markers only)
    let rx_x = 3.0 * grid.nx as f32 / 4.0;
    sources.sources.extend(bifilar_sources(rx_x, cy, cz, amplitude, frequency, false));
}

/// Populate a `ProbeSet` with three S-field probes aligned with the
/// transmitter / midpoint / receiver of the bifilar pair scenario.
///
/// The probes share the scenario's y/z center and are spaced along x at the
/// same locations where `apply_bifilar_pair_scenario` places its coils. After
/// running the simulation, the time series can be compared to measure the
/// propagation delay of the scalar wave and thus its phase velocity.
pub fn install_bifilar_pair_probes(grid: &SimulationGrid, probes: &mut ProbeSet) {
    probes.clear();

    let cy = (grid.ny / 2) as u32;
    let cz = (grid.nz / 2) as u32;
    let tx_x = (grid.nx / 4) as u32;
    let mid_x = (grid.nx / 2) as u32;
    let rx_x = (3 * grid.nx / 4) as u32;

    probes.push(Probe::new("TX (S)", [tx_x, cy, cz], ProbeField::SField));
    probes.push(Probe::new("Mid (S)", [mid_x, cy, cz], ProbeField::SField));
    probes.push(Probe::new("RX (S)", [rx_x, cy, cz], ProbeField::SField));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::grid::SimulationGrid;
    use crate::simulation::field_update::step_field_cpu;
    use crate::simulation::sources::inject_sources;
    use crate::simulation::boundaries::{apply_boundaries, BoundaryConfig};
    use crate::simulation::diagnostics::{compute_derived_fields, max_b};

    #[test]
    fn test_bifilar_sources_created() {
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();

        apply_bifilar_scenario(&mut grid, &mut sources);

        // 3 sources: 2 opposing dipoles + 1 scalar drive
        assert_eq!(sources.sources.len(), 3, "should have 3 sources (2 wires + scalar drive)");
        // Opposing dipole amplitudes
        let a0 = sources.sources[0].amplitude;
        let a1 = sources.sources[1].amplitude;
        assert!(
            (a0 + a1).abs() < 1e-10,
            "bifilar wire sources should have opposite amplitudes"
        );
        // Scalar drive source
        assert_eq!(sources.sources[2].source_type, SourceType::ScalarDrive);
    }

    #[test]
    fn test_bifilar_pair_sources_created() {
        let mut grid = SimulationGrid::new(64, 64, 64, 0.01);
        let mut sources = SourceConfig::default();

        apply_bifilar_pair_scenario(&mut grid, &mut sources);

        // 6 sources: 3 per coil (tx driven, rx passive)
        assert_eq!(sources.sources.len(), 6, "should have 6 sources (3 tx + 3 rx)");
        // Transmitter dipole pair has opposite amplitudes
        let a0 = sources.sources[0].amplitude;
        let a1 = sources.sources[1].amplitude;
        assert!(
            (a0 + a1).abs() < 1e-10,
            "transmitter pair should have opposite amplitudes"
        );
        // Receiver sources are passive (zero amplitude)
        assert_eq!(sources.sources[3].amplitude, 0.0);
        assert_eq!(sources.sources[4].amplitude, 0.0);
        assert_eq!(sources.sources[5].amplitude, 0.0);
    }

    #[test]
    fn test_bifilar_b_cancellation() {
        // Run a bifilar coil for some steps and verify that B is weaker
        // than what a single dipole of the same total amplitude would produce.
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        let boundaries = BoundaryConfig::default();

        apply_bifilar_scenario(&mut grid, &mut sources);

        let dt = grid.dt;
        for _ in 0..50 {
            let params = grid.sim_params(false);
            inject_sources(&mut grid, &sources, &params);
            step_field_cpu(&mut grid, &params, None, None);
            grid.swap_and_advance();
            apply_boundaries(&mut grid, &boundaries, dt);
        }

        let params = grid.sim_params(false);
        let derived = compute_derived_fields(&grid, &params);
        let b_bifilar = max_b(&derived);

        // Compare with a single dipole of double amplitude
        let mut grid2 = SimulationGrid::new(32, 32, 32, 0.01);
        let center = 16.0;
        let c0 = SimParams::C0;
        let freq = c0 / (16.0 * 0.01);
        let single_sources = SourceConfig {
            sources: vec![Source {
                position: [center, center, center],
                source_type: SourceType::OscillatingDipole { axis: 2 },
                amplitude: 10.0, // same total |amplitude| as bifilar pair
                frequency: freq,
                phase: 0.0,
                active: true,
            }],
        };

        for _ in 0..50 {
            let params = grid2.sim_params(false);
            inject_sources(&mut grid2, &single_sources, &params);
            step_field_cpu(&mut grid2, &params, None, None);
            grid2.swap_and_advance();
            apply_boundaries(&mut grid2, &boundaries, dt);
        }

        let params2 = grid2.sim_params(false);
        let derived2 = compute_derived_fields(&grid2, &params2);
        let b_single = max_b(&derived2);

        // Bifilar B should be significantly less than single dipole B
        // (not zero due to coarse grid — near-field doesn't perfectly cancel)
        assert!(
            b_bifilar < b_single,
            "bifilar B ({:.3e}) should be less than single dipole B ({:.3e})",
            b_bifilar,
            b_single,
        );
    }

    #[test]
    fn test_bifilar_extended_produces_s_field() {
        // In extended mode, the bifilar coil should excite the S field.
        // Use more steps to allow S to build up from the ∇·A source term.
        let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
        let mut sources = SourceConfig::default();
        let boundaries = BoundaryConfig::default();

        apply_bifilar_scenario(&mut grid, &mut sources);

        let dt = grid.dt;
        for _ in 0..300 {
            let params = grid.sim_params(true); // extended mode
            inject_sources(&mut grid, &sources, &params);
            step_field_cpu(&mut grid, &params, None, None);
            grid.swap_and_advance();
            apply_boundaries(&mut grid, &boundaries, dt);
        }

        // Check the raw s_field buffer directly (not just derived which
        // only covers interior cells)
        let s_max_raw: f32 = grid.s_read().iter()
            .map(|s| s.abs())
            .fold(0.0f32, f32::max);

        assert!(
            s_max_raw > 0.0,
            "S field should be nonzero in extended mode with bifilar source (raw max: {:.3e})",
            s_max_raw,
        );
    }
}
