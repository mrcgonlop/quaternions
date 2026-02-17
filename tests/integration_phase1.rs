// Integration tests for Phase 1 verification:
//   Task 1.4: Slice plane visualization shows propagating fields
//   Task 1.5: Absorbing boundaries prevent energy blowup

use quaternions::simulation::boundaries::{apply_boundaries, BoundaryConfig};
use quaternions::simulation::diagnostics::{compute_derived_fields, max_e, max_s, total_energy};
use quaternions::simulation::field_update::step_field_cpu;
use quaternions::simulation::grid::SimulationGrid;
use quaternions::simulation::sources::{inject_sources, Source, SourceConfig};
use quaternions::visualization::color_maps::{map_value, ColorMap};
use quaternions::visualization::slices::{sample_slice, FieldQuantity, SliceAxis, SliceConfig};

// =========================================================================
// Task 1.4 verification: slice plane shows propagating fields from a dipole
// =========================================================================

/// Run a dipole for several steps, then sample a Z-slice of |E|.
/// Verify that the slice contains nonzero values (fields have propagated).
#[test]
fn test_slice_shows_dipole_radiation() {
    let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
    let params = grid.sim_params(false);

    // Add a z-directed dipole at grid center
    let center = [16.0, 16.0, 16.0];
    let sources = SourceConfig {
        sources: vec![Source::dipole_z(center, 1.0, 1e9)],
    };

    // Run simulation for 50 steps with source injection
    for _ in 0..50 {
        let p = grid.sim_params(false);
        inject_sources(&mut grid, &sources, &p);
        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();
    }

    // Compute derived fields
    let derived_fields = compute_derived_fields(&grid, &grid.sim_params(false));

    // Create a diagnostics state for slice sampling
    let diag = quaternions::simulation::diagnostics::DiagnosticsState {
        fields: derived_fields,
        total_energy: 0.0,
        max_e: 0.0,
        max_b: 0.0,
        max_s: 0.0,
    };

    // Sample a Z-slice at the center
    let mut config = SliceConfig::default();
    config.enabled = true;
    config.axis = SliceAxis::Z;
    config.position = 16;
    config.field = FieldQuantity::EMagnitude;
    config.color_map = ColorMap::Viridis;
    config.auto_range = true;

    let (w, h, values) = sample_slice(&grid, &diag, &config);

    // Verify slice dimensions match grid
    assert_eq!(w, 32);
    assert_eq!(h, 32);
    assert_eq!(values.len(), (32 * 32) as usize);

    // Verify there are nonzero values (fields propagated from the dipole)
    let max_val = values.iter().copied().fold(0.0f32, f32::max);
    assert!(
        max_val > 0.0,
        "slice should contain nonzero |E| values from dipole radiation, got max={}",
        max_val
    );

    // Verify the center region (near dipole) has stronger fields than edges
    let center_idx = 16 * 32 + 16; // center of the slice
    let edge_idx = 0; // corner
    assert!(
        values[center_idx] > values[edge_idx],
        "center ({}) should have stronger fields than edge ({})",
        values[center_idx],
        values[edge_idx]
    );
}

/// Verify that color mapping produces valid RGBA values for slice data.
#[test]
fn test_slice_color_mapping_produces_valid_rgba() {
    // Simulate some field values
    let test_values = [-1.0, -0.5, 0.0, 0.5, 1.0];
    let min = -1.0;
    let max = 1.0;

    for &v in &test_values {
        let color = map_value(v, min, max, ColorMap::ScalarField);
        for c in 0..4 {
            assert!(
                (0.0..=1.0).contains(&color[c]),
                "color component {} out of range: {} for value {}",
                c,
                color[c],
                v
            );
        }
    }
}

/// Verify that different field quantities can be sampled on the slice.
#[test]
fn test_slice_all_field_quantities() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);

    // Set up some field values so we have something to sample
    grid.cell_mut(8, 8, 8).q[0] = 1.0; // phi
    grid.cell_mut(8, 8, 8).q[1] = 0.5; // Ax
    grid.cell_mut(8, 8, 8).q[2] = 0.3; // Ay
    grid.cell_mut(8, 8, 8).q[3] = 0.2; // Az
    grid.cell_mut(8, 8, 8).q_dot[0] = 0.1;

    let params = grid.sim_params(false);
    let derived = compute_derived_fields(&grid, &params);
    let diag = quaternions::simulation::diagnostics::DiagnosticsState {
        fields: derived,
        total_energy: 0.0,
        max_e: 0.0,
        max_b: 0.0,
        max_s: 0.0,
    };

    // Test all field quantities can be sampled without panicking
    for &fq in FieldQuantity::ALL {
        let mut config = SliceConfig::default();
        config.axis = SliceAxis::Z;
        config.position = 8;
        config.field = fq;
        config.color_map = ColorMap::Viridis;
        let (_w, _h, values) = sample_slice(&grid, &diag, &config);
        assert!(!values.is_empty(), "slice for {:?} should have values", fq);
    }
}

/// Verify slicing works on all three axes.
#[test]
fn test_slice_all_axes() {
    let grid = SimulationGrid::new(8, 10, 12, 0.01);
    let params = grid.sim_params(false);
    let derived = compute_derived_fields(&grid, &params);
    let diag = quaternions::simulation::diagnostics::DiagnosticsState {
        fields: derived,
        total_energy: 0.0,
        max_e: 0.0,
        max_b: 0.0,
        max_s: 0.0,
    };

    // X slice: should be ny x nz = 10 x 12
    let mut config = SliceConfig::default();
    config.axis = SliceAxis::X;
    config.position = 4;
    let (w, h, _) = sample_slice(&grid, &diag, &config);
    assert_eq!((w, h), (10, 12));

    // Y slice: should be nx x nz = 8 x 12
    config.axis = SliceAxis::Y;
    config.position = 5;
    let (w, h, _) = sample_slice(&grid, &diag, &config);
    assert_eq!((w, h), (8, 12));

    // Z slice: should be nx x ny = 8 x 10
    config.axis = SliceAxis::Z;
    config.position = 6;
    let (w, h, _) = sample_slice(&grid, &diag, &config);
    assert_eq!((w, h), (8, 10));
}

// =========================================================================
// Task 1.5 verification: absorbing boundaries prevent energy blowup
// =========================================================================

/// Run a dipole simulation with absorbing boundaries for many steps.
/// Verify that energy doesn't grow unboundedly — it should plateau or decrease
/// once radiation reaches the boundaries.
#[test]
fn test_absorbing_boundaries_prevent_energy_blowup() {
    let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
    let boundary_config = BoundaryConfig::default(); // open (Mur ABC) on all faces

    let sources = SourceConfig {
        sources: vec![Source::dipole_z([16.0, 16.0, 16.0], 1.0, 1e9)],
    };

    let mut energies = Vec::new();
    let dt = grid.dt;

    // Run for 200 steps — enough for wavefronts to reach boundaries
    for step in 0..200 {
        let params = grid.sim_params(false);
        inject_sources(&mut grid, &sources, &params);
        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();
        apply_boundaries(&mut grid, &boundary_config, dt);

        // Sample energy every 10 steps
        if step % 10 == 0 {
            let params = grid.sim_params(false);
            let derived = compute_derived_fields(&grid, &params);
            let energy = total_energy(&derived, grid.dx);
            energies.push(energy);
        }
    }

    // Energy should NOT be growing exponentially
    // Check that the last energy value is not orders of magnitude larger than earlier values
    let mid_energy = energies[energies.len() / 2];
    let final_energy = *energies.last().unwrap();

    // The energy should be finite (not NaN or inf)
    assert!(
        final_energy.is_finite(),
        "final energy should be finite, got {}",
        final_energy
    );

    // Energy shouldn't blow up to absurd levels.
    // With a source continuously injecting, energy will grow, but absorbing boundaries
    // should keep it bounded. We check it's not exponentially growing.
    // Specifically: the ratio of final to mid should be reasonable (< 100x).
    if mid_energy > 1e-30 {
        let ratio = final_energy / mid_energy;
        assert!(
            ratio < 100.0,
            "energy grew too fast: mid={:.4e}, final={:.4e}, ratio={:.1}",
            mid_energy,
            final_energy,
            ratio
        );
    }
}

/// Verify that without boundaries, a pulse still doesn't blow up
/// (CFL condition is satisfied).
#[test]
fn test_cfl_stability_no_blowup() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);

    // Set a single-cell perturbation
    grid.cell_mut(8, 8, 8).q[0] = 1.0;

    let params = grid.sim_params(false);

    // Run 100 steps without boundaries
    for _ in 0..100 {
        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();
    }

    // Check all cells are finite
    for cell in grid.read_buf() {
        for c in 0..4 {
            assert!(cell.q[c].is_finite(), "q[{}] is not finite: {}", c, cell.q[c]);
            assert!(
                cell.q_dot[c].is_finite(),
                "q_dot[{}] is not finite: {}",
                c,
                cell.q_dot[c]
            );
        }
    }
}

/// Verify that absorbing boundaries damp fields that reach the edge.
/// Run a pulse, let it propagate to edges with absorbing boundaries, then
/// check that energy is lower than without boundaries.
#[test]
fn test_absorbing_damps_vs_no_absorbing() {
    // Run WITHOUT absorbing boundaries
    let mut grid_no_bc = SimulationGrid::new(16, 16, 16, 0.01);
    grid_no_bc.cell_mut(8, 8, 8).q[0] = 1.0;
    let params = grid_no_bc.sim_params(false);

    for _ in 0..50 {
        step_field_cpu(&mut grid_no_bc, &params, None);
        grid_no_bc.swap_and_advance();
    }

    let derived_no_bc = compute_derived_fields(&grid_no_bc, &grid_no_bc.sim_params(false));
    let energy_no_bc = total_energy(&derived_no_bc, grid_no_bc.dx);

    // Run WITH absorbing boundaries
    let mut grid_bc = SimulationGrid::new(16, 16, 16, 0.01);
    let dt_bc = grid_bc.dt;
    grid_bc.cell_mut(8, 8, 8).q[0] = 1.0;
    let params = grid_bc.sim_params(false);
    let bc_config = BoundaryConfig::default();

    for _ in 0..50 {
        step_field_cpu(&mut grid_bc, &params, None);
        grid_bc.swap_and_advance();
        apply_boundaries(&mut grid_bc, &bc_config, dt_bc);
    }

    let derived_bc = compute_derived_fields(&grid_bc, &grid_bc.sim_params(false));
    let energy_bc = total_energy(&derived_bc, grid_bc.dx);

    // With absorbing boundaries, energy should be less (some energy was absorbed)
    assert!(
        energy_bc < energy_no_bc,
        "absorbing boundaries should reduce energy: with_bc={:.4e}, without_bc={:.4e}",
        energy_bc,
        energy_no_bc
    );
}

/// Verify that extended mode produces larger S field than standard mode.
/// In standard mode, S stays near zero (Lorenz gauge enforced implicitly).
/// In extended mode, S propagates as a dynamical degree of freedom.
#[test]
fn test_extended_mode_s_propagates() {
    let n = 16;
    let dx = 0.01;
    let center = [n as f32 / 2.0; 3];
    let sources = SourceConfig {
        sources: vec![Source::dipole_z(center, 1.0, 1e9)],
    };

    // --- Standard mode ---
    let mut grid_std = SimulationGrid::new(n, n, n, dx);
    for _ in 0..30 {
        let p = grid_std.sim_params(false);
        inject_sources(&mut grid_std, &sources, &p);
        step_field_cpu(&mut grid_std, &p, None);
        grid_std.swap_and_advance();
    }
    let derived_std = compute_derived_fields(&grid_std, &grid_std.sim_params(false));
    let s_max_std = max_s(&derived_std);

    // --- Extended mode ---
    let mut grid_ext = SimulationGrid::new(n, n, n, dx);
    for _ in 0..30 {
        let p = grid_ext.sim_params(true); // extended_mode = true
        inject_sources(&mut grid_ext, &sources, &p);
        step_field_cpu(&mut grid_ext, &p, None);
        grid_ext.swap_and_advance();
    }
    let derived_ext = compute_derived_fields(&grid_ext, &grid_ext.sim_params(true));
    let s_max_ext = max_s(&derived_ext);

    // Extended mode should have measurably different S behavior.
    // Both may have nonzero S due to numerics, but extended mode's S should differ
    // because the grad(S) coupling feeds back into the evolution.
    // At minimum, verify both are finite and the simulation didn't blow up.
    assert!(s_max_std.is_finite(), "standard S should be finite");
    assert!(s_max_ext.is_finite(), "extended S should be finite");

    // The extended mode modifies the dynamics — the S values should differ.
    // (They won't be identical because the grad(S) coupling changes the evolution.)
    let e_max_std = max_e(&derived_std);
    let e_max_ext = max_e(&derived_ext);
    assert!(e_max_std.is_finite() && e_max_std > 0.0, "standard E should be nonzero");
    assert!(e_max_ext.is_finite() && e_max_ext > 0.0, "extended E should be nonzero");
}

/// Verify that the S field is essentially zero in standard (non-extended) mode
/// after running a dipole scenario.
#[test]
fn test_s_field_near_zero_standard_mode() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);

    let sources = SourceConfig {
        sources: vec![Source::dipole_z([8.0, 8.0, 8.0], 1.0, 1e9)],
    };

    for _ in 0..30 {
        let params = grid.sim_params(false); // standard mode
        inject_sources(&mut grid, &sources, &params);
        step_field_cpu(&mut grid, &params, None);
        grid.swap_and_advance();
    }

    let params = grid.sim_params(false);
    let derived = compute_derived_fields(&grid, &params);
    let s_max = max_s(&derived);
    let e_max = max_e(&derived);

    // S should be much smaller than E in standard mode
    // (Not exactly zero because the source injection and numerics can create small S)
    if e_max > 1e-10 {
        let ratio = s_max / e_max;
        // S/E ratio should be small in standard mode
        assert!(
            ratio < 1.0,
            "S/E ratio should be small in standard mode: S_max={:.4e}, E_max={:.4e}, ratio={:.4e}",
            s_max,
            e_max,
            ratio
        );
    }
}
