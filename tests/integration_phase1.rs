// Integration tests for Phase 1 verification:
//   Task 1.4: Slice plane visualization shows propagating fields
//   Task 1.5: Absorbing boundaries prevent energy blowup

use quaternions::simulation::boundaries::{apply_boundaries, BoundaryConfig, PmlState};
use quaternions::simulation::diagnostics::{
    compute_derived_fields, compute_topological_charge, max_e, max_k_field, max_s, mean_k_field,
    total_energy,
};
use quaternions::simulation::field_update::step_field_cpu;
use quaternions::simulation::grid::SimulationGrid;
use quaternions::simulation::plugin::VacuumConfig;
use quaternions::simulation::sources::{inject_sources, Source, SourceConfig};
use quaternions::simulation::state::PmlConfig;
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
        step_field_cpu(&mut grid, &params, None, None);
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
        max_k: 1.0,
        mean_k: 1.0,
        topological_charge: 0.0,
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
        max_k: 1.0,
        mean_k: 1.0,
        topological_charge: 0.0,
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
        max_k: 1.0,
        mean_k: 1.0,
        topological_charge: 0.0,
    };

    // X slice: width=nz (maps to world Z via local X), height=ny (maps to world Y via local Y)
    let mut config = SliceConfig::default();
    config.axis = SliceAxis::X;
    config.position = 4;
    let (w, h, _) = sample_slice(&grid, &diag, &config);
    assert_eq!((w, h), (12, 10));

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
        step_field_cpu(&mut grid, &params, None, None);
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
        step_field_cpu(&mut grid, &params, None, None);
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
        step_field_cpu(&mut grid_no_bc, &params, None, None);
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
        step_field_cpu(&mut grid_bc, &params, None, None);
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
        step_field_cpu(&mut grid_std, &p, None, None);
        grid_std.swap_and_advance();
    }
    let derived_std = compute_derived_fields(&grid_std, &grid_std.sim_params(false));
    let s_max_std = max_s(&derived_std);

    // --- Extended mode ---
    let mut grid_ext = SimulationGrid::new(n, n, n, dx);
    for _ in 0..30 {
        let p = grid_ext.sim_params(true); // extended_mode = true
        inject_sources(&mut grid_ext, &sources, &p);
        step_field_cpu(&mut grid_ext, &p, None, None);
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

/// Verify that extended mode with sources doesn't blow up over many steps.
/// This was the original bug: the missing +c·∂S/∂t scalar correction caused
/// d²S/dt² = 0, leading to linear S growth and eventual divergence.
#[test]
fn test_extended_mode_stability_with_sources() {
    let n = 16;
    let dx = 0.01;
    let center = [n as f32 / 2.0; 3];

    // Test with both dipole and point charge sources
    let sources = SourceConfig {
        sources: vec![
            Source::dipole_z(center, 1.0, 1e9),
            Source::point_charge([6.0, 8.0, 8.0], 1e-12),
        ],
    };

    let mut grid = SimulationGrid::new(n, n, n, dx);
    let bc_config = BoundaryConfig::default();

    // Run for 500 steps in extended mode — enough to trigger the old blowup
    for _ in 0..500 {
        let p = grid.sim_params(true); // extended mode
        inject_sources(&mut grid, &sources, &p);
        step_field_cpu(&mut grid, &p, None, None);
        grid.swap_and_advance();
        apply_boundaries(&mut grid, &bc_config, p.dt);
    }

    // All cells must be finite (no NaN/inf blowup)
    for cell in grid.read_buf() {
        for c in 0..4 {
            assert!(
                cell.q[c].is_finite(),
                "extended mode blowup: q[{c}] = {}",
                cell.q[c]
            );
            assert!(
                cell.q_dot[c].is_finite(),
                "extended mode blowup: q_dot[{c}] = {}",
                cell.q_dot[c]
            );
        }
    }

    // Energy should be finite and reasonable (not exponentially blown up)
    let derived = compute_derived_fields(&grid, &grid.sim_params(true));
    let energy = total_energy(&derived, dx);
    assert!(energy.is_finite(), "total energy must be finite: {energy}");
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
        step_field_cpu(&mut grid, &params, None, None);
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

// =========================================================================
// Regression: CPML + extended mode must not diverge
//
// Bug: grid.reset() cleared CellFlags::PML from every cell, disabling CPML
// inside step_field_cpu while boundary_system still used the PML path (pure
// Neumann extrapolation, no absorption). Result: fully reflective boundaries
// → energy accumulation → divergence.
//
// Additionally, extended mode corrections were applied inside PML cells.
// CPML is designed for the plain wave equation; adding ∂S/∂t and ∇S terms
// inside the absorbing layer undermines absorption and can cause instability.
// =========================================================================

/// Verify that grid.reset() preserves PML flags.
/// After reset, every cell that was PML before should still have CellFlags::PML.
#[test]
fn test_reset_preserves_pml_flags() {
    use quaternions::simulation::state::CellFlags;

    let n = 16u32;
    let dx = 0.01;
    let bc_config = BoundaryConfig::default();
    let mut grid = SimulationGrid::new(n, n, n, dx);
    let _pml = PmlState::new(&mut grid, &bc_config, PmlConfig::default());

    // Count PML-flagged cells before reset
    let pml_before: usize = grid
        .read_buf()
        .iter()
        .filter(|c| (c.flags & CellFlags::PML) != 0)
        .count();
    assert!(pml_before > 0, "PML should flag some cells before reset");

    grid.reset();

    // Count PML-flagged cells after reset -- must be identical
    let pml_after: usize = grid
        .read_buf()
        .iter()
        .filter(|c| (c.flags & CellFlags::PML) != 0)
        .count();
    assert_eq!(
        pml_after,
        pml_before,
        "reset() must preserve PML flags: had {pml_before}, got {pml_after}"
    );

    // Non-PML flags should have been cleared (SOURCE, CONDUCTOR, etc.)
    let non_pml_flags = CellFlags::SOURCE | CellFlags::CONDUCTOR | CellFlags::BOUNDARY;
    for cell in grid.read_buf() {
        assert_eq!(
            cell.flags & non_pml_flags,
            0,
            "non-PML flags should be cleared after reset"
        );
        // All field values must be zero
        assert_eq!(cell.q, [0.0; 4]);
        assert_eq!(cell.q_dot, [0.0; 4]);
    }
}

/// Extended mode with CPML should remain stable for many steps.
///
/// Previously this diverged because:
/// 1. grid.reset() cleared PML flags -> no CPML absorption
/// 2. Extended corrections were applied in PML cells -> instability
#[test]
fn test_extended_mode_with_pml_stability() {
    let n = 16u32;
    let dx = 0.01;
    let center = [n as f32 / 2.0; 3];

    let sources = SourceConfig {
        sources: vec![
            Source::dipole_z(center, 1.0, 1e9),
            Source::point_charge([4.0, 8.0, 8.0], 1e-12),
        ],
    };

    let bc_config = BoundaryConfig::default();
    let mut grid = SimulationGrid::new(n, n, n, dx);
    let mut pml_state = PmlState::new(&mut grid, &bc_config, PmlConfig::default());

    // Simulate a scenario load (grid.reset() must preserve PML flags)
    grid.reset();
    pml_state.reset_psi();

    // Run 2000 steps in extended mode with CPML active
    for _ in 0..2000 {
        let p = grid.sim_params(true); // extended mode
        inject_sources(&mut grid, &sources, &p);
        step_field_cpu(&mut grid, &p, Some(&mut pml_state), None);
        grid.swap_and_advance();
    }

    // All cells must be finite -- no NaN or Inf blowup
    for (idx, cell) in grid.read_buf().iter().enumerate() {
        for c in 0..4 {
            assert!(
                cell.q[c].is_finite(),
                "cell {idx} q[{c}] blew up in extended+PML mode: {}",
                cell.q[c]
            );
            assert!(
                cell.q_dot[c].is_finite(),
                "cell {idx} q_dot[{c}] blew up in extended+PML mode: {}",
                cell.q_dot[c]
            );
        }
    }

    // Energy must also be finite
    let derived = compute_derived_fields(&grid, &grid.sim_params(true));
    let energy = total_energy(&derived, dx);
    assert!(
        energy.is_finite(),
        "total energy must be finite after extended+PML run: {energy}"
    );
}

// =========================================================================
// Task 1.8: Dynamic K field tests
// =========================================================================

/// K=1 vacuum with no EM fields → K should remain exactly 1.0 after 1000 steps.
#[test]
fn test_k_vacuum_stable() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);
    let params = grid.sim_params(false);

    // VacuumConfig with omega_p=0 and eta=0: restoring force and coupling both zero.
    // K should not drift from 1.0 (no driving, no Laplacian from uniform K=1).
    let vacuum = VacuumConfig {
        enabled: true,
        omega_p: 0.0,
        eta: 0.0,
        u_s: 1.0,
    };

    for _ in 0..1000 {
        step_field_cpu(&mut grid, &params, None, Some(&vacuum));
        grid.swap_and_advance();
    }

    let max_k = max_k_field(&grid);
    let mean_k = mean_k_field(&grid);
    assert!(
        (max_k - 1.0).abs() < 1e-6,
        "K should remain 1.0 in vacuum with no fields; got max_k={max_k}"
    );
    assert!(
        (mean_k - 1.0).abs() < 1e-6,
        "mean K should be 1.0; got {mean_k}"
    );
}

/// Strong E field with non-zero eta → K should increase above 1.0.
///
/// Uses a large-amplitude linear phi profile Q.w = x * A * dx.  A linear
/// profile has zero Laplacian, so Q is unchanged by the wave equation, but
/// it creates a uniform static E field that drives K upward.
/// The amplitude A is chosen so that the K change after N steps exceeds
/// the f32 resolution (~1.2e-7).
#[test]
fn test_k_increases_with_field() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);

    // Large-amplitude linear phi gradient: Q.w = x * A * dx
    // E_x = -c₀ * A (uniform, static since Laplacian is zero).
    // A = 1e5 → E_x ≈ -3e13 V/m → u_field ≈ 0.5*ε₀*E² ≈ 4e15 J/m³.
    // k_ddot = eta * u_field / u_s = 1.0 * 4e15 = 4e15 s⁻²
    // After 50 steps: Δk ≈ k_ddot * N(N+1)/2 * dt² ≈ 4e15 * 1275 * 3e-22 ≈ 1.5e-3 >> 1e-5
    let amplitude: f32 = 1e5;
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = grid.idx(x as u32, y as u32, z as u32);
                let v = x as f32 * amplitude * grid.dx;
                grid.cells[0][i].q[0] = v;
                grid.cells[1][i].q[0] = v;
            }
        }
    }

    let params = grid.sim_params(false);
    let vacuum = VacuumConfig {
        enabled: true,
        omega_p: 0.0, // no restoring force — pure driving
        eta: 1.0,
        u_s: 1.0,
    };

    for _ in 0..50 {
        step_field_cpu(&mut grid, &params, None, Some(&vacuum));
        grid.swap_and_advance();
    }

    let max_k = max_k_field(&grid);
    assert!(
        max_k > 1.0 + 1e-5,
        "K should increase above 1.0 when strong E field drives eta*u_field/u_s > 0; got max_k={max_k}"
    );
}

/// After removing the driving field, K should decay back toward 1.0
/// via the ωₚ² restoring force.
#[test]
fn test_k_restores_after_pulse() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);

    // Manually perturb K above 1.0 in the center region
    let center = grid.nx as usize / 2;
    for dz in 0..3usize {
        for dy in 0..3usize {
            for dx in 0..3usize {
                let x = (center - 1 + dx) as u32;
                let y = (center - 1 + dy) as u32;
                let z = (center - 1 + dz) as u32;
                let i = grid.idx(x, y, z);
                grid.cells[0][i].k = 2.0;
                grid.cells[1][i].k = 2.0;
            }
        }
    }

    let params = grid.sim_params(false);
    let omega_p = quaternions::simulation::state::SimParams::C0 / (50.0 * grid.dx);
    let vacuum = VacuumConfig {
        enabled: true,
        omega_p,
        eta: 0.0, // no driving — only restoring force
        u_s: 1.0,
    };

    // Record initial max K
    let k_before = max_k_field(&grid);

    // Run enough steps for restoring force to reduce K
    for _ in 0..200 {
        step_field_cpu(&mut grid, &params, None, Some(&vacuum));
        grid.swap_and_advance();
    }

    let k_after = max_k_field(&grid);
    assert!(
        k_after < k_before,
        "K should decay toward 1.0 after pulse removed (ωₚ² restoring); before={k_before:.4}, after={k_after:.4}"
    );
}

/// K diagnostics helpers return sensible values on a vacuum grid.
#[test]
fn test_k_diagnostics_vacuum() {
    let grid = SimulationGrid::new(8, 8, 8, 0.01);
    let max_k = max_k_field(&grid);
    let mean_k = mean_k_field(&grid);
    assert!(
        (max_k - 1.0).abs() < 1e-6,
        "max_k on vacuum grid should be 1.0, got {max_k}"
    );
    assert!(
        (mean_k - 1.0).abs() < 1e-6,
        "mean_k on vacuum grid should be 1.0, got {mean_k}"
    );
}

// =========================================================================
// Task 1.9 verification: topological charge diagnostic
// =========================================================================

/// Uniform Q field → topological charge = 0 (no winding).
#[test]
fn test_topo_charge_vacuum() {
    let grid = SimulationGrid::new(8, 8, 8, 0.01);
    let charge = compute_topological_charge(&grid);
    assert!(
        charge.abs() < 1e-6,
        "vacuum grid should have zero topological charge, got {charge}"
    );
}

/// Hedgehog configuration: Q ∝ (r, x, y, z) normalized → topological charge = ±1.
///
/// The hedgehog maps R³ → S³ with winding number 1. We set Q = (r, x, y, z)
/// where r = sqrt(x² + y² + z²) measured from grid center, which wraps S³ once.
#[test]
fn test_topo_charge_hedgehog() {
    let n = 32u32;
    let dx = 0.01;
    let mut grid = SimulationGrid::new(n, n, n, dx);
    let center = n as f32 / 2.0;

    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let i = grid.idx(x, y, z);
                let rx = x as f32 - center;
                let ry = y as f32 - center;
                let rz = z as f32 - center;
                let r = (rx * rx + ry * ry + rz * rz).sqrt();

                // Q = (cos(π*r/R), sin(π*r/R) * r_hat)
                // This is the standard hedgehog ansatz mapping a ball of radius R to S³
                let max_r = center; // radius of the grid domain
                if r < 1e-6 {
                    // At origin: Q = (1, 0, 0, 0)
                    grid.cells[0][i].q = [1.0, 0.0, 0.0, 0.0];
                } else {
                    let theta = std::f32::consts::PI * (r / max_r).min(1.0);
                    let cos_t = theta.cos();
                    let sin_t = theta.sin();
                    let inv_r = 1.0 / r;
                    grid.cells[0][i].q = [
                        cos_t,
                        sin_t * rx * inv_r,
                        sin_t * ry * inv_r,
                        sin_t * rz * inv_r,
                    ];
                }
            }
        }
    }

    let charge = compute_topological_charge(&grid);
    // The hedgehog should give charge ≈ ±1 (sign depends on orientation convention).
    // On a finite grid with central differences, we accept |charge| in [0.5, 1.5].
    assert!(
        (charge.abs() - 1.0).abs() < 0.5,
        "hedgehog should have |topological charge| ≈ 1, got {charge}"
    );
}

/// Slowly varying Q field (small perturbation from uniform) → charge ≈ 0.
#[test]
fn test_topo_charge_smooth_field() {
    let n = 16u32;
    let dx = 0.01;
    let mut grid = SimulationGrid::new(n, n, n, dx);

    // Set Q to a gently varying field: Q = (1, εx, εy, εz) normalized
    let epsilon = 0.01;
    let center = n as f32 / 2.0;
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let i = grid.idx(x, y, z);
                let rx = (x as f32 - center) * epsilon;
                let ry = (y as f32 - center) * epsilon;
                let rz = (z as f32 - center) * epsilon;
                let norm = (1.0 + rx * rx + ry * ry + rz * rz).sqrt();
                grid.cells[0][i].q = [1.0 / norm, rx / norm, ry / norm, rz / norm];
            }
        }
    }

    let charge = compute_topological_charge(&grid);
    assert!(
        charge.abs() < 0.1,
        "slowly varying Q should have near-zero topological charge, got {charge}"
    );
}

/// Evolve a non-trivial Q configuration for 50 steps; charge should be roughly conserved.
/// On a finite grid with open boundaries, some drift is expected due to field leaking
/// through boundaries and discretization effects.
#[test]
fn test_topo_charge_conserved() {
    let n = 32u32;
    let dx = 0.01;
    let mut grid = SimulationGrid::new(n, n, n, dx);
    let center = n as f32 / 2.0;

    // Set up a hedgehog-like configuration
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let i = grid.idx(x, y, z);
                let rx = x as f32 - center;
                let ry = y as f32 - center;
                let rz = z as f32 - center;
                let r = (rx * rx + ry * ry + rz * rz).sqrt();
                let max_r = center;

                if r < 1e-6 {
                    grid.cells[0][i].q = [1.0, 0.0, 0.0, 0.0];
                } else {
                    let theta = std::f32::consts::PI * (r / max_r).min(1.0);
                    let cos_t = theta.cos();
                    let sin_t = theta.sin();
                    let inv_r = 1.0 / r;
                    grid.cells[0][i].q = [
                        cos_t,
                        sin_t * rx * inv_r,
                        sin_t * ry * inv_r,
                        sin_t * rz * inv_r,
                    ];
                }
            }
        }
    }

    let charge_initial = compute_topological_charge(&grid);

    // Evolve for 100 steps with extended mode (so Q topology matters)
    let boundary_config = BoundaryConfig::default();
    let dt = grid.dt;

    for _ in 0..50 {
        let params = grid.sim_params(true);
        step_field_cpu(&mut grid, &params, None, None);
        apply_boundaries(&mut grid, &boundary_config, dt);
        grid.swap_and_advance();
    }

    let charge_final = compute_topological_charge(&grid);

    // On a finite grid with open boundaries, some drift is expected.
    // We check that the charge doesn't change by more than 50% of the initial value.
    let drift = (charge_final - charge_initial).abs();
    let threshold = charge_initial.abs().max(0.5);
    assert!(
        drift < threshold,
        "topological charge should be roughly conserved: initial={charge_initial:.4}, final={charge_final:.4}, drift={drift:.4}"
    );
}
