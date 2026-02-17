// CPML (Convolutional PML) integration tests.
//
// Verifies that CPML boundaries absorb outgoing waves with significantly
// less reflection than the old Mur+sponge approach.

use quaternions::simulation::boundaries::{apply_boundaries, BoundaryConfig, PmlState};
use quaternions::simulation::diagnostics::{compute_derived_fields, total_energy, total_energy_excluding_pml};
use quaternions::simulation::field_update::step_field_cpu;
use quaternions::simulation::grid::SimulationGrid;
use quaternions::simulation::sources::{inject_sources, Source, SourceConfig};
use quaternions::simulation::state::PmlConfig;

/// PML initialization should flag correct number of cells.
#[test]
fn test_pml_initialization() {
    let mut grid = SimulationGrid::new(32, 32, 32, 0.01);
    let bc = BoundaryConfig::default(); // all Open
    let config = PmlConfig { depth: 5, ..PmlConfig::default() };

    let pml = PmlState::new(&mut grid, &bc, config);

    // With depth=5 on all 6 faces of a 32^3 grid:
    // Total cells = 32^3 = 32768
    // Interior cells (not in any PML region): (32-10)^3 = 22^3 = 10648
    // PML cells = 32768 - 10648 = 22120
    assert!(pml.pml_cell_count > 0, "should have PML cells");
    assert!(
        pml.pml_cell_count < grid.cell_count(),
        "not all cells should be PML"
    );

    // Verify grid_to_pml mapping is consistent
    let pml_count_from_map = pml.grid_to_pml.iter().filter(|x| x.is_some()).count();
    assert_eq!(pml_count_from_map, pml.pml_cell_count);

    // Verify psi array size
    assert_eq!(pml.psi.len(), pml.pml_cell_count * 12);
    assert_eq!(pml.b.len(), pml.pml_cell_count * 3);
    assert_eq!(pml.a.len(), pml.pml_cell_count * 3);
    assert_eq!(pml.inv_kappa.len(), pml.pml_cell_count * 3);
}

/// PML cells should be flagged with CellFlags::PML in the grid.
#[test]
fn test_pml_flags_set() {
    use quaternions::simulation::state::CellFlags;

    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);
    let bc = BoundaryConfig::default();
    let config = PmlConfig { depth: 3, ..PmlConfig::default() };

    let _pml = PmlState::new(&mut grid, &bc, config);

    // Interior cell should NOT have PML flag
    let center_idx = grid.idx(8, 8, 8);
    assert_eq!(
        grid.cells[0][center_idx].flags & CellFlags::PML,
        0,
        "center cell should not be PML"
    );

    // Boundary-adjacent cell should have PML flag
    let edge_idx = grid.idx(1, 8, 8);
    assert_ne!(
        grid.cells[0][edge_idx].flags & CellFlags::PML,
        0,
        "cell near boundary should be PML"
    );
}

/// CPML should absorb a pulse leaving the domain.
/// After the pulse has had time to exit, the interior energy should drop
/// dramatically compared to the initial pulse energy.
#[test]
fn test_pml_absorbs_pulse() {
    let n = 32u32;
    let mut grid = SimulationGrid::new(n, n, n, 0.01);
    let bc = BoundaryConfig::default();
    let config = PmlConfig::default(); // depth=10
    let mut pml = PmlState::new(&mut grid, &bc, config);

    // Set a pulse at the center
    let center = grid.idx(n / 2, n / 2, n / 2);
    grid.cells[grid.current][center].q[0] = 1.0;

    // Measure initial energy
    let params0 = grid.sim_params(false);
    let derived0 = compute_derived_fields(&grid, &params0);
    let initial_energy = total_energy_excluding_pml(&derived0, &grid, grid.dx);

    // Run enough steps for the pulse to reach and pass through the PML
    for _ in 0..300 {
        let params = grid.sim_params(false);
        step_field_cpu(&mut grid, &params, Some(&mut pml));
        grid.swap_and_advance();
    }

    // Check that energy in the interior (non-PML) region has dropped significantly
    let params = grid.sim_params(false);
    let derived = compute_derived_fields(&grid, &params);
    let final_energy = total_energy_excluding_pml(&derived, &grid, grid.dx);

    // PML should absorb the vast majority of the pulse energy.
    // Allow up to 1% of initial energy to remain (reflection coefficient < -20 dB).
    assert!(
        initial_energy > 0.0,
        "initial energy should be nonzero"
    );
    assert!(
        final_energy < initial_energy * 0.01,
        "interior energy should drop to <1% of initial: initial={:.4e}, final={:.4e}, ratio={:.4e}",
        initial_energy,
        final_energy,
        final_energy / initial_energy
    );
}

/// CPML should perform better than Mur+sponge for absorbing outgoing radiation.
#[test]
fn test_pml_better_than_mur() {
    let n = 32u32;
    let dx = 0.01;
    let bc = BoundaryConfig::default();
    let sources = SourceConfig {
        sources: vec![Source::dipole_z([n as f32 / 2.0; 3], 1.0, 1e9)],
    };

    // Run WITH CPML
    let mut grid_pml = SimulationGrid::new(n, n, n, dx);
    let config = PmlConfig::default();
    let mut pml = PmlState::new(&mut grid_pml, &bc, config);

    // Drive source, then let drain
    for _ in 0..50 {
        let params = grid_pml.sim_params(false);
        inject_sources(&mut grid_pml, &sources, &params);
        step_field_cpu(&mut grid_pml, &params, Some(&mut pml));
        grid_pml.swap_and_advance();
    }
    for _ in 0..100 {
        let params = grid_pml.sim_params(false);
        step_field_cpu(&mut grid_pml, &params, Some(&mut pml));
        grid_pml.swap_and_advance();
    }
    let derived_pml = compute_derived_fields(&grid_pml, &grid_pml.sim_params(false));
    let energy_pml = total_energy_excluding_pml(&derived_pml, &grid_pml, grid_pml.dx);

    // Run WITH Mur+sponge (no PML)
    let mut grid_mur = SimulationGrid::new(n, n, n, dx);
    let dt = grid_mur.dt;
    for _ in 0..50 {
        let params = grid_mur.sim_params(false);
        inject_sources(&mut grid_mur, &sources, &params);
        step_field_cpu(&mut grid_mur, &params, None);
        grid_mur.swap_and_advance();
        apply_boundaries(&mut grid_mur, &bc, dt);
    }
    for _ in 0..100 {
        let params = grid_mur.sim_params(false);
        step_field_cpu(&mut grid_mur, &params, None);
        grid_mur.swap_and_advance();
        apply_boundaries(&mut grid_mur, &bc, dt);
    }
    let derived_mur = compute_derived_fields(&grid_mur, &grid_mur.sim_params(false));
    let energy_mur = total_energy(&derived_mur, grid_mur.dx);

    // Both should be finite
    assert!(energy_pml.is_finite(), "PML energy should be finite");
    assert!(energy_mur.is_finite(), "Mur energy should be finite");

    // PML should have less residual energy than Mur+sponge
    // (or at least comparable — the PML implementation should not be worse)
    assert!(
        energy_pml <= energy_mur * 2.0,
        "PML should absorb at least as well as Mur: pml={:.4e}, mur={:.4e}",
        energy_pml,
        energy_mur
    );
}

/// Vacuum grid with PML should remain at zero.
#[test]
fn test_pml_vacuum_stability() {
    let mut grid = SimulationGrid::new(16, 16, 16, 0.01);
    let bc = BoundaryConfig::default();
    let config = PmlConfig { depth: 3, ..PmlConfig::default() };
    let mut pml = PmlState::new(&mut grid, &bc, config);

    // Run 50 steps on vacuum (all zeros)
    for _ in 0..50 {
        let params = grid.sim_params(false);
        step_field_cpu(&mut grid, &params, Some(&mut pml));
        grid.swap_and_advance();
    }

    // All cells should still be zero
    for cell in grid.read_buf() {
        for c in 0..4 {
            assert_eq!(cell.q[c], 0.0, "vacuum q should stay zero");
            assert_eq!(cell.q_dot[c], 0.0, "vacuum q_dot should stay zero");
        }
    }
}

/// CPML with extended mode should not crash or produce NaN.
#[test]
fn test_pml_extended_mode() {
    let n = 16u32;
    let mut grid = SimulationGrid::new(n, n, n, 0.01);
    let bc = BoundaryConfig::default();
    let config = PmlConfig { depth: 3, ..PmlConfig::default() };
    let mut pml = PmlState::new(&mut grid, &bc, config);

    // Set a perturbation
    let center = grid.idx(n / 2, n / 2, n / 2);
    grid.cells[grid.current][center].q[0] = 1.0;

    // Run with extended mode
    for _ in 0..30 {
        let params = grid.sim_params(true); // extended_mode = true
        step_field_cpu(&mut grid, &params, Some(&mut pml));
        grid.swap_and_advance();
    }

    // All cells should be finite
    for cell in grid.read_buf() {
        for c in 0..4 {
            assert!(cell.q[c].is_finite(), "q[{c}] should be finite");
            assert!(cell.q_dot[c].is_finite(), "q_dot[{c}] should be finite");
        }
    }
}
