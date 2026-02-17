// Physical parameter controls — stub for future implementation.
//
// PSEUDOCODE: This panel will provide sliders and controls for:
//
// 1. Grid configuration:
//    - Grid size: dropdown (32³, 64³, 128³) — requires buffer re-creation
//    - dx: slider (0.001 .. 0.1 m) — changes domain size at fixed grid count
//    - Show computed domain size = grid_size * dx
//
// 2. Source parameters (wired to SourceConfig resource):
//    - Source type selector (dipole, point charge, bifilar coil, etc.)
//    - Frequency: slider (1 MHz .. 10 GHz)
//    - Amplitude: slider (0 .. max)
//    - Position: 3x slider (x, y, z in grid coordinates)
//    - Dipole axis: dropdown (X, Y, Z)
//
// 3. Boundary conditions (wired to BoundaryConfig resource):
//    - Type per face: dropdown (Absorbing, Conducting, Periodic)
//    - PML depth: slider (2..16 cells)
//    - Damping coefficient: slider
//
// 4. Vacuum model (wired to VacuumConfig resource):
//    - Model selector: Fixed (K=1), Puthoff PV, Custom
//    - Coupling constant α: slider (log scale)
//    - K_background: slider (0.5 .. 2.0)
//
// 5. Visualization (wired to visualization resources):
//    - Field quantity: dropdown (|E|, |B|, S, phi, energy, K)
//    - Color map: dropdown (Viridis, Plasma, Coolwarm, ScalarField)
//    - Range mode: Auto / Manual with min/max sliders
//    - Slice plane position: slider per active slice
//
// Implementation: each section is a collapsible egui::CollapsingHeader
// that reads/writes the relevant Bevy resource via system parameters.
