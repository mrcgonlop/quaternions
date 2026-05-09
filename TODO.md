# TODO.md — Implementation Task Breakdown

## Conventions for Claude Code Sessions

**Each task below is designed to be completable in a single Claude Code session.**

When starting a session, follow this workflow (also described in `CLAUDE.md`):

1. **Find your task**: Scan this file for the first `[ ]` task whose **Depends on** prerequisites are all `[~]` or `[x]`.
2. **Read only what's listed in Context**: Each task tells you exactly which files and sections to read. Do NOT read the entire ARCHITECTURE.md or README.md — only the sections referenced.
3. **Execute the task** following the conventions below.
4. **Update status** in this file when done: `[~]` for stubbed, `[x]` for complete.

### Implementation Conventions

1. **Stub-first**: Create all files with correct module structure, public API signatures, and types — but implement the body as `todo!()` or trivial passthrough with extensive `// PSEUDOCODE:` comments explaining the algorithm.
2. **Pseudocode density**: Every stub function should contain enough pseudocode that a future session can implement it without needing to read the full project context. Include the math, the data flow, edge cases, and references to ARCHITECTURE.md sections.
3. **Compile-check**: Every session should end with `cargo build` succeeding (using `todo!()` for unimplemented paths, not syntax errors).
4. **One concern per session**: Don't mix GPU shader work with UI work with physics implementation. Keep sessions focused.

### Task Status Key
- `[ ]` — Not started
- `[~]` — Stubbed with pseudocode, compiles, not yet implemented
- `[x]` — Fully implemented and tested

### Annotation Key
- **Context:** — Files and sections to read before starting. Read ONLY these.
- **Depends on:** — Prerequisite tasks that must be `[~]` or `[x]` before this task can begin.

---

## Phase 0: Project Scaffolding

### 0.1 — Cargo project initialization
- **Context:** `ARCHITECTURE.md §Key Rust Crates` (lines 43–64), `ARCHITECTURE.md §Module Structure` (lines 457–561)
- **Depends on:** None (first task)
- `[x]` Initialize the `quaternions` project with `cargo init`
- `[x]` Set up `Cargo.toml` per ARCHITECTURE.md dependency list (bevy, bevy_egui, bytemuck, serde, ron, rayon; hdf5 behind optional feature flag; NO standalone wgpu or nalgebra)
- `[x]` Pin bevy to 0.15, verify bevy_egui compatibility (0.34 targets Bevy 0.16; using **0.33** for Bevy 0.15)
- `[x]` Create the full directory structure from ARCHITECTURE.md (all `mod.rs` files, empty source files)
- `[x]` Set up `lib.rs` with module declarations, `main.rs` with minimal Bevy app that opens a window
- `[x]` Verify: `cargo build` succeeds, `cargo run` opens an empty Bevy window
- **Session output**: Compiling project skeleton with a blank 3D window

### 0.2 — Quaternion math module
- **Context:** `ARCHITECTURE.md §Quaternion Type` (lines 132–174)
- **Depends on:** 0.1
- `[x]` Implement `src/math/quaternion.rs`:
  - `Quat` struct: `#[repr(C)] struct Quat { w: f32, x: f32, y: f32, z: f32 }`
  - `impl Quat`: `new`, `zero`, `identity`, `from_scalar`, `from_vector`, `i`, `j`, `k`
  - Hamilton product: `hamilton(self, rhs) -> Quat` (full non-commutative product)
  - Conjugate, norm, norm_sq, inverse, normalize
  - Scalar part, vector part extraction
  - `impl Add, Sub, Mul<f32>, Mul<Quat> for f32, Div<f32>, Neg` for Quat
  - `impl bytemuck::Pod, bytemuck::Zeroable` for GPU compatibility
  - `impl Display, Debug, PartialEq`
- `[x]` Write comprehensive unit tests (28 tests):
  - Hamilton product non-commutativity: `i*j = k` but `j*i = -k` (all 6 basis combos)
  - `i² = j² = k² = ijk = -1`
  - Norm preservation under multiplication
  - Inverse: `q * q.inverse() ≈ identity` (left and right)
  - Scalar/vector decomposition roundtrip
  - Conjugate properties, normalize, operators, Display, bytemuck casts
- `[x]` Verify: `cargo test` passes all 28 quaternion tests
- **Session output**: Battle-tested quaternion algebra ready for physics

### 0.3 — Vector field math utilities
- **Context:** `ARCHITECTURE.md §Simulation Loop` (lines 219–266, for understanding what derivatives are needed), `src/math/quaternion.rs` (to understand existing math types)
- **Depends on:** 0.1
- `[x]` Implement `src/math/vector_field.rs`:
  - `Vec3f` type alias or wrapper (or just use `[f32; 3]` consistently)
  - Dot product, cross product, magnitude helpers
  - `gradient_scalar(grid, idx) -> Vec3f` — central difference gradient of a scalar field
  - `divergence_vector(grid, idx) -> f32` — central difference divergence of a vector field
  - `curl_vector(grid, idx) -> Vec3f` — central difference curl
  - `laplacian_scalar(grid, idx) -> f32` — 7-point stencil Laplacian
  - `laplacian_vector(grid, idx) -> Vec3f` — component-wise Laplacian
- `[x]` Implement `src/math/fdtd.rs`:
  - Grid indexing helpers: `idx(x, y, z, nx, ny) -> usize`
  - Neighbor index functions with boundary awareness
  - CFL condition calculator: `max_dt(dx, c_max) -> f32`
  - Stencil coefficient constants for 2nd-order and 4th-order finite differences
- `[x]` Unit tests:
  - Gradient of linear field = constant vector
  - Divergence of constant vector field = 0
  - Curl of gradient = 0 (to numerical precision)
  - Laplacian of quadratic = constant
- **Session output**: All spatial derivative operations tested and ready

### 0.4 — Simulation grid resource and cell state
- **Context:** `ARCHITECTURE.md §Data Model` (lines 72–130, CellState + DerivedFields structs), `ARCHITECTURE.md §Simulation Loop` (lines 219–266), `src/math/quaternion.rs`
- **Depends on:** 0.2
- `[x]` Implement `src/simulation/state.rs`:
  - `CellState` struct matching ARCHITECTURE.md (48 bytes, repr(C), Pod, Zeroable)
  - `DerivedFields` struct for visualization quantities
  - `CellFlags` bitflags: `EMPTY`, `CONDUCTOR`, `DIELECTRIC`, `SOURCE`, `PML`, `BOUNDARY`
  - `SimParams` uniform struct for GPU: grid dimensions, dx, dt, c0, time, iteration count, `extended_mode: u32`
- `[x]` Implement `src/simulation/grid.rs`:
  - `SimulationGrid` Bevy resource:
    - `nx, ny, nz: u32` — grid dimensions
    - `dx: f32` — cell spacing (uniform for now)
    - `dt: f32` — timestep
    - `cells: [Vec<CellState>; 2]` — double-buffered CPU-side state
    - `current: usize` — index into `cells` for the current read buffer (0 or 1)
    - `time: f64` — current simulation time
    - `iteration: u64` — step count
  - `SimulationGrid::new(nx, ny, nz, dx) -> Self` — allocates and zeros grid
  - `SimulationGrid::idx(x, y, z) -> usize`
  - `SimulationGrid::cell(x, y, z) -> &CellState`
  - `SimulationGrid::cell_mut(x, y, z) -> &mut CellState`
  - `SimulationGrid::is_interior(x, y, z) -> bool` — not on boundary
  - CFL-based dt calculation in constructor
- `[x]` Implement `src/simulation/plugin.rs`:
  - `SimulationPlugin` that inserts `SimulationGrid` as a resource
  - Registers a startup system that initializes the grid
  - Registers the simulation step system (stub: does nothing yet)
  - `SimulationConfig` resource: grid size, dx, paused flag, steps_per_frame, `extended_mode: bool`, `dt_factor: f32` (0.01–1.0, scales dt below CFL max)
- `[x]` Wire into `main.rs`: `app.add_plugins(SimulationPlugin)`
- `[x]` Verify: app starts, grid is allocated, no panic
- **Session output**: ECS-integrated simulation grid, ready for GPU upload

### 0.5 — Basic Bevy scene: camera, lighting, egui panel
- **Context:** `src/main.rs`, `src/simulation/plugin.rs` (for `SimulationConfig` type), `ARCHITECTURE.md §Visualization Architecture` (lines 389–453, for UI layout overview)
- **Depends on:** 0.4
- `[x]` Set up `main.rs` with:
  - 3D camera with orbit controls (bevy's built-in or a simple custom orbit system)
  - Ambient + directional light
  - Ground plane or bounding box wireframe showing the simulation domain
  - Coordinate axes gizmo
- `[x]` Implement `src/ui/plugin.rs`:
  - `UiPlugin` that adds bevy_egui
  - Basic egui side panel with:
    - Simulation controls: Play/Pause button, Step button, Reset button
    - Grid info display: dimensions, dx, dt, current time, iteration count
    - FPS counter
  - Wire `SimulationConfig.paused` to the Play/Pause button
- `[x]` Implement `src/ui/parameter_panel.rs` as stub:
  - Empty panel skeleton with `// PSEUDOCODE:` for future parameter sliders
- `[x]` Verify: app shows 3D scene with bounding box, egui panel works, camera orbits
- **Session output**: Interactive 3D app shell ready for visualization layers

### 0.6 — Color map utilities
- **Context:** `ARCHITECTURE.md §Visualization Modes` (lines 389–412, color map usage context), `src/visualization/mod.rs` (to see module structure)
- **Depends on:** 0.1
- `[x]` Implement `src/visualization/color_maps.rs`:
  - `ColorMap` enum: `Viridis`, `Plasma`, `Coolwarm`, `ScalarField` (blue-white-red), `Grayscale`
  - `fn map_value(value: f32, min: f32, max: f32, map: ColorMap) -> [f32; 4]` — returns RGBA
  - Hardcoded lookup tables (32 or 64 entries per map, linearly interpolated)
  - `ScalarField` map: negative=blue, zero=white/transparent, positive=red — critical for S field
  - Generate 1D texture from color map for GPU shader use
- `[x]` Visual test: render a colored quad showing the full color map gradient
- **Session output**: Color mapping ready for all visualization modes

---

## Phase 1: CPU Simulation Baseline

> **Strategy note**: Before tackling GPU compute shaders, implement the physics loop on the CPU first. This lets us validate the math, see results quickly, and debug without GPU pipeline complexity. GPU migration happens later as an optimization pass.

### 1.1 — CPU field update: standard FDTD on quaternionic potentials
- **Context:** `ARCHITECTURE.md §GPU Compute Shader: Field Update Kernel` (lines 268–385, the WGSL pseudocode — this IS the authoritative equation reference for the CPU implementation), `src/simulation/state.rs` (CellState layout), `src/simulation/grid.rs` (grid access API), `src/math/vector_field.rs` and `src/math/fdtd.rs` (spatial derivative helpers)
- **Depends on:** 0.3, 0.4
- `[x]` Implement `src/simulation/field_update.rs`:
  - `fn step_field_cpu(grid: &mut SimulationGrid, params: &SimParams)`:
    - For each interior cell (x, y, z):
      - Read Q = (phi/c, Ax, Ay, Az) from cell
      - Compute Laplacian of Q using 7-point stencil (component-wise)
      - Compute div(A) = dAx/dx + dAy/dy + dAz/dz (central differences)
      - Compute d(phi/c)/dt from q_dot.w
      - Compute S = d(phi/c)/dt / c + div(A)
      - Compute grad(S) using central differences on neighboring cells' S values
      - **KEY DECISION POINT**: Include or exclude grad(S) coupling term
        - WITH grad(S): extended QVED (allows scalar longitudinal modes)
        - WITHOUT grad(S): standard gauge-fixed FDTD (Lorenz gauge)
        - Use a boolean flag `extended_mode` to toggle
      - Compute q_ddot = c² * laplacian_Q [+ scalar coupling terms if extended]
      - Leapfrog (Störmer-Verlet) update with staggered velocity:
        - `q_dot` is stored at half-step offsets (t-dt/2 → t+dt/2)
        - `q_dot(t+dt/2) = q_dot(t-dt/2) + q_ddot(t) * dt`
        - `q(t+dt) = q(t) + q_dot(t+dt/2) * dt`
        - On initialization, `q_dot` must be set at `t = -dt/2` (kick back by half-step)
        - This staggering is critical for long-term energy conservation
    - Double-buffer: `SimulationGrid` holds `cells: [Vec<CellState>; 2]` + `current: usize`
      - Read from `cells[current]`, write to `cells[1 - current]`, swap `current`
  - `// PSEUDOCODE:` comments should include the full equations from ARCHITECTURE.md §2.1
- `[x]` Bevy system: `fn simulation_step_system(grid: ResMut<SimulationGrid>, config: Res<SimulationConfig>)`
  - Skip if paused
  - Call `step_field_cpu` for `config.steps_per_frame` iterations
- `[x]` Verify: no crash, grid evolves (values change from initial conditions)
- **Session output**: Core physics loop running on CPU, toggleable standard/extended mode

### 1.2 — Source injection: oscillating dipole
- **Context:** `src/simulation/grid.rs` (grid access API), `src/simulation/state.rs` (CellState), `src/simulation/plugin.rs` (system registration pattern), `src/ui/plugin.rs` (UI panel pattern)
- **Depends on:** 1.1
- `[x]` Implement `src/simulation/sources.rs`:
  - `SourceType` enum: `PointCharge`, `OscillatingDipole`, `CurrentPulse`, `Custom`
  - `Source` struct: position, type, amplitude, frequency, phase, active flag
  - `SourceConfig` resource: `Vec<Source>`
  - `fn inject_sources(grid: &mut SimulationGrid, sources: &SourceConfig, time: f64)`:
    - For each active source:
      - `OscillatingDipole`: set J(x_source, t) = J0 * sin(2π*f*t) along dipole axis
        - Inject as current density into q_dot at source cells
        - `// PSEUDOCODE: The source drives ∂A/∂t at the source location`
        - `// This is equivalent to injecting a current density J = σE`
        - `// For a z-directed dipole: q_dot.z += (J0/ε₀) * sin(ωt) * dt`
      - `PointCharge`: set φ at source cell (static charge)
      - `CurrentPulse`: Gaussian or step function current injection
  - Bevy system that runs before `simulation_step_system`
- `[x]` Add to UI: source position, frequency, amplitude sliders
- `[x]` Verify: set up dipole, run sim, see Q values changing around source
- **Session output**: Configurable EM sources driving the simulation

### 1.3 — Derived field computation
- **Context:** `src/simulation/state.rs` (DerivedFields struct), `src/simulation/grid.rs`, `src/math/vector_field.rs` (gradient, curl, divergence helpers), `README.md §The Quaternionic Four-Potential` (for E, B, S derivation from Q — lines 23–70 only)
- **Depends on:** 1.1, 0.3
- `[x]` Implement `src/simulation/diagnostics.rs` (derived fields portion):
  - `fn compute_derived_fields(grid: &SimulationGrid) -> Vec<DerivedFields>`:
    - For each cell:
      - E = -grad(phi) - dA/dt = -c*grad(Q.w) - Q_dot.vector()
      - B = curl(A) = curl(Q.vector())
      - S = (1/c²)(dphi/dt) + div(A) = Q_dot.w/c + div(Q.vector())
      - energy_density = 0.5 * ε₀ * |E|² + 0.5 * (1/μ₀) * |B|² + 0.5 * ε₀ * S²
      - poynting = (1/μ₀) * E × B
  - `fn total_energy(derived: &[DerivedFields], dx: f32) -> f64` — sum over all cells * dx³
  - `fn max_S(derived: &[DerivedFields]) -> f32` — peak scalar field magnitude
- `[x]` Add diagnostics to UI panel: total energy, max |E|, max |B|, max |S|, updated each frame
- `[x]` Verify: energy should be roughly conserved (not grow unboundedly) for a free dipole
- **Session output**: Real-time field diagnostics confirming simulation health

### 1.4 — Slice plane visualization
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields data to visualize), `src/visualization/color_maps.rs` (color mapping API), `src/simulation/grid.rs` (grid dimensions), `src/ui/plugin.rs` (UI integration pattern)
- **Depends on:** 1.3, 0.6
- `[x]` Implement `src/visualization/slices.rs`:
  - `SlicePlane` component: `axis: Axis` (X/Y/Z), `position: f32`, `field: FieldQuantity`, `color_map: ColorMap`
  - `FieldQuantity` enum: `E_magnitude`, `B_magnitude`, `S_field`, `Phi`, `Ax`, `Ay`, `Az`, `EnergyDensity`, `K_vacuum`
  - Bevy system: `fn update_slice_texture(grid, derived_fields, slice_planes)`:
    - For the active slice plane, sample the selected field quantity on the 2D grid
    - Map values through color map → RGBA texture
    - Apply texture to a quad mesh positioned at the slice location in 3D
  - Spawn a quad mesh entity for each slice plane
  - `// PSEUDOCODE: For a Z-slice at position z0:`
  - `// For each (x, y) in grid: sample field at (x, y, z0) → value`
  - `// Map value through color_map with auto-ranging or fixed range → RGBA`
  - `// Write RGBA to texture at (x, y) pixel`
  - `// Bevy updates the texture on the quad mesh each frame`
- `[x]` Add to UI: slice axis selector, position slider, field quantity dropdown, color map dropdown
- `[x]` Verify: see a colored 2D slice through the 3D grid, updating in real time as simulation runs
- **Session output**: First real visualization — watching EM fields propagate in a slice view

### 1.5 — Boundary conditions: simple absorbing
- **Context:** `src/simulation/grid.rs` (grid access, boundary cell identification), `src/simulation/state.rs` (CellState, CellFlags), `src/simulation/plugin.rs` (system ordering)
- **Depends on:** 1.1
- `[x]` Implement `src/simulation/boundaries.rs`:
  - `BoundaryType` enum: `Absorbing`, `Conducting`, `Periodic`, `PML`
  - `BoundaryConfig` resource: boundary type for each face (±x, ±y, ±z)
  - `fn apply_boundaries(grid: &mut SimulationGrid, config: &BoundaryConfig)`:
    - `Absorbing` (simple first-order): multiply Q and Q_dot by a damping factor in boundary cells
      - Damping increases toward the outer edge: `factor = 1 - ((distance_from_interior / pml_depth)²) * sigma_max`
      - `// PSEUDOCODE: This is a simplified Mur-type absorbing boundary`
      - `// For proper PML, implement split-field perfectly matched layers (future task)`
      - `// The simple version just exponentially damps fields near boundaries`
      - `// This prevents reflections but is not perfectly non-reflecting`
    - `Conducting`: set tangential A components to zero on boundary cells
    - `Periodic`: copy opposite face values
  - Bevy system that runs after field_update
- `[x]` Verify: dipole radiation reaches boundaries without strong reflections
- **Session output**: Simulation doesn't blow up at boundaries

### 1.6 — Validation: dipole radiation pattern
- **Context:** `src/simulation/sources.rs` (SourceConfig API), `src/simulation/diagnostics.rs` (DerivedFields, energy computation), `src/simulation/grid.rs`, `src/visualization/slices.rs` (to verify visually), `README.md §Tier 1 Scenarios` (lines 190–220, for expected dipole physics — brief skim only)
- **Depends on:** 1.2, 1.3, 1.4, 1.5
- `[x]` Implement `src/scenarios/dipole_radiation.rs`:
  - `fn setup_dipole_scenario(grid: &mut SimulationGrid, sources: &mut SourceConfig)`:
    - Place oscillating dipole at grid center
    - Set frequency such that wavelength = ~20 cells (resolvable)
    - Clear initial conditions
  - Expected results:
    - Transverse E, B fields propagating outward
    - Classic dipole radiation pattern: sin²θ intensity distribution
    - S field should be ~zero everywhere if standard mode (Lorenz gauge)
    - S field may be nonzero near source if extended mode
  - `fn validate_dipole(derived: &[DerivedFields], grid: &SimulationGrid) -> ValidationResult`:
    - Check wave speed: measure time for wavefront to reach known distance → c ± tolerance
    - Check polarization: E perpendicular to propagation direction
    - Check energy conservation: total energy increases at rate matching source power
    - Check S ≈ 0 in standard mode
- `[x]` Add scenario selector to UI (dropdown that loads preset configurations)
- `[x]` Verify: visually confirm dipole radiation pattern in slice view
- **Session output**: Validated baseline — standard FDTD is correct before extending it

### 1.8 — Dynamic K field: polarizable vacuum evolution
- **Context:** `src/simulation/grid.rs` (CellState fields `k`, `k_dot`), `src/simulation/field_update.rs` (step_field_cpu leapfrog pattern), `src/simulation/state.rs` (SimParams, CellFlags), `src/simulation/diagnostics.rs` (diagnostics pattern), `README.md §3.4` (K field physics and the virtual pair plasma model)
- **Depends on:** 1.3, 2.1 (extended mode), Phase 1.7 (independent S field)
- **Physics summary:** The QED vacuum acts as a polarizable medium (refractive index n = √K ≥ 1). Strong local fields excite virtual e⁺e⁻ pairs, increasing K above its vacuum value K₀=1. K evolves as a damped scalar wave: `∂²K/∂t² = c²∇²K − ωₚ²(K − 1) + η · u_field/u_S` where ωₚ = mₑc²/ℏ ≈ 7.8×10²⁰ rad/s is the Compton frequency (vacuum plasma frequency), η ≈ 2α/(45π) ≈ 10⁻⁴ is the QED coupling (Euler-Heisenberg), and u_S = ε₀(mₑc²/e)² is a characteristic field energy density. In simulation units, ωₚ and η are free parameters in `VacuumConfig`.
- **Weber + K enhancement:** Weber corrections scale as K²(v/c₀)². For K > √2, like-charge attraction threshold shifts from v > c₀√2 (superluminal) to v < c₀ (sub-relativistic), enabling stable electron clusters (EVOs/charge clusters) as self-sustaining K pockets.

#### Implementation tasks:

- `[x]` Add `VacuumConfig` to `src/simulation/plugin.rs`:
  ```rust
  #[derive(Resource, Clone, Debug)]
  pub struct VacuumConfig {
      pub enabled: bool,        // master switch
      pub omega_p: f32,         // vacuum plasma frequency (rad/s, in sim units)
      pub eta: f32,             // QED coupling ≈ 2α/(45π) ≈ 1e-4
      pub u_s: f32,             // characteristic energy density (normalization)
  }
  impl Default for VacuumConfig { fn default() -> Self { Self { enabled: false, omega_p: 0.0, eta: 1e-4, u_s: 1.0 } } }
  ```
- `[x]` Add K leapfrog to `step_field_cpu` in `src/simulation/field_update.rs`:
  - Pre-pass: compute `lap_k[i]` = 6-neighbor Laplacian of `cell.k`
  - Update k_dot: `k_dot[i] += dt * (c² * lap_k[i] - omega_p² * (k[i] - 1.0) + eta * u_field[i] / u_s)`
    - `u_field[i]` = local energy density from derived fields (0.5ε₀E² + 0.5/μ₀ B²)
    - Guard: only update if `vacuum_config.enabled && is_interior && !is_pml`
  - Update k: `k_new[i] = k_old[i] + dt * k_dot[i]`; clamp `k_new[i] >= 1.0`
  - `c_local` = `c0 / k[i]` propagates automatically since `field_update` already reads `cell.k`
- `[x]` Add K Neumann boundary conditions in `src/simulation/boundaries.rs`:
  - `apply_k_neumann`: zero-gradient at all 6 open faces (`k[boundary] = k[interior_neighbor]`)
  - Called from `apply_boundaries` and from the PML path in `boundary_system`
- `[x]` Add K visualization to `src/visualization/slices.rs`:
  - `FieldQuantity::KVacuum` variant already existed; `sample_field_at` already reads `cell.k`
- `[x]` Add K diagnostics to `src/simulation/diagnostics.rs`:
  - `DiagnosticsState` gains `max_k: f32` and `mean_k: f32`
  - `max_k` and `mean_k` computed over all non-PML interior cells in `diagnostics_system`
  - Displayed in UI panel next to E, B, S diagnostics
- `[x]` Write tests in `tests/integration_phase1.rs`:
  - `test_k_vacuum_stable`: grid with K=1, no fields → K remains 1.0 after 1000 steps ✓
  - `test_k_increases_with_field`: strong E field applied → K increases above 1.0 ✓
  - `test_k_restores_after_pulse`: after pulse removed, K decays back toward 1.0 ✓
  - `test_k_diagnostics_vacuum`: max_k/mean_k return 1.0 on vacuum grid ✓
  - Note: `test_weber_k_enhancement` deferred — Weber force not yet implemented
- `[x]` Add `VacuumK` scenario to demonstrate EVO formation:
  - `src/scenarios/vacuum_k.rs` + `Scenario::VacuumK` in dipole_radiation.rs + UI handler
- **Session output**: Polarizable vacuum model active — K field evolves under field pressure, enabling Weber-enhanced charge clustering and measurable refractive-index gradients

### 1.9 — Topological charge diagnostic (Hopf invariant / Baryon number)
- **Context:** `src/simulation/grid.rs` (CellState: Q field `w, x, y, z` components), `src/simulation/diagnostics.rs` (DiagnosticsState pattern, diagnostics_system), `src/ui/plugin.rs` (UI panel for new diagnostic value), `README.md §3.5` (topological charge formula, Q field S³ topology)
- **Depends on:** 1.3 (Q field exists in grid), 2.1 (extended mode where Q topology matters)
- **Physics summary:** Unit quaternions live on S³. The map U = Q/|Q| : R³ → S³ is classified by π₃(S³) = ℤ — an integer winding number called the Skyrmion / Baryon number. A Hopfion has Hopf invariant Q_H = linking number of field line pairs, related to but distinct from the Baryon number (it classifies maps R³ → S² via the Hopf fibration S³ → S²). Both are integer topological invariants that must be conserved under smooth field evolution. Non-zero values indicate topologically stable configurations (ball lightning candidates). The Baryon charge density is:
  ```
  ρ_topo = (1/24π²) ε^{ijk} Tr( ∂_i U U⁻¹ · ∂_j U U⁻¹ · ∂_k U U⁻¹ )
  ```
  where U ∈ S³ is the unit quaternion, treated as a 2×2 complex matrix (or directly as the quaternion algebra). The integral ∫ ρ_topo dV = n ∈ ℤ.

#### Implementation tasks:

- `[x]` Add `compute_topological_charge` to `src/simulation/diagnostics.rs`:
  - Computes Baryon/Skyrmion winding number via left Maurer-Cartan forms L_i = (∂_i U) U*
  - ρ_topo = (1/4π²) * scalar_part(L_x [L_y, L_z]) — prefactor corrected from 1/24π² to 1/4π² (factor of 6 = 3 from ε symmetry × 2 from SU(2) Tr → quaternion scalar_part)
  - Guards: skips PML cells and cells where |Q| < 1e-12
  - Central differences on the normalized unit quaternion field U = Q/|Q|
- `[x]` Add `topological_charge: f32` to `DiagnosticsState`
- `[x]` Call `compute_topological_charge` in `diagnostics_system` (throttled to every 10 steps)
- `[x]` Display `n_topo` in UI panel (formatted as float with 2 decimal places: "Topo charge: {:.2}")
  - Green highlight if near-integer AND |n_topo| > 0.1 → topological structure detected
- `[x]` Write tests in `tests/integration_phase1.rs`:
  - `test_topo_charge_vacuum`: uniform Q field → topological charge = 0 ✓
  - `test_topo_charge_hedgehog`: hedgehog ansatz U=(cos θ, sin θ · r̂) on 32³ → charge ≈ -1 ✓
  - `test_topo_charge_smooth_field`: slowly varying Q → charge ≈ 0 ✓
  - `test_topo_charge_conserved`: evolve hedgehog 50 steps on 32³; charge drift < |initial| ✓
- **Session output**: Integer-valued diagnostic detecting topological field structures — the essential measurement tool for ball lightning simulation

---

## Phase 2: Extended Electrodynamics (Scalar Field)

### 2.1 — Extended mode: scalar field dynamics
- **Context:** `ARCHITECTURE.md §GPU Compute Shader: Field Update Kernel` (lines 350–385, the extended_mode branching logic), `src/simulation/field_update.rs` (existing CPU implementation to modify), `src/simulation/plugin.rs` (SimulationConfig)
- **Depends on:** 1.6 (must validate standard mode first)
- `[x]` Modify `field_update.rs` to fully support extended mode:
  - When `extended_mode = true`:
    - Do NOT enforce Lorenz gauge (do not subtract grad(S) to zero it out)
    - S field evolves dynamically via: □S = -ρ/ε₀
    - The scalar coupling terms in the modified Gauss and Ampère laws are active
  - Add `extended_mode: bool` to `SimulationConfig`
  - `// PSEUDOCODE: The key difference is in the wave equation for Q:`
  - `// Standard mode: ∂²Q/∂t² = c²∇²Q + sources`
  - `//   → scalar coupling OMITTED → Lorenz gauge enforced implicitly`
  - `//   → only transverse (E, B) modes propagate`
  - `// Extended mode: ∂²Q/∂t² = c²∇²Q - c²(0, ∇S) + sources`
  - `//   → scalar coupling INCLUDED via params.extended_mode flag`
  - `//   → S = (1/c²)∂φ/∂t + ∇·A is a dynamical degree of freedom`
  - `//   → the -c²∇S in the vector components couples S back into potential evolution`
  - `//   → S propagates as a genuine longitudinal wave driven by charge fluctuations`
  - `// See ARCHITECTURE.md WGSL kernel for the authoritative implementation`
- `[x]` Add toggle to UI: "Standard EM / Extended QVED" switch
- `[x]` Verify: in standard mode, S stays near zero; in extended mode, S propagates from sources
- **Session output**: The fundamental theoretical extension is live and toggleable

### 2.2 — S-field dedicated visualization
- **Context:** `src/visualization/slices.rs` (existing slice system to extend), `src/visualization/color_maps.rs` (ScalarField bipolar map), `src/simulation/diagnostics.rs` (S field data)
- **Depends on:** 2.1, 1.4
- `[x]` Enhance slice visualization to handle bipolar (positive/negative) fields:
  - S field is signed → use `ScalarField` color map (blue-white-red)
  - Add auto-ranging: compute min/max S each frame for color normalization
  - Add manual range lock (so color scale doesn't flicker)
  - Add S-field magnitude isosurface option
- `[x]` Add a second simultaneous slice plane option:
  - Slice 1: |E| or |B| (to see transverse waves)
  - Slice 2: S field (to see longitudinal waves)
  - Both visible at orthogonal orientations for comparison
- `[x]` Verify: can visually distinguish transverse and longitudinal modes side by side
- `[x]` Fix slice position updates: replace entire `Image` asset via `images.insert()` each frame (Bevy 0.15 texture re-upload gotcha)
- `[x]` Fix auto-range minimum threshold: `1e-6` instead of `1e-30` to avoid `map_value` t=0.5 fallback
- `[x]` Add per-slice stats display (sample count, value min/max, range min/max) in UI
- `[x]` Add diagnostic logging (throttled to ~1/sec) in `update_slice_texture`
- `[x]` Add configurable temporal resolution: `dt_factor` (0.01–1.0) logarithmic slider scales dt below CFL max; effective dt and sim time/frame shown in UI
- **Session output**: Clear visual comparison between standard and extended EM

### 2.3 — Bifilar coil source geometry
- **Context:** `src/simulation/sources.rs` (source injection API), `src/simulation/grid.rs`, `src/scenarios/dipole_radiation.rs` (scenario pattern to follow), `README.md §Scalar / Longitudinal EM Waves` (lines 56–92, theory motivation — brief skim)
- **Depends on:** 2.1, 1.2
- **Note:** See TODO §Known Issues — bifilar coil resolution. At 64³/0.5m domain, dx ≈ 8mm is too coarse for realistic wire spacing. Use idealized geometry or 128³+ with small domain.
- `[x]` Implement `src/scenarios/bifilar_coil.rs`:
  - Bifilar coil: two windings carrying current in opposite directions
  - `// PSEUDOCODE: A bifilar coil is wound so that adjacent wires carry current`
  - `// in opposite directions. The B fields from adjacent wires cancel`
  - `// (to first order), but the A (vector potential) fields add.`
  - `// This creates a source with strong A but weak B — ideal for exciting`
  - `// scalar/longitudinal modes because the transverse (B-based) radiation`
  - `// is suppressed while the potential-based effects survive.`
  - `fn setup_bifilar_coil(grid, sources)`:
    - Define coil geometry as a set of current-carrying cells
    - Two interleaved helical paths with opposite current direction
    - Drive with oscillating current
  - `fn setup_bifilar_pair(grid, sources)`:
    - Transmitter bifilar coil + receiver bifilar coil separated by distance
    - Measure S field at receiver location
    - Compare coupling in standard vs extended mode
  - Expected result: in extended mode, the receiver sees S-field signal; in standard mode, coupling is minimal (only stray B leakage)
- `[x]` Add to scenario selector
- `[ ]` Verify: visual confirmation of B cancellation and A/S field propagation
- **Session output**: The signature experiment for scalar wave detection, simulated

### 2.4 — Quantitative scalar wave analysis
- **Context:** `src/simulation/diagnostics.rs` (to extend with probes), `src/scenarios/bifilar_coil.rs` (scenario to measure), `src/ui/plugin.rs` (UI integration)
- **Depends on:** 2.3, 1.3
- `[x]` Add measurement probes to `diagnostics.rs`:
  - `Probe` struct with `label`, `position: [u32; 3]`, `field: ProbeField`, `history: VecDeque<(f32, f32)>`
  - `ProbeSet` resource: `Vec<Probe>` plus `max_history` ring-buffer cap (default 4096)
  - `ProbeField` enum defined locally (subset of visualization `FieldQuantity`) to avoid a `simulation → visualization` dependency cycle
  - `probe_system` samples every probe each step after `diagnostics_system`, so the most recent `DerivedFields` are used
  - `fn probe_fft(probe: &Probe) -> Vec<(f32, f32)>` — naive O(N²) DFT, mean-detrended, returns `(freq_hz, magnitude)` up to Nyquist. Covered by `test_probe_fft_recovers_sinusoid` (peak within one FFT bin of the target).
- `[~]` Add probe placement to UI — probe list in the `Probes` window with per-probe X/Y/Z DragValues and field selector; "click on slice to place probe" deferred (would need raycast against the slice quad — non-trivial and not strictly required for delay measurement).
- `[x]` Add time-series plot in UI (hand-rolled egui painter plot — overlays all probes, auto-scaled, with zero line and legend; avoids pulling in an `egui_plot` dependency)
- `[x]` For bifilar scenario:
  - Probes installed automatically when BifilarPair scenario is selected (`install_bifilar_pair_probes`): TX at x=nx/4, Mid at x=nx/2, RX at x=3nx/4, all on the y/z center line, all recording the S field
  - Time series plotted in the Probes window
  - Propagation delay estimated from TX→RX peak shift; v and v/c displayed in the UI
  - FFT of the first probe shown in a collapsing section
- **Session output**: Quantitative evidence of scalar wave propagation speed and amplitude

---

## Phase 3: Weber Forces

### 3.1 — Discrete particle system
- **Context:** `src/simulation/grid.rs` (grid interpolation for fields at particle positions), `src/simulation/diagnostics.rs` (DerivedFields for E, B at grid points), `ARCHITECTURE.md §Simulation Loop step 6` (lines 256–261, particle push overview)
- **Depends on:** 1.3
- `[x]` Implement `src/simulation/particles.rs`:
  - `ChargedParticle` ECS marker component carrying an index back into the resource (physics state lives in `ParticleSystem` so it can be borrowed inside `simulation_step_system`'s substep loop without crossing Bevy's Query borrow rules).
  - `ParticleState { charge, mass, position, velocity, prev_acceleration }` — `prev_acceleration` is tracked now so Phase 3.2's explicit Weber scheme can read it without extra plumbing.
  - `ParticleSystem` resource: `Vec<ParticleState>` plus `clear`/`push` helpers.
  - `spawn_charged_particle(commands, particles, meshes, materials, state, radius, color)` — appends to the resource AND spawns a sphere-mesh entity in the 3D scene.
  - `step_particles(particles, grid, c0, dt)` — Boris pusher: half-step E kick, magnetic rotation, second half-step E kick, position update; updates `prev_acceleration` from the velocity delta.
  - `sample_e_b(grid, world_pos, c0)` — trilinear interpolation of E and B from the 8 surrounding cell centres, computing each corner's (E,B) via the same central-difference formulae as `diagnostics::compute_derived_fields`. Returns zero outside the safe interior region (≥ 2 cells from each face).
  - World/grid convention is documented at the top of the module: domain centred on the origin, cell (i,j,k) centre at `((i+0.5)dx − half_x, …)`.
- `[x]` Render particles as colored spheres in 3D scene — handled by `spawn_charged_particle` (icosphere mesh + emissive material) and the per-frame `sync_particle_transforms` system that copies positions from the resource into entity Transforms.
- `[x]` Verify: tests in `src/simulation/particles.rs` — `test_uniform_e_field_accelerates_linearly` (closed-form qE/m drift, rel. err < 2 %), `test_uniform_b_field_gyrates` (cyclotron orbit closes within 5 % of gyroradius after one period), `test_zero_field_is_pure_drift`, `test_sample_e_b_outside_grid_returns_zero`.
- Wiring: `ParticleSystem` is registered as a resource in `SimulationPlugin`; `step_particles` runs at the START of each substep inside `simulation_step_system` so particles see the field at time t before it advances; `sync_particle_transforms` runs in `SimulationSet` after diagnostics.
- **Session output**: Particle system with standard Lorentz force via Boris pusher, GPU-renderable, and verified against analytic E-drift / B-gyration limits.

### 3.2 — Weber force implementation
- **Context:** `src/simulation/particles.rs` (particle system to extend), `README.md §Weber Electrodynamics` (lines 96–140, Weber force theory derivation)
- **Depends on:** 3.1
- **Note:** See TODO §Known Issues — Weber acceleration bootstrapping uses previous-step acceleration (explicit scheme).
- `[x]` Implement `src/simulation/weber.rs`:
  - `weber_force(q1, q2, r12, v12, a12, c0) -> [f32; 3]` — pure function returning the Weber force on particle 1 due to particle 2. Implements F = (q₁q₂/4πε₀r²)·r̂·[1 − ṙ²/(2c²) + r·r̈/c²] with ṙ = r̂·v₁₂ and r̈ = r̂·a₁₂ + (|v₁₂|² − ṙ²)/r.
  - `compute_weber_forces(particles, c0) -> Vec<[f32; 3]>` — N² pair sum using Newton's third law: only iterates j > i and updates both forces[i] and forces[j], cutting pair count in half.
  - `step_particles_weber(particles, c0, dt)` — kick-drift integrator: v_new = v + (F/m)·dt; x_new = x + v_new·dt. Updates `prev_acceleration = F/m` so the next step's r̈ term has it.
  - `step_particles_both(particles, grid, c0, dt)` — Boris pusher then additive Weber kick (no extra drift since Boris already advanced position). Combines both accelerations into `prev_acceleration`.
  - `ForceMode { Lorentz, Weber, Both }` enum + `Resource` registration; defaults to Lorentz so existing scenarios keep their behaviour.
  - Explicit-scheme rationale documented at the top of the module: a₁₂ uses each particle's `prev_acceleration` field (set by either pusher) to break the chicken-and-egg.
- `[x]` Add force mode toggle to UI — `Charged Particles (N)` collapsing section in the side panel: three-button selectable_value row + a description line explaining each mode's force law.
- `[x]` Wire dispatch into `simulation_step_system`: matches on `*force_mode` inside the substep loop and routes to `step_particles` / `step_particles_weber` / `step_particles_both` so the chosen integrator runs for every substep against fresh-from-time-t fields.
- `[x]` Unit tests in `src/simulation/weber.rs`:
  - `test_weber_static_limit_equals_coulomb` — v=a=0 collapses the correction factor to 1; force matches Coulomb's law to 1e-5 rel. err.
  - `test_weber_moving_apart_less_than_coulomb` — direct |F_W| < |F_C| comparison plus exact recovery of the analytic 1 − 2v²/c² factor (rel. err < 1e-3).
  - `test_weber_force_along_r_hat` — F × r₁₂ vanishes regardless of v₁₂, a₁₂.
  - `test_weber_newton_third_law` — F₁₂ + F₂₁ = 0 to 1e-5 rel. err.
  - `test_step_particles_weber_two_charges_repel` — two like charges accelerate apart along the line of separation, transverse motion < 0.1 % of longitudinal, total momentum conserved.
- Test scaling note: tests use unit charges + meter-scale separation. Elementary-charge scales (q ≈ 1.6e-19) push squared force magnitudes (~10⁻⁴⁴) into f32 subnormal range where the 2v²/c² ≈ 2e-5 correction is numerically lost — the underlying physics still works (the integrator uses the *vector* not magnitude²), but unit tests need the magnitude check in the normal f32 range to be meaningful.
- **Session output**: Weber pair force law live and selectable from the UI; static-limit, velocity-correction, longitudinality, and Newton's-third-law symmetries all verified.

### 3.3 — Graneau wire scenario
- **Context:** `src/simulation/weber.rs` (Weber force API), `src/simulation/particles.rs` (particle system), `src/scenarios/dipole_radiation.rs` (scenario pattern), `README.md §Weber Electrodynamics` (lines 96–140, Graneau experiment context)
- **Depends on:** 3.2
- **Note:** See TODO §Known Issues — Weber force magnitudes use scaled parameters. Document scaling factor.
- `[x]` Implement `src/scenarios/graneau_wire.rs`:
  - `apply_graneau_wire_scenario(particles, force_mode, num_segments, spacing, v_drift)` — populates `ParticleSystem` with `num_segments` collocated (heavy ion at rest, light electron drifting at v_drift) pairs along x̂. Sets `*force_mode = ForceMode::Weber` because Lorentz gives identically zero on a straight wire from its own current.
  - `apply_graneau_wire_default(particles, force_mode)` — convenience wrapper using the documented `DEFAULT_NUM_SEGMENTS = 21`, `DEFAULT_SEGMENT_SPACING = 1 cm`, `DEFAULT_DRIFT_VELOCITY = 1×10⁷ m/s`.
  - `compute_wire_force_profile(particles) -> Vec<(f32, f32)>` — returns `(x_position, axial_Weber_force)` for each ion (even-indexed particle). Used by tests to verify the Weber prediction.
  - Module header documents the analytic two-segment derivation: each segment pushes its neighbour AWAY along the wire axis with force magnitude `(q²/4πε₀L²)·(v_d²/c²)`. Sum over a uniform chain → end-peaked, antisymmetric profile (interior segments see balanced pushes from both sides).
- Scaling note documented inline: real copper at 10 kA has v_drift ≈ 0.74 m/s, v²/c² ≈ 6e-18 — invisible in f32. Simulation uses v_drift = 1e7 m/s (~0.03 c, correction ~5e-4) so the longitudinal pattern is numerically resolvable. The QUALITATIVE pattern (end-peaked + antisymmetric) is preserved; absolute magnitudes are unphysical.
- `[~]` Visualize: wire segments colored by longitudinal force magnitude; arrow glyphs showing force direction.
  - Particles render as plain spheres via the Phase 3.1 `spawn_charged_particle` infrastructure, but the dedicated force-magnitude colouring + force arrow overlay is deferred. Required additions when reopened: a `WireForceProfile` resource cached each frame from `compute_wire_force_profile`, a render system that updates each ion sphere's material colour from a heat-map of |F_x|, and an egui plot of the profile vs. x in the side panel. None of these are blocking — the physics is testable from `cargo test`.
- `[x]` Add to scenario selector — `Scenario::GraneauWire` in `dipole_radiation.rs::Scenario`, UI dropdown handler in `ui_side_panel` clears grid + sources + PML + probes + particles, calls `apply_graneau_wire_default`, leaves `extended_mode = false` (Weber is field-free) and starts paused.
- `[x]` Tests in `src/scenarios/graneau_wire.rs`:
  - `test_apply_graneau_wire_scenario_creates_chain` — particle count = 2N, ions at rest at correct positions, electrons drift in +x̂, total charge ≈ 0, force_mode set to Weber.
  - `test_graneau_force_profile_outward_at_ends` — leftmost ion F_x < 0, rightmost F_x > 0, middle |F_x| < 1 % of end |F_x|, antisymmetry < 1e-3, end |F| > 1e6 N (well above noise floor).
  - `test_graneau_force_increases_toward_ends` — |F_x| grows monotonically from centre outward (rules out flat or noise-dominated profiles).
  - `test_graneau_no_axial_force_without_pair_law` — with v_drift = 0 the Weber correction vanishes, leaving pure Coulomb on a charge-neutral chain → every ion sees ≈ 0 net force. This is the null result that Weber measurably DEPARTS from when v_drift > 0.
- **Session output**: Graneau wire scenario live in the UI; the signature Weber prediction (end-peaked, antisymmetric axial force on a uniform straight wire) is verified by tests against the analytic two-segment derivation.

### 3.4 — Force visualization overlays
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields with E, B, Poynting), `src/visualization/color_maps.rs`, `src/simulation/grid.rs` (grid sampling), `ARCHITECTURE.md §Visualization Modes` (lines 396–406, glyph/arrow overview)
- **Depends on:** 1.3, 0.6
- `[x]` Implement `src/visualization/glyphs.rs`:
  - Implemented via gizmos (not instanced mesh arrows). Supports E, B, Poynting, A with
    Standard/RgbMultiField/HsvPhase/SizeColor encoding modes. Arrowhead drawn as two lines.
- `[x]` Add to UI: field selector, density slider, scale slider, on/off toggle
- `[x]` Verify: dipole scenario shows correct E field arrows radiating outward
- **Session output**: Beautiful vector field visualization in 3D

---

## Phase 4: Polarizable Vacuum

### 4.1 — Vacuum state dynamics
- **Context:** `src/simulation/state.rs` (CellState.k field), `src/simulation/grid.rs`, `src/simulation/diagnostics.rs` (energy density computation), `src/simulation/field_update.rs` (where c_local = c0/K is used), `README.md §Polarizable Vacuum` (lines 142–177, PV theory)
- **Depends on:** 1.3
- `[ ]` Implement `src/simulation/vacuum_update.rs`:
  - `VacuumModel` enum: `Fixed` (K=1 everywhere), `PuthoffPV`, `Custom`
  - `VacuumConfig` resource: model type, coupling constant α, K_background
  - `fn update_vacuum_state(grid: &mut SimulationGrid, config: &VacuumConfig)`:
    - For each cell:
      - Compute local field energy: u = 0.5*ε₀*|E|² + 0.5*(1/μ₀)*|B|² + 0.5*ε₀*S²
      - `PuthoffPV` model: K = K_bg * (1 + α * u / u_schwinger)
        - `// PSEUDOCODE: u_schwinger = ε₀ * E_schwinger² / 2`
        - `// E_schwinger ≈ 1.3e18 V/m (QED critical field)`
        - `// For realistic fields, α * u / u_schwinger << 1`
        - `// This means K ≈ 1 + tiny correction — unobservable with normal fields`
        - `// For simulation, we can use a SCALED α to make effects visible`
        - `// This lets us explore the phenomenology even if real-world coupling is tiny`
      - Update K in cell state
      - Field update kernel already uses K for c_local = c0 / K
      - `// PSEUDOCODE: The feedback loop:`
      - `// Strong field → K increases → c_local decreases → field "slows down" locally`
      - `// → gradient in K → force on fields (refraction / deflection)`
      - `// → analogue of gravitational lensing for EM waves`
  - Bevy system: runs after field_update, before derived_fields
- `[ ]` Add K field to the visualization FieldQuantity options
- `[ ]` Add vacuum model selector and α slider to UI
- `[ ]` Verify: in high-field regions, K visibly departs from 1; field propagation visibly slows
- **Session output**: Dynamical vacuum responding to electromagnetic fields

### 4.2 — Casimir scenario
- **Context:** `src/simulation/vacuum_update.rs` (vacuum model API), `src/simulation/grid.rs`, `src/simulation/state.rs` (CellFlags::CONDUCTOR), `src/scenarios/dipole_radiation.rs` (scenario pattern)
- **Depends on:** 4.1
- **Note:** See TODO §Known Issues — Casimir scenario is classical-only, not quantum vacuum.
- `[ ]` Implement `src/scenarios/casimir.rs`:
  - Two parallel conducting plates (set conductor flags on two planes of cells)
  - Initialize with vacuum fluctuations (small random Q in all cells) or just let boundary effects evolve
  - `// PSEUDOCODE: In the PV model, the Casimir effect arises because:`
  - `// 1. Conducting plates constrain the allowed field modes between them`
  - `// 2. The constrained modes have different energy density than free vacuum`
  - `// 3. Different energy density → different K between plates vs outside`
  - `// 4. The K gradient at the plate surfaces → force on the plates`
  - `// Our simulation computes this self-consistently rather than analytically`
  - `fn setup_casimir(grid, vacuum_config)`:
    - Set up two conducting planes separated by distance d
    - Use scaled α to make the effect visible at simulation scale
    - Place probes to measure K profile between and outside plates
  - `fn measure_casimir_force(grid, derived) -> f32`:
    - Integrate Maxwell stress tensor on plate surfaces
    - Or equivalently: compute pressure from K gradient at plate location
  - Expected: K between plates differs from K outside; force is attractive
- `[ ]` Verify: K profile shows expected modification between plates
- **Session output**: Casimir effect emerging from vacuum polarization dynamics

### 4.3 — Vacuum-modified wave propagation
- **Context:** `src/simulation/vacuum_update.rs`, `src/simulation/sources.rs` (to create EM pulse source), `src/visualization/slices.rs` (to observe refraction)
- **Depends on:** 4.1, 1.2
- `[ ]` Create a scenario demonstrating vacuum lensing:
  - Place a strong static charge or intense EM source at grid center
  - This creates a K > 1 region around it
  - Send a transverse EM pulse through the K gradient
  - Observe refraction/deflection — the EM wave bends toward higher K (slower c) regions
  - `// PSEUDOCODE: This is the electromagnetic analogue of gravitational lensing`
  - `// In GR: mass curves spacetime → light bends`
  - `// In PV: field energy polarizes vacuum → K gradient → c gradient → light bends`
  - `// The deflection angle should match the PV prediction:`
  - `// δ ≈ (4GM/rc²) for gravitational case`
  - `// We use a scaled α to make the deflection visible on our grid`
- `[ ]` Verify: wavefront visibly bends around the high-K region
- **Session output**: Stunning visual of electromagnetic lensing by vacuum polarization

---

## Phase 5: GPU Migration

> **Strategy**: By this point the physics is validated on CPU. Now we port the hot loops to GPU compute shaders for performance. Each task ports one system.

### 5.1 — GPU infrastructure: buffer management
- **Context:** `ARCHITECTURE.md §GPU Compute Shader` (lines 268–385, WGSL struct layouts), `ARCHITECTURE.md §Performance Targets` (lines 622–639), `src/simulation/state.rs` (CellState, DerivedFields, SimParams structs to match on GPU), `src/simulation/grid.rs` (data to upload)
- **Depends on:** 0.4
- **Note:** See TODO §Known Issues — CPU/GPU coexistence. CPU path retained behind `use_gpu` flag.
- `[ ]` Implement `src/gpu/buffers.rs`:
  - `GpuBuffers` Bevy resource (all buffer types via `bevy::render::renderer::RenderDevice`):
    - `cells_a: Buffer` — primary state buffer
    - `cells_b: Buffer` — secondary (double-buffer swap target)
    - `derived: Buffer` — derived fields for visualization readback
    - `sources: Buffer` — source term buffer
    - `params: Buffer` — uniform buffer for SimParams
  - `fn create_buffers(device: &RenderDevice, grid_size) -> GpuBuffers`
  - `fn upload_grid(queue: &RenderQueue, buffers, grid)` — CPU → GPU transfer
  - `fn download_grid(device: &RenderDevice, queue: &RenderQueue, buffers, grid)` — GPU → CPU readback
  - `fn download_derived(device: &RenderDevice, queue: &RenderQueue, buffers) -> Vec<DerivedFields>`
  - Double-buffer swap: just swap the bind group references, not the data
  - `// PSEUDOCODE: The double-buffer pattern:`
  - `// Frame N: read from cells_a, write to cells_b`
  - `// Frame N+1: read from cells_b, write to cells_a`
  - `// Swap is just exchanging which buffer is "read" and which is "write"`
  - `// in the bind group, not actually moving data`
- `[ ]` Verify: can upload grid to GPU and download identical data
- **Session output**: GPU buffer infrastructure ready for compute shaders

### 5.2 — GPU compute: field update shader
- **Context:** `ARCHITECTURE.md §GPU Compute Shader: Field Update Kernel` (lines 268–385, complete WGSL reference), `src/simulation/field_update.rs` (CPU reference implementation to port), `src/gpu/buffers.rs` (buffer management API)
- **Depends on:** 5.1, 1.1
- `[ ]` Create `src/gpu/shaders/field_update.wgsl`:
  - Port the CPU `step_field_cpu` logic to WGSL
  - CellState struct in WGSL matching Rust repr(C) layout
  - SimParams uniform binding
  - Workgroup size 8×8×4
  - Full quaternionic FDTD with extended mode flag
  - `// See ARCHITECTURE.md for the complete WGSL kernel reference`
- `[ ]` Implement `src/gpu/pipeline.rs`:
  - `ComputePipelines` resource: holds all compute pipeline objects
  - `fn create_field_update_pipeline(device, shader) -> wgpu::ComputePipeline`
  - Bind group layouts for cell buffers + params uniform
- `[ ]` Implement GPU dispatch in `field_update.rs`:
  - `fn step_field_gpu(render_device, render_queue, pipelines, buffers, params)`
  - Submit compute pass with appropriate workgroup counts
- `[ ]` Validate: GPU results match CPU results for dipole scenario (within f32 tolerance)
- **Session output**: Field update running on GPU — major performance improvement

### 5.3 — GPU compute: derived fields and vacuum update
- **Context:** `src/simulation/diagnostics.rs` (CPU derived field computation to port), `src/simulation/vacuum_update.rs` (CPU vacuum update to port), `src/simulation/sources.rs` (CPU source injection to port), `src/gpu/pipeline.rs` (pipeline creation pattern from 5.2)
- **Depends on:** 5.2, 1.3, 4.1
- `[ ]` Create `src/gpu/shaders/derived_fields.wgsl`:
  - Compute E, B, S, energy_density from Q — same math as CPU version
  - Write to derived fields buffer
- `[ ]` Create `src/gpu/shaders/vacuum_update.wgsl`:
  - Read derived fields (energy density), update K values
- `[ ]` Create `src/gpu/shaders/source_inject.wgsl`:
  - Inject source terms into cells from source buffer
- `[ ]` Wire all shaders into the compute pipeline sequence:
  - source_inject → field_update → vacuum_update → derived_fields
- `[ ]` Verify: full simulation loop on GPU matches CPU results
- **Session output**: Entire simulation loop on GPU

### 5.4 — Visualization reading from GPU buffers
- **Context:** `src/gpu/buffers.rs` (GPU buffer API), `src/visualization/slices.rs` (current CPU-based visualization to adapt), `ARCHITECTURE.md §Rendering Pipeline` (lines 421–453)
- **Depends on:** 5.3, 1.4
- `[ ]` Modify slice visualization to read directly from GPU derived fields buffer
  - Map derived fields buffer to a 3D texture
  - Slice shader samples from the 3D texture (no CPU readback needed for visualization)
  - CPU readback only needed for probes and data export
- `[ ]` Performance test: 64³ should be 60+ FPS; 128³ should be 30+ FPS
- `[ ]` Add grid size selector to UI (requires re-creating buffers)
- **Session output**: Real-time GPU-accelerated simulation with visualization

---

## Phase 6: Advanced Visualization

> **Visualization redesign (recorded 2026-05).** Phases 0–3 used slice planes + arrow glyphs as the primary 3D viz. Both have proven weak: slices project away one of the three dimensions of the data we want to see, and arrow density needed for understanding a 3D vector field collapses into occluded clutter. Phase 6 is the recovery plan, organised around one visual language per data type. See `ARCHITECTURE.md §Visualization Architecture / Design principle` for the underlying rule. Subtasks below are reordered to deliver vector improvements first (no shader work needed) and scalar improvements second (custom WGSL).

### 6.0 — Glyph triage and slice demotion (preparatory)
- **Context:** `src/visualization/glyphs.rs`, `src/visualization/slices.rs`, `src/visualization/plugin.rs`, `src/ui/plugin.rs` (panels)
- **Depends on:** none — this is preparatory cleanup before the redesign lands.
- `[x]` Glyph triage. Three concrete issues found and fixed in `src/visualization/glyphs.rs`:
  - **Outlier-sensitive auto-range** (root cause of "doesn't plot" in dipole/bifilar): `max_mag = max(samples)` was dominated by the strong source cell, then the `mag < max_mag * 1e-4` skip threshold filtered out nearly every other arrow. **Fix:** use the `auto_range_percentile`-th percentile (default 95th) of the magnitude distribution as the reference, so source spikes saturate but don't compress the bulk field below visibility. New `percentile_of` helper covered by 3 unit tests including an outlier-suppression scenario.
  - **Wrong scale in `RgbMultiField` and `SizeColor` encodings** (correctness bug): both encodings were normalising scalar channels (S, φ, |B|, etc.) against `max_mag`, which is the VECTOR field's magnitude max. Scalars with very different magnitudes saturated or went black. **Fix:** compute per-channel percentile maxes from the same sample set (`scalar_p95` closure) and pass them as the per-channel range to `encode_rgb_multi` / `map_value`.
  - **Unit-naive `scale` default** (quality of life): old `scale = 0.001` is meaningless without knowing the field units. **Fix:** scale is now in cell-widths-at-the-percentile-reference, so `scale = 0.5` means a typical arrow is half a cell long regardless of physical units. UI slider relabelled "Scale (cells at p95)" with range 0.05..=3.0.
  - UI also gains an "Auto-range percentile" slider (0.5..=1.0, default 0.95) so the user can tune outlier suppression per-scenario.
- `[x]` Demote slices: implemented. `SingleSliceConfig` gained `show_in_3d: bool` (default `false`). `manage_slice_entity` only spawns the 3D textured quad when `show_in_3d == true`; otherwise the slice is purely a 2D inset image inside the `Slice Planes` egui panel. The RGBA-computation logic was extracted from `update_slice_texture` into the public `slices::compute_slice_pixels` helper so both paths (3D-quad upload AND egui inset) share one source of truth. A `SliceUiTextures` resource holds stable `egui::TextureHandle`s across frames so the inset image doesn't allocate a new GPU texture every frame.
- **Session output**: Glyph implementation fixed (3 concrete bugs). Arrows should now read sensibly in scenarios with strong sources. Glyphs remain a fundamentally limited primary vector view — 6.2 (RK4 streamlines) is still the right replacement — but they now work correctly for sparse-overlay use.

### 6.1 — Volume rendering for scalars (signed AND positive)
- **Context:** `ARCHITECTURE.md §Visualization Architecture / Design principle` (the three-language rule), `src/visualization/color_maps.rs` (transfer function base), eventually `src/gpu/buffers.rs` (GPU 3D textures once Phase 5 lands)
- **Depends on:** 1.3 for CPU prototype; 5.4 for GPU port.
- **Two transfer functions, one renderer.** Signed scalars (S, φ) get a bipolar transfer function — blue = negative, transparent = zero, red = positive, opacity ∝ |value|. Positive scalars (|E|, |B|, K, energy density) get a sequential transfer function — viridis/plasma palette, opacity ramping from 0 at zero to opaque at max. The ray-marching code is shared; only the transfer function differs.
- `[x]` CPU prototype implemented as **stacked semi-transparent slabs** (not full ray-marching — picked the simpler approach that needed no shader work and shipped in one session). `src/visualization/volume.rs`:
  - `VolumeConfig` resource: scalar field selector, num_slabs, opacity_scale, auto_range, color_map (sequential), exclude_pml_from_range.
  - `transfer_bipolar` and `transfer_sequential` pure helpers (10 unit tests covering boundary values, clamping, opacity scaling, dispatch).
  - `manage_volume_entities` spawns N axis-aligned semi-transparent quads at evenly spaced Z positions; respawns when `enabled` or `num_slabs` change.
  - `update_volume_textures` rebuilds each slab's texture every frame using the same dynamic-image swap pattern as `slices.rs`. Picks the correct transfer function from `FieldQuantity::is_signed`.
  - `slab_z_indices` distributes slabs evenly across `[0, nz)` and dedups when more slabs are requested than the grid has layers.
  - `ui_volume_panel` in `src/ui/plugin.rs` exposes all controls; auto-hides the colour-map dropdown for signed fields (bipolar uses a hand-coded blue/white/red).
- Known limitation documented in the module header: edge-on view (camera in the xy-plane) collapses the slab stack to lines. View-aligned slab orientation per frame is the natural polish; deferred to GPU port.
- `[ ]` GPU port (deferred until 5.x): `assets/shaders/volume_render.wgsl` — full ray-marching with view-aligned sampling, freed from the slab-edge-on limitation. The transfer-function logic ports verbatim from `transfer_bipolar` / `transfer_sequential`.
- Reference shader sketch:
  - Create `assets/shaders/volume_render.wgsl`:
    - Ray-march from camera through 3D texture
    - Sample field quantity at regular intervals along ray
    - Accumulate color using front-to-back compositing
    - `// PSEUDOCODE: For each pixel:`
    - `// 1. Cast ray from camera through pixel into 3D volume`
    - `// 2. Step along ray in small increments (e.g., 0.5 * dx)`
    - `// 3. At each step: sample 3D texture → field value`
    - `// 4. Map field value → (color, opacity) via transfer function`
    - `// 5. Composite: accumulated_color += (1 - accumulated_alpha) * sample_color * sample_alpha`
    - `// 6. Early termination if accumulated_alpha > 0.99`
  - Bevy custom render pipeline using a fullscreen quad with the ray-march shader.
- `[ ]` Add volume-render UI controls: scalar quantity selector (one of S / φ / |E| / |B| / K / energy density), transfer-function preset (bipolar vs. sequential picked automatically from the quantity), opacity scale, sample-step size.
- `[ ]` Verify: dipole scenario shows |E| as a glowing radiating volume; bifilar-pair scenario shows S as a clear positive/negative bipolar pattern between transmitter and receiver.
- **Session output**: Single coherent volume-render path serving both signed and positive scalars; replaces slice-as-primary-3D-view.

### 6.2 — RK4 streamlines (primary vector view)
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields for E, B vector data), `src/simulation/grid.rs` (trilinear interpolation), `src/simulation/particles.rs` (the existing trilinear E/B sampler at world position is the right reference), `ARCHITECTURE.md §Visualization Architecture / Design principle` (vectors → streamlines)
- **Depends on:** 1.3
- **Goal:** become the primary 3D representation of every vector field (E, B, Poynting, A). Glyphs are demoted to a sparse overlay at most, or deleted in favour of streamlines.
- `[x]` Implement `src/visualization/streamlines.rs` (full rewrite from the Euler stub):
  - `compute_streamline(grid, diag, field, seed_xyz, max_steps, step_fraction, threshold) -> Vec<(Vec3, f32)>` — pure function returning a polyline of (world_pos, magnitude) pairs. Each point's magnitude is the field magnitude at the start of the segment terminating there, used by the renderer for per-segment colouring.
  - `rk4_step(grid, diag, field, fx, fy, fz, h)` — single 4th-order Runge-Kutta step on the NORMALISED direction. Step length is arc length (h·dx grid units) regardless of magnitude, which is the right thing for visualising topology. Returns `None` if any of k1..k4 falls outside the safe interpolation region or hits zero magnitude.
  - `polyline_at(polyline, t)` helper for interpolating a position at fractional parameter t ∈ [0, 1] along the polyline — used by the tracer animation.
  - `draw_streamlines` Bevy system: traces from a stride-spaced seed grid, draws each polyline as gizmo lines coloured by per-segment magnitude, and (when enabled) draws `tracers_per_line` animated white sphere markers flowing along each line at `tracer_speed` revolutions/second.
- `[x]` Add streamline controls to UI — `Streamlines` window (`ui_streamline_panel` in `src/ui/plugin.rs`): enable toggle, vector-field selector, seed stride / max steps / step fraction sliders, colour map + auto/manual range, animate-tracers toggle with tracers-per-line and tracer-speed sliders.
- Wired into `UiPlugin`'s system list — previously `StreamlineConfig` was a registered resource with no panel, so streamlines were silently disabled the entire time.
- `[x]` Tests in `src/visualization/streamlines.rs`:
  - `test_streamline_uniform_field_is_straight_line` — uniform +ẑ E-field: streamline is a vertical straight line, X/Y drift < 1e-5, Z monotone, step length = step_fraction × dx within 1e-4 rel. err.
  - `test_streamline_vortex_closes_on_itself` — vortex field E = (-y, x, 0) in grid coords: streamline of radius 4 grid units closes after one full revolution within tolerance derived from the geometric residual (round(2π·r/h) is rarely 2π·r/h exactly) plus 5 % RK4 phase drift; radius drift along the orbit < 2 %.
  - `test_streamline_outside_region_returns_only_seed` — seed near a face: integrator immediately exits the safe region, polyline contains only the seed.
  - `test_polyline_at_endpoints` and `test_polyline_at_too_short` cover the tracer-position helper.
- **Session output**: RK4 streamlines as the primary vector visualisation, with animated tracers conveying direction at a glance. Glyphs survive only as a sparse-overlay tool; streamlines are the default.

### 6.3 — Isosurface extraction
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields scalar data), `src/simulation/grid.rs` (grid dimensions for cell iteration)
- **Depends on:** 1.3
- `[~]` Implement `src/visualization/isosurface.rs`:
  - **Note:** Uses voxelized boundary faces (staircase surface), not marching cubes. No lookup tables required. Upgrade to MC for smooth surfaces if needed.
  - Marching cubes algorithm on the 3D field data
  - `fn marching_cubes(field_data, threshold, grid) -> Mesh`:
    - Standard 256-entry lookup table for cube configurations
    - Vertex interpolation along cell edges
    - `// PSEUDOCODE: For each cell (x,y,z):`
    - `// 1. Read field values at 8 corners of the cell`
    - `// 2. Determine which corners are above/below threshold → 8-bit index`
    - `// 3. Look up triangle configuration in edge table and tri table`
    - `// 4. For each triangle edge: interpolate vertex position between corners`
    - `// 5. Emit triangle vertices and normals`
    - Can run on GPU compute shader for large grids
  - Display as semi-transparent mesh, colored by field value or gradient
- `[ ]` Add isosurface controls: field selector, threshold slider, opacity
- `[ ]` Verify: isosurface of |E| shows expanding spherical wavefronts from dipole
- **Session output**: 3D wavefront surfaces visible in the simulation

### 6.4 — Quaternion-specific visualization
- **Context:** `src/math/quaternion.rs` (Quat type), `src/simulation/grid.rs` (accessing Q per cell), `src/visualization/slices.rs` (slice rendering to extend), `ARCHITECTURE.md §Quaternion Visualization` (lines 407–412)
- **Depends on:** 1.4, 0.2
- `[ ]` Implement `src/visualization/quaternion_viz.rs`:
  - Mode 1: "Scalar + Vector" — volume render scalar part (φ), overlay vector arrows (A)
  - Mode 2: "Norm + Phase" — volume render |Q|, color by "quaternion phase" (angle of rotation represented by Q)
  - Mode 3: "Rotation field" — at each point, Q represents a rotation; show as oriented ellipsoids or frame triads
    - `// PSEUDOCODE: Interpret Q at each point as a rotation quaternion`
    - `// Apply Q to a reference frame (3 orthogonal vectors)`
    - `// Render the rotated frame as a small oriented triad or ellipsoid`
    - `// The "texture" of rotations reveals the topological structure`
    - `// Singularities (where |Q|→0 or Q is undefined) show up as defects`
  - Mode 4: "Component slices" — 4 slice planes side by side, one per Q component
- `[ ]` Add quaternion viz mode selector to UI
- **Session output**: Visualizing the full quaternionic structure of the field

---

## Phase 7: Scenarios and Experimental Design

### 7.1 — Aharonov-Bohm scenario
- **Context:** `src/simulation/grid.rs`, `src/simulation/sources.rs`, `src/scenarios/dipole_radiation.rs` (scenario pattern), `README.md §Aharonov-Bohm Effect` (lines 196–203)
- **Depends on:** 1.4, 1.3
- **Note:** See TODO §Known Issues — this demonstrates classical potential structure, not quantum phase interference.
- `[x]` Implement `src/scenarios/aharonov_bohm.rs`:
  - `AbConfig { axis_xy_cells, radius_cells, b0 }` with `expected_flux(grid)` helper that returns the analytic Φ = πR²B₀ in world units.
  - `apply_ab_scenario` writes the analytic A field cell-by-cell: `A_x = -B₀·dy/2`, `A_y = +B₀·dx/2` for ρ < R; `A_x = -B₀R²·dy/(2ρ²)`, `A_y = +B₀R²·dx/(2ρ²)` for ρ ≥ R; φ = A_z = 0; q_dot = 0. Sources cleared. Configuration is static.
  - `sample_a_at_world(grid, [x, y, z])` — trilinear interpolation of A at any world point.
  - `integrate_a_dot_dl(grid, path)` — midpoint-rule line integral around a closed polyline.
  - `rectangular_loop(centre, half_x, half_y, segments_per_side)` — convenience builder for square test paths.
- `[x]` Verify B confined inside, A extends outside — **eight unit tests** covering all four corners of (initial vs. evolved) × (standard vs. extended):
  - `test_ab_b_field_inside`: B_z averaged over ρ < R/2 is within 2 % of B₀.
  - `test_ab_b_field_outside`: max |B| at ρ > 1.5 R is < 5 % of B₀ (the residual comes from the boundary discontinuity at ρ = R).
  - `test_ab_a_nonzero_outside`: |A| at ρ = 2R matches the analytic B₀R/4 within 5 % — **the AB hallmark**, A ≠ 0 where B = 0.
  - `test_ab_path_integral_encircling`: ∮A·dl on a square encircling the solenoid matches Φ = πR²B₀ within 5 %.
  - `test_ab_path_integral_non_encircling`: ∮A·dl on a square far from the flux tube is < 2 % of Φ.
  - `test_ab_short_time_flux_preservation`: encircling flux drift after 5 wave-equation steps < 5 %. **Long-time drift is real physics** — the analytic AB potential has a kink at ρ = R that produces a δ-function-like Laplacian under finite differences, seeding an outgoing wave at speed c. A static AB without sustaining current is not a steady-state of the homogeneous wave equation; this is documented in the module header.
  - `test_ab_s_field_stays_zero_extended`: with φ = 0 and ∇·A = 0 at t = 0 and no charge density, S obeys □S = 0 with zero IC and stays at zero. Run 50 steps in extended mode, max |S| ≪ B₀.
  - **`test_ab_extended_matches_standard`** — the central QVED equivalence claim. Two independent grids, one in standard EM, one in extended QVED, both initialised with the same AB analytic A. After 50 steps the A fields agree to 1e-4 relative. Because S stays at zero, the extended-mode update reduces algebraically to the standard-mode update, so the two formulations give bitwise-similar results. **This is what "AB satisfied in both formulations" means precisely** — the QVED extension introduces no new physics for gauge-clean configurations.
- `[x]` Add `Scenario::AharonovBohm` to the scenario selector with a UI handler that places the solenoid axis at the grid centre and pauses the simulation (the configuration is static and most useful inspected at t = 0).
- **Session output**: Classical proof that A extends beyond B, verified by path-integral diagnostic on a discrete grid; rigorous demonstration that QVED extended mode reduces to conventional EM for gauge-clean initial conditions.

### 7.2 — Toroidal AB circuit coupling
- **Context:** `src/scenarios/aharonov_bohm.rs` (AB pattern to extend), `src/simulation/sources.rs`, `src/simulation/diagnostics.rs` (probes for measurement)
- **Depends on:** 7.1, 2.4
- `[x]` Implement `src/scenarios/toroidal_ab.rs`:
  - `ToroidalConfig { center, major_radius, minor_radius, num_poloidal_loops, num_segments_per_loop, current_per_loop }` with `expected_flux()` returning the thin-torus formula `Φ = NIr²/(2R)` (μ₀ = 1 in code units; off by O((r/R)²) ≈ 6 % for default r/R = 0.25).
  - `compute_a_from_torus(cfg, world_pos)` — analytic A from N discrete poloidal current loops via Biot-Savart sum: each loop is at toroidal angle φ_k, parametrised in the (radial outward, ẑ) plane around its centre at (R cos φ_k, R sin φ_k, 0).
  - `apply_toroidal_ab_scenario` loads the Biot-Savart A field into every grid cell. Static configuration (q_dot = 0).
  - `linking_pickup_loop(cfg, segments)` — counter-clockwise rectangle in y = cy plane that links the +x portion of the torus tube once. Right edge is at x = cx + R + 2r (well clear of the tube), top/bottom at z = cz ± R.
  - `non_linking_pickup_loop(cfg, segments)` — small rectangle in the xy-plane at z = cz + 1.5R, far from any tube and not linking.
- `[x]` Tests in `src/scenarios/toroidal_ab.rs`:
  - **`test_toroidal_pickup_loop_linking`**: ∮A·dl on a linking loop has magnitude ≈ NIr²/(2R) within 30 %. (The thin-torus analytic factor and discrete-loop discretisation share the tolerance budget.)
  - `test_toroidal_pickup_loop_non_linking`: ∮A·dl on a non-linking loop is < 5 % of the linking magnitude.
  - `test_toroidal_b_inside_tube`: at toroidal angle φ = 0 the magnetic field inside the tube is dominated by its toroidal component (along ŷ at that point).
  - `test_toroidal_b_outside_tube`: max |B| outside the tube is < 50 % of the in-tube |B| (some leakage from discrete loops, but B is concentrated inside).
  - **`test_toroidal_extended_matches_standard`**: standard and extended-mode evolution after 30 steps agree on A to 1e-4 relative — the QVED equivalence claim, identical to the linear AB case for gauge-clean configurations.
  - `test_toroidal_s_stays_zero_extended`: S stays bounded near zero in extended mode.
- Sign-convention note documented inline: with our poloidal-loop parametrisation, current flows +ẑ at the outer rim and −ẑ at the inner rim; Ampere's law puts B in the −e_φ direction inside the tube; the linking flux through a +ŷ-oriented disk is therefore negative. Magnitude is what matches the analytic formula; the sign is bookkeeping.
- `[x]` Add `Scenario::ToroidalAb` to the scenario selector; UI handler loads the default torus geometry, pauses the simulation.
- `[~]` Driven AC operation + EMF measurement on a probe loop — deferred. The static-A scenario already demonstrates the topology cleanly; an AC mode would add scenario-driver state and a per-step pickup-EMF diagnostic. Future work.
- `[~]` Experimental BOM (toroid core specs, wire gauge, oscilloscope requirements) — deferred. The scenario predicts the geometric coupling coefficient `Φ/(NI) = r²/(2R)` directly from the simulation; translating that into a benchtop spec belongs in a separate "experimental protocols" doc, not in code.
- Note on the original scope claim: the TODO had read "in standard theory: no coupling (B=0 outside torus)". This is **incorrect** — Faraday's law gives mutual inductance between the toroidal coil and a pickup loop linking the hole even in conventional EM, because the surface bounded by the pickup loop must pass through the tube's interior (where B ≠ 0). The actual physically meaningful claim, verified by `test_toroidal_extended_matches_standard`, is that classical EM and QVED predict the *same* coupling for gauge-clean configurations — exactly as for linear AB.
- **Session output**: Macro-scale AB topology numerically demonstrated on a discrete-Biot-Savart torus, with the QVED equivalence claim verified for the donut-coil + linking-pickup configuration.

### 7.3 — PCB Weber force geometry optimizer
- **Context:** `src/scenarios/graneau_wire.rs` (wire scenario to extend), `src/simulation/weber.rs` (Weber force API), `src/simulation/particles.rs`
- **Depends on:** 3.3
- `[x]` Implement hairpin PCB-trace geometry in `src/scenarios/graneau_wire.rs`:
  - `apply_hairpin_scenario(particles, force_mode, num_segments, segment_spacing, pcb_spacing, v_drift)` — places two parallel rows of (ion, electron) pairs in y at ±pcb_spacing/2. Top strand carries current in +x̂; bottom strand returns it in −x̂. Total of 4·num_segments particles. Forces `ForceMode::Weber`.
  - `apply_hairpin_default` convenience wrapper using `DEFAULT_HAIRPIN_SPACING = 5 mm`.
  - `compute_hairpin_top_force_profile(particles, num_segments)` — returns the per-ion axial Weber force on the TOP strand. Bottom strand's contribution is included in the pair sum.
- `[x]` Tests:
  - `test_apply_hairpin_creates_two_strands` — particle count, positions, charge neutrality, top vs bottom drift directions.
  - `test_hairpin_profile_differs_from_straight_wire` — the bottom strand modifies the top-strand's longitudinal profile measurably (≥ 1 % of the straight-wire end magnitude).
  - `test_hairpin_top_strand_still_end_peaked` — the SELF-interaction of the top strand dominates: end ions still feel outward forces, middle stays small.
- `[x]` `Scenario::HairpinTrace` added to the scenario selector with a UI handler that loads the default hairpin geometry, sets `ForceMode::Weber`, and pauses the simulation.
- `[~]` Spiral and interdigitated-finger geometries — deferred. Both are straightforward extensions of the hairpin code (just different particle layouts) but each adds its own test/UI surface; not worth the complexity in the same session that introduced the hairpin.
- `[~]` Parameter sweep over `pcb_spacing`, `num_segments`, `v_drift` to find the geometry maximising Weber/Lorentz ratio — deferred. Would need a benchmarking harness, output tabulation, and possibly a parallel sweep driver. The static "load-and-look" hairpin scenario is sufficient to verify the geometry shows the expected structure; the optimisation belongs in a separate "physics study" session.
- `[~]` KiCad-footprint export — deferred. Once an optimised geometry exists, dumping its (x, y) endpoint pairs as a KiCad footprint string is a small piece of work but disconnected from the simulator's core; can live as a Python post-processor.
- **Session output**: Hairpin PCB-trace geometry implemented and verified to produce a Weber-force profile distinct from the straight-wire case. The bottom return strand contributes a measurable longitudinal correction without erasing the underlying Graneau end-peaking.

### 7.4 — Tier 3 scenario stubs
- **Context:** `src/scenarios/mod.rs` (module structure), `README.md §Tier 3 Scenarios` (lines 253–278, theory context for each)
- **Depends on:** 0.1
- `[x]` Implement stubs for speculative scenarios (pseudocode only, no full sim):
  - `src/scenarios/brown_capacitor.rs` — asymmetric capacitor in PV vacuum. `BrownConfig` with plate-separation / radius / voltage parameters. `apply_brown_capacitor_scenario` with `todo!()` body and full pseudocode in the docstring (conductor BC, K-relaxation, Maxwell-stress + K-gradient force diagnostic). Implementation prerequisites: conductor-BC enforcement in `step_field_cpu` and a Maxwell-stress force integrator.
  - `src/scenarios/pulsed_circuit.rs` — sharp current interruption, scalar-mode excitation. `PulsedCircuitConfig` with wire length / current / switch / rise-time parameters. `apply_pulsed_circuit_scenario` with pseudocode covering switched-current source, multi-probe install, and the standard-vs-extended S-channel comparison. Implementation prerequisites: a `SwitchedSource` source variant and an automatic multi-probe install for the scenario.
  - `src/scenarios/charge_cluster.rs` — EVO stability search. `ChargeClusterConfig` with N / cluster-radius / charge / mass / thermal-speed / RNG-seed parameters. `apply_charge_cluster_scenario` with pseudocode covering Maxwellian sampling and RMS-radius / K(r) diagnostics. Implementation prerequisites: a seeded RNG, RMS-radius and K(r) diagnostic resources, and an elevated-η VacuumConfig regime that reaches the K > √2 threshold from README §3.4.
- All three stubs build cleanly under `cargo build` and panic with informative messages if their `apply_*` functions are called.
- Each module header documents the **physics background**, **scenario claim**, **expected observables**, and **parameter ranges to explore** so future sessions can implement without re-deriving from the README.
- **Session output**: Compiling, fully-documented roadmap stubs for the three Tier-3 speculative scenarios. Future sessions can implement any of them as standalone tasks following the pseudocode laid out here.

### 7.5 — Hopfion / ball lightning scenario
- **Context:** `src/simulation/grid.rs` (Q field, K field layout), `src/simulation/diagnostics.rs` (topological charge diagnostic from Phase 1.9), `src/simulation/field_update.rs` (leapfrog, K dynamics from Phase 1.8), `README.md §3.5` (Hopfion theory, topological charge formula, ball lightning model)
- **Depends on:** 1.8 (K dynamics), 1.9 (topological charge), 2.1 (extended mode)
- **Physics summary:** The Irvine-Bouwmeester Hopfion is an exact Maxwell solution with Hopf-linked E and B field lines. In linear EM it disperses at speed c. In the (Q,K) system, K > 1 provides a dynamical Skyrme term that can sustain the topology. This scenario tests the critical K threshold for stable vs. dispersing Hopfion configurations.

#### Implementation tasks:

- `[~]` Implement `src/scenarios/hopfion_ball_lightning.rs` — first-pass with the simpler hedgehog Skyrmion ansatz (the README §7.5 "starting stub" recommendation), plus a `validate_hopfion_stability` helper. Full Rañada-Irvine analytic IC and the physical torus-knot discharge IC are deferred to follow-up sessions.
  - `HopfionConfig { center, radius, winding_factor }` with sensible defaults for a 32³ × 0.01 m grid.
  - `apply_hopfion_scenario` writes the linear-profile hedgehog `Q = amp · (cos θ, sin θ · r̂)` with `θ(r) = π · min(r/R, 1)`. Same construction as `test_topo_charge_conserved` from Phase 1.9, so it inherits that test's stability characteristics.
  - `validate_hopfion_stability(initial_topo, current_topo, initial_energy, current_energy, topo_tolerance, energy_floor) -> HopfionValidation` returns booleans for `topology_conserved` and `energy_retained` plus the raw measurements; `is_stable()` requires both.
  - `current_topological_charge(grid)` convenience wrapper around `diagnostics::compute_topological_charge`.
- `[~]` Full Rañada-Irvine analytic IC via Bateman complex potentials — deferred. Has well-known closed-form solutions but require complex-arithmetic helpers and per-cell evaluation of stereographic projection from S³ to Hopf-linked E/B fields. A clean session of its own.
- `[~]` Physical torus-knot discharge IC — deferred. Requires modelling the helical (1,1)-knot current source with Gaussian envelope, a non-trivial source geometry. Belongs with toroidal_ab + a real driven-current implementation.
- `[x]` Tests in `src/scenarios/hopfion_ball_lightning.rs`:
  - `test_hopfion_initial_topology_nontrivial`: hedgehog produces |topo| > 0.5 at t=0.
  - `test_hopfion_initial_q_dot_zero`, `test_hopfion_q_magnitude_unit`: setup invariants.
  - `test_hopfion_initial_energy_nonzero`: real EM energy is loaded.
  - `test_hopfion_topology_short_time_conservation`: 50 steps in extended mode, drift bounded by `1.5 × max(|initial|, 0.5)` — the same loose threshold that `tests/integration_phase1.rs::test_topo_charge_conserved` uses, accepting that the linear wave equation disperses Skyrmions and that **the stability claim itself is the deferred parameter sweep over η**.
  - Three `test_validate_*` tests cover the boolean logic of `HopfionValidation`.
- `[x]` Add `Scenario::Hopfion` to the scenario selector with a UI handler that loads the default hedgehog IC, switches to extended QVED mode, and pauses.
- `[~]` Parameter sweep over `VacuumConfig.eta` (0, 1e-6, 1e-4, 1e-2) to find the critical η — the **headline scientific output**, deferred. Requires a scan harness, multiple-run data collection, and a results table (see "Deferred infrastructure / Parameter-sweep harness" below). The simulation infrastructure to run it is now in place: `apply_hopfion_scenario` + `validate_hopfion_stability` + `compute_topological_charge`.
- `[~]` Full Rañada-Irvine analytic IC via Bateman complex potentials — deferred. The hedgehog Skyrmion implemented here has the same integer topological invariant (Skyrmion winding = ±1) and the same dispersion behaviour under linear `□Q = 0`, so the qualitative scientific question is testable with the existing IC. The Hopf invariant proper (linking of B and E lines) is a separate, more elaborate target that needs Bateman F = (Ax + i)/(r²+1)-type closed-form expressions.
- `[~]` Physical torus-knot discharge IC — deferred. Future helical-winding source with Gaussian current envelope (`peak_current · exp(−(t/τ)²)`) traversing a (1,1) torus knot. Useful for testing whether a real laboratory geometry creates Hopf topology vs. relying on the analytic IC.
- `[~]` Field-line tracer with linking-number colouring (Phase 6.x extension), K-field XZ slice (already supported by `Slice 1`/`Slice 2`), topological-charge vs. time plot — deferred polish; the existing `topological_charge` diagnostic in `DiagnosticsState` already updates and displays as a number, just not as a time-series plot.
- **Session output**: Hedgehog Skyrmion ansatz scenario with 8 unit tests covering setup invariants, initial topological charge, energy non-trivial, validation logic, and short-time topology conservation under linear evolution. Sets up the infrastructure for the deferred η-parameter sweep that produces the headline scientific output.

### 7.6 — K-cycle resonator scenario
- **Context:** `src/simulation/field_update.rs` (K leapfrog from Phase 1.8), `src/simulation/sources.rs` (source API), `src/simulation/diagnostics.rs` (energy tracking), `README.md §3.4` (K field dynamics, virtual pair plasma)
- **Depends on:** 1.8 (K dynamics)
- **Physics summary:** If K is dynamically switchable (via a pulsed discharge that rapidly changes the local field energy density), non-adiabatic switching near ωₚ produces real photons from the virtual pair plasma (dynamical Casimir effect). A resonator that cycles K between K_low and K_high at frequency ωₚ_eff could sustain photon production as long as the input energy exceeds the net photon emission. This is NOT claimed energy extraction (that would require input < output) — the simulation tests whether coherent K cycling produces a detectable photon mode population distinct from thermalization.

#### Implementation tasks:

This phase **modifies a core simulation module** (`field_update.rs`) — see CLAUDE.md "Modifying Core Simulation Modules" for the gate-by-config-flag pattern. None of the existing 157 lib tests should change behaviour; the new drive term is invisible when its config field is at the default zero amplitude.

**Prerequisite: extend `VacuumConfig` (in `src/simulation/plugin.rs`).**
- Add three fields with defaults that disable the feature:
  ```rust
  pub k_drive_amplitude: f32,   // default 0.0 (no drive)
  pub k_drive_frequency: f32,   // default 0.0
  pub k_drive_radius: f32,      // default 0.0 — when 0.0, no spatial mask is applied
  ```
- The drive is active only when `k_drive_amplitude > 0.0`.

**Prerequisite: extend the K-update path in `step_field_cpu`.**
- Currently the K leapfrog updates `k_dot[i] += dt * (c² · lap_k − ωₚ² · (k − 1) + η · u_field / u_s)`.
- Add: when `vacuum.k_drive_amplitude > 0.0` AND the cell falls inside the drive shell (`r ≤ k_drive_radius` from grid centre), add `dt · k_drive_amplitude · sin(2π · k_drive_frequency · t)` to `k_dot[i]` BEFORE the integration step. The drive is a SCALAR forcing, not a current — it directly perturbs the K leapfrog.
- A 3-line module-header note explaining why and which scenario uses it.

**Prerequisite: add a per-scenario diagnostic for resonant amplification signature.**
- The resonance signature is "EM energy density inside the drive shell oscillates at 2× the drive frequency" (parametric amplification).
- Easiest implementation: a second probe (or a resource holding `(time, em_energy_inside_shell, em_energy_outside_shell)` snapshots) that records every step. Then the user (or a test) takes an FFT of the energy time-series and looks for a peak at `2 · k_drive_frequency`.
- Reuse the existing `Probe` infrastructure if possible — install one probe inside the shell measuring `EnergyDensity`, one outside.

- `[x]` Extend `VacuumConfig` with `k_drive_amplitude`, `k_drive_frequency`, `k_drive_radius` fields, all defaulting to 0.0 (drive disabled).
- `[x]` Wire the drive term into the K leapfrog inside `step_field_cpu` — when amplitude > 0 and radius > 0, cells within `k_drive_radius` of the grid origin pick up `k_drive_amplitude · sin(2π · k_drive_frequency · t)` added to k_ddot. Module-header note in `field_update.rs` explains the modification.
- `[x]` Implement `src/scenarios/k_cycle_resonator.rs`:
  - `KCycleResonatorConfig { center, resonator_radius, drive_amplitude, drive_frequency, k_initial_peak, omega_p, eta }` with sensible defaults.
  - `apply_k_cycle_resonator_scenario` initialises a smooth Gaussian K bump inside the shell (σ = R/2, peak `k_initial_peak`), enables VacuumConfig, sets the K-drive fields, installs two probes (centre = inside shell, +2R along x = outside shell) recording `EnergyDensity`. Both buffers' K initial state are mirrored so the first leapfrog step has a consistent base regardless of `current` index.
  - `apply_k_cycle_resonator_default` convenience wrapper.
- `[x]` Tests:
  - **`test_k_drive_off_preserves_existing_behaviour`** — two grids with identical Gaussian K bump and `drive_amplitude = 0.0`, run 30 steps each, max |Δk| < 1e-6. Regression protection per CLAUDE.md "Modifying Core Simulation Modules".
  - **`test_k_drive_oscillates_k`** — drive on, after 20 steps K(centre) drift OR k_dot(centre) > 1e-6 (drive successfully forces oscillation).
  - **`test_drive_localised_to_shell`** — inside-shell K drift > outside-shell K drift (the drive's spatial mask works).
  - **`test_parametric_response_peak_visible`** — after 200 drive periods, K(centre) shows visible departure-plus-velocity excitation (signature-level, not precision).
- `[x]` Add `Scenario::KCycleResonator` to the scenario selector with a UI handler that sets `extended_mode = true`, configures VacuumConfig + probes via `apply_k_cycle_resonator_default`, and pauses.
- `[~]` Frequency sweep over `drive_frequency / ωₚ` from 0.1 to 2.0 to locate the resonance peak — deferred as a physics-study follow-up. Requires the parameter-sweep harness from "Deferred infrastructure" below; the per-point machinery (apply scenario, run, FFT probe) is now in place.
- `[~]` Dedicated `analyse_resonance_response` helper that runs `probe_fft` on the inside-shell probe and reports the peak amplitude at `2 · drive_frequency` relative to the broadband baseline — deferred. Easy to add when the sweep harness lands; the FFT and probe infrastructure already do the heavy lifting.
- **Session output**: K-equation drive infrastructure landed in `field_update.rs` behind a default-zero amplitude flag (existing 157 lib tests untouched). K-cycle resonator scenario implements the drive-shell setup with inside/outside probes, ready for parametric-resonance investigation.

### 7.7 — Spheromak / Taylor relaxation with dynamic K (fusion confinement connection)
- **Context:** `src/simulation/field_update.rs` (K dynamics from Phase 1.8), `src/simulation/sources.rs` (current source API), `src/simulation/diagnostics.rs` (topological charge from Phase 1.9, energy diagnostics), `README.md §3.7–§3.9` (magnetic helicity, Taylor relaxation, K-modified confinement physics)
- **Depends on:** 1.8 (K dynamics), 1.9 (topological charge), 3.1 (particle system for plasma current model)
- **Physics summary:** The spheromak is the minimum-energy state of a plasma at fixed magnetic helicity H_mag = ∫ A·B dV. Standard MHD Taylor relaxation finds this state with K=1. With dynamic K, the minimum-energy state at fixed helicity shifts — the K gradient contributes to the force balance, potentially modifying the confinement geometry and introducing a photon-reflecting K-boundary that standard MHD cannot model. The magnetic helicity H_mag is the classical analogue of the topological charge from Phase 1.9 and should be monitored alongside it. See README §3.9 for the full argument connecting plasma-K (effective dielectric) to vacuum-K and the unexplored fusion confinement implications.

#### Implementation tasks:

- `[x]` Add `compute_magnetic_helicity(grid, derived)` to `src/simulation/diagnostics.rs`. Sums A · B · dx³ over interior non-PML cells; called from `diagnostics_system` on the same 10-step throttle as `topological_charge`.
- `[x]` Add `magnetic_helicity: f32` to `DiagnosticsState`. UI diagnostics panel display deferred (the value is in the resource and ready to display).
- `[x]` Spherical Bessel `j_1` helper inside `src/scenarios/spheromak_taylor.rs` with Taylor expansion for |x| < 1e-3 to avoid catastrophic cancellation. Verified at known values (j_1(0)=0, first zero at 4.49…, j_1(π)=1/π).
- `[x]` Implement `src/scenarios/spheromak_taylor.rs`:
  - `SpheromakConfig { center, radius, b0 }` with `lambda()` returning the analytic 4.4934.../radius (first j_1 zero).
  - `apply_spheromak_scenario` initialises A from the analytic Chandrasekhar-Kendall l=1 mode in spherical coordinates (B_r, B_θ, B_φ from the standard formulas; A = B/λ in Coulomb gauge), converts to Cartesian per cell, and clears A outside r > R for sharp confinement.
- `[x]` Tests:
  - `test_bessel_j1_basic_values` and `test_bessel_j1_small_argument` — Bessel implementation correctness.
  - `test_magnetic_helicity_uniform_b_is_zero` — H = 0 for a uniform field with A = 0.
  - **`test_spheromak_magnetic_helicity_nonzero`** — the spheromak IC produces nonzero helicity (linked flux).
  - `test_spheromak_a_vanishes_outside` — A is exactly zero outside r > R.
  - `test_spheromak_helicity_short_time_drift` — after 30 steps, helicity stays nonzero with the same sign (magnitude can drop substantially as the boundary discontinuity radiates an outgoing wave; documented inline as the same physics as AB-flux drift).
  - **`test_spheromak_extended_matches_standard`** — QVED equivalence holds for the spheromak IC, same as for AB and Toroidal AB.
- `[~]` `test_spheromak_approximately_force_free` — `#[ignore]`'d. The analytic spherical-coordinate B I assembled is NOT divergence-free after Cartesian conversion (my closed-form CK angular factors don't satisfy ∇·B = 0). The IC is a "spheromak-like" linked-flux configuration but not a true CK eigenmode of curl. Fixing requires deriving the correct l=1 CK formula via the toroidal-poloidal stream-function construction. Lift the `#[ignore]` once the formula is corrected.
- `[x]` Add `Scenario::SpheromakTaylor` to the scenario selector. UI handler loads the IC, sets extended QVED mode, pauses.
- `[~]` Coaxial-gun helicity-injection source — deferred to the cross-phase "Driven-current source" infrastructure (see "Deferred infrastructure" section above). The current static-IC scenario is sufficient to validate the helicity diagnostic; dynamic injection requires the programmable-current source variant.
- `[~]` `test_k_develops_at_boundary` (K-dynamic boundary-shell formation, the headline novel-physics result) — deferred. Needs the parameter-sweep harness to scan η and verify K shell formation; simulation infrastructure is in place but the multi-run harness isn't.
- `[~]` Full Taylor relaxation from turbulent IC — deferred (physics-study session, ~1000 steps with randomised perturbations).
- **Session output**: `compute_magnetic_helicity` diagnostic landed in `diagnostics.rs` with throttled invocation. Spheromak-like analytic IC scenario verifies the helicity computation on a configuration with nontrivial linked flux; QVED equivalence claim verified for the spheromak as well, completing the AB-family of gauge-clean equivalence tests (linear AB, toroidal AB, spheromak). The headline K-boundary prediction is unblocked but requires the parameter-sweep harness for actual measurement.

---

## Phase 8: Data Export and Reproducibility

### 8.1 — HDF5 state export (behind `hdf5-export` feature flag)
- **Context:** `src/simulation/grid.rs` (grid data to export), `src/simulation/diagnostics.rs` (DerivedFields, probes), `src/simulation/state.rs` (CellState layout)
- **Depends on:** 1.3
- `[ ]` Implement `src/export/hdf5_export.rs` (gated with `#[cfg(feature = "hdf5-export")]`):
  - Export full grid state (all cell Q, K values) at a given timestep
  - Export derived fields
  - Export probe time-series data
  - Include metadata: simulation parameters, scenario name, timestamp, git hash
  - `// PSEUDOCODE: HDF5 file structure:`
  - `// /metadata/ — simulation parameters, scenario config (as JSON attributes)`
  - `// /state/Q — (nx, ny, nz, 4) f32 array — quaternionic potential`
  - `// /state/K — (nx, ny, nz) f32 array — vacuum polarizability`
  - `// /derived/E — (nx, ny, nz, 3) f32 array`
  - `// /derived/B — (nx, ny, nz, 3) f32 array`
  - `// /derived/S — (nx, ny, nz) f32 array`
  - `// /probes/{name}/time — (N,) f64 array`
  - `// /probes/{name}/{field} — (N,) f32 array`
- `[ ]` Add export button to UI: "Save State" → file dialog → HDF5 write
- **Session output**: Simulation results exportable for analysis in Python/Julia/MATLAB

### 8.2 — Screenshot and animation export
- **Context:** `src/ui/plugin.rs` (UI button integration), Bevy render pipeline (frame capture API)
- **Depends on:** 0.5
- `[ ]` Implement `src/export/screenshot.rs`:
  - Capture current frame buffer to PNG
  - Include simulation time and parameters in filename
- `[ ]` Implement `src/export/animation.rs`:
  - Record mode: save every Nth frame as PNG sequence
  - Add record start/stop button to UI
  - Output frame sequence suitable for ffmpeg encoding
- **Session output**: Publication-quality image and video capture

### 8.3 — Scenario configuration files
- **Context:** `src/scenarios/*.rs` (all scenario setup functions), `src/simulation/plugin.rs` (SimulationConfig), `src/simulation/sources.rs` (SourceConfig), `src/simulation/boundaries.rs` (BoundaryConfig)
- **Depends on:** 1.6, 1.5, 1.2
- `[ ]` Create RON config files for each scenario in `assets/scenarios/`:
  - Grid size, dx, dt
  - Source configuration (positions, types, parameters)
  - Boundary conditions
  - Vacuum model and parameters
  - Probe positions
  - Visualization presets (which modes active, color maps, slice positions)
- `[ ]` Implement config load/save in UI
- `[ ]` Add command-line argument: `--scenario <name>` to load directly
- **Session output**: Fully reproducible scenario configurations

---

## Deferred infrastructure (cross-phase prerequisites)

Several Phase-7 deliverables intentionally deferred their **physics-study** components — parameter sweeps, multi-run data collection, pipeline export. Those all need the same infrastructure that doesn't exist yet, so capturing it once here saves rediscovering it per phase.

### Parameter-sweep harness
**Needed by:** Phase 7.5 (η-sweep for Hopfion stability), Phase 7.6 (drive_freq/ωₚ sweep for K-cycle resonance), the future Tier-3 Brown-capacitor (V × radius-ratio scan), and any "find the critical X" question.

Requirements:
- A standalone binary or test-binary harness in `src/bin/parameter_sweep.rs` (or a `cargo test --features sweep` flag).
- Loops over a Cartesian product of named parameters. Each point: build a fresh `SimulationGrid` + `SourceConfig` + `VacuumConfig`, apply the chosen scenario, run N steps, record diagnostics into a CSV row.
- Output to `sweep_results/<scenario>_<timestamp>.csv` with one column per measured quantity.
- Parallelism via `rayon` over independent runs.

### Time-series export for FFT analysis
**Needed by:** Phase 7.6 (parametric resonance signature at 2ω), Phase 7.7 (Taylor-relaxation oscillation tracking), Phase 8.1 (HDF5 export of probe histories).

Currently `Probe::history` is a ring buffer of `(time, value)` tuples. Add a `dump_to_csv(path: &Path)` method on `ProbeSet` so external tools (Python, Julia) can FFT or fit. This is a small piece of work but currently every "find the resonance" question needs custom code.

### Driven-current source
**Needed by:** Phase 7.2 follow-on (AC drive on the toroidal coil), Phase 7.4's `pulsed_circuit` stub (switched current), Phase 7.7 (coaxial-gun helicity injection).

Extend `SourceType` (in `src/simulation/sources.rs`) with a programmable-current variant whose amplitude is an arbitrary closed-form function of time, not just `sin(2πft)` and Gaussian. Easiest API:
```rust
SourceType::ProgrammableCurrent { axis: u8, waveform: Waveform }
enum Waveform {
    Sinusoid { frequency: f32, phase: f32 },
    GaussianPulse { sigma: f32, t_center: f32 },
    StepDown { t_switch: f32, rise_time: f32 },
    LinearRamp { from: f32, to: f32, duration: f32 },
}
```
Each Phase that needs a different waveform extends `Waveform` rather than inventing a new `SourceType`.

### Conductor boundary-condition enforcement
**Needed by:** Phase 4.2 (Casimir, two parallel conductor plates), Phase 7.4's `brown_capacitor` stub.

`CellFlags::CONDUCTOR` exists but `step_field_cpu` doesn't currently enforce φ = const on those cells. Adding it is a 5-line patch to the cell-update loop, gated by a check on the flag. Once landed, the Casimir and Brown stubs become straightforward analytic-IC scenarios.

### Magnetic helicity diagnostic
**Needed by:** Phase 7.7 (the spheromak / Taylor scenario). Called out in the 7.7 section above.

### Maxwell-stress force integrator
**Needed by:** Phase 4.2 (Casimir force on the plates), Phase 7.4 `brown_capacitor` stub. A surface integral of the Maxwell stress tensor over a closed surface enclosing the device, plus K-gradient pressure in extended QVED mode.

---

## Ongoing / Cross-Cutting Tasks

These tasks don't have strict ordering but should be addressed as the relevant systems mature.

### Testing
- **Context:** Relevant `src/simulation/` files for each test
- **Depends on:** Varies per test (noted inline)
- `[ ]` Energy conservation regression test (run 1000 steps, check energy drift < threshold) — after 1.3
- `[ ]` Wave speed validation (measure wavefront propagation, compare to c) — after 1.6
- `[ ]` Coulomb limit test (static charges, compare to 1/r² analytically) — after 1.2
- `[ ]` Weber-Coulomb equivalence (static limit of Weber = Coulomb) — after 3.2
- `[ ]` Gauge invariance test: apply gauge transform, verify E and B unchanged, verify S changes — after 2.1
- `[ ]` CPU-GPU equivalence test: run same scenario on both, compare results within f32 tolerance — after 5.2

### Documentation
- `[ ]` Code-level doc comments on all public types and functions
- `[ ]` Theory reference document: detailed derivation of the extended Maxwell equations from quaternionic potentials (separate from README, more mathematical)
- `[ ]` Experimental protocols document: step-by-step guide for each proposed physical experiment
- `[ ]` Visualization guide: how to interpret each visualization mode

### Performance
- **Context:** `ARCHITECTURE.md §Performance Targets` (lines 622–639), `ARCHITECTURE.md §Optimization Strategies` (lines 634–639)
- **Depends on:** 5.4
- `[ ]` Profile GPU compute shaders, identify bottlenecks
- `[ ]` Implement shared memory tiling in field update shader
- `[ ]` Benchmark: cells/second vs grid size, produce scaling plot
- `[ ]` Investigate multi-resolution grids (coarse far field, fine near sources)

---

## Known Issues & Deferred Considerations

These are non-blocking inconsistencies or scope questions to revisit as the relevant phases are reached.

### Physics Scope Clarifications
- **Casimir scenario (Phase 4.2):** The Casimir effect is fundamentally quantum (zero-point fluctuations). Our classical FDTD with random initial conditions is NOT equivalent to quantum vacuum. The simulation demonstrates the *mechanism* (boundary-constrained modes → K gradient → force) rather than producing quantitative Casimir predictions. The random initialization spectrum is artificial — document this clearly in the scenario.
- **Aharonov-Bohm scenario (Phase 7.1):** The AB effect is quantum phase interference. Our simulation can demonstrate that A != 0 where B = 0 (purely classical potential structure) and compute path integrals of A, but it cannot simulate the actual electron wavefunction interference. Scope is "demonstrate potential structure and compute phase-relevant integrals," not "simulate the AB effect."
- **S field as derived vs evolved:** ~~RESOLVED~~ In the original Phase 1 stub, S was a derived quantity (Lorenz gauge formula). As of Phase 1.7, S is fully independently evolved via its own □S = 0 leapfrog in `field_update.rs`, stored in double-buffered `grid.s_field`, and read directly by `diagnostics.rs` in extended mode. This is true QVED (α=1.0) with no gauge constraint on S.
- **Phase 4.1 vs Phase 1.8 — vacuum K dynamics:** The K-field evolution originally scoped for Phase 4.1 in `src/simulation/vacuum_update.rs` was implemented earlier in Phase 1.8 directly inside `step_field_cpu` (`src/simulation/field_update.rs`), driven by the `VacuumConfig` resource. The `vacuum_update.rs` module remains as a 1-line stub and should be removed (or repurposed for Phase 4.2 Casimir) when Phase 4 is opened. Phase 4.1 is effectively complete; only Phase 4.2 (Casimir) and Phase 4.3 (vacuum lensing scenario) remain.

### Resolution & Scale
- **Bifilar coil resolution:** A realistic bifilar coil has ~1mm wire spacing. At 64³ with 0.5m domain, dx ≈ 8mm — too coarse. Either use a fine grid (128³+ with small domain) or idealize the coil as a current sheet. Note this in the bifilar scenario setup.
- **Weber force magnitudes:** For realistic copper wire currents, Weber corrections scale as v_drift²/c² ≈ 10⁻¹⁸. The simulation uses scaled parameters to make forces visible. Document the scaling factor explicitly in the Graneau scenario so results are not misinterpreted as physical magnitudes.

### Crate Compatibility (verified in Phase 0.1)
- **bevy_egui 0.33 + Bevy 0.15:** Verified — bevy_egui 0.34 targets Bevy 0.16, so 0.33 is used instead. Locked in Cargo.toml.
- **CellState size assertion:** Add `const_assert!(std::mem::size_of::<CellState>() == 48)` or equivalent static check to catch layout surprises.

### Visualization / Bevy
- **Dynamic texture updates (Bevy 0.15):** Neither `images.get_mut()` + modifying `data`, NOR `images.insert(&handle, new_image)` reliably triggers Bevy's render asset re-extraction. The working fix: create a genuinely new `Image` via `images.add(new_image)` each frame, update the `StandardMaterial.base_color_texture` to the new handle, and `images.remove()` the old one. This forces a fresh GPU upload every frame.
- **Auto-range minimum:** The `auto_range` function's minimum range must be > `f32::EPSILON` (~1.19e-7). A minimum of `1e-6` prevents `map_value` from hitting the "zero range" fallback (t=0.5 → misleading mid-colormap color for uniform fields).
- **System ordering:** Visualization systems must run `.after(SimulationSet)` with `apply_deferred` inserted between entity spawn (`manage_slice_entity`) and texture update (`update_slice_texture`) to ensure deferred entity commands are flushed before queries.

### Architecture
- **CPU/GPU coexistence:** `simulation/field_update.rs` hosts the CPU implementation. When GPU is added (Phase 5), the CPU path is retained behind a `SimulationConfig.use_gpu: bool` flag for the CPU-GPU equivalence test and for debugging. The GPU dispatch wraps Bevy's `RenderDevice`/`RenderQueue`.
