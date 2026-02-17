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
- `[ ]` Implement `src/scenarios/bifilar_coil.rs`:
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
- `[ ]` Add to scenario selector
- `[ ]` Verify: visual confirmation of B cancellation and A/S field propagation
- **Session output**: The signature experiment for scalar wave detection, simulated

### 2.4 — Quantitative scalar wave analysis
- **Context:** `src/simulation/diagnostics.rs` (to extend with probes), `src/scenarios/bifilar_coil.rs` (scenario to measure), `src/ui/plugin.rs` (UI integration)
- **Depends on:** 2.3, 1.3
- `[ ]` Add measurement probes to `diagnostics.rs`:
  - `Probe` struct: position, records time-series of selected field quantities
  - `ProbeSet` resource: `Vec<Probe>`
  - Bevy system: each step, sample field at probe locations, append to time series
  - `fn probe_fft(probe: &Probe, field: FieldQuantity) -> Vec<(f32, f32)>` — frequency spectrum
- `[ ]` Add probe placement to UI (click on slice to place probe)
- `[ ]` Add time-series plot in UI (simple egui plot of probe data vs time)
- `[ ]` For bifilar scenario:
  - Place probes at transmitter, midpoint, and receiver
  - Plot S vs time at each location
  - Measure propagation delay → compute S-wave speed
  - Compare with c (should be equal in uniform vacuum)
- **Session output**: Quantitative evidence of scalar wave propagation speed and amplitude

---

## Phase 3: Weber Forces

### 3.1 — Discrete particle system
- **Context:** `src/simulation/grid.rs` (grid interpolation for fields at particle positions), `src/simulation/diagnostics.rs` (DerivedFields for E, B at grid points), `ARCHITECTURE.md §Simulation Loop step 6` (lines 256–261, particle push overview)
- **Depends on:** 1.3
- `[ ]` Implement `src/simulation/particles.rs`:
  - `ChargedParticle` Bevy component: charge (q), mass (m), position (Vec3), velocity (Vec3)
  - `ParticleSystem` resource: collection of particles, integration method
  - `fn spawn_particle(commands, charge, mass, position, velocity)` — creates ECS entity
  - `fn update_particles_lorentz(particles, grid)`:
    - For each particle: interpolate E, B at particle position from grid
    - Lorentz force: F = q(E + v × B)
    - Boris pusher for velocity update (standard stable EM particle pusher)
    - Position update: x += v * dt
    - `// PSEUDOCODE: Boris pusher algorithm:`
    - `// 1. Half-step electric acceleration: v⁻ = v + (q*dt/2m)*E`
    - `// 2. Rotation by B: t = (q*dt/2m)*B; s = 2t/(1+|t|²)`
    - `//    v' = v⁻ + v⁻×t; v⁺ = v⁻ + v'×s`
    - `// 3. Second half-step electric: v_new = v⁺ + (q*dt/2m)*E`
    - `// 4. Position: x_new = x + v_new * dt`
- `[ ]` Render particles as colored spheres in 3D scene
- `[ ]` Verify: single particle in uniform E field accelerates linearly; particle in uniform B field gyrates
- **Session output**: Particle system with standard Lorentz force, visualized

### 3.2 — Weber force implementation
- **Context:** `src/simulation/particles.rs` (particle system to extend), `README.md §Weber Electrodynamics` (lines 96–140, Weber force theory derivation)
- **Depends on:** 3.1
- **Note:** See TODO §Known Issues — Weber acceleration bootstrapping uses previous-step acceleration (explicit scheme).
- `[ ]` Implement `src/simulation/weber.rs`:
  - `fn weber_force(q1, q2, r12, v12, a12) -> Vec3`:
    - `r = |r12|; r_hat = r12 / r`
    - `r_dot = dot(v12, r_hat)` — relative radial velocity
    - `r_ddot = dot(a12, r_hat) + (|v12|² - r_dot²) / r` — relative radial acceleration
    - `F = (q1*q2 / 4πε₀r²) * r_hat * [1 - r_dot²/(2c²) + r*r_ddot/c²]`
    - `// PSEUDOCODE: The three terms:`
    - `// Term 1: [1] — Coulomb repulsion/attraction`
    - `// Term 2: [-r_dot²/(2c²)] — velocity correction (always reduces force)`
    - `// Term 3: [+r*r_ddot/c²] — acceleration correction (can attract or repel)`
    - `// The acceleration term is what produces longitudinal forces along current`
    - `// For steady current in a wire, the drift velocity creates a systematic`
    - `// r_dot between adjacent segments that produces a net longitudinal force`
    - `// NOTE: The acceleration a12 requires knowing current-step acceleration,`
    - `// which depends on the force we're computing (chicken-and-egg).`
    - `// Resolution: use previous-step acceleration (explicit scheme). This is`
    - `// first-order accurate and stable for v << c. If instabilities appear,`
    - `// switch to a predictor-corrector: predict a(t+dt) from a(t), compute`
    - `// force, then correct. For the Graneau scenario (steady drift), previous-`
    - `// step acceleration is sufficient since a changes slowly.`
  - `fn update_particles_weber(particles, grid)`:
    - For each particle pair: compute Weber force (N² algorithm for now)
    - Sum forces on each particle
    - Same integration as Lorentz but with Weber force instead
    - Option to compute BOTH and show comparison
  - `ForceMode` enum: `Lorentz`, `Weber`, `Both` (shows side-by-side)
- `[ ]` Add force mode toggle to UI
- `[ ]` Unit tests:
  - Two static charges: Weber force = Coulomb force (velocity terms vanish)
  - Two charges moving apart: Weber force < Coulomb (velocity correction)
  - Verify force is always along the line connecting charges (longitudinal)
- **Session output**: Weber force calculation with tests proving correctness in known limits

### 3.3 — Graneau wire scenario
- **Context:** `src/simulation/weber.rs` (Weber force API), `src/simulation/particles.rs` (particle system), `src/scenarios/dipole_radiation.rs` (scenario pattern), `README.md §Weber Electrodynamics` (lines 96–140, Graneau experiment context)
- **Depends on:** 3.2
- **Note:** See TODO §Known Issues — Weber force magnitudes use scaled parameters. Document scaling factor.
- `[ ]` Implement `src/scenarios/graneau_wire.rs`:
  - Model a straight wire as a chain of N charged segments
  - Each segment: positive ions (stationary) + electrons (drifting with velocity v_drift)
  - Apply high pulse current → v_drift increases suddenly
  - `fn setup_graneau_wire(commands, particles, grid)`:
    - Create chain of positive charges (fixed) and negative charges (mobile) along a line
    - Set electron drift velocity corresponding to desired current I
    - `// PSEUDOCODE: For a wire of cross-section A with n electron density:`
    - `// I = n * e * v_drift * A, so v_drift = I / (n * e * A)`
    - `// For copper: n ≈ 8.5e28 /m³`
    - `// For I = 10 kA, A = 1 mm²: v_drift ≈ 0.74 m/s`
    - `// Weber correction terms scale as v_drift²/c² ≈ 6e-18 — tiny!`
    - `// But the COLLECTIVE effect of ~10²³ charge pairs can sum to measurable force`
    - `// The simulation uses scaled parameters to make forces visible`
  - `fn compute_wire_force_profile(particles) -> Vec<(f32, Vec3)>`:
    - For each segment, compute total Weber force (sum over all other segments)
    - Return force magnitude along wire axis vs position
    - Compare with Lorentz prediction (should be zero for straight wire — Lorentz gives no longitudinal force on a straight wire!)
  - Expected result: Weber predicts nonzero longitudinal force that varies along wire length; Lorentz predicts zero
- `[ ]` Visualize: wire segments colored by longitudinal force magnitude; arrow glyphs showing force direction
- `[ ]` Add to scenario selector
- **Session output**: Dramatic visual comparison — Weber predicts wire-breaking forces that Lorentz cannot explain

### 3.4 — Force visualization overlays
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields with E, B, Poynting), `src/visualization/color_maps.rs`, `src/simulation/grid.rs` (grid sampling), `ARCHITECTURE.md §Visualization Modes` (lines 396–406, glyph/arrow overview)
- **Depends on:** 1.3, 0.6
- `[ ]` Implement `src/visualization/glyphs.rs`:
  - Arrow mesh generation: cone + cylinder
  - `GlyphField` component: which vector field to display, sampling density, scale factor
  - Bevy system: sample vector field on a subgrid, spawn/update instanced arrow meshes
  - Color arrows by magnitude using selected color map
  - Fields to support: E, B, A (vector potential), Poynting, Weber force, Lorentz force
  - `// PSEUDOCODE: For each sample point (x,y,z) on the subgrid:`
  - `// 1. Read vector field value V at (x,y,z)`
  - `// 2. Compute arrow length = |V| * scale_factor`
  - `// 3. Compute arrow orientation = V / |V|`
  - `// 4. Compute arrow color = color_map(|V|, min_val, max_val)`
  - `// 5. Set instance transform: translate to (x,y,z), rotate to align with V, scale by length`
- `[ ]` Add to UI: field selector, density slider, scale slider, on/off toggle
- `[ ]` Verify: dipole scenario shows correct E field arrows radiating outward
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

### 6.1 — Volume rendering
- **Context:** `ARCHITECTURE.md §Visualization Modes` (lines 381–385, volume render description), `src/visualization/color_maps.rs` (transfer function base), `src/gpu/buffers.rs` (3D texture from derived fields)
- **Depends on:** 5.4
- `[ ]` Implement `src/visualization/volume.rs`:
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
  - Transfer function: maps field value → (color, opacity)
    - Configurable via egui curve editor or presets
    - For S field: make zero transparent, positive/negative opaque with different colors
    - For |E|: transparent at low values, opaque at high values
  - Bevy custom render pipeline using a fullscreen quad with the ray-march shader
- `[ ]` Add volume render toggle and transfer function controls to UI
- `[ ]` Verify: 3D radiation pattern visible as glowing volume
- **Session output**: Gorgeous 3D visualization of propagating electromagnetic fields

### 6.2 — Streamlines
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields for E, B vector data), `src/simulation/grid.rs` (trilinear interpolation needed), `ARCHITECTURE.md §Visualization Modes` (lines 397–400, streamline description)
- **Depends on:** 1.3
- `[ ]` Implement `src/visualization/streamlines.rs`:
  - `StreamlineConfig`: seed points, field to trace, integration steps, step size
  - `fn compute_streamlines(derived_fields, config) -> Vec<Vec<Vec3>>`:
    - From each seed point, integrate using 4th-order Runge-Kutta
    - `// PSEUDOCODE: RK4 streamline integration:`
    - `// p = seed_point`
    - `// for step in 0..max_steps:`
    - `//   k1 = field_at(p) * h`
    - `//   k2 = field_at(p + k1/2) * h`
    - `//   k3 = field_at(p + k2/2) * h`
    - `//   k4 = field_at(p + k3) * h`
    - `//   p_new = p + (k1 + 2*k2 + 2*k3 + k4) / 6`
    - `//   append p_new to streamline`
    - `//   stop if |field_at(p_new)| < threshold or p_new outside grid`
    - Uses trilinear interpolation of field values between grid points
  - Render streamlines as polyline meshes, colored by field magnitude
  - Optional: animated particles moving along streamlines
- `[ ]` Add streamline controls to UI: field selector, seed placement, density
- `[ ]` Verify: dipole scenario shows classic dipole field line pattern
- **Session output**: Electric and magnetic field lines in 3D

### 6.3 — Isosurface extraction
- **Context:** `src/simulation/diagnostics.rs` (DerivedFields scalar data), `src/simulation/grid.rs` (grid dimensions for cell iteration)
- **Depends on:** 1.3
- `[ ]` Implement `src/visualization/isosurface.rs`:
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
- `[ ]` Implement `src/scenarios/aharonov_bohm.rs`:
  - Solenoid geometry: cylindrical region with uniform B inside, zero B outside
  - Set up by initializing A as the known analytic solution:
    - Inside (r < R): A_θ = B*r/2
    - Outside (r > R): A_θ = B*R²/(2r)
  - Verify B is confined inside solenoid, A extends outside
  - Compute ∮A·dl for paths encircling the solenoid vs paths that don't
  - Visualize: B shown as volume inside solenoid; A shown as streamlines outside
  - `// PSEUDOCODE: The AB phase is Φ = (q/ℏ) ∮ A·dl`
  - `// For a path encircling the solenoid: Φ = (q/ℏ) * B * π * R²`
  - `// For a path NOT encircling: Φ = 0`
  - `// The simulation should show A ≠ 0 outside where B = 0`
  - `// This is the visual proof that potentials carry physical information`
- **Session output**: Visual proof that A extends beyond B — potentials are physical

### 7.2 — Toroidal AB circuit coupling
- **Context:** `src/scenarios/aharonov_bohm.rs` (AB pattern to extend), `src/simulation/sources.rs`, `src/simulation/diagnostics.rs` (probes for measurement)
- **Depends on:** 7.1, 2.4
- `[ ]` Implement `src/scenarios/toroidal_ab.rs`:
  - Toroidal solenoid (B confined entirely inside torus)
  - External pickup loop threaded through the torus hole
  - Drive AC current through toroid
  - Measure induced signal in pickup loop
  - In standard theory: no coupling (B=0 outside torus)
  - In QVED: A-field extends outside, may produce coupling
  - `// PSEUDOCODE: This is the macro-scale AB experiment.`
  - `// If we measure voltage in the pickup loop with NO B-field leakage,`
  - `// it proves macroscopic AB coupling via potentials.`
  - `// The simulation predicts the signal magnitude for given geometry and current.`
  - `// This directly translates to an experiment: wind a toroid, add a pickup loop,`
  - `// drive with AC, look for signal with shielded sensitive measurement.`
  - Output: predicted coupling coefficient vs frequency, geometry parameters
- `[ ]` Include experimental BOM: toroid core specs, wire gauge, oscilloscope requirements
- **Session output**: Experimentally testable prediction with component list

### 7.3 — PCB Weber force geometry optimizer
- **Context:** `src/scenarios/graneau_wire.rs` (wire scenario to extend), `src/simulation/weber.rs` (Weber force API), `src/simulation/particles.rs`
- **Depends on:** 3.3
- `[ ]` Implement `src/scenarios/graneau_wire.rs` extended version:
  - Instead of straight wire, model PCB trace geometries:
    - Hairpin turns (current doubles back on itself)
    - Spiral traces
    - Interdigitated fingers
  - For each geometry: compute Weber force map under pulsed current
  - Find the geometry that maximizes longitudinal force per unit current
  - `// PSEUDOCODE: The Weber longitudinal force between two current elements is:`
  - `// dF = (μ₀/4π) * I² * dl₁·dl₂ / r² * [correction terms]`
  - `// The correction produces a force ALONG the current direction`
  - `// that doesn't exist in the Lorentz/Biot-Savart formulation.`
  - `// For a hairpin geometry, the two parallel traces create strong`
  - `// longitudinal forces because the current elements are close and antiparallel.`
  - `// The simulation can search over trace width, spacing, and shape`
  - `// to find the geometry that maximizes the anomalous force.`
  - Output: optimal PCB geometry with predicted force magnitude vs pulse current
- `[ ]` Export optimized geometry as KiCad footprint (or at least dimensions)
- **Session output**: Ready-to-fabricate PCB design for Weber force measurement

### 7.4 — Tier 3 scenario stubs
- **Context:** `src/scenarios/mod.rs` (module structure), `README.md §Tier 3 Scenarios` (lines 253–278, theory context for each)
- **Depends on:** 0.1
- `[ ]` Implement stubs for speculative scenarios (pseudocode only, no full sim):
  - `src/scenarios/brown_capacitor.rs` — asymmetric capacitor in PV vacuum
  - `src/scenarios/pulsed_circuit.rs` — sharp current interruption, scalar mode excitation
  - `src/scenarios/charge_cluster.rs` — EVO stability search
  - Each file: scenario description comments, setup function signature, expected observables, parameter ranges to explore, relevant references
  - Full implementation deferred until framework is validated
- **Session output**: Roadmap stubs for all speculative experiments

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
- **S field as derived vs evolved:** The README states S obeys its own wave equation □S = -ρ/ε₀. In our implementation, S is NOT independently evolved — it is a *derived* quantity computed from Q at each step. The wave-like behavior of S emerges from the coupled Q evolution when extended mode is active. This is physically equivalent but numerically different from solving a separate S PDE. If numerical issues arise (S not propagating cleanly), consider adding explicit S evolution as an alternative.

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
