# CLAUDE.md — Session Instructions for Claude Code

## Project Identity

**Quaternions** — a real-time 3D simulator for Quaternionic Vacuum Electrodynamics (QVED), built in Rust with Bevy + wgpu.

The project explores extended electrodynamics: quaternionic potentials without gauge fixing, Weber forces, and polarizable vacuum models. The physics is unconventional but the engineering is rigorous.

## Session Workflow

### Step 1: Find your task

Open `TODO.md` and find the first task with status `[ ]` (not started). Each task has:

- **Context:** — the specific files and sections to read before starting. Read ONLY these.
- **Depends on:** — prerequisite tasks that must be `[x]` or `[~]` before this task can begin.

If the dependencies aren't met, skip to the next `[ ]` task.

### Step 2: Read only what the task specifies

Do NOT read the entire ARCHITECTURE.md or README.md upfront. Each task's **Context** line tells you exactly which sections are relevant. This keeps context focused and prevents wasting tokens on irrelevant material.

- `ARCHITECTURE.md §Key Rust Crates` means read only that section
- `ARCHITECTURE.md §Module Structure` means read only the directory tree
- `src/math/quaternion.rs` means read the existing file

### Step 3: Execute the task

Follow the conventions in TODO.md's header:

1. **Stub-first:** Create files with correct signatures, types, and `todo!()` bodies with dense `// PSEUDOCODE:` comments
2. **Compile-check:** Session must end with `cargo build` succeeding
3. **One concern per session:** Don't mix unrelated work
4. **Test what the task says to test:** Run `cargo test` when the task includes test items

### Step 4: Update status and documentation

After completing a task:

1. **Update TODO.md status:**
   - `[~]` if stubbed with pseudocode but not fully implemented
   - `[x]` if fully implemented and tested

2. **Update project documents if the task produced new knowledge:**
   - `ARCHITECTURE.md` — fix any version numbers, API signatures, or technical details that turned out differently than planned (e.g., crate version incompatibilities discovered during implementation)
   - `TODO.md §Known Issues` — resolve or annotate any known issues addressed by the task
   - `CLAUDE.md` — if a new project rule or convention was discovered, add it to the relevant section
   - Do NOT update documents speculatively — only when the task produced a concrete correction or discovery

## Project Rules

### Dependencies
- **NO standalone wgpu crate** — use Bevy's internal `RenderDevice`/`RenderQueue`
- **NO nalgebra** — custom `Quat` type and `[f32; 3]` cover all math
- **HDF5 is optional** — behind `--features hdf5-export` flag
- Pin bevy to 0.15; bevy_egui pinned to **0.33** (0.34+ targets Bevy 0.16)

### Code Style
- `#[repr(C)]` on all GPU-shared structs
- `bytemuck::Pod + Zeroable` on all GPU-transferable types
- Leapfrog integration uses **staggered half-step velocity** (Störmer-Verlet), not symplectic Euler
- Double-buffering: `cells: [Vec<CellState>; 2]` with `current: usize` swap index
- `extended_mode` flag controls standard vs QVED physics — never hardcode one mode

### Naming
- The project is called **quaternions** (not "qved" or "QVED")
- "QVED" refers only to the physics theory/mode, not the project name

### Scenarios

These conventions emerged from Phases 3–7 and apply to every new scenario.

- **Code units, not SI.** Analytic-IC scenarios (Aharonov-Bohm, Toroidal AB, Hopfion, Hairpin, Graneau wire) use code units with μ₀ = 1, ε₀ = 1, dropping physical prefactors. The simulator's wave-equation kernel is unit-agnostic; trying to mix SI charges with code-units potentials produces silent f32 subnormal underflow (see `memory/f32_subnormal_pitfall.md`). Scaling factor explanations live in the scenario docstring.
- **Apply functions are static loaders, not drivers.** `apply_*_scenario(grid, sources, config)` resets the grid, writes the analytic Q field cell-by-cell, and clears sources. Driven (time-varying) scenarios use the existing `Source` infrastructure, NOT extra simulation steps inside the apply function.
- **Tier-3 / speculative stubs.** Scenarios in `src/scenarios/{brown_capacitor, pulsed_circuit, charge_cluster}.rs` are compiling pseudocode-only stubs: a `*Config` struct with `Default`, an `apply_*_scenario` function with `todo!()` body, and a module-level docstring covering physics background, expected observables, parameter ranges, and concrete implementation prerequisites. Future sessions implement them without re-deriving from the README.
- **Wire to the UI selector.** Every active scenario gets a `Scenario::*` variant in `src/scenarios/dipole_radiation.rs` and a UI handler arm in `src/ui/plugin.rs::ui_side_panel`. The handler calls `apply_*_scenario`, resets PML state, sets `extended_mode` and `paused` appropriately, and clears any other resources the scenario doesn't use (`particle_system`, `probes`).

### Testing Patterns

- **`test_*_extended_matches_standard`** for any scenario with gauge-clean initial conditions (φ = 0, ∇·A = 0 in the bulk). Run the same setup on two independent grids — one in standard mode, one in extended QVED — for 30–50 steps. Assert max relative difference in A is < 1e-4. This is the canonical "QVED reduces to conventional EM for gauge-clean configs" verification and ships with AB, Toroidal AB, and (extended to) Hopfion.
- **Topology / line-integral / flux tests should compare magnitudes**, not signed values, unless the sign convention is explicitly documented. Sign depends on path orientation, current-direction parametrisation, and right-hand-rule choices that easily flip between equivalent formulations. The Toroidal AB scenario hit this — the linking flux came out negative because of how the poloidal loops were parametrised, even though the magnitude matched the analytic formula.
- **Loose thresholds for topological diagnostics under linear wave evolution.** Hedgehog Skyrmion configurations are NOT steady states of `□Q = 0`; they disperse over a few wave-traversal times. The Phase 1.9 `test_topo_charge_conserved` uses `drift < max(|initial|, 0.5)` as the criterion. Phase 7.5's hedgehog test follows the same precedent. The headline scientific claim ("K-coupling sustains topology") is the **deferred parameter sweep** over `VacuumConfig.eta`, not the linear-wave-equation conservation test.

### Modifying Core Simulation Modules

Scenarios in Phases 7.6+ (K-cycle resonator) and 7.7 (Spheromak / Taylor with K) require modifying `step_field_cpu` (`src/simulation/field_update.rs`) and `diagnostics.rs` (adding `compute_magnetic_helicity`). The pattern:

1. **Add new behaviour behind a config flag**, default off. New `*Config` resource fields, NOT new mutually-exclusive modes. The K-drive term lives next to the existing `omega_p` / `eta` in `VacuumConfig`. The magnetic-helicity diagnostic lives in `DiagnosticsState` next to `topological_charge`.
2. **Existing tests must continue to pass without modification.** A new feature with `default = false` is invisible to existing scenarios.
3. **New tests gate-enable the feature.** Tests construct a `*Config` with the new flag/field set explicitly; assert the new behaviour. Both with-flag and without-flag execution paths get coverage.
4. **Document the modification in the module header** (a 3–5 line note explaining what was added and which scenario uses it).

### Build Commands
```bash
cargo build                              # Debug build
cargo build --release                    # Release build
cargo test                               # Run all tests
cargo run --release                      # Run with default scenario
cargo build --release --features hdf5-export  # With HDF5 support
```

## File Map

| File | Purpose | When to Read |
|------|---------|-------------|
| `README.md` | Physics theory and motivation | Only if task involves physics equations or scenario design |
| `ARCHITECTURE.md` | Technical implementation guide | Only the sections your task's **Context** line specifies |
| `TODO.md` | Task breakdown with session instructions | Always — find your task here |
| `src/` | Implementation source | Read existing files when your task depends on them |
