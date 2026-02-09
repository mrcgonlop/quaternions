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
