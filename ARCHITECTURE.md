# ARCHITECTURE.md — Technical Implementation Guide

## Quaternions — Quaternionic Vacuum Electrodynamics Simulator

---

## Design Philosophy

**Visualization-first, physics-accurate, computationally efficient.**

The simulator is designed so that every computational construct has a corresponding visual representation. The user should be able to see the quaternionic fields evolving in real time, manipulate source configurations interactively, and immediately observe how changes propagate through the vacuum. The visualization is not a post-processing step — it is the primary interface to the physics.

This means:
- The simulation loop and the render loop are tightly coupled but independently configurable (simulation can run faster or slower than rendering)
- Every physical quantity (Q, S, E, B, K, forces) has a dedicated visualization mode
- 3D spatial structure is always visible; higher-dimensional data (quaternion components, vacuum state) is encoded in color, opacity, glyphs, or selectable overlays
- Interactive parameter adjustment without stopping the simulation

---

## Technology Stack

### Core: Rust + Bevy + wgpu

**Rust** is the implementation language for the entire project. Reasons:

- Memory safety without garbage collection — critical for real-time simulation where allocation pauses cause frame drops
- Zero-cost abstractions — the quaternion math compiles to the same SIMD instructions you'd write in C
- Excellent WebGPU support through `wgpu` — the simulation kernel runs on GPU compute shaders
- The type system enforces physical correctness: quaternion types, field types, and index types are distinct and cannot be accidentally mixed
- Cargo ecosystem provides everything needed without dependency hell

**Bevy** is the application framework. It provides:

- Entity Component System (ECS) architecture — maps naturally to the simulation domain
- Real-time 3D rendering with PBR pipeline
- Camera controls, input handling, UI framework (bevy_egui)
- Built-in parallelism for systems
- Plugin architecture for modular simulation components

**wgpu** handles GPU compute. The heavy numerical work (field updates, force computations) runs in WGSL compute shaders on the GPU. **We access wgpu through Bevy's internal render resources** (`RenderDevice`, `RenderQueue`, `RenderAdapter`) rather than pulling in a standalone `wgpu` crate. This avoids version conflicts since Bevy bundles its own wgpu internally — using a separate `wgpu` dependency would cause type mismatches. All GPU buffer creation and compute dispatch goes through Bevy's render world.

### Key Rust Crates

```toml
[dependencies]
bevy = "0.15"                    # Application framework, ECS, rendering, AND wgpu access
bevy_egui = "0.33"               # Immediate-mode UI for parameter controls (0.33 targets Bevy 0.15; 0.34+ targets Bevy 0.16)
# NOTE: No standalone wgpu dependency — use Bevy's internal wgpu via RenderDevice/RenderQueue
bytemuck = "1.14"                # Zero-copy casting for GPU buffer data
serde = { version = "1", features = ["derive"] }  # Serialization for configs/saves
ron = "0.8"                      # Rusty Object Notation for config files
rayon = "1.8"                    # CPU parallelism for non-GPU tasks

[features]
default = []
hdf5-export = ["hdf5"]          # Optional: requires libhdf5 system library

[dependencies.hdf5]
version = "0.8"
optional = true                  # HDF5 is painful on Windows; gate behind feature flag
```

> **Note on nalgebra:** Removed as a direct dependency. The custom `Quat` type and `[f32; 3]` vectors cover all needed math. If matrix operations become necessary later (e.g., stress tensors), nalgebra can be re-added at that point.

### Alternative Visualization Approaches (considered and rejected)

| Option | Pros | Cons | Verdict |
|--------|------|------|---------|
| Rust + Bevy + wgpu | Type safety, GPU compute, real-time 3D, single language | Steeper initial setup | **Selected** |
| Three.js + WebGPU | Fast prototyping, browser-based | Weak numerics, no native quaternion math, JS performance ceiling | Good for demo/web viewer, not for core sim |
| Julia + GLMakie | Beautiful math notation, fast prototyping | Limited rendering control, no GPU compute shaders, GC pauses | Good for prototyping equations, not for real-time sim |
| Python + VTK/PyVista | Huge ecosystem | Far too slow for real-time; would require C++ backend anyway | No |
| C++ + Vulkan | Maximum performance | Development velocity far too slow; unsafe | No |

**Note:** A Three.js web viewer could be added later as a secondary visualization target, reading simulation data exported from the Rust core. This would allow sharing results in a browser without requiring the full Rust toolchain.

---

## Data Model

### Simulation Grid

The simulation domain is a regular 3D Cartesian grid. Each cell stores:

```rust
/// Physical state at a single grid point
#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct CellState {
    // Quaternionic potential: Q = (phi/c, Ax, Ay, Az)
    q: [f32; 4],

    // Time derivative of Q (for leapfrog integration)
    q_dot: [f32; 4],

    // Vacuum polarizability index K(x,t)
    // K = 1.0 is standard vacuum
    // K > 1.0 is polarized vacuum (higher permittivity)
    k: f32,

    // Time derivative of K
    k_dot: f32,

    // Material flags (conductor, dielectric, source, boundary)
    flags: u32,

    // Padding for GPU alignment (vec4 alignment)
    _pad: f32,
}
// Total: 48 bytes per cell, well-aligned for GPU
```

### Derived Quantities (computed on demand for visualization)

```rust
/// Derived electromagnetic fields at a grid point
struct DerivedFields {
    // Electric field: E = -grad(phi) - dA/dt
    e: [f32; 3],

    // Magnetic field: B = curl(A)
    b: [f32; 3],

    // Scalar field: S = (1/c²)(dphi/dt) + div(A)
    // This is the KEY quantity — zero in standard EM, nonzero in QVED
    s: f32,

    // Energy density: standard EM + scalar field + vacuum polarization
    energy_density: f32,

    // Poynting vector (energy flow)
    poynting: [f32; 3],

    // Local speed of light: c_local = c0 / K
    c_local: f32,
}
```

### Quaternion Type

```rust
/// Hamilton quaternion: q = w + xi + yj + zk
#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct Quat {
    w: f32,  // scalar part (phi/c for EM potential)
    x: f32,  // i component (Ax)
    y: f32,  // j component (Ay)
    z: f32,  // k component (Az)
}

impl Quat {
    /// Hamilton product (non-commutative)
    fn mul(self, rhs: Quat) -> Quat {
        Quat {
            w: self.w*rhs.w - self.x*rhs.x - self.y*rhs.y - self.z*rhs.z,
            x: self.w*rhs.x + self.x*rhs.w + self.y*rhs.z - self.z*rhs.y,
            y: self.w*rhs.y - self.x*rhs.z + self.y*rhs.w + self.z*rhs.x,
            z: self.w*rhs.z + self.x*rhs.y - self.y*rhs.x + self.z*rhs.w,
        }
    }

    /// Conjugate: q* = w - xi - yj - zk
    fn conj(self) -> Quat {
        Quat { w: self.w, x: -self.x, y: -self.y, z: -self.z }
    }

    /// Norm squared: |q|² = w² + x² + y² + z²
    fn norm_sq(self) -> f32 {
        self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z
    }

    /// Scalar part
    fn scalar(self) -> f32 { self.w }

    /// Vector part
    fn vector(self) -> [f32; 3] { [self.x, self.y, self.z] }
}
```

On the GPU (WGSL), this maps directly to `vec4<f32>` with a custom Hamilton product function.

---

## Simulation Architecture

### System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        Bevy Application                         │
│                                                                 │
│  ┌──────────────┐  ┌──────────────┐  ┌───────────────────────┐ │
│  │  UI Systems   │  │  Render      │  │  Simulation Plugin    │ │
│  │  (egui)       │  │  Systems     │  │                       │ │
│  │              │  │              │  │  ┌─────────────────┐  │ │
│  │ • Parameters │  │ • Field viz  │  │  │ GPU Compute     │  │ │
│  │ • Mode select│  │ • Glyph viz  │  │  │                 │  │ │
│  │ • Source edit │  │ • Slice viz  │  │  │ • Field update  │  │ │
│  │ • Run/pause  │  │ • Isosurface │  │  │ • Source inject  │  │ │
│  │ • Export     │  │ • Streamline │  │  │ • Boundary cond │  │ │
│  │              │  │              │  │  │ • Derived fields│  │ │
│  └──────────────┘  └──────────────┘  │  │ • Weber forces  │  │ │
│                                       │  └─────────────────┘  │ │
│                                       │                       │ │
│                                       │  ┌─────────────────┐  │ │
│                                       │  │ CPU Systems     │  │ │
│                                       │  │                 │  │ │
│                                       │  │ • Particle push │  │ │
│                                       │  │ • Diagnostics   │  │ │
│                                       │  │ • Data export   │  │ │
│                                       │  └─────────────────┘  │ │
│                                       └───────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                              │
                    ┌─────────┴─────────┐
                    │    GPU Buffers     │
                    │                   │
                    │  cell_state[]     │ ← Primary simulation state
                    │  derived_fields[] │ ← Computed for visualization
                    │  source_terms[]   │ ← Current/charge sources
                    │  particle_buf[]   │ ← Discrete charge particles
                    └───────────────────┘
```

### Simulation Loop (per timestep)

```
1. INJECT SOURCES
   - Write current density J(x,t) and charge density ρ(x,t) into source buffer
   - For pulsed circuits: update source waveforms based on circuit model

2. UPDATE QUATERNIONIC POTENTIAL (GPU compute shader)
   For each grid cell (x,y,z):
     - Compute spatial derivatives of Q using finite differences
       ∂Q/∂x, ∂Q/∂y, ∂Q/∂z, ∇²Q (Laplacian)
     - Compute vacuum-modified wave equation:
       ∂²Q/∂t² = (c₀²/K²) ∇²Q - (c₀²/K²) ∇(∇·A + K²/c₀² ∂φ/∂t) + source terms
     - NOTE: the ∇(∇·A + ...) term is NOT set to zero — this is what
       distinguishes QVED from standard FDTD EM simulation
     - Update Q and Q_dot using leapfrog (Verlet) integration

3. UPDATE VACUUM STATE (GPU compute shader)
   For each grid cell:
     - Compute local field energy density u(x) from Q
     - Update K(x) based on vacuum polarization model:
       K(x) = 1 + α · u(x) / u_critical
       where α is the vacuum polarizability coupling constant
       and u_critical is the Schwinger limit (~4.6 × 10¹⁸ V/m)
     - For gravitational sources: K(x) = exp(2GM/rc₀²)

4. APPLY BOUNDARY CONDITIONS (GPU compute shader)
   - Absorbing boundaries (PML - Perfectly Matched Layer)
   - Conducting boundaries (Q tangential components → 0)
   - Periodic boundaries (for infinite structure simulation)
   - Custom boundaries (for specific experimental geometries)

5. COMPUTE DERIVED FIELDS (GPU compute shader, only for visible cells)
   - E, B, S from Q using finite differences
   - Energy density, Poynting vector
   - For Weber force mode: compute inter-particle forces

6. PARTICLE PUSH (CPU or GPU)
   - For discrete charge simulations:
     - Compute force on each particle from local fields
     - Include Weber correction terms for particle-particle interactions
     - Update positions and velocities (Boris pusher or similar)

7. VISUALIZATION UPDATE (Bevy render systems)
   - Read derived fields from GPU buffer
   - Update visualization meshes/textures
   - Render frame
```

### GPU Compute Shader: Field Update Kernel

The core simulation kernel in WGSL (WebGPU Shading Language):

```wgsl
// Core field update compute shader (simplified)

struct CellState {
    q: vec4<f32>,      // quaternionic potential (phi/c, Ax, Ay, Az)
    q_dot: vec4<f32>,  // time derivative
    k: f32,            // vacuum polarizability
    k_dot: f32,
    flags: u32,
    _pad: f32,
};

@group(0) @binding(0) var<storage, read> cells_in: array<CellState>;
@group(0) @binding(1) var<storage, read_write> cells_out: array<CellState>;
@group(0) @binding(2) var<uniform> params: SimParams;

fn idx(x: u32, y: u32, z: u32) -> u32 {
    return x + y * params.nx + z * params.nx * params.ny;
}

// Hamilton product of two quaternions stored as vec4
fn qmul(a: vec4<f32>, b: vec4<f32>) -> vec4<f32> {
    return vec4<f32>(
        a.x*b.x - a.y*b.y - a.z*b.z - a.w*b.w,
        a.x*b.y + a.y*b.x + a.z*b.w - a.w*b.z,
        a.x*b.z - a.y*b.w + a.z*b.x + a.w*b.y,
        a.x*b.w + a.y*b.z - a.z*b.y + a.w*b.x,
    );
}

@compute @workgroup_size(8, 8, 4)
fn update_field(@builtin(global_invocation_id) gid: vec3<u32>) {
    let x = gid.x; let y = gid.y; let z = gid.z;
    if (x >= params.nx || y >= params.ny || z >= params.nz) { return; }

    let i = idx(x, y, z);
    let cell = cells_in[i];
    if ((cell.flags & BOUNDARY_FLAG) != 0u) { return; }

    let dx = params.dx;
    let dt = params.dt;
    let c0 = params.c0;
    let K = cell.k;
    let c_local = c0 / K;

    // Fetch neighbors for finite differences
    let qxp = cells_in[idx(x+1, y, z)].q;
    let qxm = cells_in[idx(x-1, y, z)].q;
    let qyp = cells_in[idx(x, y+1, z)].q;
    let qym = cells_in[idx(x, y-1, z)].q;
    let qzp = cells_in[idx(x, y, z+1)].q;
    let qzm = cells_in[idx(x, y, z-1)].q;

    // Laplacian of Q (component-wise, standard 7-point stencil)
    let lap_q = (qxp + qxm + qyp + qym + qzp + qzm - 6.0 * cell.q) / (dx * dx);

    // Divergence of A (vector part of Q): dAx/dx + dAy/dy + dAz/dz
    let div_a = (qxp.y - qxm.y + qyp.z - qym.z + qzp.w - qzm.w) / (2.0 * dx);

    // d(phi/c)/dt from q_dot
    let dphi_dt_over_c = cell.q_dot.x;

    // Scalar field: S = (1/c²)(dphi/dt) + div(A) = dphi_dt_over_c / c + div_a
    let S = dphi_dt_over_c / c_local + div_a;

    // Gradient of S
    let Sxp = compute_S_at(idx(x+1, y, z));
    let Sxm = compute_S_at(idx(x-1, y, z));
    let Syp = compute_S_at(idx(x, y+1, z));
    let Sym = compute_S_at(idx(x, y-1, z));
    let Szp = compute_S_at(idx(x, y, z+1));
    let Szm = compute_S_at(idx(x, y, z-1));
    let grad_S = vec3<f32>(
        (Sxp - Sxm) / (2.0 * dx),
        (Syp - Sym) / (2.0 * dx),
        (Szp - Szm) / (2.0 * dx),
    );

    // Extended wave equation for Q:
    // d²Q/dt² = c_local² * lap_Q + source_terms                    (STANDARD MODE)
    // d²Q/dt² = c_local² * lap_Q - c_local² * (0, grad_S) + src   (EXTENDED MODE)
    //
    // The -c_local² * grad_S coupling in the vector components is what allows
    // the scalar field S to propagate as a genuine longitudinal wave.
    // In standard mode, this term is OMITTED, which implicitly enforces
    // Lorenz gauge and restricts solutions to transverse modes only.
    // In extended mode, this term is INCLUDED, coupling S back into
    // the potential evolution and enabling scalar/longitudinal dynamics.
    let scalar_coupling = vec4<f32>(0.0, grad_S.x, grad_S.y, grad_S.z);

    var q_ddot = c_local * c_local * lap_q;
    if (params.extended_mode != 0u) {
        // EXTENDED QVED: include scalar coupling → longitudinal modes propagate
        q_ddot = q_ddot - c_local * c_local * scalar_coupling;
    }
    // STANDARD: omit scalar_coupling → Lorenz gauge enforced implicitly

    // Leapfrog (Störmer-Verlet) integration — velocity is staggered at half-steps
    // q_dot is stored at t - dt/2; after update it represents t + dt/2
    // q is stored at t; after update it represents t + dt
    // This staggering is critical for energy conservation over long runs.
    // Initialization must set q_dot at t = -dt/2 (half-step offset).
    var new_q_dot = cell.q_dot + q_ddot * dt;   // q_dot(t+dt/2) = q_dot(t-dt/2) + a(t)*dt
    var new_q = cell.q + new_q_dot * dt;         // q(t+dt) = q(t) + q_dot(t+dt/2)*dt

    cells_out[i].q = new_q;
    cells_out[i].q_dot = new_q_dot;
    cells_out[i].k = cell.k; // Updated separately in vacuum update pass
}
```

This is a simplified version — the production shader would include PML boundary handling, source injection, and vacuum coupling terms. The `params.extended_mode` flag (a `u32` in `SimParams`) controls whether the scalar coupling is active.

> **SimParams must include:** `extended_mode: u32` (0 = standard, 1 = extended QVED). This is a uniform value set from the CPU-side `SimulationConfig.extended_mode: bool`.

---

## Visualization Architecture

### Visualization Modes

The simulator supports multiple simultaneous visualization modes, toggled and configured through the egui panel:

#### 1. Volume Rendering (Primary Mode)
- 3D scalar field displayed as semi-transparent colored volume
- Selectable quantity: |E|, |B|, S (scalar field), K (vacuum state), energy density, |Q|
- Configurable color map and opacity transfer function
- GPU ray-marching through 3D texture

#### 2. Slice Planes
- Orthogonal slice planes (XY, XZ, YZ) at adjustable positions
- Color-mapped scalar field on each slice
- Can display vector fields (E, B, Poynting) as overlaid arrows
- Multiple slices can be active simultaneously

#### 3. Vector Field Glyphs
- Arrows or cones showing field direction and magnitude at sampled points
- Selectable field: E, B, A (vector potential), Poynting vector, force
- Adjustable sampling density and scale

#### 4. Streamlines
- 3D curves following field lines of E, B, A, or Poynting vector
- Seeded from user-placed points or automatically from sources
- Animated particles flowing along streamlines to show direction and magnitude

#### 5. Isosurfaces
- 3D surfaces of constant value for any scalar quantity
- Adjustable threshold
- Useful for visualizing wavefronts (isosurface of |E| or S) and vacuum structure (isosurface of K)

#### 6. Quaternion Visualization
- The full quaternion Q has 4 components at each point — more than can be shown with a single 3D visualization
- Approach: display the vector part (A) as arrows and encode the scalar part (φ) as color on the arrows or as a background volume
- Alternative: separate the quaternion into modulus |Q| (volume render) and "phase" (encoded as hue)
- The Hamilton product structure can be visualized by showing how Q rotates a reference quaternion at each point — this gives a rotation field that can be rendered as oriented ellipsoids

#### 7. Source/Geometry Display
- Conductors shown as solid metallic surfaces
- Current paths shown as glowing tubes with animated flow direction
- Charge distributions shown as colored spheres (+ red, - blue)
- Coil geometries rendered with wire-frame detail
- PCB traces shown as flat conductors on board surface

### Rendering Pipeline

```
Simulation GPU Buffers
        │
        ▼
┌────────────────────┐
│ Derived Field Pass  │  GPU compute: Q → E, B, S, energy, K
│ (compute shader)    │  Writes to 3D textures for visualization
└────────────────────┘
        │
        ▼
┌────────────────────┐
│ Visualization Pass  │  Bevy render systems
│                     │
│  ├─ Volume Render   │  Ray-march through 3D texture
│  ├─ Slice Render    │  Textured quads with field data
│  ├─ Glyph Render    │  Instanced arrow meshes
│  ├─ Streamlines     │  GPU-computed curves → line meshes
│  ├─ Isosurfaces     │  Marching cubes on GPU → triangle mesh
│  └─ Geometry        │  Standard Bevy mesh rendering
└────────────────────┘
        │
        ▼
┌────────────────────┐
│ UI Overlay          │  bevy_egui panel
│                     │
│  ├─ Parameter panel │  Sliders for all physical parameters
│  ├─ Vis controls    │  Mode selection, color maps, thresholds
│  ├─ Source editor   │  Place/configure sources interactively
│  ├─ Diagnostics     │  Total energy, max S, divergence checks
│  └─ Export controls │  Screenshot, data export, animation record
└────────────────────┘
```

---

## Module Structure

```
quaternions/
├── Cargo.toml
├── README.md
├── ARCHITECTURE.md
├── src/
│   ├── main.rs                     # Bevy app setup, plugin registration
│   ├── lib.rs                      # Re-exports, top-level types
│   │
│   ├── math/
│   │   ├── mod.rs
│   │   ├── quaternion.rs           # Quat type with Hamilton algebra
│   │   ├── vector_field.rs         # 3D vector field operations
│   │   └── fdtd.rs                 # Finite difference helpers (stencils, etc.)
│   │
│   ├── simulation/
│   │   ├── mod.rs
│   │   ├── plugin.rs               # Bevy plugin: registers all sim systems
│   │   ├── grid.rs                 # Grid resource: dimensions, spacing, double-buffered cells
│   │   ├── state.rs                # CellState, DerivedFields definitions
│   │   ├── field_update.rs         # Field evolution (CPU impl first, GPU dispatch added later)
│   │   ├── vacuum_update.rs        # GPU compute for vacuum polarization K(x,t)
│   │   ├── sources.rs              # Source injection: currents, charges, waveforms
│   │   ├── boundaries.rs           # PML, conducting, periodic boundaries
│   │   ├── weber.rs                # Weber force computation for discrete charges
│   │   ├── particles.rs            # Charged particle system (ECS entities + pusher)
│   │   └── diagnostics.rs          # Energy conservation, divergence checks
│   │
│   ├── gpu/
│   │   ├── mod.rs
│   │   ├── pipeline.rs             # Compute pipeline setup and management
│   │   ├── buffers.rs              # GPU buffer creation and synchronization
│   │   └── shaders/
│   │       ├── field_update.wgsl   # Core quaternionic FDTD kernel
│   │       ├── vacuum_update.wgsl  # Vacuum state evolution kernel
│   │       ├── derived_fields.wgsl # E, B, S computation from Q
│   │       ├── weber_force.wgsl    # Weber force between particles
│   │       ├── pml_boundary.wgsl   # Absorbing boundary implementation
│   │       └── source_inject.wgsl  # Source term injection
│   │
│   ├── visualization/
│   │   ├── mod.rs
│   │   ├── plugin.rs               # Bevy plugin: registers all viz systems
│   │   ├── volume.rs               # Volume rendering (ray-march shader)
│   │   ├── slices.rs               # Slice plane rendering
│   │   ├── glyphs.rs               # Vector field glyph instancing
│   │   ├── streamlines.rs          # Streamline computation and rendering
│   │   ├── isosurface.rs           # Marching cubes isosurface extraction
│   │   ├── geometry.rs             # Source/conductor geometry rendering
│   │   ├── color_maps.rs           # Scientific color maps (viridis, plasma, etc.)
│   │   └── quaternion_viz.rs       # Special quaternion visualization modes
│   │
│   ├── ui/
│   │   ├── mod.rs
│   │   ├── plugin.rs               # Bevy plugin: egui panels
│   │   ├── parameter_panel.rs      # Physical parameter controls
│   │   ├── visualization_panel.rs  # Viz mode selection and configuration
│   │   ├── source_editor.rs        # Interactive source placement
│   │   ├── diagnostics_panel.rs    # Real-time simulation diagnostics
│   │   └── export_panel.rs         # Data export and screenshot controls
│   │
│   ├── scenarios/
│   │   ├── mod.rs
│   │   ├── aharonov_bohm.rs        # Solenoid geometry, phase computation
│   │   ├── casimir.rs              # Parallel plate vacuum polarization
│   │   ├── graneau_wire.rs         # High-current conductor with Weber forces
│   │   ├── bifilar_coil.rs         # Scalar wave transmitter/receiver
│   │   ├── toroidal_ab.rs          # Toroidal coil AB-effect circuit
│   │   ├── brown_capacitor.rs      # Asymmetric capacitor in PV vacuum
│   │   ├── pulsed_circuit.rs       # Sharp current interruption scalar excitation
│   │   ├── charge_cluster.rs       # EVO/charge cluster stability search
│   │   └── dipole_radiation.rs     # Standard dipole (validation benchmark)
│   │
│   └── export/
│       ├── mod.rs
│       ├── hdf5_export.rs          # Full state export to HDF5 (behind hdf5-export feature)
│       ├── screenshot.rs           # Image capture
│       └── animation.rs            # Frame sequence recording
│
├── assets/
│   ├── shaders/                    # Bevy render shaders (volume, slice, etc.)
│   │   ├── volume_render.wgsl
│   │   ├── slice_render.wgsl
│   │   ├── glyph_render.wgsl
│   │   └── streamline_render.wgsl
│   ├── color_maps/                 # Color map textures (1D lookup)
│   │   ├── viridis.png
│   │   ├── plasma.png
│   │   ├── coolwarm.png
│   │   └── scalar_field.png        # Custom: blue(neg) → white(zero) → red(pos)
│   └── scenarios/                  # Preset scenario configurations
│       ├── aharonov_bohm.ron
│       ├── graneau_wire.ron
│       └── bifilar_coil.ron
│
└── tests/
    ├── conservation.rs             # Energy conservation tests
    ├── wave_propagation.rs         # Transverse wave speed = c validation
    ├── coulomb_limit.rs            # Static limit recovers Coulomb's law
    ├── weber_limit.rs              # Near-field limit recovers Weber force
    └── gauge_invariance.rs         # Gauge transform doesn't change E, B
                                    # (but DOES change S — this is the test)
```

---

## Implementation Phases

### Phase 0: Scaffolding (Week 1-2)
- Set up Bevy application with camera, basic 3D scene, egui panel
- Implement Quat type with full Hamilton algebra and tests
- Create GPU buffer infrastructure (grid allocation, double-buffering)
- Implement basic slice plane visualization with dummy data
- **Milestone:** Spinning quaternion field displayed on a color-mapped slice plane

### Phase 1: Standard FDTD Baseline (Week 3-5)
- Implement standard FDTD electromagnetic update on GPU (Yee algorithm equivalent for potentials)
- Implement dipole source injection
- Implement PML absorbing boundaries
- Add volume rendering for field magnitude
- Validate: dipole radiation pattern, wave speed = c, energy conservation
- **Milestone:** Beautiful real-time dipole radiation in 3D, visually verified correct

### Phase 2: Extended Electrodynamics (Week 6-8)
- Remove Lorenz gauge enforcement from field update kernel
- Add scalar field S computation and visualization
- Implement bifilar coil source geometry (B-canceling, A-surviving)
- Implement S-field specific visualization (separate color map, dedicated slice)
- Test: drive oscillating source, observe S propagation vs transverse wave propagation
- **Milestone:** First visualization of scalar longitudinal wave propagating from a bifilar source

### Phase 3: Weber Forces and Particles (Week 9-11)
- Implement discrete charged particle system (ECS entities)
- Implement Weber force computation between particle pairs
- Compare Lorentz vs Weber force predictions for current-carrying geometries
- Implement Graneau wire scenario: model wire as chain of charged segments
- Visualize force distribution along wire
- **Milestone:** Side-by-side comparison of Lorentz vs Weber force on a current-carrying conductor

### Phase 4: Polarizable Vacuum (Week 12-14)
- Implement vacuum state K(x,t) as dynamical variable
- Implement vacuum polarization response to field energy
- Add K visualization (volume render or isosurface)
- Test Casimir scenario: plates modify K between them
- Test feedback: strong field → K change → field propagation change
- **Milestone:** Casimir-like vacuum polarization profile between conducting plates

### Phase 5: Scenario Library (Week 15-18)
- Implement all Tier 1 and Tier 2 scenarios from README
- Add scenario selection UI with descriptions and preset parameters
- Implement measurement tools: point probes, line integrals, surface integrals
- Add data export for quantitative analysis
- **Milestone:** Complete interactive exploration of all non-speculative scenarios

### Phase 6: Advanced Features and Speculative Scenarios (Week 19+)
- Implement Tier 3 scenarios (Brown capacitor, pulsed circuit, charge clusters)
- Add topological charge computation
- Implement multi-resolution grids for large-scale + fine-detail simulation
- Optimize compute shaders for larger grids
- Web viewer (Three.js) for sharing results

---

## Performance Targets

| Grid Size | Cells | GPU Memory | Target Framerate | Use Case |
|-----------|-------|------------|-----------------|----------|
| 64³ | 262K | ~12 MB | 60+ FPS | Real-time interactive exploration |
| 128³ | 2.1M | ~96 MB | 30+ FPS | Detailed single-phenomenon study |
| 256³ | 16.8M | ~768 MB | 5-10 FPS | High-resolution visualization |
| 512³ | 134M | ~6 GB | Offline | Publication-quality renders |

These targets assume a modern discrete GPU (RTX 3060+ or equivalent). The 48 bytes/cell storage means even 128³ fits comfortably in VRAM with room for derived field buffers and visualization textures.

### Optimization Strategies
- **Double-buffered** state arrays (read from one, write to other, swap)
- **Lazy derived field computation**: only compute E, B, S for cells that are currently visible (slice planes, volume render sample points)
- **Adaptive timestep**: CFL condition determines maximum stable dt; scale with smallest dx and fastest c_local
- **Workgroup size tuning**: 8×8×4 = 256 threads per workgroup matches most GPU architectures
- **Shared memory tiling**: load cell neighborhoods into workgroup shared memory to reduce global memory reads

---

## Experimental Interface

The simulator is designed to bridge directly to physical experiments. Each scenario includes:

### Circuit/Geometry Specification
- Conductor paths defined as 3D curves (coils, PCB traces, wire geometries)
- Parameterized by physical dimensions (mm), wire gauge, and material properties
- Exportable as KiCad footprints for PCB scenarios

### Measurement Points
- Virtual probes placed at specific locations to record field values over time
- Equivalent to where you would place a physical sensor (B-field probe, antenna, voltmeter)
- Time-series data exportable as CSV for comparison with experimental measurements

### Predicted Observables
- For each scenario, the simulator outputs specific quantitative predictions:
  - "At this point, the S field amplitude should be X V/m for input current Y A"
  - "The Weber longitudinal force on this segment should be Z N for pulse current W A"
  - "The voltage induced in the secondary coil by A-field coupling should be V μV"
- These predictions are what the physical experiment would measure to validate or falsify the theory

### Bill of Materials Generator
- For experimental scenarios, the simulator can output a component list:
  - Coil specifications (turns, diameter, wire gauge, inductance)
  - Capacitor values and voltage ratings
  - PCB dimensions and trace widths
  - Required instrumentation (oscilloscope bandwidth, current probe sensitivity)

---

## Data Export Formats

- **HDF5**: Full simulation state snapshots for offline analysis (Python, Julia, MATLAB compatible)
- **CSV**: Time-series data from virtual probes
- **PNG/EXR**: High-resolution screenshots and frame sequences
- **RON**: Scenario configuration files (human-readable, version-controllable)
- **KiCad**: PCB trace geometries for experimental scenarios
- **STL**: 3D-printable coil forms and experimental fixtures

---

## Getting Started

### Prerequisites
- Rust toolchain (stable, 1.75+)
- GPU with Vulkan, Metal, or DX12 support
- For HDF5 export (optional): `cargo build --features hdf5-export` (requires libhdf5 system library)

### Build and Run
```bash
# Build (release mode recommended for simulation performance)
cargo build --release

# Run with default scenario (dipole radiation)
cargo run --release

# Run specific scenario
cargo run --release -- --scenario bifilar_coil

# Run with custom grid size
cargo run --release -- --grid 128 --scenario graneau_wire

# Build with HDF5 data export support (requires libhdf5 installed)
cargo build --release --features hdf5-export
```

### First Steps
1. Start with the dipole radiation scenario to verify the visualization works
2. Switch to the Aharonov-Bohm scenario to see potentials extending beyond field regions
3. Open the bifilar coil scenario and toggle the S-field visualization to see the scalar mode
4. Compare Lorentz vs Weber forces in the Graneau wire scenario
5. Enable vacuum polarization in any scenario and observe K(x) responding to field intensity