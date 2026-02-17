// CellState, DerivedFields, CellFlags, SimParams definitions
// See ARCHITECTURE.md §Data Model (lines 80–137)

use bytemuck::{Pod, Zeroable};

/// Physical state at a single grid point.
///
/// Layout matches ARCHITECTURE.md exactly: 48 bytes, GPU-aligned.
/// Q = (phi/c, Ax, Ay, Az) is the quaternionic four-potential.
/// q_dot stores dQ/dt at half-step offsets for Störmer-Verlet integration.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CellState {
    /// Quaternionic potential: Q = (phi/c, Ax, Ay, Az)
    pub q: [f32; 4],

    /// Time derivative of Q (stored at half-step for leapfrog)
    pub q_dot: [f32; 4],

    /// Vacuum polarizability index K(x,t). K=1.0 is standard vacuum.
    pub k: f32,

    /// Time derivative of K
    pub k_dot: f32,

    /// Material/boundary flags (see CellFlags)
    pub flags: u32,

    /// Padding for GPU vec4 alignment
    pub _pad: f32,
}

// SAFETY: CellState is #[repr(C)] with only f32 and u32 fields, no padding gaps,
// no invalid bit patterns. Total = 4*4 + 4*4 + 4 + 4 + 4 + 4 = 48 bytes.
unsafe impl Zeroable for CellState {}
unsafe impl Pod for CellState {}

// Compile-time size check per ARCHITECTURE.md / Known Issues §CellState size assertion
const _: () = assert!(std::mem::size_of::<CellState>() == 48);

impl CellState {
    /// Create a default vacuum cell (all fields zero, K=1, no flags).
    #[inline]
    pub fn vacuum() -> Self {
        Self {
            q: [0.0; 4],
            q_dot: [0.0; 4],
            k: 1.0,
            k_dot: 0.0,
            flags: 0,
            _pad: 0.0,
        }
    }
}

impl Default for CellState {
    fn default() -> Self {
        Self::vacuum()
    }
}

/// Derived electromagnetic fields at a grid point, computed on demand for visualization.
/// Not stored in the simulation grid — recomputed each frame from CellState.
#[repr(C)]
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct DerivedFields {
    /// Electric field: E = -grad(phi) - dA/dt = -c*grad(Q.w) - Q_dot.vector()
    pub e: [f32; 3],

    /// Magnetic field: B = curl(A) = curl(Q.vector())
    pub b: [f32; 3],

    /// Scalar field: S = (1/c²)(dphi/dt) + div(A) = Q_dot.w/c + div(Q.vector())
    /// Zero in standard EM (Lorenz gauge); nonzero in QVED extended mode.
    pub s: f32,

    /// Energy density: 0.5*ε₀*|E|² + 0.5*(1/μ₀)*|B|² + 0.5*ε₀*S²
    pub energy_density: f32,

    /// Poynting vector (energy flow): (1/μ₀) * E × B
    pub poynting: [f32; 3],

    /// Local speed of light: c_local = c0 / K
    pub c_local: f32,
}

// SAFETY: DerivedFields is #[repr(C)] with only f32 fields, no padding issues.
unsafe impl Zeroable for DerivedFields {}
unsafe impl Pod for DerivedFields {}

/// Bitflags for cell material/boundary classification.
/// Stored as u32 in CellState.flags.
pub struct CellFlags;

impl CellFlags {
    pub const EMPTY: u32 = 0;
    pub const CONDUCTOR: u32 = 1 << 0;
    pub const DIELECTRIC: u32 = 1 << 1;
    pub const SOURCE: u32 = 1 << 2;
    pub const PML: u32 = 1 << 3;
    pub const BOUNDARY: u32 = 1 << 4;
}

/// Simulation parameters uniform struct, matching GPU uniform buffer layout.
/// Passed to compute shaders and used by CPU simulation loop.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SimParams {
    /// Grid dimensions
    pub nx: u32,
    pub ny: u32,
    pub nz: u32,

    /// Grid spacing (uniform, meters)
    pub dx: f32,

    /// Timestep (seconds)
    pub dt: f32,

    /// Speed of light in vacuum (m/s)
    pub c0: f32,

    /// Current simulation time (seconds)
    pub time: f32,

    /// Current iteration count
    pub iteration: u32,

    /// Extended mode flag: 0 = standard EM (Lorenz gauge), 1 = QVED extended
    pub extended_mode: u32,

    /// Padding to align to 16 bytes (GPU uniform alignment)
    pub _pad: [u32; 3],
}

unsafe impl Zeroable for SimParams {}
unsafe impl Pod for SimParams {}

impl SimParams {
    /// Speed of light in vacuum (m/s).
    pub const C0: f32 = 2.998e8;

    /// Create SimParams from grid configuration.
    pub fn new(nx: u32, ny: u32, nz: u32, dx: f32, dt: f32, extended_mode: bool) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            c0: Self::C0,
            time: 0.0,
            iteration: 0,
            extended_mode: if extended_mode { 1 } else { 0 },
            _pad: [0; 3],
        }
    }
}

/// CPML (Convolutional PML) configuration parameters.
///
/// Controls the absorbing boundary layer that replaces Mur+sponge.
/// Default values provide ~-60 dB reflection for typical simulations.
#[derive(Clone, Debug)]
pub struct PmlConfig {
    /// Number of PML cells on each Open face (default 10).
    pub depth: usize,
    /// Polynomial grading order for sigma profile (default 3.0).
    pub grading_order: f32,
    /// Kappa maximum for coordinate stretching (default 1.0 = no stretching).
    pub kappa_max: f32,
    /// Alpha maximum for CFS-PML (default 0.0). Stabilizes low-frequency damping.
    pub alpha_max: f32,
}

impl Default for PmlConfig {
    fn default() -> Self {
        Self {
            depth: 10,
            grading_order: 3.0,
            kappa_max: 1.0,
            alpha_max: 0.0,
        }
    }
}

impl PmlConfig {
    /// Compute the optimal sigma_max for the scalar wave equation PML.
    ///
    /// For d²Q/dt² = c²∇²Q with CFS-PML, the optimal conductivity is:
    ///   sigma_max = (p+1) * c / (2 * d_pml) * (-ln R)
    /// where d_pml = depth * dx is the physical PML thickness and R is
    /// the target reflection coefficient (1e-6 ≈ -120 dB).
    ///
    /// This ensures sigma * dt is O(1) so the CPML exponential decay
    /// b = exp(-(sigma/kappa + alpha) * dt) provides real absorption.
    pub fn sigma_max(&self, dx: f32) -> f32 {
        let c0 = SimParams::C0;
        let d_pml = self.depth as f32 * dx;
        let r_target: f32 = 1e-6;
        (self.grading_order + 1.0) * c0 * (-r_target.ln()) / (2.0 * d_pml)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cellstate_size() {
        assert_eq!(std::mem::size_of::<CellState>(), 48);
    }

    #[test]
    fn test_cellstate_vacuum_defaults() {
        let cell = CellState::vacuum();
        assert_eq!(cell.q, [0.0; 4]);
        assert_eq!(cell.q_dot, [0.0; 4]);
        assert_eq!(cell.k, 1.0);
        assert_eq!(cell.k_dot, 0.0);
        assert_eq!(cell.flags, CellFlags::EMPTY);
    }

    #[test]
    fn test_cellstate_bytemuck() {
        let cell = CellState::vacuum();
        let bytes: &[u8] = bytemuck::bytes_of(&cell);
        assert_eq!(bytes.len(), 48);
        let roundtrip: &CellState = bytemuck::from_bytes(bytes);
        assert_eq!(*roundtrip, cell);
    }

    #[test]
    fn test_cell_flags_distinct() {
        // All flags should be distinct powers of 2
        let flags = [
            CellFlags::CONDUCTOR,
            CellFlags::DIELECTRIC,
            CellFlags::SOURCE,
            CellFlags::PML,
            CellFlags::BOUNDARY,
        ];
        for i in 0..flags.len() {
            for j in (i + 1)..flags.len() {
                assert_eq!(flags[i] & flags[j], 0, "flags {} and {} overlap", i, j);
            }
        }
    }

    #[test]
    fn test_simparams_new() {
        let p = SimParams::new(64, 64, 64, 0.01, 1e-12, true);
        assert_eq!(p.nx, 64);
        assert_eq!(p.extended_mode, 1);
    }
}
