// Finite difference helpers: grid indexing, stencils, CFL condition
// Provides the low-level grid infrastructure used by all spatial derivative operations.

// --- Grid indexing ---

/// Convert 3D grid coordinates (x, y, z) to a flat array index.
///
/// Layout: x varies fastest (row-major in x), then y, then z.
/// idx = x + y * nx + z * nx * ny
#[inline]
pub fn idx(x: usize, y: usize, z: usize, nx: usize, ny: usize) -> usize {
    x + y * nx + z * nx * ny
}

/// Convert a flat array index back to 3D grid coordinates (x, y, z).
#[inline]
pub fn idx_to_xyz(i: usize, nx: usize, ny: usize) -> (usize, usize, usize) {
    let z = i / (nx * ny);
    let rem = i % (nx * ny);
    let y = rem / nx;
    let x = rem % nx;
    (x, y, z)
}

// --- Neighbor index functions ---

/// Neighbor indices for cell at (x, y, z), with boundary awareness.
/// Returns (x_minus, x_plus, y_minus, y_plus, z_minus, z_plus) as flat indices.
/// Clamps to grid boundaries (Neumann-style: repeat edge values).
#[inline]
pub fn neighbors_clamped(
    x: usize,
    y: usize,
    z: usize,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (usize, usize, usize, usize, usize, usize) {
    let xm = if x > 0 { x - 1 } else { 0 };
    let xp = if x + 1 < nx { x + 1 } else { nx - 1 };
    let ym = if y > 0 { y - 1 } else { 0 };
    let yp = if y + 1 < ny { y + 1 } else { ny - 1 };
    let zm = if z > 0 { z - 1 } else { 0 };
    let zp = if z + 1 < nz { z + 1 } else { nz - 1 };

    (
        idx(xm, y, z, nx, ny),
        idx(xp, y, z, nx, ny),
        idx(x, ym, z, nx, ny),
        idx(x, yp, z, nx, ny),
        idx(x, y, zm, nx, ny),
        idx(x, y, zp, nx, ny),
    )
}

/// Check if a cell is strictly interior (not on any face of the grid).
/// Interior cells are safe for central-difference stencils without boundary checks.
#[inline]
pub fn is_interior(x: usize, y: usize, z: usize, nx: usize, ny: usize, nz: usize) -> bool {
    x >= 1 && x + 1 < nx && y >= 1 && y + 1 < ny && z >= 1 && z + 1 < nz
}

// --- CFL condition ---

/// Compute the maximum stable timestep for explicit FDTD on a 3D grid.
///
/// CFL condition for 3D wave equation: dt <= dx / (c_max * sqrt(3))
///
/// `dx`: grid spacing (uniform).
/// `c_max`: maximum wave speed in the domain (typically c0 for vacuum).
///
/// Returns the CFL-limited timestep with a safety factor of 0.9.
#[inline]
pub fn max_dt(dx: f32, c_max: f32) -> f32 {
    // In 3D, the Courant number must satisfy: c * dt / dx <= 1/sqrt(3)
    // So dt <= dx / (c * sqrt(3))
    let cfl_limit = dx / (c_max * 3.0_f32.sqrt());
    cfl_limit * 0.9 // 90% safety factor
}

// --- Stencil coefficients ---

/// 2nd-order central difference stencil coefficients for first derivative.
/// df/dx ≈ (-1/2 * f[i-1] + 0 * f[i] + 1/2 * f[i+1]) / dx
pub const STENCIL_1ST_ORDER_2: [f32; 3] = [-0.5, 0.0, 0.5];

/// 2nd-order central difference stencil coefficients for second derivative (Laplacian component).
/// d²f/dx² ≈ (1 * f[i-1] - 2 * f[i] + 1 * f[i+1]) / dx²
pub const STENCIL_2ND_ORDER_2: [f32; 3] = [1.0, -2.0, 1.0];

/// 4th-order central difference stencil coefficients for first derivative.
/// df/dx ≈ (1/12 * f[i-2] - 2/3 * f[i-1] + 0 * f[i] + 2/3 * f[i+1] - 1/12 * f[i+2]) / dx
pub const STENCIL_1ST_ORDER_4: [f32; 5] = [1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0];

/// 4th-order central difference stencil coefficients for second derivative.
/// d²f/dx² ≈ (-1/12 * f[i-2] + 4/3 * f[i-1] - 5/2 * f[i] + 4/3 * f[i+1] - 1/12 * f[i+2]) / dx²
pub const STENCIL_2ND_ORDER_4: [f32; 5] = [
    -1.0 / 12.0,
    4.0 / 3.0,
    -5.0 / 2.0,
    4.0 / 3.0,
    -1.0 / 12.0,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_idx_roundtrip() {
        let (nx, ny, _nz) = (16, 16, 16);
        for z in 0..16 {
            for y in 0..16 {
                for x in 0..16 {
                    let i = idx(x, y, z, nx, ny);
                    let (rx, ry, rz) = idx_to_xyz(i, nx, ny);
                    assert_eq!((x, y, z), (rx, ry, rz), "roundtrip failed for ({x},{y},{z})");
                }
            }
        }
    }

    #[test]
    fn test_is_interior() {
        let (nx, ny, nz) = (8, 8, 8);
        // Corner — not interior
        assert!(!is_interior(0, 0, 0, nx, ny, nz));
        // Edge — not interior
        assert!(!is_interior(7, 4, 4, nx, ny, nz));
        assert!(!is_interior(4, 0, 4, nx, ny, nz));
        // Interior
        assert!(is_interior(1, 1, 1, nx, ny, nz));
        assert!(is_interior(4, 4, 4, nx, ny, nz));
        assert!(is_interior(6, 6, 6, nx, ny, nz));
    }

    #[test]
    fn test_neighbors_clamped_interior() {
        let (nx, ny, nz) = (8, 8, 8);
        let (xm, xp, ym, yp, zm, zp) = neighbors_clamped(4, 4, 4, nx, ny, nz);
        assert_eq!(xm, idx(3, 4, 4, nx, ny));
        assert_eq!(xp, idx(5, 4, 4, nx, ny));
        assert_eq!(ym, idx(4, 3, 4, nx, ny));
        assert_eq!(yp, idx(4, 5, 4, nx, ny));
        assert_eq!(zm, idx(4, 4, 3, nx, ny));
        assert_eq!(zp, idx(4, 4, 5, nx, ny));
    }

    #[test]
    fn test_neighbors_clamped_boundary() {
        let (nx, ny, nz) = (8, 8, 8);
        let (xm, xp, ym, yp, zm, zp) = neighbors_clamped(0, 0, 0, nx, ny, nz);
        // At origin, minus neighbors clamp to 0
        assert_eq!(xm, idx(0, 0, 0, nx, ny));
        assert_eq!(xp, idx(1, 0, 0, nx, ny));
        assert_eq!(ym, idx(0, 0, 0, nx, ny));
        assert_eq!(yp, idx(0, 1, 0, nx, ny));
        assert_eq!(zm, idx(0, 0, 0, nx, ny));
        assert_eq!(zp, idx(0, 0, 1, nx, ny));
    }

    #[test]
    fn test_cfl_condition() {
        let dx = 0.01; // 1 cm
        let c = 3.0e8;  // speed of light
        let dt = max_dt(dx, c);
        // dt should be < dx / (c * sqrt(3)) and positive
        let cfl_limit = dx / (c * 3.0_f32.sqrt());
        assert!(dt > 0.0);
        assert!(dt < cfl_limit, "dt should be below CFL limit");
        assert!((dt - cfl_limit * 0.9).abs() < 1e-25, "dt should be 0.9 * CFL limit");
    }

    #[test]
    fn test_stencil_coefficients_sum() {
        // First derivative stencils should sum to 0 (antisymmetric)
        let sum_1st_2: f32 = STENCIL_1ST_ORDER_2.iter().sum();
        assert!(sum_1st_2.abs() < 1e-6, "1st order 2nd-accuracy stencil should sum to 0");

        let sum_1st_4: f32 = STENCIL_1ST_ORDER_4.iter().sum();
        assert!(sum_1st_4.abs() < 1e-6, "1st order 4th-accuracy stencil should sum to 0");

        // Second derivative stencils should also sum to 0 (property of Laplacian stencils)
        let sum_2nd_2: f32 = STENCIL_2ND_ORDER_2.iter().sum();
        assert!(sum_2nd_2.abs() < 1e-6, "2nd order 2nd-accuracy stencil should sum to 0");

        let sum_2nd_4: f32 = STENCIL_2ND_ORDER_4.iter().sum();
        assert!(sum_2nd_4.abs() < 1e-6, "2nd order 4th-accuracy stencil should sum to 0");
    }
}
