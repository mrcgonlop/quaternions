// 3D vector field operations: gradient, divergence, curl, Laplacian
// Used by the FDTD simulation loop for spatial derivatives of quaternionic potentials.
// All operations use central finite differences on a uniform grid.

/// Type alias for 3D vectors stored as [f32; 3] — [x, y, z].
pub type Vec3f = [f32; 3];

// --- Basic vector operations ---

/// Dot product of two 3D vectors.
#[inline]
pub fn dot(a: Vec3f, b: Vec3f) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3D vectors.
#[inline]
pub fn cross(a: Vec3f, b: Vec3f) -> Vec3f {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean magnitude of a 3D vector.
#[inline]
pub fn magnitude(v: Vec3f) -> f32 {
    dot(v, v).sqrt()
}

/// Squared magnitude of a 3D vector (avoids sqrt).
#[inline]
pub fn magnitude_sq(v: Vec3f) -> f32 {
    dot(v, v)
}

/// Component-wise addition of two 3D vectors.
#[inline]
pub fn add(a: Vec3f, b: Vec3f) -> Vec3f {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Component-wise subtraction of two 3D vectors.
#[inline]
pub fn sub(a: Vec3f, b: Vec3f) -> Vec3f {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiplication of a 3D vector.
#[inline]
pub fn scale(v: Vec3f, s: f32) -> Vec3f {
    [v[0] * s, v[1] * s, v[2] * s]
}

// --- Finite-difference spatial derivative operations ---
//
// All functions below operate on flat arrays indexed by a grid helper function.
// The grid is assumed uniform with spacing `dx` in all directions.
// Central differences: df/dx ≈ (f[i+1] - f[i-1]) / (2*dx)
//
// Grid indexing convention:
//   idx(x, y, z) = x + y * nx + z * nx * ny
// where nx, ny are the grid dimensions along X and Y.
//
// Callers must ensure (x, y, z) are interior points (1..n-2) to avoid
// out-of-bounds access. Boundary handling is done at a higher level.

/// Gradient of a scalar field using 2nd-order central differences.
///
/// grad(f) = (df/dx, df/dy, df/dz)
///
/// `scalar_field`: flat array of scalar values indexed by grid position.
/// `idx`: linear index of the cell at (x, y, z).
/// `nx`, `ny`: grid dimensions used to compute neighbor offsets.
/// `inv_2dx`: precomputed 1.0 / (2.0 * dx).
#[inline]
pub fn gradient_scalar(
    scalar_field: &[f32],
    idx: usize,
    nx: usize,
    ny: usize,
    inv_2dx: f32,
) -> Vec3f {
    // Neighbor offsets in the flat array:
    //   +x: idx + 1,          -x: idx - 1
    //   +y: idx + nx,         -y: idx - nx
    //   +z: idx + nx*ny,      -z: idx - nx*ny
    let stride_y = nx;
    let stride_z = nx * ny;

    [
        (scalar_field[idx + 1] - scalar_field[idx - 1]) * inv_2dx,
        (scalar_field[idx + stride_y] - scalar_field[idx - stride_y]) * inv_2dx,
        (scalar_field[idx + stride_z] - scalar_field[idx - stride_z]) * inv_2dx,
    ]
}

/// Divergence of a vector field using 2nd-order central differences.
///
/// div(V) = dVx/dx + dVy/dy + dVz/dz
///
/// `vx`, `vy`, `vz`: flat arrays for each component of the vector field.
/// `idx`: linear index of the cell at (x, y, z).
/// `nx`, `ny`: grid dimensions.
/// `inv_2dx`: precomputed 1.0 / (2.0 * dx).
#[inline]
pub fn divergence_vector(
    vx: &[f32],
    vy: &[f32],
    vz: &[f32],
    idx: usize,
    nx: usize,
    ny: usize,
    inv_2dx: f32,
) -> f32 {
    let stride_y = nx;
    let stride_z = nx * ny;

    // dVx/dx + dVy/dy + dVz/dz using central differences
    (vx[idx + 1] - vx[idx - 1]) * inv_2dx
        + (vy[idx + stride_y] - vy[idx - stride_y]) * inv_2dx
        + (vz[idx + stride_z] - vz[idx - stride_z]) * inv_2dx
}

/// Curl of a vector field using 2nd-order central differences.
///
/// curl(V) = (dVz/dy - dVy/dz, dVx/dz - dVz/dx, dVy/dx - dVx/dy)
///
/// `vx`, `vy`, `vz`: flat arrays for each component of the vector field.
/// `idx`: linear index of the cell at (x, y, z).
/// `nx`, `ny`: grid dimensions.
/// `inv_2dx`: precomputed 1.0 / (2.0 * dx).
#[inline]
pub fn curl_vector(
    vx: &[f32],
    vy: &[f32],
    vz: &[f32],
    idx: usize,
    nx: usize,
    ny: usize,
    inv_2dx: f32,
) -> Vec3f {
    let stride_y = nx;
    let stride_z = nx * ny;

    [
        // (dVz/dy - dVy/dz)
        (vz[idx + stride_y] - vz[idx - stride_y]) * inv_2dx
            - (vy[idx + stride_z] - vy[idx - stride_z]) * inv_2dx,
        // (dVx/dz - dVz/dx)
        (vx[idx + stride_z] - vx[idx - stride_z]) * inv_2dx
            - (vz[idx + 1] - vz[idx - 1]) * inv_2dx,
        // (dVy/dx - dVx/dy)
        (vy[idx + 1] - vy[idx - 1]) * inv_2dx
            - (vx[idx + stride_y] - vx[idx - stride_y]) * inv_2dx,
    ]
}

/// 7-point stencil Laplacian of a scalar field.
///
/// ∇²f ≈ (f[x+1] + f[x-1] + f[y+1] + f[y-1] + f[z+1] + f[z-1] - 6*f[center]) / dx²
///
/// `scalar_field`: flat array of scalar values.
/// `idx`: linear index of the cell at (x, y, z).
/// `nx`, `ny`: grid dimensions.
/// `inv_dx2`: precomputed 1.0 / (dx * dx).
#[inline]
pub fn laplacian_scalar(
    scalar_field: &[f32],
    idx: usize,
    nx: usize,
    ny: usize,
    inv_dx2: f32,
) -> f32 {
    let stride_y = nx;
    let stride_z = nx * ny;

    (scalar_field[idx + 1]
        + scalar_field[idx - 1]
        + scalar_field[idx + stride_y]
        + scalar_field[idx - stride_y]
        + scalar_field[idx + stride_z]
        + scalar_field[idx - stride_z]
        - 6.0 * scalar_field[idx])
        * inv_dx2
}

/// Component-wise 7-point stencil Laplacian of a vector field.
///
/// ∇²V = (∇²Vx, ∇²Vy, ∇²Vz)
///
/// `vx`, `vy`, `vz`: flat arrays for each component.
/// `idx`: linear index of the cell at (x, y, z).
/// `nx`, `ny`: grid dimensions.
/// `inv_dx2`: precomputed 1.0 / (dx * dx).
#[inline]
pub fn laplacian_vector(
    vx: &[f32],
    vy: &[f32],
    vz: &[f32],
    idx: usize,
    nx: usize,
    ny: usize,
    inv_dx2: f32,
) -> Vec3f {
    [
        laplacian_scalar(vx, idx, nx, ny, inv_dx2),
        laplacian_scalar(vy, idx, nx, ny, inv_dx2),
        laplacian_scalar(vz, idx, nx, ny, inv_dx2),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::fdtd;

    const EPS: f32 = 1e-4;

    /// Create a test grid and fill a scalar field with a linear function f(x,y,z) = ax + by + cz.
    /// The gradient of this field should be (a, b, c) everywhere (interior).
    #[test]
    fn test_gradient_of_linear_field() {
        let (nx, ny, nz) = (8, 8, 8);
        let dx = 0.1;
        let inv_2dx = 1.0 / (2.0 * dx);
        let n = nx * ny * nz;

        let (a, b, c) = (2.0_f32, -3.0_f32, 1.5_f32);
        let mut field = vec![0.0f32; n];

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    field[idx] = a * x as f32 * dx + b * y as f32 * dx + c * z as f32 * dx;
                }
            }
        }

        // Check gradient at an interior point
        let idx = fdtd::idx(4, 4, 4, nx, ny);
        let grad = gradient_scalar(&field, idx, nx, ny, inv_2dx);

        assert!((grad[0] - a).abs() < EPS, "dfdx: expected {a}, got {}", grad[0]);
        assert!((grad[1] - b).abs() < EPS, "dfdy: expected {b}, got {}", grad[1]);
        assert!((grad[2] - c).abs() < EPS, "dfdz: expected {c}, got {}", grad[2]);
    }

    /// Divergence of a constant vector field should be zero.
    #[test]
    fn test_divergence_of_constant_field() {
        let (nx, ny, nz) = (8, 8, 8);
        let dx = 0.1;
        let inv_2dx = 1.0 / (2.0 * dx);
        let n = nx * ny * nz;

        // Constant vector field V = (3, -2, 7)
        let vx = vec![3.0f32; n];
        let vy = vec![-2.0f32; n];
        let vz = vec![7.0f32; n];

        let idx = fdtd::idx(4, 4, 4, nx, ny);
        let div = divergence_vector(&vx, &vy, &vz, idx, nx, ny, inv_2dx);

        assert!(div.abs() < EPS, "div of constant field should be 0, got {div}");
    }

    /// Curl of a gradient should be zero (to numerical precision).
    /// Take f(x,y,z) = x² + y² + z², so grad(f) = (2x, 2y, 2z).
    /// curl(grad(f)) should be zero.
    #[test]
    fn test_curl_of_gradient_is_zero() {
        let (nx, ny, nz) = (8, 8, 8);
        let dx = 0.1;
        let inv_2dx = 1.0 / (2.0 * dx);
        let n = nx * ny * nz;

        // grad(f) = (2x, 2y, 2z) for f = x² + y² + z²
        let mut vx = vec![0.0f32; n];
        let mut vy = vec![0.0f32; n];
        let mut vz = vec![0.0f32; n];

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    vx[idx] = 2.0 * x as f32 * dx;
                    vy[idx] = 2.0 * y as f32 * dx;
                    vz[idx] = 2.0 * z as f32 * dx;
                }
            }
        }

        let idx = fdtd::idx(4, 4, 4, nx, ny);
        let c = curl_vector(&vx, &vy, &vz, idx, nx, ny, inv_2dx);

        assert!(c[0].abs() < EPS, "curl_x should be 0, got {}", c[0]);
        assert!(c[1].abs() < EPS, "curl_y should be 0, got {}", c[1]);
        assert!(c[2].abs() < EPS, "curl_z should be 0, got {}", c[2]);
    }

    /// Laplacian of a quadratic f(x,y,z) = x² + y² + z² should be constant = 2 + 2 + 2 = 6.
    #[test]
    fn test_laplacian_of_quadratic() {
        let (nx, ny, nz) = (8, 8, 8);
        let dx = 0.1;
        let inv_dx2 = 1.0 / (dx * dx);
        let n = nx * ny * nz;

        let mut field = vec![0.0f32; n];

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    let fx = x as f32 * dx;
                    let fy = y as f32 * dx;
                    let fz = z as f32 * dx;
                    field[idx] = fx * fx + fy * fy + fz * fz;
                }
            }
        }

        let idx = fdtd::idx(4, 4, 4, nx, ny);
        let lap = laplacian_scalar(&field, idx, nx, ny, inv_dx2);

        // ∇²(x²+y²+z²) = 2 + 2 + 2 = 6
        assert!((lap - 6.0).abs() < EPS, "laplacian should be 6, got {lap}");
    }

    // --- Basic vector operation tests ---

    #[test]
    fn test_dot_product() {
        assert!((dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < EPS);
    }

    #[test]
    fn test_cross_product() {
        // i × j = k
        let result = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((result[0]).abs() < EPS);
        assert!((result[1]).abs() < EPS);
        assert!((result[2] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_magnitude() {
        assert!((magnitude([3.0, 4.0, 0.0]) - 5.0).abs() < EPS);
    }
}
