// Quat type with Hamilton algebra
// See ARCHITECTURE.md §Quaternion Type (lines 140–174)

use bytemuck::{Pod, Zeroable};
use std::fmt;
use std::ops::{Add, Sub, Mul, Div, Neg};

/// Hamilton quaternion: q = w + xi + yj + zk
///
/// In the EM context:
/// - w = scalar part (φ/c for the electromagnetic four-potential)
/// - (x, y, z) = vector part (Ax, Ay, Az)
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Quat {
    pub w: f32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

// SAFETY: Quat is #[repr(C)] with four f32 fields, no padding, no invalid bit patterns.
unsafe impl Zeroable for Quat {}
unsafe impl Pod for Quat {}

impl Quat {
    /// Create a quaternion from four components.
    #[inline]
    pub fn new(w: f32, x: f32, y: f32, z: f32) -> Self {
        Self { w, x, y, z }
    }

    /// Zero quaternion (additive identity).
    #[inline]
    pub fn zero() -> Self {
        Self { w: 0.0, x: 0.0, y: 0.0, z: 0.0 }
    }

    /// Multiplicative identity: 1 + 0i + 0j + 0k.
    #[inline]
    pub fn identity() -> Self {
        Self { w: 1.0, x: 0.0, y: 0.0, z: 0.0 }
    }

    /// Pure scalar quaternion: s + 0i + 0j + 0k.
    #[inline]
    pub fn from_scalar(s: f32) -> Self {
        Self { w: s, x: 0.0, y: 0.0, z: 0.0 }
    }

    /// Pure vector quaternion: 0 + v[0]i + v[1]j + v[2]k.
    #[inline]
    pub fn from_vector(v: [f32; 3]) -> Self {
        Self { w: 0.0, x: v[0], y: v[1], z: v[2] }
    }

    /// Basis element i = (0, 1, 0, 0).
    #[inline]
    pub fn i() -> Self {
        Self { w: 0.0, x: 1.0, y: 0.0, z: 0.0 }
    }

    /// Basis element j = (0, 0, 1, 0).
    #[inline]
    pub fn j() -> Self {
        Self { w: 0.0, x: 0.0, y: 1.0, z: 0.0 }
    }

    /// Basis element k = (0, 0, 0, 1).
    #[inline]
    pub fn k() -> Self {
        Self { w: 0.0, x: 0.0, y: 0.0, z: 1.0 }
    }

    /// Hamilton product (non-commutative quaternion multiplication).
    ///
    /// (a + bi + cj + dk)(e + fi + gj + hk) =
    ///   (ae - bf - cg - dh)
    /// + (af + be + ch - dg)i
    /// + (ag - bh + ce + df)j
    /// + (ah + bg - cf + de)k
    #[inline]
    pub fn hamilton(self, rhs: Quat) -> Quat {
        Quat {
            w: self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            x: self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
        }
    }

    /// Conjugate: q* = w - xi - yj - zk.
    #[inline]
    pub fn conj(self) -> Quat {
        Quat { w: self.w, x: -self.x, y: -self.y, z: -self.z }
    }

    /// Norm squared: |q|² = w² + x² + y² + z².
    #[inline]
    pub fn norm_sq(self) -> f32 {
        self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Norm (magnitude): |q| = sqrt(w² + x² + y² + z²).
    #[inline]
    pub fn norm(self) -> f32 {
        self.norm_sq().sqrt()
    }

    /// Multiplicative inverse: q⁻¹ = q* / |q|².
    /// Returns None for zero quaternions.
    #[inline]
    pub fn inverse(self) -> Option<Quat> {
        let n2 = self.norm_sq();
        if n2 < f32::EPSILON * f32::EPSILON {
            None
        } else {
            let inv_n2 = 1.0 / n2;
            Some(Quat {
                w: self.w * inv_n2,
                x: -self.x * inv_n2,
                y: -self.y * inv_n2,
                z: -self.z * inv_n2,
            })
        }
    }

    /// Normalize to unit quaternion. Returns None for zero quaternions.
    #[inline]
    pub fn normalize(self) -> Option<Quat> {
        let n = self.norm();
        if n < f32::EPSILON {
            None
        } else {
            let inv_n = 1.0 / n;
            Some(Quat {
                w: self.w * inv_n,
                x: self.x * inv_n,
                y: self.y * inv_n,
                z: self.z * inv_n,
            })
        }
    }

    /// Extract scalar part: w.
    #[inline]
    pub fn scalar(self) -> f32 {
        self.w
    }

    /// Extract vector part: [x, y, z].
    #[inline]
    pub fn vector(self) -> [f32; 3] {
        [self.x, self.y, self.z]
    }
}

// --- Operator overloads ---

impl Add for Quat {
    type Output = Quat;
    #[inline]
    fn add(self, rhs: Quat) -> Quat {
        Quat {
            w: self.w + rhs.w,
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Quat {
    type Output = Quat;
    #[inline]
    fn sub(self, rhs: Quat) -> Quat {
        Quat {
            w: self.w - rhs.w,
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f32> for Quat {
    type Output = Quat;
    #[inline]
    fn mul(self, s: f32) -> Quat {
        Quat {
            w: self.w * s,
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

impl Mul<Quat> for f32 {
    type Output = Quat;
    #[inline]
    fn mul(self, q: Quat) -> Quat {
        q * self
    }
}

impl Div<f32> for Quat {
    type Output = Quat;
    #[inline]
    fn div(self, s: f32) -> Quat {
        let inv = 1.0 / s;
        Quat {
            w: self.w * inv,
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

impl Neg for Quat {
    type Output = Quat;
    #[inline]
    fn neg(self) -> Quat {
        Quat {
            w: -self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl fmt::Display for Quat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}i + {}j + {}k)", self.w, self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f32 = 1e-6;

    fn approx_eq(a: Quat, b: Quat) -> bool {
        (a.w - b.w).abs() < EPS
            && (a.x - b.x).abs() < EPS
            && (a.y - b.y).abs() < EPS
            && (a.z - b.z).abs() < EPS
    }

    // --- Hamilton product basis rules ---

    #[test]
    fn test_i_times_j_eq_k() {
        let result = Quat::i().hamilton(Quat::j());
        assert!(approx_eq(result, Quat::k()), "i*j should be k, got {result}");
    }

    #[test]
    fn test_j_times_i_eq_neg_k() {
        // Non-commutativity: j * i = -k
        let result = Quat::j().hamilton(Quat::i());
        assert!(approx_eq(result, -Quat::k()), "j*i should be -k, got {result}");
    }

    #[test]
    fn test_j_times_k_eq_i() {
        let result = Quat::j().hamilton(Quat::k());
        assert!(approx_eq(result, Quat::i()), "j*k should be i, got {result}");
    }

    #[test]
    fn test_k_times_j_eq_neg_i() {
        let result = Quat::k().hamilton(Quat::j());
        assert!(approx_eq(result, -Quat::i()), "k*j should be -i, got {result}");
    }

    #[test]
    fn test_k_times_i_eq_j() {
        let result = Quat::k().hamilton(Quat::i());
        assert!(approx_eq(result, Quat::j()), "k*i should be j, got {result}");
    }

    #[test]
    fn test_i_times_k_eq_neg_j() {
        let result = Quat::i().hamilton(Quat::k());
        assert!(approx_eq(result, -Quat::j()), "i*k should be -j, got {result}");
    }

    // --- i² = j² = k² = ijk = -1 ---

    #[test]
    fn test_i_squared_eq_neg_one() {
        let result = Quat::i().hamilton(Quat::i());
        assert!(approx_eq(result, Quat::from_scalar(-1.0)), "i² should be -1, got {result}");
    }

    #[test]
    fn test_j_squared_eq_neg_one() {
        let result = Quat::j().hamilton(Quat::j());
        assert!(approx_eq(result, Quat::from_scalar(-1.0)), "j² should be -1, got {result}");
    }

    #[test]
    fn test_k_squared_eq_neg_one() {
        let result = Quat::k().hamilton(Quat::k());
        assert!(approx_eq(result, Quat::from_scalar(-1.0)), "k² should be -1, got {result}");
    }

    #[test]
    fn test_ijk_eq_neg_one() {
        let ij = Quat::i().hamilton(Quat::j());
        let ijk = ij.hamilton(Quat::k());
        assert!(approx_eq(ijk, Quat::from_scalar(-1.0)), "ijk should be -1, got {ijk}");
    }

    // --- Norm preservation under multiplication ---

    #[test]
    fn test_norm_preservation() {
        let a = Quat::new(1.0, 2.0, 3.0, 4.0);
        let b = Quat::new(-0.5, 1.5, -2.5, 0.7);
        let product = a.hamilton(b);
        let expected_norm = a.norm() * b.norm();
        let actual_norm = product.norm();
        assert!(
            (expected_norm - actual_norm).abs() < EPS * expected_norm,
            "|a*b| should be |a|*|b|: expected {expected_norm}, got {actual_norm}"
        );
    }

    // --- Inverse ---

    #[test]
    fn test_inverse_gives_identity() {
        let q = Quat::new(1.0, 2.0, -3.0, 0.5);
        let q_inv = q.inverse().expect("non-zero quaternion should have inverse");
        let product = q.hamilton(q_inv);
        assert!(
            approx_eq(product, Quat::identity()),
            "q * q⁻¹ should be identity, got {product}"
        );
    }

    #[test]
    fn test_right_inverse() {
        let q = Quat::new(-2.0, 0.3, 1.7, -0.9);
        let q_inv = q.inverse().expect("non-zero quaternion should have inverse");
        let product = q_inv.hamilton(q);
        assert!(
            approx_eq(product, Quat::identity()),
            "q⁻¹ * q should be identity, got {product}"
        );
    }

    #[test]
    fn test_zero_has_no_inverse() {
        assert!(Quat::zero().inverse().is_none());
    }

    // --- Scalar/vector decomposition roundtrip ---

    #[test]
    fn test_scalar_vector_roundtrip() {
        let q = Quat::new(3.0, -1.0, 2.0, 0.5);
        let reconstructed = Quat::from_scalar(q.scalar()) + Quat::from_vector(q.vector());
        assert!(approx_eq(q, reconstructed), "scalar + vector should reconstruct q");
    }

    // --- Conjugate properties ---

    #[test]
    fn test_conjugate_conjugate_is_identity() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        assert!(approx_eq(q.conj().conj(), q));
    }

    #[test]
    fn test_q_times_conj_is_norm_sq() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        let result = q.hamilton(q.conj());
        let expected = Quat::from_scalar(q.norm_sq());
        assert!(approx_eq(result, expected), "q*q* should be |q|², got {result}");
    }

    // --- Normalize ---

    #[test]
    fn test_normalize() {
        let q = Quat::new(3.0, 4.0, 0.0, 0.0);
        let n = q.normalize().expect("non-zero");
        assert!((n.norm() - 1.0).abs() < EPS, "normalized should have unit norm");
    }

    #[test]
    fn test_normalize_zero_returns_none() {
        assert!(Quat::zero().normalize().is_none());
    }

    // --- Operator overloads ---

    #[test]
    fn test_add() {
        let a = Quat::new(1.0, 2.0, 3.0, 4.0);
        let b = Quat::new(0.5, -1.0, 1.5, -2.0);
        assert!(approx_eq(a + b, Quat::new(1.5, 1.0, 4.5, 2.0)));
    }

    #[test]
    fn test_sub() {
        let a = Quat::new(1.0, 2.0, 3.0, 4.0);
        let b = Quat::new(0.5, -1.0, 1.5, -2.0);
        assert!(approx_eq(a - b, Quat::new(0.5, 3.0, 1.5, 6.0)));
    }

    #[test]
    fn test_scalar_mul() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        assert!(approx_eq(q * 2.0, Quat::new(2.0, 4.0, 6.0, 8.0)));
    }

    #[test]
    fn test_scalar_mul_commutative() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        assert!(approx_eq(q * 3.0, 3.0 * q));
    }

    #[test]
    fn test_scalar_div() {
        let q = Quat::new(2.0, 4.0, 6.0, 8.0);
        assert!(approx_eq(q / 2.0, Quat::new(1.0, 2.0, 3.0, 4.0)));
    }

    #[test]
    fn test_neg() {
        let q = Quat::new(1.0, -2.0, 3.0, -4.0);
        assert!(approx_eq(-q, Quat::new(-1.0, 2.0, -3.0, 4.0)));
    }

    // --- Display ---

    #[test]
    fn test_display() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        let s = format!("{q}");
        assert!(s.contains("1") && s.contains("2") && s.contains("3") && s.contains("4"));
    }

    // --- Bytemuck ---

    #[test]
    fn test_bytemuck_cast() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0);
        let bytes: &[u8] = bytemuck::bytes_of(&q);
        assert_eq!(bytes.len(), 16); // 4 x f32
        let roundtrip: &Quat = bytemuck::from_bytes(bytes);
        assert_eq!(*roundtrip, q);
    }

    #[test]
    fn test_bytemuck_slice_cast() {
        let quats = vec![Quat::new(1.0, 0.0, 0.0, 0.0), Quat::new(0.0, 1.0, 0.0, 0.0)];
        let bytes: &[u8] = bytemuck::cast_slice(&quats);
        assert_eq!(bytes.len(), 32); // 2 * 16
        let roundtrip: &[Quat] = bytemuck::cast_slice(bytes);
        assert_eq!(roundtrip, &quats[..]);
    }
}
