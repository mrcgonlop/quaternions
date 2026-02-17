// Scientific color maps for field visualization.
// Each map is a 32-entry RGBA lookup table with linear interpolation.

use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

/// Available color maps for field visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorMap {
    /// Perceptually uniform sequential (dark purple → yellow). Good for magnitudes.
    Viridis,
    /// Perceptually uniform sequential (dark purple → bright yellow). High contrast.
    Plasma,
    /// Diverging (blue → white → red). Good for signed quantities near zero.
    Coolwarm,
    /// Diverging bipolar (blue → white/transparent → red). Critical for S field.
    /// Zero maps to transparent white, making it ideal for overlay visualization.
    ScalarField,
    /// Simple grayscale (black → white).
    Grayscale,
}

impl ColorMap {
    /// All available color maps, for UI dropdowns.
    pub const ALL: &'static [ColorMap] = &[
        ColorMap::Viridis,
        ColorMap::Plasma,
        ColorMap::Coolwarm,
        ColorMap::ScalarField,
        ColorMap::Grayscale,
    ];

    /// Human-readable name for UI display.
    pub fn name(&self) -> &'static str {
        match self {
            ColorMap::Viridis => "Viridis",
            ColorMap::Plasma => "Plasma",
            ColorMap::Coolwarm => "Coolwarm",
            ColorMap::ScalarField => "Scalar Field",
            ColorMap::Grayscale => "Grayscale",
        }
    }
}

/// Map a scalar value to an RGBA color using the specified color map.
///
/// `value` is clamped and normalized to [0, 1] using `min` and `max`.
/// Returns [r, g, b, a] with each component in [0.0, 1.0].
pub fn map_value(value: f32, min: f32, max: f32, map: ColorMap) -> [f32; 4] {
    let range = max - min;
    let t = if range.abs() < f32::EPSILON {
        0.5
    } else {
        ((value - min) / range).clamp(0.0, 1.0)
    };

    let lut = lut_for(map);
    sample_lut(lut, t)
}

/// Map a scalar value through the ScalarField bipolar map.
///
/// Convenience function for signed fields like S.
/// Maps: negative → blue, zero → transparent white, positive → red.
/// `abs_max` is the symmetric range: value is mapped over [-abs_max, +abs_max].
pub fn map_bipolar(value: f32, abs_max: f32) -> [f32; 4] {
    map_value(value, -abs_max, abs_max, ColorMap::ScalarField)
}

/// Sample a LUT with linear interpolation at position t in [0, 1].
fn sample_lut(lut: &[[f32; 4]], t: f32) -> [f32; 4] {
    let n = lut.len();
    let pos = t * (n - 1) as f32;
    let idx = (pos as usize).min(n - 2);
    let frac = pos - idx as f32;

    let c0 = lut[idx];
    let c1 = lut[idx + 1];

    [
        c0[0] + (c1[0] - c0[0]) * frac,
        c0[1] + (c1[1] - c0[1]) * frac,
        c0[2] + (c1[2] - c0[2]) * frac,
        c0[3] + (c1[3] - c0[3]) * frac,
    ]
}

/// Get the LUT for a given color map.
fn lut_for(map: ColorMap) -> &'static [[f32; 4]] {
    match map {
        ColorMap::Viridis => &LUT_VIRIDIS,
        ColorMap::Plasma => &LUT_PLASMA,
        ColorMap::Coolwarm => &LUT_COOLWARM,
        ColorMap::ScalarField => &LUT_SCALAR_FIELD,
        ColorMap::Grayscale => &LUT_GRAYSCALE,
    }
}

/// Generate a 1D Bevy Image (256 pixels) for use as a GPU texture.
///
/// The texture is RGBA8 (sRGB), 256x1 pixels, suitable for sampling
/// in shaders with `textureSample(colormap_tex, sampler, vec2(t, 0.5))`.
pub fn generate_texture(map: ColorMap) -> Image {
    let width = 256u32;
    let mut data = Vec::with_capacity((width * 4) as usize);

    for i in 0..width {
        let t = i as f32 / (width - 1) as f32;
        let lut = lut_for(map);
        let c = sample_lut(lut, t);
        data.push((c[0] * 255.0) as u8);
        data.push((c[1] * 255.0) as u8);
        data.push((c[2] * 255.0) as u8);
        data.push((c[3] * 255.0) as u8);
    }

    Image::new(
        Extent3d {
            width,
            height: 1,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        data,
        TextureFormat::Rgba8UnormSrgb,
        default(),
    )
}

// =============================================================================
// Lookup Tables — 32 entries each, [R, G, B, A] in [0.0, 1.0]
// =============================================================================

/// Viridis: perceptually uniform sequential colormap (matplotlib).
/// 32 samples from the canonical 256-entry table.
#[rustfmt::skip]
const LUT_VIRIDIS: [[f32; 4]; 32] = [
    [0.267, 0.004, 0.329, 1.0],
    [0.277, 0.050, 0.375, 1.0],
    [0.282, 0.094, 0.417, 1.0],
    [0.278, 0.137, 0.451, 1.0],
    [0.270, 0.176, 0.477, 1.0],
    [0.258, 0.213, 0.497, 1.0],
    [0.244, 0.249, 0.510, 1.0],
    [0.228, 0.283, 0.518, 1.0],
    [0.213, 0.316, 0.521, 1.0],
    [0.197, 0.347, 0.522, 1.0],
    [0.183, 0.377, 0.519, 1.0],
    [0.170, 0.407, 0.514, 1.0],
    [0.157, 0.435, 0.506, 1.0],
    [0.146, 0.462, 0.497, 1.0],
    [0.137, 0.489, 0.485, 1.0],
    [0.130, 0.515, 0.470, 1.0],
    [0.131, 0.540, 0.452, 1.0],
    [0.140, 0.565, 0.432, 1.0],
    [0.159, 0.589, 0.408, 1.0],
    [0.190, 0.612, 0.380, 1.0],
    [0.231, 0.634, 0.349, 1.0],
    [0.282, 0.654, 0.314, 1.0],
    [0.340, 0.673, 0.277, 1.0],
    [0.404, 0.690, 0.237, 1.0],
    [0.472, 0.706, 0.193, 1.0],
    [0.545, 0.720, 0.147, 1.0],
    [0.621, 0.731, 0.104, 1.0],
    [0.699, 0.739, 0.073, 1.0],
    [0.777, 0.744, 0.064, 1.0],
    [0.853, 0.745, 0.093, 1.0],
    [0.924, 0.741, 0.154, 1.0],
    [0.993, 0.906, 0.144, 1.0],
];

/// Plasma: perceptually uniform sequential colormap (matplotlib).
#[rustfmt::skip]
const LUT_PLASMA: [[f32; 4]; 32] = [
    [0.050, 0.030, 0.528, 1.0],
    [0.107, 0.024, 0.578, 1.0],
    [0.160, 0.016, 0.620, 1.0],
    [0.211, 0.009, 0.653, 1.0],
    [0.261, 0.005, 0.678, 1.0],
    [0.310, 0.008, 0.693, 1.0],
    [0.356, 0.019, 0.700, 1.0],
    [0.401, 0.039, 0.698, 1.0],
    [0.443, 0.062, 0.689, 1.0],
    [0.483, 0.087, 0.674, 1.0],
    [0.521, 0.112, 0.654, 1.0],
    [0.556, 0.137, 0.631, 1.0],
    [0.589, 0.162, 0.606, 1.0],
    [0.620, 0.188, 0.579, 1.0],
    [0.649, 0.213, 0.551, 1.0],
    [0.676, 0.239, 0.522, 1.0],
    [0.702, 0.265, 0.493, 1.0],
    [0.726, 0.291, 0.463, 1.0],
    [0.749, 0.319, 0.432, 1.0],
    [0.770, 0.347, 0.401, 1.0],
    [0.791, 0.376, 0.369, 1.0],
    [0.810, 0.407, 0.337, 1.0],
    [0.828, 0.439, 0.304, 1.0],
    [0.846, 0.474, 0.269, 1.0],
    [0.863, 0.511, 0.233, 1.0],
    [0.879, 0.551, 0.195, 1.0],
    [0.895, 0.595, 0.153, 1.0],
    [0.910, 0.642, 0.110, 1.0],
    [0.925, 0.693, 0.069, 1.0],
    [0.939, 0.748, 0.035, 1.0],
    [0.951, 0.808, 0.029, 1.0],
    [0.940, 0.975, 0.131, 1.0],
];

/// Coolwarm: diverging colormap (blue → white → red).
#[rustfmt::skip]
const LUT_COOLWARM: [[f32; 4]; 32] = [
    [0.230, 0.299, 0.754, 1.0],
    [0.266, 0.330, 0.777, 1.0],
    [0.303, 0.362, 0.799, 1.0],
    [0.342, 0.395, 0.819, 1.0],
    [0.383, 0.429, 0.838, 1.0],
    [0.424, 0.464, 0.855, 1.0],
    [0.467, 0.499, 0.870, 1.0],
    [0.510, 0.535, 0.884, 1.0],
    [0.554, 0.571, 0.896, 1.0],
    [0.598, 0.607, 0.906, 1.0],
    [0.642, 0.643, 0.914, 1.0],
    [0.686, 0.679, 0.920, 1.0],
    [0.729, 0.714, 0.924, 1.0],
    [0.771, 0.748, 0.925, 1.0],
    [0.812, 0.781, 0.924, 1.0],
    [0.850, 0.812, 0.920, 1.0],
    [0.884, 0.812, 0.883, 1.0],
    [0.899, 0.784, 0.838, 1.0],
    [0.910, 0.752, 0.790, 1.0],
    [0.917, 0.718, 0.741, 1.0],
    [0.921, 0.681, 0.690, 1.0],
    [0.920, 0.642, 0.639, 1.0],
    [0.916, 0.601, 0.588, 1.0],
    [0.908, 0.558, 0.537, 1.0],
    [0.897, 0.513, 0.487, 1.0],
    [0.882, 0.467, 0.438, 1.0],
    [0.863, 0.419, 0.390, 1.0],
    [0.841, 0.370, 0.343, 1.0],
    [0.816, 0.319, 0.298, 1.0],
    [0.788, 0.267, 0.255, 1.0],
    [0.757, 0.213, 0.214, 1.0],
    [0.706, 0.016, 0.150, 1.0],
];

/// ScalarField: bipolar diverging map for signed fields (S field).
/// negative (blue) → zero (transparent white) → positive (red).
/// Alpha ramps from 1.0 at extremes to 0.0 at zero — makes zero values transparent.
#[rustfmt::skip]
const LUT_SCALAR_FIELD: [[f32; 4]; 32] = [
    [0.019, 0.188, 0.760, 1.00],
    [0.071, 0.232, 0.793, 0.97],
    [0.129, 0.280, 0.824, 0.94],
    [0.190, 0.331, 0.851, 0.90],
    [0.254, 0.385, 0.876, 0.87],
    [0.320, 0.441, 0.897, 0.83],
    [0.388, 0.499, 0.916, 0.77],
    [0.456, 0.558, 0.931, 0.70],
    [0.525, 0.617, 0.944, 0.63],
    [0.593, 0.676, 0.954, 0.53],
    [0.660, 0.735, 0.963, 0.43],
    [0.726, 0.793, 0.971, 0.33],
    [0.789, 0.849, 0.978, 0.23],
    [0.849, 0.903, 0.985, 0.13],
    [0.907, 0.952, 0.992, 0.06],
    [0.960, 0.980, 0.998, 0.00],
    [0.998, 0.960, 0.960, 0.00],
    [0.992, 0.930, 0.907, 0.06],
    [0.985, 0.890, 0.849, 0.13],
    [0.978, 0.840, 0.789, 0.23],
    [0.971, 0.782, 0.726, 0.33],
    [0.963, 0.720, 0.660, 0.43],
    [0.954, 0.655, 0.593, 0.53],
    [0.944, 0.590, 0.525, 0.63],
    [0.931, 0.524, 0.456, 0.70],
    [0.916, 0.459, 0.388, 0.77],
    [0.897, 0.396, 0.320, 0.83],
    [0.876, 0.336, 0.254, 0.87],
    [0.851, 0.280, 0.190, 0.90],
    [0.824, 0.228, 0.129, 0.94],
    [0.793, 0.180, 0.071, 0.97],
    [0.760, 0.137, 0.019, 1.00],
];

// =============================================================================
// Multi-encoding color modes
// =============================================================================

/// Color encoding mode for visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorEncoding {
    /// Single-field scalar color map (existing behavior).
    Standard,
    /// Map 3 independent field quantities to R, G, B channels.
    RgbMultiField,
    /// Hue = field direction angle, Value = magnitude.
    HsvPhase,
    /// Glyph-only: size from one quantity, color from another.
    SizeColor,
}

impl ColorEncoding {
    pub const ALL: &'static [ColorEncoding] = &[
        ColorEncoding::Standard,
        ColorEncoding::RgbMultiField,
        ColorEncoding::HsvPhase,
        ColorEncoding::SizeColor,
    ];

    /// Encodings valid for 2D slice visualization (excludes SizeColor).
    pub const SLICE_MODES: &'static [ColorEncoding] = &[
        ColorEncoding::Standard,
        ColorEncoding::RgbMultiField,
        ColorEncoding::HsvPhase,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            ColorEncoding::Standard => "Standard",
            ColorEncoding::RgbMultiField => "RGB Multi-Field",
            ColorEncoding::HsvPhase => "HSV Phase",
            ColorEncoding::SizeColor => "Size + Color",
        }
    }
}

/// Which 2D plane to compute the phase angle from for HSV encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhasePlane {
    /// atan2(field_y, field_x)
    XY,
    /// atan2(field_z, field_x)
    XZ,
    /// atan2(field_z, field_y)
    YZ,
}

impl PhasePlane {
    pub const ALL: &'static [PhasePlane] = &[PhasePlane::XY, PhasePlane::XZ, PhasePlane::YZ];

    pub fn name(&self) -> &'static str {
        match self {
            PhasePlane::XY => "XY",
            PhasePlane::XZ => "XZ",
            PhasePlane::YZ => "YZ",
        }
    }
}

/// Map three independently-ranged scalar values to RGB channels.
///
/// Each channel is normalized to [0, 1] within its own (min, max) range.
/// Returns [R, G, B, 1.0].
pub fn encode_rgb_multi(
    r_val: f32,
    r_min: f32,
    r_max: f32,
    g_val: f32,
    g_min: f32,
    g_max: f32,
    b_val: f32,
    b_min: f32,
    b_max: f32,
) -> [f32; 4] {
    let normalize = |v: f32, lo: f32, hi: f32| -> f32 {
        let range = hi - lo;
        if range.abs() < 1e-6 {
            0.5
        } else {
            ((v - lo) / range).clamp(0.0, 1.0)
        }
    };
    [
        normalize(r_val, r_min, r_max),
        normalize(g_val, g_min, g_max),
        normalize(b_val, b_min, b_max),
        1.0,
    ]
}

/// HSV phase encoding: hue from field vector direction, value from magnitude.
///
/// `field_vec` is a 3-component field vector (e.g., E, B, or Poynting).
/// `plane` selects which two components determine the hue angle.
/// `max_mag` sets the maximum magnitude for brightness normalization.
///
/// Returns [R, G, B, A] with A=1.0.
pub fn encode_hsv_phase(field_vec: [f32; 3], plane: PhasePlane, max_mag: f32) -> [f32; 4] {
    let (a, b) = match plane {
        PhasePlane::XY => (field_vec[0], field_vec[1]),
        PhasePlane::XZ => (field_vec[0], field_vec[2]),
        PhasePlane::YZ => (field_vec[1], field_vec[2]),
    };

    // Hue from direction angle, normalized to [0, 1]
    let hue = b.atan2(a); // -pi..pi
    let hue_normalized = (hue + std::f32::consts::PI) / (2.0 * std::f32::consts::PI);

    // Value from magnitude
    let mag = (field_vec[0] * field_vec[0]
        + field_vec[1] * field_vec[1]
        + field_vec[2] * field_vec[2])
        .sqrt();
    let value = (mag / max_mag.max(1e-6)).clamp(0.0, 1.0);

    hsv_to_rgb(hue_normalized, 0.9, value)
}

/// Convert HSV (all components in [0, 1]) to RGBA.
pub fn hsv_to_rgb(h: f32, s: f32, v: f32) -> [f32; 4] {
    let h6 = h * 6.0;
    let c = v * s;
    let x = c * (1.0 - (h6 % 2.0 - 1.0).abs());
    let m = v - c;
    let (r, g, b) = match h6 as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    [r + m, g + m, b + m, 1.0]
}

// =============================================================================
// Lookup Tables — 32 entries each, [R, G, B, A] in [0.0, 1.0]
// =============================================================================

/// Grayscale: simple linear black → white.
#[rustfmt::skip]
const LUT_GRAYSCALE: [[f32; 4]; 32] = [
    [0.000, 0.000, 0.000, 1.0],
    [0.032, 0.032, 0.032, 1.0],
    [0.065, 0.065, 0.065, 1.0],
    [0.097, 0.097, 0.097, 1.0],
    [0.129, 0.129, 0.129, 1.0],
    [0.161, 0.161, 0.161, 1.0],
    [0.194, 0.194, 0.194, 1.0],
    [0.226, 0.226, 0.226, 1.0],
    [0.258, 0.258, 0.258, 1.0],
    [0.290, 0.290, 0.290, 1.0],
    [0.323, 0.323, 0.323, 1.0],
    [0.355, 0.355, 0.355, 1.0],
    [0.387, 0.387, 0.387, 1.0],
    [0.419, 0.419, 0.419, 1.0],
    [0.452, 0.452, 0.452, 1.0],
    [0.484, 0.484, 0.484, 1.0],
    [0.516, 0.516, 0.516, 1.0],
    [0.548, 0.548, 0.548, 1.0],
    [0.581, 0.581, 0.581, 1.0],
    [0.613, 0.613, 0.613, 1.0],
    [0.645, 0.645, 0.645, 1.0],
    [0.677, 0.677, 0.677, 1.0],
    [0.710, 0.710, 0.710, 1.0],
    [0.742, 0.742, 0.742, 1.0],
    [0.774, 0.774, 0.774, 1.0],
    [0.806, 0.806, 0.806, 1.0],
    [0.839, 0.839, 0.839, 1.0],
    [0.871, 0.871, 0.871, 1.0],
    [0.903, 0.903, 0.903, 1.0],
    [0.935, 0.935, 0.935, 1.0],
    [0.968, 0.968, 0.968, 1.0],
    [1.000, 1.000, 1.000, 1.0],
];

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f32 = 0.02; // LUT interpolation tolerance

    #[test]
    fn test_map_value_min_max() {
        // At min → should be first LUT entry
        let c = map_value(0.0, 0.0, 1.0, ColorMap::Grayscale);
        assert!(c[0] < 0.01 && c[1] < 0.01 && c[2] < 0.01, "min should be black");
        assert!((c[3] - 1.0).abs() < EPS, "alpha should be 1.0");

        // At max → should be last LUT entry
        let c = map_value(1.0, 0.0, 1.0, ColorMap::Grayscale);
        assert!(c[0] > 0.99 && c[1] > 0.99 && c[2] > 0.99, "max should be white");
    }

    #[test]
    fn test_map_value_midpoint_grayscale() {
        let c = map_value(0.5, 0.0, 1.0, ColorMap::Grayscale);
        assert!((c[0] - 0.5).abs() < EPS, "mid gray R: got {}", c[0]);
        assert!((c[1] - 0.5).abs() < EPS, "mid gray G: got {}", c[1]);
        assert!((c[2] - 0.5).abs() < EPS, "mid gray B: got {}", c[2]);
    }

    #[test]
    fn test_map_value_clamping() {
        // Values outside range should clamp
        let below = map_value(-10.0, 0.0, 1.0, ColorMap::Grayscale);
        let at_min = map_value(0.0, 0.0, 1.0, ColorMap::Grayscale);
        assert_eq!(below, at_min, "below-min should clamp to min");

        let above = map_value(100.0, 0.0, 1.0, ColorMap::Grayscale);
        let at_max = map_value(1.0, 0.0, 1.0, ColorMap::Grayscale);
        assert_eq!(above, at_max, "above-max should clamp to max");
    }

    #[test]
    fn test_scalar_field_bipolar() {
        // Zero should be nearly transparent
        let zero = map_bipolar(0.0, 1.0);
        assert!(zero[3] < 0.05, "zero should be transparent, got alpha={}", zero[3]);

        // Positive extreme should be red-ish with high alpha
        let pos = map_bipolar(1.0, 1.0);
        assert!(pos[0] > 0.5, "positive should be red, got R={}", pos[0]);
        assert!(pos[3] > 0.9, "positive should be opaque, got alpha={}", pos[3]);

        // Negative extreme should be blue-ish with high alpha
        let neg = map_bipolar(-1.0, 1.0);
        assert!(neg[2] > 0.5, "negative should be blue, got B={}", neg[2]);
        assert!(neg[3] > 0.9, "negative should be opaque, got alpha={}", neg[3]);
    }

    #[test]
    fn test_viridis_endpoints() {
        let start = map_value(0.0, 0.0, 1.0, ColorMap::Viridis);
        // Viridis starts dark purple
        assert!(start[0] < 0.3 && start[2] > 0.3, "viridis start should be purplish");

        let end = map_value(1.0, 0.0, 1.0, ColorMap::Viridis);
        // Viridis ends bright yellow
        assert!(end[0] > 0.9 && end[1] > 0.8, "viridis end should be yellowish");
    }

    #[test]
    fn test_generate_texture() {
        let img = generate_texture(ColorMap::Viridis);
        assert_eq!(img.size().x, 256);
        assert_eq!(img.size().y, 1);
        assert_eq!(img.data.len(), 256 * 4); // RGBA8
    }

    #[test]
    fn test_equal_min_max() {
        // When min == max, should return midpoint color (t=0.5)
        let c = map_value(5.0, 5.0, 5.0, ColorMap::Grayscale);
        assert!((c[0] - 0.5).abs() < EPS, "equal range should give midpoint");
    }

    #[test]
    fn test_all_colormaps_have_32_entries() {
        for &map in ColorMap::ALL {
            let lut = lut_for(map);
            assert_eq!(lut.len(), 32, "{:?} LUT should have 32 entries", map);
        }
    }

    #[test]
    fn test_all_colormaps_rgba_in_range() {
        for &map in ColorMap::ALL {
            let lut = lut_for(map);
            for (i, entry) in lut.iter().enumerate() {
                for (ch, &v) in entry.iter().enumerate() {
                    assert!(
                        (0.0..=1.0).contains(&v),
                        "{:?} LUT[{i}][{ch}] = {v} out of [0,1]",
                        map
                    );
                }
            }
        }
    }

    #[test]
    fn test_encode_rgb_multi_basic() {
        let c = encode_rgb_multi(1.0, 0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.0, 1.0);
        assert!((c[0] - 1.0).abs() < 0.01, "R should be 1.0");
        assert!((c[1] - 0.5).abs() < 0.01, "G should be 0.5");
        assert!((c[2] - 0.0).abs() < 0.01, "B should be 0.0");
        assert_eq!(c[3], 1.0);
    }

    #[test]
    fn test_encode_rgb_multi_zero_range() {
        let c = encode_rgb_multi(5.0, 5.0, 5.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!((c[0] - 0.5).abs() < 0.01, "zero-range should give 0.5");
        assert!((c[1] - 0.5).abs() < 0.01);
        assert!((c[2] - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_hsv_to_rgb_red() {
        // Hue=0 should be red-ish
        let c = hsv_to_rgb(0.0, 1.0, 1.0);
        assert!(c[0] > 0.9, "h=0 should be red, got R={}", c[0]);
        assert!(c[1] < 0.1, "h=0 G should be low");
        assert!(c[2] < 0.1, "h=0 B should be low");
    }

    #[test]
    fn test_hsv_to_rgb_green() {
        // Hue=1/3 should be green-ish
        let c = hsv_to_rgb(1.0 / 3.0, 1.0, 1.0);
        assert!(c[1] > 0.9, "h=1/3 should be green, got G={}", c[1]);
    }

    #[test]
    fn test_hsv_to_rgb_zero_value() {
        // Value=0 should be black
        let c = hsv_to_rgb(0.5, 1.0, 0.0);
        assert!(c[0] < 0.01 && c[1] < 0.01 && c[2] < 0.01);
    }

    #[test]
    fn test_encode_hsv_phase_xy() {
        // Field pointing along +X: atan2(0, 1) = 0
        let c = encode_hsv_phase([1.0, 0.0, 0.0], PhasePlane::XY, 1.0);
        assert!(c[0] + c[1] + c[2] > 0.1, "should have some color");
        assert_eq!(c[3], 1.0);

        // Zero field should be dark (value=0)
        let c = encode_hsv_phase([0.0, 0.0, 0.0], PhasePlane::XY, 1.0);
        assert!(c[0] < 0.1 && c[1] < 0.1 && c[2] < 0.1);
    }

    #[test]
    fn test_encode_hsv_phase_magnitude_scales() {
        let weak = encode_hsv_phase([0.1, 0.0, 0.0], PhasePlane::XY, 1.0);
        let strong = encode_hsv_phase([1.0, 0.0, 0.0], PhasePlane::XY, 1.0);
        let weak_brightness = weak[0] + weak[1] + weak[2];
        let strong_brightness = strong[0] + strong[1] + strong[2];
        assert!(
            strong_brightness > weak_brightness,
            "stronger field should be brighter"
        );
    }
}
