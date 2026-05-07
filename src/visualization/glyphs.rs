// 3D vector field glyph visualization using Bevy gizmos.
//
// Renders arrow glyphs at subsampled grid points showing vector field direction
// and magnitude. Supports E, B, Poynting, and vector potential (A) fields.
// Color and size can be driven by multiple encoding modes.

use bevy::prelude::*;

use crate::math::fdtd;
use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;
use crate::visualization::color_maps::{self, ColorEncoding, ColorMap, PhasePlane};
use crate::visualization::slices::FieldQuantity;

/// Which vector field to visualize as glyphs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GlyphField {
    Electric,
    Magnetic,
    Poynting,
    VectorPotential,
}

impl GlyphField {
    pub const ALL: &'static [GlyphField] = &[
        GlyphField::Electric,
        GlyphField::Magnetic,
        GlyphField::Poynting,
        GlyphField::VectorPotential,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            GlyphField::Electric => "E field",
            GlyphField::Magnetic => "B field",
            GlyphField::Poynting => "Poynting",
            GlyphField::VectorPotential => "A (potential)",
        }
    }
}

/// Configuration for RGB multi-field glyph coloring.
#[derive(Clone, Debug)]
pub struct GlyphRgbConfig {
    pub r_field: FieldQuantity,
    pub g_field: FieldQuantity,
    pub b_field: FieldQuantity,
}

impl Default for GlyphRgbConfig {
    fn default() -> Self {
        Self {
            r_field: FieldQuantity::EMagnitude,
            g_field: FieldQuantity::BMagnitude,
            b_field: FieldQuantity::SField,
        }
    }
}

/// Configuration for HSV phase glyph coloring.
#[derive(Clone, Debug)]
pub struct GlyphHsvConfig {
    pub plane: PhasePlane,
}

impl Default for GlyphHsvConfig {
    fn default() -> Self {
        Self {
            plane: PhasePlane::XY,
        }
    }
}

/// Configuration for Size+Color glyph encoding.
#[derive(Clone, Debug)]
pub struct GlyphSizeColorConfig {
    pub size_field: FieldQuantity,
    pub color_field: FieldQuantity,
    pub color_map: ColorMap,
}

impl Default for GlyphSizeColorConfig {
    fn default() -> Self {
        Self {
            size_field: FieldQuantity::EMagnitude,
            color_field: FieldQuantity::BMagnitude,
            color_map: ColorMap::Viridis,
        }
    }
}

/// Configuration resource for 3D vector field glyph visualization.
#[derive(Resource)]
pub struct GlyphConfig {
    pub enabled: bool,
    pub field: GlyphField,
    /// Sample every Nth cell along each axis.
    pub stride: u32,
    /// Length of an arrow at the auto-range reference magnitude, in units of
    /// `grid.dx` (cell widths). With `scale = 0.5` an arrow at the 95th-
    /// percentile magnitude is half a cell long; rare outliers grow
    /// proportionally larger but stay bounded relative to the grid.
    pub scale: f32,
    /// Color encoding mode.
    pub color_encoding: ColorEncoding,
    /// Color map used when encoding is Standard.
    pub color_map: ColorMap,
    /// Auto-range for color mapping AND length normalization. Uses the
    /// `auto_range_percentile`-th percentile rather than the max so a single
    /// source-cell outlier doesn't squash every other arrow into invisibility.
    pub auto_range: bool,
    /// Percentile (0..1) used as the auto-range reference. Default 0.95
    /// keeps 95% of the field amplitude in the readable range and clips the
    /// top 5% (typically near sources) to a uniform "saturated" bin.
    pub auto_range_percentile: f32,
    pub manual_min: f32,
    pub manual_max: f32,
    /// RGB multi-field config.
    pub rgb_config: GlyphRgbConfig,
    /// HSV phase config.
    pub hsv_config: GlyphHsvConfig,
    /// Size+Color config.
    pub size_color_config: GlyphSizeColorConfig,
}

impl Default for GlyphConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            field: GlyphField::Electric,
            stride: 4,
            // dx-relative: 0.5 cells at the percentile reference.
            scale: 0.5,
            color_encoding: ColorEncoding::Standard,
            color_map: ColorMap::Viridis,
            auto_range: true,
            auto_range_percentile: 0.95,
            manual_min: 0.0,
            manual_max: 1.0,
            rgb_config: GlyphRgbConfig::default(),
            hsv_config: GlyphHsvConfig::default(),
            size_color_config: GlyphSizeColorConfig::default(),
        }
    }
}

/// Compute the `p`-th percentile (0..1) of a slice of magnitudes.
///
/// Sorts a copy of the input and returns the value at the percentile index.
/// Returns 0.0 for empty input. Used by the auto-range path so that strong
/// near-source spikes don't compress the rest of the field into invisibility.
pub fn percentile_of(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f32> = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p = p.clamp(0.0, 1.0);
    let idx = ((p * (v.len() - 1) as f32).round() as usize).min(v.len() - 1);
    v[idx]
}

/// Sample the selected vector field at a given flat grid index.
#[inline]
fn sample_vector_at(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    idx: usize,
    field: GlyphField,
) -> [f32; 3] {
    match field {
        GlyphField::Electric => diag.fields[idx].e,
        GlyphField::Magnetic => diag.fields[idx].b,
        GlyphField::Poynting => diag.fields[idx].poynting,
        GlyphField::VectorPotential => {
            let cell = &grid.read_buf()[idx];
            [cell.q[1], cell.q[2], cell.q[3]]
        }
    }
}

/// Sample a scalar field quantity at a given flat grid index (for encoding).
#[inline]
fn sample_scalar_at(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    idx: usize,
    field: FieldQuantity,
) -> f32 {
    match field {
        FieldQuantity::EMagnitude => {
            let e = diag.fields[idx].e;
            (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt()
        }
        FieldQuantity::BMagnitude => {
            let b = diag.fields[idx].b;
            (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt()
        }
        FieldQuantity::SField => diag.fields[idx].s,
        FieldQuantity::Phi => grid.read_buf()[idx].q[0],
        FieldQuantity::Ax => grid.read_buf()[idx].q[1],
        FieldQuantity::Ay => grid.read_buf()[idx].q[2],
        FieldQuantity::Az => grid.read_buf()[idx].q[3],
        FieldQuantity::EnergyDensity => diag.fields[idx].energy_density,
        FieldQuantity::KVacuum => grid.read_buf()[idx].k,
    }
}

/// Convert grid coordinates to world position (domain centered at origin).
#[inline]
fn grid_to_world(x: usize, y: usize, z: usize, grid: &SimulationGrid) -> Vec3 {
    Vec3::new(
        (x as f32 + 0.5) * grid.dx - grid.nx as f32 * grid.dx * 0.5,
        (y as f32 + 0.5) * grid.dx - grid.ny as f32 * grid.dx * 0.5,
        (z as f32 + 0.5) * grid.dx - grid.nz as f32 * grid.dx * 0.5,
    )
}

/// Bevy system: draw vector field glyphs using gizmos.
pub fn draw_glyph_arrows(
    mut gizmos: Gizmos,
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    config: Res<GlyphConfig>,
) {
    if !config.enabled {
        return;
    }
    let Some(grid) = grid else { return };
    if diag.fields.is_empty() {
        return;
    }

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let stride = config.stride.max(1) as usize;

    // First pass: collect samples, vectors, and magnitudes for auto-ranging
    let mut samples: Vec<(usize, usize, usize, usize, [f32; 3], f32)> = Vec::new();

    for z in (0..nz).step_by(stride) {
        for y in (0..ny).step_by(stride) {
            for x in (0..nx).step_by(stride) {
                let idx = fdtd::idx(x, y, z, nx, ny);
                let vec = sample_vector_at(&grid, &diag, idx, config.field);
                let mag = (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt();
                samples.push((x, y, z, idx, vec, mag));
            }
        }
    }

    // Auto-range uses the configurable percentile (default 95th) of the
    // sample magnitudes, NOT the absolute max. With `max_mag = max(samples)`
    // a single near-source outlier (where |E| is orders of magnitude above
    // the rest of the field) compresses every other arrow below the visibility
    // threshold and the visualization looks empty. Percentile clipping keeps
    // the bulk of the field readable while letting strong spots saturate.
    let max_mag = if config.auto_range {
        let mags: Vec<f32> = samples.iter().map(|&(_, _, _, _, _, m)| m).collect();
        percentile_of(&mags, config.auto_range_percentile).max(1e-6)
    } else {
        config.manual_max.max(1e-6)
    };
    let min_mag = if config.auto_range {
        0.0
    } else {
        config.manual_min
    };

    // For RgbMultiField and SizeColor we encode SCALAR fields whose magnitudes
    // are unrelated to the vector-field max. Reusing `max_mag` (a vector
    // magnitude) for these scalars caused channels to saturate or go black —
    // a real correctness bug. Compute per-channel percentile maxes from the
    // same sample set, lazily, only for the encodings that need them.
    let scalar_p95 = |field: FieldQuantity| -> f32 {
        let vals: Vec<f32> = samples
            .iter()
            .map(|&(_, _, _, idx, _, _)| sample_scalar_at(&grid, &diag, idx, field).abs())
            .collect();
        percentile_of(&vals, config.auto_range_percentile).max(1e-6)
    };
    let (r_max, g_max, b_max) = match config.color_encoding {
        ColorEncoding::RgbMultiField => (
            scalar_p95(config.rgb_config.r_field),
            scalar_p95(config.rgb_config.g_field),
            scalar_p95(config.rgb_config.b_field),
        ),
        _ => (1.0, 1.0, 1.0),
    };
    let (size_max, color_max) = match config.color_encoding {
        ColorEncoding::SizeColor => (
            scalar_p95(config.size_color_config.size_field),
            scalar_p95(config.size_color_config.color_field),
        ),
        _ => (1.0, 1.0),
    };

    // Arrow lengths are normalised to grid.dx so the visualization is unit-
    // agnostic: `config.scale` directly means "cell widths at the reference
    // magnitude". An arrow at the percentile reference is `scale * dx` long;
    // rare outliers grow proportionally larger.
    let target_unit = grid.dx * config.scale;

    // Second pass: draw arrows
    for &(x, y, z, idx, vec, mag) in &samples {
        if mag < max_mag * 1e-4 {
            continue;
        }

        let dir = Vec3::new(vec[0], vec[1], vec[2]) / mag;
        let pos = grid_to_world(x, y, z, &grid);

        let (arrow_length, color) = match config.color_encoding {
            ColorEncoding::Standard => {
                let length = (mag / max_mag) * target_unit;
                let rgba = color_maps::map_value(mag, min_mag, max_mag, config.color_map);
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
            ColorEncoding::RgbMultiField => {
                let length = (mag / max_mag) * target_unit;
                let r_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.r_field);
                let g_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.g_field);
                let b_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.b_field);
                let rgba = color_maps::encode_rgb_multi(
                    r_val, 0.0, r_max, g_val, 0.0, g_max, b_val, 0.0, b_max,
                );
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
            ColorEncoding::HsvPhase => {
                let length = (mag / max_mag) * target_unit;
                let rgba =
                    color_maps::encode_hsv_phase(vec, config.hsv_config.plane, max_mag);
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
            ColorEncoding::SizeColor => {
                let size_val =
                    sample_scalar_at(&grid, &diag, idx, config.size_color_config.size_field);
                let color_val = sample_scalar_at(
                    &grid,
                    &diag,
                    idx,
                    config.size_color_config.color_field,
                );
                let length = (size_val.abs() / size_max) * target_unit;
                let rgba = color_maps::map_value(
                    color_val.abs(),
                    0.0,
                    color_max,
                    config.size_color_config.color_map,
                );
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
        };

        if arrow_length < 1e-10 {
            continue;
        }

        let end = pos + dir * arrow_length;

        // Draw arrow shaft
        gizmos.line(pos, end, color);

        // Draw arrowhead
        let head_len = arrow_length * 0.25;
        let perp = if dir.dot(Vec3::Y).abs() < 0.9 {
            dir.cross(Vec3::Y).normalize()
        } else {
            dir.cross(Vec3::X).normalize()
        };
        gizmos.line(end, end - dir * head_len + perp * head_len * 0.3, color);
        gizmos.line(end, end - dir * head_len - perp * head_len * 0.3, color);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_percentile_of_basic() {
        let v: Vec<f32> = (1..=100).map(|i| i as f32).collect();
        // 50th percentile of 1..=100 (101 values? no, 100 values) is the
        // value at idx = round(0.5 * 99) = 50, i.e. v[50] = 51.0.
        assert!((percentile_of(&v, 0.5) - 51.0).abs() < 0.5);
        // 95th percentile: idx = round(0.95 * 99) = 94, v[94] = 95.0.
        assert!((percentile_of(&v, 0.95) - 95.0).abs() < 0.5);
        // 100th percentile = max = 100.0.
        assert_eq!(percentile_of(&v, 1.0), 100.0);
        // 0th percentile = min = 1.0.
        assert_eq!(percentile_of(&v, 0.0), 1.0);
    }

    #[test]
    fn test_percentile_of_empty_returns_zero() {
        assert_eq!(percentile_of(&[], 0.95), 0.0);
    }

    #[test]
    fn test_percentile_clipping_suppresses_outlier() {
        // 99 small values and 1 huge spike (mimics a near-source cell next
        // to a sea of bulk-radiation cells in the dipole scenario).
        let mut v: Vec<f32> = vec![1.0; 99];
        v.push(1.0e6);
        // The simple max would give 1e6, dominating the auto-range.
        let max = v.iter().cloned().fold(0.0f32, f32::max);
        assert_eq!(max, 1.0e6);
        // The 95th percentile sits in the bulk → 1.0, not the spike.
        let p95 = percentile_of(&v, 0.95);
        assert!(p95 < 10.0, "p95 should be bulk, got {p95}");
    }
}
