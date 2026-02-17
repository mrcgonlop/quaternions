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
    /// Scale factor for arrow length (world units per field unit).
    pub scale: f32,
    /// Color encoding mode.
    pub color_encoding: ColorEncoding,
    /// Color map used when encoding is Standard.
    pub color_map: ColorMap,
    /// Auto-range for color mapping.
    pub auto_range: bool,
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
            scale: 0.001,
            color_encoding: ColorEncoding::Standard,
            color_map: ColorMap::Viridis,
            auto_range: true,
            manual_min: 0.0,
            manual_max: 1.0,
            rgb_config: GlyphRgbConfig::default(),
            hsv_config: GlyphHsvConfig::default(),
            size_color_config: GlyphSizeColorConfig::default(),
        }
    }
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

    // First pass: collect samples and magnitudes for auto-ranging
    let mut samples: Vec<(usize, usize, usize, usize, f32)> = Vec::new();

    for z in (0..nz).step_by(stride) {
        for y in (0..ny).step_by(stride) {
            for x in (0..nx).step_by(stride) {
                let idx = fdtd::idx(x, y, z, nx, ny);
                let vec = sample_vector_at(&grid, &diag, idx, config.field);
                let mag = (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt();
                samples.push((x, y, z, idx, mag));
            }
        }
    }

    let max_mag = if config.auto_range {
        samples
            .iter()
            .map(|&(_, _, _, _, m)| m)
            .fold(0.0f32, f32::max)
            .max(1e-6)
    } else {
        config.manual_max.max(1e-6)
    };
    let min_mag = if config.auto_range {
        0.0
    } else {
        config.manual_min
    };

    // Second pass: draw arrows
    for &(x, y, z, idx, mag) in &samples {
        if mag < max_mag * 1e-4 {
            continue;
        }

        let vec = sample_vector_at(&grid, &diag, idx, config.field);
        let dir = Vec3::new(vec[0], vec[1], vec[2]) / mag;
        let pos = grid_to_world(x, y, z, &grid);

        let (arrow_length, color) = match config.color_encoding {
            ColorEncoding::Standard => {
                let length = mag * config.scale;
                let rgba = color_maps::map_value(mag, min_mag, max_mag, config.color_map);
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
            ColorEncoding::RgbMultiField => {
                let length = mag * config.scale;
                let r_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.r_field);
                let g_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.g_field);
                let b_val = sample_scalar_at(&grid, &diag, idx, config.rgb_config.b_field);
                let rgba = color_maps::encode_rgb_multi(
                    r_val, 0.0, max_mag, g_val, 0.0, max_mag, b_val, 0.0, max_mag,
                );
                (length, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]))
            }
            ColorEncoding::HsvPhase => {
                let length = mag * config.scale;
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
                let length = size_val.abs() * config.scale;
                let rgba = color_maps::map_value(
                    color_val.abs(),
                    0.0,
                    max_mag,
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
