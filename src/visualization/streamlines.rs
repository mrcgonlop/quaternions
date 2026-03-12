// Streamline computation and rendering.
//
// Integrates vector field lines from uniformly seeded starting points using
// Euler integration with trilinear interpolation. Each streamline follows the
// chosen vector field until it leaves the domain or exhausts max_steps.
//
// Rendered as connected line segments via Bevy gizmos. Color encodes local
// field magnitude using a configurable color map.

use bevy::prelude::*;

use crate::math::fdtd;
use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;
use crate::visualization::color_maps::{self, ColorMap};
use crate::visualization::glyphs::GlyphField;

/// Configuration for streamline visualization.
#[derive(Resource)]
pub struct StreamlineConfig {
    pub enabled: bool,
    /// Which vector field to trace streamlines along.
    pub field: GlyphField,
    /// Seed one streamline every `seed_stride` cells along each axis.
    pub seed_stride: u32,
    /// Maximum number of Euler integration steps per streamline.
    pub max_steps: u32,
    /// Step size as a fraction of grid spacing dx.
    pub step_fraction: f32,
    /// Color map for field magnitude encoding.
    pub color_map: ColorMap,
    /// If true, auto-scale the color range to the maximum field magnitude.
    pub auto_range: bool,
    /// Manual maximum for color range (used when auto_range is false).
    pub manual_max: f32,
}

impl Default for StreamlineConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            field: GlyphField::Electric,
            seed_stride: 6,
            max_steps: 64,
            step_fraction: 0.5,
            color_map: ColorMap::Plasma,
            auto_range: true,
            manual_max: 1.0,
        }
    }
}

/// Convert fractional grid coordinates to world position.
///
/// fx=0 corresponds to the center of cell (0,*,*).
#[inline]
fn frac_to_world(fx: f32, fy: f32, fz: f32, grid: &SimulationGrid) -> Vec3 {
    Vec3::new(
        (fx + 0.5) * grid.dx - grid.nx as f32 * grid.dx * 0.5,
        (fy + 0.5) * grid.dx - grid.ny as f32 * grid.dx * 0.5,
        (fz + 0.5) * grid.dx - grid.nz as f32 * grid.dx * 0.5,
    )
}

/// Trilinearly interpolate a vector field at fractional grid coordinates.
///
/// Returns [0, 0, 0] when both ix and ix+1 (and equivalents for y, z) are not
/// within the grid — i.e., outside the safe interpolation region.
fn trilinear_sample(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    field: GlyphField,
    fx: f32,
    fy: f32,
    fz: f32,
) -> [f32; 3] {
    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;

    if fx < 0.0
        || fx >= (nx - 1) as f32
        || fy < 0.0
        || fy >= (ny - 1) as f32
        || fz < 0.0
        || fz >= (nz - 1) as f32
    {
        return [0.0; 3];
    }

    let ix = fx as usize;
    let iy = fy as usize;
    let iz = fz as usize;
    let tx = fx - ix as f32;
    let ty = fy - iy as f32;
    let tz = fz - iz as f32;

    let w = [
        (1.0 - tx) * (1.0 - ty) * (1.0 - tz),
        tx * (1.0 - ty) * (1.0 - tz),
        (1.0 - tx) * ty * (1.0 - tz),
        tx * ty * (1.0 - tz),
        (1.0 - tx) * (1.0 - ty) * tz,
        tx * (1.0 - ty) * tz,
        (1.0 - tx) * ty * tz,
        tx * ty * tz,
    ];
    let cx = [ix, ix + 1, ix, ix + 1, ix, ix + 1, ix, ix + 1];
    let cy = [iy, iy, iy + 1, iy + 1, iy, iy, iy + 1, iy + 1];
    let cz = [iz, iz, iz, iz, iz + 1, iz + 1, iz + 1, iz + 1];

    let mut result = [0.0f32; 3];
    for i in 0..8 {
        let idx = fdtd::idx(cx[i], cy[i], cz[i], nx, ny);
        let v = sample_vec_at(grid, diag, idx, field);
        result[0] += w[i] * v[0];
        result[1] += w[i] * v[1];
        result[2] += w[i] * v[2];
    }
    result
}

/// Sample a vector field at a flat grid index.
#[inline]
fn sample_vec_at(
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

/// Bevy system: trace streamlines from subsampled seed points and draw with gizmos.
pub fn draw_streamlines(
    mut gizmos: Gizmos,
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    config: Res<StreamlineConfig>,
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
    let stride = config.seed_stride.max(1) as usize;
    let step_fraction = config.step_fraction.max(0.01);

    // Auto-range: scan sampled cells for maximum field magnitude.
    let max_mag = if config.auto_range {
        let mut m = 1e-6f32;
        for z in (0..nz).step_by(stride) {
            for y in (0..ny).step_by(stride) {
                for x in (0..nx).step_by(stride) {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    let v = sample_vec_at(&grid, &diag, idx, config.field);
                    let mag2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                    m = m.max(mag2.sqrt());
                }
            }
        }
        m
    } else {
        config.manual_max.max(1e-6)
    };

    let threshold = max_mag * 1e-4;

    // Trace one streamline from each interior seed point.
    for sz in (1..nz.saturating_sub(1)).step_by(stride) {
        for sy in (1..ny.saturating_sub(1)).step_by(stride) {
            for sx in (1..nx.saturating_sub(1)).step_by(stride) {
                let mut fx = sx as f32;
                let mut fy = sy as f32;
                let mut fz = sz as f32;
                let mut prev = frac_to_world(fx, fy, fz, &grid);

                for _ in 0..config.max_steps {
                    let v = trilinear_sample(&grid, &diag, config.field, fx, fy, fz);
                    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    if mag < threshold {
                        break;
                    }

                    let inv_mag = 1.0 / mag;
                    fx += v[0] * inv_mag * step_fraction;
                    fy += v[1] * inv_mag * step_fraction;
                    fz += v[2] * inv_mag * step_fraction;

                    // Stop if we've left the safe interpolation region (2 cells from edge).
                    if fx < 1.0
                        || fx >= (nx - 2) as f32
                        || fy < 1.0
                        || fy >= (ny - 2) as f32
                        || fz < 1.0
                        || fz >= (nz - 2) as f32
                    {
                        break;
                    }

                    let cur = frac_to_world(fx, fy, fz, &grid);
                    let rgba = color_maps::map_value(mag, 0.0, max_mag, config.color_map);
                    gizmos.line(
                        prev,
                        cur,
                        Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]),
                    );
                    prev = cur;
                }
            }
        }
    }
}
