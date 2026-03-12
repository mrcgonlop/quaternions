// Isosurface extraction and rendering.
//
// Extracts a scalar isosurface as a set of boundary faces between voxels that
// are above and below the isovalue, rendered as wireframe quads via Bevy gizmos.
//
// This is a voxelized (staircase) approach: for each pair of adjacent cells
// where one is inside (>= isovalue) and the other is outside (< isovalue), the
// shared face is drawn as four line segments. The result resembles Minecraft-style
// voxel boundaries — not smooth, but always correct and requires no lookup tables.
//
// Use `voxel_stride` > 1 to subsample the grid for better performance.

use bevy::prelude::*;

use crate::math::fdtd;
use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;

/// Scalar field to extract the isosurface from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IsoField {
    EMagnitude,
    BMagnitude,
    SField,
    EnergyDensity,
    Phi,
    KVacuum,
}

impl IsoField {
    pub const ALL: &'static [IsoField] = &[
        IsoField::EMagnitude,
        IsoField::BMagnitude,
        IsoField::SField,
        IsoField::EnergyDensity,
        IsoField::Phi,
        IsoField::KVacuum,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            IsoField::EMagnitude => "|E|",
            IsoField::BMagnitude => "|B|",
            IsoField::SField => "S field",
            IsoField::EnergyDensity => "Energy density",
            IsoField::Phi => "phi/c",
            IsoField::KVacuum => "K (vacuum index)",
        }
    }
}

/// Configuration for isosurface visualization.
#[derive(Resource)]
pub struct IsosurfaceConfig {
    pub enabled: bool,
    /// Scalar field to extract the isosurface from.
    pub field: IsoField,
    /// Isosurface threshold value; cells with scalar >= isovalue are "inside".
    pub isovalue: f32,
    /// Only check every `voxel_stride`-th cell in each axis (performance control).
    pub voxel_stride: u32,
    /// Color of the wireframe faces.
    pub color: Color,
    /// If true, auto-set isovalue to 50% of the field maximum each frame.
    pub auto_isovalue: bool,
}

impl Default for IsosurfaceConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            field: IsoField::EMagnitude,
            isovalue: 1.0,
            voxel_stride: 2,
            color: Color::srgba(0.2, 0.8, 1.0, 0.7),
            auto_isovalue: false,
        }
    }
}

/// Sample a scalar field at a flat grid index.
#[inline]
fn sample_scalar(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    idx: usize,
    field: IsoField,
) -> f32 {
    match field {
        IsoField::EMagnitude => {
            let e = diag.fields[idx].e;
            (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt()
        }
        IsoField::BMagnitude => {
            let b = diag.fields[idx].b;
            (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt()
        }
        IsoField::SField => diag.fields[idx].s,
        IsoField::EnergyDensity => diag.fields[idx].energy_density,
        IsoField::Phi => grid.read_buf()[idx].q[0],
        IsoField::KVacuum => grid.read_buf()[idx].k,
    }
}

/// Convert cell-corner coordinates to world position.
///
/// Cell corners are at multiples of dx. Cell center (ix, iy, iz) is bounded by
/// corners (ix, iy, iz) and (ix+1, iy+1, iz+1). Corner ix maps to:
///   world_x = ix * dx - nx * dx / 2
#[inline]
fn corner_to_world(ix: usize, iy: usize, iz: usize, grid: &SimulationGrid) -> Vec3 {
    Vec3::new(
        ix as f32 * grid.dx - grid.nx as f32 * grid.dx * 0.5,
        iy as f32 * grid.dx - grid.ny as f32 * grid.dx * 0.5,
        iz as f32 * grid.dx - grid.nz as f32 * grid.dx * 0.5,
    )
}

/// Draw a quad wireframe (4 edges) for a boundary face.
#[inline]
fn draw_face(gizmos: &mut Gizmos, a: Vec3, b: Vec3, c: Vec3, d: Vec3, color: Color) {
    gizmos.line(a, b, color);
    gizmos.line(b, c, color);
    gizmos.line(c, d, color);
    gizmos.line(d, a, color);
}

/// Bevy system: extract and draw the isosurface as boundary face wireframes.
pub fn draw_isosurface(
    mut gizmos: Gizmos,
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    config: Res<IsosurfaceConfig>,
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
    let stride = config.voxel_stride.max(1) as usize;

    // Determine isovalue (auto = 50% of field maximum).
    let isovalue = if config.auto_isovalue {
        let mut max_val = 1e-6f32;
        for z in (0..nz).step_by(stride) {
            for y in (0..ny).step_by(stride) {
                for x in (0..nx).step_by(stride) {
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    let v = sample_scalar(&grid, &diag, idx, config.field);
                    max_val = max_val.max(v.abs());
                }
            }
        }
        max_val * 0.5
    } else {
        config.isovalue
    };

    let color = config.color;

    // Scan every cell at the stride resolution. For each of the three axis
    // directions, check if this cell and its neighbor differ in inside/outside
    // status. If so, draw the shared boundary face as a wireframe quad.
    for z in (0..nz).step_by(stride) {
        for y in (0..ny).step_by(stride) {
            for x in (0..nx).step_by(stride) {
                let i0 = fdtd::idx(x, y, z, nx, ny);
                let v0 = sample_scalar(&grid, &diag, i0, config.field);
                let inside0 = v0 >= isovalue;

                // +X face: between (x, y, z) and (x+1, y, z).
                // Face corners in cell-corner space: (x+1, y..y+1, z..z+1).
                if x + 1 < nx {
                    let i1 = fdtd::idx(x + 1, y, z, nx, ny);
                    let inside1 = sample_scalar(&grid, &diag, i1, config.field) >= isovalue;
                    if inside0 != inside1 {
                        let a = corner_to_world(x + 1, y, z, &grid);
                        let b = corner_to_world(x + 1, y + 1, z, &grid);
                        let c = corner_to_world(x + 1, y + 1, z + 1, &grid);
                        let d = corner_to_world(x + 1, y, z + 1, &grid);
                        draw_face(&mut gizmos, a, b, c, d, color);
                    }
                }

                // +Y face: between (x, y, z) and (x, y+1, z).
                // Face corners in cell-corner space: (x..x+1, y+1, z..z+1).
                if y + 1 < ny {
                    let i1 = fdtd::idx(x, y + 1, z, nx, ny);
                    let inside1 = sample_scalar(&grid, &diag, i1, config.field) >= isovalue;
                    if inside0 != inside1 {
                        let a = corner_to_world(x, y + 1, z, &grid);
                        let b = corner_to_world(x + 1, y + 1, z, &grid);
                        let c = corner_to_world(x + 1, y + 1, z + 1, &grid);
                        let d = corner_to_world(x, y + 1, z + 1, &grid);
                        draw_face(&mut gizmos, a, b, c, d, color);
                    }
                }

                // +Z face: between (x, y, z) and (x, y, z+1).
                // Face corners in cell-corner space: (x..x+1, y..y+1, z+1).
                if z + 1 < nz {
                    let i1 = fdtd::idx(x, y, z + 1, nx, ny);
                    let inside1 = sample_scalar(&grid, &diag, i1, config.field) >= isovalue;
                    if inside0 != inside1 {
                        let a = corner_to_world(x, y, z + 1, &grid);
                        let b = corner_to_world(x + 1, y, z + 1, &grid);
                        let c = corner_to_world(x + 1, y + 1, z + 1, &grid);
                        let d = corner_to_world(x, y + 1, z + 1, &grid);
                        draw_face(&mut gizmos, a, b, c, d, color);
                    }
                }
            }
        }
    }
}
