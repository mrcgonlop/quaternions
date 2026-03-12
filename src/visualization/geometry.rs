// Source and conductor geometry rendering.
//
// Draws visual markers at grid cells that are flagged as SOURCE, CONDUCTOR, or
// BOUNDARY. Uses Bevy gizmos — no mesh spawning required.
//
// SOURCE cells:    3-axis cross + wireframe box outline (yellow-gold).
// CONDUCTOR cells: wireframe box outline (lavender).
// BOUNDARY cells:  small wireframe box (gray, disabled by default).

use bevy::prelude::*;

use crate::math::fdtd;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::CellFlags;

/// Configuration for geometry visualization.
#[derive(Resource)]
pub struct GeometryConfig {
    pub enabled: bool,
    /// Draw markers at SOURCE-flagged cells.
    pub show_sources: bool,
    /// Draw outlines at CONDUCTOR-flagged cells.
    pub show_conductors: bool,
    /// Draw outlines at BOUNDARY-flagged cells (noisy; off by default).
    pub show_boundary: bool,
    pub source_color: Color,
    pub conductor_color: Color,
    pub boundary_color: Color,
}

impl Default for GeometryConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            show_sources: true,
            show_conductors: true,
            show_boundary: false,
            source_color: Color::srgb(1.0, 0.8, 0.0),
            conductor_color: Color::srgb(0.6, 0.6, 0.9),
            boundary_color: Color::srgba(0.4, 0.4, 0.4, 0.3),
        }
    }
}

/// Convert grid cell indices to world position (cell center, domain centered at origin).
#[inline]
fn grid_to_world(x: usize, y: usize, z: usize, grid: &SimulationGrid) -> Vec3 {
    Vec3::new(
        (x as f32 + 0.5) * grid.dx - grid.nx as f32 * grid.dx * 0.5,
        (y as f32 + 0.5) * grid.dx - grid.ny as f32 * grid.dx * 0.5,
        (z as f32 + 0.5) * grid.dx - grid.nz as f32 * grid.dx * 0.5,
    )
}

/// Draw a wireframe box (12 edges) centered at `center` with side length `size`.
fn draw_wire_box(gizmos: &mut Gizmos, center: Vec3, size: f32, color: Color) {
    let h = size * 0.5;
    let corners = [
        center + Vec3::new(-h, -h, -h),
        center + Vec3::new(h, -h, -h),
        center + Vec3::new(h, h, -h),
        center + Vec3::new(-h, h, -h),
        center + Vec3::new(-h, -h, h),
        center + Vec3::new(h, -h, h),
        center + Vec3::new(h, h, h),
        center + Vec3::new(-h, h, h),
    ];
    let edges: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    for (a, b) in edges {
        gizmos.line(corners[a], corners[b], color);
    }
}

/// Draw a 3-axis cross marker centered at `center` with arm length `size`.
fn draw_cross(gizmos: &mut Gizmos, center: Vec3, size: f32, color: Color) {
    let h = size * 0.5;
    gizmos.line(center - Vec3::X * h, center + Vec3::X * h, color);
    gizmos.line(center - Vec3::Y * h, center + Vec3::Y * h, color);
    gizmos.line(center - Vec3::Z * h, center + Vec3::Z * h, color);
}

/// Bevy system: draw source, conductor, and boundary cell markers.
pub fn draw_geometry(
    mut gizmos: Gizmos,
    grid: Option<Res<SimulationGrid>>,
    config: Res<GeometryConfig>,
) {
    if !config.enabled {
        return;
    }
    let Some(grid) = grid else { return };

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;
    let nz = grid.nz as usize;
    let dx = grid.dx;
    let buf = grid.read_buf();

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = fdtd::idx(x, y, z, nx, ny);
                let flags = buf[i].flags;

                if config.show_sources && (flags & CellFlags::SOURCE) != 0 {
                    let pos = grid_to_world(x, y, z, &grid);
                    draw_cross(&mut gizmos, pos, dx * 1.5, config.source_color);
                    draw_wire_box(&mut gizmos, pos, dx, config.source_color);
                }

                if config.show_conductors && (flags & CellFlags::CONDUCTOR) != 0 {
                    let pos = grid_to_world(x, y, z, &grid);
                    draw_wire_box(&mut gizmos, pos, dx, config.conductor_color);
                }

                if config.show_boundary && (flags & CellFlags::BOUNDARY) != 0 {
                    let pos = grid_to_world(x, y, z, &grid);
                    draw_wire_box(&mut gizmos, pos, dx * 0.5, config.boundary_color);
                }
            }
        }
    }
}
