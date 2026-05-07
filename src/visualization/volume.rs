// Volume rendering — primary 3D scalar view (Phase 6.1, CPU prototype).
//
// Renders a chosen scalar field as a stack of N axis-aligned semi-transparent
// quads (slabs). Each slab is a 2D textured quad at a fixed Z position; the
// texture is the field's value on that Z layer mapped through a transfer
// function. Bevy's standard alpha-blending pipeline composites the slabs back-
// to-front automatically when the camera looks roughly down ±Z.
//
// Two transfer functions, picked automatically from the field's signedness
// (see `FieldQuantity::is_signed`):
//   * **Bipolar** for signed scalars (S, φ, A_x, A_y, A_z): blue at −max,
//     transparent at zero, red at +max. Opacity ramps with |value|/max.
//   * **Sequential** for positive scalars (|E|, |B|, K, energy density): a
//     selectable perceptual palette (default viridis); opacity ramps from
//     0 at zero to opaque near max.
//
// Known limitation (acceptable for the CPU prototype): when the camera is
// nearly edge-on to the slab stack (line-of-sight in the xy-plane), the
// slabs collapse to lines and the volume effectively disappears. View-
// aligned slabs (re-orienting per-frame to face the camera) is a polish
// task tracked in Phase 6.1 GPU-port notes.

use bevy::prelude::*;
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::state::CellFlags;
use crate::visualization::color_maps::{self, ColorMap};
use crate::visualization::slices::FieldQuantity;

/// Per-render-frame configuration for the volume view.
#[derive(Resource)]
pub struct VolumeConfig {
    pub enabled: bool,
    /// Scalar quantity to render. Sign determines the transfer function.
    pub field: FieldQuantity,
    /// Number of slabs along Z. Fewer = faster, more aliasing; more = slower,
    /// smoother. Capped at the grid's `nz` so we never duplicate slabs.
    pub num_slabs: u32,
    /// Multiplier on per-pixel alpha. Lower for "see-through" rendering of
    /// dense volumes; higher to make weak fields more visible.
    pub opacity_scale: f32,
    /// Auto-scale the transfer function range to the live field max.
    pub auto_range: bool,
    /// Manual range upper bound (used when `auto_range = false`). For signed
    /// fields the bipolar transfer function uses [-manual_max, +manual_max].
    pub manual_max: f32,
    /// Color map used by the SEQUENTIAL transfer function (positive fields).
    /// Bipolar uses a hand-coded blue/white/red mapping that ignores this.
    pub color_map: ColorMap,
    /// Skip PML cells when computing the auto-range — strong CPML auxiliary
    /// fields can otherwise dominate the colour scale.
    pub exclude_pml_from_range: bool,
    // Internal — track whether the entity stack needs to be respawned.
    last_enabled: bool,
    last_num_slabs: u32,
}

impl Default for VolumeConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            field: FieldQuantity::EMagnitude,
            num_slabs: 16,
            opacity_scale: 0.6,
            auto_range: true,
            manual_max: 1.0,
            color_map: ColorMap::Viridis,
            exclude_pml_from_range: true,
            last_enabled: false,
            last_num_slabs: 0,
        }
    }
}

/// Marker component on volume slab entities. Stores the slab's Z grid index
/// so the texture-update system knows which Z layer to sample.
#[derive(Component, Clone, Copy)]
pub struct VolumeSlab {
    pub z_index: u32,
}

/// Holds the dynamic image handle for a slab so the update system can swap it.
#[derive(Component)]
pub struct VolumeSlabTexture(pub Handle<Image>);

// ---------------------------------------------------------------------------
// Transfer functions — pure helpers, fully testable.
// ---------------------------------------------------------------------------

/// Bipolar transfer function for signed scalars. Maps `value ∈ [−max, +max]`
/// to RGBA. Negative values shade toward blue, positive toward red, with the
/// linear opacity `|value|/max · opacity_scale`. Returns RGBA in 0..=1.
pub fn transfer_bipolar(value: f32, max: f32, opacity_scale: f32) -> [f32; 4] {
    if max <= 0.0 {
        return [0.0, 0.0, 0.0, 0.0];
    }
    let t = (value / max).clamp(-1.0, 1.0);
    let alpha = (t.abs() * opacity_scale).clamp(0.0, 1.0);
    if t >= 0.0 {
        // White → red as t goes 0 → 1.
        let r = 1.0;
        let g = 1.0 - t;
        let b = 1.0 - t;
        [r, g, b, alpha]
    } else {
        // White → blue as t goes 0 → −1.
        let s = -t;
        let r = 1.0 - s;
        let g = 1.0 - s;
        let b = 1.0;
        [r, g, b, alpha]
    }
}

/// Sequential transfer function for positive scalars. Maps `value ∈ [0, max]`
/// through the supplied colour map; opacity ramps with the magnitude.
pub fn transfer_sequential(
    value: f32,
    max: f32,
    opacity_scale: f32,
    color_map: ColorMap,
) -> [f32; 4] {
    if max <= 0.0 {
        return [0.0, 0.0, 0.0, 0.0];
    }
    let v = value.max(0.0);
    let rgba = color_maps::map_value(v, 0.0, max, color_map);
    let t = (v / max).clamp(0.0, 1.0);
    let alpha = (t * opacity_scale).clamp(0.0, 1.0);
    [rgba[0], rgba[1], rgba[2], alpha]
}

/// Pick the transfer function and produce RGBA for one sample.
fn transfer_for(
    field: FieldQuantity,
    value: f32,
    max: f32,
    opacity_scale: f32,
    color_map: ColorMap,
) -> [f32; 4] {
    if field.is_signed() {
        transfer_bipolar(value, max, opacity_scale)
    } else {
        transfer_sequential(value, max, opacity_scale, color_map)
    }
}

// ---------------------------------------------------------------------------
// Field sampling — duplicated locally to avoid leaking the `slices` private
// helpers. Kept in sync with `slices::sample_field_at`.
// ---------------------------------------------------------------------------

fn sample_field_at(
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

/// Compute the auto-range bound (positive scalar) for the chosen field over
/// every interior cell of the grid, optionally excluding PML cells.
fn compute_max(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    field: FieldQuantity,
    exclude_pml: bool,
) -> f32 {
    let cells = grid.read_buf();
    let mut m = 0.0f32;
    for (i, cell) in cells.iter().enumerate() {
        if exclude_pml && (cell.flags & CellFlags::PML) != 0 {
            continue;
        }
        let v = sample_field_at(grid, diag, i, field).abs();
        if v > m {
            m = v;
        }
    }
    m.max(1e-6)
}

// ---------------------------------------------------------------------------
// Slab geometry / placement
// ---------------------------------------------------------------------------

/// World-Z position of slab `i` of `n` slabs in the volume.
fn slab_world_z(grid: &SimulationGrid, z_index: u32) -> f32 {
    (z_index as f32 + 0.5) * grid.dx - grid.nz as f32 * grid.dx * 0.5
}

/// Choose `n` slab Z indices spread across `[0, nz)` so the volume is
/// covered approximately uniformly. For `n >= nz` returns one slab per
/// layer; for `n < nz` picks evenly spaced layers.
fn slab_z_indices(nz: u32, n: u32) -> Vec<u32> {
    let n = n.max(1).min(nz.max(1));
    if n == 1 {
        return vec![nz / 2];
    }
    let mut out = Vec::with_capacity(n as usize);
    for i in 0..n {
        // Even distribution including both endpoints when n > 1.
        let z = ((i as f32 / (n - 1) as f32) * (nz - 1) as f32).round() as u32;
        out.push(z.min(nz - 1));
    }
    // Ensure no duplicates if nz is small relative to n.
    out.dedup();
    out
}

// ---------------------------------------------------------------------------
// Bevy systems
// ---------------------------------------------------------------------------

/// Spawn slab entities when enabled toggles on or `num_slabs` changes;
/// despawn when disabled.
pub fn manage_volume_entities(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut images: ResMut<Assets<Image>>,
    mut config: ResMut<VolumeConfig>,
    grid: Option<Res<SimulationGrid>>,
    existing: Query<Entity, With<VolumeSlab>>,
) {
    let Some(grid) = grid else { return };

    let needs_respawn =
        config.enabled != config.last_enabled || config.num_slabs != config.last_num_slabs;
    if !needs_respawn {
        return;
    }

    // Tear down any existing slabs.
    for entity in existing.iter() {
        commands.entity(entity).despawn_recursive();
    }

    config.last_enabled = config.enabled;
    config.last_num_slabs = config.num_slabs;

    if !config.enabled {
        return;
    }

    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let z_indices = slab_z_indices(nz, config.num_slabs);
    let pixel_count = (nx * ny) as usize;
    let placeholder = vec![0u8; pixel_count * 4];

    // Slab quad spans the whole xy domain.
    let mesh = meshes.add(Rectangle::new(
        nx as f32 * grid.dx,
        ny as f32 * grid.dx,
    ));

    for &z in &z_indices {
        let image = Image::new(
            Extent3d {
                width: nx,
                height: ny,
                depth_or_array_layers: 1,
            },
            TextureDimension::D2,
            placeholder.clone(),
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
        );
        let image_handle = images.add(image);
        let material = materials.add(StandardMaterial {
            base_color_texture: Some(image_handle.clone()),
            unlit: true,
            alpha_mode: AlphaMode::Blend,
            double_sided: true,
            cull_mode: None,
            ..default()
        });

        let zw = slab_world_z(&grid, z);
        // Bevy's `Rectangle` mesh sits in the XY plane (normal +Z) by
        // default. That's exactly what we want for a Z-stacked slab: the
        // quad lives at (x, y, zw) for x ∈ [−nx·dx/2, +nx·dx/2] and y
        // similarly. NO rotation — applying one (e.g., the previous
        // `from_rotation_x(π/2)`) flips every slab into the XZ plane,
        // making them all overlap and producing the "single slab" look.
        let transform = Transform::from_translation(Vec3::new(0.0, 0.0, zw));

        commands.spawn((
            Mesh3d(mesh.clone()),
            MeshMaterial3d(material),
            transform,
            VolumeSlab { z_index: z },
            VolumeSlabTexture(image_handle),
        ));
    }
}

/// Each frame: re-sample the chosen field at every slab's Z layer and rebuild
/// the slab textures. Mirrors the Bevy 0.15 dynamic-texture pattern from
/// `slices.rs` — we create a fresh `Image` and swap the material's handle,
/// because `images.get_mut` and `images.insert` don't reliably re-extract.
pub fn update_volume_textures(
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    config: Res<VolumeConfig>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut query: Query<(
        &VolumeSlab,
        &mut VolumeSlabTexture,
        &MeshMaterial3d<StandardMaterial>,
    )>,
) {
    if !config.enabled {
        return;
    }
    let Some(grid) = grid else { return };
    if diag.fields.is_empty() {
        return;
    }

    let max = if config.auto_range {
        compute_max(&grid, &diag, config.field, config.exclude_pml_from_range)
    } else {
        config.manual_max.max(1e-6)
    };

    let nx = grid.nx as usize;
    let ny = grid.ny as usize;

    for (slab, mut slab_tex, mat_handle) in query.iter_mut() {
        let z = slab.z_index as usize;
        if z >= grid.nz as usize {
            continue;
        }

        let mut data = Vec::with_capacity(nx * ny * 4);
        // Iterate y from top → bottom so the texture's V axis aligns with
        // the world's −Y, matching the slice convention.
        for y_iter in 0..ny {
            let y = ny - 1 - y_iter;
            for x in 0..nx {
                let idx = grid.idx(x as u32, y as u32, z as u32);
                let v = sample_field_at(&grid, &diag, idx, config.field);
                let rgba = transfer_for(
                    config.field,
                    v,
                    max,
                    config.opacity_scale,
                    config.color_map,
                );
                data.push((rgba[0] * 255.0) as u8);
                data.push((rgba[1] * 255.0) as u8);
                data.push((rgba[2] * 255.0) as u8);
                data.push((rgba[3] * 255.0) as u8);
            }
        }

        let new_image = Image::new(
            Extent3d {
                width: grid.nx,
                height: grid.ny,
                depth_or_array_layers: 1,
            },
            TextureDimension::D2,
            data,
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
        );
        let new_handle = images.add(new_image);

        if let Some(material) = materials.get_mut(&mat_handle.0) {
            material.base_color_texture = Some(new_handle.clone());
        }

        // Drop the old image asset to free memory.
        images.remove(&slab_tex.0);
        slab_tex.0 = new_handle;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Bipolar transfer function characterisation -----------------------------

    #[test]
    fn test_bipolar_zero_is_transparent() {
        let rgba = transfer_bipolar(0.0, 1.0, 1.0);
        assert_eq!(rgba[3], 0.0);
    }

    #[test]
    fn test_bipolar_positive_max_is_red() {
        let rgba = transfer_bipolar(1.0, 1.0, 1.0);
        assert!((rgba[0] - 1.0).abs() < 1e-6);
        assert!(rgba[1] < 1e-6, "g should be 0, got {}", rgba[1]);
        assert!(rgba[2] < 1e-6, "b should be 0, got {}", rgba[2]);
        assert!((rgba[3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bipolar_negative_max_is_blue() {
        let rgba = transfer_bipolar(-1.0, 1.0, 1.0);
        assert!(rgba[0] < 1e-6, "r should be 0, got {}", rgba[0]);
        assert!(rgba[1] < 1e-6, "g should be 0, got {}", rgba[1]);
        assert!((rgba[2] - 1.0).abs() < 1e-6);
        assert!((rgba[3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bipolar_clamps_above_max() {
        let r1 = transfer_bipolar(1.0, 1.0, 1.0);
        let r2 = transfer_bipolar(5.0, 1.0, 1.0);
        for i in 0..4 {
            assert!((r1[i] - r2[i]).abs() < 1e-6);
        }
    }

    #[test]
    fn test_bipolar_opacity_scale_lowers_alpha() {
        let full = transfer_bipolar(1.0, 1.0, 1.0);
        let half = transfer_bipolar(1.0, 1.0, 0.5);
        assert!((full[3] - 1.0).abs() < 1e-6);
        assert!((half[3] - 0.5).abs() < 1e-6);
    }

    // Sequential transfer function ------------------------------------------

    #[test]
    fn test_sequential_zero_is_transparent() {
        let rgba = transfer_sequential(0.0, 1.0, 1.0, ColorMap::Viridis);
        assert_eq!(rgba[3], 0.0);
    }

    #[test]
    fn test_sequential_max_is_opaque() {
        let rgba = transfer_sequential(1.0, 1.0, 1.0, ColorMap::Viridis);
        assert!((rgba[3] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_sequential_negative_treated_as_zero() {
        let rgba = transfer_sequential(-0.5, 1.0, 1.0, ColorMap::Viridis);
        assert!(rgba[3].abs() < 1e-6);
    }

    #[test]
    fn test_sequential_alpha_linear_in_value() {
        let half = transfer_sequential(0.5, 1.0, 1.0, ColorMap::Viridis);
        assert!(
            (half[3] - 0.5).abs() < 1e-6,
            "alpha should be 0.5 at half-max, got {}",
            half[3]
        );
    }

    // Auto-dispatch ---------------------------------------------------------

    #[test]
    fn test_transfer_dispatch_signed_uses_bipolar() {
        let rgba_signed = transfer_for(
            FieldQuantity::SField,
            -1.0,
            1.0,
            1.0,
            ColorMap::Viridis,
        );
        let bipolar = transfer_bipolar(-1.0, 1.0, 1.0);
        for i in 0..4 {
            assert!((rgba_signed[i] - bipolar[i]).abs() < 1e-6);
        }
    }

    #[test]
    fn test_transfer_dispatch_positive_uses_sequential() {
        let rgba_pos = transfer_for(
            FieldQuantity::EMagnitude,
            0.5,
            1.0,
            1.0,
            ColorMap::Viridis,
        );
        let seq = transfer_sequential(0.5, 1.0, 1.0, ColorMap::Viridis);
        for i in 0..4 {
            assert!((rgba_pos[i] - seq[i]).abs() < 1e-6);
        }
    }

    // Slab placement --------------------------------------------------------

    #[test]
    fn test_slab_z_indices_distribute_evenly() {
        let zs = slab_z_indices(32, 4);
        assert_eq!(zs.first(), Some(&0));
        assert_eq!(zs.last(), Some(&31));
        for w in zs.windows(2) {
            assert!(w[1] > w[0], "indices not increasing: {:?}", zs);
        }
    }

    #[test]
    fn test_slab_z_indices_dedups_when_n_too_large() {
        let zs = slab_z_indices(4, 8);
        assert!(zs.len() <= 4, "got {} slabs for nz=4", zs.len());
        for w in zs.windows(2) {
            assert_ne!(w[0], w[1]);
        }
    }

    #[test]
    fn test_slab_world_z_centres() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        let z0 = slab_world_z(&grid, 0);
        let half = 8.0 * 0.01 * 0.5;
        assert!((z0 - (-half + 0.5 * 0.01)).abs() < 1e-6);
        let zlast = slab_world_z(&grid, 7);
        assert!((zlast - (half - 0.5 * 0.01)).abs() < 1e-6);
    }
}
