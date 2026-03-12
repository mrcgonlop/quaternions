// Slice plane visualization: renders 2D color-mapped cross-sections through the 3D grid.
//
// Supports two simultaneous slice planes for comparing different field quantities
// (e.g., |E| on one slice and S field on another at orthogonal orientations).
//
// A textured quad is positioned in the 3D scene at each slice location.
// Each frame the textures are updated by sampling chosen field quantities
// from DiagnosticsState (or raw Q from the grid) and mapping through color maps.

use bevy::prelude::*;
use bevy::render::render_asset::RenderAssetUsages;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};

use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;
use crate::visualization::color_maps::{self, ColorEncoding, ColorMap, PhasePlane};
use crate::visualization::glyphs::GlyphField;

/// Number of simultaneous slice planes supported.
pub const NUM_SLICES: usize = 2;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Which axis the slice is perpendicular to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SliceAxis {
    X, // slice is a YZ plane at some x
    Y, // slice is an XZ plane at some y
    Z, // slice is an XY plane at some z
}

impl SliceAxis {
    pub const ALL: &'static [SliceAxis] = &[SliceAxis::X, SliceAxis::Y, SliceAxis::Z];

    pub fn name(&self) -> &'static str {
        match self {
            SliceAxis::X => "X (YZ plane)",
            SliceAxis::Y => "Y (XZ plane)",
            SliceAxis::Z => "Z (XY plane)",
        }
    }

    /// Maximum grid index along this axis.
    pub fn max_index(&self, grid: &SimulationGrid) -> u32 {
        match self {
            SliceAxis::X => grid.nx.saturating_sub(1),
            SliceAxis::Y => grid.ny.saturating_sub(1),
            SliceAxis::Z => grid.nz.saturating_sub(1),
        }
    }
}

/// Which scalar field quantity to visualize on the slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FieldQuantity {
    EMagnitude,
    BMagnitude,
    SField,
    Phi,
    Ax,
    Ay,
    Az,
    EnergyDensity,
    KVacuum,
}

impl FieldQuantity {
    pub const ALL: &'static [FieldQuantity] = &[
        FieldQuantity::EMagnitude,
        FieldQuantity::BMagnitude,
        FieldQuantity::SField,
        FieldQuantity::Phi,
        FieldQuantity::Ax,
        FieldQuantity::Ay,
        FieldQuantity::Az,
        FieldQuantity::EnergyDensity,
        FieldQuantity::KVacuum,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            FieldQuantity::EMagnitude => "|E|",
            FieldQuantity::BMagnitude => "|B|",
            FieldQuantity::SField => "S field",
            FieldQuantity::Phi => "Phi (Q.w)",
            FieldQuantity::Ax => "Ax (Q.x)",
            FieldQuantity::Ay => "Ay (Q.y)",
            FieldQuantity::Az => "Az (Q.z)",
            FieldQuantity::EnergyDensity => "Energy density",
            FieldQuantity::KVacuum => "K (vacuum)",
        }
    }

    /// Whether this quantity is signed (can be negative).
    pub fn is_signed(&self) -> bool {
        matches!(
            self,
            FieldQuantity::SField
                | FieldQuantity::Phi
                | FieldQuantity::Ax
                | FieldQuantity::Ay
                | FieldQuantity::Az
        )
    }

    /// Suggested default color map for this field quantity.
    pub fn default_color_map(&self) -> ColorMap {
        if self.is_signed() {
            ColorMap::ScalarField
        } else {
            ColorMap::Viridis
        }
    }
}

/// Configuration for a single slice plane.
#[derive(Clone, Debug)]
pub struct SingleSliceConfig {
    pub enabled: bool,
    pub axis: SliceAxis,
    /// Grid-coordinate position along the slice axis (integer index).
    pub position: u32,
    pub field: FieldQuantity,
    pub color_map: ColorMap,
    pub auto_range: bool,
    pub manual_min: f32,
    pub manual_max: f32,
    /// Color encoding mode (Standard uses existing single-field color map).
    pub color_encoding: ColorEncoding,
    /// RGB multi-field: which field maps to each channel.
    pub rgb_r_field: FieldQuantity,
    pub rgb_g_field: FieldQuantity,
    pub rgb_b_field: FieldQuantity,
    /// HSV phase: which vector field and plane.
    pub hsv_vector_field: GlyphField,
    pub hsv_phase_plane: PhasePlane,
}

impl Default for SingleSliceConfig {
    fn default() -> Self {
        Self::default_primary()
    }
}

impl SingleSliceConfig {
    /// Primary slice default: Z axis, |E|, Viridis.
    pub fn default_primary() -> Self {
        Self {
            enabled: true,
            axis: SliceAxis::Z,
            position: 0,
            field: FieldQuantity::EMagnitude,
            color_map: ColorMap::Viridis,
            auto_range: true,
            manual_min: 0.0,
            manual_max: 1.0,
            color_encoding: ColorEncoding::Standard,
            rgb_r_field: FieldQuantity::EMagnitude,
            rgb_g_field: FieldQuantity::BMagnitude,
            rgb_b_field: FieldQuantity::SField,
            hsv_vector_field: GlyphField::Electric,
            hsv_phase_plane: PhasePlane::XY,
        }
    }

    /// Secondary slice default: X axis, S field, ScalarField bipolar map, disabled.
    pub fn default_secondary() -> Self {
        Self {
            enabled: false,
            axis: SliceAxis::X,
            position: 0,
            field: FieldQuantity::SField,
            color_map: ColorMap::ScalarField,
            auto_range: true,
            manual_min: -1.0,
            manual_max: 1.0,
            color_encoding: ColorEncoding::Standard,
            rgb_r_field: FieldQuantity::EMagnitude,
            rgb_g_field: FieldQuantity::BMagnitude,
            rgb_b_field: FieldQuantity::SField,
            hsv_vector_field: GlyphField::Electric,
            hsv_phase_plane: PhasePlane::XY,
        }
    }
}

// Keep backward-compatible type alias for tests and external use.
pub type SliceConfig = SingleSliceConfig;

/// Resource holding configuration for all slice planes.
#[derive(Resource)]
pub struct SliceConfigs {
    pub slices: [SingleSliceConfig; NUM_SLICES],
    /// Track previous state per-slice for respawn logic.
    prev_axis: [Option<SliceAxis>; NUM_SLICES],
    prev_enabled: [bool; NUM_SLICES],
}

impl Default for SliceConfigs {
    fn default() -> Self {
        Self {
            slices: [
                SingleSliceConfig::default_primary(),
                SingleSliceConfig::default_secondary(),
            ],
            prev_axis: [None; NUM_SLICES],
            prev_enabled: [false; NUM_SLICES],
        }
    }
}

/// Marker component on slice quad entities.
#[derive(Component)]
pub struct SlicePlane;

/// Identifies which slice index (0 or 1) an entity belongs to.
#[derive(Component, Clone, Copy, PartialEq, Eq)]
pub struct SliceIndex(pub usize);

/// Component storing the handle to the slice's dynamic texture.
#[derive(Component)]
pub struct SliceTexture(pub Handle<Image>);

/// Per-slice statistics computed each frame for UI display.
#[derive(Clone, Debug, Default)]
pub struct SliceStats {
    /// Number of values sampled on the slice.
    pub sample_count: usize,
    /// Minimum sampled value.
    pub value_min: f32,
    /// Maximum sampled value.
    pub value_max: f32,
    /// Color mapping range used (after auto-range or manual).
    pub range_min: f32,
    pub range_max: f32,
}

/// Resource holding per-slice statistics for UI display.
#[derive(Resource, Default)]
pub struct AllSliceStats {
    pub stats: [SliceStats; NUM_SLICES],
}

// ---------------------------------------------------------------------------
// Field sampling
// ---------------------------------------------------------------------------

/// Extract a single scalar value for the chosen field quantity at flat index `idx`.
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

/// Sample the entire 2D slice and return (width, height, flat value array).
///
/// The returned array is in row-major order: for v in 0..height, for u in 0..width.
/// (u, v) map to grid axes depending on the slice axis.
pub fn sample_slice(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    config: &SingleSliceConfig,
) -> (u32, u32, Vec<f32>) {
    let (width, height) = slice_dimensions(grid, config.axis);
    let pos = config.position.min(config.axis.max_index(grid));

    let mut values = Vec::with_capacity((width * height) as usize);

    // Iterate v in reverse: image row 0 (UV v=0) is at the TOP of the quad
    // (positive local Y = positive world Y/Z), so the highest grid index
    // must come first to avoid a vertical flip.
    for v in (0..height).rev() {
        for u in 0..width {
            let (x, y, z) = match config.axis {
                // u→Z, v→Y (matches transform: local X→world Z, local Y→world Y)
                SliceAxis::X => (pos, v, u),
                SliceAxis::Y => (u, pos, v),
                SliceAxis::Z => (u, v, pos),
            };
            let idx = grid.idx(x, y, z);
            values.push(sample_field_at(grid, diag, idx, config.field));
        }
    }

    (width, height, values)
}

/// Sample a vector field at a flat grid index.
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

/// Return (width, height) of the 2D slice for a given axis.
///
/// Width maps to the quad's local X (texture u), height to local Y (texture v).
/// After the slice transform, these must correspond to the correct world axes:
///   X-slice: quad local X → world Z, local Y → world Y  ⇒  width=nz, height=ny
///   Y-slice: quad local X → world X, local Y → world Z  ⇒  width=nx, height=nz
///   Z-slice: quad local X → world X, local Y → world Y  ⇒  width=nx, height=ny
fn slice_dimensions(grid: &SimulationGrid, axis: SliceAxis) -> (u32, u32) {
    match axis {
        SliceAxis::X => (grid.nz, grid.ny),
        SliceAxis::Y => (grid.nx, grid.nz),
        SliceAxis::Z => (grid.nx, grid.ny),
    }
}

// ---------------------------------------------------------------------------
// Bevy systems
// ---------------------------------------------------------------------------

/// Spawn or despawn slice quad entities when config.enabled or axis changes.
pub fn manage_slice_entity(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut images: ResMut<Assets<Image>>,
    mut configs: ResMut<SliceConfigs>,
    grid: Option<Res<SimulationGrid>>,
    existing: Query<(Entity, &SliceIndex), With<SlicePlane>>,
) {
    let Some(grid) = grid else { return };

    for si in 0..NUM_SLICES {
        let cfg = &configs.slices[si];
        let needs_respawn = cfg.enabled != configs.prev_enabled[si]
            || (cfg.enabled && configs.prev_axis[si] != Some(cfg.axis));

        if !needs_respawn {
            continue;
        }

        // Despawn existing entity for this slice index
        for (entity, idx) in existing.iter() {
            if idx.0 == si {
                commands.entity(entity).despawn_recursive();
            }
        }

        configs.prev_enabled[si] = configs.slices[si].enabled;
        configs.prev_axis[si] = Some(configs.slices[si].axis);

        if !configs.slices[si].enabled {
            continue;
        }

        let axis = configs.slices[si].axis;
        let position = configs.slices[si].position;

        // Create texture sized to the slice
        let (w, h) = slice_dimensions(&grid, axis);
        let pixel_count = (w * h) as usize;
        let data = vec![128u8; pixel_count * 4]; // mid-gray placeholder

        let image = Image::new(
            Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            TextureDimension::D2,
            data,
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
        );
        let image_handle = images.add(image);

        let mesh = meshes.add(Rectangle::new(1.0, 1.0));

        let material = materials.add(StandardMaterial {
            base_color_texture: Some(image_handle.clone()),
            unlit: true,
            alpha_mode: AlphaMode::Blend,
            double_sided: true,
            cull_mode: None,
            ..default()
        });

        commands.spawn((
            Mesh3d(mesh),
            MeshMaterial3d(material),
            compute_slice_transform(&grid, axis, position),
            SlicePlane,
            SliceIndex(si),
            SliceTexture(image_handle),
        ));
    }
}

/// Update slice textures each frame from field data.
///
/// Creates a new Image asset each frame and swaps the material's texture handle.
/// This is necessary because Bevy 0.15 does not reliably re-extract modified
/// Image assets to the render world via `get_mut` or `insert`.
pub fn update_slice_texture(
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    configs: Res<SliceConfigs>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut stats: ResMut<AllSliceStats>,
    mut query: Query<
        (&SliceIndex, &mut SliceTexture, &MeshMaterial3d<StandardMaterial>),
        With<SlicePlane>,
    >,
    mut frame_count: Local<u64>,
) {
    let Some(grid) = grid else { return };
    if diag.fields.is_empty() {
        return;
    }

    *frame_count += 1;
    let do_log = *frame_count % 60 == 1; // log once per ~second

    for (slice_idx, mut tex, mat_handle) in query.iter_mut() {
        let si = slice_idx.0;
        if si >= NUM_SLICES {
            continue;
        }
        let cfg = &configs.slices[si];
        if !cfg.enabled {
            continue;
        }

        let (width, height) = slice_dimensions(&grid, cfg.axis);
        let pos = cfg.position.min(cfg.axis.max_index(&grid));
        let pixel_count = (width * height) as usize;

        // Build RGBA pixel data based on encoding mode
        let (rgba, sample_count, data_min, data_max, min_val, max_val) =
            match cfg.color_encoding {
                ColorEncoding::Standard | ColorEncoding::SizeColor => {
                    // Standard single-field encoding (SizeColor not applicable to slices, treat as Standard)
                    let (_w, _h, values) = sample_slice(&grid, &diag, cfg);
                    let dmin = values.iter().copied().fold(f32::INFINITY, f32::min);
                    let dmax = values.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                    let (mn, mx) = if cfg.auto_range {
                        auto_range(&values, cfg.field.is_signed())
                    } else {
                        (cfg.manual_min, cfg.manual_max)
                    };
                    let mut pix = Vec::with_capacity(values.len() * 4);
                    for &v in &values {
                        let color = color_maps::map_value(v, mn, mx, cfg.color_map);
                        pix.push((color[0] * 255.0) as u8);
                        pix.push((color[1] * 255.0) as u8);
                        pix.push((color[2] * 255.0) as u8);
                        pix.push((color[3] * 255.0) as u8);
                    }
                    (pix, values.len(), dmin, dmax, mn, mx)
                }
                ColorEncoding::RgbMultiField => {
                    // Sample 3 independent fields and map to RGB
                    let mut r_vals = Vec::with_capacity(pixel_count);
                    let mut g_vals = Vec::with_capacity(pixel_count);
                    let mut b_vals = Vec::with_capacity(pixel_count);
                    for v in (0..height).rev() {
                        for u in 0..width {
                            let (x, y, z) = match cfg.axis {
                                SliceAxis::X => (pos, v, u),
                                SliceAxis::Y => (u, pos, v),
                                SliceAxis::Z => (u, v, pos),
                            };
                            let idx = grid.idx(x, y, z);
                            r_vals.push(sample_field_at(&grid, &diag, idx, cfg.rgb_r_field));
                            g_vals.push(sample_field_at(&grid, &diag, idx, cfg.rgb_g_field));
                            b_vals.push(sample_field_at(&grid, &diag, idx, cfg.rgb_b_field));
                        }
                    }
                    let r_abs_max = r_vals.iter().copied().fold(0.0f32, |a, b| a.max(b.abs())).max(1e-6);
                    let g_abs_max = g_vals.iter().copied().fold(0.0f32, |a, b| a.max(b.abs())).max(1e-6);
                    let b_abs_max = b_vals.iter().copied().fold(0.0f32, |a, b| a.max(b.abs())).max(1e-6);
                    // Use symmetric range [-max, +max] for signed fields so
                    // negative values map to 0.0 and positive to 1.0 (zero = 0.5).
                    // For unsigned fields, min=0 is fine since values are non-negative.
                    let r_signed = cfg.rgb_r_field.is_signed();
                    let g_signed = cfg.rgb_g_field.is_signed();
                    let b_signed = cfg.rgb_b_field.is_signed();
                    let r_min = if r_signed { -r_abs_max } else { 0.0 };
                    let g_min = if g_signed { -g_abs_max } else { 0.0 };
                    let b_min = if b_signed { -b_abs_max } else { 0.0 };
                    let mut pix = Vec::with_capacity(pixel_count * 4);
                    for i in 0..pixel_count {
                        let color = color_maps::encode_rgb_multi(
                            r_vals[i], r_min, r_abs_max,
                            g_vals[i], g_min, g_abs_max,
                            b_vals[i], b_min, b_abs_max,
                        );
                        pix.push((color[0] * 255.0) as u8);
                        pix.push((color[1] * 255.0) as u8);
                        pix.push((color[2] * 255.0) as u8);
                        pix.push((color[3] * 255.0) as u8);
                    }
                    (pix, pixel_count, 0.0, r_abs_max.max(g_abs_max).max(b_abs_max), 0.0, 1.0)
                }
                ColorEncoding::HsvPhase => {
                    // Sample vector field, encode as HSV
                    let mut max_mag = 0.0f32;
                    let mut vectors = Vec::with_capacity(pixel_count);
                    for v in (0..height).rev() {
                        for u in 0..width {
                            let (x, y, z) = match cfg.axis {
                                SliceAxis::X => (pos, v, u),
                                SliceAxis::Y => (u, pos, v),
                                SliceAxis::Z => (u, v, pos),
                            };
                            let idx = grid.idx(x, y, z);
                            let vec = sample_vector_at(&grid, &diag, idx, cfg.hsv_vector_field);
                            let mag = (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]).sqrt();
                            max_mag = max_mag.max(mag);
                            vectors.push(vec);
                        }
                    }
                    max_mag = max_mag.max(1e-6);
                    let mut pix = Vec::with_capacity(pixel_count * 4);
                    for vec in &vectors {
                        let color = color_maps::encode_hsv_phase(*vec, cfg.hsv_phase_plane, max_mag);
                        pix.push((color[0] * 255.0) as u8);
                        pix.push((color[1] * 255.0) as u8);
                        pix.push((color[2] * 255.0) as u8);
                        pix.push((color[3] * 255.0) as u8);
                    }
                    (pix, pixel_count, 0.0, max_mag, 0.0, max_mag)
                }
            };

        // Update stats for UI display
        stats.stats[si] = SliceStats {
            sample_count,
            value_min: data_min,
            value_max: data_max,
            range_min: min_val,
            range_max: max_val,
        };

        if do_log {
            info!(
                "Slice {si} [{:?} pos={}]: encoding={:?}, field={:?}, samples={}, value=[{:.3e}, {:.3e}], range=[{:.3e}, {:.3e}]",
                cfg.axis, cfg.position, cfg.color_encoding, cfg.field,
                sample_count, data_min, data_max, min_val, max_val
            );
        }

        // Create a genuinely new Image asset with a new handle.
        // Bevy 0.15 does not reliably re-extract modified images to the GPU,
        // so we must create a fresh asset and swap the material's texture reference.
        let new_image = Image::new(
            Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            TextureDimension::D2,
            rgba,
            TextureFormat::Rgba8UnormSrgb,
            RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
        );
        let new_handle = images.add(new_image);

        // Swap the material's texture to the new image
        if let Some(material) = materials.get_mut(&mat_handle.0) {
            material.base_color_texture = Some(new_handle.clone());
        }

        // Remove the old image asset to prevent memory leak
        images.remove(&tex.0);

        // Update the component to track the new handle
        tex.0 = new_handle;
    }
}

/// Update slice quads' positions and orientations when position sliders move.
pub fn update_slice_transform(
    grid: Option<Res<SimulationGrid>>,
    configs: Res<SliceConfigs>,
    mut query: Query<(&SliceIndex, &mut Transform), With<SlicePlane>>,
) {
    let Some(grid) = grid else { return };

    for (slice_idx, mut t) in query.iter_mut() {
        let si = slice_idx.0;
        if si >= NUM_SLICES {
            continue;
        }
        let cfg = &configs.slices[si];
        if !cfg.enabled {
            continue;
        }
        *t = compute_slice_transform(&grid, cfg.axis, cfg.position);
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Compute the 3D transform for the slice quad given axis and grid position.
fn compute_slice_transform(
    grid: &SimulationGrid,
    axis: SliceAxis,
    position: u32,
) -> Transform {
    let domain_x = grid.nx as f32 * grid.dx;
    let domain_y = grid.ny as f32 * grid.dx;
    let domain_z = grid.nz as f32 * grid.dx;

    // World-space offset: grid centered at origin
    // Position in world = (index + 0.5) * dx - domain/2
    let world_pos = |idx: u32, dim_size: u32| -> f32 {
        (idx as f32 + 0.5) * grid.dx - (dim_size as f32 * grid.dx * 0.5)
    };

    match axis {
        SliceAxis::X => {
            let x = world_pos(position, grid.nx);
            Transform {
                translation: Vec3::new(x, 0.0, 0.0),
                rotation: Quat::from_rotation_y(std::f32::consts::FRAC_PI_2),
                scale: Vec3::new(domain_z, domain_y, 1.0),
            }
        }
        SliceAxis::Y => {
            let y = world_pos(position, grid.ny);
            Transform {
                translation: Vec3::new(0.0, y, 0.0),
                rotation: Quat::from_rotation_x(-std::f32::consts::FRAC_PI_2),
                scale: Vec3::new(domain_x, domain_z, 1.0),
            }
        }
        SliceAxis::Z => {
            let z = world_pos(position, grid.nz);
            Transform {
                translation: Vec3::new(0.0, 0.0, z),
                rotation: Quat::IDENTITY,
                scale: Vec3::new(domain_x, domain_y, 1.0),
            }
        }
    }
}

/// Compute auto-range (min, max) for color mapping.
/// For signed fields, uses symmetric range around zero.
/// For unsigned fields, uses [0, max] (or small epsilon if all zero).
pub fn auto_range(values: &[f32], signed: bool) -> (f32, f32) {
    if values.is_empty() {
        return (0.0, 1.0);
    }

    // Minimum range must be > f32::EPSILON so that map_value doesn't
    // fall back to t=0.5 (which produces misleading mid-colormap color).
    let min_range = 1e-6_f32;

    if signed {
        let abs_max = values
            .iter()
            .map(|v| v.abs())
            .fold(0.0f32, f32::max)
            .max(min_range);
        (-abs_max, abs_max)
    } else {
        let max = values.iter().copied().fold(0.0f32, f32::max).max(min_range);
        (0.0, max)
    }
}

// ---------------------------------------------------------------------------
// Initialization helper
// ---------------------------------------------------------------------------

/// Set the initial slice positions to the grid center once the grid is available.
pub fn initialize_slice_position(
    grid: Option<Res<SimulationGrid>>,
    mut configs: ResMut<SliceConfigs>,
    mut done: Local<bool>,
) {
    if *done {
        return;
    }
    let Some(grid) = grid else { return };
    for cfg in configs.slices.iter_mut() {
        cfg.position = cfg.axis.max_index(&grid) / 2;
    }
    *done = true;
}
