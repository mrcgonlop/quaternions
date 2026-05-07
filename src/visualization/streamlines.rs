// Streamline computation and rendering — primary 3D vector view.
//
// Integrates field lines from uniformly seeded starting points using
// 4th-order Runge-Kutta on the NORMALISED field direction (so step length is
// arc length, independent of magnitude — the right thing for "show me the
// topology"). Each streamline follows the chosen vector field until it
// leaves the safe interpolation region or the field magnitude drops below
// a relative threshold.
//
// Rendered as connected line segments via Bevy gizmos, coloured by local
// field magnitude using a configurable color map. Optional animated tracer
// dots flow along each streamline to convey direction at a glance.
//
// Phase 6.2: replaces the previous Euler-stepping implementation; together
// with the future demotion of slices to a 2D inset and volume rendering
// for scalars (Phase 6.1), this becomes the primary 3D view of every
// vector field (E, B, Poynting, A).

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
    /// Maximum number of RK4 integration steps per streamline.
    pub max_steps: u32,
    /// Step size as a fraction of grid spacing dx. Each RK4 step advances
    /// the streamline by `step_fraction * dx` in arc length.
    pub step_fraction: f32,
    /// Color map for field magnitude encoding.
    pub color_map: ColorMap,
    /// If true, auto-scale the color range to the maximum sampled field
    /// magnitude. Auto-range uses the simple max here (not percentile) —
    /// streamlines self-truncate when field is small, so a single source
    /// outlier doesn't suppress the rest of the picture the way it does
    /// for arrow glyphs.
    pub auto_range: bool,
    /// Manual maximum for color range (used when auto_range is false).
    pub manual_max: f32,
    /// If true, animate tracer dots flowing along each streamline.
    pub animate_tracers: bool,
    /// How many tracers ride each streamline simultaneously.
    pub tracers_per_line: u32,
    /// Tracer flow speed in "fraction of streamline per second". 0.5 means
    /// a tracer traverses half the streamline each second.
    pub tracer_speed: f32,
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
            animate_tracers: true,
            tracers_per_line: 2,
            tracer_speed: 0.5,
        }
    }
}

/// Convert fractional grid coordinates to world position. Cell (0, 0, 0)'s
/// centre lives at world `(-half_x + dx/2, -half_y + dx/2, -half_z + dx/2)`.
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
/// Returns `[0, 0, 0]` when the sample point is outside the safe
/// interpolation region (≤ 1 cell from any face).
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

/// Sample the normalised field DIRECTION at fractional grid coords.
/// Returns `(direction, magnitude)`. Returns `(zeros, 0.0)` outside the
/// safe interpolation region or when the local field is zero.
fn sample_direction(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    field: GlyphField,
    fx: f32,
    fy: f32,
    fz: f32,
) -> ([f32; 3], f32) {
    let v = trilinear_sample(grid, diag, field, fx, fy, fz);
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag < 1e-30 {
        return ([0.0; 3], 0.0);
    }
    let inv = 1.0 / mag;
    ([v[0] * inv, v[1] * inv, v[2] * inv], mag)
}

/// Single 4th-order Runge-Kutta step of arc length `h` (in grid units —
/// fraction of dx) along the normalised field direction.
///
/// Returns `Some((fx_new, fy_new, fz_new, mag_at_start))` on success, or
/// `None` if any intermediate sample falls outside the safe interpolation
/// region (which means the streamline would step out of the domain).
fn rk4_step(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    field: GlyphField,
    fx: f32,
    fy: f32,
    fz: f32,
    h: f32,
) -> Option<(f32, f32, f32, f32)> {
    let (k1, mag) = sample_direction(grid, diag, field, fx, fy, fz);
    if mag == 0.0 {
        return None;
    }
    let (k2, m2) = sample_direction(
        grid,
        diag,
        field,
        fx + 0.5 * h * k1[0],
        fy + 0.5 * h * k1[1],
        fz + 0.5 * h * k1[2],
    );
    if m2 == 0.0 {
        return None;
    }
    let (k3, m3) = sample_direction(
        grid,
        diag,
        field,
        fx + 0.5 * h * k2[0],
        fy + 0.5 * h * k2[1],
        fz + 0.5 * h * k2[2],
    );
    if m3 == 0.0 {
        return None;
    }
    let (k4, m4) = sample_direction(
        grid,
        diag,
        field,
        fx + h * k3[0],
        fy + h * k3[1],
        fz + h * k3[2],
    );
    if m4 == 0.0 {
        return None;
    }

    let dx = (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
    let dy = (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
    let dz = (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0;

    Some((fx + h * dx, fy + h * dy, fz + h * dz, mag))
}

/// Trace a single streamline from a fractional-grid seed point.
///
/// Returns a polyline as `Vec<(world_position, local_magnitude)>`. The
/// magnitude entry is the field magnitude at the START of the segment
/// terminating at that point — used by the renderer for per-segment
/// magnitude colouring.
///
/// Stops on any of:
///   * `max_steps` exceeded.
///   * Position leaves the safe interpolation region (≤ 1 cell from a face).
///   * Local field magnitude drops below `mag_threshold`.
pub fn compute_streamline(
    grid: &SimulationGrid,
    diag: &DiagnosticsState,
    field: GlyphField,
    seed_fx: f32,
    seed_fy: f32,
    seed_fz: f32,
    max_steps: usize,
    step_fraction: f32,
    mag_threshold: f32,
) -> Vec<(Vec3, f32)> {
    let nx = grid.nx as f32;
    let ny = grid.ny as f32;
    let nz = grid.nz as f32;
    let h = step_fraction.max(0.01);

    let mut points = Vec::with_capacity(max_steps + 1);
    let mut fx = seed_fx;
    let mut fy = seed_fy;
    let mut fz = seed_fz;
    points.push((frac_to_world(fx, fy, fz, grid), 0.0));

    for _ in 0..max_steps {
        let Some((fx_new, fy_new, fz_new, mag)) =
            rk4_step(grid, diag, field, fx, fy, fz, h)
        else {
            break;
        };
        if mag < mag_threshold {
            break;
        }
        // Stay inside the safe interpolation region.
        if fx_new < 1.0
            || fx_new >= nx - 2.0
            || fy_new < 1.0
            || fy_new >= ny - 2.0
            || fz_new < 1.0
            || fz_new >= nz - 2.0
        {
            break;
        }
        fx = fx_new;
        fy = fy_new;
        fz = fz_new;
        points.push((frac_to_world(fx, fy, fz, grid), mag));
    }

    points
}

/// Linearly interpolate position at fractional parameter `t ∈ [0, 1]` along
/// a polyline. Returns `None` for a polyline shorter than two points.
fn polyline_at(polyline: &[(Vec3, f32)], t: f32) -> Option<Vec3> {
    let n = polyline.len();
    if n < 2 {
        return None;
    }
    let segments = (n - 1) as f32;
    let s = t.clamp(0.0, 0.999_999) * segments;
    let seg_idx = s as usize;
    let seg_t = s - seg_idx as f32;
    let p0 = polyline[seg_idx].0;
    let p1 = polyline[seg_idx + 1].0;
    Some(p0.lerp(p1, seg_t))
}

/// Bevy system: trace streamlines from subsampled seed points, draw the
/// polyline, and (optionally) draw animated tracer dots.
pub fn draw_streamlines(
    mut gizmos: Gizmos,
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    config: Res<StreamlineConfig>,
    time: Res<Time>,
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
    let max_steps = config.max_steps as usize;

    // Auto-range over every sampled point — streamlines self-truncate when
    // the field is small, so the auto-range max is normally driven by the
    // peak interesting region. Outliers compress the colour scale but don't
    // hide line geometry the way they hide arrows in glyph view.
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
    let elapsed = time.elapsed_secs();
    let tracer_radius = grid.dx * 0.18;

    for sz in (1..nz.saturating_sub(1)).step_by(stride) {
        for sy in (1..ny.saturating_sub(1)).step_by(stride) {
            for sx in (1..nx.saturating_sub(1)).step_by(stride) {
                let polyline = compute_streamline(
                    &grid,
                    &diag,
                    config.field,
                    sx as f32,
                    sy as f32,
                    sz as f32,
                    max_steps,
                    step_fraction,
                    threshold,
                );

                if polyline.len() < 2 {
                    continue;
                }

                // Polyline body — colour each segment by the magnitude
                // recorded at its start point.
                for i in 1..polyline.len() {
                    let (a, _) = polyline[i - 1];
                    let (b, mag) = polyline[i];
                    let rgba = color_maps::map_value(mag, 0.0, max_mag, config.color_map);
                    gizmos.line(a, b, Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]));
                }

                // Animated tracer dots, evenly spaced along the polyline.
                if config.animate_tracers && config.tracers_per_line > 0 {
                    let n = config.tracers_per_line as f32;
                    let phase = (elapsed * config.tracer_speed).fract();
                    for k in 0..config.tracers_per_line {
                        let t = ((k as f32 / n) + phase).fract();
                        if let Some(pos) = polyline_at(&polyline, t) {
                            // Bright, opaque dot at the tracer position.
                            gizmos.sphere(
                                Isometry3d::from_translation(pos),
                                tracer_radius,
                                Color::WHITE,
                            );
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::state::DerivedFields;

    /// Build a `DiagnosticsState` whose every cell has the given E vector.
    /// Lets tests drive `compute_streamline` against analytic expectations
    /// without running the full simulation.
    fn uniform_e_field_diag(n: usize, e: [f32; 3]) -> DiagnosticsState {
        let derived = DerivedFields {
            e,
            ..Default::default()
        };
        DiagnosticsState {
            fields: vec![derived; n],
            ..Default::default()
        }
    }

    /// Build a `DiagnosticsState` with E = (-y, x, 0) in grid coordinates,
    /// centred on the grid centre. Field lines are circles in the xy plane;
    /// streamlines should close on themselves to RK4 precision.
    fn vortex_e_field_diag(grid: &SimulationGrid) -> DiagnosticsState {
        let nx = grid.nx as usize;
        let ny = grid.ny as usize;
        let nz = grid.nz as usize;
        let cx = grid.nx as f32 * 0.5;
        let cy = grid.ny as f32 * 0.5;
        let mut fields = vec![DerivedFields::default(); nx * ny * nz];
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dx = x as f32 - cx;
                    let dy = y as f32 - cy;
                    let idx = fdtd::idx(x, y, z, nx, ny);
                    fields[idx].e = [-dy, dx, 0.0];
                }
            }
        }
        DiagnosticsState {
            fields,
            ..Default::default()
        }
    }

    /// In a uniform +ẑ field the streamline must be a vertical straight
    /// line, with each step advancing the same arc length (= step_fraction
    /// × dx) since the integrator works on the NORMALISED direction.
    #[test]
    fn test_streamline_uniform_field_is_straight_line() {
        let grid = SimulationGrid::new(16, 16, 16, 0.01);
        let diag = uniform_e_field_diag(16 * 16 * 16, [0.0, 0.0, 1.0]);
        let max_mag = 1.0;
        let threshold = max_mag * 1e-4;

        let polyline = compute_streamline(
            &grid,
            &diag,
            GlyphField::Electric,
            8.0,
            8.0,
            4.0,
            50,
            0.5,
            threshold,
        );

        assert!(polyline.len() > 5, "expected several integration steps");

        // No drift in x or y.
        let x0 = polyline[0].0.x;
        let y0 = polyline[0].0.y;
        for &(p, _) in &polyline {
            assert!((p.x - x0).abs() < 1e-5, "x drifted: {} vs {x0}", p.x);
            assert!((p.y - y0).abs() < 1e-5, "y drifted: {} vs {y0}", p.y);
        }

        // Strictly monotone in z.
        for i in 1..polyline.len() {
            assert!(
                polyline[i].0.z > polyline[i - 1].0.z,
                "z not monotone at step {i}"
            );
        }

        // Step length should equal step_fraction * dx for a unit-length
        // normalised direction.
        let expected_step = 0.5 * grid.dx;
        let dz = polyline[1].0.z - polyline[0].0.z;
        let rel_err = ((dz - expected_step) / expected_step).abs();
        assert!(
            rel_err < 1e-4,
            "step length: expected {expected_step}, got {dz} (rel_err = {rel_err})"
        );
    }

    /// In a vortex field E = (−y, x, 0), streamlines are circles. A circle
    /// of radius r > 0 should close: after enough RK4 steps to traverse
    /// 2π·r in arc length, the streamline returns near the seed.
    #[test]
    fn test_streamline_vortex_closes_on_itself() {
        let grid = SimulationGrid::new(32, 32, 16, 0.01);
        let diag = vortex_e_field_diag(&grid);

        // Seed off-centre so we have a well-defined radius.
        let seed_x: f32 = 20.0;
        let seed_y: f32 = 16.0;
        let seed_z: f32 = 8.0;
        let radius_grid: f32 = (seed_x - 16.0_f32).hypot(seed_y - 16.0_f32); // = 4.0

        // Arc length around the circle = 2π·r grid units. With step_fraction
        // = h, one full revolution takes round(2π·r/h) steps. We pick
        // `max_steps` to integrate as close to a single full revolution as
        // possible — going past it would put the endpoint partway around a
        // second revolution, which is fine for a topology test but bad for
        // a closure test.
        let h: f32 = 0.5;
        let max_steps = (2.0 * std::f32::consts::PI * radius_grid / h).round() as usize;
        // Magnitude on this circle is r = 4 grid units; threshold scales
        // with grid coords here since the vortex field is unitless.
        let threshold = 1.0e-3;

        let polyline = compute_streamline(
            &grid,
            &diag,
            GlyphField::Electric,
            seed_x,
            seed_y,
            seed_z,
            max_steps,
            h,
            threshold,
        );

        assert!(polyline.len() > 30, "expected a long streamline");

        // After max_steps, the endpoint should be near the seed. Exact
        // closure is impossible because round(2π·r/h) is rarely 2π·r/h
        // exactly — the residual angular gap is δ = (2π·r − max_steps·h)/r
        // and the chord across that gap is 2·r·sin(δ/2). Tolerance below
        // accommodates this geometric residual plus RK4 phase drift.
        let seed_world = polyline[0].0;
        let end_world = polyline.last().unwrap().0;
        let close_err = (end_world - seed_world).length();
        let radius_world = radius_grid * grid.dx;
        let arc_residual =
            2.0 * std::f32::consts::PI * radius_grid - max_steps as f32 * h;
        let chord_residual_grid =
            2.0 * radius_grid * (arc_residual / (2.0 * radius_grid)).abs().sin();
        let chord_residual_world = chord_residual_grid * grid.dx;
        // 5% of radius on top of the geometric chord — RK4 phase drift.
        let tol = chord_residual_world + 0.05 * radius_world;
        assert!(
            close_err < tol,
            "vortex streamline failed to close: error = {close_err:.3e}, radius = {radius_world:.3e}, geometric residual = {chord_residual_world:.3e}, tol = {tol:.3e}"
        );

        // Every point should remain near the same radius (no spiral). The
        // vortex centre in WORLD coordinates is the centre-of-grid cell at
        // (16, 16, *), which `frac_to_world` places at +dx/2 offset from the
        // domain origin (cell-centre indexing). Compute it explicitly so the
        // test isn't tied to the grid-size arithmetic.
        let centre_world = frac_to_world(16.0, 16.0, seed_z, &grid);
        for (i, &(p, _)) in polyline.iter().enumerate() {
            let dx = p.x - centre_world.x;
            let dy = p.y - centre_world.y;
            let r = (dx * dx + dy * dy).sqrt();
            let dr = (r - radius_world).abs() / radius_world;
            assert!(
                dr < 0.02,
                "step {i}: radius drifted by {dr:.3e} relative (r = {r:.3e}, expected {radius_world:.3e})"
            );
        }
    }

    /// Streamline starting outside the safe region returns just the seed.
    #[test]
    fn test_streamline_outside_region_returns_only_seed() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        let diag = uniform_e_field_diag(8 * 8 * 8, [0.0, 0.0, 1.0]);
        // Seed exactly on a face — RK4 step will leave the safe region
        // immediately.
        let polyline = compute_streamline(
            &grid,
            &diag,
            GlyphField::Electric,
            0.5,
            0.5,
            0.5,
            10,
            0.5,
            1e-6,
        );
        assert_eq!(polyline.len(), 1, "expected only the seed point");
    }

    /// Polyline parameter interpolation: t = 0 → first point, t = 1 → last.
    #[test]
    fn test_polyline_at_endpoints() {
        let pl = vec![
            (Vec3::new(0.0, 0.0, 0.0), 1.0),
            (Vec3::new(1.0, 0.0, 0.0), 1.0),
            (Vec3::new(2.0, 0.0, 0.0), 1.0),
        ];
        let p0 = polyline_at(&pl, 0.0).unwrap();
        let p1 = polyline_at(&pl, 0.999_999).unwrap();
        assert!((p0 - Vec3::ZERO).length() < 1e-5);
        assert!((p1 - Vec3::new(2.0, 0.0, 0.0)).length() < 1e-3);
        // Midpoint
        let pm = polyline_at(&pl, 0.5).unwrap();
        assert!((pm - Vec3::new(1.0, 0.0, 0.0)).length() < 1e-5);
    }

    #[test]
    fn test_polyline_at_too_short() {
        let pl = vec![(Vec3::ZERO, 0.0)];
        assert!(polyline_at(&pl, 0.5).is_none());
        let empty: Vec<(Vec3, f32)> = vec![];
        assert!(polyline_at(&empty, 0.5).is_none());
    }
}
