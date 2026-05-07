// UiPlugin: egui-based control panels for the simulator.

use bevy::prelude::*;
use bevy::diagnostic::{DiagnosticsStore, FrameTimeDiagnosticsPlugin};
use bevy_egui::{egui, EguiContexts, EguiPlugin};

use crate::scenarios::bifilar_coil;
use crate::scenarios::dipole_radiation::{self, Scenario};
use crate::scenarios::graneau_wire;
use crate::scenarios::vacuum_k;
use crate::simulation::diagnostics::{
    probe_fft, DiagnosticsState, Probe, ProbeField, ProbeSet,
};
use crate::simulation::grid::SimulationGrid;
use crate::simulation::particles::ParticleSystem;
use crate::simulation::plugin::{SimulationConfig, VacuumConfig};
use crate::simulation::sources::{Source, SourceConfig, SourceType};
use crate::simulation::weber::ForceMode;
use crate::visualization::color_maps::{ColorEncoding, ColorMap, PhasePlane};
use crate::visualization::glyphs::{GlyphConfig, GlyphField};
use crate::visualization::slices::{
    self, AllSliceStats, FieldQuantity, SliceAxis, SliceConfigs, NUM_SLICES,
};
use crate::visualization::streamlines::StreamlineConfig;
use crate::visualization::volume::VolumeConfig;

/// Tracks which scenario is currently selected in the UI.
#[derive(Resource, Default)]
struct SelectedScenario {
    current: Option<Scenario>,
}

/// Stable egui texture handles for the slice 2D inset images.
///
/// Egui requires a `TextureHandle` whose lifetime spans frames; recreating
/// it via `ctx.load_texture` every frame allocates a new GPU texture each
/// time. We keep one handle per slice and call `set` on it instead.
#[derive(Resource, Default)]
struct SliceUiTextures {
    handles: [Option<egui::TextureHandle>; NUM_SLICES],
}

/// Plugin that adds bevy_egui and the simulation control panels.
pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin)
            .add_plugins(FrameTimeDiagnosticsPlugin::default())
            .init_resource::<SelectedScenario>()
            .init_resource::<SliceUiTextures>()
            .add_systems(
                Update,
                (
                    ui_side_panel,
                    ui_source_panel,
                    ui_diagnostics_panel,
                    ui_slice_panel,
                    ui_glyph_panel,
                    ui_streamline_panel,
                    ui_volume_panel,
                    ui_probe_panel,
                ),
            );
    }
}

/// Main egui side panel: simulation controls, grid info, FPS counter.
fn ui_side_panel(
    mut contexts: EguiContexts,
    mut config: ResMut<SimulationConfig>,
    mut grid: Option<ResMut<SimulationGrid>>,
    diagnostics: Res<DiagnosticsStore>,
    mut sources: ResMut<SourceConfig>,
    mut selected_scenario: ResMut<SelectedScenario>,
    mut pml: Option<ResMut<crate::simulation::boundaries::PmlState>>,
    mut vacuum_config: ResMut<VacuumConfig>,
    mut probes: ResMut<ProbeSet>,
    mut force_mode: ResMut<ForceMode>,
    mut particle_system: ResMut<ParticleSystem>,
) {
    let ctx = contexts.ctx_mut();

    egui::SidePanel::left("sim_control_panel")
        .resizable(true)
        .default_width(240.0)
        .show(ctx, |ui| {
            ui.heading("Simulation");
            ui.separator();

            // --- Play / Pause / Step controls ---
            ui.horizontal(|ui| {
                let label = if config.paused { "Play" } else { "Pause" };
                if ui.button(label).clicked() {
                    config.paused = !config.paused;
                }
                if ui.button("Step").clicked() {
                    config.step_requested = true;
                }
                if ui.button("Reset").clicked() {
                    if let Some(ref mut grid) = grid {
                        grid.reset();
                    }
                    // Clear stale CPML auxiliary fields so the new run starts clean.
                    if let Some(ref mut pml_state) = pml {
                        pml_state.reset_psi();
                    }
                    sources.sources.clear();
                    probes.clear();
                    particle_system.clear();
                    selected_scenario.current = None;
                    config.paused = true;
                }
            });

            ui.separator();

            // --- Mode toggle ---
            ui.horizontal(|ui| {
                ui.label("Mode:");
                ui.selectable_value(&mut config.extended_mode, false, "Standard EM");
                ui.selectable_value(&mut config.extended_mode, true, "QVED Extended");
            });

            // --- Particle force law selector ---
            ui.collapsing(
                format!("Charged Particles ({})", particle_system.particles.len()),
                |ui| {
                    ui.label("Force law:");
                    ui.horizontal(|ui| {
                        for &mode in ForceMode::ALL {
                            ui.selectable_value(force_mode.as_mut(), mode, mode.name());
                        }
                    });
                    ui.label(match *force_mode {
                        ForceMode::Lorentz => {
                            "F = q(E + v × B), Boris pusher, fields from grid."
                        }
                        ForceMode::Weber => {
                            "F = (q₁q₂/4πε₀r²)·r̂·[1 − ṙ²/2c² + r·r̈/c²]. \
                             Direct pair sum, no field coupling."
                        }
                        ForceMode::Both => {
                            "Lorentz + Weber applied additively each substep."
                        }
                    });
                },
            );

            ui.separator();

            // --- Steps per frame ---
            ui.horizontal(|ui| {
                ui.label("Steps/frame:");
                ui.add(egui::DragValue::new(&mut config.steps_per_frame).range(1..=1000));
            });

            // --- Timestep factor ---
            ui.horizontal(|ui| {
                ui.label("dt factor:");
                ui.add(egui::Slider::new(&mut config.dt_factor, 0.01..=1.0).logarithmic(true));
            });

            if let Some(ref grid) = grid {
                let effective_dt = grid.dt * config.dt_factor;
                ui.label(format!("Effective dt: {:.3e} s", effective_dt));
                ui.label(format!(
                    "Sim time/frame: {:.3e} s",
                    effective_dt * config.steps_per_frame as f32
                ));
            }

            ui.separator();

            // --- Grid info ---
            ui.heading("Grid Info");
            if let Some(ref grid) = grid {
                egui::Grid::new("grid_info").show(ui, |ui| {
                    ui.label("Dimensions:");
                    ui.label(format!("{}x{}x{}", grid.nx, grid.ny, grid.nz));
                    ui.end_row();

                    ui.label("dx:");
                    ui.label(format!("{:.4} m", grid.dx));
                    ui.end_row();

                    ui.label("dt:");
                    ui.label(format!("{:.4e} s", grid.dt));
                    ui.end_row();

                    ui.label("Time:");
                    ui.label(format!("{:.6e} s", grid.time));
                    ui.end_row();

                    ui.label("Iteration:");
                    ui.label(format!("{}", grid.iteration));
                    ui.end_row();

                    ui.label("Cells:");
                    ui.label(format!("{}", grid.cell_count()));
                    ui.end_row();
                });
            } else {
                ui.label("Grid not initialized");
            }

            ui.separator();

            // --- FPS counter ---
            ui.heading("Performance");
            if let Some(fps) = diagnostics.get(&FrameTimeDiagnosticsPlugin::FPS) {
                if let Some(value) = fps.smoothed() {
                    ui.label(format!("FPS: {:.0}", value));
                }
            }
            if let Some(frame_time) = diagnostics.get(&FrameTimeDiagnosticsPlugin::FRAME_TIME) {
                if let Some(value) = frame_time.smoothed() {
                    ui.label(format!("Frame: {:.2} ms", value));
                }
            }

            ui.separator();

            // --- Scenario selector ---
            ui.heading("Scenarios");
            let current_name = selected_scenario
                .current
                .map_or("(none)", |s| s.name());
            egui::ComboBox::from_id_salt("scenario_selector")
                .selected_text(current_name)
                .show_ui(ui, |ui| {
                    if ui
                        .selectable_label(selected_scenario.current.is_none(), "(none)")
                        .clicked()
                    {
                        selected_scenario.current = None;
                    }
                    for &scenario in Scenario::ALL {
                        if ui
                            .selectable_label(
                                selected_scenario.current == Some(scenario),
                                scenario.name(),
                            )
                            .clicked()
                        {
                            selected_scenario.current = Some(scenario);
                            // Apply scenario configuration; grid.reset() inside preserves
                            // PML flags, but stale CPML psi must be zeroed separately.
                            if let Some(ref mut grid) = grid {
                                // Scenario switch invalidates any previously recorded
                                // probe history — clear by default. BifilarPair
                                // re-installs a default probe set below. Likewise
                                // clear particles; GraneauWire repopulates them.
                                probes.clear();
                                particle_system.clear();
                                match scenario {
                                    Scenario::DipoleRadiation => {
                                        dipole_radiation::apply_dipole_scenario(grid, &mut sources);
                                        if let Some(ref mut pml_state) = pml {
                                            pml_state.reset_psi();
                                        }
                                        config.paused = true;
                                    }
                                    Scenario::VacuumK => {
                                        vacuum_k::apply_vacuum_k_scenario(
                                            grid,
                                            &mut sources,
                                            &mut vacuum_config,
                                        );
                                        if let Some(ref mut pml_state) = pml {
                                            pml_state.reset_psi();
                                        }
                                        config.extended_mode = false;
                                        config.paused = true;
                                    }
                                    Scenario::BifilarCoil => {
                                        bifilar_coil::apply_bifilar_scenario(
                                            grid,
                                            &mut sources,
                                        );
                                        if let Some(ref mut pml_state) = pml {
                                            pml_state.reset_psi();
                                        }
                                        config.extended_mode = true; // best viewed in extended mode
                                        config.paused = true;
                                    }
                                    Scenario::BifilarPair => {
                                        bifilar_coil::apply_bifilar_pair_scenario(
                                            grid,
                                            &mut sources,
                                        );
                                        bifilar_coil::install_bifilar_pair_probes(
                                            grid,
                                            &mut probes,
                                        );
                                        if let Some(ref mut pml_state) = pml {
                                            pml_state.reset_psi();
                                        }
                                        config.extended_mode = true;
                                        config.paused = true;
                                    }
                                    Scenario::GraneauWire => {
                                        // Graneau is a particle-only scenario:
                                        // no field sources, no PML transients,
                                        // pure Weber pair force on the chain.
                                        grid.reset();
                                        sources.sources.clear();
                                        if let Some(ref mut pml_state) = pml {
                                            pml_state.reset_psi();
                                        }
                                        graneau_wire::apply_graneau_wire_default(
                                            &mut particle_system,
                                            &mut force_mode,
                                        );
                                        config.extended_mode = false;
                                        config.paused = true;
                                    }
                                }
                            }
                        }
                    }
                });

            ui.separator();

            // --- Parameter panel placeholder ---
            ui.heading("Parameters");
            ui.label("(parameter sliders coming in future tasks)");

            // PSEUDOCODE: Future parameter sliders will include:
            // - Grid size selector (requires re-creating buffers)
            // - dx slider (changes resolution vs domain size trade-off)
            // - Boundary condition type selector
            // - Vacuum model parameters (coupling constant α)
            // - Visualization mode toggles
        });
}

/// Diagnostics panel: total energy, max |E|, max |B|, max |S|.
fn ui_diagnostics_panel(
    mut contexts: EguiContexts,
    diag: Res<DiagnosticsState>,
    pml: Option<Res<crate::simulation::boundaries::PmlState>>,
) {
    let ctx = contexts.ctx_mut();

    egui::Window::new("Diagnostics")
        .default_pos([260.0, 300.0])
        .default_width(220.0)
        .show(ctx, |ui| {
            egui::Grid::new("diag_grid").show(ui, |ui| {
                ui.label("Total Energy:");
                ui.label(format!("{:.4e} J", diag.total_energy));
                ui.end_row();

                ui.label("Max |E|:");
                ui.label(format!("{:.4e} V/m", diag.max_e));
                ui.end_row();

                ui.label("Max |B|:");
                ui.label(format!("{:.4e} T", diag.max_b));
                ui.end_row();

                ui.label("Max |S|:");
                ui.label(format!("{:.4e}", diag.max_s));
                ui.end_row();

                ui.label("Max K:");
                ui.label(format!("{:.4}", diag.max_k));
                ui.end_row();

                ui.label("Mean K:");
                ui.label(format!("{:.4}", diag.mean_k));
                ui.end_row();

                ui.label("Topo charge:");
                let n_topo = diag.topological_charge;
                let near_integer = (n_topo.round() - n_topo).abs() < 0.1;
                let label = format!("{:.2}", n_topo);
                if near_integer && n_topo.abs() > 0.1 {
                    ui.colored_label(egui::Color32::from_rgb(100, 255, 100), label);
                } else {
                    ui.label(label);
                }
                ui.end_row();
            });

            if let Some(pml) = &pml {
                ui.separator();
                ui.label(format!(
                    "PML: depth={}, {} cells",
                    pml.config.depth, pml.pml_cell_count
                ));
            }
        });
}

/// Source configuration panel: add/remove sources, adjust parameters.
fn ui_source_panel(
    mut contexts: EguiContexts,
    mut sources: ResMut<SourceConfig>,
    grid: Option<Res<SimulationGrid>>,
) {
    let Some(ref grid) = grid else { return };
    let ctx = contexts.ctx_mut();

    egui::Window::new("Sources")
        .default_pos([260.0, 10.0])
        .default_width(260.0)
        .show(ctx, |ui| {
            // Quick-add buttons
            ui.horizontal(|ui| {
                if ui.button("Add Dipole (Z)").clicked() {
                    let center = [grid.nx as f32 / 2.0, grid.ny as f32 / 2.0, grid.nz as f32 / 2.0];
                    sources.sources.push(Source::dipole_z(center, 1.0, 1e9));
                }
                if ui.button("Add Charge").clicked() {
                    let center = [grid.nx as f32 / 2.0, grid.ny as f32 / 2.0, grid.nz as f32 / 2.0];
                    sources.sources.push(Source::point_charge(center, 1e-10));
                }
            });

            ui.separator();

            // Per-source controls
            let mut to_remove = None;
            for (idx, source) in sources.sources.iter_mut().enumerate() {
                ui.push_id(idx, |ui| {
                    ui.horizontal(|ui| {
                        let label = match &source.source_type {
                            SourceType::PointCharge => "Charge",
                            SourceType::OscillatingDipole { .. } => "Dipole",
                            SourceType::CurrentPulse { .. } => "Pulse",
                            SourceType::ScalarDrive => "Scalar",
                        };
                        ui.checkbox(&mut source.active, label);
                        if ui.small_button("X").clicked() {
                            to_remove = Some(idx);
                        }
                    });

                    egui::Grid::new(format!("source_{idx}")).show(ui, |ui| {
                        ui.label("X:");
                        ui.add(egui::DragValue::new(&mut source.position[0]).range(0.0..=(grid.nx as f32 - 1.0)).speed(0.1));
                        ui.label("Y:");
                        ui.add(egui::DragValue::new(&mut source.position[1]).range(0.0..=(grid.ny as f32 - 1.0)).speed(0.1));
                        ui.label("Z:");
                        ui.add(egui::DragValue::new(&mut source.position[2]).range(0.0..=(grid.nz as f32 - 1.0)).speed(0.1));
                        ui.end_row();

                        ui.label("Amplitude:");
                        ui.add(egui::DragValue::new(&mut source.amplitude).speed(0.01));
                        ui.end_row();

                        if matches!(source.source_type, SourceType::OscillatingDipole { .. } | SourceType::CurrentPulse { .. }) {
                            ui.label("Freq (Hz):");
                            ui.add(egui::DragValue::new(&mut source.frequency).speed(1e6));
                            ui.end_row();
                        }
                    });

                    ui.separator();
                });
            }

            if let Some(idx) = to_remove {
                sources.sources.remove(idx);
            }
        });
}

/// Slice plane configuration panel — supports two simultaneous slices.
fn ui_slice_panel(
    mut contexts: EguiContexts,
    mut configs: ResMut<SliceConfigs>,
    grid: Option<Res<SimulationGrid>>,
    diag: Res<DiagnosticsState>,
    stats: Res<AllSliceStats>,
    mut ui_textures: ResMut<SliceUiTextures>,
) {
    let ctx = contexts.ctx_mut();

    egui::Window::new("Slice Planes")
        .default_pos([800.0, 10.0])
        .default_width(240.0)
        .show(ctx, |ui| {
            let names = ["Slice 1 (Primary)", "Slice 2 (S-field)"];
            for si in 0..NUM_SLICES {
                let id_prefix = format!("slice_{si}");

                egui::CollapsingHeader::new(names[si])
                    .default_open(si == 0)
                    .show(ui, |ui| {
                        let cfg = &mut configs.slices[si];

                        ui.checkbox(&mut cfg.enabled, "Enabled");

                        if !cfg.enabled {
                            return;
                        }

                        ui.separator();

                        // Axis selector
                        ui.label("Slice axis:");
                        ui.horizontal(|ui| {
                            for &axis in SliceAxis::ALL {
                                ui.selectable_value(&mut cfg.axis, axis, axis.name());
                            }
                        });

                        // Position slider
                        if let Some(ref grid) = grid {
                            let max = cfg.axis.max_index(grid);
                            ui.horizontal(|ui| {
                                ui.label("Position:");
                                let mut pos = cfg.position as f32;
                                ui.add(egui::Slider::new(&mut pos, 0.0..=max as f32).integer());
                                cfg.position = pos as u32;
                            });
                        }

                        ui.separator();

                        // Field quantity selector
                        let prev_field = cfg.field;
                        ui.label("Field:");
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_field"))
                            .selected_text(cfg.field.name())
                            .show_ui(ui, |ui| {
                                for &fq in FieldQuantity::ALL {
                                    ui.selectable_value(&mut cfg.field, fq, fq.name());
                                }
                            });

                        // Auto-switch color map when field changes to/from signed
                        if cfg.field != prev_field {
                            cfg.color_map = cfg.field.default_color_map();
                        }

                        // Color encoding mode
                        ui.label("Encoding:");
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_enc"))
                            .selected_text(cfg.color_encoding.name())
                            .show_ui(ui, |ui| {
                                for &enc in ColorEncoding::SLICE_MODES {
                                    ui.selectable_value(&mut cfg.color_encoding, enc, enc.name());
                                }
                            });

                        match cfg.color_encoding {
                            ColorEncoding::Standard | ColorEncoding::SizeColor => {
                                // Color map selector (standard mode)
                                ui.label("Color map:");
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_cmap"))
                                    .selected_text(cfg.color_map.name())
                                    .show_ui(ui, |ui| {
                                        for &cm in ColorMap::ALL {
                                            ui.selectable_value(&mut cfg.color_map, cm, cm.name());
                                        }
                                    });
                            }
                            ColorEncoding::RgbMultiField => {
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_rgb_r"))
                                    .selected_text(format!("R: {}", cfg.rgb_r_field.name()))
                                    .show_ui(ui, |ui| {
                                        for &fq in FieldQuantity::ALL {
                                            ui.selectable_value(&mut cfg.rgb_r_field, fq, fq.name());
                                        }
                                    });
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_rgb_g"))
                                    .selected_text(format!("G: {}", cfg.rgb_g_field.name()))
                                    .show_ui(ui, |ui| {
                                        for &fq in FieldQuantity::ALL {
                                            ui.selectable_value(&mut cfg.rgb_g_field, fq, fq.name());
                                        }
                                    });
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_rgb_b"))
                                    .selected_text(format!("B: {}", cfg.rgb_b_field.name()))
                                    .show_ui(ui, |ui| {
                                        for &fq in FieldQuantity::ALL {
                                            ui.selectable_value(&mut cfg.rgb_b_field, fq, fq.name());
                                        }
                                    });
                            }
                            ColorEncoding::HsvPhase => {
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_hsv_field"))
                                    .selected_text(format!("Vector: {:?}", cfg.hsv_vector_field))
                                    .show_ui(ui, |ui| {
                                        for &gf in GlyphField::ALL {
                                            ui.selectable_value(&mut cfg.hsv_vector_field, gf, gf.name());
                                        }
                                    });
                                egui::ComboBox::from_id_salt(format!("{id_prefix}_hsv_plane"))
                                    .selected_text(format!("Plane: {}", cfg.hsv_phase_plane.name()))
                                    .show_ui(ui, |ui| {
                                        for &pp in PhasePlane::ALL {
                                            ui.selectable_value(&mut cfg.hsv_phase_plane, pp, pp.name());
                                        }
                                    });
                            }
                        }

                        ui.separator();

                        // Range controls
                        ui.checkbox(&mut cfg.auto_range, "Auto range");
                        if !cfg.auto_range {
                            ui.horizontal(|ui| {
                                ui.label("Min:");
                                ui.add(egui::DragValue::new(&mut cfg.manual_min).speed(0.01));
                            });
                            ui.horizontal(|ui| {
                                ui.label("Max:");
                                ui.add(egui::DragValue::new(&mut cfg.manual_max).speed(0.01));
                            });
                        }

                        // Show-in-3D toggle (default off — slice is now a
                        // 2D inset only). The 3D quad competes visually with
                        // the volume render and is rarely useful in 3D.
                        ui.separator();
                        ui.checkbox(&mut cfg.show_in_3d, "Show 3D quad in scene");

                        // 2D inset image — the live slice rendered into the
                        // panel. Computed every frame from the same helper
                        // the 3D-quad path uses, so they stay in sync.
                        ui.separator();
                        if let Some(ref grid_ref) = grid {
                            if !diag.fields.is_empty() {
                                let cfg_snapshot = cfg.clone();
                                let (rgba, w, h, _) = slices::compute_slice_pixels(
                                    grid_ref,
                                    &diag,
                                    &cfg_snapshot,
                                );
                                let color_image = egui::ColorImage::from_rgba_unmultiplied(
                                    [w as usize, h as usize],
                                    &rgba,
                                );
                                match &mut ui_textures.handles[si] {
                                    Some(handle) => {
                                        handle.set(color_image, egui::TextureOptions::LINEAR);
                                    }
                                    slot @ None => {
                                        *slot = Some(ctx.load_texture(
                                            format!("slice_inset_{si}"),
                                            color_image,
                                            egui::TextureOptions::LINEAR,
                                        ));
                                    }
                                }
                                if let Some(handle) = &ui_textures.handles[si] {
                                    let inset_size = egui::vec2(180.0, 180.0);
                                    ui.image((handle.id(), inset_size));
                                }
                            }
                        }

                        // Stats display
                        ui.separator();
                        let st = &stats.stats[si];
                        if st.sample_count > 0 {
                            ui.label(format!("Samples: {}", st.sample_count));
                            ui.label(format!("Value: [{:.4e}, {:.4e}]", st.value_min, st.value_max));
                            ui.label(format!("Range: [{:.4e}, {:.4e}]", st.range_min, st.range_max));
                        } else {
                            ui.label("No data");
                        }
                    });

                if si < NUM_SLICES - 1 {
                    ui.separator();
                }
            }
        });
}

/// Glyph configuration panel: vector field arrows in 3D.
fn ui_glyph_panel(mut contexts: EguiContexts, mut config: ResMut<GlyphConfig>) {
    let ctx = contexts.ctx_mut();

    egui::Window::new("Vector Glyphs")
        .default_pos([800.0, 300.0])
        .default_width(240.0)
        .show(ctx, |ui| {
            ui.checkbox(&mut config.enabled, "Enabled");

            if !config.enabled {
                return;
            }

            // Vector field selector
            let current_name = config.field.name();
            egui::ComboBox::from_label("Field")
                .selected_text(current_name)
                .show_ui(ui, |ui| {
                    for &gf in GlyphField::ALL {
                        ui.selectable_value(&mut config.field, gf, gf.name());
                    }
                });

            // Stride
            ui.add(
                egui::Slider::new(&mut config.stride, 1..=8).text("Stride"),
            );

            // Scale (cell widths at the auto-range reference magnitude).
            // 0.5 means a "typical" arrow is half a cell long; outliers grow
            // proportionally but stay bounded relative to the grid.
            ui.add(
                egui::Slider::new(&mut config.scale, 0.05..=3.0)
                    .text("Scale (cells at p95)"),
            );

            // Auto-range percentile — clips the top (1 - p) fraction of the
            // magnitude distribution so a single source spike doesn't compress
            // every other arrow into invisibility.
            ui.add(
                egui::Slider::new(&mut config.auto_range_percentile, 0.5..=1.0)
                    .text("Auto-range percentile"),
            );

            ui.separator();

            // Color encoding mode
            let enc_name = config.color_encoding.name();
            egui::ComboBox::from_label("Encoding")
                .selected_text(enc_name)
                .show_ui(ui, |ui| {
                    for &enc in ColorEncoding::ALL {
                        ui.selectable_value(&mut config.color_encoding, enc, enc.name());
                    }
                });

            match config.color_encoding {
                ColorEncoding::Standard => {
                    // Color map
                    let map_name = config.color_map.name();
                    egui::ComboBox::from_label("Color Map")
                        .selected_text(map_name)
                        .show_ui(ui, |ui| {
                            for &cm in ColorMap::ALL {
                                ui.selectable_value(&mut config.color_map, cm, cm.name());
                            }
                        });
                }
                ColorEncoding::RgbMultiField => {
                    let r_name = config.rgb_config.r_field.name();
                    egui::ComboBox::from_label("R channel")
                        .selected_text(r_name)
                        .show_ui(ui, |ui| {
                            for &fq in FieldQuantity::ALL {
                                ui.selectable_value(&mut config.rgb_config.r_field, fq, fq.name());
                            }
                        });
                    let g_name = config.rgb_config.g_field.name();
                    egui::ComboBox::from_label("G channel")
                        .selected_text(g_name)
                        .show_ui(ui, |ui| {
                            for &fq in FieldQuantity::ALL {
                                ui.selectable_value(&mut config.rgb_config.g_field, fq, fq.name());
                            }
                        });
                    let b_name = config.rgb_config.b_field.name();
                    egui::ComboBox::from_label("B channel")
                        .selected_text(b_name)
                        .show_ui(ui, |ui| {
                            for &fq in FieldQuantity::ALL {
                                ui.selectable_value(&mut config.rgb_config.b_field, fq, fq.name());
                            }
                        });
                }
                ColorEncoding::HsvPhase => {
                    let plane_name = config.hsv_config.plane.name();
                    egui::ComboBox::from_label("Phase Plane")
                        .selected_text(plane_name)
                        .show_ui(ui, |ui| {
                            for &pp in PhasePlane::ALL {
                                ui.selectable_value(&mut config.hsv_config.plane, pp, pp.name());
                            }
                        });
                }
                ColorEncoding::SizeColor => {
                    let size_name = config.size_color_config.size_field.name();
                    egui::ComboBox::from_label("Size field")
                        .selected_text(size_name)
                        .show_ui(ui, |ui| {
                            for &fq in FieldQuantity::ALL {
                                ui.selectable_value(
                                    &mut config.size_color_config.size_field,
                                    fq,
                                    fq.name(),
                                );
                            }
                        });

                    let color_name = config.size_color_config.color_field.name();
                    egui::ComboBox::from_label("Color field")
                        .selected_text(color_name)
                        .show_ui(ui, |ui| {
                            for &fq in FieldQuantity::ALL {
                                ui.selectable_value(
                                    &mut config.size_color_config.color_field,
                                    fq,
                                    fq.name(),
                                );
                            }
                        });

                    let cm_name = config.size_color_config.color_map.name();
                    egui::ComboBox::from_label("Color Map")
                        .selected_text(cm_name)
                        .show_ui(ui, |ui| {
                            for &cm in ColorMap::ALL {
                                ui.selectable_value(
                                    &mut config.size_color_config.color_map,
                                    cm,
                                    cm.name(),
                                );
                            }
                        });
                }
            }

            ui.separator();
            ui.checkbox(&mut config.auto_range, "Auto range");
            if !config.auto_range {
                ui.add(egui::DragValue::new(&mut config.manual_min).prefix("Min: ").speed(0.01));
                ui.add(egui::DragValue::new(&mut config.manual_max).prefix("Max: ").speed(0.01));
            }
        });
}

/// Streamline panel: RK4-traced field lines with optional animated tracers.
///
/// Phase 6.2 — primary 3D vector view, replaces arrow glyphs as the default
/// way to read a vector field. Each seed point traces forward along the
/// field direction using 4th-order Runge-Kutta; the polyline is coloured by
/// local magnitude. Tracer dots flow along the line over time so direction
/// is readable at a glance.
fn ui_streamline_panel(
    mut contexts: EguiContexts,
    mut config: ResMut<StreamlineConfig>,
) {
    let ctx = contexts.ctx_mut();

    egui::Window::new("Streamlines")
        .default_pos([800.0, 540.0])
        .default_width(240.0)
        .show(ctx, |ui| {
            ui.checkbox(&mut config.enabled, "Enabled");
            if !config.enabled {
                return;
            }

            // Vector field selector (same set as glyphs).
            let field_name = config.field.name();
            egui::ComboBox::from_label("Field")
                .selected_text(field_name)
                .show_ui(ui, |ui| {
                    for &gf in GlyphField::ALL {
                        ui.selectable_value(&mut config.field, gf, gf.name());
                    }
                });

            ui.add(
                egui::Slider::new(&mut config.seed_stride, 2..=12).text("Seed stride"),
            );
            ui.add(
                egui::Slider::new(&mut config.max_steps, 16..=512).text("Max steps"),
            );
            ui.add(
                egui::Slider::new(&mut config.step_fraction, 0.1..=1.0)
                    .text("Step (×dx)"),
            );

            ui.separator();

            // Colour map.
            let map_name = config.color_map.name();
            egui::ComboBox::from_label("Color map")
                .selected_text(map_name)
                .show_ui(ui, |ui| {
                    for &cm in ColorMap::ALL {
                        ui.selectable_value(&mut config.color_map, cm, cm.name());
                    }
                });
            ui.checkbox(&mut config.auto_range, "Auto range");
            if !config.auto_range {
                ui.add(
                    egui::DragValue::new(&mut config.manual_max)
                        .prefix("Max: ")
                        .speed(0.01),
                );
            }

            ui.separator();

            // Animated tracer controls.
            ui.checkbox(&mut config.animate_tracers, "Animated tracers");
            if config.animate_tracers {
                ui.add(
                    egui::Slider::new(&mut config.tracers_per_line, 0..=6)
                        .text("Tracers / line"),
                );
                ui.add(
                    egui::Slider::new(&mut config.tracer_speed, 0.05..=2.0)
                        .text("Tracer speed"),
                );
            }
        });
}

/// Volume rendering panel — primary 3D scalar view (Phase 6.1).
///
/// Stack of axis-aligned semi-transparent slabs. Signed scalars (S, φ, A_*)
/// auto-render with a bipolar blue/white/red transfer function; positive
/// scalars (|E|, |B|, K, energy density) use a configurable sequential
/// palette. Best viewed roughly along ±Z; edge-on views collapse the slabs.
fn ui_volume_panel(
    mut contexts: EguiContexts,
    mut config: ResMut<VolumeConfig>,
) {
    let ctx = contexts.ctx_mut();

    egui::Window::new("Volume")
        .default_pos([800.0, 60.0])
        .default_width(240.0)
        .show(ctx, |ui| {
            ui.checkbox(&mut config.enabled, "Enabled");
            if !config.enabled {
                return;
            }

            // Scalar field selector (full FieldQuantity list — sign auto-
            // dispatches to the correct transfer function).
            let field_name = config.field.name();
            egui::ComboBox::from_label("Field")
                .selected_text(field_name)
                .show_ui(ui, |ui| {
                    for &fq in FieldQuantity::ALL {
                        ui.selectable_value(&mut config.field, fq, fq.name());
                    }
                });
            ui.label(if config.field.is_signed() {
                "Transfer: bipolar (blue / clear / red)"
            } else {
                "Transfer: sequential (palette + opacity ramp)"
            });

            ui.add(
                egui::Slider::new(&mut config.num_slabs, 4..=64).text("Slabs"),
            );
            ui.add(
                egui::Slider::new(&mut config.opacity_scale, 0.05..=2.0)
                    .text("Opacity scale"),
            );

            ui.separator();

            // Sequential-only controls (the bipolar transfer function uses
            // a hand-coded palette that ignores ColorMap).
            if !config.field.is_signed() {
                let map_name = config.color_map.name();
                egui::ComboBox::from_label("Color map")
                    .selected_text(map_name)
                    .show_ui(ui, |ui| {
                        for &cm in ColorMap::ALL {
                            ui.selectable_value(&mut config.color_map, cm, cm.name());
                        }
                    });
            }

            ui.checkbox(&mut config.auto_range, "Auto range");
            if !config.auto_range {
                ui.add(
                    egui::DragValue::new(&mut config.manual_max)
                        .prefix("Max: ")
                        .speed(0.01),
                );
            }
            ui.checkbox(
                &mut config.exclude_pml_from_range,
                "Exclude PML from range",
            );
        });
}

/// Probe panel: add/remove probes, view time series and optional FFT.
///
/// Probes are point samplers that record a selected scalar field at a fixed
/// grid location each simulation step. Used for quantitative analysis like
/// measuring scalar-wave propagation delay in the bifilar pair scenario.
fn ui_probe_panel(
    mut contexts: EguiContexts,
    mut probes: ResMut<ProbeSet>,
    grid: Option<Res<SimulationGrid>>,
) {
    let Some(ref grid) = grid else { return };
    let ctx = contexts.ctx_mut();

    egui::Window::new("Probes")
        .default_pos([260.0, 500.0])
        .default_width(360.0)
        .default_height(380.0)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("Add @ center").clicked() {
                    let cx = grid.nx / 2;
                    let cy = grid.ny / 2;
                    let cz = grid.nz / 2;
                    let label = format!("Probe {}", probes.probes.len());
                    probes.push(Probe::new(label, [cx, cy, cz], ProbeField::SField));
                }
                if ui.button("Clear histories").clicked() {
                    probes.clear_histories();
                }
                if ui.button("Remove all").clicked() {
                    probes.clear();
                }
            });

            ui.label(format!(
                "max_history: {}  (samples / probe)",
                probes.max_history
            ));
            ui.separator();

            if probes.probes.is_empty() {
                ui.label("No probes. Select the Bifilar Pair scenario to install defaults, or 'Add @ center'.");
                return;
            }

            // --- Per-probe controls ---
            let mut to_remove = None;
            let nx_max = grid.nx.saturating_sub(1);
            let ny_max = grid.ny.saturating_sub(1);
            let nz_max = grid.nz.saturating_sub(1);

            for (idx, probe) in probes.probes.iter_mut().enumerate() {
                ui.push_id(idx, |ui| {
                    ui.horizontal(|ui| {
                        ui.label(format!("{}:", probe.label));
                        ui.label(format!(
                            "{} @ [{}, {}, {}] → {:.3e}",
                            probe.field.name(),
                            probe.position[0],
                            probe.position[1],
                            probe.position[2],
                            probe.latest(),
                        ));
                        if ui.small_button("X").clicked() {
                            to_remove = Some(idx);
                        }
                    });

                    egui::Grid::new("probe_pos").show(ui, |ui| {
                        ui.label("X:");
                        ui.add(egui::DragValue::new(&mut probe.position[0]).range(0..=nx_max).speed(1));
                        ui.label("Y:");
                        ui.add(egui::DragValue::new(&mut probe.position[1]).range(0..=ny_max).speed(1));
                        ui.label("Z:");
                        ui.add(egui::DragValue::new(&mut probe.position[2]).range(0..=nz_max).speed(1));
                        ui.end_row();

                        ui.label("Field:");
                        egui::ComboBox::from_id_salt("probe_field")
                            .selected_text(probe.field.name())
                            .show_ui(ui, |ui| {
                                for &pf in ProbeField::ALL {
                                    ui.selectable_value(&mut probe.field, pf, pf.name());
                                }
                            });
                        ui.end_row();
                    });

                    ui.separator();
                });
            }

            if let Some(idx) = to_remove {
                probes.probes.remove(idx);
            }

            // --- Time series plot ---
            ui.separator();
            ui.label("Time series (all probes overlaid):");
            draw_probe_time_series(ui, &probes.probes);

            // --- Propagation delay (bifilar 3-probe helper) ---
            if probes.probes.len() >= 3 {
                ui.separator();
                if let Some((delay_s, speed)) = estimate_propagation_delay(&probes.probes, grid.dx)
                {
                    ui.label(format!(
                        "Peak shift TX→RX: Δt = {:.3e} s  ⇒  v ≈ {:.3e} m/s",
                        delay_s, speed,
                    ));
                    ui.label(format!(
                        "c (reference): {:.3e} m/s   ratio v/c ≈ {:.3}",
                        crate::simulation::state::SimParams::C0,
                        speed / crate::simulation::state::SimParams::C0,
                    ));
                }
            }

            // --- Optional FFT of first probe ---
            if let Some(first) = probes.probes.first() {
                ui.separator();
                ui.collapsing(format!("FFT of {}", first.label), |ui| {
                    let spectrum = probe_fft(first);
                    if spectrum.is_empty() {
                        ui.label("(need ≥ 2 samples with monotonic time)");
                    } else {
                        draw_probe_spectrum(ui, &spectrum);
                    }
                });
            }
        });
}

/// Draw an overlaid line plot of all probe time series using the egui
/// painter. Each probe gets its own color and auto-scales to the global
/// (min,max) bounding box. Deliberately lightweight to avoid adding an
/// `egui_plot` dependency.
fn draw_probe_time_series(ui: &mut egui::Ui, probes: &[Probe]) {
    let desired_size = egui::vec2(ui.available_width(), 140.0);
    let (rect, _resp) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
    let painter = ui.painter_at(rect);

    // Frame
    painter.rect_stroke(
        rect,
        0.0,
        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY),
        egui::StrokeKind::Outside,
    );

    // Compute global bounds across all probes with data
    let mut t_min = f32::INFINITY;
    let mut t_max = f32::NEG_INFINITY;
    let mut v_min = f32::INFINITY;
    let mut v_max = f32::NEG_INFINITY;
    for p in probes {
        for &(t, v) in &p.history {
            if t < t_min { t_min = t; }
            if t > t_max { t_max = t; }
            if v < v_min { v_min = v; }
            if v > v_max { v_max = v; }
        }
    }
    if !t_min.is_finite() || t_max <= t_min {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "(no samples yet — run the simulation)",
            egui::FontId::monospace(11.0),
            egui::Color32::GRAY,
        );
        return;
    }
    // Guard against zero-variance series
    if (v_max - v_min).abs() < 1e-30 {
        v_min -= 1.0;
        v_max += 1.0;
    }

    let colors = [
        egui::Color32::from_rgb(255, 180, 90),
        egui::Color32::from_rgb(120, 220, 255),
        egui::Color32::from_rgb(160, 255, 160),
        egui::Color32::from_rgb(255, 120, 200),
        egui::Color32::from_rgb(220, 220, 120),
    ];

    // Zero line
    if v_min < 0.0 && v_max > 0.0 {
        let y0 = remap(0.0, v_min, v_max, rect.bottom(), rect.top());
        painter.line_segment(
            [
                egui::pos2(rect.left(), y0),
                egui::pos2(rect.right(), y0),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_rgb(70, 70, 70)),
        );
    }

    for (i, p) in probes.iter().enumerate() {
        if p.history.len() < 2 { continue; }
        let color = colors[i % colors.len()];
        let points: Vec<egui::Pos2> = p
            .history
            .iter()
            .map(|&(t, v)| {
                egui::pos2(
                    remap(t, t_min, t_max, rect.left(), rect.right()),
                    remap(v, v_min, v_max, rect.bottom(), rect.top()),
                )
            })
            .collect();
        painter.add(egui::Shape::line(points, egui::Stroke::new(1.5, color)));
    }

    // Legend + axes labels
    let font = egui::FontId::monospace(10.0);
    painter.text(
        rect.left_top() + egui::vec2(4.0, 2.0),
        egui::Align2::LEFT_TOP,
        format!("v:[{:.2e}, {:.2e}]", v_min, v_max),
        font.clone(),
        egui::Color32::LIGHT_GRAY,
    );
    painter.text(
        rect.right_bottom() + egui::vec2(-4.0, -2.0),
        egui::Align2::RIGHT_BOTTOM,
        format!("t:[{:.2e}, {:.2e}] s", t_min, t_max),
        font.clone(),
        egui::Color32::LIGHT_GRAY,
    );
    for (i, p) in probes.iter().enumerate() {
        let color = colors[i % colors.len()];
        painter.text(
            rect.left_top() + egui::vec2(4.0, 16.0 + 12.0 * i as f32),
            egui::Align2::LEFT_TOP,
            format!("■ {}", p.label),
            font.clone(),
            color,
        );
    }
}

/// Draw a magnitude spectrum as a thin bar chart. Highlights the peak bin.
fn draw_probe_spectrum(ui: &mut egui::Ui, spectrum: &[(f32, f32)]) {
    let desired_size = egui::vec2(ui.available_width(), 100.0);
    let (rect, _resp) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_stroke(
        rect,
        0.0,
        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY),
        egui::StrokeKind::Outside,
    );

    if spectrum.is_empty() { return; }

    let f_max = spectrum.last().map(|(f, _)| *f).unwrap_or(1.0);
    let m_max = spectrum.iter().map(|(_, m)| *m).fold(0.0f32, f32::max).max(1e-30);

    let mut peak_f = 0.0f32;
    let mut peak_m = 0.0f32;
    for &(f, m) in spectrum {
        if m > peak_m { peak_m = m; peak_f = f; }
        let x = remap(f, 0.0, f_max, rect.left(), rect.right());
        let y = remap(m, 0.0, m_max, rect.bottom(), rect.top());
        painter.line_segment(
            [egui::pos2(x, rect.bottom()), egui::pos2(x, y)],
            egui::Stroke::new(1.0, egui::Color32::from_rgb(120, 220, 255)),
        );
    }

    let font = egui::FontId::monospace(10.0);
    painter.text(
        rect.left_top() + egui::vec2(4.0, 2.0),
        egui::Align2::LEFT_TOP,
        format!("peak: {:.3e} Hz", peak_f),
        font.clone(),
        egui::Color32::LIGHT_GRAY,
    );
    painter.text(
        rect.right_bottom() + egui::vec2(-4.0, -2.0),
        egui::Align2::RIGHT_BOTTOM,
        format!("f_max: {:.3e} Hz", f_max),
        font,
        egui::Color32::LIGHT_GRAY,
    );
}

/// Linear remap of `x` from `[a, b]` to `[c, d]`, with `a != b` guaranteed
/// by callers.
fn remap(x: f32, a: f32, b: f32, c: f32, d: f32) -> f32 {
    let t = (x - a) / (b - a);
    c + t * (d - c)
}

/// Estimate scalar-wave propagation delay from the first three probes.
///
/// Assumes probe[0] is at the transmitter and probe[2] is at the receiver.
/// Returns (Δt, v) where Δt is the time lag of the receiver's peak-|amplitude|
/// sample relative to the transmitter's, and v is |Δx|/Δt using the positions
/// of probe[0] and probe[2] scaled by grid dx.
fn estimate_propagation_delay(probes: &[Probe], dx: f32) -> Option<(f32, f32)> {
    if probes.len() < 3 { return None; }
    let tx = &probes[0];
    let rx = &probes[2];
    let tx_peak = peak_time(&tx.history)?;
    let rx_peak = peak_time(&rx.history)?;
    let delay = rx_peak - tx_peak;
    if delay <= 0.0 { return None; }
    let d = distance_cells(tx.position, rx.position) as f32 * dx;
    Some((delay, d / delay))
}

fn peak_time(history: &VecDequeF32) -> Option<f32> {
    let mut best_t = None;
    let mut best_abs = 0.0f32;
    for &(t, v) in history {
        let a = v.abs();
        if a > best_abs {
            best_abs = a;
            best_t = Some(t);
        }
    }
    best_t
}

type VecDequeF32 = std::collections::VecDeque<(f32, f32)>;

fn distance_cells(a: [u32; 3], b: [u32; 3]) -> u32 {
    let dx = a[0].max(b[0]) - a[0].min(b[0]);
    let dy = a[1].max(b[1]) - a[1].min(b[1]);
    let dz = a[2].max(b[2]) - a[2].min(b[2]);
    // Euclidean rounded to nearest cell — probes are axis-aligned in the
    // bifilar scenario, so this collapses to |dx|.
    (((dx * dx + dy * dy + dz * dz) as f32).sqrt()).round() as u32
}
