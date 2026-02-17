// UiPlugin: egui-based control panels for the simulator.

use bevy::prelude::*;
use bevy::diagnostic::{DiagnosticsStore, FrameTimeDiagnosticsPlugin};
use bevy_egui::{egui, EguiContexts, EguiPlugin};

use crate::scenarios::dipole_radiation::{self, Scenario};
use crate::simulation::diagnostics::DiagnosticsState;
use crate::simulation::grid::SimulationGrid;
use crate::simulation::plugin::SimulationConfig;
use crate::simulation::sources::{Source, SourceConfig, SourceType};
use crate::visualization::color_maps::{ColorEncoding, ColorMap, PhasePlane};
use crate::visualization::glyphs::{GlyphConfig, GlyphField};
use crate::visualization::slices::{AllSliceStats, FieldQuantity, SliceAxis, SliceConfigs, NUM_SLICES};

/// Tracks which scenario is currently selected in the UI.
#[derive(Resource, Default)]
struct SelectedScenario {
    current: Option<Scenario>,
}

/// Plugin that adds bevy_egui and the simulation control panels.
pub struct UiPlugin;

impl Plugin for UiPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(EguiPlugin)
            .add_plugins(FrameTimeDiagnosticsPlugin::default())
            .init_resource::<SelectedScenario>()
            .add_systems(Update, (ui_side_panel, ui_source_panel, ui_diagnostics_panel, ui_slice_panel, ui_glyph_panel));
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
                    // Single-step: unpause for one frame, then re-pause
                    // (handled by the simulation step system checking a one-shot flag)
                    config.paused = false;
                    // TODO: implement single-step logic in simulation_step_system
                }
                if ui.button("Reset").clicked() {
                    if let Some(ref mut grid) = grid {
                        grid.reset();
                    }
                    sources.sources.clear();
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
                            // Apply scenario configuration
                            if let Some(ref mut grid) = grid {
                                match scenario {
                                    Scenario::DipoleRadiation => {
                                        dipole_radiation::apply_dipole_scenario(grid, &mut sources);
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
    stats: Res<AllSliceStats>,
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

            // Scale (logarithmic)
            ui.add(
                egui::Slider::new(&mut config.scale, 1e-6..=1.0)
                    .logarithmic(true)
                    .text("Scale"),
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
