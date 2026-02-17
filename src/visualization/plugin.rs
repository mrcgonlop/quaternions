// Bevy plugin: registers all visualization systems.
// Runs AFTER SimulationSet so that DiagnosticsState is up-to-date.

use bevy::prelude::*;

use crate::simulation::plugin::SimulationSet;
use crate::visualization::glyphs::{self, GlyphConfig};
use crate::visualization::slices::{self, AllSliceStats, SliceConfigs};

pub struct VisualizationPlugin;

impl Plugin for VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SliceConfigs>()
            .init_resource::<AllSliceStats>()
            .init_resource::<GlyphConfig>()
            .add_systems(
                Update,
                (
                    slices::initialize_slice_position,
                    slices::manage_slice_entity,
                    // apply_deferred ensures entities spawned by manage_slice_entity
                    // are available for queries in update_slice_texture.
                    apply_deferred,
                    slices::update_slice_texture,
                    slices::update_slice_transform,
                    glyphs::draw_glyph_arrows,
                )
                    .chain()
                    .after(SimulationSet),
            );
    }
}
