use bevy::prelude::*;
use bevy::input::mouse::{MouseMotion, MouseWheel};
use quaternions::simulation::plugin::SimulationPlugin;
use quaternions::ui::plugin::UiPlugin;
use quaternions::visualization::plugin::VisualizationPlugin;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Quaternions — QVED Simulator".into(),
                resolution: (1280.0, 720.0).into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(SimulationPlugin)
        .add_plugins(UiPlugin)
        .add_plugins(VisualizationPlugin)
        .add_systems(Startup, setup_scene)
        .add_systems(Update, (orbit_camera_system, draw_scene_gizmos))
        .run();
}

/// Orbit camera state, stored as a resource.
#[derive(Resource)]
struct OrbitCamera {
    /// Distance from the target point.
    distance: f32,
    /// Vertical angle (radians, clamped to avoid gimbal lock).
    pitch: f32,
    /// Horizontal angle (radians).
    yaw: f32,
    /// The point the camera orbits around.
    target: Vec3,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            distance: 0.5,
            pitch: 0.6,
            yaw: 0.8,
            target: Vec3::ZERO,
        }
    }
}

/// Marker component for the orbit camera entity.
#[derive(Component)]
struct MainCamera;

fn setup_scene(mut commands: Commands) {
    // Insert orbit camera resource
    commands.insert_resource(OrbitCamera::default());

    // 3D camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(5.0, 5.0, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        MainCamera,
    ));

    // Directional light
    commands.spawn((
        DirectionalLight {
            illuminance: 10000.0,
            ..default()
        },
        Transform::from_rotation(Quat::from_euler(EulerRot::XYZ, -0.5, 0.5, 0.0)),
    ));

    // Ambient light so we can see things from all angles
    commands.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 200.0,
    });
}

/// Simple orbit camera: right-mouse-drag to orbit, scroll to zoom.
fn orbit_camera_system(
    mouse_button: Res<ButtonInput<MouseButton>>,
    mut mouse_motion: EventReader<MouseMotion>,
    mut mouse_wheel: EventReader<MouseWheel>,
    mut orbit: ResMut<OrbitCamera>,
    mut query: Query<&mut Transform, With<MainCamera>>,
) {
    // Accumulate mouse motion while right button is held
    let mut delta = Vec2::ZERO;
    for event in mouse_motion.read() {
        if mouse_button.pressed(MouseButton::Right) {
            delta += event.delta;
        }
    }

    // Apply rotation
    if delta.length_squared() > 0.0 {
        orbit.yaw -= delta.x * 0.005;
        orbit.pitch -= delta.y * 0.005;
        orbit.pitch = orbit.pitch.clamp(
            -std::f32::consts::FRAC_PI_2 + 0.05,
            std::f32::consts::FRAC_PI_2 - 0.05,
        );
    }

    // Scroll to zoom
    for event in mouse_wheel.read() {
        orbit.distance *= 1.0 - event.y * 0.1;
        orbit.distance = orbit.distance.clamp(0.05, 100.0);
    }

    // Update camera transform
    if let Ok(mut transform) = query.get_single_mut() {
        let rotation = Quat::from_euler(EulerRot::YXZ, orbit.yaw, orbit.pitch, 0.0);
        let offset = rotation * Vec3::new(0.0, 0.0, orbit.distance);
        transform.translation = orbit.target + offset;
        transform.look_at(orbit.target, Vec3::Y);
    }
}

/// Draw bounding box wireframe and coordinate axes using gizmos.
fn draw_scene_gizmos(
    mut gizmos: Gizmos,
    grid: Option<Res<quaternions::simulation::grid::SimulationGrid>>,
) {
    // Coordinate axes at origin
    let axis_len = 0.05;
    gizmos.line(Vec3::ZERO, Vec3::X * axis_len, Color::srgb(1.0, 0.2, 0.2)); // X = red
    gizmos.line(Vec3::ZERO, Vec3::Y * axis_len, Color::srgb(0.2, 1.0, 0.2)); // Y = green
    gizmos.line(Vec3::ZERO, Vec3::Z * axis_len, Color::srgb(0.3, 0.3, 1.0)); // Z = blue

    // Bounding box showing simulation domain
    if let Some(grid) = grid {
        let half = Vec3::new(
            grid.nx as f32 * grid.dx * 0.5,
            grid.ny as f32 * grid.dx * 0.5,
            grid.nz as f32 * grid.dx * 0.5,
        );
        let color = Color::srgba(0.5, 0.5, 0.5, 0.5);

        // Draw 12 edges of the bounding box centered at origin
        let corners = [
            Vec3::new(-half.x, -half.y, -half.z),
            Vec3::new(half.x, -half.y, -half.z),
            Vec3::new(half.x, half.y, -half.z),
            Vec3::new(-half.x, half.y, -half.z),
            Vec3::new(-half.x, -half.y, half.z),
            Vec3::new(half.x, -half.y, half.z),
            Vec3::new(half.x, half.y, half.z),
            Vec3::new(-half.x, half.y, half.z),
        ];

        // Bottom face
        gizmos.line(corners[0], corners[1], color);
        gizmos.line(corners[1], corners[2], color);
        gizmos.line(corners[2], corners[3], color);
        gizmos.line(corners[3], corners[0], color);
        // Top face
        gizmos.line(corners[4], corners[5], color);
        gizmos.line(corners[5], corners[6], color);
        gizmos.line(corners[6], corners[7], color);
        gizmos.line(corners[7], corners[4], color);
        // Vertical edges
        gizmos.line(corners[0], corners[4], color);
        gizmos.line(corners[1], corners[5], color);
        gizmos.line(corners[2], corners[6], color);
        gizmos.line(corners[3], corners[7], color);
    }
}
