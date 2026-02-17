// SimulationGrid: Bevy resource holding the double-buffered simulation state.
// See ARCHITECTURE.md §Data Model and §Simulation Loop.

use bevy::prelude::*;

use crate::math::fdtd;
use crate::simulation::state::{CellState, SimParams};

/// The primary simulation grid resource.
///
/// Holds double-buffered cell state for the entire 3D domain.
/// Read from `cells[current]`, write to `cells[1 - current]`, then swap.
#[derive(Resource)]
pub struct SimulationGrid {
    /// Grid dimensions along each axis.
    pub nx: u32,
    pub ny: u32,
    pub nz: u32,

    /// Uniform cell spacing (meters).
    pub dx: f32,

    /// Timestep (seconds), computed from CFL condition.
    pub dt: f32,

    /// Double-buffered cell state. Index with `current` for the read buffer.
    pub cells: [Vec<CellState>; 2],

    /// Index into `cells` for the current read buffer (0 or 1).
    pub current: usize,

    /// Current simulation time (seconds).
    pub time: f64,

    /// Step count.
    pub iteration: u64,
}

impl SimulationGrid {
    /// Allocate and initialize a simulation grid with vacuum cells.
    ///
    /// `dt` is computed from the CFL condition: dt = 0.9 * dx / (c0 * sqrt(3)).
    pub fn new(nx: u32, ny: u32, nz: u32, dx: f32) -> Self {
        let n = (nx as usize) * (ny as usize) * (nz as usize);
        let dt = fdtd::max_dt(dx, SimParams::C0);

        Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            cells: [
                vec![CellState::vacuum(); n],
                vec![CellState::vacuum(); n],
            ],
            current: 0,
            time: 0.0,
            iteration: 0,
        }
    }

    /// Total number of cells in the grid.
    #[inline]
    pub fn cell_count(&self) -> usize {
        (self.nx as usize) * (self.ny as usize) * (self.nz as usize)
    }

    /// Convert 3D coordinates to flat index.
    #[inline]
    pub fn idx(&self, x: u32, y: u32, z: u32) -> usize {
        fdtd::idx(x as usize, y as usize, z as usize, self.nx as usize, self.ny as usize)
    }

    /// Read-only access to a cell in the current (read) buffer.
    #[inline]
    pub fn cell(&self, x: u32, y: u32, z: u32) -> &CellState {
        let i = self.idx(x, y, z);
        &self.cells[self.current][i]
    }

    /// Mutable access to a cell in the current (read) buffer.
    /// Use for source injection or initialization — NOT for the update step
    /// (which should write to `cells[1 - current]`).
    #[inline]
    pub fn cell_mut(&mut self, x: u32, y: u32, z: u32) -> &mut CellState {
        let i = self.idx(x, y, z);
        &mut self.cells[self.current][i]
    }

    /// Check if (x, y, z) is an interior cell (not on any boundary face).
    #[inline]
    pub fn is_interior(&self, x: u32, y: u32, z: u32) -> bool {
        fdtd::is_interior(
            x as usize,
            y as usize,
            z as usize,
            self.nx as usize,
            self.ny as usize,
            self.nz as usize,
        )
    }

    /// Reference to the current read buffer.
    #[inline]
    pub fn read_buf(&self) -> &[CellState] {
        &self.cells[self.current]
    }

    /// Reference to the write buffer (opposite of current).
    #[inline]
    pub fn write_buf(&self) -> &[CellState] {
        &self.cells[1 - self.current]
    }

    /// Mutable reference to the write buffer.
    #[inline]
    pub fn write_buf_mut(&mut self) -> &mut [CellState] {
        &mut self.cells[1 - self.current]
    }

    /// Swap read/write buffers and advance time.
    pub fn swap_and_advance(&mut self) {
        self.current = 1 - self.current;
        self.time += self.dt as f64;
        self.iteration += 1;
    }

    /// Swap read/write buffers and advance time by a custom dt.
    pub fn swap_and_advance_with_dt(&mut self, dt: f32) {
        self.current = 1 - self.current;
        self.time += dt as f64;
        self.iteration += 1;
    }

    /// Build a SimParams snapshot of current grid state.
    pub fn sim_params(&self, extended_mode: bool) -> SimParams {
        SimParams::new(
            self.nx,
            self.ny,
            self.nz,
            self.dx,
            self.dt,
            extended_mode,
        )
    }

    /// Reset the grid to vacuum: zero all fields, reset time and iteration.
    pub fn reset(&mut self) {
        for cell in self.cells[0].iter_mut() {
            *cell = CellState::vacuum();
        }
        for cell in self.cells[1].iter_mut() {
            *cell = CellState::vacuum();
        }
        self.current = 0;
        self.time = 0.0;
        self.iteration = 0;
    }

    /// Build a SimParams snapshot with a custom dt (for dt_factor scaling).
    pub fn sim_params_with_dt(&self, extended_mode: bool, dt: f32) -> SimParams {
        SimParams::new(
            self.nx,
            self.ny,
            self.nz,
            self.dx,
            dt,
            extended_mode,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_creation() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        assert_eq!(grid.cell_count(), 512);
        assert_eq!(grid.nx, 8);
        assert!(grid.dt > 0.0);
        assert_eq!(grid.current, 0);
        assert_eq!(grid.time, 0.0);
        assert_eq!(grid.iteration, 0);
    }

    #[test]
    fn test_grid_cell_access() {
        let grid = SimulationGrid::new(4, 4, 4, 0.1);
        let cell = grid.cell(2, 2, 2);
        assert_eq!(cell.k, 1.0); // vacuum default
        assert_eq!(cell.q, [0.0; 4]);
    }

    #[test]
    fn test_grid_cell_mut() {
        let mut grid = SimulationGrid::new(4, 4, 4, 0.1);
        grid.cell_mut(1, 1, 1).q[0] = 42.0;
        assert_eq!(grid.cell(1, 1, 1).q[0], 42.0);
    }

    #[test]
    fn test_grid_is_interior() {
        let grid = SimulationGrid::new(8, 8, 8, 0.01);
        assert!(!grid.is_interior(0, 0, 0));
        assert!(!grid.is_interior(7, 4, 4));
        assert!(grid.is_interior(4, 4, 4));
        assert!(grid.is_interior(1, 1, 1));
    }

    #[test]
    fn test_grid_swap() {
        let mut grid = SimulationGrid::new(4, 4, 4, 0.1);
        assert_eq!(grid.current, 0);
        let dt = grid.dt as f64;
        grid.swap_and_advance();
        assert_eq!(grid.current, 1);
        assert_eq!(grid.iteration, 1);
        assert!((grid.time - dt).abs() < 1e-15);
        grid.swap_and_advance();
        assert_eq!(grid.current, 0);
        assert_eq!(grid.iteration, 2);
    }

    #[test]
    fn test_grid_double_buffer_isolation() {
        let mut grid = SimulationGrid::new(4, 4, 4, 0.1);
        // Write to current buffer
        grid.cell_mut(1, 1, 1).q[0] = 99.0;
        // Write buffer should still be zeros
        let write_idx = grid.idx(1, 1, 1);
        assert_eq!(grid.write_buf()[write_idx].q[0], 0.0);
    }
}
