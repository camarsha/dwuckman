// Deals with everything related to wave functions
pub struct WaveFunction {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub start_idx: usize,
}

impl WaveFunction {
    pub fn new(grid: &[f64]) -> WaveFunction {
        let grid_size = grid.len();
        WaveFunction {
            re: vec![0.0; grid_size],
            im: vec![0.0; grid_size],
            start_idx: 0,
        }
    }

    pub fn setup(&mut self, h: f64, l: f64) {
        // Simple implementation of a l dependent advance
        // of the starting integration point and step
        self.start_idx = (f64::sqrt((l * (l + 1.0)) / 12.0) + 1.0) as usize;
        self.re[self.start_idx] = h.powi(l as i32 + 1);
    }
}
