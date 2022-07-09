use num::complex::Complex;

pub fn woods_saxon(x: &[f64], V: f64, r: f64, a: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
    for (i, ele) in x.iter().enumerate() {
        result[i] = -V / (1.0 + f64::exp((ele - r) / a));
    }
    result
}

pub fn coulomb(x: &[f64], z1: f64, z2: f64, rc: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
    let e2 = 1.439976;
    for (i, &ele) in x.iter().enumerate() {
        if ele < rc {
            result[i] = (z1 * z2 * e2) / (2.0 * rc) * (3.0 - (x[i]).powi(2) / rc.powi(2));
        } else {
            result[i] = (z1 * z2 * e2) / x[i];
        }
    }
    result
}

pub fn centrifugal(x: &[f64], l: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
    let l2 = (l * (l + 1.0));
    for (i, &ele) in x.iter().enumerate() {
        result[i] = l2 / ele.powi(2);
    }
    result
}

pub fn add_pot(v1: &mut [f64], v2: &[f64]) {
    for i in 0..v1.len() {
        v1[i] += v2[i];
    }
}

// Hold the form factor (i.e sum of the potentials) for a given l.
pub struct FormFactor {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub grid: Vec<f64>,
}

impl FormFactor {
    pub fn new(grid: &[f64]) -> FormFactor {
        // grid of radii that each potential will be evaluated on
        let grid_size = grid.len();
        FormFactor {
            re: vec![0.0; grid_size],
            im: vec![0.0; grid_size],
            grid: grid.to_vec(),
        }
    }

    pub fn add_woods_saxon(&mut self, V: f64, r: f64, a: f64, re: bool) {
        let temp: Vec<f64> = woods_saxon(&self.grid, V, r, a);
        if re {
            add_pot(self.re.as_mut_slice(), temp.as_slice());
        } else {
            add_pot(self.im.as_mut_slice(), temp.as_slice());
        }
    }

    pub fn add_coulomb(&mut self, z1: f64, z2: f64, rc: f64) {
        let temp: Vec<f64> = coulomb(&self.grid, z1, z2, rc);
        add_pot(self.re.as_mut_slice(), temp.as_slice());
    }

    pub fn scale(&mut self, k: f64, mu: f64) {
        // apply the proper scaling to the l-independent parts
        let hbar: f64 = 197.3269804;
        let c: f64 = -(2.0 * mu) / hbar.powi(2);

        for i in 0..self.re.len() {
            self.re[i] = -k.powi(2) - self.re[i] * c;
            self.im[i] = -c * self.im[i];
        }
    }

    pub fn add_centrifugal(grid: &[f64], l: f64, re: &[f64]) -> Vec<f64> {
        let mut temp: Vec<f64> = centrifugal(grid, l);
        for i in 0..temp.len() {
            temp[i] += re[i];
        }
        temp
    }
}
