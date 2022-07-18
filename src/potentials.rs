use crate::constants::*;

/// A struct to hold all currently defined potential parameters.
/// The alternative is to start passing 20+ arguments to functions.
pub struct PotParams {}

/// Woods-Saxon form factor. Given a positive V, returns -f(V, r, a).
pub fn woods_saxon(x: &[f64], V: f64, r: f64, a: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
    for (i, ele) in x.iter().enumerate() {
        result[i] = -V / (1.0 + f64::exp((ele - r) / a));
    }
    result
}

/// Derivative Woods-Saxon. Given a positive V, returns -f(V, r, a).
/// Matches convention that is written -4a*df/dr
pub fn der_woods_saxon(x: &[f64], V: f64, r: f64, a: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
    for (i, ele) in x.iter().enumerate() {
        result[i] = -4.0 * (V * f64::exp((ele - r) / a)) / (1.0 + f64::exp((ele - r) / a)).powi(2);
    }
    result
}

pub fn coulomb(x: &[f64], z1: f64, z2: f64, rc: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];
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
    let l2 = l * (l + 1.0);
    for (i, &ele) in x.iter().enumerate() {
        result[i] = l2 / ele.powi(2);
    }
    result
}

pub fn spin_orbit(x: &[f64], V: f64, r: f64, a: f64, l: f64, s: f64, mu: f64) -> Vec<f64> {
    let mut result = vec![0.0; x.len()];

    let j = l + s; // s is assumed to have a sign
    let j_term = j * (j + 1.0);
    let l_term = l * (l + 1.0);
    let s_term = s.abs() * (s.abs() + 1.0); // s needs to be positive now
    let spin_term = (j_term - l_term - s_term) / 2.0;
    //scaling
    let c: f64 = -(2.0 * mu) / hbar.powi(2);

    for (i, &ele) in x.iter().enumerate() {
        result[i] = 2.0 * spin_term * (V / ele * f64::exp((ele - r) / a))
            / (1.0 + f64::exp((ele - r) / a)).powi(2)
            * c;
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
    pub mu: f64,
    pub k: f64,
    pub eta: f64,
    pub V_so: f64,
    pub r_so: f64,
    pub a_so: f64,
}

impl FormFactor {
    pub fn new(grid: &[f64], mu: f64, k: f64, eta: f64) -> FormFactor {
        // grid of radii that each potential will be evaluated on
        let grid_size = grid.len();
        FormFactor {
            re: vec![0.0; grid_size],
            im: vec![0.0; grid_size],
            grid: grid.to_vec(),
            mu,
            k,
            eta,
            V_so: 0.0,
            r_so: 0.0,
            a_so: 0.0,
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

    pub fn add_der_woods_saxon(&mut self, V: f64, r: f64, a: f64, re: bool) {
        let temp: Vec<f64> = der_woods_saxon(&self.grid, V, r, a);
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

    /// This one just simply initializes the spin orbit parameters
    pub fn add_spin_orbit(&mut self, V: f64, r: f64, a: f64) {
        self.V_so = V;
        self.r_so = r;
        self.a_so = a;
    }

    pub fn scale(&mut self, k: f64, mu: f64) {
        // apply the proper scaling to the l-independent parts
        let c: f64 = -(2.0 * mu) / hbar.powi(2);

        for i in 0..self.re.len() {
            self.re[i] = -k.powi(2) - self.re[i] * c;
            self.im[i] = -c * self.im[i];
        }
    }

    pub fn update_centrifugal(&self, l: f64) -> Vec<f64> {
        let mut temp: Vec<f64> = centrifugal(self.grid.as_slice(), l);
        for i in 0..temp.len() {
            temp[i] += self.re[i];
        }
        temp
    }

    pub fn update_spin_orbit(&self, l: f64, s: f64) -> Vec<f64> {
        let mut temp: Vec<f64> = spin_orbit(
            self.grid.as_slice(),
            l,
            s,
            self.V_so,
            self.r_so,
            self.a_so,
            self.mu,
        );
        for i in 0..temp.len() {
            // need to properly scale this as well
            temp[i] += self.re[i];
        }
        temp
    }
}
