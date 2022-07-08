use num::complex::Complex;
use rgsl::bessel::{jl, yl};
use rgsl::coulomb::wave_FG_e;
use rgsl::{Result, Value};
use std::f64::consts::E;

pub fn coulomb_functions(rho: f64, eta: f64, l: f64) -> Vec<f64> {
    let mut exp_F = 0.0_f64;
    let mut exp_G = 0.0_f64;
    let mut overlow: Value = Value::OverFlow;

    // Get the coulomb functions at rho
    let (overflow, mut F, mut Fp, mut G, mut Gp) =
        wave_FG_e(eta, rho, l, 0, &mut exp_F, &mut exp_G);

    // deal with potential overflow
    if !overflow.is_success() {
        F.val = F.val * exp_F.powf(E);
        G.val = G.val * exp_G.powf(E);
        Fp.val = Fp.val * exp_F.powf(E);
        Gp.val = Gp.val * exp_F.powf(E);
        println! {"Overflow in coulomb wave functions!"}
    }

    vec![F.val, Fp.val, G.val, Gp.val]
}

pub fn spherical_bessel_functions(rho: f64, l: f64) -> Vec<f64> {
    /* returns the values and values of the derivatives of the spherical Bessel
    (i.e regular and irregular, sometimes called Neumann functions) at the
    values of rho and l.

    */

    let mut exp_F = 0.0_f64;
    let mut exp_G = 0.0_f64;
    let mut overlow: Value = Value::OverFlow;

    let int_l = l as i32;
    let l2 = int_l + 1;

    // regular
    let j = jl(int_l, rho);
    let j2 = jl(l2, rho);
    // from the Wronskian
    let dj = (l + 1.0) / rho * j - j2;

    // irregular
    let y = yl(int_l, rho);
    let y2 = yl(l2, rho);
    // from the Wronskian
    let dy = (l + 1.0) / rho * y - y2;

    vec![j, dj, y, dy]
}

pub fn in_out_coulomb_functions(rho: f64, eta: f64, l: f64) -> Vec<Complex<f64>> {
    // Convert coulomb function values to the associated hankel functions

    let r: Vec<f64>;

    r = coulomb_functions(rho, eta, l);

    let mut H_plus = Complex::new(r[0], r[2]);
    let mut H_minus = Complex::new(r[0], -r[2]);
    let mut dH_plus = Complex::new(r[1], r[3]);
    let mut dH_minus = Complex::new(r[1], -r[3]);

    vec![H_plus, dH_plus, H_minus, dH_minus]
}

pub fn in_out_bessel_functions(rho: f64, l: f64) -> Vec<Complex<f64>> {
    // Form the incoming and outgoing linear independent solutions
    let r: Vec<f64>;

    r = spherical_bessel_functions(rho, l);

    let mut u_plus = Complex::new(r[0], r[2]);
    let mut u_minus = Complex::new(r[0], -r[2]);
    let mut du_plus = Complex::new(r[1], r[3]);
    let mut du_minus = Complex::new(r[1], -r[3]);

    vec![u_plus, du_plus, u_minus, du_minus]
}

pub fn phase_shift(
    phi_R: Complex<f64>,
    phi_Rh: Complex<f64>,
    rho_R: f64,
    rho_Rh: f64,
    eta: f64,
    l: f64,
) -> Complex<f64> {
    /* Calculate the nuclear phase shift using the psuedo-Wronskian
    suggested by Carl.
     */

    let fun_R: Vec<f64>; // asymptotic functions at R
    let fun_Rh: Vec<f64>; // asymptotic functions at Rh

    // calculate the appropriate functions for either a charged or neutral particle
    match (eta.ceil() as i32) {
        0 => {
            // r
            fun_R = spherical_bessel_functions(rho_R, l);
            // r + h
            fun_Rh = spherical_bessel_functions(rho_Rh, l);
        }
        _ => {
            // r
            fun_R = coulomb_functions(rho_R, eta, l);
            // r + h
            fun_Rh = coulomb_functions(rho_Rh, eta, l);
        }
    }

    // Do the matching
    let num = (phi_R * fun_Rh[0]) - (phi_Rh * fun_R[0]);
    let denom = (phi_R * fun_Rh[2]) - (phi_Rh * fun_R[2]);

    (num / denom).atan()
}

pub fn backward_difference(
    phi: Complex<f64>,
    phi1: Complex<f64>,
    phi2: Complex<f64>,
    h: f64,
) -> Complex<f64> {
    (3.0 * phi - 4.0 * phi1 + phi2) / (2.0 * h)
}

pub fn phase_shift_derivative(
    phi: Complex<f64>,
    phi1: Complex<f64>,
    phi2: Complex<f64>,
    r_match: f64,
    k: f64,
    eta: f64,
    l: f64,
    h: f64,
) -> Complex<f64> {
    let fun_R: Vec<Complex<f64>>; // asymptotic functions at R
    let rho_R = r_match * k;
    let d_phi = backward_difference(phi, phi1, phi2, h);

    // calculate the appropriate functions for either a charged or neutral particle
    match (eta.ceil() as i32) {
        0 => {
            // r
            fun_R = in_out_bessel_functions(rho_R, l);
        }
        _ => {
            // r

            fun_R = in_out_coulomb_functions(rho_R, eta, l);
        }
    }

    let R_matrix: Complex<f64> = (1.0 / r_match) * (phi / d_phi);
    let S_matrix: Complex<f64> =
        (fun_R[2] - (rho_R * fun_R[3] * R_matrix)) / (fun_R[0] - (rho_R * fun_R[1] * R_matrix));
    S_matrix.ln() / (2.0 * Complex::new(0.0, 1.0))
}
