use num::complex::Complex;
use std::f64::consts::PI;
mod cross_section;
mod integrate;
mod matching;
mod potentials;
mod wave_function;
use potentials::FormFactor;
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;
use wave_function::WaveFunction;
extern crate rayon;
use rayon::prelude::*;

fn main() {
    // constants for the calculation
    let u_to_MeV = 931.49410242;
    let hbar: f64 = 197.3269804;
    let e_charge = 1.602176634e-19;
    let e2 = 1.439976;
    let E_lab = 30.0;

    // reaction constants
    let m1 = 4.0 * u_to_MeV;
    let m2 = 12.0 * u_to_MeV;
    let z1 = 2.0;
    let z2 = 6.0;
    let E_cm = E_lab * (m2 / (m1 + m2));
    let mu = (m1 * m2) / (m1 + m2);
    let k = f64::sqrt((2.0 * mu * E_cm) / hbar.powi(2));
    let eta = ((z1 * z2) * e2) * (mu / (hbar.powi(2) * k));

    // partial waves
    let num_l: usize = 40;
    let ell: Vec<f64> = (0..num_l).map(|x| x as f64).collect();

    // setup the grid
    let h = 0.1; // units of fm
    let r_match = 40.0; //
    let points = (r_match / h) as usize + 1;
    let mut r_grid: Vec<f64> = (0..points).map(|x| h * x as f64).collect();

    r_grid[0] = 1e-10;

    let mut ff: FormFactor = FormFactor::new(&r_grid);
    let a13 = (12.0_f64).powf(1.0 / 3.0);

    let V = 150.0;
    let r = 1.25 * a13;
    let a = 0.65;
    let rc = 1.3 * a13;

    // setup the potentials
    ff.add_woods_saxon(V, r, a, true);
    ff.add_woods_saxon(10.0, r, 0.35, false);
    ff.add_coulomb(z1, z2, rc);
    ff.scale(k, mu);

    let mut total_ff: Vec<FormFactor> = Vec::new();

    // set up the angles for the cross sections
    let mut angles: Vec<f64> = (0..180).map(|x| x as f64).collect();
    angles[0] = 1e-4;
    let angles_rad: Vec<f64> = angles.iter().map(|x| x * PI / 180.0).collect();

    // Amplitude Sum
    //let A: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];

    // matching radii
    let Rh_idx = r_grid.len() - 1;
    let R_idx = r_grid.len() - 2;
    let rho_R = r_grid[R_idx] * k;
    let rho_Rh = r_grid[Rh_idx] * k;

    // short hand for complex numbers
    let j: Complex<f64> = Complex::i();

    // unpack the struct now that all of the
    // l-independent potential has been calculated
    let FormFactor { re, im, grid } = ff;

    // parallel partial wave loop
    let ell_ampl: Vec<Vec<Complex<f64>>> = ell
        .into_par_iter()
        .map(|l| {
            // create wave function
            let mut phi = WaveFunction::new(&r_grid);

            // starting values for integration
            phi.setup(h, l);

            // add centrifugal term
            let re = FormFactor::add_centrifugal(&r_grid, l, re.as_slice());

            // special case for l=1, see Melkanoff
            if l as i32 == 1 {
                phi.re[phi.start_idx - 1] = 2.0 / re[0];
            }

            integrate::fox_goodwin_coupled(
                h,
                &re[..],
                &im[..],
                phi.re.as_mut_slice(),
                phi.im.as_mut_slice(),
                phi.start_idx,
            );

            // values for matching using Psuedo-Wronskian

            let phih = Complex::new(phi.re[R_idx + 1], phi.im[R_idx + 1]);
            let phi0 = Complex::new(phi.re[R_idx], phi.im[R_idx]);
            let phi1 = Complex::new(phi.re[R_idx - 1], phi.im[R_idx - 1]);
            let phi2 = Complex::new(phi.re[R_idx - 2], phi.im[R_idx - 2]);

            let phase_shift = matching::phase_shift(phi0, phih, rho_R, rho_Rh, eta, l);

            //      println!("{} {} {}", l, phase_shift.re, phase_shift.im);
            //Calculate and print S matrix
            let S_mat: Complex<f64> = (2.0 * j * phase_shift).exp();

            let a = cross_section::spin_zero_amp(&angles_rad, phase_shift, l, eta, k);
            a
        })
        .collect();

    // sum the amplitudes
    let mut A: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];

    for ele in ell_ampl.iter() {
        A = A.iter().zip(ele.iter()).map(|(&x, &y)| x + y).collect();
    }

    let sigma: Vec<f64> = cross_section::diff_cross_section(&angles_rad, eta, k, A.as_slice());

    // I like ratio to rutherford
    let R = cross_section::rutherford_cs(&angles_rad, eta, k);
    let RR: Vec<f64> = sigma.iter().zip(R.iter()).map(|(&x, &y)| x / y).collect();
    // output to text file
    for (sig, ang) in sigma.iter().zip(angles) {
        println!("{} {}", ang, sig);
    }
}
