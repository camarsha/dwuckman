use crate::{
    cross_section, integrate, matching, potentials::FormFactor, wave_function::WaveFunction,
};
use num::Complex;
use rayon::prelude::*;

// Module to help reduce the length of the main loop

pub fn setup_grid(r_match: f64, h: f64) -> Vec<f64> {
    /*
    the grid has a few hidden items
    1) the first point is always 1e-10 fm

    2) later on, matching for the phase shifts
    we will use a pseudo-Wronskian method. This entails
    going to the point R_match + h. Lets go ahead and put
    the additional point in.

     */

    // number of integration points on the grid always
    // will always be greater than r_match, and always have
    // one additional point
    let points = (r_match / h).ceil() as usize + 1;

    let r_grid: Vec<f64> = (0..points)
        .map(|p| {
            let r: f64 = if p == 0 { 1e-10 } else { (p as f64) * h };
            r
        })
        .collect();

    r_grid
}

/// sets up the form factor based on the potentials that are entered
pub fn setup_form_factor(
    r_grid: &[f64],
    V: f64,
    r: f64,
    a: f64,
    W: f64,
    r_i: f64,
    a_i: f64,
    V_so: f64,
    r_so: f64,
    a_so: f64,
    z1: f64,
    z2: f64,
    r_c: f64,
    mu: f64,
    k: f64,
    eta: f64,
) -> FormFactor {
    let mut ff: FormFactor = FormFactor::new(r_grid, mu, k, eta);
    // setup the potentials
    if V != 0.0 {
        ff.add_woods_saxon(V, r, a, true);
    };
    if W != 0.0 {
        ff.add_woods_saxon(W, r_i, a_i, false);
    };
    if z1 != 0.0 {
        ff.add_coulomb(z1, z2, r_c);
    };
    if V_so != 0.0 {
        ff.add_spin_orbit(V_so, r_so, a_so);
    }
    ff.scale(k, mu);

    ff
}

/// Returns a tuple of the index and rho values at the matching radius
pub fn match_points(r_grid: &[f64], k: f64) -> (usize, f64, f64) {
    let r_idx = r_grid.len() - 2;
    let rho_r = r_grid[r_idx] * k;
    let rho_rh = r_grid[r_idx + 1] * k;
    (r_idx, rho_r, rho_rh)
}

/// performs integration and returns the total scattering amplitudes, loop is parallel
pub fn partial_waves_par(
    r_grid: &[f64],
    ff: FormFactor,
    angles: &[f64],
    num_l: i32,
    h: f64,
) -> Vec<Complex<f64>> {
    // convert l values to f64 for calculations
    let ell: Vec<f64> = (0..num_l).map(|x| x as f64).collect();

    // the grid point match at
    let (r_idx, rho_r, rho_rh) = match_points(r_grid, ff.k);

    // parallel partial wave loop
    let ell_ampl: Vec<Vec<Complex<f64>>> = ell
        .into_par_iter()
        .map(|l| {
            // create wave function
            let mut phi = WaveFunction::new(r_grid);

            // starting values for integration
            phi.setup(h, l);

            // add centrifugal term
            let re_l = ff.update_centrifugal(ff.re.as_slice(), l);

            // special case for l=1, see Melkanoff
            if l as i32 == 1 {
                phi.re[phi.start_idx - 1] = 2.0 / re_l[0];
            }

            integrate::fox_goodwin_coupled(
                h,
                re_l.as_slice(),
                ff.im.as_slice(),
                phi.re.as_mut_slice(),
                phi.im.as_mut_slice(),
                phi.start_idx,
            );

            // values for matching using Psuedo-Wronskian

            let phi_r = Complex::new(phi.re[r_idx], phi.im[r_idx]);
            let phi_rh = Complex::new(phi.re[r_idx + 1], phi.im[r_idx + 1]);

            let phase_shift = matching::phase_shift(phi_r, phi_rh, rho_r, rho_rh, ff.eta, l);
            cross_section::spin_zero_amp(angles, phase_shift, l, ff.eta, ff.k)
        })
        .collect();

    // implementation needs work
    let mut total_ampl: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];

    // sum across the ell components
    for ele in ell_ampl.iter() {
        total_ampl = total_ampl
            .iter()
            .zip(ele.iter())
            .map(|(&x, &y)| x + y)
            .collect();
    }

    total_ampl
}

/// performs integration and returns the total scattering amplitudes
pub fn partial_waves(
    r_grid: &[f64],
    ff: FormFactor,
    angles: &[f64],
    num_l: i32,
    h: f64,
) -> Vec<Complex<f64>> {
    // convert l values to f64 for calculations
    let ell: Vec<f64> = (0..num_l).map(|x| x as f64).collect();

    // the grid point match at
    let (r_idx, rho_r, rho_rh) = match_points(r_grid, ff.k);

    let mut total_ampl: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];

    // parallel partial wave loop
    for l in ell.into_iter() {
        // create wave function
        let mut phi = WaveFunction::new(r_grid);

        // starting values for integration
        phi.setup(h, l);

        // add centrifugal term, changed to make it more clear, old way worked
        // due to different scope of the for loop
        let re_l = ff.update_centrifugal(ff.re.as_slice(), l);

        // special case for l=1, see Melkanoff
        if l as i32 == 1 {
            phi.re[phi.start_idx - 1] = 2.0 / re_l[0];
        }

        integrate::fox_goodwin_coupled(
            h,
            re_l.as_slice(),
            ff.im.as_slice(),
            phi.re.as_mut_slice(),
            phi.im.as_mut_slice(),
            phi.start_idx,
        );

        // values for matching using Psuedo-Wronskian

        let phi_r = Complex::new(phi.re[r_idx], phi.im[r_idx]);
        let phi_rh = Complex::new(phi.re[r_idx + 1], phi.im[r_idx + 1]);

        let phase_shift = matching::phase_shift(phi_r, phi_rh, rho_r, rho_rh, ff.eta, l);
        println!("{} {}", phase_shift.re, phase_shift.im);

        total_ampl = total_ampl
            .iter()
            .zip(cross_section::spin_zero_amp(angles, phase_shift, l, ff.eta, ff.k).into_iter())
            .map(|(x, y)| x + y)
            .collect();
    }

    total_ampl
}

// performs integration for spin one half projectiles
pub fn partial_waves_half_par(
    r_grid: &[f64],
    ff: FormFactor,
    angles: &[f64],
    num_l: i32,
    h: f64,
) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    // check if there is spin-orbit, if not, just do the zero case
    if ff.V_so == 0.0 {
        let a_theta = partial_waves_par(r_grid, ff, angles, num_l, h);
        let b_theta = vec![Complex::new(0.0, 0.0); angles.len()];
        return (a_theta, b_theta);
    };

    // convert l values to f64 for calculations
    let ell: Vec<f64> = (0..num_l).map(|x| x as f64).collect();

    // the grid point match at
    let (r_idx, rho_r, rho_rh) = match_points(r_grid, ff.k);

    // parallel wave loop with spin up and down
    let (a_vec, b_vec): (Vec<Vec<Complex<f64>>>, Vec<Vec<Complex<f64>>>) = ell
        .into_par_iter()
        .map(|l| {
            // just hard code the two spins for now
            let mut phi_m = WaveFunction::new(r_grid);
            let mut phi_p = WaveFunction::new(r_grid);
            // starting values for integration
            phi_m.setup(h, l);
            phi_p.setup(h, l);

            // add centrifugal and spin_orbit term
            let re_l = ff.update_centrifugal(ff.re.as_slice(), l);
            let re_j_m = ff.update_spin_orbit(&re_l, l, -0.5);
            let re_j_p = ff.update_spin_orbit(&re_l, l, 0.5);

            // special case for l=1, see Melkanoff
            if l as i32 == 1 {
                phi_m.re[phi_m.start_idx - 1] = 2.0 / re_j_m[0];
                phi_p.re[phi_p.start_idx - 1] = 2.0 / re_j_p[0];
            }

            integrate::fox_goodwin_coupled(
                h,
                re_j_m.as_slice(),
                ff.im.as_slice(),
                phi_m.re.as_mut_slice(),
                phi_m.im.as_mut_slice(),
                phi_m.start_idx,
            );

            integrate::fox_goodwin_coupled(
                h,
                re_j_p.as_slice(),
                ff.im.as_slice(),
                phi_p.re.as_mut_slice(),
                phi_p.im.as_mut_slice(),
                phi_p.start_idx,
            );

            // values for matching using Psuedo-Wronskian

            let phi_r_m = Complex::new(phi_m.re[r_idx], phi_m.im[r_idx]);
            let phi_rh_m = Complex::new(phi_m.re[r_idx + 1], phi_m.im[r_idx + 1]);

            let phi_r_p = Complex::new(phi_p.re[r_idx], phi_p.im[r_idx]);
            let phi_rh_p = Complex::new(phi_p.re[r_idx + 1], phi_p.im[r_idx + 1]);

            // push to the vector
            let ps_m = matching::phase_shift(phi_r_m, phi_rh_m, rho_r, rho_rh, ff.eta, l);
            let ps_p = matching::phase_shift(phi_r_p, phi_rh_p, rho_r, rho_rh, ff.eta, l);
            println!("{} {} {} {}", ps_m.re, ps_m.im, ps_p.re, ps_m.im);
            cross_section::spin_half_ampl(angles, ps_m, ps_p, l, ff.eta, ff.k)
        })
        .collect();

    // implementation needs work
    let mut a_total: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];
    let mut b_total: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); angles.len()];

    // sums for a and b
    for ele in a_vec {
        a_total = a_total
            .iter()
            .zip(ele.iter())
            .map(|(&x, &y)| x + y)
            .collect();
    }

    for ele in b_vec {
        b_total = b_total
            .iter()
            .zip(ele.iter())
            .map(|(&x, &y)| x + y)
            .collect();
    }

    (a_total, b_total)
}
