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
    if r_c != 0.0 {
        ff.add_coulomb(z1, z2, r_c);
    };
    ff.scale(k, mu);

    ff
}

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
            // two phase shifts per l now.
            let mut phase_shifts: Vec<Complex<f64>> = vec![];
            // iterate over the spins
            for ele in [-0.5, 0.5] {
                let mut phi = WaveFunction::new(r_grid);
                // starting values for integration
                phi.setup(h, l);

                // add centrifugal term
                let re_l = ff.update_centrifugal(ff.re.as_slice(), l);
                let re_j = ff.update_spin_orbit(re_l.as_slice(), l, ele);

                // special case for l=1, see Melkanoff
                if l as i32 == 1 {
                    phi.re[phi.start_idx - 1] = 2.0 / re_j[0];
                }

                integrate::fox_goodwin_coupled(
                    h,
                    re_j.as_slice(),
                    ff.im.as_slice(),
                    phi.re.as_mut_slice(),
                    phi.im.as_mut_slice(),
                    phi.start_idx,
                );

                // values for matching using Psuedo-Wronskian

                let phi_r = Complex::new(phi.re[r_idx], phi.im[r_idx]);
                let phi_rh = Complex::new(phi.re[r_idx + 1], phi.im[r_idx + 1]);

                // push to the vector
                phase_shifts.push(matching::phase_shift(
                    phi_r, phi_rh, rho_r, rho_rh, ff.eta, l,
                ));
            }

            cross_section::spin_half_ampl(angles, phase_shifts[0], phase_shifts[1], l, ff.eta, ff.k)
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
