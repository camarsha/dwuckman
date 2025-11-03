mod calculation;
mod constants;
mod cross_section;
mod integrate;
mod matching;
mod potentials;
mod wave_function;
use constants::*;
use potentials::FormFactor;
use pyo3::prelude::*;
use std::f64::consts::PI;

#[inline(always)]
fn deg_to_rad(angles: &[f64]) -> Vec<f64> {
    // check and convert angles
    let rad_angles: Vec<f64> = angles
        .iter()
        .map(|&x| {
            let y: f64 = if x < 1e-4 {
                1e-2 * PI / 180.0 // 1e-4 is the smallest angle we will consider
            } else {
                x * PI / 180.0
            };
            y
        })
        .collect();
    rad_angles
}

#[inline(always)]
fn mass_constants(m1: f64, m2: f64) -> (f64, f64, f64, f64, f64) {
    let a1 = m1.round();
    let a2 = m2.round();
    // convert to MeV
    let m1_mev = m1 * U_TO_MEV;
    let m2_mev = m2 * U_TO_MEV;
    // Scale the radii
    let a13 = a2.powf(1.0 / 3.0);
    (m1_mev, m2_mev, a1, a2, a13)
}

#[inline(always)]
fn com_energy(energy_lab: f64, m1: f64, m2: f64) -> f64 {
    energy_lab * (m2 / (m1 + m2))
}

#[inline(always)]
/// Calculate the reduced mass. IT IS EXPECT THE MASSES ARE IN MeV!
fn reduced_mass(m1: f64, m2: f64) -> f64 {
    (m1 * m2) / (m1 + m2)
}

#[inline(always)]
fn wave_number(energy_com: f64, mu: f64) -> f64 {
    f64::sqrt((2.0 * mu * energy_com) / HBAR.powi(2))
}

#[inline(always)]
fn coulomb_wf_eta(z1: f64, z2: f64, mu: f64, k: f64) -> f64 {
    ((z1 * z2) * E2) * (mu / (HBAR.powi(2) * k))
}

#[inline(always)]
fn common_spin_zero(
    m1: f64,
    z1: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    pot_params: Vec<f64>,
    _partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
) -> (Vec<f64>, Vec<f64>, FormFactor, f64, f64, f64) {
    let (m1_mev, m2_mev, _a1, _a2, a13) = mass_constants(m1, m2);

    // energy constants
    let energy_com = com_energy(energy_lab, m1, m2);
    let mu = reduced_mass(m1_mev, m2_mev);
    let k = wave_number(energy_com, mu);
    let eta = coulomb_wf_eta(z1, z2, mu, k);

    // check and convert angles
    let angles: Vec<f64> = deg_to_rad(&angles);

    // setup the grid and the potentials
    let r_grid: Vec<f64> = calculation::setup_grid(r_match, dr);
    let ff: FormFactor = calculation::setup_form_factor(
        r_grid.as_slice(),
        pot_params.as_slice(),
        a13,
        z1,
        z2,
        mu,
        k,
        eta,
    );
    (angles, r_grid, ff, mu, eta, k)
}
#[pyfunction]
fn phase_shift_spin_zero(
    m1: f64,
    z1: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    pot_params: Vec<f64>,
    partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (_angles, r_grid, ff, _mu, _eta, _k) = common_spin_zero(
        m1,
        z1,
        m2,
        z2,
        energy_lab,
        pot_params,
        partial_waves,
        angles,
        r_match,
        dr,
    );
    let ps = calculation::calc_phase_shifts(r_grid.as_slice(), ff, partial_waves, dr);
    let l = ps.iter().map(|x| x.l).collect();
    let re = ps.iter().map(|x| x.val.re).collect();
    let im = ps.iter().map(|x| x.val.im).collect();
    (l, re, im)
}

#[pyfunction]
fn wave_function_spin_zero(
    m1: f64,
    z1: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    pot_params: Vec<f64>,
    ell: i32,
    partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (_angles, r_grid, ff, _mu, _eta, _k) = common_spin_zero(
        m1,
        z1,
        m2,
        z2,
        energy_lab,
        pot_params,
        partial_waves,
        angles,
        r_match,
        dr,
    );
    let wf = calculation::calc_wave_function(&r_grid, ff, ell, dr);
    (r_grid, wf.re, wf.im)
}

/// Elastic scattering for spin zero particles also returns rutherford.
///fn spin_zero(
///     m1: f64,
///     z1: f64,
///     m2: f64,
///     z2: f64,
///     energy_lab: f64,
///     V: f64,
///     r: f64,
///     a: f64,
///     W: f64,
///     r_i: f64,
///     a_i: f64,
///     r_c: f64,
///     partial_waves: i32,
///     angles: Vec<f64>,
///     r_match: f64,
///     dr: f64,
/// ) -> (Vec<f64>, Vec<f64>)
#[allow(non_snake_case, clippy::too_many_arguments)]
#[pyfunction]
fn spin_zero(
    m1: f64,
    z1: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    pot_params: Vec<f64>,
    partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
) -> (f64, Vec<f64>, Vec<f64>) {
    // common parameters that are needed.
    let (angles, r_grid, ff, _mu, eta, k) = common_spin_zero(
        m1,
        z1,
        m2,
        z2,
        energy_lab,
        pot_params,
        partial_waves,
        angles,
        r_match,
        dr,
    );
    // calculate the scattering amplitude note that ff will be moved
    let ps = calculation::calc_phase_shifts(r_grid.as_slice(), ff, partial_waves, dr);
    let mel_coeff = cross_section::melkanoff_coeff(ps.as_slice());
    let tot_cs = cross_section::cross_section_spin_zero(ps.as_slice(), mel_coeff.as_slice(), k);
    let diff_cs: Vec<f64> = cross_section::diff_cross_spin_zero(
        angles.as_slice(),
        ps.as_slice(),
        mel_coeff.as_slice(),
        k,
        eta,
    );
    let ruth: Vec<f64> = cross_section::rutherford_cs(angles.as_slice(), k, eta);
    (tot_cs, diff_cs, ruth)
}

/// #[pyfunction]
/// fn spin_half(
///     a1: f64,
///     m1: f64,
///     z1: f64,
///     a2: f64,
///     m2: f64,
///     z2: f64,
///     energy_lab: f64,
///     V: f64,
///     r: f64,
///     a: f64,
///     W: f64,
///     r_i: f64,
///     a_i: f64,
///     V_so: f64,
///     r_so: f64,
///     a_so: f64,
///     r_c: f64,
///     partial_waves: i32,
///     angles: Vec<f64>,
///     r_match: f64,
///     dr: f64,
/// ) -> (Vec<f64>, Vec<f64>, Vec<f64>)

// #[pyfunction]
// fn spin_half(
//     a1: f64,
//     m1: f64,
//     z1: f64,
//     a2: f64,
//     m2: f64,
//     z2: f64,
//     energy_lab: f64,
//     V: f64,
//     r: f64,
//     a: f64,
//     W: f64,
//     r_i: f64,
//     a_i: f64,
//     V_so: f64,
//     r_so: f64,
//     a_so: f64,
//     r_c: f64,
//     partial_waves: i32,
//     angles: Vec<f64>,
//     r_match: f64,
//     dr: f64,
// ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
//     // reaction constants

//     // convert to MeV
//     let m1 = m1 * u_to_MeV;
//     let m2 = m2 * u_to_MeV;

//     // Scale the radii
//     let a13 = a2.powf(1.0 / 3.0);
//     let r = r * a13;
//     let r_i = r_i * a13;
//     let r_so = r_so * a13;
//     let r_c = r_c * a13;

//     let energy_com = energy_lab * (m2 / (m1 + m2));
//     let mu = (m1 * m2) / (m1 + m2);
//     let k = f64::sqrt((2.0 * mu * energy_com) / hbar.powi(2));
//     let eta = ((z1 * z2) * e2) * (mu / (hbar.powi(2) * k));

//     // check and convert angles
//     let angles: Vec<f64> = deg_to_rad(&angles);

//     // setup the grid and the potentials
//     let r_grid: Vec<f64> = calculation::setup_grid(r_match, dr);
//     let ff: FormFactor = calculation::setup_form_factor(
//         r_grid.as_slice(),
//         V,
//         r,
//         a,
//         W,
//         r_i,
//         a_i,
//         V_so,
//         r_so,
//         a_so,
//         z1,
//         z2,
//         r_c,
//         mu,
//         k,
//         eta,
//     );

//     // calculate the scattering amplitude not that ff will be moved
//     let (a_theta, b_theta): (Vec<Complex<f64>>, Vec<Complex<f64>>) =
//         calculation::partial_waves_half_par(r_grid.as_slice(), ff, &angles, partial_waves, dr);

//     // cross section in mb
//     //    let sigma: Vec<f64> =
//     let (sigma, pol) = cross_section::all_observables(&angles, &a_theta, &b_theta, k, eta);
//     let ruth: Vec<f64> = cross_section::rutherford_cs(&angles, k, eta);
//     (sigma, pol, ruth)
// }

/// Lightweight optical model used in Python written in Rust.
#[pymodule]
fn dwuckman(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(spin_zero, m)?)?;
    m.add_function(wrap_pyfunction!(phase_shift_spin_zero, m)?)?;
    m.add_function(wrap_pyfunction!(wave_function_spin_zero, m)?)?;
    Ok(())
}
