mod calculation;
mod constants;
mod cross_section;
mod integrate;
mod matching;
mod potentials;
mod wave_function;
use constants::*;
use num::complex::Complex;
use potentials::FormFactor;
use pyo3::prelude::*;
use std::f64::consts::PI;

fn deg_to_rad(angles: &[f64]) -> Vec<f64> {
    // check and convert angles
    let rad_angles: Vec<f64> = angles
        .into_iter()
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

/// Elastic scattering for spin zero particles also returns rutherford.
///fn spin_zero(
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
///     r_c: f64,
///     partial_waves: i32,
///     angles: Vec<f64>,
///     r_match: f64,
///     dr: f64,
/// ) -> (Vec<f64>, Vec<f64>)
#[pyfunction]
fn spin_zero(
    a1: f64,
    m1: f64,
    z1: f64,
    a2: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    V: f64,
    r: f64,
    a: f64,
    W: f64,
    r_i: f64,
    a_i: f64,
    r_c: f64,
    partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
    par: bool,
) -> (Vec<f64>, Vec<f64>) {
    // reaction constants

    // convert to MeV
    let m1 = m1 * u_to_MeV;
    let m2 = m2 * u_to_MeV;

    // Scale the radii
    let a13 = a2.powf(1.0 / 3.0);
    let r = r * a13;
    let r_i = r_i * a13;
    let r_c = r_c * a13;

    let energy_com = energy_lab * (m2 / (m1 + m2));
    let mu = (m1 * m2) / (m1 + m2);
    let k = f64::sqrt((2.0 * mu * energy_com) / hbar.powi(2));
    let eta = ((z1 * z2) * e2) * (mu / (hbar.powi(2) * k));

    // check and convert angles
    let angles: Vec<f64> = deg_to_rad(&angles);

    // setup the grid and the potentials
    let r_grid: Vec<f64> = calculation::setup_grid(r_match, dr);
    let ff: FormFactor = calculation::setup_form_factor(
        r_grid.as_slice(),
        V,
        r,
        a,
        W,
        r_i,
        a_i,
        0.0,
        0.0,
        0.0,
        z1,
        z2,
        r_c,
        mu,
        k,
        eta,
    );

    // calculate the scattering amplitude not that ff will be moved
    let total_ampl: Vec<Complex<f64>> = if par {
        calculation::partial_waves_par(r_grid.as_slice(), ff, angles.as_slice(), partial_waves, dr)
    } else {
        calculation::partial_waves(r_grid.as_slice(), ff, angles.as_slice(), partial_waves, dr)
    };

    // cross section in mb
    let sigma: Vec<f64> =
        cross_section::diff_cross_section(angles.as_slice(), total_ampl.as_slice(), k, eta);
    let ruth: Vec<f64> = cross_section::rutherford_cs(angles.as_slice(), k, eta);
    (sigma, ruth)
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

#[pyfunction]
fn spin_half(
    a1: f64,
    m1: f64,
    z1: f64,
    a2: f64,
    m2: f64,
    z2: f64,
    energy_lab: f64,
    V: f64,
    r: f64,
    a: f64,
    W: f64,
    r_i: f64,
    a_i: f64,
    V_so: f64,
    r_so: f64,
    a_so: f64,
    r_c: f64,
    partial_waves: i32,
    angles: Vec<f64>,
    r_match: f64,
    dr: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // reaction constants

    // convert to MeV
    let m1 = m1 * u_to_MeV;
    let m2 = m2 * u_to_MeV;

    // Scale the radii
    let a13 = a2.powf(1.0 / 3.0);
    let r = r * a13;
    let r_i = r_i * a13;
    let r_so = r_so * a13;
    let r_c = r_c * a13;

    let energy_com = energy_lab * (m2 / (m1 + m2));
    let mu = (m1 * m2) / (m1 + m2);
    let k = f64::sqrt((2.0 * mu * energy_com) / hbar.powi(2));
    let eta = ((z1 * z2) * e2) * (mu / (hbar.powi(2) * k));

    // check and convert angles
    let angles: Vec<f64> = deg_to_rad(&angles);

    // setup the grid and the potentials
    let r_grid: Vec<f64> = calculation::setup_grid(r_match, dr);
    let ff: FormFactor = calculation::setup_form_factor(
        r_grid.as_slice(),
        V,
        r,
        a,
        W,
        r_i,
        a_i,
        V_so,
        r_so,
        a_so,
        z1,
        z2,
        r_c,
        mu,
        k,
        eta,
    );

    // calculate the scattering amplitude not that ff will be moved
    let (a_theta, b_theta): (Vec<Complex<f64>>, Vec<Complex<f64>>) =
        calculation::partial_waves_half_par(r_grid.as_slice(), ff, &angles, partial_waves, dr);

    // cross section in mb
    //    let sigma: Vec<f64> =
    let (sigma, pol) = cross_section::all_observables(&angles, &a_theta, &b_theta, k, eta);
    let ruth: Vec<f64> = cross_section::rutherford_cs(&angles, k, eta);
    (sigma, pol, ruth)
}

/// Lightweight optical model used in Python written in Rust.
#[pymodule]
fn dwuckman(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(spin_zero, m)?)?;
    m.add_function(wrap_pyfunction!(spin_half, m)?)?;
    Ok(())
}
