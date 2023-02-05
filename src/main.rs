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
) -> Vec<matching::PhaseShift> {
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

    // calculate the scattering amplitude note that ff will be moved
    calculation::calc_phase_shifts(r_grid.as_slice(), ff, partial_waves, dr)
}

fn main() {
    let angles = (0..180).map(|x| x as f64).collect();

    println!(
        "{:?}",
        spin_zero(
            4.0, 4.0, 2.0, 12.0, 12.0, 6.0, 20.0, 140.0, 1.25, 0.65, 10.0, 1.15, 0.83, 1.35, 80,
            angles, 40.0, 0.1, true,
        )
    )
}
