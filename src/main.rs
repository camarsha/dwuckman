mod calculation;
mod constants;
mod cross_section;
mod integrate;
mod matching;
mod potentials;
mod wave_function;
use constants::*;
//use core::slice::SlicePattern;
use num::complex::Complex;
use potentials::FormFactor;
use pyo3::prelude::*;
use std::f64::consts::PI;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

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
) -> (f64, Vec<f64>, Vec<f64>) {
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

    println!("k = {} ; eta = {}", k, eta);

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

fn main() {
    // output file
    let mut out_file = BufWriter::new(File::create("calc.csv").unwrap());
    writeln!(out_file, "theta,sigma,ruth,tot").unwrap();

    let angles: Vec<f64> = (0..180).map(|x| x as f64).collect();
    //    let nangles: usize = angles.len();

    let (tot, diff, ruth) = spin_zero(
        4.0,
        4.0,
        2.0,
        88.0,
        88.0,
        38.0,
        21.0,
        140.0,
        1.25,
        0.65,
        10.0,
        1.15,
        0.83,
        1.3,
        40,
        angles.clone(),
        50.0,
        0.1,
        true,
    );

    println!("Total Cross Section = {}", tot);

    for i in 0..angles.len() {
        writeln!(out_file, "{},{},{},{}", angles[i], diff[i], ruth[i], tot).unwrap();
    }
}
