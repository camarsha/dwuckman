use crate::matching;
use num::complex::{Complex, Complex64};
use rgsl::gamma_beta::gamma::lngamma_complex_e;
use rgsl::legendre::associated_polynomials::legendre_Plm;
use rgsl::legendre::polynomials::legendre_Pl;
use std::f64::consts::PI;

pub fn coulomb_phase_shift(l: f64, eta: f64) -> f64 {
    let zr = l + 1.0;
    let zi = eta;

    let (_lnr, arg) = lngamma_complex_e(zr, zi).unwrap_or_else(|_| {
        panic!(
            "Failed calculation of colomb phase shift for l={:.3}, eta={:.3}",
            l, eta
        );
    });
    // We just care about arg
    arg.val
}

pub fn coulomb_ampl(angles: &[f64], k: f64, eta: f64) -> Vec<Complex<f64>> {
    // This is the point-Coulomb scattering amplitude, i.e it will give you the Rutherford cross section.
    let sig_0 = coulomb_phase_shift(0.0, eta); // l = 0 phase shift
    let mut f = vec![Complex::new(0.0, 0.0); angles.len()];

    for (i, &a) in angles.iter().enumerate() {
        let denom = Complex::new(2.0 * k * (a / 2.0).sin().powi(2), 0.0);
        let num = Complex::new(0.0, -((a / 2.0).sin().powi(2).ln()) * eta + (2.0 * sig_0));
        f[i] = (-eta / denom) * num.exp();
    }

    f
}

pub fn rutherford_cs(angles: &[f64], k: f64, eta: f64) -> Vec<f64> {
    // Calculate the Rutherford differential cross section in mb/sr
    let f: Vec<Complex<f64>> = coulomb_ampl(angles, k, eta);
    let ruth: Vec<f64> = f.iter().map(|x| 10.0 * x.norm_sqr()).collect();
    ruth
}

/// Calculate C_l = -i/2 * [exp(2 * i * delta_l) - 1]
pub fn melkanoff_coeff(phase_shifts: &[matching::PhaseShift]) -> Vec<Complex64> {
    phase_shifts
        .iter()
        .map(|ps| (-Complex::i() / 2.0) * ((2.0_f64 * Complex::i() * ps.val).exp() - 1.0))
        .collect()
}

/// Calculate the total reaction cross section for spin zero particles
/// Update: This now is used to check for convergence as well
/// the method for checking convergence.
pub fn cross_section_spin_zero(
    phase_shifts: &[matching::PhaseShift],
    mel_coeff: &[Complex64],
    k: f64,
    absend: f64,
) -> (usize, f64) {
    let coeff = 10.0 * (4.0 * PI) / (k.powi(2));
    let cs_l: Vec<f64> = mel_coeff
        .iter()
        .zip(phase_shifts.iter())
        .map(|(mc, ps)| coeff * (((2.0 * ps.l) + 1.0) * (mc.im - mc.norm_sqr())))
        .collect();
    // convergence check
    let max_l = matching::converged_values(&cs_l, absend);
    // sum
    let result = cs_l[..max_l].iter().sum();
    (max_l, result)
}

pub fn cross_section_spin_zero_no_check(
    phase_shifts: &[matching::PhaseShift],
    mel_coeff: &[Complex64],
    k: f64,
) -> f64 {
    let coeff = 10.0 * (4.0 * PI) / (k.powi(2));
    let cs_l: Vec<f64> = mel_coeff
        .iter()
        .zip(phase_shifts.iter())
        .map(|(mc, ps)| coeff * (((2.0 * ps.l) + 1.0) * (mc.im - mc.norm_sqr())))
        .collect();
    cs_l.iter().sum()
}

pub fn diff_cross_spin_zero(
    angles: &[f64],
    phase_shifts: &[matching::PhaseShift],
    mel_coeff: &[Complex64],
    k: f64,
    eta: f64,
) -> Vec<f64> {
    // Coulomb part
    let f_coul = coulomb_ampl(angles, k, eta);

    // Do the l sum
    let mut f_nuc = vec![Complex::new(0.0, 0.0); angles.len()];
    let cos_angle: Vec<f64> = angles.iter().map(|&x| f64::cos(x)).collect();
    for (&ps, &mc) in phase_shifts.iter().zip(mel_coeff.iter()) {
        let int_l = ps.l as i32;
        let coul_ps: f64 = coulomb_phase_shift(ps.l, eta); //coulomb phase shift
        let coul_term: Complex<f64> = (2.0 * coul_ps * Complex::i()).exp();

        for (i, &x) in cos_angle.iter().enumerate() {
            f_nuc[i] += mc * coul_term * 1.0 / k * (2.0 * ps.l + 1.0) * legendre_Pl(int_l, x);
        }
    }

    // Now we derive the total amplitudes and cross section
    // Sum the coulomb and nuclear terms, take the squared norm, and multiply by 10
    // to get mb/sr.
    f_nuc
        .iter()
        .zip(f_coul.iter())
        .map(|(&x, &y)| (x + y).norm_sqr() * 10.0)
        .collect()
}

/// Now we have l - 1/2 and l + 1/2.
pub fn _spin_half_ampl(
    angles: &[f64],
    phase_shift_minus: Complex<f64>,
    phase_shift_plus: Complex<f64>,
    l: f64,
    k: f64,
    eta: f64,
) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    let coul_ps: f64 = coulomb_phase_shift(l, eta); //coulomb phase shift
    let coul_term: Complex<f64> = (2.0 * coul_ps * Complex::i()).exp();
    let int_l = l as i32;
    // regular polynomial terms
    let pl: Vec<f64> = angles
        .iter()
        .map(|x| 1.0 / k * legendre_Pl(int_l, f64::cos(*x)))
        .collect();

    // associated polynomial terms, exception for l = 0
    let pl_1: Vec<f64> = if l > 0.0 {
        angles
            .iter()
            .map(|x| -1.0 / k * legendre_Plm(int_l, 1, f64::cos(*x))) // TO DO: Check the phase convention
            .collect()
    } else {
        vec![0.0; angles.len()]
    };

    // We have two coeff this time
    let c_minus: Complex<f64> =
        (-Complex::i() / 2.0) * ((2.0_f64 * Complex::i() * phase_shift_minus).exp() - 1.0);
    let c_plus: Complex<f64> =
        (-Complex::i() / 2.0) * ((2.0_f64 * Complex::i() * phase_shift_plus).exp() - 1.0);

    let a_theta: Vec<Complex<f64>> = pl
        .iter()
        .map(|x| *x * coul_term * ((l + 1.0) * c_plus + l * c_minus))
        .collect();

    let b_theta: Vec<Complex<f64>> = pl_1
        .iter()
        .map(|x| *x * Complex::i() * coul_term * (c_plus - c_minus))
        .collect();

    (a_theta, b_theta)
}

pub fn _diff_cross_section(angles: &[f64], f_nuc: &[Complex<f64>], k: f64, eta: f64) -> Vec<f64> {
    // Coulomb part
    let f_coul = coulomb_ampl(angles, k, eta);
    // total is nuclear + coulomb
    let f: Vec<Complex<f64>> = f_nuc
        .iter()
        .zip(f_coul.iter())
        .map(|(&x, &y)| x + y)
        .collect();
    // modulus and factor of 10 gets you to mb/sr
    let cs = f.iter().map(|&x| 10.0 * x.norm_sqr()).collect();
    cs
}

/// This calculates differential cross section and first order analyzing power
pub fn _all_observables(
    angles: &[f64],
    a_nuc: &[Complex<f64>],
    b_nuc: &[Complex<f64>],
    k: f64,
    eta: f64,
) -> (Vec<f64>, Vec<f64>) {
    // Coulomb part
    let f_coul = coulomb_ampl(angles, k, eta);
    // A(theta) has the nuclear + coulomb part
    let a_nuc_coul: Vec<Complex<f64>> = a_nuc
        .iter()
        .zip(f_coul.iter())
        .map(|(&x, &y)| x + y)
        .collect();
    // |A|^2 + |B|^2 for cross section
    // modulus and factor of 10 gets you to mb/sr
    let cs: Vec<f64> = a_nuc_coul
        .iter()
        .zip(b_nuc.iter())
        .map(|(&a, &b)| 10.0 * (a.norm_sqr() + b.norm_sqr()))
        .collect();
    // (A*)B + A (B*) / |A|^2 + |B|^2, and I know exactly what I did
    // I am still trying to confirm, but I think A_y = sqrt(2) * iT11
    // Eq. 3.15, Polarization Phenomena in Physics.
    let anal_power: Vec<f64> = a_nuc_coul
        .iter()
        .zip(b_nuc.iter())
        .zip(cs.iter())
        .map(|((&a, &b), &denom)| {
            let term1 = a.conj() * b;
            let term2 = a * b.conj();
            (term1 + term2).re / (denom / 10.0)
        })
        .collect();
    (cs, anal_power)
}
