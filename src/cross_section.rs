use num::complex::Complex;
use rgsl::gamma_beta::gamma::lngamma_complex_e;
use rgsl::legendre::polynomials::legendre_Pl;
use rgsl::{Result, Value};

pub fn coulomb_phase_shift(l: f64, eta: f64) -> f64 {
    let zr = l + 1.0;
    let zi = eta;
    let mut overflow: Value = Value::OverFlow;

    let (overflow, mut lnr, mut arg) = lngamma_complex_e(zr, zi);
    // We just care about arg
    arg.val
}

pub fn coulomb_ampl(angles: &[f64], eta: f64, k: f64) -> Vec<Complex<f64>> {
    // This is the point-Coulomb scattering amplitude, i.e it will give you the Rutherford cross section.
    let sig_0 = coulomb_phase_shift(0.0, eta); // l = 0 phase shift
    let mut f = vec![Complex::new(0.0, 0.0); angles.len()];

    for (i, &a) in angles.iter().enumerate() {
        let denom = Complex::new(2.0 * k * (a / 2.0).sin().powi(2), 0.0);
        let num = Complex::new(
            0.0,
            -1.0 * ((a / 2.0).sin().powi(2).ln()) * eta + (2.0 * sig_0),
        );
        f[i] = (-eta / denom) * num.exp();
    }

    f
}

pub fn rutherford_cs(angles: &[f64], eta: f64, k: f64) -> Vec<f64> {
    // Calculate the Rutherford differential cross section in mb/sr
    let f: Vec<Complex<f64>> = coulomb_ampl(angles, eta, k);
    let dR: Vec<f64> = f.iter().map(|x| x.norm_sqr()).collect();
    dR
}

// calculate spin 0 amplitude

pub fn spin_zero_amp(
    angles: &[f64],
    phase_shift: Complex<f64>,
    l: f64,
    eta: f64,
    k: f64,
) -> Vec<Complex<f64>> {
    let coul_ps: f64 = coulomb_phase_shift(l, eta); //coulomb phase shift
    let coul_term: Complex<f64> = (2.0 * coul_ps * Complex::i()).exp();
    let int_l = l as i32;
    let Pl: Vec<f64> = angles
        .iter()
        .map(|x| 1.0 / k * (2.0 * l + 1.0) * legendre_Pl(int_l, f64::cos(*x)))
        .collect();
    let C: Complex<f64> =
        (-Complex::i() / 2.0) * ((2.0_f64 * Complex::i() * phase_shift).exp() - 1.0);
    Pl.iter().map(|x| *x * coul_term * C).collect()
}

pub fn diff_cross_section(angles: &[f64], eta: f64, k: f64, f_nuc: &[Complex<f64>]) -> Vec<f64> {
    // Coulomb part
    let f_coul = coulomb_ampl(&angles, eta, k);
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
