use num::complex::{Complex, Complex64};
use rgsl::bessel::{jl, yl};
use rgsl::coulomb::wave_FG_e;
use rgsl::{Result, Value};
use std::cmp::Ordering;
use std::f64::consts::E;

/* We define all of the phase shift structs here. There is one for each supported spin value
Each of them implements PartialEq and PartialOrd so that we can sort after the parallel partial wave loop.
*/

// we need to be able to sort values based on l, s, and j so we need to derive these traits
#[derive(Debug, Clone, Copy)]
pub struct PhaseShift {
    pub val: Complex64,
    pub l: f64,
}

impl PhaseShift {
    pub fn new(val: Complex64, l: f64) -> Self {
        PhaseShift { val, l }
    }
}

impl PartialEq for PhaseShift {
    fn eq(&self, other: &Self) -> bool {
        self.l == other.l
    }
}

impl PartialOrd for PhaseShift {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.l.partial_cmp(&other.l)
    }
}

/*
Now we have the functions that actual calculate the phase shifts
*/

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
        //        println! {"Overflow in coulomb wave functions!"}
    }

    vec![F.val, Fp.val, G.val, Gp.val]
}

pub fn phase_shift(
    phi_R: Complex64,
    phi_Rh: Complex64,
    rho_R: f64,
    rho_Rh: f64,
    eta: f64,
    l: f64,
) -> PhaseShift {
    /* Calculate the nuclear phase shift using the psuedo-Wronskian
    suggested by Carl.
     */

    // asymptotic functions at R
    let fun_R: Vec<f64> = coulomb_functions(rho_R, eta, l);
    // asymptotic functions at Rh
    let fun_Rh: Vec<f64> = coulomb_functions(rho_Rh, eta, l);

    // Do the matching
    let num = (phi_R * fun_Rh[0]) - (phi_Rh * fun_R[0]);
    let denom = (phi_R * fun_Rh[2]) - (phi_Rh * fun_R[2]);

    PhaseShift::new(-1.0 * (num / denom).atan(), l)
}
