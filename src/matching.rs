use num::complex::{Complex, Complex64};
use rgsl::coulomb::wave_FG_e;
use std::cmp::Ordering;

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
#[allow(non_snake_case, unused_mut)]
pub fn coulomb_functions(rho: f64, eta: f64, l: f64) -> Vec<f64> {
    let mut exp_F = 0.0_f64;
    let mut exp_G = 0.0_f64;

    // Get the coulomb functions at rho
    let (mut F, mut Fp, mut G, mut Gp) = wave_FG_e(eta, rho, l, 0, &mut exp_F, &mut exp_G)
        .unwrap_or_else(|_| {
            panic!(
                "Overflow in coulomb wavefunctions: rho={:.3} eta={:.3} l={}",
                rho, eta, l as i32
            );
        });

    vec![F.val, Fp.val, G.val, Gp.val]
}

#[allow(non_snake_case)]
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

#[inline(always)]
pub fn s_matrix(phase_shift: Complex64) -> Complex64 {
    (2.0_f64 * Complex::i() * phase_shift).exp()
}

/// Given a vector of PhaseShift structs that are already ordered according to
/// l_i > l_{i - 1}, return a new vector that only has the strictly increasing
/// real S-matrix values.
pub fn converged_values(phase_shifts: &[PhaseShift]) -> Vec<PhaseShift> {
    let mut stop_l: usize = phase_shifts.len();
    let mut begin_check = false;
    for (i, &ele) in phase_shifts.iter().enumerate() {
        let re = s_matrix(ele.val).re;
        if begin_check && (re < s_matrix(phase_shifts[i - 1].val).re) || (re >= 1.0) {
            stop_l = i; // stopping index is exclusive
            break;
        }

        // We first wait until the phase shift crosses a threshold value of 0.99
        // to start checking for convergence
        if re > 0.99 {
            begin_check = true
        };
    }
    if !begin_check {
        // This means that the values never crossed 0.99, raise an error.
        panic!(
            "Non-convergence in phase shifts! Last value: {:.5}",
            s_matrix(phase_shifts.last().unwrap().val).re
        );
    }
    phase_shifts[..stop_l].to_vec()
}
