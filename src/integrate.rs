pub fn fox_goodwin(h: f64, q: &[f64], phi: &mut [f64], start_idx: usize) {
    // fox goodwin algorithm for a real wave function and potential. Used for testing.

    let g = h.powi(2) / 12.0;
    // three step algorithm
    let mut y1 = 0.0_f64;
    let mut y2 = 0.0_f64;
    let mut y3 = 0.0_f64;
    let end = phi.len();
    let start = start_idx + 1;

    for i in start..end {
        y1 = (1.0 - (g * q[i]));
        y2 = (2.0 + (10.0 * g * q[i - 1]));
        y3 = (1.0 - (g * q[i - 2]));
        phi[i] = 1.0 / y1 * (y2 * phi[i - 1] - y3 * phi[i - 2]);
        // check if renormalization is needed, taken from ECIS
        if f64::abs(phi[i]) > 1e15 {
            //and if so do it
            for j in 0..i {
                phi[i] = phi[i] * 1e-30
            }
        };
    }
}

pub fn fox_goodwin_coupled(
    h: f64,
    q_r: &[f64],
    q_i: &[f64],
    phi_r: &mut [f64],
    phi_i: &mut [f64],
    start_idx: usize,
) {
    // Just plug in complex potentials and wave functions and follow
    // mathematica blindly - Caleb Marshall, scholar

    let g = h.powi(2) / 12.0;

    let (mut ar1, mut ar2, mut ar3) = (0.0, 0.0, 0.0);
    let (mut ai1, mut ai2, mut ai3) = (0.0, 0.0, 0.0);
    let (mut br, mut bi) = (0.0, 0.0);
    let (mut cr, mut ci) = (0.0, 0.0);
    let mut det = 0.0;

    let end = phi_r.len();
    let start = start_idx + 1;

    for i in start..end {
        ar3 = 1.0 - q_r[i] * g;
        ar2 = 1.0 - q_r[i - 1] * g;
        ar1 = 1.0 - q_r[i - 2] * g;

        ai3 = -q_i[i] * g;
        ai2 = -q_i[i - 1] * g;
        ai1 = -q_i[i - 2] * g;

        det = ar3.powi(2) + ai3.powi(2);

        br = 12.0 - 10.0 * ar2;
        bi = -10.0 * ai2;

        cr = br * phi_r[i - 1] - bi * phi_i[i - 1] - ar1 * phi_r[i - 2] + ai1 * phi_i[i - 2];
        ci = bi * phi_r[i - 1] + br * phi_i[i - 1] - ai1 * phi_r[i - 2] - ar1 * phi_i[i - 2];

        phi_r[i] = (cr * ar3 + ci * ai3) / det;
        phi_i[i] = (ci * ar3 - cr * ai3) / det;

        // check if renormalization is needed
        if f64::abs(phi_r[i]) > 1e15 {
            //and if so do it
            for j in 0..i {
                phi_r[i] = phi_r[i] * 1e-30
            }
        }
    }
}
