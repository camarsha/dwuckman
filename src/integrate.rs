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

/// This solves Z'' = A * Z for Z and Z complex using Crowell-Fox-Goodwin method.
/// The equation is Z_i = (12.0 - 10.0 * C_{i-1})Z_{i - 1} - C_{i-2}Z_{i-2} * C^{-1}_i.
/// This is a matrix equation since both the potential and wave function are complex.
pub fn fox_goodwin_coupled(
    h: f64,
    q_r: &[f64],
    q_i: &[f64],
    phi_r: &mut [f64],
    phi_i: &mut [f64],
    start_idx: usize,
) {
    let g = h.powi(2) / 12.0;

    let end = phi_r.len();
    let start = start_idx + 1;
    for i in start..end {
        // real terms in front of the wave function
        let cr3 = 1.0 - q_r[i] * g;
        let cr2 = 1.0 - q_r[i - 1] * g;
        let cr1 = 1.0 - q_r[i - 2] * g;

        // imaginary terms in front of the wave function.
        let ci3 = -q_i[i] * g;
        let ci2 = -q_i[i - 1] * g;
        let ci1 = -q_i[i - 2] * g;

        // for the i-1 terms we need the coefficent to be
        // 2 + 5/6h^2 * A not 1 - 1/12h^2 * A, so this transformation
        // gives us the right form.
        let cr2 = 12.0 - 10.0 * cr2;
        let ci2 = -10.0 * ci2;

        // These are the real and imaginary parts of the numerator.
        // Basically just two sets of (cr + i ci ) * (phi_r + i phi_i).
        let real =
            cr2 * phi_r[i - 1] - ci2 * phi_i[i - 1] - cr1 * phi_r[i - 2] + ci1 * phi_i[i - 2];
        let im = ci2 * phi_r[i - 1] + cr2 * phi_i[i - 1] - ci1 * phi_r[i - 2] - cr1 * phi_i[i - 2];

        // to finish up we have one final matrix multiplication between the two terms above
        // and the A^{-1}_i term on bottom.
        let det = cr3.powi(2) + ci3.powi(2);

        // (cr3  ci3)   (Rl)   1/
        // (-ci3 cr3) * (Im) * det
        phi_r[i] = (real * cr3 + im * ci3) / det;
        phi_i[i] = (im * cr3 - real * ci3) / det;
    }
}
