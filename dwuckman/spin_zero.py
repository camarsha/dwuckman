from .dwuckman import spin_zero, wave_function_spin_zero, phase_shift_spin_zero
import numpy as np


class SpinZero:
    def __init__(
        self, m1, z1, m2, z2, n_partial_waves, r_match, dr, absend=1e-3
    ) -> None:
        self.m1 = m1
        self.m2 = m2
        self.z1 = z1
        self.z2 = z2
        self.n_partial_waves = n_partial_waves
        self.r_match = r_match
        self.dr = dr
        self.absend = absend

    def wave_function(self, energy_lab, ell, pot_params):
        r, re, im = wave_function_spin_zero(
            self.m1,
            self.z1,
            self.m2,
            self.z2,
            energy_lab,
            pot_params,
            ell,
            self.n_partial_waves,
            [10.0],  # just a dummy angle
            self.r_match,
            self.dr,
        )
        return np.asarray(r), np.asarray(re), np.asarray(im)

    def phase_shifts(self, energy_lab, pot_params):
        l, re, im = phase_shift_spin_zero(
            self.m1,
            self.z1,
            self.m2,
            self.z2,
            energy_lab,
            pot_params,
            self.n_partial_waves,
            [10.0],  # just a dummy angle
            self.r_match,
            self.dr,
        )
        return np.asarray(l), np.asarray(re), np.asarray(im)

    def s_matrix(self, energy_lab, pot_params):
        l, re, im = self.phase_shifts(energy_lab, pot_params)
        s = np.exp(2.0j * (re + (1j * im)))
        return (np.asarray(l), s.real, s.imag)

    def cross_section(self, energy_lab, angles, pot_params):
        tot, diff, ruth = spin_zero(
            self.m1,
            self.z1,
            self.m2,
            self.z2,
            energy_lab,
            pot_params,
            self.n_partial_waves,
            angles,
            self.r_match,
            self.dr,
            self.absend,
        )
        return tot, np.asarray(diff), np.asarray(ruth)
