from .dwuckman import spin_zero, wave_function_spin_zero, phase_shift_spin_zero
import numpy as np


class SpinZero:
    def __init__(self, m1, z1, m2, z2, n_partial_waves, r_match, dr) -> None:
        self.m1 = m1
        self.m2 = m2
        self.z1 = z1
        self.z2 = z2
        self.n_partial_waves = n_partial_waves
        self.r_match = r_match
        self.dr = dr

    def wave_function(self, energy_lab, ell, pot_params):
        return wave_function_spin_zero(
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

    def phase_shifts(self, energy_lab, pot_params):
        return phase_shift_spin_zero(
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
        )
        return tot, np.asarray(diff), np.asarray(ruth)
