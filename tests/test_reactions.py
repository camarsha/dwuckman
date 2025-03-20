import dwuckman
import subprocess
import numpy as np
import pytest


class ECISInterface:
    def __init__(
        self, input_filename: str, output_filename="output.txt"
    ) -> None:
        self.input_filename = input_filename
        self.output_filename = output_filename

    def _omp_line(self, param1: float, param2: float, param3: float) -> str:
        return f"{param1:10.5f}{param2:10.5f}{param3:10.5f}\n"

    def _masses_spins_ecis(
        self,
        e_lab: float,
        spin: float,
        proj_mass: float,
        res_mass: float,
        z_proj: float,
        z_res: float,
    ) -> str:
        z_prod = z_proj * z_res
        return f" 0.00 1 1+{e_lab:10.5f}{spin:10.5f}{proj_mass:10.5f}{res_mass:10.5f}{z_prod:10.5f}\n"

    def _ecis_calc_parameters(self) -> str:
        return "    1   20    1    1\n"

    def _100_flags_for_ecis(self) -> str:
        # Produces the header for the ecis file
        header = (
            "Spherical optical model                                                 \n"
            + "FFFFFFFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF\n"
            + "FFFFFFFFTFFFTTTFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF\n"
        )
        return header

    def _i_dont_know_lines(self) -> str:
        return (
            "   0.00000   0.00000              1.e-10    1.e-10    1.e-30\n\n"
        )

    def _end_of_ecis_input(self) -> str:
        return "fin\n"

    def run_ecis(self):
        subprocess.run(
            f"ecis < {self.input_filename} > {self.output_filename}",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return self

    def read_ecis_output(self):
        with open(self.output_filename, "r") as f:
            read_in_line = False
            angle = []
            cs = []
            ratio_to_ruth = []
            for line in f:
                if read_in_line:
                    line = line.split()
                    angle.append(float(line[0]))
                    cs.append(float(line[1]))
                    ratio_to_ruth.append(float(line[2]))
                    # end of the list
                    if angle[-1] == 180.0:
                        break
                else:
                    if "angle   cross-section      c. s./ruther." in line:
                        read_in_line = True
            # remove the file after a read
            return (np.array(angle), np.array(cs), np.array(ratio_to_ruth))

    def write_ecis_input(
        self,
        a1: float,
        z1: float,
        a2: float,
        z2: float,
        e_lab: float,
        v: float,
        r: float,
        a: float,
        w: float,
        r_i: float,
        a_i: float,
        r_c: float,
    ):
        """Update the ecis input based on the given parameters.

        :param a1:
        :param z1:
        :param a2:
        :param z2:
        :param e_lab:
        :param v:
        :param r:
        :param a:
        :param w:
        :param r_i:
        :param a_i:
        :param r_c:
        :returns:

        """

        with open(self.input_filename, "w+") as f:
            f.write(self._100_flags_for_ecis())
            f.write(self._ecis_calc_parameters())
            f.write(self._i_dont_know_lines())
            f.write(self._masses_spins_ecis(e_lab, 0.0, a1, a2, z1, z2))
            f.write(self._omp_line(v, r, a))
            f.write(self._omp_line(w, r_i, a_i))
            f.write(self._omp_line(0.0, 1.2, 0.5))  # real surface
            f.write(self._omp_line(0.0, 1.2, 0.5))  # imaginary surface
            f.write(self._omp_line(0.0, 1.0, 0.6))  # spin orbit
            f.write(self._omp_line(0.0, 1.0, 0.6))  # imaginary spin orbit
            f.write(self._omp_line(r_c, 0.0, 0.0))
            f.write(self._omp_line(0.0, 0.0, 0.0))  # ???? Just zeros always
            f.write(self._omp_line(0.00001, 1.0, 180.0))  # angles
            f.write(self._end_of_ecis_input())
        return self


def dwuckman_to_ecis_comp(a1, z1, a2, z2, e_lab, v, r, a, w, r_i, a_i, r_c):
    ecis = ECISInterface("test_ecis.input")
    angles, cs_ecis, csr_ecis = (
        ecis.write_ecis_input(a1, z1, a2, z2, e_lab, v, r, a, w, r_i, a_i, r_c)
        .run_ecis()
        .read_ecis_output()
    )
    tot, cs_dm, csr_dm = dwuckman.spin_zero(
        a1,
        z1,
        a2,
        z2,
        e_lab,
        v,
        r,
        a,
        w,
        r_i,
        a_i,
        0.0,
        0.0,
        0.0,
        r_c,
        60,
        angles[1:],
        40.0,
        0.01,
    )
    return (
        angles[1:],
        cs_ecis[1:],
        csr_ecis[1:],
        cs_dm,
        csr_dm,
        np.abs(np.array(cs_dm) - cs_ecis[1:]) / cs_ecis[1:],
        np.abs(np.array(csr_dm) - csr_ecis[1:]) / csr_ecis[1:],
    )


def test_dwuckman_ecis_1():
    _, _, _, _, _, diff, diffr = dwuckman_to_ecis_comp(
        4.0, 2.0, 86.0, 36.0, 10.0, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 1.3
    )
    assert np.nanmax(diff) < 0.03


def test_dwuckman_ecis_2():
    _, _, _, _, _, diff, diffr = dwuckman_to_ecis_comp(
        4.0, 2.0, 86.0, 36.0, 15.0, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 1.3
    )
    assert np.nanmax(diff) < 0.03


def test_dwuckman_ecis_3():
    _, _, _, _, _, diff, diffr = dwuckman_to_ecis_comp(
        4.0, 2.0, 86.0, 36.0, 20.0, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 1.3
    )
    assert np.nanmax(diff) < 0.03


def test_dwuckman_ecis_4():
    _, _, _, _, _, diff, diffr = dwuckman_to_ecis_comp(
        4.0, 2.0, 86.0, 36.0, 40.0, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 1.3
    )
    assert np.nanmax(diff) < 0.03


def test_dwuckman_ecis_5():
    _, _, _, _, _, diff, diffr = dwuckman_to_ecis_comp(
        4.0, 2.0, 20.0, 10.0, 20.0, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 1.3
    )
    assert np.nanmax(diff) < 0.03
