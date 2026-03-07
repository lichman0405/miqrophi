"""Regression tests for Level 1 coincidence-function peaks."""

import math

from miqrophi import Lattice2D, SUBSTRATE_DB, coincidence


class TestLevel1Regression:
    def test_square_on_square_peak_family(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        result = coincidence.compute(SUBSTRATE_DB["Au_100"], hkust1, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 5
        assert math.isclose(result.theta_peaks[0], 270.0, abs_tol=0.2)
        assert math.isclose(result.theta_peaks[1], 90.0, abs_tol=0.2)
        assert math.isclose(result.theta_peaks[3], 50.76, abs_tol=0.2)
        assert math.isclose(result.phi_peaks[0], 1.0, rel_tol=1e-9)
        assert result.phi_peaks[3] > 0.94

    def test_hexagonal_domains_are_nearly_degenerate(self):
        hex_mof = Lattice2D(4.92, 4.92, 120.0, "Hex-MOF(2xa)")
        result = coincidence.compute(SUBSTRATE_DB["Au_111"], hex_mof, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 6
        expected = [90.0, 270.0, 329.94, 149.94, 30.06, 210.06]
        for got, want in zip(result.theta_peaks[:6], expected):
            assert math.isclose(got, want, abs_tol=0.2)
        assert min(result.phi_peaks[:6]) > 0.9999
