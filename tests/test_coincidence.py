"""Regression tests for Level 1 coincidence-function peaks."""

import math

from miqrophi import Lattice2D, SUBSTRATE_DB, coincidence


class TestLevel1Regression:
    def test_square_on_square_peak_family(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        result = coincidence.compute(SUBSTRATE_DB["Au_100"], hkust1, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 5
        # stable lexsort: within equal-phi peaks, ascending theta order
        # -> 90° (smaller) comes before 270°
        assert math.isclose(result.theta_peaks[0], 90.0, abs_tol=0.2)
        assert math.isclose(result.theta_peaks[1], 270.0, abs_tol=0.2)
        # 50.76° and 140.76° are 4-fold-symmetry-equivalent peaks; which one
        # appears at index 3 depends on floating-point tie-breaking in the sort.
        _EQUIV_FAMILY = [50.76, 140.76]
        assert any(math.isclose(result.theta_peaks[3], ex, abs_tol=0.2)
                   for ex in _EQUIV_FAMILY)
        assert math.isclose(result.phi_peaks[0], 1.0, rel_tol=1e-9)
        assert result.phi_peaks[3] > 0.94

    def test_hexagonal_domains_are_nearly_degenerate(self):
        hex_mof = Lattice2D(4.92, 4.92, 120.0, "Hex-MOF(2xa)")
        result = coincidence.compute(SUBSTRATE_DB["Au_111"], hex_mof, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 6
        # stable lexsort: within equal-phi peaks, ascending theta order
        # -> 90° < 270°; 149.94° < 329.94°; 30.06° < 210.06°
        expected = [90.0, 270.0, 149.94, 329.94, 30.06, 210.06]
        for got, want in zip(result.theta_peaks[:6], expected):
            assert math.isclose(got, want, abs_tol=0.2)
        assert min(result.phi_peaks[:6]) > 0.9999
