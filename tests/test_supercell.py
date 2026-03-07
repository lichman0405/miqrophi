"""Regression tests for Level 2 supercell matching."""

import math

from miqrophi import Lattice2D, SUBSTRATE_DB, coincidence, supercell


class TestLevel2Regression:
    def test_commensurate_square_match_at_rotated_peak(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        l1 = coincidence.compute(SUBSTRATE_DB["Au_100"], hkust1, G_cutoff=8.0, sigma=0.4)

        matches = supercell.find_matches(
            SUBSTRATE_DB["Au_100"],
            hkust1,
            theta_deg=l1.theta_peaks[3],
            eta_tol=0.04,
        )

        assert len(matches) >= 1
        best = matches[0]
        # theta_peaks[3] may be either symmetry-equivalent peak (50.76° or 140.76°)
        _EQUIV = [50.76, 140.76]
        assert any(math.isclose(best.theta_deg, ex, abs_tol=0.2) for ex in _EQUIV)
        assert math.isclose(best.eta, 0.002141, rel_tol=1e-3)
        assert math.isclose(best.area, 5876.2, rel_tol=1e-3)
        # N (substrate supercell) is orientation-independent
        assert best.N.tolist() == [[17, 8], [-8, 17]]

    def test_axis_aligned_peak_gives_small_supercell_solution(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        matches = supercell.find_matches(
            SUBSTRATE_DB["Au_100"],
            hkust1,
            theta_deg=270.0,
            eta_tol=0.04,
        )

        assert len(matches) == 1
        best = matches[0]
        assert math.isclose(best.eta, 0.020169, rel_tol=1e-3)
        assert math.isclose(best.area, 1348.4, rel_tol=1e-3)
        assert best.M.tolist() == [[0, 2], [-2, 0]]
        assert best.N.tolist() == [[-9, 0], [0, -9]]
