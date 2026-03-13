"""Regression tests for Level 1 coincidence-function peaks."""

import math

from miqrophi import Lattice2D, SUBSTRATE_DB, coincidence


class TestLevel1Regression:
    def test_square_on_square_peak_family(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        result = coincidence.compute(SUBSTRATE_DB["Au_100"], hkust1, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 5
        # stable lexsort: within equal-phi peaks, ascending theta order
        # -> 90° (smaller) must appear before 270° in the output list.
        # 180° may also appear as a symmetry-equivalent max-phi peak depending on
        # platform float precision; we check relative ordering rather than fixed indices.
        assert math.isclose(result.phi_peaks[0], 1.0, rel_tol=1e-9)
        idx_90  = next(i for i, t in enumerate(result.theta_peaks)
                       if math.isclose(t, 90.0, abs_tol=0.2))
        idx_270 = next(i for i, t in enumerate(result.theta_peaks)
                       if math.isclose(t, 270.0, abs_tol=0.2))
        assert idx_90 < idx_270, "90° must appear before 270° (ascending-theta tiebreak)"
        # Secondary family: 50.76° and 140.76° are 4-fold-symmetry equivalents.
        _EQUIV_FAMILY = [50.76, 140.76]
        sec_idx = next(
            i for i, t in enumerate(result.theta_peaks)
            if any(math.isclose(t, ex, abs_tol=0.2) for ex in _EQUIV_FAMILY)
        )
        assert result.phi_peaks[sec_idx] > 0.94

    def test_hexagonal_domains_are_nearly_degenerate(self):
        hex_mof = Lattice2D(4.92, 4.92, 120.0, "Hex-MOF(2xa)")
        result = coincidence.compute(SUBSTRATE_DB["Au_111"], hex_mof, G_cutoff=8.0, sigma=0.4)

        assert len(result.theta_peaks) >= 6
        # stable lexsort: within equal-phi peaks, ascending theta order.
        # The three peak pairs are nearly degenerate; their inter-pair order may
        # vary across platforms.  We verify only the intra-pair property:
        # within each pair the smaller theta always comes before the larger one.
        expected_pairs = [(90.0, 270.0), (149.94, 329.94), (30.06, 210.06)]
        for small, large in expected_pairs:
            idx_small = next(i for i, t in enumerate(result.theta_peaks)
                             if math.isclose(t, small, abs_tol=0.2))
            idx_large = next(i for i, t in enumerate(result.theta_peaks)
                             if math.isclose(t, large, abs_tol=0.2))
            assert idx_small < idx_large, (
                f"{small}° must appear before {large}° (ascending-theta tiebreak)"
            )
        assert min(result.phi_peaks[:6]) > 0.9999
