"""End-to-end regression tests for the public matcher API."""

import math

from miqrophi import EpitaxyMatcher, Lattice2D, MatcherConfig, SUBSTRATE_DB


class TestEpitaxyMatcherRegression:
    def test_hkust1_on_au100_best_match_is_stable(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        matcher = EpitaxyMatcher(MatcherConfig(sigma=0.4, eta_tol=0.04))

        df = matcher.run(SUBSTRATE_DB["Au_100"], hkust1, verbose=False)

        assert df is not None
        assert not df.empty
        row = df.iloc[0]
        assert math.isclose(row["theta (deg)"], 50.76, abs_tol=0.2)
        assert math.isclose(row["eta"], 0.002141, rel_tol=1e-3)
        assert math.isclose(row["area (A2)"], 5876.2, rel_tol=1e-3)
        assert row["M"] == [[1, 4], [-4, 1]]
        assert row["N"] == [[17, 8], [-8, 17]]
        assert bool(row["L0_feasible"]) is True

    def test_hkust1_on_au111_returns_approximate_matches_only(self):
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        matcher = EpitaxyMatcher(MatcherConfig(sigma=0.4, eta_tol=0.04))

        df = matcher.run(SUBSTRATE_DB["Au_111"], hkust1, verbose=False)

        assert df is not None
        assert not df.empty
        assert set(df["L0_feasible"].unique()) == {False}
        row = df.iloc[0]
        assert math.isclose(row["theta (deg)"], 165.06, abs_tol=0.2)
        assert math.isclose(row["eta"], 0.001781, rel_tol=1e-3)
        assert math.isclose(row["area (A2)"], 7987.7, rel_tol=1e-3)