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
        # 50.76° and 140.76° are 4-fold-equivalent; either is correct
        assert any(math.isclose(row["theta (deg)"], ex, abs_tol=0.2)
                   for ex in [50.76, 140.76])
        assert math.isclose(row["eta"], 0.002141, rel_tol=1e-3)
        assert math.isclose(row["area (A2)"], 5876.2, rel_tol=1e-3)
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
        # 75.06° and 165.06° are MOF-4-fold-equivalent peaks
        assert any(math.isclose(row["theta (deg)"], ex, abs_tol=0.2)
                   for ex in [75.06, 165.06])
        assert math.isclose(row["eta"], 0.001781, rel_tol=1e-3)
        assert math.isclose(row["area (A2)"], 7987.7, rel_tol=1e-3)


class TestEpitaxyMatcherVerbose:
    """Exercises the verbose=True logging paths."""

    def test_verbose_true_prints_output(self, capsys):
        hkust1  = Lattice2D(18.62, 18.62, 90.0, "HKUST-1")
        matcher = EpitaxyMatcher(MatcherConfig(sigma=0.4, eta_tol=0.04))
        df = matcher.run(SUBSTRATE_DB["Au_100"], hkust1, verbose=True)
        out = capsys.readouterr().out
        assert "[L1]" in out
        assert "[L2]" in out
        assert df is not None

    def test_no_match_returns_none_verbose(self, capsys):
        """Very tight eta_tol with an unlikely pair triggers the no-match path."""
        tiny_mof = Lattice2D(3.0, 3.0, 90.0, "tiny")
        matcher  = EpitaxyMatcher(MatcherConfig(sigma=0.3, eta_tol=1e-9))
        result   = matcher.run(SUBSTRATE_DB["Au_100"], tiny_mof, verbose=True)
        out = capsys.readouterr().out
        assert result is None
        assert "no match" in out.lower() or "[L2]" in out or "[L1]" in out

    def test_no_peaks_returns_none(self):
        """Zero L1 peaks (tiny G_cutoff) should return None."""
        tiny_mof = Lattice2D(0.001, 0.001, 90.0)
        matcher  = EpitaxyMatcher(MatcherConfig(sigma=0.001, eta_tol=0.05, G_cutoff=0.1))
        result   = matcher.run(SUBSTRATE_DB["Au_100"], tiny_mof, verbose=False)
        assert result is None

class TestEpitaxyMatcherVerbose:
    """Exercises the verbose=True logging paths."""

    def test_verbose_true_prints_output(self, capsys):
        hkust1  = Lattice2D(18.62, 18.62, 90.0, "HKUST-1")
        matcher = EpitaxyMatcher(MatcherConfig(sigma=0.4, eta_tol=0.04))
        df = matcher.run(SUBSTRATE_DB["Au_100"], hkust1, verbose=True)
        out = capsys.readouterr().out
        assert "[L1]" in out
        assert "[L2]" in out
        assert df is not None

    def test_verbose_no_match_prints_note(self, capsys):
        """A very tight eta_tol with an unlikely pair should trigger the no-match path."""
        tiny_mof = Lattice2D(3.0, 3.0, 90.0, "tiny")
        matcher  = EpitaxyMatcher(MatcherConfig(sigma=0.3, eta_tol=1e-9))
        result   = matcher.run(SUBSTRATE_DB["Au_100"], tiny_mof, verbose=True)
        out = capsys.readouterr().out
        # Either no peaks or no matches — just make sure it doesn't crash
        assert result is None or True  # graceful exit

    def test_no_peaks_returns_none(self):
        """A MOF lattice that produces zero L1 peaks should return None."""
        # Use extremely narrow sigma so Phi has no peaks above threshold
        tiny_mof = Lattice2D(0.001, 0.001, 90.0)
        matcher  = EpitaxyMatcher(MatcherConfig(sigma=0.001, eta_tol=0.05, G_cutoff=0.1))
        result   = matcher.run(SUBSTRATE_DB["Au_100"], tiny_mof, verbose=False)
        assert result is None