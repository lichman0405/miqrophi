"""Tests for miqrophi.visualize — all plotting and PDF-generation functions."""
import math
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # non-interactive backend; must be set before any plt import
import matplotlib.pyplot as plt
import numpy as np
import pytest

from miqrophi import Lattice2D, SUBSTRATE_DB, coincidence, supercell
from miqrophi.visualize import (
    _rot,
    _real_points_in_box,
    _recip_points_in_disk,
    plot_phi_curve,
    plot_lattice_overlay,
    plot_leed_pattern,
    plot_strain_ellipse,
    plot_match_card,
    generate_pdf_report,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def _square_pair():
    """Au(100) + HKUST-1 as a well-tested commensurate pair."""
    lat_sub = SUBSTRATE_DB["Au_100"]
    lat_mof = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
    return lat_sub, lat_mof


@pytest.fixture(scope="module")
def _l1(_square_pair):
    lat_sub, lat_mof = _square_pair
    return coincidence.compute(lat_sub, lat_mof, G_cutoff=8.0, sigma=0.4)


@pytest.fixture(scope="module")
def _match(_square_pair, _l1):
    lat_sub, lat_mof = _square_pair
    matches = supercell.find_matches(
        lat_sub, lat_mof,
        theta_deg=_l1.theta_peaks[0],
        eta_tol=0.10,
    )
    assert matches, "Expected at least one match for fixture setup"
    return matches[0]


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

class TestRotMatrix:
    def test_identity_at_zero(self):
        R = _rot(0.0)
        assert np.allclose(R, np.eye(2), atol=1e-12)

    def test_90_degrees(self):
        R = _rot(90.0)
        v = R @ np.array([1.0, 0.0])
        assert np.allclose(v, [0.0, 1.0], atol=1e-12)

    def test_180_degrees(self):
        R = _rot(180.0)
        v = R @ np.array([1.0, 0.0])
        assert np.allclose(v, [-1.0, 0.0], atol=1e-10)

    def test_returns_2x2(self):
        assert _rot(45.0).shape == (2, 2)


class TestRealPointsInBox:
    def test_square_lattice_has_gamma_point(self):
        A = np.eye(2) * 3.0
        pts = _real_points_in_box(A, half_box=1.5)
        assert any(np.allclose(p, [0, 0]) for p in pts)

    def test_count_scales_with_box(self):
        A = np.eye(2) * 1.0
        pts_small = _real_points_in_box(A, half_box=1.5)
        pts_large = _real_points_in_box(A, half_box=3.5)
        assert len(pts_large) > len(pts_small)

    def test_all_points_within_box(self):
        A = np.eye(2) * 2.5
        half = 5.0
        pts = _real_points_in_box(A, half_box=half)
        assert np.all(np.abs(pts) <= half + 1e-10)

    def test_empty_lattice_returns_empty(self):
        # very tiny box, only origin should fit
        A = np.eye(2) * 100.0
        pts = _real_points_in_box(A, half_box=0.01)
        assert len(pts) == 1  # only the origin


class TestRecipPointsInDisk:
    def test_excludes_gamma(self):
        B = np.eye(2) * (2 * np.pi)
        pts = _recip_points_in_disk(B, G_max=20.0)
        assert not any(np.allclose(p, [0, 0]) for p in pts)

    def test_all_within_radius(self):
        B = np.eye(2) * (2 * np.pi / 4.0)
        G_max = 4.0
        pts = _recip_points_in_disk(B, G_max=G_max)
        norms = np.linalg.norm(pts, axis=1)
        assert np.all(norms < G_max + 1e-10)

    def test_larger_radius_more_points(self):
        B = np.eye(2) * (2 * np.pi / 4.0)
        pts_small = _recip_points_in_disk(B, G_max=2.0)
        pts_large = _recip_points_in_disk(B, G_max=8.0)
        assert len(pts_large) > len(pts_small)


# ---------------------------------------------------------------------------
# plot_phi_curve
# ---------------------------------------------------------------------------

class TestPlotPhiCurve:
    def test_returns_axes_no_ax(self, _l1):
        from matplotlib.axes import Axes
        ax = plot_phi_curve(_l1)
        assert isinstance(ax, Axes)
        plt.close("all")

    def test_returns_axes_with_ax(self, _l1):
        from matplotlib.axes import Axes
        fig, ax_in = plt.subplots()
        ax_out = plot_phi_curve(_l1, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_custom_title(self, _l1):
        ax = plot_phi_curve(_l1, title="My title")
        assert ax.get_title() == "My title"
        plt.close("all")

    def test_xlim_matches_theta_grid(self, _l1):
        ax = plot_phi_curve(_l1)
        xlim = ax.get_xlim()
        assert math.isclose(xlim[0], _l1.theta_grid[0], abs_tol=1.0)
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_lattice_overlay
# ---------------------------------------------------------------------------

class TestPlotLatticeOverlay:
    def test_returns_axes(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_lattice_overlay(lat_sub, lat_mof, _match)
        from matplotlib.axes import Axes
        assert isinstance(ax, Axes)
        plt.close("all")

    def test_accepts_existing_ax(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        fig, ax_in = plt.subplots()
        ax_out = plot_lattice_overlay(lat_sub, lat_mof, _match, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_custom_title(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_lattice_overlay(lat_sub, lat_mof, _match, title="Test overlay")
        assert "Test overlay" in ax.get_title()
        plt.close("all")

    def test_aspect_equal(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_lattice_overlay(lat_sub, lat_mof, _match)
        # get_aspect() may return 'equal' or 1.0 depending on matplotlib version
        assert ax.get_aspect() in ("equal", 1.0)
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_leed_pattern
# ---------------------------------------------------------------------------

class TestPlotLeedPattern:
    def test_returns_axes(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_leed_pattern(lat_sub, lat_mof, _match)
        from matplotlib.axes import Axes
        assert isinstance(ax, Axes)
        plt.close("all")

    def test_accepts_existing_ax(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        fig, ax_in = plt.subplots()
        ax_out = plot_leed_pattern(lat_sub, lat_mof, _match, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_custom_g_max(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_leed_pattern(lat_sub, lat_mof, _match, G_max=2.0)
        assert isinstance(ax, object)
        plt.close("all")

    def test_custom_title(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_leed_pattern(lat_sub, lat_mof, _match, title="LEED test")
        assert "LEED test" in ax.get_title()
        plt.close("all")

    def test_aspect_equal(self, _square_pair, _match):
        lat_sub, lat_mof = _square_pair
        ax = plot_leed_pattern(lat_sub, lat_mof, _match)
        assert ax.get_aspect() in ("equal", 1.0)
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_strain_ellipse
# ---------------------------------------------------------------------------

class TestPlotStrainEllipse:
    def test_returns_axes(self, _match):
        from matplotlib.axes import Axes
        ax = plot_strain_ellipse(_match)
        assert isinstance(ax, Axes)
        plt.close("all")

    def test_accepts_existing_ax(self, _match):
        fig, ax_in = plt.subplots()
        ax_out = plot_strain_ellipse(_match, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_custom_title(self, _match):
        ax = plot_strain_ellipse(_match, title="Strain test")
        assert "Strain test" in ax.get_title()
        plt.close("all")

    def test_aspect_equal(self, _match):
        ax = plot_strain_ellipse(_match)
        assert ax.get_aspect() in ("equal", 1.0)
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_match_card
# ---------------------------------------------------------------------------

class TestPlotMatchCard:
    def test_returns_figure(self, _square_pair, _l1, _match):
        from matplotlib.figure import Figure
        lat_sub, lat_mof = _square_pair
        fig = plot_match_card(lat_sub, lat_mof, _l1, _match)
        assert isinstance(fig, Figure)
        plt.close("all")

    def test_saves_png(self, _square_pair, _l1, _match, tmp_path):
        lat_sub, lat_mof = _square_pair
        out = tmp_path / "card.png"
        fig = plot_match_card(lat_sub, lat_mof, _l1, _match,
                              save_path=str(out), dpi=72)
        plt.close("all")
        assert out.exists()
        assert out.stat().st_size > 1000

    def test_has_four_axes(self, _square_pair, _l1, _match):
        lat_sub, lat_mof = _square_pair
        fig = plot_match_card(lat_sub, lat_mof, _l1, _match)
        assert len(fig.axes) == 4
        plt.close("all")


# ---------------------------------------------------------------------------
# generate_pdf_report
# ---------------------------------------------------------------------------

class TestGeneratePdfReport:
    def test_creates_pdf(self, _square_pair, _l1, _match, tmp_path):
        lat_sub, lat_mof = _square_pair
        matches = supercell.find_matches(
            lat_sub, lat_mof,
            theta_deg=_l1.theta_peaks[0],
            eta_tol=0.10,
        )
        out = tmp_path / "report.pdf"
        path = generate_pdf_report(
            lat_sub, lat_mof, _l1, matches,
            title="Test Report",
            save_path=str(out),
            dpi=72,
        )
        assert Path(path).exists()
        assert Path(path).stat().st_size > 5000

    def test_returns_absolute_path(self, _square_pair, _l1, _match, tmp_path):
        lat_sub, lat_mof = _square_pair
        matches = [_match]
        out = tmp_path / "sub" / "r.pdf"
        path = generate_pdf_report(
            lat_sub, lat_mof, _l1, matches,
            save_path=str(out), dpi=72,
        )
        assert Path(path).is_absolute()

    def test_empty_matches_raises(self, _square_pair, _l1, tmp_path):
        lat_sub, lat_mof = _square_pair
        with pytest.raises(ValueError, match="empty"):
            generate_pdf_report(
                lat_sub, lat_mof, _l1, [],
                save_path=str(tmp_path / "r.pdf"),
            )

    def test_top_n_respected(self, _square_pair, _l1, tmp_path):
        """top_n=1 should still produce a valid PDF."""
        lat_sub, lat_mof = _square_pair
        matches = supercell.find_matches(
            lat_sub, lat_mof,
            theta_deg=_l1.theta_peaks[0],
            eta_tol=0.10,
        )
        out = tmp_path / "top1.pdf"
        path = generate_pdf_report(
            lat_sub, lat_mof, _l1, matches,
            save_path=str(out), top_n=1, dpi=72,
        )
        assert Path(path).exists()


# ---------------------------------------------------------------------------
# Hexagonal pair (broader coverage: more coincident LEED spots)
# ---------------------------------------------------------------------------

class TestHexagonalPair:
    """Run the four panel functions on an Au(111) / hex-MOF pair."""

    def setup_method(self):
        self.lat_sub = SUBSTRATE_DB["Au_111"]
        self.lat_mof = Lattice2D(4.92, 4.92, 120.0, "Hex-MOF")
        self.l1 = coincidence.compute(
            self.lat_sub, self.lat_mof, G_cutoff=8.0, sigma=0.4
        )
        self.matches = supercell.find_matches(
            self.lat_sub, self.lat_mof,
            theta_deg=self.l1.theta_peaks[0],
            eta_tol=0.10,
        )

    def test_phi_curve_hex(self):
        ax = plot_phi_curve(self.l1)
        assert ax is not None
        plt.close("all")

    def test_leed_hex(self):
        if not self.matches:
            pytest.skip("No matches found for hex pair")
        ax = plot_leed_pattern(self.lat_sub, self.lat_mof, self.matches[0])
        assert ax is not None
        plt.close("all")

    def test_overlay_hex(self):
        if not self.matches:
            pytest.skip("No matches found for hex pair")
        ax = plot_lattice_overlay(self.lat_sub, self.lat_mof, self.matches[0])
        assert ax is not None
        plt.close("all")

    def test_strain_hex(self):
        if not self.matches:
            pytest.skip("No matches found for hex pair")
        ax = plot_strain_ellipse(self.matches[0])
        assert ax is not None
        plt.close("all")


# ---------------------------------------------------------------------------
# _safe_stem helper
# ---------------------------------------------------------------------------

class TestSafeStem:
    def test_valid_name_unchanged(self):
        from miqrophi.visualize import _safe_stem
        assert _safe_stem("report") == "report"

    def test_spaces_replaced(self):
        from miqrophi.visualize import _safe_stem
        result = _safe_stem("my report title")
        assert " " not in result

    def test_special_chars_replaced(self):
        from miqrophi.visualize import _safe_stem
        result = _safe_stem("Au/111:test")
        assert "/" not in result
        assert ":" not in result

    def test_empty_string_fallback(self):
        from miqrophi.visualize import _safe_stem
        assert _safe_stem("") == "report"

    def test_only_invalid_chars_fallback(self):
        from miqrophi.visualize import _safe_stem
        assert _safe_stem("///") == "report"
