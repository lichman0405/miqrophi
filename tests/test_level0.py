"""Tests for Level 0 discriminant check."""
import pytest

from miqrocal.lattice import Lattice2D, SUBSTRATE_DB
from miqrocal import level0


class TestLevel0Feasibility:
    def test_square_on_square_feasible(self):
        lat_sq = Lattice2D(18.62, 18.62, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_sq)
        assert result.feasible is True

    def test_hexagonal_on_hexagonal_feasible(self):
        lat_hex = Lattice2D(4.92, 4.92, 120.0)
        result = level0.check(SUBSTRATE_DB["Au_111"], lat_hex)
        assert result.feasible is True

    def test_square_on_hexagonal_forbidden(self):
        lat_sq = Lattice2D(18.62, 18.62, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_111"], lat_sq)
        assert result.feasible is False

    def test_hexagonal_on_square_forbidden(self):
        lat_hex = Lattice2D(4.92, 4.92, 120.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_hex)
        assert result.feasible is False

    def test_different_scale_square_still_feasible(self):
        """Scaling should not affect feasibility (discriminant is scale-independent)."""
        lat_large = Lattice2D(100.0, 100.0, 90.0)
        lat_small = Lattice2D(1.0, 1.0, 90.0)
        result = level0.check(lat_small, lat_large)
        assert result.feasible is True

    def test_rectangular_on_square_feasible(self):
        """Rectangular (a≠b, γ=90°) on square: same quadratic field."""
        lat_rect = Lattice2D(5.0, 10.0, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_rect)
        assert result.feasible is True


class TestLevel0Result:
    def test_delta_square(self):
        import math
        lat_sq = Lattice2D(4.0, 4.0, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_sq)
        assert abs(result.delta_sub - (-4.0)) < 1e-6
        assert abs(result.delta_mof - (-4.0)) < 1e-6

    def test_delta_hexagonal(self):
        import math
        lat_hex = Lattice2D(2.88, 2.88, 120.0)
        result = level0.check(SUBSTRATE_DB["Au_111"], lat_hex)
        assert abs(result.delta_sub - (-3.0)) < 1e-4
        assert abs(result.delta_mof - (-3.0)) < 1e-4

    def test_ratio_same_type_is_one(self):
        lat_sq = Lattice2D(10.0, 10.0, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_sq)
        assert abs(result.ratio - 1.0) < 1e-6

    def test_L0_feasible_appears_in_str(self):
        lat_sq = Lattice2D(18.62, 18.62, 90.0)
        result = level0.check(SUBSTRATE_DB["Au_100"], lat_sq)
        # result should have a non-empty message
        assert len(result.message) > 0


class TestLevel0KnownCases:
    def test_HKUST1_Au100_feasible(self):
        """HKUST-1(010) on Au(100): both square → feasible."""
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        result = level0.check(SUBSTRATE_DB["Au_100"], hkust1)
        assert result.feasible is True

    def test_HKUST1_Au111_forbidden(self):
        """HKUST-1(010) on Au(111): square vs hexagonal → forbidden."""
        hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
        result = level0.check(SUBSTRATE_DB["Au_111"], hkust1)
        assert result.feasible is False
