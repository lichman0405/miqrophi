"""Tests for Lattice2D dataclass and SUBSTRATE_DB."""
import math

import numpy as np
import pytest

from miqrocal.lattice import Lattice2D, SUBSTRATE_DB


class TestLattice2D:
    def test_square_matrix(self):
        lat = Lattice2D(4.0, 4.0, 90.0)
        np.testing.assert_allclose(lat.A, [[4.0, 0.0], [0.0, 4.0]], atol=1e-10)

    def test_square_area(self):
        lat = Lattice2D(4.0, 4.0, 90.0)
        assert math.isclose(lat.omega, 16.0, rel_tol=1e-10)

    def test_hexagonal_row1_x(self):
        lat = Lattice2D(2.88, 2.88, 120.0)
        assert math.isclose(lat.A[0, 0], 2.88, rel_tol=1e-10)
        assert math.isclose(lat.A[0, 1], 0.0,  abs_tol=1e-10)

    def test_hexagonal_row2_y(self):
        lat = Lattice2D(2.88, 2.88, 120.0)
        expected_y = 2.88 * math.sin(math.radians(120))
        assert math.isclose(lat.A[1, 1], expected_y, rel_tol=1e-6)

    def test_hexagonal_area(self):
        lat = Lattice2D(2.88, 2.88, 120.0)
        expected = 2.88 ** 2 * math.sin(math.radians(120))
        assert math.isclose(lat.omega, expected, rel_tol=1e-6)

    def test_recip_matrix_identity(self):
        """A · Bᵀ = 2π I for any lattice."""
        for gamma in (60.0, 90.0, 120.0):
            lat = Lattice2D(3.0, 4.0, gamma)
            B = lat.recip_matrix()
            product = lat.A @ B.T
            np.testing.assert_allclose(
                product, 2 * math.pi * np.eye(2), atol=1e-8,
                err_msg=f"Failed for gamma={gamma}",
            )

    def test_label_default_empty(self):
        lat = Lattice2D(3.0, 3.0, 90.0)
        assert lat.label == ""

    def test_label_set(self):
        lat = Lattice2D(3.0, 3.0, 90.0, "Au(100)")
        assert lat.label == "Au(100)"

    def test_oblique_area(self):
        lat = Lattice2D(3.0, 4.0, 70.0)
        expected = 3.0 * 4.0 * math.sin(math.radians(70.0))
        assert math.isclose(lat.omega, expected, rel_tol=1e-6)


class TestSubstrateDB:
    def test_total_count(self):
        assert len(SUBSTRATE_DB) == 11

    def test_expected_keys(self):
        expected = {
            "Au_111", "Au_100", "Ag_111", "Cu_111", "Pt_111",
            "Graphene", "HOPG", "ITO_111", "MgO_001", "SrTiO3_001", "Si_001",
        }
        assert set(SUBSTRATE_DB.keys()) == expected

    def test_hexagonal_gamma(self):
        for key in ("Au_111", "Ag_111", "Cu_111", "Pt_111", "Graphene", "HOPG", "ITO_111"):
            assert SUBSTRATE_DB[key].gamma_deg == 120.0, f"{key} should be hexagonal"

    def test_square_gamma(self):
        for key in ("Au_100", "MgO_001", "SrTiO3_001", "Si_001"):
            assert SUBSTRATE_DB[key].gamma_deg == 90.0, f"{key} should be square"

    def test_graphene_hopg_same_lattice(self):
        g  = SUBSTRATE_DB["Graphene"]
        hg = SUBSTRATE_DB["HOPG"]
        assert math.isclose(g.a, hg.a, rel_tol=1e-10)

    def test_all_entries_are_lattice2d(self):
        for key, val in SUBSTRATE_DB.items():
            assert isinstance(val, Lattice2D), f"{key} is not a Lattice2D"

    def test_positive_areas(self):
        for key, lat in SUBSTRATE_DB.items():
            assert lat.omega > 0, f"{key} has non-positive area"
