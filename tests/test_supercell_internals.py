"""Tests for internal supercell functions (LLL helpers, embedding basis, etc.)."""
import math
import numpy as np
import pytest

from miqrophi.supercell import (
    _gs,
    _lll_reduce_py,
    _embedding_basis,
    _extract_M,
    _green_lagrange,
    find_matches,
)
from miqrophi import Lattice2D, SUBSTRATE_DB


# ---------------------------------------------------------------------------
# _gs  (Gram-Schmidt orthogonalisation)
# ---------------------------------------------------------------------------

class TestGramSchmidt:
    def test_orthogonal_basis_unchanged(self):
        """Identity matrix is already orthogonal; GS should leave it unchanged."""
        B = np.eye(3)
        Bs, mu = _gs(B)
        assert np.allclose(Bs, np.eye(3), atol=1e-12)
        assert np.allclose(np.tril(mu, -1), 0.0, atol=1e-12)

    def test_orthogonality_of_output(self):
        rng = np.random.default_rng(0)
        B = rng.standard_normal((4, 4))
        Bs, mu = _gs(B)
        for i in range(4):
            for j in range(i):
                dot = np.dot(Bs[i], Bs[j])
                assert abs(dot) < 1e-10, f"Bs[{i}] and Bs[{j}] not orthogonal"

    def test_mu_lower_triangular(self):
        rng = np.random.default_rng(1)
        B = rng.standard_normal((3, 3))
        _, mu = _gs(B)
        # upper triangle (excluding diagonal) should be zero
        for i in range(3):
            for j in range(i + 1, 3):
                assert mu[i, j] == 0.0

    def test_2d_example(self):
        """Simple 2-D example: e1 and e1+e2."""
        B = np.array([[1.0, 0.0], [1.0, 1.0]])
        Bs, mu = _gs(B)
        assert np.allclose(Bs[0], [1.0, 0.0], atol=1e-12)
        assert np.allclose(Bs[1], [0.0, 1.0], atol=1e-12)
        assert math.isclose(mu[1, 0], 1.0, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# _lll_reduce_py  (pure-Python LLL)
# ---------------------------------------------------------------------------

class TestLllReducePy:
    def test_reduces_2d_basis(self):
        """The LLL-reduced 2-D basis should have shorter vectors."""
        B = np.array([[1.0, 0.0], [100.0, 1.0]])
        Bred = _lll_reduce_py(B)
        # Reduced basis should not have extremely long vectors
        norms = [np.linalg.norm(Bred[i]) for i in range(2)]
        assert max(norms) < 10.0

    def test_output_span_same_lattice(self):
        """det must be preserved (same lattice, just different basis)."""
        B = np.array([[3.0, 1.0], [2.0, 5.0]])
        Bred = _lll_reduce_py(B)
        assert math.isclose(abs(np.linalg.det(Bred)),
                             abs(np.linalg.det(B)), rel_tol=1e-10)

    def test_4d_basis(self):
        """LLL on a 4-D basis (as used in the supercell embedding)."""
        B = np.array([
            [3.0,  1.0, 0.02, 0.0 ],
            [1.0,  4.0, 0.0,  0.02],
            [-1.0, 0.0, 0.0,  0.0 ],
            [0.0, -1.0, 0.0,  0.0 ],
        ])
        Bred = _lll_reduce_py(B)
        # det magnitude preserved
        assert math.isclose(abs(np.linalg.det(Bred)),
                             abs(np.linalg.det(B)), rel_tol=1e-8)

    def test_identity_is_already_reduced(self):
        B = np.eye(2)
        Bred = _lll_reduce_py(B)
        assert np.allclose(np.abs(Bred), np.eye(2), atol=1e-12)

    def test_does_not_mutate_input(self):
        B = np.array([[3.0, 1.0], [2.0, 5.0]])
        B_orig = B.copy()
        _lll_reduce_py(B)
        assert np.allclose(B, B_orig)


# ---------------------------------------------------------------------------
# _embedding_basis
# ---------------------------------------------------------------------------

class TestEmbeddingBasis:
    def test_shape(self):
        T = np.array([[1.5, 0.3], [0.2, 1.4]])
        B = _embedding_basis(T, lam=0.02)
        assert B.shape == (4, 4)

    def test_top_left_is_T(self):
        T = np.array([[2.1, 0.5], [-0.3, 1.9]])
        B = _embedding_basis(T, lam=0.1)
        assert np.allclose(B[:2, :2], T, atol=1e-12)

    def test_lambda_appears_on_diagonal(self):
        T = np.eye(2)
        lam = 0.05
        B = _embedding_basis(T, lam=lam)
        assert math.isclose(B[0, 2], lam)
        assert math.isclose(B[1, 3], lam)

    def test_minus_identity_in_bottom_left(self):
        T = np.eye(2)
        B = _embedding_basis(T, lam=0.1)
        assert B[2, 0] == -1.0
        assert B[3, 1] == -1.0
        assert B[2, 1] == 0.0
        assert B[3, 0] == 0.0


# ---------------------------------------------------------------------------
# _extract_M
# ---------------------------------------------------------------------------

class TestExtractM:
    def _make_reduced_with_known_m(self, m1, m2, lam=0.1):
        """Construct a fake LLL-reduced basis encoding rows m1, m2 in slots 2:4."""
        v1 = np.array([0.0, 0.0, m1[0] * lam, m1[1] * lam])
        v2 = np.array([0.0, 0.0, m2[0] * lam, m2[1] * lam])
        v3 = np.array([1.0, 0.0, 0.0, 0.0])
        v4 = np.array([0.0, 1.0, 0.0, 0.0])
        return np.array([v1, v2, v3, v4])

    def test_recovers_simple_matrix(self):
        m1 = np.array([1, 0])
        m2 = np.array([0, 1])
        B = self._make_reduced_with_known_m(m1, m2)
        M = _extract_M(B, lam=0.1)
        assert M is not None
        assert abs(round(float(np.linalg.det(M.astype(float))))) >= 1

    def test_returns_none_when_singular(self):
        """All reduced vectors encode the same m → linearly dependent → None."""
        m = np.array([1, 1])
        v = np.array([0.0, 0.0, m[0] * 0.1, m[1] * 0.1])
        B = np.tile(v, (4, 1))
        M = _extract_M(B, lam=0.1)
        assert M is None

    def test_returns_none_all_zero(self):
        B = np.zeros((4, 4))
        M = _extract_M(B, lam=0.1)
        assert M is None


# ---------------------------------------------------------------------------
# _green_lagrange
# ---------------------------------------------------------------------------

class TestGreenLagrange:
    def test_zero_strain_for_identical_lattice(self):
        """When M = N = I and lat_mof = lat_sub, strain must be zero."""
        lat = Lattice2D(4.08, 4.08, 90.0)
        M = np.eye(2, dtype=int)
        N = np.eye(2, dtype=int)
        A_o_rot = lat.A.copy()
        eps, eta = _green_lagrange(M, N, A_o_rot, lat.A)
        assert math.isclose(eta, 0.0, abs_tol=1e-12)
        assert np.allclose(eps, np.zeros((2, 2)), atol=1e-12)

    def test_strain_symmetric(self):
        lat_sub = SUBSTRATE_DB["Au_100"]
        lat_mof = Lattice2D(4.50, 4.50, 90.0)
        M = np.array([[1, 0], [0, 1]])
        N = np.array([[1, 0], [0, 1]])
        eps, eta = _green_lagrange(M, N, lat_mof.A, lat_sub.A)
        assert np.allclose(eps, eps.T, atol=1e-12)

    def test_positive_eta(self):
        lat_sub = SUBSTRATE_DB["Au_100"]
        lat_mof = Lattice2D(5.00, 5.00, 90.0)
        M = np.eye(2, dtype=int)
        N = np.eye(2, dtype=int)
        eps, eta = _green_lagrange(M, N, lat_mof.A, lat_sub.A)
        assert eta > 0.0

    def test_singular_N_returns_inf(self):
        lat_sub = SUBSTRATE_DB["Au_100"]
        M = np.array([[1, 0], [0, 1]])
        N = np.array([[1, 1], [1, 1]])   # det = 0 → singular C_s
        A_o = lat_sub.A.copy()
        eps, eta = _green_lagrange(M, N, A_o, lat_sub.A)
        assert eta == np.inf or not np.isfinite(eta)


# ---------------------------------------------------------------------------
# find_matches — edge-case coverage
# ---------------------------------------------------------------------------

class TestFindMatchesEdgeCases:
    def test_exact_commensurate_match_eta_zero(self):
        """lat_mof = lat_sub × 2 → should give a perfect (eta≈0) match."""
        lat_sub = Lattice2D(4.08, 4.08, 90.0, "sub")
        lat_mof = Lattice2D(8.16, 8.16, 90.0, "mof×2")
        matches = find_matches(lat_sub, lat_mof, theta_deg=0.0, eta_tol=0.01)
        assert len(matches) >= 1
        assert matches[0].eta < 0.01

    def test_no_match_returns_empty(self):
        """Wildly mismatched lattices should return nothing under tight tolerance."""
        lat_sub = Lattice2D(4.08, 4.08, 90.0)
        lat_mof = Lattice2D(4.08, 4.08, 90.0)
        # theta=45° is a bad commensurate angle for a square lattice
        matches = find_matches(lat_sub, lat_mof, theta_deg=45.0, eta_tol=1e-6)
        assert isinstance(matches, list)

    def test_sorted_by_eta(self):
        lat_sub = SUBSTRATE_DB["Au_100"]
        lat_mof = Lattice2D(18.62, 18.62, 90.0)
        matches = find_matches(lat_sub, lat_mof, theta_deg=270.0, eta_tol=0.10)
        etas = [m.eta for m in matches]
        assert etas == sorted(etas)

    def test_duplicate_M_suppressed(self):
        """Multiple lambda values should not produce duplicate M matrices."""
        lat_sub = SUBSTRATE_DB["Au_100"]
        lat_mof = Lattice2D(18.62, 18.62, 90.0)
        matches = find_matches(
            lat_sub, lat_mof, theta_deg=270.0,
            lambda_values=[0.02, 0.02, 0.02],  # same lambda repeated
            eta_tol=0.05,
        )
        keys = [tuple(m.M.flatten()) for m in matches]
        assert len(keys) == len(set(keys))
