"""
Level 2 — Floating-point LLL lattice reduction and strain tensor.

For each candidate rotation angle theta from Level 1, a 4-D embedding
lattice is constructed whose shortest vectors correspond to the optimal
integer supercell matrices (M, N).  The epitaxial strain is quantified
by the Green-Lagrange tensor eps = 0.5 * (F^T F - I), where
F = C_o . C_s^{-1} maps the substrate supercell onto the MOF supercell.

Embedding lattice basis (rows are basis vectors, lambda is penalty weight):

    B = [[ T[0,0],  T[0,1],  lam,  0  ],
         [ T[1,0],  T[1,1],   0,  lam ],
         [  -1,       0,       0,   0  ],
         [   0,      -1,       0,   0  ]]

Short vectors minimise  ||T.m - n||^2 + lam^2 * ||m||^2,
i.e. the trade-off between strain (residual) and supercell size.

Reference: docs/lattice_matching_theory.md, Section 5
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .lattice import Lattice2D


@dataclass
class MatchResult:
    """A single epitaxial supercell match produced by Level 2."""

    theta_deg:   float       # rotation angle used (degrees)
    M:           np.ndarray  # (2, 2) integer MOF supercell matrix
    N:           np.ndarray  # (2, 2) integer substrate supercell matrix
    strain:      np.ndarray  # (2, 2) Green-Lagrange strain tensor
    eta:         float       # Frobenius norm of strain, ||eps||_F
    area:        float       # substrate supercell area (Angstrom^2)
    lambda_used: float       # embedding parameter lambda that found this result


_DEFAULT_LAMBDAS: list[float] = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _lll_reduce(B: np.ndarray, delta: float = 0.75) -> np.ndarray:
    """
    Floating-point LLL lattice basis reduction (Gram-Schmidt variant).

    Parameters
    ----------
    B     : (d, d) lattice basis, rows are basis vectors
    delta : Lovasz parameter in (0.25, 1); default 0.75
    """
    B = B.copy().astype(float)
    d = len(B)

    def _gs(B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        Bs = np.zeros_like(B)
        mu = np.zeros((d, d))
        for i in range(d):
            Bs[i] = B[i].copy()
            for j in range(i):
                mu[i, j] = np.dot(B[i], Bs[j]) / np.dot(Bs[j], Bs[j])
                Bs[i] -= mu[i, j] * Bs[j]
        return Bs, mu

    k = 1
    while k < d:
        Bs, mu = _gs(B)
        for j in range(k - 1, -1, -1):
            if abs(mu[k, j]) > 0.5:
                B[k] -= round(float(mu[k, j])) * B[j]
                Bs, mu = _gs(B)
        lovasz = np.dot(Bs[k], Bs[k]) >= (delta - mu[k, k-1] ** 2) * np.dot(Bs[k-1], Bs[k-1])
        if lovasz:
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]
            k = max(k - 1, 1)
    return B


def _embedding_basis(T: np.ndarray, lam: float) -> np.ndarray:
    """Build the 4-D embedding lattice basis for transfer matrix T and weight lam."""
    return np.array(
        [[ T[0, 0],  T[0, 1],  lam,  0.0],
         [ T[1, 0],  T[1, 1],  0.0,  lam],
         [-1.0,      0.0,       0.0,  0.0],
         [ 0.0,     -1.0,       0.0,  0.0]],
        dtype=float,
    )


def _extract_M(reduced_B: np.ndarray, lam: float) -> Optional[np.ndarray]:
    """
    Extract the 2x2 integer MOF supercell matrix M from the LLL-reduced basis.

    Each reduced basis vector v satisfies v[2:4] = lam * (m1, m2), so
    m = round(v[2:4] / lam).  Two linearly independent rows are required
    (det(M) != 0).
    """
    candidates: list[np.ndarray] = []
    for v in reduced_B:
        m = np.round(v[2:4] / lam).astype(int)
        if not np.all(m == 0):
            candidates.append(m)
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            M = np.array([candidates[i], candidates[j]])
            if abs(round(float(np.linalg.det(M.astype(float))))) >= 1:
                return M
    return None


def _green_lagrange(
    M: np.ndarray,
    N: np.ndarray,
    A_o_rot: np.ndarray,
    A_s: np.ndarray,
) -> tuple[np.ndarray, float]:
    """
    Compute the Green-Lagrange strain tensor eps = 0.5 * (F^T F - I) and
    its Frobenius norm eta = ||eps||_F.

    F = C_o . C_s^{-1} is the deformation gradient that maps the substrate
    supercell C_s = N . A_s onto the (rotated) MOF supercell C_o = M . A_o_rot.
    The Green-Lagrange form is rotation-invariant, so equivalent orientation
    domains yield the same eta.
    """
    C_o = M.astype(float) @ A_o_rot
    C_s = N.astype(float) @ A_s
    try:
        F = C_o @ np.linalg.inv(C_s)
    except np.linalg.LinAlgError:
        inf2 = np.full((2, 2), np.inf)
        return inf2, np.inf
    eps = 0.5 * (F.T @ F - np.eye(2))
    return eps, float(np.linalg.norm(eps, "fro"))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def find_matches(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    theta_deg: float,
    *,
    lambda_values: list[float] = _DEFAULT_LAMBDAS,
    eta_tol: float = 0.05,
) -> list[MatchResult]:
    """
    Find integer supercell pairs (M, N) for the given rotation angle theta.

    Multiple lambda values are scanned to expose the strain-vs-supercell-size
    Pareto front.  Results with eta < eta_tol are returned, sorted by eta.
    """
    th = np.radians(theta_deg)
    R  = np.array([[np.cos(th), -np.sin(th)],
                   [np.sin(th),  np.cos(th)]])

    A_o_rot = R @ lat_mof.A                          # rotated MOF lattice matrix
    T       = A_o_rot @ np.linalg.inv(lat_sub.A)    # transfer matrix

    results: list[MatchResult] = []
    seen:    set[tuple]        = set()

    for lam in lambda_values:
        B4     = _embedding_basis(T, lam)
        B4_red = _lll_reduce(B4)
        M      = _extract_M(B4_red, lam)
        if M is None:
            continue

        key = tuple(M.flatten())
        if key in seen:
            continue
        seen.add(key)

        N = np.round(M.astype(float) @ T).astype(int)
        if abs(round(float(np.linalg.det(N.astype(float))))) < 1:
            continue

        eps, eta = _green_lagrange(M, N, A_o_rot, lat_sub.A)
        if eta > eta_tol:
            continue

        area = abs(float(np.linalg.det(N.astype(float) @ lat_sub.A)))
        results.append(MatchResult(
            theta_deg=theta_deg,
            M=M, N=N,
            strain=eps, eta=eta,
            area=area, lambda_used=lam,
        ))

    results.sort(key=lambda r: r.eta)
    return results
