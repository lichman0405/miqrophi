"""
Coincidence filter — reciprocal-space overlap function Phi(theta).

Constructs the weighted overlap integral between two reciprocal-lattice
point sets as a function of the rotation angle theta:

    Phi(theta) = sum_{Gs, Gm} w_s * w_m
                 * exp(-|Gs - R_theta * Gm|^2 / (2 * sigma^2))

Weights w = exp(-(|G|/G_cutoff)^2) implement a Debye-Waller-like
suppression of high-frequency reflections.

Phi(theta) is expanded as a Fourier series via the Jacobi-Anger identity

    exp(z cos(phi)) = sum_n  I_n(z) exp(i*n*phi)

and all significant peaks are located by FFT-based evaluation.

Reference: docs/lattice_matching_theory.md, Section 4
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks
from scipy.special import iv as bessel_iv

from .lattice import Lattice2D


@dataclass
class CoincidenceResult:
    """Result of the Level-1 coincidence-function computation."""

    theta_peaks: list[float]  # peak positions (degrees), sorted by height desc
    phi_peaks:   list[float]  # normalised peak heights in [0, 1]
    phi_curve:   np.ndarray   # full normalised Phi(theta) curve, shape (K,)
    theta_grid:  np.ndarray   # corresponding angles (degrees), shape (K,)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _recip_points(
    lat: Lattice2D,
    G_cutoff: float,
    max_hk: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate the set of reciprocal-lattice points {G_hk} with |G| < G_cutoff.

    Returns
    -------
    G_vecs  : (N, 2) Cartesian coordinates in Angstrom^-1
    weights : (N,)   Debye-Waller weights exp(-(|G|/G_cutoff)^2)
    """
    B = lat.recip_matrix()
    vecs, wts = [], []
    for h in range(-max_hk, max_hk + 1):
        for k in range(-max_hk, max_hk + 1):
            if h == 0 and k == 0:
                continue
            G  = h * B[0] + k * B[1]
            Gn = float(np.linalg.norm(G))
            if Gn < G_cutoff:
                vecs.append(G)
                wts.append(np.exp(-(Gn / G_cutoff) ** 2))
    return np.array(vecs, dtype=float), np.array(wts, dtype=float)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    *,
    G_cutoff:        float = 8.0,
    sigma:           float = 0.3,
    n_max:           int   = 12,
    K:               int   = 2000,
    peak_threshold:  float = 0.15,
) -> CoincidenceResult:
    """
    Compute the coincidence function Phi(theta) and return all significant peaks.

    Parameters
    ----------
    G_cutoff        : reciprocal-space truncation radius (Angstrom^-1)
    sigma           : Gaussian width controlling angular tolerance (Angstrom^-1)
    n_max           : Fourier truncation order (n in [-n_max, n_max])
    K               : number of theta samples in [0, 360 deg)
    peak_threshold  : minimum relative peak height to report (fraction of max)
    """
    Gs, ws = _recip_points(lat_sub, G_cutoff)
    Gm, wm = _recip_points(lat_mof, G_cutoff)

    Gs_norm = np.linalg.norm(Gs, axis=1)          # (Ns,)
    Gm_norm = np.linalg.norm(Gm, axis=1)          # (Nm,)
    Gs_ang  = np.arctan2(Gs[:, 1], Gs[:, 0])      # (Ns,)
    Gm_ang  = np.arctan2(Gm[:, 1], Gm[:, 0])      # (Nm,)

    # Pairwise quantities, shape (Ns, Nm)
    A_ij     = np.outer(ws, wm) * np.exp(
        -(Gs_norm[:, None] ** 2 + Gm_norm[None, :] ** 2) / (2.0 * sigma ** 2)
    )
    z_ij     = np.outer(Gs_norm, Gm_norm) / sigma ** 2
    alpha_ij = Gs_ang[:, None] - Gm_ang[None, :]  # fixed inter-set angle

    # Fourier coefficients c_n via Jacobi-Anger expansion
    n_arr = np.arange(-n_max, n_max + 1)
    c = np.array(
        [np.sum(A_ij * bessel_iv(n, z_ij) * np.exp(-1j * n * alpha_ij))
         for n in n_arr],
        dtype=complex,
    )

    # Reconstruct Phi(theta) on a uniform grid (direct Fourier synthesis)
    theta_grid = np.linspace(0.0, 2.0 * np.pi, K, endpoint=False)
    Phi = np.real(
        np.einsum("n,nk->k", c,
                  np.exp(1j * n_arr[:, None] * theta_grid[None, :]))
    )
    phi_max  = Phi.max()
    Phi_norm = Phi / phi_max if phi_max > 0.0 else Phi

    # Peak detection
    peaks_idx, _ = find_peaks(
        Phi_norm,
        height=peak_threshold,
        distance=K // 36,   # minimum separation ~10 deg
    )
    order        = np.argsort(Phi_norm[peaks_idx])[::-1]
    theta_peaks  = np.degrees(theta_grid[peaks_idx[order]]).tolist()
    phi_peaks    = Phi_norm[peaks_idx[order]].tolist()

    return CoincidenceResult(
        theta_peaks=theta_peaks,
        phi_peaks=phi_peaks,
        phi_curve=Phi_norm,
        theta_grid=np.degrees(theta_grid),
    )
