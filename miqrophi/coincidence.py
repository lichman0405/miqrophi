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
from functools import lru_cache

import numpy as np
from scipy.signal import find_peaks
from scipy.special import iv as bessel_iv

from .lattice import Lattice2D


@dataclass
class CoincidenceResult:
    """Result of the Level-1 coincidence-function computation."""

    theta_peaks: list[float]  # peak positions (degrees), sorted by Φ desc then θ asc
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

    Implementation note
    -------------------
    Fully vectorised with NumPy meshgrid — no Python-level loops, replacing
    the previous O(max_hk^2) list-append loop.
    """
    B = lat.recip_matrix()
    idx = np.arange(-max_hk, max_hk + 1)
    hh, kk = np.meshgrid(idx, idx, indexing='ij')
    hh, kk = hh.ravel(), kk.ravel()
    # exclude the Gamma point (0, 0)
    nonzero = (hh != 0) | (kk != 0)
    hh, kk = hh[nonzero], kk[nonzero]
    G_vecs  = hh[:, None] * B[0] + kk[:, None] * B[1]   # (N_cand, 2)
    G_norms = np.linalg.norm(G_vecs, axis=1)             # (N_cand,)
    inside  = G_norms < G_cutoff
    G_vecs  = G_vecs[inside]
    G_norms = G_norms[inside]
    weights = np.exp(-(G_norms / G_cutoff) ** 2)
    return G_vecs, weights


@lru_cache(maxsize=128)
def _recip_points_cached(
    a: float, b: float, gamma_deg: float,
    G_cutoff: float, max_hk: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Cached wrapper around _recip_points.  In batch runs the same substrate
    (or MOF surface) is queried many times — caching avoids redundant work.
    The returned arrays are read-only and must not be mutated.
    """
    lat = Lattice2D(a, b, gamma_deg)
    G_vecs, weights = _recip_points(lat, G_cutoff, max_hk)
    G_vecs.flags.writeable = False
    weights.flags.writeable = False
    return G_vecs, weights


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
    Gs, ws = _recip_points_cached(lat_sub.a, lat_sub.b, lat_sub.gamma_deg, G_cutoff)
    Gm, wm = _recip_points_cached(lat_mof.a, lat_mof.b, lat_mof.gamma_deg, G_cutoff)

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
    # ── sparse magnitude pre-filter ─────────────────────────────────────────
    # For large z = |Gs||Gm|/σ², the maximum contribution of a pair is:
    #   A_ij * I_n(z_ij) ~ exp(-( |Gs| - |Gm| )^2 / (2*sigma^2))
    # Pairs with ||Gs|-|Gm|| > _MAG_SIGMA_MULT * sigma are negligible and
    # skipped — the Bessel evaluation is only done for the surviving pairs.
    _MAG_SIGMA_MULT = 5.0
    mag_diff = np.abs(Gs_norm[:, None] - Gm_norm[None, :])  # (Ns, Nm)
    sig_s, sig_m = np.where(mag_diff < _MAG_SIGMA_MULT * sigma)

    n_arr = np.arange(-n_max, n_max + 1)                            # (2*n_max+1,)
    if len(sig_s) == 0:
        c = np.zeros(len(n_arr), dtype=complex)
    else:
        A_flat = A_ij[sig_s, sig_m]                                 # (Npairs,)
        z_flat = z_ij[sig_s, sig_m]
        α_flat = alpha_ij[sig_s, sig_m]

        # Exploit I_n(z) = I_{-n}(z): compute only non-negative orders.
        # This halves the Bessel function evaluations vs computing all
        # 2*n_max+1 orders independently.
        n_pos  = np.arange(0, n_max + 1)                            # (n_max+1,)
        I_pos  = bessel_iv(n_pos[:, None], z_flat[None, :])        # (n_max+1, Npairs)
        ph_pos = np.exp(-1j * n_pos[:, None] * α_flat[None, :])    # (n_max+1, Npairs)
        # c_pos[k] = Σ_p A_p · I_k(z_p) · exp(-i·k·α_p)
        c_pos  = (I_pos * ph_pos) @ A_flat                         # (n_max+1,)
        # Negative orders: c_{-k} = conj(c_k)  (because A_flat is real)
        c = np.concatenate([np.conj(c_pos[1:][::-1]), c_pos])      # (2*n_max+1,)

    # Reconstruct Phi(theta) via zero-padded IFFT.
    # Phi(theta_k) = sum_n c_n exp(i*n*theta_k), theta_k = 2*pi*k/K
    # Place c_n into the DFT layout recognised by numpy.fft.ifft:
    #   positive n -> indices 0 .. n_max
    #   negative n -> indices K-n_max .. K-1
    theta_grid = np.linspace(0.0, 2.0 * np.pi, K, endpoint=False)
    c_fft = np.zeros(K, dtype=complex)
    c_fft[:n_max + 1]  = c[n_max:]   # c_0 .. c_{n_max}
    c_fft[K - n_max:]  = c[:n_max]   # c_{-n_max} .. c_{-1}
    Phi = np.real(np.fft.ifft(c_fft)) * K
    phi_max  = Phi.max()
    Phi_norm = Phi / phi_max if phi_max > 0.0 else Phi

    # Peak detection
    peaks_idx, _ = find_peaks(
        Phi_norm,
        height=peak_threshold,
        distance=K // 36,   # minimum separation ~10 deg
    )
    # Sort primarily by descending phi height, secondarily by ascending theta
    # (stable secondary key) so that equal-phi peaks are always in the same
    # deterministic order regardless of NumPy version or floating-point jitter.
    order = np.lexsort((np.degrees(theta_grid[peaks_idx]),
                        -Phi_norm[peaks_idx]))
    theta_peaks  = np.degrees(theta_grid[peaks_idx[order]]).tolist()
    phi_peaks    = Phi_norm[peaks_idx[order]].tolist()

    return CoincidenceResult(
        theta_peaks=theta_peaks,
        phi_peaks=phi_peaks,
        phi_curve=Phi_norm,
        theta_grid=np.degrees(theta_grid),
    )
