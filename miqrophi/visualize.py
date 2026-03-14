"""
Visualisation utilities for miqrophi epitaxy matching results.

Five plotting functions, one per plan:

    plot_phi_curve       — Level-1 coincidence function Φ(θ)           [Plan B]
    plot_lattice_overlay — real-space supercell overlay                 [Plan A]
    plot_leed_pattern    — simulated LEED diffraction pattern           [Plan C]
    plot_strain_ellipse  — principal-strain ellipse and tensor summary  [Plan D]
    plot_match_card      — 2×2 summary figure combining all four panels [Plan E]

All single-panel functions accept an optional ``ax`` keyword argument so they
can be embedded in larger user-defined figures.

References
----------
docs/lattice_matching_theory.md, Sections 3-5
"""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .coincidence import CoincidenceResult
from .lattice import Lattice2D
from .supercell import MatchResult

if TYPE_CHECKING:
    from matplotlib.animation import FuncAnimation

# ---------------------------------------------------------------------------
# Colour palette  (Tableau-10, colour-blind-friendly)
# ---------------------------------------------------------------------------
_C_SUB  = "#1f77b4"   # blue   — substrate lattice
_C_MOF  = "#d62728"   # red    — MOF overlayer
_C_CELL = "#2ca02c"   # green  — supercell boundary / selected angle
_C_PHI  = "#9467bd"   # purple — coincidence function curve
_C_PEAK = "#ff7f0e"   # orange — Phi peaks
_C_COIN = "#2ca02c"   # green  — coincident LEED reflections


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

def _rot(theta_deg: float) -> np.ndarray:
    """2D counter-clockwise rotation matrix for *theta_deg* degrees."""
    t = np.radians(theta_deg)
    return np.array([[np.cos(t), -np.sin(t)],
                     [np.sin(t),  np.cos(t)]])


def _real_points_in_box(
    A: np.ndarray,
    half_box: float,
    n_max: int = 30,
) -> np.ndarray:
    """
    Return Cartesian coordinates of lattice points P = h·A[0] + k·A[1]
    for |h|, |k| ≤ n_max that fall within the square [-half_box, half_box]².
    """
    idx = np.arange(-n_max, n_max + 1)
    hh, kk = np.meshgrid(idx, idx, indexing='ij')
    pts = hh.ravel()[:, None] * A[0] + kk.ravel()[:, None] * A[1]  # (N, 2)
    mask = (np.abs(pts[:, 0]) <= half_box) & (np.abs(pts[:, 1]) <= half_box)
    pts = pts[mask]
    return pts if len(pts) > 0 else np.empty((0, 2))


def _recip_points_in_disk(
    B: np.ndarray,
    G_max: float,
    n_max: int = 20,
) -> np.ndarray:
    """
    Return reciprocal lattice points G = h·B[0] + k·B[1] with |G| < G_max,
    excluding the Gamma point (0, 0).
    """
    idx = np.arange(-n_max, n_max + 1)
    hh, kk = np.meshgrid(idx, idx, indexing='ij')
    hh, kk = hh.ravel(), kk.ravel()
    nonzero = (hh != 0) | (kk != 0)
    hh, kk = hh[nonzero], kk[nonzero]
    pts = hh[:, None] * B[0] + kk[:, None] * B[1]   # (N, 2)
    pts = pts[np.linalg.norm(pts, axis=1) < G_max]
    return pts if len(pts) > 0 else np.empty((0, 2))


# ---------------------------------------------------------------------------
# Plan B — Φ(θ) coincidence function
# ---------------------------------------------------------------------------

def plot_phi_curve(
    l1: CoincidenceResult,
    *,
    ax: Optional[Axes] = None,
    title: str = "Coincidence Function Φ(θ)",
) -> Axes:
    """
    Plot the normalised Level-1 coincidence function Φ(θ) with significant
    peaks annotated.

    The curve is the output of :func:`miqrophi.coincidence.compute`.  Peaks
    correspond to rotation angles where the reciprocal-lattice point sets
    show maximum overlap, i.e. candidate epitaxial orientations.

    Parameters
    ----------
    l1    : CoincidenceResult from ``level1.compute()``
    ax    : existing Axes to draw into; a new figure is created if *None*
    title : panel title
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 3))

    ax.plot(l1.theta_grid, l1.phi_curve, color=_C_PHI, lw=1.2, zorder=2)
    ax.fill_between(l1.theta_grid, l1.phi_curve, alpha=0.12, color=_C_PHI)

    for th, ph in zip(l1.theta_peaks[:5], l1.phi_peaks[:5]):
        ax.axvline(th, color=_C_PEAK, lw=0.9, ls="--", alpha=0.7, zorder=1)
        ax.text(
            th, min(ph + 0.04, 1.08),
            f"{th:.1f}°",
            ha="center", va="bottom", fontsize=7, color=_C_PEAK,
        )

    ax.set_xlabel("Rotation angle θ (deg)")
    ax.set_ylabel("Normalised Φ(θ)")
    ax.set_title(title, fontsize=9)
    ax.set_xlim(l1.theta_grid[0], l1.theta_grid[-1])
    ax.set_ylim(0, 1.18)
    ax.grid(True, lw=0.4, alpha=0.45)
    return ax


# ---------------------------------------------------------------------------
# Plan A — real-space lattice overlay
# ---------------------------------------------------------------------------

def plot_lattice_overlay(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    match: MatchResult,
    *,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
) -> Axes:
    """
    Draw real-space substrate and (rotated) MOF lattice points together with
    the coincidence supercell boundary.

    The supercell parallelogram is defined by the integer matrix N acting on
    the substrate basis vectors: C_s = N · A_s.

    Parameters
    ----------
    lat_sub, lat_mof : substrate and MOF ``Lattice2D`` objects
    match            : ``MatchResult`` from ``level2.find_matches()``
    ax               : existing Axes; a new figure is created if *None*
    title            : panel title (auto-generated if *None*)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    R     = _rot(match.theta_deg)
    A_s   = lat_sub.A.copy()
    A_mof = (R @ lat_mof.A.T).T          # rotate MOF basis rows by R

    # Supercell vectors: rows of N @ A_s
    C_s    = match.N.astype(float) @ A_s
    v1, v2 = C_s[0], C_s[1]
    span   = max(
        np.linalg.norm(v1),
        np.linalg.norm(v2),
        np.linalg.norm(v1 + v2),
    ) * 1.35
    span = max(span, 5.0)               # ensure a minimum visible window

    pts_s = _real_points_in_box(A_s,   span)
    pts_m = _real_points_in_box(A_mof, span)

    ax.scatter(pts_s[:, 0], pts_s[:, 1],
               s=18, c=_C_SUB, zorder=3, alpha=0.80,
               label=f"{lat_sub.label or 'Substrate'}")
    ax.scatter(pts_m[:, 0], pts_m[:, 1],
               s=18, c=_C_MOF, zorder=3, alpha=0.80, marker="^",
               label=f"{lat_mof.label or 'Overlayer'} (θ={match.theta_deg:.1f}°)")

    # Supercell parallelogram
    corners = np.array([[0, 0], v1, v1 + v2, v2, [0, 0]])
    ax.plot(corners[:, 0], corners[:, 1],
            "-", color=_C_CELL, lw=1.8, zorder=4,
            label=f"Supercell  ({match.area:.0f} Å²)")

    ax.scatter([0], [0], s=40, c="black", zorder=5)

    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_aspect("equal")
    ax.axhline(0, lw=0.3, c="grey")
    ax.axvline(0, lw=0.3, c="grey")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_title(title or f"Real-space overlay  η = {match.eta:.5f}", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, lw=0.3, alpha=0.35)
    return ax


# ---------------------------------------------------------------------------
# Plan C — simulated LEED pattern
# ---------------------------------------------------------------------------

def plot_leed_pattern(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    match: MatchResult,
    *,
    ax: Optional[Axes] = None,
    G_max: Optional[float] = None,
    title: Optional[str] = None,
) -> Axes:
    """
    Plot a simulated LEED diffraction pattern in reciprocal space.

    Substrate spots are drawn as blue circles; rotated MOF spots as red
    triangles.  Points that coincide within a proximity tolerance are
    highlighted as green stars — these are the superstructure reflections
    visible in experiment.

    Parameters
    ----------
    G_max : reciprocal-space cutoff radius (Å⁻¹).
            Defaults to 3.2 × |b₁| of the substrate.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    R     = _rot(match.theta_deg)
    B_s   = lat_sub.recip_matrix()
    # Reciprocal vectors transform identically to real-space vectors under R
    B_mof = (R @ lat_mof.recip_matrix().T).T

    if G_max is None:
        G_max = 3.2 * float(np.linalg.norm(B_s[0]))

    Gs = _recip_points_in_disk(B_s,   G_max)
    Gm = _recip_points_in_disk(B_mof, G_max)

    # Coincidence: |Gs_i − Gm_j| < tol
    tol = 0.18 * float(np.linalg.norm(B_s[0]))
    coin_s: list[int] = []
    coin_m: list[int] = []
    if len(Gs) and len(Gm):
        dist = np.linalg.norm(
            Gs[:, None, :] - Gm[None, :, :], axis=-1
        )                                 # shape (N_s, N_m)
        ci, cj = np.where(dist < tol)
        coin_s, coin_m = list(ci), list(cj)

    # Substrate spots
    if len(Gs):
        mask_s = np.ones(len(Gs), dtype=bool)
        if coin_s:
            mask_s[coin_s] = False
        ax.scatter(
            Gs[mask_s, 0], Gs[mask_s, 1],
            s=28, c=_C_SUB, zorder=3, alpha=0.85,
            label=f"{lat_sub.label or 'Substrate'}",
        )

    # MOF spots
    if len(Gm):
        mask_m = np.ones(len(Gm), dtype=bool)
        if coin_m:
            mask_m[coin_m] = False
        ax.scatter(
            Gm[mask_m, 0], Gm[mask_m, 1],
            s=28, c=_C_MOF, marker="^", zorder=3, alpha=0.85,
            label=f"{lat_mof.label or 'Overlayer'} (rotated)",
        )

    # Superstructure (coincident) spots
    if coin_s:
        ax.scatter(
            Gs[coin_s, 0], Gs[coin_s, 1],
            s=90, c=_C_COIN, marker="*", zorder=5,
            label=f"Superstructure ({len(coin_s)} spots)",
        )

    # (00) beam
    ax.scatter([0], [0], s=90, c="black", zorder=6, label="(00)")

    ax.set_aspect("equal")
    ax.set_xlabel(r"$k_x$  (Å$^{-1}$)")
    ax.set_ylabel(r"$k_y$  (Å$^{-1}$)")
    ax.set_title(title or f"LEED pattern  θ = {match.theta_deg:.1f}°", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, lw=0.3, alpha=0.35)
    return ax


# ---------------------------------------------------------------------------
# Plan D — strain tensor visualisation
# ---------------------------------------------------------------------------

def plot_strain_ellipse(
    match: MatchResult,
    *,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
) -> Axes:
    """
    Visualise the Green-Lagrange strain tensor ε.

    A reference unit circle (grey dashed) is compared with the deformed
    shape (I + ε)·circle (red solid).  Arrows indicate the two principal
    strain directions with their eigenvalues.  A summary box lists all
    tensor components.

    The first-order mapping p → (I + ε)p is used for visualisation; this
    is accurate whenever η = ‖ε‖_F ≪ 1, which is the physically relevant
    regime.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))

    eps = match.strain

    # Reference unit circle
    phi_arr = np.linspace(0, 2 * np.pi, 360)
    circle  = np.array([np.cos(phi_arr), np.sin(phi_arr)])   # (2, 360)
    ax.plot(circle[0], circle[1], ":", color="grey", lw=1.0,
            zorder=1, label="Reference (unstrained)")

    # Deformed ellipse
    deformed = (np.eye(2) + eps) @ circle
    ax.plot(deformed[0], deformed[1], "-", color=_C_MOF, lw=1.8,
            zorder=2, label="Strained shape")

    # Principal strains via eigendecomposition of symmetric ε
    eigvals, eigvecs = np.linalg.eigh(eps)   # ascending order
    for val, vec in zip(eigvals, eigvecs.T):
        tip    = (1.0 + val) * vec
        offset = vec * 0.20
        ax.annotate(
            "", xy=tip, xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=_C_PEAK, lw=2.0),
            zorder=3,
        )
        ax.text(
            tip[0] + offset[0], tip[1] + offset[1],
            f"ε = {val:+.4f}",
            ha="center", va="center", fontsize=7.5, color=_C_PEAK,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75),
        )

    ax.axhline(0, lw=0.4, c="grey")
    ax.axvline(0, lw=0.4, c="grey")

    r = max(1.30, float(np.abs(deformed).max()) * 1.25)
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    summary = (
        f"ε₁₁ = {eps[0, 0]:+.5f}\n"
        f"ε₂₂ = {eps[1, 1]:+.5f}\n"
        f"ε₁₂ = {eps[0, 1]:+.5f}\n"
        f"η   = {match.eta:.5f}"
    )
    ax.text(
        0.04, 0.04, summary,
        transform=ax.transAxes,
        fontsize=7.5, va="bottom", family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="grey", alpha=0.88),
    )

    ax.set_title(title or "Green-Lagrange strain tensor", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, lw=0.3, alpha=0.35)
    return ax


# ---------------------------------------------------------------------------
# Plan E — match card (2 × 2 combined summary figure)
# ---------------------------------------------------------------------------

def plot_match_card(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    l1: CoincidenceResult,
    match: MatchResult,
    *,
    save_path: Optional[str] = None,
    dpi: int = 150,
) -> Figure:
    """
    Generate a 2×2 "match card" figure combining all four visualisation panels:

        ┌─────────────────────┬─────────────────────┐
        │  Φ(θ) curve         │  Real-space overlay  │
        │  (Plan B)           │  (Plan A)            │
        ├─────────────────────┼─────────────────────┤
        │  LEED pattern       │  Strain ellipse      │
        │  (Plan C)           │  (Plan D)            │
        └─────────────────────┴─────────────────────┘

    Parameters
    ----------
    lat_sub, lat_mof : substrate and MOF ``Lattice2D`` objects
    l1               : ``CoincidenceResult`` for the Φ curve panel
    match            : ``MatchResult`` to visualise in all four panels
    save_path        : if provided, save the figure here (PNG / SVG / PDF)
    dpi              : pixel density when saving to a raster format

    Returns
    -------
    ``matplotlib.figure.Figure`` — caller is responsible for ``plt.close(fig)``
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    fig.suptitle(
        f"{lat_mof.label or 'Overlayer'}  on  "
        f"{lat_sub.label or 'Substrate'}\n"
        f"θ = {match.theta_deg:.2f}°    "
        f"η = {match.eta:.5f}    "
        f"Area = {match.area:.0f} Å²",
        fontsize=11,
        fontweight="bold",
    )

    # Top-left: Phi curve with the chosen angle highlighted
    plot_phi_curve(l1, ax=axes[0, 0],
                   title="Level-1: Coincidence Function Φ(θ)")
    axes[0, 0].axvline(
        match.theta_deg, color=_C_CELL, lw=1.8, zorder=5,
        label=f"Selected  θ = {match.theta_deg:.1f}°",
    )
    axes[0, 0].legend(fontsize=7, loc="upper right")

    # Top-right: real-space overlay
    plot_lattice_overlay(lat_sub, lat_mof, match,
                         ax=axes[0, 1],
                         title="Level-2: Real-space lattice overlay")

    # Bottom-left: LEED simulation
    plot_leed_pattern(lat_sub, lat_mof, match,
                      ax=axes[1, 0],
                      title="Simulated LEED pattern")

    # Bottom-right: strain ellipse
    plot_strain_ellipse(match,
                        ax=axes[1, 1],
                        title="Green-Lagrange strain tensor")

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

    if save_path is not None:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    return fig


# ---------------------------------------------------------------------------
# PDF report  (two-page: text summary + match-card figure)
# ---------------------------------------------------------------------------

_INVALID_CHARS = re.compile(r'[\\/:*?"<>|\[\]\s]+')


def _safe_stem(text: str) -> str:
    """Return a filename-safe stem (Windows-compatible: no / \\ : * ? \" < > | [ ])."""
    return _INVALID_CHARS.sub('_', text).strip('_') or 'report'


_PLACEHOLDER = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{title}</title>
<style>
  body{{font-family:Arial,Helvetica,sans-serif;margin:40px;color:#222;background:#fafafa}}
  h1{{font-size:1.6em;border-bottom:2px solid #1f77b4;padding-bottom:6px}}
  h2{{font-size:1.15em;color:#1f77b4;margin-top:32px}}
  table{{border-collapse:collapse;width:100%;margin-top:8px;font-size:0.9em}}
  th{{background:#1f77b4;color:#fff;padding:6px 10px;text-align:left}}
  td{{padding:5px 10px;border-bottom:1px solid #ddd}}
  tr:nth-child(even){{background:#f0f4fa}}
  .best{{background:#d4edda !important;font-weight:bold}}
  .card{{margin-top:18px;text-align:center}}
  .card img{{max-width:100%;border:1px solid #ccc;border-radius:4px;box-shadow:2px 2px 6px #bbb}}
  .footer{{margin-top:40px;font-size:0.75em;color:#888}}
</style>
</head>
<body>
<h1>{title}</h1>
<p><b>Generated:</b> {timestamp}</p>

<h2>Lattice parameters</h2>
<table>
<tr><th></th><th>a (Å)</th><th>b (Å)</th><th>γ (°)</th></tr>
<tr><td><b>Substrate</b> — {sub_label}</td>
    <td>{sub_a:.3f}</td><td>{sub_b:.3f}</td><td>{sub_g:.1f}</td></tr>
<tr><td><b>Overlayer</b> — {mof_label}</td>
    <td>{mof_a:.3f}</td><td>{mof_b:.3f}</td><td>{mof_g:.1f}</td></tr>
</table>

<h2>Top matches (ranked by strain index η)</h2>
{match_table}

<h2>Best match — visualisation</h2>
<div class="card"><img src="data:image/png;base64,{img_b64}" alt="match card"/></div>

<div class="footer">Generated by miqrophi · {timestamp}</div>
</body>
</html>"""


def generate_pdf_report(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    l1: CoincidenceResult,
    matches: "list[MatchResult]",
    *,
    title: str = "Epitaxy Matching Report",
    save_path: str = "output/report.pdf",
    dpi: int = 150,
    top_n: int = 10,
) -> str:
    """
    Generate a two-page PDF report using matplotlib's built-in PDF backend
    (no extra dependencies required).

    * **Page 1** — text summary: lattice parameters + ranked match table.
    * **Page 2** — match-card figure (Φ curve, real-space overlay, LEED,
      strain ellipse).

    Parameters
    ----------
    lat_sub, lat_mof : substrate and MOF ``Lattice2D`` objects
    l1               : ``CoincidenceResult`` (used for the Φ(θ) panel on page 2)
    matches          : list of ``MatchResult`` objects, pre-sorted by η
    title            : report title shown on page 1
    save_path        : output path for the ``.pdf`` file
    dpi              : figure resolution (affects rasterised elements)
    top_n            : number of matches to include in the summary table

    Returns
    -------
    Absolute path to the saved PDF file.
    """
    from datetime import datetime

    from matplotlib.backends.backend_pdf import PdfPages

    if not matches:
        raise ValueError("matches list is empty — nothing to report")

    best = matches[0]
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)

    with PdfPages(save_path) as pdf:
        # ── Page 1: text summary ──────────────────────────────────────────
        fig = plt.figure(figsize=(11, 8.5))
        fig.patch.set_facecolor("white")

        fig.text(0.5, 0.965, title,
                 ha="center", va="top", fontsize=15, fontweight="bold")
        fig.text(0.5, 0.933,
                 f"Generated: {timestamp}    |    miqrophi epitaxy matcher",
                 ha="center", va="top", fontsize=9, color="#666")

        # Lattice-parameter table
        ax1 = fig.add_axes([0.06, 0.775, 0.88, 0.135])
        ax1.axis("off")
        ax1.set_title("Lattice Parameters", loc="left",
                      fontsize=11, fontweight="bold", color="#1f77b4", pad=5)
        lat_rows = [
            [f"Substrate  —  {lat_sub.label or 'Substrate'}",
             f"{lat_sub.a:.3f}", f"{lat_sub.b:.3f}",
             f"{lat_sub.gamma_deg:.1f}°"],
            [f"Overlayer  —  {lat_mof.label or 'Overlayer'}",
             f"{lat_mof.a:.3f}", f"{lat_mof.b:.3f}",
             f"{lat_mof.gamma_deg:.1f}°"],
        ]
        tbl1 = ax1.table(
            cellText=lat_rows,
            colLabels=["", "a (Å)", "b (Å)", "γ"],
            bbox=[0, 0, 1, 1],
            cellLoc="center",
        )
        tbl1.auto_set_font_size(False)
        tbl1.set_fontsize(10)
        for (row, col), cell in tbl1.get_celld().items():
            if row == 0:
                cell.set_facecolor("#1f77b4")
                cell.set_text_props(color="white", fontweight="bold")
            if col == 0 and row > 0:
                cell.set_text_props(ha="left")

        # Match table
        n_show = min(top_n, len(matches))
        tbl_height = min(0.60, (n_show + 1) * 0.048 + 0.04)
        ax2 = fig.add_axes([0.06, 0.755 - tbl_height, 0.88, tbl_height])
        ax2.axis("off")
        ax2.set_title("Top Matches — ranked by strain index η", loc="left",
                      fontsize=11, fontweight="bold", color="#1f77b4", pad=5)
        match_rows = []
        for i, m in enumerate(matches[:n_show]):
            rank_str = f"\u2605 {i + 1}" if i == 0 else str(i + 1)
            match_rows.append([
                rank_str,
                f"{m.theta_deg:.2f}",
                f"{m.eta:.5f}",
                f"{m.strain[0, 0]:+.4f}",
                f"{m.strain[1, 1]:+.4f}",
                f"{m.strain[0, 1]:+.4f}",
                f"{m.area:.0f}",
            ])
        tbl2 = ax2.table(
            cellText=match_rows,
            colLabels=["Rank", "θ (°)", "η",
                       "ε₁₁", "ε₂₂", "ε₁₂", "Area (Å²)"],
            bbox=[0, 0, 1, 1],
            cellLoc="center",
        )
        tbl2.auto_set_font_size(False)
        tbl2.set_fontsize(9)
        for (row, col), cell in tbl2.get_celld().items():
            if row == 0:
                cell.set_facecolor("#1f77b4")
                cell.set_text_props(color="white", fontweight="bold")
            elif row == 1:
                cell.set_facecolor("#d4edda")

        fig.text(0.5, 0.02,
                 f"Generated by miqrophi  ·  {timestamp}",
                 ha="center", va="bottom", fontsize=8, color="#aaa")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── Page 2: match-card figure ─────────────────────────────────────
        fig_card = plot_match_card(lat_sub, lat_mof, l1, best,
                                   save_path=None, dpi=dpi)
        pdf.savefig(fig_card, bbox_inches="tight")
        plt.close(fig_card)

    return os.path.abspath(save_path)


# ---------------------------------------------------------------------------
# Algorithm animation — coincidence search
# ---------------------------------------------------------------------------

def animate_coincidence_search(  # pragma: no cover
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    *,
    G_cutoff: float = 8.0,
    sigma: float = 0.3,
    n_frames: int = 360,
    interval: int = 25,
    save_path: Optional[str] = None,
    fps: int = 30,
    dpi: int = 120,
) -> "FuncAnimation":
    """
    Animate the coincidence search: show the MOF reciprocal lattice rotating
    over the substrate lattice while the Φ(θ) curve is drawn in real time.

    Each frame corresponds to one rotation angle θ stepping from 0° to 360°.

    Left panel  — reciprocal space: substrate points (blue, fixed) and MOF
                  points (red, rotating).  Point size encodes the
                  Debye-Waller weight; overlapping pairs flash green.
    Right panel — Φ(θ) curve accumulated frame-by-frame; detected peaks are
                  annotated with an orange vertical marker.

    Parameters
    ----------
    lat_sub, lat_mof : substrate and MOF ``Lattice2D`` objects
    G_cutoff         : reciprocal truncation radius (Å⁻¹)
    sigma            : Gaussian width for coincidence score (Å⁻¹)
    n_frames         : number of animation frames (= angular resolution)
    interval         : delay between frames in milliseconds
    save_path        : if given, save the animation to this path (.gif or .mp4)
    fps              : frames per second when saving
    dpi              : resolution when saving
    """
    import matplotlib
    import matplotlib.animation as animation
    from scipy.signal import find_peaks

    # ── Pre-compute the full Φ(θ) curve (needed for peak detection) ──────
    from .coincidence import compute as _coincidence_compute
    full_result = _coincidence_compute(
        lat_sub, lat_mof,
        G_cutoff=G_cutoff, sigma=sigma, K=n_frames,
    )
    theta_grid = full_result.theta_grid   # (n_frames,) degrees
    phi_curve  = full_result.phi_curve    # (n_frames,) normalised

    # Peaks on the full curve (for annotation)
    peak_indices, _ = find_peaks(phi_curve, height=0.15, distance=n_frames // 36)
    peak_thetas     = theta_grid[peak_indices]
    peak_phis       = phi_curve[peak_indices]

    # ── Static reciprocal-lattice point sets ─────────────────────────────
    Gs_all = _recip_points_in_disk(lat_sub.recip_matrix(), G_cutoff)  # (Ns, 2)
    Gm_all = _recip_points_in_disk(lat_mof.recip_matrix(), G_cutoff)  # (Nm, 2)

    # Debye-Waller sizes for scatter
    def _dw_sizes(vecs: np.ndarray, base: float = 80.0) -> np.ndarray:
        norms = np.linalg.norm(vecs, axis=1)
        return base * np.exp(-2.0 * (norms / G_cutoff) ** 2) + 8.0

    s_sub = _dw_sizes(Gs_all)

    # ── Figure layout ─────────────────────────────────────────────────────
    matplotlib.use("Agg")
    fig, (ax_recip, ax_phi) = plt.subplots(
        1, 2, figsize=(12, 5.5),
        gridspec_kw={"width_ratios": [1, 1.3]},
    )
    fig.patch.set_facecolor("#0d1117")
    for ax in (ax_recip, ax_phi):
        ax.set_facecolor("#161b22")
        for spine in ax.spines.values():
            spine.set_color("#30363d")
        ax.tick_params(colors="#8b949e", labelsize=8)
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")

    # Left panel — reciprocal space
    lim = G_cutoff * 1.05
    ax_recip.set_xlim(-lim, lim)
    ax_recip.set_ylim(-lim, lim)
    ax_recip.set_aspect("equal")
    ax_recip.set_xlabel("$G_x$  (Å⁻¹)", fontsize=9)
    ax_recip.set_ylabel("$G_y$  (Å⁻¹)", fontsize=9)
    ax_recip.set_title(
        f"Reciprocal space  ·  {lat_sub.label or 'substrate'} (blue)  +  "
        f"{lat_mof.label or 'MOF'} (red, rotating)",
        color="#c9d1d9", fontsize=9, pad=6,
    )
    ax_recip.axhline(0, color="#30363d", lw=0.5)
    ax_recip.axvline(0, color="#30363d", lw=0.5)

    ax_recip.scatter(
        Gs_all[:, 0], Gs_all[:, 1],
        s=s_sub, c=_C_SUB, alpha=0.75, zorder=3, label="Substrate",
    )
    sc_mof = ax_recip.scatter(
        Gm_all[:, 0], Gm_all[:, 1],
        s=_dw_sizes(Gm_all), c=_C_MOF, alpha=0.75, zorder=4, label="MOF",
    )
    sc_coin = ax_recip.scatter(
        [], [], s=120, c=_C_COIN, marker="*", zorder=5, label="Coincident",
    )
    ax_recip.legend(
        loc="upper right", fontsize=7,
        facecolor="#161b22", edgecolor="#30363d", labelcolor="#c9d1d9",
    )
    theta_label = ax_recip.text(
        0.03, 0.04, "", transform=ax_recip.transAxes,
        color="#f0e68c", fontsize=10, fontweight="bold",
    )

    # Right panel — Φ(θ) curve
    ax_phi.set_xlim(0, 360)
    ax_phi.set_ylim(-0.05, 1.12)
    ax_phi.set_xlabel("θ  (degrees)", fontsize=9)
    ax_phi.set_ylabel("Φ(θ)  (normalised)", fontsize=9)
    ax_phi.set_title(
        "Coincidence function  Φ(θ)  —  accumulated in real time",
        color="#c9d1d9", fontsize=9, pad=6,
    )
    ax_phi.axhline(0, color="#30363d", lw=0.5)

    line_phi, = ax_phi.plot([], [], color=_C_PHI, lw=1.5, zorder=3)
    vline = ax_phi.axvline(x=0, color="#f0e68c", lw=1.0, ls="--", alpha=0.7)
    peak_markers = []   # will grow as peaks come into view

    fig.tight_layout(pad=2.0)

    # ── Animation update function ──────────────────────────────────────────
    def _update(frame_idx: int):
        theta_deg = theta_grid[frame_idx]
        theta_rad = np.radians(theta_deg)

        # Rotate MOF reciprocal lattice
        R = np.array([[ np.cos(theta_rad), -np.sin(theta_rad)],
                      [ np.sin(theta_rad),  np.cos(theta_rad)]])
        Gm_rot = Gm_all @ R.T
        sc_mof.set_offsets(Gm_rot)
        sc_mof.set_sizes(_dw_sizes(Gm_rot))

        # Highlight coincident pairs: |Gs - Gm_rot| < 2*sigma
        dists = np.linalg.norm(
            Gs_all[:, None, :] - Gm_rot[None, :, :], axis=2
        )  # (Ns, Nm)
        coincident_mask = dists < 2.0 * sigma
        coin_pts_sub = Gs_all[coincident_mask.any(axis=1)]
        coin_pts_mof = Gm_rot[coincident_mask.any(axis=0)]
        coin_all = np.vstack([coin_pts_sub, coin_pts_mof]) if (
            len(coin_pts_sub) and len(coin_pts_mof)
        ) else np.empty((0, 2))
        sc_coin.set_offsets(coin_all if len(coin_all) else np.empty((0, 2)))

        # Theta label
        theta_label.set_text(f"θ = {theta_deg:.1f}°")

        # Accumulate Φ curve up to current frame
        line_phi.set_data(theta_grid[:frame_idx + 1], phi_curve[:frame_idx + 1])

        # Moving cursor
        vline.set_xdata([theta_deg, theta_deg])

        # Add peak markers as they scroll into view
        for i, (pt, pp) in enumerate(zip(peak_thetas, peak_phis)):
            if pt <= theta_deg and len(peak_markers) <= i:
                vl = ax_phi.axvline(x=pt, color=_C_PEAK, lw=1.2, ls=":", alpha=0.9)
                txt = ax_phi.text(
                    pt + 2, pp + 0.04,
                    f"{pt:.1f}°",
                    color=_C_PEAK, fontsize=7, va="bottom",
                )
                peak_markers.append((vl, txt))

        return sc_mof, sc_coin, line_phi, vline, theta_label

    anim = animation.FuncAnimation(
        fig, _update,
        frames=n_frames,
        interval=interval,
        blit=False,
    )

    if save_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        suffix = os.path.splitext(save_path)[1].lower()
        if suffix == ".gif":
            writer = animation.PillowWriter(fps=fps)
        else:
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
        anim.save(save_path, writer=writer, dpi=dpi)

    return anim
