"""Top-level epitaxy matching pipeline (Levels 0 -> 1 -> 2)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from .lattice import Lattice2D
from . import level0, level1, level2


@dataclass
class MatcherConfig:
    """Configuration parameters for EpitaxyMatcher."""

    G_cutoff:      float       = 8.0   # reciprocal-space truncation radius (Angstrom^-1)
    sigma:         float       = 0.3   # Phi(theta) Gaussian peak width (Angstrom^-1)
    eta_tol:       float       = 0.05  # Green-Lagrange strain Frobenius norm tolerance
    top_theta:     int         = 8     # max Level-1 peak angles forwarded to Level 2
    lambda_values: list[float] = field(
        default_factory=lambda: [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
    )


class EpitaxyMatcher:
    """Three-level hierarchical epitaxy matching pipeline."""

    def __init__(self, config: Optional[MatcherConfig] = None) -> None:
        self.config = config or MatcherConfig()

    def run(
        self,
        lat_sub: Lattice2D,
        lat_mof: Lattice2D,
        verbose: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Execute the full Level 0 -> 1 -> 2 pipeline.

        Returns
        -------
        pd.DataFrame of matches sorted by eta, or None if no match satisfies
        the strain tolerance.

        Column description
        ------------------
        theta (deg)  : rotation angle of the MOF supercell relative to substrate
        M, N         : 2x2 integer supercell matrices (lists)
        eta          : Frobenius norm of the Green-Lagrange strain tensor
        eps_11/22/12 : individual strain tensor components
        area (A2)    : substrate supercell area in Angstrom^2
        L0_feasible  : whether Level 0 judges exact commensurability algebraically
                       feasible
        """
        cfg  = self.config
        _log = print if verbose else (lambda *a, **kw: None)

        _log(f"\n{'─' * 60}")
        _log(f"  {lat_mof.label or 'overlayer'}  on  {lat_sub.label or 'substrate'}")
        _log(f"{'─' * 60}")

        # ── Level 0 ────────────────────────────────────────────────
        l0 = level0.check(lat_sub, lat_mof)
        _log(f"  [L0] delta_sub={l0.delta_sub:.4f}  delta_mof={l0.delta_mof:.4f}"
             f"  ratio={l0.ratio:.4f}  ->  {l0.message}")

        # ── Level 1 ────────────────────────────────────────────────
        l1 = level1.compute(
            lat_sub, lat_mof,
            G_cutoff=cfg.G_cutoff,
            sigma=cfg.sigma,
        )
        peak_summary = (
            f"  top: theta={l1.theta_peaks[0]:.1f} deg (Phi={l1.phi_peaks[0]:.3f})"
            if l1.theta_peaks else ""
        )
        _log(f"  [L1] {len(l1.theta_peaks)} peaks found{peak_summary}")

        if not l1.theta_peaks:
            _log("  -> no coincidence peaks; epitaxial growth not feasible")
            return None

        # ── Level 2 ────────────────────────────────────────────────
        all_matches: list[level2.MatchResult] = []
        for th in l1.theta_peaks[: cfg.top_theta]:
            matches = level2.find_matches(
                lat_sub, lat_mof,
                theta_deg=th,
                lambda_values=cfg.lambda_values,
                eta_tol=cfg.eta_tol,
            )
            all_matches.extend(matches)
            if matches and verbose:
                b = matches[0]
                _log(f"  [L2] theta={th:7.2f} deg  eta={b.eta:.5f}"
                     f"  eps=({b.strain[0,0]:+.4f}, {b.strain[1,1]:+.4f},"
                     f" {b.strain[0,1]:+.4f})  area={b.area:.1f} A^2")

        if not all_matches:
            note = "  (consistent with L0 symmetry prohibition)" if not l0.feasible else ""
            _log(f"\n  -> no matches with eta < {cfg.eta_tol}{note}")
            return None

        rows = [
            {
                "theta (deg)": round(m.theta_deg, 2),
                "M":           m.M.tolist(),
                "N":           m.N.tolist(),
                "eta":         round(m.eta, 6),
                "eps_11":      round(m.strain[0, 0], 5),
                "eps_22":      round(m.strain[1, 1], 5),
                "eps_12":      round(m.strain[0, 1], 5),
                "area (A2)":   round(m.area, 1),
                "L0_feasible": l0.feasible,
            }
            for m in all_matches
        ]
        df = (
            pd.DataFrame(rows)
            .sort_values("eta")
            .drop_duplicates(subset=["eta", "area (A2)"])
            .reset_index(drop=True)
        )

        tag = "commensurate" if l0.feasible else "approximate only (L0 forbidden)"
        _log(f"\n  Best match: eta={df['eta'].iloc[0]:.5f}  [{tag}]")
        return df
