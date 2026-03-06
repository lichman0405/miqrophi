"""
Batch epitaxy matching: screen many CIF files against multiple substrates.

Public API
----------
BatchConfig : configuration dataclass
batch_run   : run the full pipeline on a list / glob of CIF files,
              write per-pair outputs (PNG, PDF) and a summary CSV,
              return a summary DataFrame.

Output directory structure
--------------------------
Each call to :func:`batch_run` creates a **timestamped run directory** inside
the configured ``output_dir``::

    {output_dir}/
      run_{YYYY-MM-DD_HH-MM-SS}[_{tag}]/   ← one directory per run
        summary.csv                          ← all (MOF × substrate) results
        {mof}_{hkl}__{substrate}/            ← one sub-directory per pair
          match_card.png
          report.pdf
        ...

The sub-directory name is derived entirely from the CIF filename, the surface
Miller indices, and the substrate key — it contains **no spaces, colons, or
other Windows-unsafe characters** and can be used as a stable identifier to
trace results back to their inputs.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")   # non-interactive; safe to import before plt
import matplotlib.pyplot as plt
import pandas as pd

from .lattice import SUBSTRATE_DB
from .matcher import EpitaxyMatcher, MatcherConfig
from .cif_parser import best_surface_lattice, read_cell
from . import level0 as _level0
from . import level1 as _level1
from . import level2 as _level2
from .visualize import plot_match_card, generate_pdf_report

_INVALID = re.compile(r'[\\/:*?"<>|\[\]\s]+')


def _safe(text: str) -> str:
    """Return a Windows-safe, traceable filename stem."""
    return _INVALID.sub("_", text).strip("_") or "unnamed"


# ---------------------------------------------------------------------------
# BatchConfig
# ---------------------------------------------------------------------------

@dataclass
class BatchConfig:
    """
    Configuration for :func:`batch_run`.

    Parameters
    ----------
    substrates : ``SUBSTRATE_DB`` keys to screen against.
        ``None`` (default) uses **all** entries in the database.
    n_faces : number of BFDH-ranked faces to try per CIF (default 1 = top
        face only; increase to 3–5 to explore alternative orientations).
    outputs : controls which files are written for each matched pair.
        Accepted values: ``"csv"``, ``"pdf"``, ``"png"``.
        Pass an empty set to suppress all file output (in-memory only).
    output_dir : parent directory for all run output (default ``"output"``).
    run_tag : optional human-readable suffix appended to the run directory
        name, e.g. ``"mof_screen_v1"``.  May contain letters, digits, and
        hyphens; other characters are sanitised automatically.
    matcher_cfg : :class:`MatcherConfig` applied to every pair.
        ``None`` uses library defaults (σ=0.3, η_tol=0.05).
    top_n_report : number of match rows included in each PDF report (default 10).
    verbose : print per-pair progress to stdout (default True).
    """

    substrates:   Optional[list[str]]     = None
    n_faces:      int                      = 1
    outputs:      set[str]                 = field(
        default_factory=lambda: {"csv", "pdf", "png"}
    )
    output_dir:   str                      = "output"
    run_tag:      str                      = ""
    matcher_cfg:  Optional[MatcherConfig]  = None
    top_n_report: int                      = 10
    verbose:      bool                     = True


# ---------------------------------------------------------------------------
# batch_run
# ---------------------------------------------------------------------------

def batch_run(
    cif_paths: "list[str | Path] | str",
    *,
    config: Optional[BatchConfig] = None,
) -> pd.DataFrame:
    """
    Run the three-level epitaxy matching pipeline on multiple CIF files.

    For each CIF the BFDH rule selects up to ``config.n_faces`` surfaces,
    which are then matched against every substrate in ``config.substrates``
    (or the full :data:`~miqrocal.SUBSTRATE_DB` if *None*).

    Parameters
    ----------
    cif_paths : a glob pattern string (e.g. ``"data/*.cif"``) **or** an
        explicit list of CIF file paths.
    config : :class:`BatchConfig`; library defaults used when *None*.

    Returns
    -------
    ``pd.DataFrame`` — one row per *(MOF face × substrate)* pair tried.

    Summary CSV columns
    -------------------
    +---------------+--------------------------------------------------------+
    | Column        | Description                                            |
    +===============+========================================================+
    | cif_file      | basename of the source CIF                             |
    | mof_name      | compound name from CIF (or file stem)                  |
    | hkl           | Miller indices, e.g. ``"(0,0,1)"``                    |
    | d_hkl_A       | interplanar spacing (Å)                                |
    | mof_a/b/gamma | 2D surface lattice parameters                          |
    | substrate     | substrate key from SUBSTRATE_DB                        |
    | match_found   | True if any match with η < η_tol                      |
    | theta_deg     | rotation angle of best match (NaN if no match)         |
    | eta           | strain index η of best match (NaN if no match)         |
    | eps_11/22/12  | Green–Lagrange strain components                       |
    | area_A2       | substrate supercell area (Å²)                          |
    | L0_feasible   | algebraic commensurability flag                        |
    | png_path      | path to match-card PNG relative to ``output_dir``     |
    | pdf_path      | path to PDF report relative to ``output_dir``         |
    +---------------+--------------------------------------------------------+
    """
    cfg  = config or BatchConfig()
    mcfg = cfg.matcher_cfg or MatcherConfig()

    # ── resolve CIF list ────────────────────────────────────────────────────
    if isinstance(cif_paths, str):
        import glob as _glob
        paths = sorted(_glob.glob(cif_paths, recursive=True))
        if not paths:
            raise FileNotFoundError(f"No CIFs matched glob pattern: {cif_paths!r}")
    else:
        paths = [str(p) for p in cif_paths]

    # ── substrate map ────────────────────────────────────────────────────────
    keys = cfg.substrates if cfg.substrates is not None else list(SUBSTRATE_DB.keys())
    subs = {k: SUBSTRATE_DB[k] for k in keys if k in SUBSTRATE_DB}
    if not subs:
        raise ValueError(f"No valid substrate keys found in SUBSTRATE_DB: {keys}")

    # ── timestamped run directory ────────────────────────────────────────────
    ts       = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    dir_name = f"run_{ts}_{_safe(cfg.run_tag)}" if cfg.run_tag else f"run_{ts}"
    run_dir  = Path(cfg.output_dir) / dir_name
    run_dir.mkdir(parents=True, exist_ok=True)

    _log = print if cfg.verbose else (lambda *a, **kw: None)
    _log(f"\n{'═' * 62}")
    _log(f"  miqrocal batch_run")
    _log(f"  {len(paths)} CIF(s)  ×  {len(subs)} substrate(s)"
         f"  ×  up to {cfg.n_faces} face(s)  =  up to"
         f" {len(paths) * len(subs) * cfg.n_faces} pairs")
    _log(f"  Run directory: {run_dir.resolve()}")
    _log(f"  Outputs:       {sorted(cfg.outputs) or ['none (in-memory only)']}")
    _log(f"{'═' * 62}")

    rows: list[dict] = []

    for cif in paths:
        cif_path = Path(cif)
        try:
            cell     = read_cell(cif_path)
            mof_name = str(cell.get("name", cif_path.stem))
        except Exception:
            mof_name = cif_path.stem

        _log(f"\n  ┌─ {cif_path.name}  ({mof_name})")

        try:
            faces = best_surface_lattice(cif_path, n_faces=cfg.n_faces)
        except Exception as exc:
            _log(f"  │  WARN: BFDH failed — {exc}")
            continue

        for hkl, d_hkl, lat_mof in faces:
            hkl_str = f"({','.join(str(x) for x in hkl)})"
            _log(f"  │  face {hkl_str}  d={d_hkl:.2f} Å  "
                 f"a={lat_mof.a:.2f} b={lat_mof.b:.2f} γ={lat_mof.gamma_deg:.1f}°")

            # per-face sub-directory prefix (shared across substrates)
            hkl_compact = "".join(str(x) for x in hkl).replace("-", "m")
            face_stem   = f"{_safe(cif_path.stem)}_{hkl_compact}"

            for sub_key, lat_sub in subs.items():

                # default row — filled in regardless of match
                row: dict = {
                    "cif_file":    cif_path.name,
                    "mof_name":    mof_name,
                    "hkl":         hkl_str,
                    "d_hkl_A":     round(d_hkl, 3),
                    "mof_a":       round(lat_mof.a, 3),
                    "mof_b":       round(lat_mof.b, 3),
                    "mof_gamma":   round(lat_mof.gamma_deg, 2),
                    "substrate":   sub_key,
                    "match_found": False,
                    "theta_deg":   float("nan"),
                    "eta":         float("nan"),
                    "eps_11":      float("nan"),
                    "eps_22":      float("nan"),
                    "eps_12":      float("nan"),
                    "area_A2":     float("nan"),
                    "L0_feasible": False,
                    "png_path":    "",
                    "pdf_path":    "",
                }

                # ── Level 1 + 2 (keep raw objects for visualisation) ────────
                try:
                    l1_res = _level1.compute(
                        lat_sub, lat_mof,
                        G_cutoff=mcfg.G_cutoff,
                        sigma=mcfg.sigma,
                    )
                    all_matches: list[_level2.MatchResult] = []
                    for th in l1_res.theta_peaks[: mcfg.top_theta]:
                        all_matches.extend(
                            _level2.find_matches(
                                lat_sub, lat_mof,
                                theta_deg=th,
                                lambda_values=mcfg.lambda_values,
                                eta_tol=mcfg.eta_tol,
                            )
                        )
                    all_matches.sort(key=lambda m: m.eta)
                except Exception as exc:
                    _log(f"  │    {sub_key:16s}  ERROR: {exc}")
                    rows.append(row)
                    continue

                if not all_matches:
                    _log(f"  │    {sub_key:16s}  — no match  (η_tol={mcfg.eta_tol})")
                    rows.append(row)
                    continue

                # ── match found ─────────────────────────────────────────────
                best  = all_matches[0]
                l0    = _level0.check(lat_sub, lat_mof)
                row.update({
                    "match_found": True,
                    "theta_deg":   round(best.theta_deg, 2),
                    "eta":         round(best.eta, 6),
                    "eps_11":      round(best.strain[0, 0], 5),
                    "eps_22":      round(best.strain[1, 1], 5),
                    "eps_12":      round(best.strain[0, 1], 5),
                    "area_A2":     round(best.area, 1),
                    "L0_feasible": l0.feasible,
                })

                # per-pair sub-directory: {mof_stem}_{hkl}__{substrate_key}
                pair_dir = run_dir / f"{face_stem}__{sub_key}"
                if cfg.outputs & {"png", "pdf"}:
                    pair_dir.mkdir(exist_ok=True)

                if "png" in cfg.outputs:
                    png_path = pair_dir / "match_card.png"
                    try:
                        fig = plot_match_card(
                            lat_sub, lat_mof, l1_res, best,
                            save_path=str(png_path), dpi=150,
                        )
                        plt.close(fig)
                        row["png_path"] = str(png_path.relative_to(
                            Path(cfg.output_dir).resolve()
                        ))
                    except Exception as exc:
                        _log(f"  │    {sub_key:16s}  WARN PNG: {exc}")

                if "pdf" in cfg.outputs:
                    pdf_path = pair_dir / "report.pdf"
                    try:
                        generate_pdf_report(
                            lat_sub, lat_mof, l1_res, all_matches,
                            title=f"{mof_name} {hkl_str}  /  {sub_key}",
                            save_path=str(pdf_path),
                            top_n=cfg.top_n_report,
                        )
                        row["pdf_path"] = str(pdf_path.relative_to(
                            Path(cfg.output_dir).resolve()
                        ))
                    except Exception as exc:
                        _log(f"  │    {sub_key:16s}  WARN PDF: {exc}")

                _log(
                    f"  │    {sub_key:16s}  η={best.eta:.5f}"
                    f"  θ={best.theta_deg:.1f}°"
                    f"  L0={'✓' if l0.feasible else '✗'}"
                )
                rows.append(row)

        _log(f"  └─ done")

    # ── summary ─────────────────────────────────────────────────────────────
    df = pd.DataFrame(rows) if rows else pd.DataFrame()

    if "csv" in cfg.outputs:
        csv_path = run_dir / "summary.csv"
        df.to_csv(csv_path, index=False)
        _log(f"\n  Summary CSV  → {csv_path.resolve()}")

    n_match = int(df["match_found"].sum()) if not df.empty and "match_found" in df else 0
    _log(f"\n{'═' * 62}")
    _log(f"  Done.  {len(df)} pairs evaluated,  {n_match} match(es) found.")
    _log(f"  Run directory: {run_dir.resolve()}")
    _log(f"{'═' * 62}\n")

    return df
