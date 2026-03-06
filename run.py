"""
Demonstration entry point for the miqrocal epitaxy matching pipeline.

Runs the three-level matching pipeline on three test cases, then generates
a match-card visualisation (PNG) for each commensurate result.

Usage
-----
    python run.py

Outputs are written to the ``output/`` directory (created automatically).
"""

import os
import re

import matplotlib
matplotlib.use("Agg")          # non-interactive backend; no display required
import matplotlib.pyplot as plt

from miqrocal import (
    EpitaxyMatcher, MatcherConfig, Lattice2D, SUBSTRATE_DB,
    plot_match_card, generate_pdf_report,
)
from miqrocal import level1, level2


def _run_text(
    label: str,
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    matcher: EpitaxyMatcher,
) -> None:
    """Print tabular results to stdout."""
    print(f"\n{'=' * 60}")
    print(f"  {label}")
    df = matcher.run(lat_sub, lat_mof)
    if df is not None:
        cols = ["theta (deg)", "eta", "eps_11", "eps_22", "eps_12",
                "area (A2)", "L0_feasible"]
        print(df[cols].head(5).to_string(index=False))


def _run_visual(
    label: str,
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    cfg: MatcherConfig,
    out_dir: str = "output",
) -> None:
    """
    Run Level 1 + Level 2 directly to obtain raw result objects, then
    save a match-card figure for the best match found.
    """
    l1 = level1.compute(
        lat_sub, lat_mof,
        G_cutoff=cfg.G_cutoff,
        sigma=cfg.sigma,
    )
    if not l1.theta_peaks:
        print(f"  [vis] no peaks — skipping figure for '{label}'")
        return

    all_matches: list[level2.MatchResult] = []
    for th in l1.theta_peaks[: cfg.top_theta]:
        all_matches.extend(
            level2.find_matches(
                lat_sub, lat_mof,
                theta_deg=th,
                lambda_values=cfg.lambda_values,
                eta_tol=cfg.eta_tol,
            )
        )

    if not all_matches:
        print(f"  [vis] no matches — skipping figure for '{label}'")
        return

    best = min(all_matches, key=lambda m: m.eta)

    os.makedirs(out_dir, exist_ok=True)
    safe  = re.sub(r'[\\/:*?"<>|\[\]\s]+', '_', label).strip('_')
    path  = os.path.join(out_dir, f"{safe}.png")
    fig   = plot_match_card(lat_sub, lat_mof, l1, best, save_path=path, dpi=150)
    plt.close(fig)
    print(f"  [vis] saved → {path}")


def _run_report(
    label: str,
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    cfg: MatcherConfig,
    out_dir: str = "output",
) -> None:
    """
    Run Level 1 + Level 2, then write a two-page PDF report:
    page 1 is the ranked match table, page 2 is the match-card figure.
    """
    l1 = level1.compute(
        lat_sub, lat_mof,
        G_cutoff=cfg.G_cutoff,
        sigma=cfg.sigma,
    )
    if not l1.theta_peaks:
        print(f"  [report] no peaks — skipping '{label}'")
        return

    all_matches: list[level2.MatchResult] = []
    for th in l1.theta_peaks[: cfg.top_theta]:
        all_matches.extend(
            level2.find_matches(
                lat_sub, lat_mof,
                theta_deg=th,
                lambda_values=cfg.lambda_values,
                eta_tol=cfg.eta_tol,
            )
        )

    if not all_matches:
        print(f"  [report] no matches — skipping '{label}'")
        return

    all_matches.sort(key=lambda m: m.eta)

    os.makedirs(out_dir, exist_ok=True)
    safe = re.sub(r'[\\/:*?"<>|\[\]\s]+', '_', label).strip('_')
    path = os.path.join(out_dir, f"{safe}_report.pdf")
    saved = generate_pdf_report(
        lat_sub, lat_mof, l1, all_matches,
        title=label,
        save_path=path,
    )
    print(f"  [report] saved → {saved}")


def main() -> None:
    cfg     = MatcherConfig(sigma=0.4, eta_tol=0.04)
    matcher = EpitaxyMatcher(cfg)

    hkust1  = Lattice2D(18.62, 18.62,  90.0, "HKUST-1(010)")
    hex_mof = Lattice2D( 4.92,  4.92, 120.0, "Hex-MOF(2xa)")

    cases = [
        ("HKUST-1(010) / Au(111)  [expected: L0 forbidden]",
         SUBSTRATE_DB["Au_111"], hkust1),
        ("HKUST-1(010) / Au(100)  [expected: commensurate]",
         SUBSTRATE_DB["Au_100"], hkust1),
        ("Hexagonal MOF / Au(111)  [expected: multiple rotation domains]",
         SUBSTRATE_DB["Au_111"], hex_mof),
    ]

    # ── Text output ────────────────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("  Text results")
    for label, lat_sub, lat_mof in cases:
        _run_text(label, lat_sub, lat_mof, matcher)

    # ── PNG match-card figures ──────────────────────────────────────────────
    print("\n" + "─" * 60)
    print("  Generating match-card figures  (output/ directory)")
    for label, lat_sub, lat_mof in cases:
        _run_visual(label, lat_sub, lat_mof, cfg)

    # ── HTML reports (figure + text table, self-contained) ─────────────────
    print("\n" + "─" * 60)
    print("  Generating PDF reports  (output/ directory)")
    for label, lat_sub, lat_mof in cases:
        _run_report(label, lat_sub, lat_mof, cfg)


if __name__ == "__main__":
    main()
