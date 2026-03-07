"""Demonstration entry point for the miqrocal epitaxy matching pipeline.

Runs three representative lattice-matching cases, prints the top tabular
results, and writes PNG/PDF outputs to the ``output/`` directory.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from miqrophi import (
    EpitaxyMatcher,
    Lattice2D,
    MatcherConfig,
    SUBSTRATE_DB,
    animate_coincidence_search,
    generate_pdf_report,
    coincidence,
    supercell,
    plot_match_card,
)

_SAFE = re.compile(r'[\\/:*?"<>|\[\]\s]+')


def _safe_stem(label: str) -> str:
    return _SAFE.sub("_", label).strip("_")


def _collect_matches(
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    cfg: MatcherConfig,
) -> tuple[coincidence.CoincidenceResult, list[supercell.MatchResult]]:
    l1 = coincidence.compute(
        lat_sub,
        lat_mof,
        G_cutoff=cfg.G_cutoff,
        sigma=cfg.sigma,
    )
    matches: list[supercell.MatchResult] = []
    for theta_deg in l1.theta_peaks[: cfg.top_theta]:
        matches.extend(
            supercell.find_matches(
                lat_sub,
                lat_mof,
                theta_deg=theta_deg,
                lambda_values=cfg.lambda_values,
                eta_tol=cfg.eta_tol,
            )
        )
    matches.sort(key=lambda match: match.eta)
    return l1, matches


def _print_results(
    label: str,
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    matcher: EpitaxyMatcher,
) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {label}")
    df = matcher.run(lat_sub, lat_mof, verbose=False)
    if df is None or df.empty:
        print("  no matches within configured eta_tol")
        return
    cols = [
        "theta (deg)",
        "eta",
        "eps_11",
        "eps_22",
        "eps_12",
        "area (A2)",
        "L0_feasible",
    ]
    print(df[cols].head(5).to_string(index=False))


def _write_outputs(
    label: str,
    lat_sub: Lattice2D,
    lat_mof: Lattice2D,
    cfg: MatcherConfig,
    output_dir: Path,
) -> None:
    l1, matches = _collect_matches(lat_sub, lat_mof, cfg)
    if not l1.theta_peaks:
        print(f"  [skip] no Level-1 peaks for '{label}'")
        return
    if not matches:
        print(f"  [skip] no Level-2 matches for '{label}'")
        return

    safe_stem = _safe_stem(label)
    png_path = output_dir / f"{safe_stem}.png"
    pdf_path = output_dir / f"{safe_stem}_report.pdf"

    fig = plot_match_card(lat_sub, lat_mof, l1, matches[0], save_path=str(png_path), dpi=150)
    plt.close(fig)
    saved_pdf = generate_pdf_report(
        lat_sub,
        lat_mof,
        l1,
        matches,
        title=label,
        save_path=str(pdf_path),
    )
    print(f"  [png] {png_path}")
    print(f"  [pdf] {saved_pdf}")

    gif_path = output_dir / f"{safe_stem}_animation.gif"
    animate_coincidence_search(
        lat_sub, lat_mof,
        G_cutoff=cfg.G_cutoff,
        sigma=cfg.sigma,
        n_frames=360,
        interval=25,
        save_path=str(gif_path),
        fps=30,
        dpi=100,
    )
    print(f"  [gif] {gif_path}")


def main() -> None:
    cfg = MatcherConfig(sigma=0.4, eta_tol=0.04)
    matcher = EpitaxyMatcher(cfg)
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    hkust1 = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")
    hex_mof = Lattice2D(4.92, 4.92, 120.0, "Hex-MOF(2xa)")

    cases = [
        (
            "HKUST-1(010) / Au(111)  [expected: L0 forbidden]",
            SUBSTRATE_DB["Au_111"],
            hkust1,
        ),
        (
            "HKUST-1(010) / Au(100)  [expected: commensurate]",
            SUBSTRATE_DB["Au_100"],
            hkust1,
        ),
        (
            "Hexagonal MOF / Au(111)  [expected: multiple rotation domains]",
            SUBSTRATE_DB["Au_111"],
            hex_mof,
        ),
    ]

    print("\n" + "-" * 60)
    print("  miqrocal demo")
    for label, lat_sub, lat_mof in cases:
        _print_results(label, lat_sub, lat_mof, matcher)

    print("\n" + "-" * 60)
    print(f"  Writing figures and reports to {output_dir.resolve()}")
    for label, lat_sub, lat_mof in cases:
        _write_outputs(label, lat_sub, lat_mof, cfg, output_dir)


if __name__ == "__main__":
    main()
