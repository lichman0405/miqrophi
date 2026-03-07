"""
miqrocal.cli — Command-line interface for the miqrocal package.

Sub-commands
------------
substrates   List all built-in substrate keys and their lattice parameters.
run          Run matching for specific CIF file(s) against chosen substrates.
batch        Screen a glob of CIF files against multiple substrates.

Usage
-----
    miqrocal substrates
    miqrocal run FILE.cif [FILE2.cif ...] [options]
    miqrocal batch "examples/*.cif"  [options]

Common options
--------------
    --substrates KEY ...    substrate keys from SUBSTRATE_DB (default: all)
    --n-faces N             BFDH faces to consider per CIF (default: 3)
    --output-dir DIR        directory for PNG/PDF/CSV output  (default: output)
    --tag TAG               run-tag appended to output directory name
    --n-jobs N              parallel workers: 1=serial, -1=all CPUs (default: 1)
    --outputs csv,pdf,png   comma-separated output types (default: csv,pdf,png)
    --eta-tol FLOAT         Green-Lagrange strain tolerance (default: 0.05)
    --sigma FLOAT           Phi(theta) peak width in Angstrom^-1 (default: 0.3)
    --quiet                 suppress per-pair progress output
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="miqrocal",
        description="2D epitaxial lattice matching from CIF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # ── substrates ────────────────────────────────────────────────────────
    sub.add_parser(
        "substrates",
        help="List all built-in substrate keys and lattice parameters.",
    )

    # ── shared run/batch options ──────────────────────────────────────────
    _common = argparse.ArgumentParser(add_help=False)
    _common.add_argument(
        "--substrates", nargs="+", metavar="KEY",
        help="Substrate keys from SUBSTRATE_DB (default: all built-in).",
    )
    _common.add_argument("--n-faces",    type=int,   default=3,       metavar="N")
    _common.add_argument("--output-dir", type=str,   default="output", metavar="DIR")
    _common.add_argument("--tag",        type=str,   default="",       metavar="TAG")
    _common.add_argument("--n-jobs",     type=int,   default=1,        metavar="N")
    _common.add_argument(
        "--outputs", type=str, default="csv,pdf,png",
        help="Comma-separated subset of {csv,pdf,png} (default: csv,pdf,png).",
    )
    _common.add_argument("--eta-tol", type=float, default=0.05, metavar="FLOAT")
    _common.add_argument("--sigma",   type=float, default=0.3,  metavar="FLOAT")
    _common.add_argument("--quiet",   action="store_true")

    # ── run ───────────────────────────────────────────────────────────────
    p_run = sub.add_parser(
        "run",
        parents=[_common],
        help="Run matching for one or more CIF files.",
    )
    p_run.add_argument(
        "cif_files", nargs="+", metavar="FILE.cif",
        help="One or more CIF file paths.",
    )

    # ── batch ─────────────────────────────────────────────────────────────
    p_batch = sub.add_parser(
        "batch",
        parents=[_common],
        help="Screen a glob pattern of CIF files.",
    )
    p_batch.add_argument(
        "glob", metavar="GLOB",
        help='Glob pattern, e.g.  "examples/*.cif"  (quote to avoid shell expansion).',
    )

    return parser


def _resolve_substrates(keys: list[str] | None) -> dict:
    from miqrophi import SUBSTRATE_DB
    if keys is None:
        return dict(SUBSTRATE_DB)
    missing = [k for k in keys if k not in SUBSTRATE_DB]
    if missing:
        sys.exit(
            f"[miqrocal] Unknown substrate key(s): {missing}\n"
            f"Run  miqrocal substrates  to list valid keys."
        )
    return {k: SUBSTRATE_DB[k] for k in keys}


def _parse_outputs(raw: str) -> set[str]:
    parts = {p.strip().lower() for p in raw.split(",")}
    valid = {"csv", "pdf", "png"}
    bad   = parts - valid
    if bad:
        sys.exit(f"[miqrocal] Unknown output type(s): {bad}; valid: {valid}")
    return parts


def _print_summary(df: pd.DataFrame) -> None:
    matched = df[df["match_found"]] if "match_found" in df.columns else df
    total   = len(df)
    print(f"\n  {len(matched)} / {total} pairs matched.")
    if not matched.empty:
        cols = [c for c in ["mof_name", "hkl", "substrate",
                             "theta_deg", "eta", "L0_feasible"] if c in matched.columns]
        print(matched[cols].sort_values("eta").head(10).to_string(index=False))


# ---------------------------------------------------------------------------
# sub-command handlers
# ---------------------------------------------------------------------------

def _cmd_substrates(_args: argparse.Namespace) -> int:
    from miqrophi import SUBSTRATE_DB
    print(f"\n{'Key':<20} {'a (Å)':>8} {'b (Å)':>8} {'γ (°)':>8}")
    print("─" * 48)
    for key, lat in SUBSTRATE_DB.items():
        print(f"{key:<20} {lat.a:>8.3f} {lat.b:>8.3f} {lat.gamma_deg:>8.1f}")
    print()
    return 0


def _cmd_run_or_batch(args: argparse.Namespace, glob_or_files) -> int:
    from miqrophi import BatchConfig, MatcherConfig, batch_run

    subs    = _resolve_substrates(args.substrates)
    outputs = _parse_outputs(args.outputs)
    mcfg    = MatcherConfig(eta_tol=args.eta_tol, sigma=args.sigma)
    cfg     = BatchConfig(
        substrates   = list(subs.keys()),
        n_faces      = args.n_faces,
        output_dir   = args.output_dir,
        run_tag      = args.tag,
        n_jobs       = args.n_jobs,
        outputs      = outputs,
        matcher_cfg  = mcfg,
        verbose      = not args.quiet,
    )

    df = batch_run(glob_or_files, config=cfg)
    _print_summary(df)
    return 0


# ---------------------------------------------------------------------------
# main entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    """
    CLI entry point registered as ``miqrocal`` by pyproject.toml.

    Returns an integer exit code so it can be tested without calling
    :func:`sys.exit` directly.
    """
    parser = _build_parser()
    args   = parser.parse_args(argv)

    if args.command == "substrates":
        return _cmd_substrates(args)

    if args.command == "run":
        return _cmd_run_or_batch(args, args.cif_files)

    if args.command == "batch":
        return _cmd_run_or_batch(args, args.glob)

    # no sub-command given
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
