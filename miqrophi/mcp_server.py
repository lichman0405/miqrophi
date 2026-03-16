"""
miqrophi MCP Server

Exposes 2D epitaxial lattice matching as MCP tools for AI agents:

  - list_substrates   : enumerate all built-in substrate lattice entries
  - parse_cif_surface : extract 2D surface lattice parameters from a CIF
  - analyze_epitaxy   : run the full Level 0 → 1 → 2 epitaxy pipeline
"""
from __future__ import annotations

import os
import tempfile
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # non-interactive backend; must precede pyplot import
import matplotlib.pyplot as plt  # noqa: E402
from mcp.server.fastmcp import FastMCP  # noqa: E402

from . import coincidence, discriminant, supercell  # noqa: E402
from .cif_parser import best_surface_lattice, read_cell, surface_lattice  # noqa: E402
from .lattice import SUBSTRATE_DB, Lattice2D  # noqa: E402
from .visualize import (  # noqa: E402
    animate_coincidence_search,
    generate_pdf_report,
    plot_match_card,
)

# ---------------------------------------------------------------------------
# Server instance
# ---------------------------------------------------------------------------

mcp = FastMCP(
    "miqrophi",
    instructions=(
        "Tools for 2D epitaxial lattice matching between crystalline "
        "overlayers (MOFs, COFs, molecular films) and single-crystal "
        "substrates. Determines optimal rotation angles, supercell "
        "matrices, and epitaxial strain via a three-level hierarchical "
        "screening pipeline (algebraic pre-filter → reciprocal-space "
        "coincidence → real-space supercell optimisation)."
    ),
)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_cif_path(
    cif_path: Optional[str],
    cif_content: Optional[str],
) -> str:
    """Return a filesystem path to a CIF file, writing a temp file if needed."""
    if cif_path is not None:
        if not os.path.isfile(cif_path):
            raise FileNotFoundError(f"File not found: {cif_path}")
        return cif_path
    if cif_content is not None:
        with tempfile.NamedTemporaryFile(
            suffix=".cif", mode="w", delete=False, encoding="utf-8",
        ) as fh:
            fh.write(cif_content)
            return fh.name
    raise ValueError("Provide either cif_path or cif_content.")


def _lattice_info(lat: Lattice2D) -> dict:
    return {
        "a": round(lat.a, 4),
        "b": round(lat.b, 4),
        "gamma_deg": round(lat.gamma_deg, 2),
        "label": lat.label,
    }


def _resolve_substrate(
    substrate_name: Optional[str],
    substrate_a: Optional[float],
    substrate_b: Optional[float],
    substrate_gamma_deg: Optional[float],
) -> Lattice2D:
    """Build a substrate Lattice2D from either a DB key or explicit params."""
    if substrate_name is not None:
        if substrate_name not in SUBSTRATE_DB:
            raise ValueError(
                f"Unknown substrate '{substrate_name}'. "
                f"Valid keys: {list(SUBSTRATE_DB.keys())}"
            )
        return SUBSTRATE_DB[substrate_name]
    if substrate_a is not None and substrate_b is not None and substrate_gamma_deg is not None:
        return Lattice2D(substrate_a, substrate_b, substrate_gamma_deg, label="custom-substrate")
    raise ValueError(
        "Provide either substrate_name (e.g. 'Au_111') or all three of "
        "substrate_a, substrate_b, substrate_gamma_deg."
    )


def _resolve_overlayer(
    cif_path: Optional[str],
    cif_content: Optional[str],
    hkl: Optional[list[int]],
    overlayer_a: Optional[float],
    overlayer_b: Optional[float],
    overlayer_gamma_deg: Optional[float],
) -> Lattice2D:
    """Build an overlayer Lattice2D from CIF or explicit params."""
    if overlayer_a is not None and overlayer_b is not None and overlayer_gamma_deg is not None:
        return Lattice2D(overlayer_a, overlayer_b, overlayer_gamma_deg, label="custom-overlayer")
    if cif_path is not None or cif_content is not None:
        real_path = _resolve_cif_path(cif_path, cif_content)
        is_temp = cif_content is not None
        try:
            if hkl is not None:
                return surface_lattice(real_path, hkl=tuple(hkl))
            else:
                faces = best_surface_lattice(real_path, n_faces=1)
                if not faces:
                    raise ValueError("No valid BFDH surface found in CIF.")
                return faces[0][2]
        finally:
            if is_temp:
                os.unlink(real_path)
    raise ValueError(
        "Provide CIF input (cif_path or cif_content) or all three of "
        "overlayer_a, overlayer_b, overlayer_gamma_deg."
    )


# ---------------------------------------------------------------------------
# MCP Tools
# ---------------------------------------------------------------------------

@mcp.tool()
def list_substrates() -> dict:
    """
    List all built-in substrate lattice entries available for epitaxy matching.

    Returns a dictionary mapping substrate keys (e.g. 'Au_111', 'Cu_111',
    'HOPG') to their 2D lattice parameters (a, b in Angstrom; gamma in
    degrees) and human-readable labels.
    """
    entries = {}
    for key, lat in SUBSTRATE_DB.items():
        entries[key] = _lattice_info(lat)
    return {"substrates": entries, "count": len(entries)}


@mcp.tool()
def parse_cif_surface(
    cif_path: Optional[str] = None,
    cif_content: Optional[str] = None,
    hkl: Optional[list[int]] = None,
    n_faces: int = 5,
) -> dict:
    """
    Extract 2D surface lattice parameters from a CIF crystal structure file.

    When hkl is given (e.g. [0,0,1]), returns the lattice for that specific
    surface plane.  When hkl is omitted, uses the BFDH morphology rule to
    automatically rank the most prominent crystal faces and returns the
    top n_faces surfaces.

    Parameters
    ----------
    cif_path    : Absolute path to a CIF file on disk.
    cif_content : Raw CIF text (alternative to cif_path; useful when
                  chaining with mofstructure remove_guest output).
    hkl         : Miller indices [h, k, l] for a specific surface plane.
                  If omitted, BFDH automatic face selection is used.
    n_faces     : Number of BFDH faces to return (default 5, ignored when
                  hkl is specified).

    Returns
    -------
    dict with cell_params (3D unit cell) and faces list, each containing
    hkl, interplanar spacing d_hkl, and 2D lattice parameters (a, b, gamma).
    """
    real_path = _resolve_cif_path(cif_path, cif_content)
    is_temp = cif_content is not None
    try:
        cell = read_cell(real_path)
        cell_params = {
            k: cell[k] for k in ("a", "b", "c", "alpha", "beta", "gamma")
            if k in cell
        }
        if "name" in cell:
            cell_params["name"] = str(cell["name"])

        faces_out = []
        if hkl is not None:
            hkl_tuple = tuple(hkl)
            lat = surface_lattice(real_path, hkl=hkl_tuple)
            faces_out.append({
                "hkl": list(hkl_tuple),
                **_lattice_info(lat),
            })
        else:
            faces = best_surface_lattice(real_path, n_faces=n_faces)
            for hkl_t, d_hkl, lat in faces:
                faces_out.append({
                    "hkl": list(hkl_t),
                    "d_hkl": round(d_hkl, 3),
                    **_lattice_info(lat),
                })

        return {"cell_params": cell_params, "faces": faces_out}
    finally:
        if is_temp:
            os.unlink(real_path)


@mcp.tool()
def analyze_epitaxy(
    substrate_name: Optional[str] = None,
    substrate_a: Optional[float] = None,
    substrate_b: Optional[float] = None,
    substrate_gamma_deg: Optional[float] = None,
    cif_path: Optional[str] = None,
    cif_content: Optional[str] = None,
    hkl: Optional[list[int]] = None,
    overlayer_a: Optional[float] = None,
    overlayer_b: Optional[float] = None,
    overlayer_gamma_deg: Optional[float] = None,
    sigma: float = 0.3,
    eta_tol: float = 0.05,
    top_n: int = 10,
    output_dir: Optional[str] = None,
) -> dict:
    """
    Run the full three-level epitaxy matching pipeline and return structured
    results including algebraic feasibility, coincidence peaks, and ranked
    supercell matches with strain tensors.

    The substrate can be specified by name (e.g. 'Au_111') or by explicit
    lattice parameters (substrate_a, substrate_b, substrate_gamma_deg).

    The overlayer can be specified as:
      - A CIF file (cif_path or cif_content), optionally with hkl to select
        a specific surface plane (otherwise the best BFDH face is used).
      - Explicit lattice parameters (overlayer_a, overlayer_b,
        overlayer_gamma_deg).

    Parameters
    ----------
    substrate_name      : Built-in substrate key (e.g. 'Au_111', 'Cu_111').
    substrate_a/b/gamma : Explicit substrate lattice (Angstrom, degrees).
    cif_path            : Path to overlayer CIF file.
    cif_content         : Raw CIF text for the overlayer.
    hkl                 : Miller indices [h,k,l] of the overlayer surface.
    overlayer_a/b/gamma : Explicit overlayer lattice (Angstrom, degrees).
    sigma               : Gaussian width for coincidence function (default 0.3).
    eta_tol             : Maximum Green-Lagrange strain norm (default 0.05).
    top_n               : Max number of matches to return (default 10).
    output_dir          : If provided, generate match_card.png and report.pdf
                          in this directory.

    Returns
    -------
    dict with keys: substrate, overlayer, level0, level1, level2, and
    optionally files (paths to generated PNG/PDF).
    """
    lat_sub = _resolve_substrate(
        substrate_name, substrate_a, substrate_b, substrate_gamma_deg,
    )
    lat_mof = _resolve_overlayer(
        cif_path, cif_content, hkl,
        overlayer_a, overlayer_b, overlayer_gamma_deg,
    )

    # ── Level 0: algebraic discriminant ───────────────────────────────
    l0 = discriminant.check(lat_sub, lat_mof)
    level0 = {
        "feasible": l0.feasible,
        "delta_sub": round(l0.delta_sub, 4),
        "delta_mof": round(l0.delta_mof, 4),
        "ratio": round(l0.ratio, 4),
        "message": l0.message,
    }

    # ── Level 1: coincidence function peaks ───────────────────────────
    l1 = coincidence.compute(lat_sub, lat_mof, sigma=sigma)
    level1 = {
        "n_peaks": len(l1.theta_peaks),
        "peaks": [
            {"theta_deg": round(t, 2), "phi": round(p, 6)}
            for t, p in zip(l1.theta_peaks[:top_n], l1.phi_peaks[:top_n])
        ],
    }

    # ── Level 2: supercell matching ───────────────────────────────────
    all_matches: list[supercell.MatchResult] = []
    for th in l1.theta_peaks[:8]:
        all_matches.extend(
            supercell.find_matches(lat_sub, lat_mof, theta_deg=th, eta_tol=eta_tol)
        )
    all_matches.sort(key=lambda m: m.eta)

    # Deduplicate by (eta, area) rounded
    seen = set()
    unique: list[supercell.MatchResult] = []
    for m in all_matches:
        key = (round(m.eta, 6), round(m.area, 1))
        if key not in seen:
            seen.add(key)
            unique.append(m)
    all_matches = unique[:top_n]

    level2 = {
        "n_matches": len(all_matches),
        "matches": [
            {
                "theta_deg": round(m.theta_deg, 2),
                "M": m.M.tolist(),
                "N": m.N.tolist(),
                "eta": round(m.eta, 6),
                "eps_11": round(m.strain[0, 0], 5),
                "eps_22": round(m.strain[1, 1], 5),
                "eps_12": round(m.strain[0, 1], 5),
                "area_A2": round(m.area, 1),
            }
            for m in all_matches
        ],
    }

    # ── Optional: generate visualisation outputs ──────────────────────
    files = None
    if output_dir is not None and all_matches:
        os.makedirs(output_dir, exist_ok=True)
        best = all_matches[0]
        files = {}

        png_path = os.path.join(output_dir, "match_card.png")
        try:
            fig = plot_match_card(lat_sub, lat_mof, l1, best, save_path=png_path, dpi=150)
            plt.close(fig)
            files["png_path"] = os.path.abspath(png_path)
        except Exception:
            pass

        pdf_path = os.path.join(output_dir, "report.pdf")
        try:
            generate_pdf_report(
                lat_sub, lat_mof, l1, all_matches,
                title=f"{lat_mof.label}  on  {lat_sub.label}",
                save_path=pdf_path,
                top_n=top_n,
            )
            files["pdf_path"] = os.path.abspath(pdf_path)
        except Exception:
            pass

    result = {
        "substrate": _lattice_info(lat_sub),
        "overlayer": _lattice_info(lat_mof),
        "level0": level0,
        "level1": level1,
        "level2": level2,
    }
    if files:
        result["files"] = files
    return result


@mcp.tool()
def generate_coincidence_animation(
    substrate_name: Optional[str] = None,
    substrate_a: Optional[float] = None,
    substrate_b: Optional[float] = None,
    substrate_gamma_deg: Optional[float] = None,
    cif_path: Optional[str] = None,
    cif_content: Optional[str] = None,
    hkl: Optional[list[int]] = None,
    overlayer_a: Optional[float] = None,
    overlayer_b: Optional[float] = None,
    overlayer_gamma_deg: Optional[float] = None,
    sigma: float = 0.3,
    n_frames: int = 360,
    fps: int = 30,
    output_dir: Optional[str] = None,
) -> dict:
    """
    Generate a GIF animation of the coincidence lattice search process.

    The animation shows the MOF reciprocal lattice rotating over the
    substrate reciprocal lattice (left panel), while the coincidence
    function Phi(theta) is drawn in real time (right panel).

    Substrate and overlayer inputs follow the same convention as
    analyze_epitaxy.

    Parameters
    ----------
    substrate_name      : Built-in substrate key (e.g. 'Au_111', 'Cu_111').
    substrate_a/b/gamma : Explicit substrate lattice (Angstrom, degrees).
    cif_path            : Path to overlayer CIF file.
    cif_content         : Raw CIF text for the overlayer.
    hkl                 : Miller indices [h,k,l] of the overlayer surface.
    overlayer_a/b/gamma : Explicit overlayer lattice (Angstrom, degrees).
    sigma               : Gaussian width for coincidence score (default 0.3).
    n_frames            : Number of animation frames / angular steps (default 360).
    fps                 : Frames per second in the output GIF (default 30).
    output_dir          : Directory to save the GIF. If omitted, a temp dir is used.

    Returns
    -------
    dict with gif_path (absolute path to the saved .gif file) and metadata.
    """
    lat_sub = _resolve_substrate(
        substrate_name, substrate_a, substrate_b, substrate_gamma_deg,
    )
    lat_mof = _resolve_overlayer(
        cif_path, cif_content, hkl,
        overlayer_a, overlayer_b, overlayer_gamma_deg,
    )

    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix="miqrophi_anim_")
    else:
        os.makedirs(output_dir, exist_ok=True)

    gif_path = os.path.join(output_dir, "coincidence_search.gif")
    animate_coincidence_search(
        lat_sub, lat_mof,
        sigma=sigma,
        n_frames=n_frames,
        save_path=gif_path,
        fps=fps,
    )
    plt.close("all")

    return {
        "gif_path": os.path.abspath(gif_path),
        "substrate": _lattice_info(lat_sub),
        "overlayer": _lattice_info(lat_mof),
        "n_frames": n_frames,
        "fps": fps,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Run the MCP server over stdio (for featherflow / Claude Desktop)."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
