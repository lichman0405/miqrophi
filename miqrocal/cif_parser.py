"""
CIF parser — extract a 2D surface unit cell from a 3D crystal CIF file.

Given a CIF file and a set of Miller indices (h, k, l), this module:

  1. Reads the six unit-cell parameters (a, b, c, alpha, beta, gamma) from
     the CIF ``_cell_length_*`` and ``_cell_angle_*`` keys.
  2. Builds the 3x3 Cartesian lattice matrix A (rows = real-space basis
     vectors a1, a2, a3).
  3. Enumerates all short integer-combination lattice vectors that lie in
     the (hkl) surface plane (condition: h·m1 + k·m2 + l·m3 = 0).
  4. Returns the two shortest linearly independent in-plane vectors as a
     ``Lattice2D`` object.

The parser is intentionally lightweight — it handles only the cell-parameter
block and requires no third-party crystallography library (pymatgen, ase, etc.).

Reference
---------
docs/lattice_matching_theory.md, Section 2 (lattice geometry).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np

from .lattice import Lattice2D


# ---------------------------------------------------------------------------
# CIF key → value extraction
# ---------------------------------------------------------------------------

_CIF_KEYS = {
    "a":     re.compile(r"_cell_length_a\s+([\d.]+)", re.IGNORECASE),
    "b":     re.compile(r"_cell_length_b\s+([\d.]+)", re.IGNORECASE),
    "c":     re.compile(r"_cell_length_c\s+([\d.]+)", re.IGNORECASE),
    "alpha": re.compile(r"_cell_angle_alpha\s+([\d.]+)", re.IGNORECASE),
    "beta":  re.compile(r"_cell_angle_beta\s+([\d.]+)", re.IGNORECASE),
    "gamma": re.compile(r"_cell_angle_gamma\s+([\d.]+)", re.IGNORECASE),
    "name":  re.compile(
        r"(?:_chemical_name_common|_chemical_formula_sum)\s+['\"]?([^\n'\"]+)['\"]?",
        re.IGNORECASE,
    ),
}

# Space-group parsing — used by BFDH systematic-absence filter
_SG_SYMBOL_PAT = re.compile(
    r"(?:_symmetry_space_group_name_H-M|_space_group_name_H-M_alt)"
    r"\s+['\"]?\s*([A-Za-z])",
    re.IGNORECASE,
)
_SG_NUMBER_PAT = re.compile(
    r"(?:_symmetry_Int_Tables_number|_space_group_IT_number)\s+(\d+)",
    re.IGNORECASE,
)

# Space-group number → lattice-centering type.  Only non-P entries listed;
# everything else defaults to 'P'.
_SG_LATTICE: dict[int, str] = {
    # C-centred monoclinic
    **{n: "C" for n in [5, 8, 9, 12, 13, 14, 15]},
    # C/A orthorhombic
    **{n: "C" for n in [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]},
    **{n: "A" for n in [38, 39, 40, 41]},
    # I orthorhombic
    **{n: "I" for n in [23, 24, 44, 45, 46, 71, 72, 73, 74]},
    # F orthorhombic
    **{n: "F" for n in [22, 42, 43, 69, 70]},
    # I tetragonal
    **{n: "I" for n in [79, 80, 82, 87, 88, 97, 98,
                         107, 108, 109, 110, 119, 120, 121, 122,
                         139, 140, 141, 142]},
    # R trigonal (hexagonal setting)
    **{n: "R" for n in [146, 148, 155, 160, 161, 166, 167]},
    # I cubic
    **{n: "I" for n in [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]},
    # F cubic
    **{n: "F" for n in [196, 202, 203, 209, 210, 216, 225, 226, 227, 228]},
}


def _lattice_type(cif_text: str) -> str:
    """
    Extract the lattice-centering type (P, I, F, A, B, C, R) from CIF text.

    First tries the Hermann–Mauguin symbol (first character); falls back to
    the space-group number; returns ``'P'`` if neither is present.
    """
    m = _SG_SYMBOL_PAT.search(cif_text)
    if m:
        lt = m.group(1).upper()
        if lt in {"P", "I", "F", "A", "B", "C", "R", "H"}:
            return lt
    m = _SG_NUMBER_PAT.search(cif_text)
    if m:
        return _SG_LATTICE.get(int(m.group(1)), "P")
    return "P"


def _is_allowed(h: int, k: int, l: int, lattice_type: str) -> bool:
    """
    Return True if the reflection (h, k, l) is systematically allowed for
    the given lattice-centering type.

    Rules applied (general conditions from ITA):
    * P — all allowed
    * I — h+k+l = 2n
    * F — h, k, l all odd or all even (unmixed parities)
    * A — k+l = 2n
    * B — h+l = 2n
    * C — h+k = 2n
    * R — −h+k+l = 3n  (hexagonal setting, obverse)
    """
    lt = lattice_type.upper()
    if lt == "P":
        return True
    if lt == "I":
        return (h + k + l) % 2 == 0
    if lt == "F":
        parities = {h % 2, k % 2, l % 2}
        return len(parities) == 1          # all same parity
    if lt == "A":
        return (k + l) % 2 == 0
    if lt == "B":
        return (h + l) % 2 == 0
    if lt == "C":
        return (h + k) % 2 == 0
    if lt == "R":
        return (-h + k + l) % 3 == 0
    return True                            # H or unknown: allow all


def _parse_cell(cif_text: str) -> dict[str, float | str]:
    """
    Extract cell parameters and optional compound name from CIF text.

    Returns a dict with float keys ``a, b, c, alpha, beta, gamma``
    and optional string key ``name``.

    Raises
    ------
    ValueError if any required cell parameter is missing.
    """
    result: dict[str, float | str] = {}
    for key, pat in _CIF_KEYS.items():
        m = pat.search(cif_text)
        if m:
            result[key] = float(m.group(1)) if key != "name" else m.group(1).strip()

    required = {"a", "b", "c", "alpha", "beta", "gamma"}
    missing  = required - result.keys()
    if missing:
        raise ValueError(f"CIF is missing required cell parameter(s): {missing}")

    return result


# ---------------------------------------------------------------------------
# 3D lattice matrix
# ---------------------------------------------------------------------------

def _cartesian_matrix(
    a: float, b: float, c: float,
    alpha_deg: float, beta_deg: float, gamma_deg: float,
) -> np.ndarray:
    """
    Build the 3x3 Cartesian lattice matrix.  Rows are real-space basis
    vectors (a1, a2, a3) using the standard crystallographic convention:

        a1 = a * [1, 0, 0]
        a2 = b * [cos(gamma), sin(gamma), 0]
        a3 = c * [cos(beta), (cos(alpha) - cos(beta)*cos(gamma))/sin(gamma), v33]
    """
    al = np.radians(alpha_deg)
    be = np.radians(beta_deg)
    ga = np.radians(gamma_deg)

    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([b * np.cos(ga), b * np.sin(ga), 0.0])

    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz_sq = c ** 2 - cx ** 2 - cy ** 2
    cz = float(np.sqrt(max(cz_sq, 0.0)))
    a3 = np.array([cx, cy, cz])

    return np.array([a1, a2, a3], dtype=float)


# ---------------------------------------------------------------------------
# 2D surface lattice from Miller indices
# ---------------------------------------------------------------------------

def _in_plane_vectors(
    A3: np.ndarray,
    hkl: tuple[int, int, int],
    n_max: int = 6,
) -> list[np.ndarray]:
    """
    Return all Cartesian lattice vectors v = m1·a1 + m2·a2 + m3·a3 that lie
    in the (hkl) plane, i.e. satisfy  h·m1 + k·m2 + l·m3 = 0, for integer
    combinations with |mi| ≤ n_max (excluding the zero vector).

    The list is sorted by vector length (ascending).
    """
    h, k, l = hkl
    vecs: list[tuple[float, np.ndarray]] = []

    for m1 in range(-n_max, n_max + 1):
        for m2 in range(-n_max, n_max + 1):
            for m3 in range(-n_max, n_max + 1):
                if m1 == 0 and m2 == 0 and m3 == 0:
                    continue
                if h * m1 + k * m2 + l * m3 != 0:
                    continue
                v = m1 * A3[0] + m2 * A3[1] + m3 * A3[2]
                vecs.append((float(np.linalg.norm(v)), v))

    vecs.sort(key=lambda x: x[0])
    return [v for _, v in vecs]


def _pick_two_shortest(vecs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """
    From a list of in-plane vectors (sorted by length), return the two
    shortest linearly independent ones (non-parallel, det != 0).
    """
    v1 = vecs[0]
    for v2 in vecs[1:]:
        # 2D det via cross product magnitude / (|v1||v2|) — reject near-parallel
        cross = abs(float(v1[0] * v2[1] - v1[1] * v2[0])) if len(v1) >= 2 else 0.0
        if len(v1) == 3:
            cross = float(np.linalg.norm(np.cross(v1, v2)))
        if cross > 1e-4 * np.linalg.norm(v1) * np.linalg.norm(v2):
            return v1, v2
    raise ValueError(
        "Could not find two linearly independent in-plane vectors. "
        "Try increasing n_max in surface_lattice()."
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def read_cell(cif_path: str | Path) -> dict[str, float | str]:
    """
    Parse a CIF file and return its unit-cell parameters as a plain dict.

    Returns
    -------
    dict with float keys ``a, b, c, alpha, beta, gamma`` (Å and degrees)
    and optional string key ``name``.

    Example
    -------
    >>> cell = read_cell("examples/HKUST-1.cif")
    >>> cell
    {'a': 26.343, 'b': 26.343, 'c': 26.343,
     'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0,
     'name': 'HKUST-1'}
    """
    text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    return _parse_cell(text)


def surface_lattice(
    cif_path: str | Path,
    hkl: tuple[int, int, int] = (0, 0, 1),
    *,
    label: Optional[str] = None,
    n_max: int = 6,
) -> Lattice2D:
    """
    Extract the 2D surface unit cell for the (hkl) plane from a CIF file.

    The two shortest linearly independent lattice vectors that lie in the
    (hkl) cutting plane are found by integer enumeration.  Their lengths
    and included angle define the returned ``Lattice2D``.

    Parameters
    ----------
    cif_path : path to a CIF file (only cell-parameter block is required)
    hkl      : Miller indices of the surface plane (default: (0,0,1))
    label    : human-readable identifier; auto-generated from the CIF name
               and Miller indices if *None*
    n_max    : maximum absolute value of integer coefficients to search
               (increase for very high-index surfaces)

    Returns
    -------
    ``Lattice2D`` with the 2D surface unit cell parameters

    Raises
    ------
    ValueError  if required cell parameters are missing, or if two
                linearly independent in-plane vectors cannot be found.

    Example
    -------
    >>> from miqrocal.cif_parser import surface_lattice
    >>> lat = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
    >>> lat
    Lattice2D(a=26.343, b=26.343, gamma_deg=90.0, label='HKUST-1 (001)')

    >>> lat_111 = surface_lattice("examples/HKUST-1.cif", hkl=(1, 1, 1))
    >>> round(lat_111.a, 3)
    18.628
    """
    cell = read_cell(cif_path)

    A3 = _cartesian_matrix(
        float(cell["a"]), float(cell["b"]), float(cell["c"]),
        float(cell["alpha"]), float(cell["beta"]), float(cell["gamma"]),
    )

    in_plane = _in_plane_vectors(A3, hkl, n_max=n_max)
    if len(in_plane) < 2:
        raise ValueError(
            f"Fewer than 2 in-plane vectors found for hkl={hkl}. "
            "Try increasing n_max."
        )

    v1, v2 = _pick_two_shortest(in_plane)

    a_2d     = float(np.linalg.norm(v1))
    b_2d     = float(np.linalg.norm(v2))
    cos_g    = float(np.dot(v1, v2) / (a_2d * b_2d))
    # numerical guard
    cos_g    = max(-1.0, min(1.0, cos_g))
    gamma_2d = float(np.degrees(np.arccos(cos_g)))

    if label is None:
        name = str(cell.get("name", Path(cif_path).stem))
        h, k, l = hkl
        label = f"{name} ({h}{k}{l})"

    return Lattice2D(a=a_2d, b=b_2d, gamma_deg=gamma_2d, label=label)


# ---------------------------------------------------------------------------
# B2 — BFDH automatic face selection
# ---------------------------------------------------------------------------

def _recip_3d(A3: np.ndarray) -> np.ndarray:
    """3x3 reciprocal lattice matrix.  Rows are b1, b2, b3  (A · B^T = 2π I)."""
    return 2.0 * np.pi * np.linalg.inv(A3).T


def _canonical_hkl(h: int, k: int, l: int) -> bool:
    """
    Return True if (h, k, l) is the canonical representative of its Friedel pair.

    The canonical form has its first non-zero index positive, ensuring that
    (hkl) and (-h,-k,-l) are not both included in the face list.
    """
    for x in (h, k, l):
        if x > 0:
            return True
        if x < 0:
            return False
    return False   # (0, 0, 0) — excluded separately


def bfdh_faces(
    cif_path: str | Path,
    *,
    n_faces: int = 5,
    hkl_max: int = 3,
) -> list[tuple[tuple[int, int, int], float]]:
    """
    Rank crystal faces by interplanar spacing using the BFDH morphology rule.

    The **Bravais–Friedel–Donnay–Harker (BFDH)** rule states that the faces
    with the largest interplanar spacings d_hkl grow slowest and are therefore
    the most prominent on the crystal morphology.  These are also the most
    likely surfaces to be exposed in a thin-film epitaxy experiment.

    The interplanar spacing is computed from the reciprocal lattice vector:

        d_hkl = 2π / |G_hkl|    where G_hkl = h·b1 + k·b2 + l·b3

    **Note on systematic absences**: lattice-centering extinction rules are
    applied automatically.  The space group is read from the CIF
    (``_symmetry_space_group_name_H-M`` or ``_symmetry_Int_Tables_number``);
    if absent the structure is assumed primitive (P).  Screw-axis and
    glide-plane extinctions (which depend on zone-specific conditions) are
    *not* applied — for BFDH morphology prediction the centering condition
    alone provides a large improvement over the uncorrected list.

    Parameters
    ----------
    cif_path : path to a CIF file
    n_faces  : number of top faces to return (default 5)
    hkl_max  : maximum absolute value of h, k, l to enumerate (default 3)

    Returns
    -------
    List of ``((h, k, l), d_hkl)`` tuples sorted by d_hkl descending.

    Example
    -------
    >>> # HKUST-1 (Fm-3m, F-centering): {100} family absent → {111} family first
    >>> faces = bfdh_faces("examples/HKUST-1.cif", n_faces=3)
    >>> [(hkl, round(d, 2)) for hkl, d in faces]
    [((1, 1, 1), 15.21), ((1, 1, -1), 15.21), ((1, -1, 1), 15.21)]
    """
    cell = read_cell(cif_path)
    cif_text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    lat_type = _lattice_type(cif_text)

    A3   = _cartesian_matrix(
        float(cell["a"]), float(cell["b"]), float(cell["c"]),
        float(cell["alpha"]), float(cell["beta"]), float(cell["gamma"]),
    )
    B3 = _recip_3d(A3)

    faces: list[tuple[float, tuple[int, int, int]]] = []
    for h in range(-hkl_max, hkl_max + 1):
        for k in range(-hkl_max, hkl_max + 1):
            for l in range(-hkl_max, hkl_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                if not _canonical_hkl(h, k, l):
                    continue
                if not _is_allowed(h, k, l, lat_type):
                    continue
                G   = h * B3[0] + k * B3[1] + l * B3[2]
                Gn  = float(np.linalg.norm(G))
                d   = 2.0 * np.pi / Gn
                faces.append((d, (h, k, l)))

    faces.sort(reverse=True)
    return [(hkl, d) for d, hkl in faces[:n_faces]]


def best_surface_lattice(
    cif_path: str | Path,
    *,
    n_faces: int = 5,
    hkl_max: int = 3,
    n_max: int = 6,
) -> list[tuple[tuple[int, int, int], float, Lattice2D]]:
    """
    Automatically select the most morphologically important surface planes
    from a CIF file and return their 2D lattice objects.

    Combines :func:`bfdh_faces` (BFDH ranking) with
    :func:`surface_lattice` (2D cell extraction) in one call.

    Parameters
    ----------
    cif_path : path to a CIF file
    n_faces  : number of top BFDH faces to evaluate (default 5)
    hkl_max  : maximum |hkl| index for BFDH enumeration (default 3)
    n_max    : integer-vector search range for in-plane vector enumeration

    Returns
    -------
    List of ``((h, k, l), d_hkl, Lattice2D)`` tuples, sorted by d_hkl
    descending (most prominent face first).  Faces for which 2D lattice
    extraction fails are silently skipped.

    Example
    -------
    >>> # HKUST-1 (Fm-3m, F-centering): {100}/{010}/{001} absent → {111} family leads
    >>> results = best_surface_lattice("examples/HKUST-1.cif", n_faces=3)
    >>> for hkl, d, lat in results:
    ...     print(hkl, f"d={d:.2f}A", f"a={lat.a:.2f} b={lat.b:.2f} g={lat.gamma_deg:.1f}")
    (1, 1, 1)   d=15.21A  a=37.26  b=37.26  g=60.0
    (1, 1, -1)  d=15.21A  a=37.26  b=37.26  g=60.0
    (1, -1, 1)  d=15.21A  a=37.26  b=37.26  g=60.0
    """
    ranked = bfdh_faces(cif_path, n_faces=n_faces, hkl_max=hkl_max)
    cell   = read_cell(cif_path)
    name   = str(cell.get("name", Path(cif_path).stem))

    results: list[tuple[tuple[int, int, int], float, Lattice2D]] = []
    for hkl, d in ranked:
        h, k, l = hkl
        try:
            lat = surface_lattice(
                cif_path,
                hkl=hkl,
                label=f"{name} ({h}{k}{l})",
                n_max=n_max,
            )
            results.append((hkl, d, lat))
        except ValueError:
            pass   # skip degenerate faces

    return results
