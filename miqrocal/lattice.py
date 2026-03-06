"""2D lattice representation and predefined substrate database."""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass
class Lattice2D:
    """
    A 2D periodic lattice defined by two basis vector lengths and their
    included angle.  The first basis vector is aligned with the x-axis.

    The 2x2 lattice matrix A stores basis vectors as rows:
        A[0] = a1  (along x)
        A[1] = a2  (at gamma_deg from x)

    Attributes
    ----------
    a, b       : basis vector lengths (Angstrom)
    gamma_deg  : angle between basis vectors (degrees)
    label      : human-readable identifier
    A          : (2, 2) Cartesian lattice matrix (set by __post_init__)
    omega      : unit cell area in Angstrom^2  (set by __post_init__)
    """

    a:         float
    b:         float
    gamma_deg: float
    label:     str = ""

    def __post_init__(self) -> None:
        g = np.radians(self.gamma_deg)
        self.A = np.array(
            [[self.a,             0.0          ],
             [self.b * np.cos(g), self.b * np.sin(g)]],
            dtype=float,
        )
        self.omega: float = abs(float(np.linalg.det(self.A)))

    def recip_matrix(self) -> np.ndarray:
        """
        Reciprocal lattice matrix B (rows are reciprocal vectors), defined by
            A . B^T = 2*pi * I
        """
        return 2.0 * np.pi * np.linalg.inv(self.A).T


# ---------------------------------------------------------------------------
# Predefined substrate database
# ---------------------------------------------------------------------------

SUBSTRATE_DB: dict[str, Lattice2D] = {
    # ── Noble-metal single crystals ────────────────────────────────────────
    "Au_111":      Lattice2D(2.88,  2.88,  120.0, "Au(111)"),       # FCC, hex
    "Au_100":      Lattice2D(4.08,  4.08,   90.0, "Au(100)"),       # FCC, square
    "Ag_111":      Lattice2D(2.89,  2.89,  120.0, "Ag(111)"),       # FCC, hex
    "Cu_111":      Lattice2D(2.56,  2.56,  120.0, "Cu(111)"),       # FCC, hex
    "Pt_111":      Lattice2D(2.775, 2.775, 120.0, "Pt(111)"),       # FCC, hex
    # ── Carbon substrates ───────────────────────────────────────────────────
    "Graphene":    Lattice2D(2.46,  2.46,  120.0, "Graphene"),      # monolayer
    "HOPG":        Lattice2D(2.46,  2.46,  120.0, "HOPG(0001)"),    # graphite basal plane
    # ── Oxide substrates ────────────────────────────────────────────────────
    "ITO_111":     Lattice2D(4.13,  4.13,  120.0, "ITO(111)"),
    "MgO_001":     Lattice2D(4.211, 4.211,  90.0, "MgO(001)"),      # rock salt, square
    "SrTiO3_001":  Lattice2D(3.905, 3.905,  90.0, "SrTiO3(001)"),   # perovskite, square
    # ── Semiconductor substrates ────────────────────────────────────────────
    "Si_001":      Lattice2D(5.431, 5.431,  90.0, "Si(001)"),       # diamond cubic
}
