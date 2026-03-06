"""
miqrocal — 2D epitaxial lattice matching via three-level hierarchical screening.

Public API
----------
EpitaxyMatcher   : main pipeline class
MatcherConfig    : configuration dataclass
Lattice2D        : 2D lattice representation
SUBSTRATE_DB     : predefined substrate lattice database

CIF input (no third-party crystallography library required)
------------------------------------------------------------
read_cell            : parse unit-cell parameters from a CIF file
surface_lattice      : extract Lattice2D for a given (hkl) surface from a CIF
bfdh_faces           : rank faces by interplanar spacing (BFDH morphology rule)
best_surface_lattice : auto-select top BFDH faces and return their Lattice2D objects

Visualisation (requires matplotlib)
------------------------------------
plot_phi_curve       : Level-1 coincidence function Phi(theta)
plot_lattice_overlay : real-space supercell overlay
plot_leed_pattern    : simulated LEED diffraction pattern
plot_strain_ellipse  : principal-strain ellipse and tensor summary
plot_match_card      : 2x2 combined match-card figure
generate_pdf_report  : two-page PDF report (text summary + match-card figure)
"""

from .lattice     import Lattice2D, SUBSTRATE_DB
from .matcher     import EpitaxyMatcher, MatcherConfig
from .cif_parser  import read_cell, surface_lattice, bfdh_faces, best_surface_lattice
from .visualize   import (
    plot_phi_curve,
    plot_lattice_overlay,
    plot_leed_pattern,
    plot_strain_ellipse,
    plot_match_card,
    generate_pdf_report,
)

__all__ = [
    "Lattice2D", "SUBSTRATE_DB",
    "EpitaxyMatcher", "MatcherConfig",
    "read_cell", "surface_lattice", "bfdh_faces", "best_surface_lattice",
    "plot_phi_curve", "plot_lattice_overlay", "plot_leed_pattern",
    "plot_strain_ellipse", "plot_match_card", "generate_pdf_report",
]
