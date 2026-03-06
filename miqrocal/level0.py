"""
Level 0 — Quadratic form discriminant check.

Determines whether two 2D lattices can be exactly commensurate based on
the algebraic structure of their associated binary quadratic forms.
Exact commensurability requires both lattices to belong to the same
imaginary quadratic field, which is tested by checking whether the ratio
of their normalised discriminants is a rational perfect square.

Reference: docs/lattice_matching_theory.md, Section 3
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from .lattice import Lattice2D


@dataclass(frozen=True)
class Level0Result:
    """Result of the Level-0 symmetry pre-check."""

    feasible:  bool   # True if exact commensurability is algebraically allowed
    delta_sub: float  # normalised discriminant of the substrate lattice
    delta_mof: float  # normalised discriminant of the overlayer lattice
    ratio:     float  # delta_mof / delta_sub
    message:   str


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _normalised_discriminant(gamma_deg: float) -> float:
    """
    Return the normalised discriminant  delta = -4 * sin^2(gamma)  of the
    canonical binary quadratic form for a 2D lattice with angle gamma.

    The a^2 * b^2 dimensional factor is removed so that delta depends only
    on lattice *type*, not on scale:
        square     (gamma = 90 deg)  -> delta = -4
        hexagonal  (gamma = 120 deg) -> delta = -3
    """
    return -4.0 * np.sin(np.radians(gamma_deg)) ** 2


def _is_rational_perfect_square(x: float, max_denom: int = 10_000) -> bool:
    """
    Return True if x is a rational perfect square, i.e. x = (p/q)^2 for
    integers p, q.  Uses exact integer arithmetic to avoid floating-point
    false positives.
    """
    if x < 0.0:
        return False
    frac = Fraction(x).limit_denominator(max_denom)
    p, q = frac.numerator, frac.denominator

    def _is_square(n: int) -> bool:
        s = int(round(n ** 0.5))
        return s * s == n

    return _is_square(p) and _is_square(q)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def check(lat_sub: Lattice2D, lat_mof: Lattice2D) -> Level0Result:
    """
    Run the Level-0 symmetry pre-check.

    The two lattices are algebraically compatible for exact commensurability
    when delta_mof / delta_sub is a rational perfect square (i.e. both
    lattices belong to the same imaginary quadratic field Q(sqrt(delta))).

    Examples
    --------
    Square on square   -> ratio = (-4)/(-4) = 1   -> feasible
    Hex    on hex      -> ratio = (-3)/(-3) = 1   -> feasible
    Square on hex      -> ratio = (-4)/(-3) = 4/3 -> forbidden
    """
    d_s = _normalised_discriminant(lat_sub.gamma_deg)
    d_o = _normalised_discriminant(lat_mof.gamma_deg)
    ratio = d_o / d_s
    ok = _is_rational_perfect_square(ratio)
    msg = (
        "exact commensurability algebraically feasible"
        if ok else
        "symmetry-forbidden: discriminant ratio is not a rational perfect square"
    )
    return Level0Result(
        feasible=ok, delta_sub=d_s, delta_mof=d_o, ratio=ratio, message=msg
    )
