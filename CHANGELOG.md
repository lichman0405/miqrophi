# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version numbers follow [Semantic Versioning](https://semver.org/).

---

## [0.1.0] — 2026-03-06

### Added

- `miqrocal` package with modular three-level architecture:
  - `level0.py` — algebraic symmetry pre-check via quadratic-form
    discriminant ratio; rejects square/hexagonal mismatches in O(1) without
    any numerical tolerance
  - `level1.py` — reciprocal-space coincidence function Φ(θ) expanded as a
    Fourier series via the Jacobi–Anger identity; FFT-based peak finding
    identifies optimal rotation angles in O(N_G² + K log K)
  - `level2.py` — 4-D floating-point LLL lattice reduction finds optimal
    integer supercell matrices (M, N) for each candidate angle; mismatch
    quantified by the rotation-invariant Green–Lagrange strain tensor ε
  - `matcher.py` — `EpitaxyMatcher` pipeline class with `MatcherConfig`
    dataclass for centralised parameter management
  - `lattice.py` — `Lattice2D` dataclass and `SUBSTRATE_DB` (Graphene,
    Au(111), Au(100), ITO(111), Cu(111), Ag(111))
- `run.py` — demonstration entry point with three validation test cases
- `docs/lattice_matching_theory.md` — full mathematical derivation covering
  all three levels, comparison table, and references
- `README.md` — installation, quick start, API reference, worked examples,
  and known limitations
- `pyproject.toml` — PEP 517/518 package metadata

### Design decisions

- Green–Lagrange tensor ε replaces ad-hoc length/angle percentage checks;
  it is rotation-invariant and each component has explicit physical meaning
- LLL `lambda` parameter scan exposes the strain-vs-supercell-size Pareto
  front; no artificial `max_index` cutoff required
- Level 0 uses exact rational arithmetic (`fractions.Fraction`) to avoid
  floating-point false positives in the perfect-square test
- Pure Python with numpy/scipy/pandas; no compiled extensions required

### Validated test cases

| System | Expected | Result |
|--------|----------|--------|
| HKUST-1(010) / Au(111) | L0 forbidden (square vs hexagonal) | ✓ `L0_feasible=False` |
| HKUST-1(010) / Au(100) | Commensurate (square vs square) | ✓ η ≈ 0.0021, θ ≈ 51° |
| Hex-MOF / Au(111) | Multiple equivalent rotation domains | ✓ C6-symmetric peak pattern |

---

## [Unreleased]

### Planned

- CIF file parser to extract 2D surface lattice parameters from 3D crystal
  structures given a Miller index (hkl)
- Julia rewrite of Level 1 and Level 2 hot paths for high-throughput
  screening (no GIL, native multi-threading)
- `batch_run()` method in `EpitaxyMatcher` for screening many MOF/substrate
  pairs with optional multiprocessing
- Surface reconstruction database (Au(111) 22×√3, Au(100) 5×20, etc.)
- Visualisation: Φ(θ) curve plot and coincidence lattice overlay diagram
