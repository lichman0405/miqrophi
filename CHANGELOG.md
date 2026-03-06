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

## [0.3.0] — 2026-03-06

### Added

- **Screw-axis / glide-plane extinctions** (`cif_parser.py`): `_SG_ZONE_CONDS`
  dictionary (ITA Vol A, Table 2.2.13) covers ~55 of the most common MOF space
  groups (P2₁/c SG†14, P2₁2₁2₁ SG†19, Pnma SG†62, Pbca SG†61, Fd-3m SG†227,
  Ia-3d SG†230, etc.).  Zone-specific helper functions `_sg_number`,
  `_in_zone`, `_rule_ok`, `_is_zone_allowed` integrated into `bfdh_faces()`.
- **Parallel batch processing** (`batch.py`): `BatchConfig.n_jobs` parameter;
  `1` = serial (default), `-1` = all CPUs, `N > 1` = N workers via
  `concurrent.futures.ProcessPoolExecutor`.  Single-CIF logic extracted to
  `_run_single_cif()` plus module-level picklable `_worker()` function.
- **Command-line interface** (`miqrocal/cli.py`): `miqrocal substrates`,
  `miqrocal run FILE ...`, `miqrocal batch "glob/**/*.cif"` sub-commands;
  supports `--substrates`, `--n-faces`, `--output-dir`, `--tag`, `--n-jobs`,
  `--outputs`, `--eta-tol`, `--sigma`, `--quiet`.
- **Test suite** (`tests/`): 4 test modules covering
  `Lattice2D`/`SUBSTRATE_DB`, Level-0 discriminant, CIF centering + screw/glide
  zone conditions, and `batch_run` (serial vs parallel row-count equality).
- **`__version__ = "0.3.0"`** added to `miqrocal/__init__.py`.
- `run.py` simplified to use `EpitaxyMatcher` directly with `Lattice2D` inputs
  (no CIF required for the demo).

### Changed

- `bfdh_faces()` docstring updated to reflect that screw/glide extinctions are
  now applied in addition to centering conditions.
- `pyproject.toml`: repository URL corrected to `lichman0405/miqrophi`;
  `[project.scripts]` entry point `miqrocal = "miqrocal.cli:main"` added;
  version bumped to 0.3.0.
- `README.md`: clone URL corrected; File Structure section updated to include
  `batch.py`, `cli.py`, `tests/`.

---

## [0.2.0] — 2026-03-06

### Added

- **`batch_run()` / `BatchConfig`** (`batch.py`): single-call batch pipeline
  that screens a glob of CIF files against multiple substrates and writes
  per-pair match-card PNG, PDF report, and a summary CSV with traceable paths.
- **BFDH face selection with centering extinctions** (`cif_parser.py`):
  `bfdh_faces()` and `best_surface_lattice()` now apply F/I/A/B/C/R
  lattice-centering systematic-absence rules when ranking faces.
- **`SUBSTRATE_DB` expanded** to 11 entries: Au(111), Au(100), Ag(111),
  Cu(111), Pt(111), Graphene, HOPG(0001), ITO(111), MgO(001),
  SrTiO₃(001), Si(001).
- **PDF report** (`visualize.py`): `generate_pdf_report()` writes a two-page
  PDF (ranked match table + match-card figure) with no external dependencies
  beyond matplotlib.
- `README.md` extended with Quick Start sections B/C, API tables, worked
  examples, and visualisation function reference.

---

## [Unreleased]

### Planned

- Julia rewrite of Level 1 and Level 2 hot paths for high-throughput
  screening (no GIL, native multi-threading)
- Surface reconstruction database (Au(111) 22×√3, Au(100) 5×20, etc.)
- DFT surface-energy weighting to replace the BFDH approximation
- GUI / Jupyter widget for interactive face selection and result inspection
