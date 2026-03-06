# miqrocal

**miqrocal** is a Python library for predicting the epitaxial compatibility
between MOF (metal-organic framework) thin films and crystalline substrates
using a three-level hierarchical lattice-matching algorithm.

MOF thin-film epitaxy has emerged as a route to integrate porous, functional
materials with device-grade substrates \[1, 2\].  Identifying epitaxia-compatible
substrate–overlayer pairs is the first bottleneck in any growth campaign.

The classical approach to this problem is the brute-force supercell enumeration
introduced by **Zur & McGill (1984)** \[3\], which iterates over all integer
matrix pairs (M, N) and scales as O(N_m² · N_s²).  miqrocal replaces that
exhaustive search with a three-level hierarchy:

| Level | Method | Key reference | Cost |
|-------|--------|---------------|------|
| 0 | Quadratic-form discriminant check | Number theory (classical) | O(1) |
| 1 | Reciprocal-space coincidence function Φ(θ) via Jacobi–Anger expansion \[4\] | Watson 1944 \[5\] | O(N_G² + K log K) |
| 2 | Floating-point LLL lattice reduction \[6\] → Green–Lagrange strain tensor | Lenstra, Lenstra & Lovász 1982 | O(d⁴ log B) |

Surface-exposure probabilities are estimated using the **BFDH morphological
rule** \[7\] (Donnay & Harker 1937), which ranks faces by interplanar spacing
$d_{hkl}$ without requiring a force-field calculation.

The benchmark example CIF for HKUST-1 is based on the crystal structure first
reported by **Chui et al. (1999)** \[8\].

For a complete mathematical and physicochemical derivation see
[`docs/lattice_matching_theory.md`](docs/lattice_matching_theory.md).

---

## Requirements

- Python ≥ 3.10
- numpy
- scipy
- pandas

No compiled extensions required.

---

## Installation

```bash
git clone https://github.com/lichman0405/miqrophi.git
cd miqrophi
pip install -e .
```

---

## Quick Start

### Option A — automatic face selection from CIF (recommended)

```python
from miqrocal import best_surface_lattice, EpitaxyMatcher, SUBSTRATE_DB, generate_pdf_report
from miqrocal import level1, level2

# Let BFDH rule pick the most prominent face automatically
results = best_surface_lattice("examples/0000[Co][sql]2[ASR]2.cif", n_faces=3)
hkl, d_hkl, lat_mof = results[0]   # largest d_hkl → most likely exposed face
print(hkl, f"d={d_hkl:.2f} Å")    # -> (0,0,1)  d=10.62 Å

# Run matching against all substrates and find the best one
from miqrocal import EpitaxyMatcher, MatcherConfig, SUBSTRATE_DB
cfg     = MatcherConfig(sigma=0.4, eta_tol=0.05)
matcher = EpitaxyMatcher(cfg)

best_eta, best_key, best_df = 1.0, None, None
for key, lat_sub in SUBSTRATE_DB.items():
    df = matcher.run(lat_sub, lat_mof, verbose=False)
    if df is not None and df["eta"].iloc[0] < best_eta:
        best_eta, best_key, best_df = df["eta"].iloc[0], key, df
print(f"Best substrate: {best_key}  η={best_eta:.4f}")
```

### Option B — from a specific CIF face

```python
from miqrocal import EpitaxyMatcher, MatcherConfig, SUBSTRATE_DB
from miqrocal.cif_parser import surface_lattice

lat_mof = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
# -> Lattice2D(a=26.343, b=26.343, gamma_deg=90.0)

cfg     = MatcherConfig(sigma=0.4, eta_tol=0.04)
matcher = EpitaxyMatcher(cfg)
df      = matcher.run(lat_sub=SUBSTRATE_DB["Au_100"], lat_mof=lat_mof)
print(df)
```

### Option C — from known lattice parameters

```python
from miqrocal import EpitaxyMatcher, MatcherConfig, Lattice2D, SUBSTRATE_DB

lat_mof = Lattice2D(a=26.343, b=26.343, gamma_deg=90.0, label="HKUST-1(001)")
cfg     = MatcherConfig(sigma=0.4, eta_tol=0.04)
matcher = EpitaxyMatcher(cfg)
df      = matcher.run(lat_sub=SUBSTRATE_DB["Au_100"], lat_mof=lat_mof)
print(df)
```

Sample output:

```
 theta (deg)      eta  eps_11  eps_22  eps_12  area (A2)  L0_feasible
       50.76 0.002141 0.00151 0.00151    -0.0     5876.2         True
       90.00 0.020169 0.01426 0.01426     0.0     1348.4         True
```

### Generating a PDF report

```python
from miqrocal import generate_pdf_report
# (obtain l1 and sorted matches via level1.compute / level2.find_matches)
generate_pdf_report(lat_sub, lat_mof, l1, matches,
                    title="HKUST-1 / Au(100)",
                    save_path="output/report.pdf")
```

The PDF contains two pages: a text summary table (lattice params + ranked
matches) and the 2×2 match-card figure.  No external dependencies beyond
matplotlib are required.

Run the bundled demonstration (saves PNG + PDF to `output/`):

```bash
python run.py
```

For a full explanation of all inputs, outputs, and how to interpret the figures,
see [docs/usage_guide.md](docs/usage_guide.md).

---

## API Reference

### `Lattice2D`

```python
Lattice2D(a: float, b: float, gamma_deg: float, label: str = "")
```

Represents a 2D periodic lattice.  The first basis vector lies along the
x-axis.  After construction:

| Attribute | Type | Description |
|-----------|------|-------------|
| `A` | `(2, 2) ndarray` | Cartesian lattice matrix (rows = basis vectors) |
| `omega` | `float` | Unit-cell area (Å²) |

Convenience method `recip_matrix()` returns the 2×2 reciprocal lattice
matrix B satisfying A · Bᵀ = 2π I.

### `SUBSTRATE_DB`

Pre-defined substrate lattices (bulk-truncated, ideal surfaces):

| Key | Description | a (Å) | γ (°) |
|-----|-------------|--------|--------|
| `"Au_111"` | Au(111) | 2.880 | 120 |
| `"Au_100"` | Au(100) | 4.080 | 90 |
| `"Ag_111"` | Ag(111) | 2.890 | 120 |
| `"Cu_111"` | Cu(111) | 2.560 | 120 |
| `"Pt_111"` | Pt(111) | 2.775 | 120 |
| `"Graphene"` | Graphene (1×1) | 2.460 | 120 |
| `"HOPG"` | HOPG(0001) | 2.460 | 120 |
| `"ITO_111"` | ITO(111) | 4.130 | 120 |
| `"MgO_001"` | MgO(001) | 4.211 | 90 |
| `"SrTiO3_001"` | SrTiO₃(001) | 3.905 | 90 |
| `"Si_001"` | Si(001) | 5.431 | 90 |

### `MatcherConfig`

```python
MatcherConfig(
    G_cutoff      = 8.0,   # reciprocal-space truncation radius (Å⁻¹)
    sigma         = 0.3,   # Φ(θ) Gaussian peak width (Å⁻¹); controls angular tolerance
    eta_tol       = 0.05,  # Green–Lagrange strain Frobenius norm tolerance
    top_theta     = 8,     # number of Level-1 peaks forwarded to Level 2
    lambda_values = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0],  # LLL embedding weights
)
```

**Parameter guide:**

`sigma` — Larger values make Φ(θ) peaks broader and more peaks are found.  A
good starting value is ~0.1–0.5 × (2π / a_sub).

`eta_tol` — The Frobenius norm of the Green–Lagrange strain tensor.  Matches
with η < 0.02 are generally considered physically feasible.  Increasing this
threshold to 0.05–0.10 reveals approximate (strained) matches.

`lambda_values` — Each value in this list produces one LLL solution.  Small λ
(0.02) favours low-strain large supercells; large λ (1.0) favours small
supercells with potentially higher strain.  The union of results covers the
Pareto front.

### `EpitaxyMatcher`

```python
matcher = EpitaxyMatcher(config: MatcherConfig = None)
df      = matcher.run(lat_sub, lat_mof, verbose=True) -> pd.DataFrame | None
```

Returns a `DataFrame` sorted by `eta` or `None` when no match satisfies
`eta_tol`.  The `L0_feasible` column distinguishes true commensurate matches
from approximate (symmetry-forbidden) ones.

---

## Output Columns

| Column | Unit | Description |
|--------|------|-------------|
| `theta (deg)` | ° | Rotation angle of MOF supercell relative to substrate |
| `M` | — | 2×2 integer MOF supercell matrix (list of lists) |
| `N` | — | 2×2 integer substrate supercell matrix |
| `eta` | — | ‖ε‖_F, Frobenius norm of Green–Lagrange strain tensor |
| `eps_11` | — | Normal strain along first supercell basis vector |
| `eps_22` | — | Normal strain along second supercell basis vector |
| `eps_12` | — | Shear strain (related to misfit dislocation density) |
| `area (A2)` | Å² | Substrate supercell area |
| `L0_feasible` | bool | True if Level 0 permits exact commensurability |

---

## Worked Examples

### Example 1: False-positive rejection — square MOF on hexagonal substrate

HKUST-1 (square, γ = 90°) on Au(111) (hexagonal, γ = 120°).  Level 0
immediately reports a symmetry prohibition (discriminant ratio 4/3, not a
perfect square).  The algorithm continues to Level 1/2 to bound the
approximate match, flagging all results with `L0_feasible = False`.

```python
df = matcher.run(SUBSTRATE_DB["Au_111"], hkust1)
# L0_feasible is False for all rows
# best approximate match has large area (>7000 Å²), confirming impracticality
```

### Example 2: Commensurate match — square MOF on square substrate

HKUST-1 on Au(100) (both square).  Level 0 passes (ratio = 1).  The best
match is at θ ≈ 51° with η ≈ 0.0021, corresponding to a (1, 4) / (−4, 1)
rotated supercell with biaxial strain ε₁₁ = ε₂₂ ≈ 0.15%.

```python
df = matcher.run(SUBSTRATE_DB["Au_100"], hkust1)
# L0_feasible = True, eta ~ 0.002
```

### Example 3: Custom substrate or MOF

```python
my_substrate = Lattice2D(a=3.15, b=3.15, gamma_deg=120.0, label="Custom hex")
my_mof       = Lattice2D(a=12.3, b=15.6, gamma_deg=90.0,  label="MOF-X(100)")
df = EpitaxyMatcher().run(my_substrate, my_mof)
```

---

## Limitations

The current implementation covers **geometric commensurability only**.
A complete epitaxial feasibility assessment additionally requires:

1. **Surface termination**: which crystal face is exposed is approximated here
   by the BFDH rule (`best_surface_lattice()`); DFT surface-energy calculations
   are needed for definitive assignments.
2. **Surface reconstruction**: real substrates (e.g. Au(111) 22×√3) differ
   from bulk-truncated parameters in `SUBSTRATE_DB`.
3. **Chemical compatibility**: geometric matching does not guarantee chemical
   adhesion.  Functional-group / substrate affinity must be assessed separately.

---

## Visualisation

`miqrocal.visualize` provides five functions for communicating results to
experimentalists.  All single-panel functions accept an optional `ax` keyword
so they can be embedded in custom figures.

| Function | Panel | Content |
|---|---|---|
| `plot_phi_curve(l1)` | B | Φ(θ) curve with peak annotations |
| `plot_lattice_overlay(lat_sub, lat_mof, match)` | A | Real-space substrate + MOF lattice points + supercell parallelogram |
| `plot_leed_pattern(lat_sub, lat_mof, match)` | C | Simulated LEED: substrate (blue), MOF (red), superstructure spots (green ★) |
| `plot_strain_ellipse(match)` | D | Unit circle vs. deformed ellipse, principal-strain arrows, ε component table |
| `plot_match_card(lat_sub, lat_mof, l1, match)` | E | 2×2 figure combining all four panels |
| `generate_pdf_report(lat_sub, lat_mof, l1, matches)` | — | Two-page PDF: text summary (p.1) + match-card figure (p.2) |

### Quick usage

```python
import matplotlib.pyplot as plt
from miqrocal import level1, level2, SUBSTRATE_DB, Lattice2D
from miqrocal import plot_match_card

lat_sub = SUBSTRATE_DB["Au_100"]
lat_mof = Lattice2D(18.62, 18.62, 90.0, "HKUST-1(010)")

l1      = level1.compute(lat_sub, lat_mof, G_cutoff=8.0, sigma=0.4)
matches = level2.find_matches(lat_sub, lat_mof,
                              theta_deg=l1.theta_peaks[0],
                              eta_tol=0.04)
best    = min(matches, key=lambda m: m.eta)

fig = plot_match_card(lat_sub, lat_mof, l1, best, save_path="match.png")
plt.close(fig)
```

Running `python run.py` automatically saves match-card PNGs to the `output/`
directory for all three demo cases.

---

## File Structure

```
miqrocal/
  __init__.py      public API  (__version__)
  lattice.py       Lattice2D dataclass and SUBSTRATE_DB
  level0.py        quadratic-form discriminant check
  level1.py        reciprocal-space coincidence function Phi(theta)
  level2.py        LLL lattice reduction and Green-Lagrange strain tensor
  matcher.py       EpitaxyMatcher pipeline and MatcherConfig
  batch.py         BatchConfig, batch_run (serial + parallel)
  cli.py           miqrocal CLI entry point
  visualize.py     five plotting functions (Plans A-E)
tests/
  test_lattice.py  unit tests for Lattice2D and SUBSTRATE_DB
  test_level0.py   unit tests for Level-0 algebraic check
  test_cif_parser.py  tests for centering/screw-glide extinctions
  test_batch.py    tests for batch_run (serial + parallel)
docs/
  lattice_matching_theory.md   full mathematical derivation
output/            generated figures (created at runtime)
run.py             demonstration entry point
pyproject.toml     package metadata
CHANGELOG.md       version history
LICENSE            MIT license
```

---

## References

We gratefully acknowledge the foundational works upon which miqrocal is built:

\[1\] I. Stassen, N. Burtch, A. Talin, P. Falcaro, M. Allendorf, and R. Ameloot,
"An updated roadmap for the integration of metal–organic frameworks with
electronic devices and chemical sensors,"
*Chem. Soc. Rev.* **46**, 3185–3241 (2017).
https://doi.org/10.1039/C7CS00122C

\[2\] O. Shekhah, J. Liu, R. A. Fischer, and C. Wöll,
"MOF thin films: existing and future applications,"
*Chem. Soc. Rev.* **40**, 1081–1106 (2011).
https://doi.org/10.1039/C0CS00147C

\[3\] A. Zur and T. C. McGill,
"Lattice match: An application to heteroepitaxy,"
*J. Appl. Phys.* **55**, 378–386 (1984).
https://doi.org/10.1063/1.333084
— *The seminal brute-force supercell enumeration method that miqrocal supersedes.*

\[4\] The Jacobi–Anger expansion  $e^{z\cos\varphi} = \sum_n I_n(z)\,e^{in\varphi}$
converts the two-dimensional coincidence integral into an efficiently
evaluable Fourier series; $I_n$ are modified Bessel functions of the first kind.

\[5\] G. N. Watson,
*A Treatise on the Theory of Bessel Functions*, 2nd ed.,
Cambridge University Press (1944), §2.22.
— *Standard reference for Bessel-function identities used in Level 1.*

\[6\] A. K. Lenstra, H. W. Lenstra Jr., and L. Lovász,
"Factoring polynomials with rational coefficients,"
*Math. Ann.* **261**, 515–534 (1982).
https://doi.org/10.1007/BF01457454
— *The LLL algorithm: guarantees a near-shortest basis in polynomial time.
  miqrocal uses this to find the optimal integer supercell matrix pair (M, N)
  without exhaustive enumeration.*

\[7\] J. D. H. Donnay and D. Harker,
"A new law of crystal morphology extending the law of Bravais,"
*Am. Mineral.* **22**, 446–467 (1937).
— *The Donnay–Harker extension of the Bravais–Friedel morphological rule;
  the basis of the BFDH face-ranking used in `bfdh_faces()` and
  `best_surface_lattice()`.*

\[8\] S. S.-Y. Chui, S. M.-F. Lo, J. P. H. Charmant, A. G. Orpen, and I. D. Williams,
"A Chemically Functionalizable Nanoporous Material [Cu₃(TMA)₂(H₂O)₃]ₙ,"
*Science* **283**, 1148–1150 (1999).
https://doi.org/10.1126/science.283.5405.1148
— *First report of HKUST-1 (MOF-199); crystal structure used in
  `examples/HKUST-1.cif`.*

---

## License

MIT License — see [LICENSE](LICENSE).
