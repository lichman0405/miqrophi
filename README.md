# miqrophi

**miqrophi** (`pip install miqrophi`) is a Python library for predicting the
epitaxial compatibility between MOF (metal-organic framework) thin films and
crystalline substrates using a three-level hierarchical lattice-matching
algorithm.

MOF thin-film epitaxy has emerged as a route to integrate porous, functional
materials with device-grade substrates \[1, 2\].  Identifying epitaxially
compatible substrate–overlayer pairs is the first bottleneck in any growth
campaign.

The classical approach to this problem is the brute-force supercell enumeration
introduced by **Zur & McGill (1984)** \[3\], which iterates over all integer
matrix pairs (M, N) and scales as O(N_m² · N_s²).  miqrophi replaces that
exhaustive search with a three-level hierarchy:

| Level | Module | Method | Cost |
|-------|--------|--------|------|
| 0 | `discriminant` | Quadratic-form discriminant check | O(1) |
| 1 | `coincidence` | Reciprocal-space coincidence function Φ(θ) via Jacobi–Anger + FFT | O(N_G² + K log K) |
| 2 | `supercell` | Floating-point LLL lattice reduction → Green–Lagrange strain tensor | O(d⁴ log B) |

Surface-exposure probabilities are estimated using the **BFDH morphological
rule** \[7\] (Donnay & Harker 1937), which ranks faces by interplanar spacing
$d_{hkl}$ without requiring a force-field calculation.

For a complete mathematical derivation see
[`docs/lattice_matching_theory.md`](docs/lattice_matching_theory.md).

---

## Performance

| Metric | Value |
|--------|-------|
| Single pair (L0 + L1 + L2), warm | **~80 ms** |
| 1 MOF × 11 built-in substrates | **~900 ms** |
| L1 vectorised speedup vs. original | **4.4×** |
| L2 Numba JIT speedup (warm) | **~17×** |

Level 1 uses a sparse reciprocal-space pre-filter and exploits
$I_{-n}(z) = I_n(z)$ symmetry to halve Bessel-function evaluations.
Level 2 uses an optional Numba-compiled LLL kernel (`pip install miqrophi[fast]`).
Both fall back gracefully to pure NumPy / SciPy when Numba is absent.

---

## Requirements

- Python ≥ 3.10 (tested on 3.10 – 3.13)
- `numpy >= 1.24`
- `scipy >= 1.10`
- `pandas >= 2.0`
- `matplotlib >= 3.7`

No compiled extensions required for the base install.

---

## Installation

```bash
# base install
pip install miqrophi

# with Numba JIT acceleration for the LLL kernel (~17× speedup on L2)
pip install "miqrophi[fast]"

# development install from source
git clone https://github.com/lichman0405/miqrophi.git
cd miqrophi
pip install -e ".[fast,dev]"
```

---

## Quick Start

### Option A — automatic face selection from CIF (recommended)

```python
from miqrophi.cif_parser import best_surface_lattice
from miqrophi import EpitaxyMatcher, MatcherConfig, SUBSTRATE_DB

# BFDH rule picks the most prominent face automatically
results = best_surface_lattice("examples/HKUST-1.cif", n_faces=3)
hkl, d_hkl, lat_mof = results[0]          # largest d_hkl -> most likely face
print(hkl, f"d={d_hkl:.2f} Å")            # -> (1, 1, 1)  d=15.21 Å

# Match against all built-in substrates, pick the best
matcher = EpitaxyMatcher(MatcherConfig(eta_tol=0.05))
best_eta, best_key = 1.0, None
for key, lat_sub in SUBSTRATE_DB.items():
    df = matcher.run(lat_sub, lat_mof, verbose=False)
    if df is not None and df["eta"].iloc[0] < best_eta:
        best_eta, best_key = df["eta"].iloc[0], key
print(f"Best substrate: {best_key}  η={best_eta:.4f}")
# -> Best substrate: Graphene  η=0.0002
```

### Option B — from a specific CIF face

```python
from miqrophi import EpitaxyMatcher, MatcherConfig, SUBSTRATE_DB
from miqrophi.cif_parser import surface_lattice

lat_mof = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
# -> Lattice2D(a=26.343, b=26.343, gamma_deg=90.0)

matcher = EpitaxyMatcher(MatcherConfig(eta_tol=0.04))
df = matcher.run(lat_sub=SUBSTRATE_DB["Au_100"], lat_mof=lat_mof)
print(df)
```

### Option C — from known lattice parameters

```python
from miqrophi import EpitaxyMatcher, MatcherConfig, Lattice2D, SUBSTRATE_DB

lat_mof = Lattice2D(a=26.343, b=26.343, gamma_deg=90.0, label="HKUST-1(001)")
matcher = EpitaxyMatcher(MatcherConfig(eta_tol=0.04))
df = matcher.run(lat_sub=SUBSTRATE_DB["Au_100"], lat_mof=lat_mof)
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
from miqrophi import generate_pdf_report
from miqrophi import coincidence, supercell, SUBSTRATE_DB, Lattice2D

lat_sub = SUBSTRATE_DB["Au_100"]
lat_mof = Lattice2D(26.343, 26.343, 90.0, "HKUST-1(001)")

l1      = coincidence.compute(lat_sub, lat_mof)
matches = supercell.find_matches(lat_sub, lat_mof, theta_deg=l1.theta_peaks[0])

generate_pdf_report(lat_sub, lat_mof, l1, matches,
                    title="HKUST-1 / Au(100)",
                    save_path="output/report.pdf")
```

Run the bundled demonstration (saves PNG + PDF to `output/`):

```bash
python run.py
```

---

## Batch Processing

```python
from miqrophi import batch_run, BatchConfig

df = batch_run(
    "examples/*.cif",
    config=BatchConfig(
        substrates=["Au_111", "Au_100", "Graphene"],
        n_faces=2,
        outputs={"csv", "pdf"},
        eta_tol=0.05,           # passed through to MatcherConfig
    ),
)
print(df[["cif_file", "substrate", "eta"]].head())
```

A timestamped run directory is always created under `output_dir`;
pass `outputs=set()` for in-memory-only operation (no CSV/PNG/PDF written).

---

## API Reference

### `Lattice2D`

```python
Lattice2D(a: float, b: float, gamma_deg: float, label: str = "")
```

Represents a 2D periodic lattice.  The first basis vector lies along the x-axis.

| Attribute | Type | Description |
|-----------|------|-------------|
| `A` | `(2, 2) ndarray` | Cartesian lattice matrix (rows = basis vectors) |
| `omega` | `float` | Unit-cell area (Å²) |

`recip_matrix()` returns the 2×2 reciprocal lattice matrix B satisfying A · Bᵀ = 2π I.

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
    sigma         = 0.3,   # Φ(θ) Gaussian peak width (Å⁻¹)
    eta_tol       = 0.05,  # Green–Lagrange strain Frobenius norm tolerance
    top_theta     = 8,     # Level-1 peaks forwarded to Level 2
    lambda_values = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0],
)
```

`sigma` — Larger values broaden Φ(θ) peaks; ~0.1–0.5 × (2π / a_sub) is a
good starting range.

`eta_tol` — Matches with η < 0.02 are generally considered physically
feasible; 0.05–0.10 reveals strained approximate matches.

`lambda_values` — Each value traces one point on the strain-vs-supercell-size
Pareto front.  Small λ (0.02) favours low-strain large supercells; large λ
(1.0) favours compact supercells.

### `EpitaxyMatcher`

```python
matcher = EpitaxyMatcher(config: MatcherConfig = None)
df      = matcher.run(lat_sub, lat_mof, verbose=True) -> pd.DataFrame | None
```

Returns a `DataFrame` sorted by `eta`, or `None` when no match satisfies
`eta_tol`.

---

## Output Columns

| Column | Unit | Description |
|--------|------|-------------|
| `theta (deg)` | ° | Rotation angle of MOF supercell relative to substrate |
| `M` | — | 2×2 integer MOF supercell matrix |
| `N` | — | 2×2 integer substrate supercell matrix |
| `eta` | — | ‖ε‖_F, Frobenius norm of Green–Lagrange strain tensor |
| `eps_11` | — | Normal strain along first supercell basis vector |
| `eps_22` | — | Normal strain along second supercell basis vector |
| `eps_12` | — | Shear strain |
| `area (A2)` | Å² | Substrate supercell area |
| `L0_feasible` | bool | True if Level 0 permits exact commensurability |

---

## Worked Examples

### Example 1: False-positive rejection — square MOF on hexagonal substrate

HKUST-1 (square, γ = 90°) on Au(111) (hexagonal, γ = 120°).  Level 0
immediately reports a symmetry prohibition (discriminant ratio 4/3, not a
perfect square).  Results are returned but flagged `L0_feasible = False`.

```python
from miqrophi import EpitaxyMatcher, SUBSTRATE_DB, Lattice2D

hkust1  = Lattice2D(26.343, 26.343, 90.0, "HKUST-1(001)")
df = EpitaxyMatcher().run(SUBSTRATE_DB["Au_111"], hkust1)
# L0_feasible is False for all rows
```

### Example 2: Commensurate match — square MOF on square substrate

HKUST-1 on Au(100) (both square).  Level 0 passes (ratio = 1).  Best match
at θ ≈ 51°, η ≈ 0.0021, biaxial strain ε₁₁ = ε₂₂ ≈ 0.15%.

```python
df = EpitaxyMatcher().run(SUBSTRATE_DB["Au_100"], hkust1)
# L0_feasible = True, eta ~ 0.002
```

### Example 3: Custom substrate or MOF

```python
from miqrophi import EpitaxyMatcher, Lattice2D

my_sub = Lattice2D(3.15, 3.15, 120.0, "Custom hex")
my_mof = Lattice2D(12.3, 15.6,  90.0, "MOF-X(100)")
df = EpitaxyMatcher().run(my_sub, my_mof)
```



## Limitations

The current implementation covers **geometric commensurability only**:

1. **Surface termination** is approximated by the BFDH rule; DFT surface-energy
   calculations are needed for definitive assignments.
2. **Surface reconstruction** (e.g. Au(111) 22×√3) is not modelled.
3. **Chemical compatibility** (adhesion, functional-group affinity) must be
   assessed separately.

---

## Visualisation

`miqrophi.visualize` provides six functions:

| Function | Content |
|---|---|
| `plot_phi_curve(l1)` | Φ(θ) curve with peak annotations |
| `plot_lattice_overlay(lat_sub, lat_mof, match)` | Real-space supercell overlay |
| `plot_leed_pattern(lat_sub, lat_mof, match)` | Simulated LEED pattern |
| `plot_strain_ellipse(match)` | Principal-strain ellipse + component table |
| `plot_match_card(lat_sub, lat_mof, l1, match)` | 2×2 combined figure |
| `generate_pdf_report(...)` | Two-page PDF: text summary + match card |

```python
from miqrophi import plot_match_card, coincidence, supercell, SUBSTRATE_DB, Lattice2D
import matplotlib.pyplot as plt

lat_sub = SUBSTRATE_DB["Au_100"]
lat_mof = Lattice2D(26.343, 26.343, 90.0, "HKUST-1(001)")

l1      = coincidence.compute(lat_sub, lat_mof)
matches = supercell.find_matches(lat_sub, lat_mof, theta_deg=l1.theta_peaks[0])

fig = plot_match_card(lat_sub, lat_mof, l1, matches[0], save_path="match.png")
plt.close(fig)
```

---

## File Structure

```
miqrophi/               Python package (import as: import miqrophi)
  __init__.py           public API exports
  lattice.py            Lattice2D dataclass and SUBSTRATE_DB
  discriminant.py       Level 0 — quadratic-form symmetry pre-filter
  coincidence.py        Level 1 — reciprocal-space coincidence function Φ(θ)
  supercell.py          Level 2 — LLL lattice reduction + Green–Lagrange strain
  matcher.py            EpitaxyMatcher pipeline and MatcherConfig
  batch.py              BatchConfig, batch_run (serial + parallel)
  cif_parser.py         CIF reader, BFDH face ranking, surface lattice extraction
  cli.py                miqrophi CLI entry point
  visualize.py          plotting and PDF report generation
tests/
  test_lattice.py
  test_discriminant.py
  test_cif_parser.py
  test_batch.py
docs/
  lattice_matching_theory.md   full mathematical derivation
  usage_guide.md               user guide
examples/                      sample CIF files
output/                        generated figures (created at runtime)
run.py                         demonstration entry point
pyproject.toml
CHANGELOG.md
LICENSE
```

---

## CLI

```bash
# list all built-in substrates
miqrophi substrates

# match a single CIF against specific substrates
miqrophi run examples/HKUST-1.cif --substrates Au_111 Au_100 Graphene

# batch screen a folder of CIFs
miqrophi batch "examples/*.cif" --n-faces 2 --eta-tol 0.05 --outputs csv,pdf
```

---

## References

\[1\] Shekhah, O. et al. *Chem. Soc. Rev.* **2011**, 40, 1081.  
\[2\] Falcaro, P. et al. *Nat. Chem.* **2016**, 8, 925.  
\[3\] Zur, A.; McGill, T. C. *J. Appl. Phys.* **1984**, 55, 378.  
\[4\] Watson, G. N. *A Treatise on the Theory of Bessel Functions*, 2nd ed., 1944.  
\[5\] Lenstra, A. K.; Lenstra, H. W.; Lovász, L. *Math. Ann.* **1982**, 261, 515.  
\[6\] Donnay, J. D. H.; Harker, D. *Am. Mineral.* **1937**, 22, 446.  
\[7\] Chui, S. S.-Y. et al. *Science* **1999**, 283, 1148.
