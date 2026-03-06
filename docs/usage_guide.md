# miqrocal Usage Guide

This guide explains, step by step, how to go from a crystal structure file
(CIF) all the way to an experimentally interpretable epitaxy prediction, using
HKUST-1 on Au(100) as a worked example.

---

## Table of Contents

1. [What the algorithm actually does](#1-what-the-algorithm-actually-does)
2. [How to describe a lattice (inputs)](#2-how-to-describe-a-lattice-inputs)
3. [Running the pipeline](#3-running-the-pipeline)
4. [Reading the output table](#4-reading-the-output-table)
5. [Reading the match-card figure](#5-reading-the-match-card-figure)
6. [Starting from a CIF file](#6-starting-from-a-cif-file)
7. [Common questions from experimentalists](#7-common-questions-from-experimentalists)
8. [Substrate database](#8-substrate-database)
9. [Choosing a surface plane (Miller indices)](#9-choosing-a-surface-plane-miller-indices)
10. [Batch pipeline (`batch_run`)](#10-batch-pipeline-batch_run)

---

## 1. What the algorithm actually does

The algorithm answers one question:

> If I put MOF crystal X on substrate Y, at what rotation angle θ do the
> lattices match well enough for a coherent thin film to grow?

It does this in three stages, each acting as a filter:

```
CIF file                        substrate DB
   │                                │
   ▼                                ▼
 Lattice2D ──────────────────► Lattice2D
                                    │
                           ┌────────┴────────┐
                           │   Level 0        │  O(1)
                           │   Discriminant   │  Algebraic prohibition check
                           │   check          │  (can the two symmetries ever
                           └────────┬────────┘   match exactly?)
                                    │
                           ┌────────┴────────┐
                           │   Level 1        │  O(N² + K log K)
                           │   Φ(θ) curve     │  Which rotation angles θ give
                           │                  │  maximum reciprocal-space overlap?
                           └────────┬────────┘
                                    │  candidate θ values
                           ┌────────┴────────┐
                           │   Level 2        │  O(d⁴ log B)
                           │   LLL reduction  │  Find the actual integer supercell
                           │   + strain ε     │  matrices M, N and measure strain
                           └────────┬────────┘
                                    │
                              Output table
                              + match-card figure
```

---

## 2. How to describe a lattice (inputs)

A 2D surface lattice is fully specified by three numbers:

| Parameter | Symbol | Unit | Meaning |
|-----------|--------|------|---------|
| `a`       | a      | Å    | Length of first basis vector |
| `b`       | b      | Å    | Length of second basis vector |
| `gamma_deg` | γ  | deg  | Angle between the two basis vectors |

The first basis vector always points along the x-axis.  The second is at
angle γ from it.

### Common lattice types

| System | γ | Example |
|--------|----|---------|
| Square / rectangular | 90° | Au(100), Si(001) |
| Hexagonal | 120° | Au(111), graphene |
| Oblique | any | low-symmetry MOFs |

### How to find a, b, γ for your MOF

**Option A — from a CIF file (recommended)**

```python
from miqrocal.cif_parser import surface_lattice

# Which face of the crystal is exposed?
# For HKUST-1 grown on Au(100), experimentally the (001) face is observed.
lat_mof = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
print(lat_mof)
# Lattice2D(a=26.343, b=26.343, gamma_deg=90.0, label='... (001)')
```

**Option B — from published lattice parameters**

```python
from miqrocal import Lattice2D

# HKUST-1 (001) surface: a = b = 26.343 Å, γ = 90°
lat_mof = Lattice2D(a=26.343, b=26.343, gamma_deg=90.0, label="HKUST-1(001)")
```

**Option C — from a VESTA or Mercury screenshot**

Measure the two repeat distances along the surface unit cell edges and the
angle between them using VESTA's bond/length tool.

---

## 3. Running the pipeline

```python
from miqrocal import EpitaxyMatcher, MatcherConfig, Lattice2D, SUBSTRATE_DB
from miqrocal.cif_parser import surface_lattice

# ── Define the MOF surface ──────────────────────────────────────────────────
lat_mof = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))

# ── Choose the substrate ────────────────────────────────────────────────────
lat_sub = SUBSTRATE_DB["Au_100"]   # Au(100): a = b = 4.08 Å, γ = 90°

# ── Configure the matcher ───────────────────────────────────────────────────
cfg = MatcherConfig(
    sigma   = 0.4,    # peak width in Phi(theta); larger = more forgiving
    eta_tol = 0.04,   # strain tolerance; 0.04 means ≤ 4% strain accepted
)

# ── Run ─────────────────────────────────────────────────────────────────────
matcher = EpitaxyMatcher(cfg)
df      = matcher.run(lat_sub, lat_mof)
print(df)
```

---

## 4. Reading the output table

```
 theta (deg)      eta  eps_11  eps_22  eps_12  area (A2)  L0_feasible
       50.76 0.002141 0.00151 0.00151    -0.0     5876.2         True
       90.00 0.020169 0.01426 0.01426     0.0     1348.4         True
```

| Column | What it means | How to use it |
|--------|---------------|---------------|
| `theta (deg)` | Rotation angle of the MOF unit cell relative to the substrate | Set sample azimuth to this angle in the deposition chamber |
| `eta` | Frobenius norm of the strain tensor ‖ε‖_F (dimensionless) | Rule of thumb: η < 0.02 → low strain; η < 0.05 → still feasible |
| `eps_11` | Normal strain along first supercell axis | Positive = tensile, negative = compressive |
| `eps_22` | Normal strain along second supercell axis | Same convention as eps_11 |
| `eps_12` | Shear strain | Non-zero → unit cell shape changes, misfit dislocations expected |
| `area (A2)` | Substrate supercell area in Å² | Larger → fewer coincidence points per unit area → harder to nucleate |
| `L0_feasible` | True if exact commensurability is algebraically possible | False = approximate match only; real film will always have misfit dislocations |

### Practical decision guide

```
L0_feasible = True  AND  eta < 0.02  AND  area < 2000 Å²
       └─────── ideal case: expect coherent thin film growth ───────┘

L0_feasible = True  AND  0.02 < eta < 0.05
       └─── near-coherent: growth possible, some misfit dislocations ──┘

L0_feasible = False  OR  eta > 0.05
       └──── use as approximate prediction only; likely polycrystalline ┘
```

---

## 5. Reading the match-card figure

```python
import matplotlib.pyplot as plt
from miqrocal import level1, level2, plot_match_card, SUBSTRATE_DB
from miqrocal.cif_parser import surface_lattice

lat_sub = SUBSTRATE_DB["Au_100"]
lat_mof = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
cfg     = MatcherConfig(sigma=0.4, eta_tol=0.04)

l1      = level1.compute(lat_sub, lat_mof, G_cutoff=8.0, sigma=cfg.sigma)
matches = level2.find_matches(lat_sub, lat_mof,
                              theta_deg=l1.theta_peaks[0],
                              eta_tol=cfg.eta_tol)
best    = min(matches, key=lambda m: m.eta)

fig = plot_match_card(lat_sub, lat_mof, l1, best, save_path="match.png")
plt.close(fig)
```

The 2×2 figure contains:

```
┌─────────────────────────┬─────────────────────────┐
│  TOP-LEFT: Φ(θ) curve   │  TOP-RIGHT: Real-space  │
│                         │  lattice overlay         │
│  Purple curve = overlap │                          │
│  function vs angle.     │  Blue dots = substrate   │
│  Peaks = candidate θ.   │  Red triangles = MOF     │
│  Green line = selected. │  Green box = supercell   │
├─────────────────────────┼─────────────────────────┤
│  BOTTOM-LEFT: LEED      │  BOTTOM-RIGHT: Strain   │
│  simulation             │  ellipse                 │
│                         │                          │
│  Blue circles = subst.  │  Grey dashed = reference │
│  Red triangles = MOF    │  Red solid = strained    │
│  Green stars = super-   │  Arrows = principal      │
│  structure reflections  │  strain directions       │
└─────────────────────────┴─────────────────────────┘
```

**What to tell your experimentalist colleague:**

- The **LEED panel** shows exactly what pattern to expect on the LEED screen,
  including the extra superstructure spots (green stars).  Count and measure
  the green-star positions to verify the prediction experimentally.

- The **real-space overlay** shows how many substrate unit cells are covered by
  one MOF unit cell.  A small, compact parallelogram = easier nucleation.

- The **strain ellipse** — if it looks like a circle, strain is small and
  isotropic (best case).  If one arrow is much longer than the other, the film
  is stretched more in one direction (check eps_11 vs eps_22).

---

## 6. Starting from a CIF file

### Step 1 — obtain a CIF

Sources for MOF CIF files:

| Source | URL | Notes |
|--------|-----|-------|
| Cambridge Structural Database (CSD) | https://www.ccdc.cam.ac.uk/ | Largest collection; subscription required for full access |
| Crystallography Open Database (COD) | https://www.crystallography.net/cod/ | Free, CC0 licensed |
| Materials Project | https://materialsproject.org/ | Free API; computed structures |
| CoRE MOF database | https://zenodo.org/record/3370144 | ~10 000 curated MOF CIFs |
| MOFXDB | https://mof.tech.northwestern.edu/ | Northwestern curated set |

### Step 2 — choose a surface plane (Miller indices)

The (hkl) plane to expose is determined by the growth chemistry.  Common
choices:

| MOF family | Typical facet | Reason |
|-----------|---------------|--------|
| HKUST-1 (cubic, Fm-3m) | (111) | F-centering: `{100}`/`{010}`/`{001}` all absent; `{111}` is the true BFDH top face |
| UiO-66 (cubic, Fm-3m) | (001), (111) | (001) smallest cell; (111) for hex. subs. |
| MIL-53 (orthorhombic) | (100) | Pore aperture exposed |
| ZIF-8 (cubic, I-43m) | (001) | Flat square face |

If you are unsure, check the BFDH (Bravais–Friedel–Donnay–Harker) morphology
prediction in VESTA or Mercury.

### Step 3 — extract the 2D lattice and run

```python
from miqrocal import EpitaxyMatcher, MatcherConfig, SUBSTRATE_DB
from miqrocal.cif_parser import surface_lattice, read_cell

# Inspect CIF first
cell = read_cell("my_mof.cif")
print(f"a={cell['a']:.3f}  b={cell['b']:.3f}  c={cell['c']:.3f} Å")
print(f"α={cell['alpha']:.1f}  β={cell['beta']:.1f}  γ={cell['gamma']:.1f} °")

# Extract 2D surface lattice for the (001) face
lat_mof = surface_lattice("my_mof.cif", hkl=(0, 0, 1), label="My-MOF(001)")

# Run matching against all common substrates
cfg     = MatcherConfig(sigma=0.4, eta_tol=0.05)
matcher = EpitaxyMatcher(cfg)

for name, lat_sub in SUBSTRATE_DB.items():
    df = matcher.run(lat_sub, lat_mof, verbose=False)
    if df is not None:
        best = df.iloc[0]
        print(f"{name:15s}  θ={best['theta (deg)']:6.2f}°  "
              f"η={best['eta']:.4f}  L0={best['L0_feasible']}")
```

---

## 7. Common questions from experimentalists

**Q: I got `L0_feasible = False`. Can the film still grow?**

Yes — `L0_feasible = False` means that *perfect* (zero-strain) commensurability
is algebraically impossible due to the symmetry mismatch between the two
lattices.  The film will still form provided η is small enough (< ~0.05), but
it will always contain misfit dislocations.  This is normal for most real
epitaxial systems.

---

**Q: There are many rows with similar θ values. Which one is correct?**

Multiple rows differing by ~90° (or 60° for hexagonal) are symmetry-equivalent
rotation domains.  They give the same physical result.  Report only the
smallest unique θ in (0°, 90°) or (0°, 60°) for square or hexagonal
substrates, respectively.

Rows with substantially different θ values represent genuinely different
orientational relationships.  In experiment, you may observe multiple
coexisting domains; their relative abundance depends on kinetics, not just
thermodynamics.

---

**Q: How many substrate unit cells does one MOF unit cell cover?**

This is `round(area (A2) / lat_sub.omega)`.  For HKUST-1/Au(100) at θ = 50.76°
with area = 5876 Å² and Au(100) cell area = 4.08² = 16.6 Å²:

    N_cells ≈ 5876 / 16.6 ≈ 354  substrate cells per MOF supercell

---

**Q: I increased eta_tol but got no results. What's wrong?**

Either the lattice symmetry is incompatible (check `L0_feasible` from Level 0
— it appears in the verbose output even when no match is found), or the
`sigma` parameter is too narrow.  Try:

```python
cfg = MatcherConfig(sigma=0.6, eta_tol=0.08)
```

---

**Q: Can I use this for non-MOF overlayers (organic molecules, 2D materials)?**

Yes.  Any periodic 2D overlayer with known lattice constants works.  Just
provide the correct `a, b, gamma_deg`.  For example, for a PTCDA monolayer
(a = 14.86 Å, b = 11.96 Å, γ = 90°):

```python
ptcda = Lattice2D(14.86, 11.96, 90.0, "PTCDA(102)")
```

---

## 8. Substrate database

The built-in database covers the most common epitaxy substrates:

| Key | Material | a (Å) | b (Å) | γ (°) | System |
|-----|----------|-------|-------|--------|--------|
| **Noble-metal single crystals** |||||
| `"Au_111"` | Au(111) | 2.880 | 2.880 | 120 | hexagonal |
| `"Au_100"` | Au(100) | 4.080 | 4.080 |  90 | square |
| `"Ag_111"` | Ag(111) | 2.890 | 2.890 | 120 | hexagonal |
| `"Cu_111"` | Cu(111) | 2.560 | 2.560 | 120 | hexagonal |
| `"Pt_111"` | Pt(111) | 2.775 | 2.775 | 120 | hexagonal |
| **Carbon substrates** |||||
| `"Graphene"` | Graphene | 2.46 | 2.46 | 120 | hexagonal |
| `"HOPG"` | HOPG(0001) graphite | 2.46 | 2.46 | 120 | hexagonal |
| **Oxide substrates** |||||
| `"ITO_111"` | ITO(111) | 4.130 | 4.130 | 120 | hexagonal |
| `"MgO_001"` | MgO(001) | 4.211 | 4.211 |  90 | square |
| `"SrTiO3_001"` | SrTiO₃(001) | 3.905 | 3.905 |  90 | square |
| **Semiconductor substrates** |||||
| `"Si_001"` | Si(001) | 5.431 | 5.431 |  90 | square |

For a substrate not in the database, create it manually:

```python
substrate = Lattice2D(a=3.524, b=3.524, gamma_deg=120.0, label="Ni(111)")
```

---

## 9. Choosing a surface plane (Miller indices)

### Automatic selection: `best_surface_lattice()` (recommended)

For any CIF, let the BFDH rule decide automatically:

```python
from miqrocal import best_surface_lattice

# Returns top-5 faces ranked by interplanar spacing (most prominent first)
results = best_surface_lattice("my_mof.cif", n_faces=5)
for hkl, d_hkl, lat in results:
    print(f"{str(hkl):<12}  d={d_hkl:.2f} Å  "
          f"a={lat.a:.2f}  b={lat.b:.2f}  γ={lat.gamma_deg:.1f}°")
```

The face with the largest $d_{hkl}$ is usually the one exposed in
experiment.  Use this to narrow down which faces to pass to the matcher.

**Example output — Co-sql MOF:**
```
(0, 0, 1)   d=10.62 Å   a=5.96  b=7.07  γ=92.8°   ← layer plane, top choice
(0, 1, 0)   d= 7.06 Å   a=5.96  b=10.67 γ=95.2°
(1, 0, 0)   d= 5.93 Å   a=7.07  b=10.67 γ=90.0°
```

**Example output — HKUST-1 (Fm-3m, F-centering applied):**
```
(1, 1, 1)   d=15.21 Å   a=37.26  b=37.26  γ=60.0°  ← {111} family: first BFDH face after F-centering
(1, 1, -1)  d=15.21 Å   a=37.26  b=37.26  γ=60.0°  ← symmetry-equivalent variant
(1, -1, 1)  d=15.21 Å   a=37.26  b=37.26  γ=60.0°  ← symmetry-equivalent variant
```

> **Why not `(1,0,0)`?**  HKUST-1 belongs to space group Fm̄3m (No. 225).
> For F-centering, h, k, l must all be odd or all be even.  The index set
> (1, 0, 0) is mixed parity → absent.  Same for (0,1,0) and (0,0,1).
> The true first-allowed face is (1,1,1) with d = 26.34/√3 ≈ 15.21 Å,
> consistent with STM/AFM observations of HKUST-1 on Au(111).

### Manual selection

If you know which face is exposed (e.g. from literature), specify hkl explicitly:

```python
from miqrocal.cif_parser import surface_lattice

# Square (001) face: natural choice for Au(100)
lat_001 = surface_lattice("examples/HKUST-1.cif", hkl=(0, 0, 1))
# a=26.343  b=26.343  γ=90°

# (110) face: rectangular, exposed on some stepped surfaces
lat_110 = surface_lattice("examples/HKUST-1.cif", hkl=(1, 1, 0))
# a=26.343  b=37.255  γ=90°

# (111) face: hexagonal, natural choice for Au(111)
lat_111 = surface_lattice("examples/HKUST-1.cif", hkl=(1, 1, 1))
# a=37.255  b=37.255  γ=60°
```

If the output cell seems too large (> 50 Å), the crystal may have a small
primitive cell that can be used instead.  Check `read_cell()` output for the
raw 3D cell size.

---

## 10. Batch pipeline (`batch_run`)

For screening many CIF files against multiple substrates in one call, use
`batch_run()` instead of writing a manual loop.

### Basic usage

```python
from miqrocal import batch_run, BatchConfig, MatcherConfig

df = batch_run(
    "data/*.cif",                          # glob pattern or list of paths
    config=BatchConfig(
        n_faces     = 3,                   # BFDH top-3 faces per MOF
        substrates  = ["Au_111", "Cu_111", "Au_100"],
        run_tag     = "screen_v1",         # appended to the output directory name
        outputs     = {"csv", "pdf", "png"},
        matcher_cfg = MatcherConfig(eta_tol=0.04),
    ),
)

# Filter and display matches
matched = df[df.match_found].sort_values("eta")
print(matched[["mof_name", "hkl", "substrate", "theta_deg", "eta"]].to_string())
```

### `BatchConfig` fields

| Field | Default | Description |
|-------|---------|-------------|
| `substrates` | `None` | List of `SUBSTRATE_DB` keys to screen against; `None` = all 11 substrates |
| `n_faces` | `1` | BFDH top-N faces per MOF |
| `outputs` | `{"csv","pdf","png"}` | Controls which files are written (see below) |
| `output_dir` | `"output"` | Root directory for all run output |
| `run_tag` | `""` | Human-readable suffix appended to the run directory name |
| `matcher_cfg` | `None` | Full `MatcherConfig` for fine-grained tuning; library defaults used when `None` |
| `top_n_report` | `10` | Rows shown in each PDF report |
| `verbose` | `True` | Print per-pair progress to stdout |

### Output directory structure

Each `batch_run()` call creates a self-contained, **timestamped run directory**
inside `output_dir` so that successive screening runs never overwrite each other:

```
output/
  run_2026-03-06_12-00-00_screen_v1/    ← one directory per batch_run() call
    summary.csv                          ← one row per (MOF face × substrate) pair
    HKUST-1_111__Au_111/                 ← {mof_stem}_{hkl}__{substrate}
      match_card.png
      report.pdf
    HKUST-1_111__Cu_111/
      match_card.png
      report.pdf
    Co_sql_001__Au_100/
      match_card.png
      report.pdf
    ...
```

All directory and file names use only alphanumeric characters, underscores, and
hyphens — no spaces, colons, or Windows-unsafe characters — so paths are safe
and **traceable** across all platforms.

### Controlling output content

Pass any subset of `{"csv", "pdf", "png"}` to `outputs`:

```python
# In-memory only, no files written (fastest for exploratory analysis)
df = batch_run("data/*.cif", config=BatchConfig(outputs=set()))

# Only the summary CSV (no per-pair graphics)
df = batch_run("data/*.cif", config=BatchConfig(outputs={"csv"}))

# Full output: CSV + PNG match cards + PDF reports
df = batch_run("data/*.cif", config=BatchConfig(outputs={"csv", "pdf", "png"}))
```

### Summary CSV columns

| Column | Meaning |
|--------|---------|
| `cif_file` | Basename of the source CIF |
| `mof_name` | Compound name from CIF (or file stem) |
| `hkl` | Miller indices, e.g. `"(1,1,1)"` |
| `d_hkl_A` | Interplanar spacing (Å) |
| `mof_a` / `mof_b` / `mof_gamma` | 2D surface lattice parameters |
| `substrate` | Substrate key from `SUBSTRATE_DB` |
| `match_found` | `True` if any match with η < η_tol |
| `theta_deg` | Rotation angle of best match (`NaN` if none) |
| `eta` | Strain index η of best match |
| `eps_11` / `eps_22` / `eps_12` | Green–Lagrange strain components |
| `area_A2` | Substrate supercell area (Å²) |
| `L0_feasible` | Algebraic commensurability flag |
| `png_path` | Path to PNG relative to `output_dir` |
| `pdf_path` | Path to PDF relative to `output_dir` |

### Note on systematic absences

`bfdh_faces()` applies **lattice-centering extinction rules** automatically.
The space group symbol is read from `_symmetry_space_group_name_H-M` (or
`_symmetry_Int_Tables_number` as fallback) to determine the Bravais centering
type.  The general centering conditions are then applied during face
enumeration:

| Centering | Condition for reflection to be **present** |
|-----------|-----------------------------------------|
| P | all *hkl* |
| I | *h* + *k* + *l* = 2*n* |
| F | *h*, *k*, *l* all odd **or** all even |
| A | *k* + *l* = 2*n* |
| B | *h* + *l* = 2*n* |
| C | *h* + *k* = 2*n* |
| R | −*h* + *k* + *l* = 3*n* |

Screw-axis and glide-plane extinctions are not applied (centering alone is
sufficient to correct all common cubic and non-primitive MOF space groups).
If the space group cannot be parsed, the structure is treated as P (no
extinction) to remain compatible with minimal-metadata CIFs.
