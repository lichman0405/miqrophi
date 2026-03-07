# CIF to 2D Surface Lattice: Crystallographic Algorithm

> This document derives the mathematical foundations of `miqrophi.cif_parser`
> — how a 3D crystal structure in CIF format is projected onto a 2D surface
> unit cell for a given set of Miller indices (hkl).

---

## 1  Problem Statement

A CIF file encodes a three-dimensional periodic crystal via six **unit-cell
parameters**: the lengths a, b, c (Å) of the three basis vectors and the
angles α, β, γ (°) between them.  The surface exposed in an epitaxy
experiment is the (hkl) crystal face — a family of parallel planes whose
periodically repeating 2D structure is what contacts the substrate.

The goal is to extract the 2D lattice parameters (a₂D, b₂D, γ₂D) of this
surface from the bulk CIF, without any external crystallography library.

---

## 2  3D Cartesian Lattice Matrix

### 2.1  Standard Crystallographic Convention

From the unit-cell parameters (a, b, c, α, β, γ), the three basis vectors are
placed in Cartesian space according to the **right-handed standard setting**
adopted by IUCr:

$$\mathbf{a}_1 = \begin{pmatrix} a \\ 0 \\ 0 \end{pmatrix}, \qquad
  \mathbf{a}_2 = \begin{pmatrix} b\cos\gamma \\ b\sin\gamma \\ 0 \end{pmatrix}, \qquad
  \mathbf{a}_3 = \begin{pmatrix} c_x \\ c_y \\ c_z \end{pmatrix}$$

where the components of **a**₃ are determined by the three remaining metric
conditions $\mathbf{a}_i \cdot \mathbf{a}_j = \text{const}$:

$$c_x = c\cos\beta$$

$$c_y = c\,\frac{\cos\alpha - \cos\beta\cos\gamma}{\sin\gamma}$$

$$c_z = \sqrt{c^2 - c_x^2 - c_y^2}$$

**Derivation of c_y.**  The condition $\mathbf{a}_2 \cdot \mathbf{a}_3 = bc\cos\alpha$ gives

$$b\cos\gamma \cdot c_x + b\sin\gamma \cdot c_y = bc\cos\alpha$$

Substituting $c_x = c\cos\beta$ and dividing by $b$ yields the formula above.

### 2.2  Why This Convention?

Placing **a**₁ along x and **a**₂ in the xy-plane leaves **a**₃ as the most
general vector.  This choice:

- Uniquely determines the matrix from (a,b,c,α,β,γ) without additional
  degrees of freedom (no arbitrary rotation of the whole crystal).
- Equals the convention used in VESTA, CrystalMaker, and most visualisation
  software, ensuring consistent orientations.
- Makes the 2×2 surface matrix a simple sub-block when hkl = (0,0,1).

The resulting **3×3 Cartesian lattice matrix** is

$$A_3 = \begin{pmatrix} — \mathbf{a}_1 — \\ — \mathbf{a}_2 — \\ — \mathbf{a}_3 — \end{pmatrix}$$

with rows equal to the Cartesian basis vectors.  The cell volume is
$V = |\det A_3|$.

---

## 3  Reciprocal Lattice and the In-Plane Condition

### 3.1  3D Reciprocal Basis

The **reciprocal lattice basis** vectors $\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3$
are uniquely defined by the biorthogonality relation:

$$\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\,\delta_{ij}$$

In matrix form, if $A_3$ has rows $\mathbf{a}_i$, the reciprocal matrix $B_3$
(rows $\mathbf{b}_j$) satisfies

$$A_3 \, B_3^{\,T} = 2\pi\,I \quad\Longrightarrow\quad B_3 = 2\pi\,(A_3^{-1})^T$$

### 3.2  Miller Planes and the Normal Vector

The (hkl) family of planes is defined by all lattice points $\mathbf{r}$ with

$$h\,m_1 + k\,m_2 + l\,m_3 = n, \quad n \in \mathbb{Z}$$

where $\mathbf{r} = m_1\mathbf{a}_1 + m_2\mathbf{a}_2 + m_3\mathbf{a}_3$.
The **reciprocal lattice vector**

$$\mathbf{G}_{hkl} = h\,\mathbf{b}_1 + k\,\mathbf{b}_2 + l\,\mathbf{b}_3$$

is **normal** to the (hkl) plane family.  Its magnitude gives the
interplanar spacing:

$$\boxed{d_{hkl} = \frac{2\pi}{|\mathbf{G}_{hkl}|}}$$

### 3.3  In-Plane Vector Condition

A lattice translation $\mathbf{v} = m_1\mathbf{a}_1 + m_2\mathbf{a}_2 + m_3\mathbf{a}_3$
lies **within** the (hkl) surface plane if and only if it is perpendicular to
the plane normal $\mathbf{G}_{hkl}$:

$$\mathbf{G}_{hkl} \cdot \mathbf{v} = 0$$

Expanding via biorthogonality $\mathbf{b}_i \cdot \mathbf{a}_j = 2\pi\,\delta_{ij}$:

$$\mathbf{G}_{hkl}\cdot\mathbf{v}
  = (h\mathbf{b}_1 + k\mathbf{b}_2 + l\mathbf{b}_3)
    \cdot (m_1\mathbf{a}_1 + m_2\mathbf{a}_2 + m_3\mathbf{a}_3)
  = 2\pi(h m_1 + k m_2 + l m_3)$$

Setting this to zero:

$$\boxed{h\,m_1 + k\,m_2 + l\,m_3 = 0}$$

This single linear Diophantine equation in three integers defines the 2D
sublattice of bulk lattice vectors that lie within the (hkl) surface plane.
The implementation enumerates all $(m_1, m_2, m_3)$ with $|m_i| \leq N_{\max}$
that satisfy this condition, computes the Cartesian vector
$\mathbf{v} = m_1\mathbf{a}_1 + m_2\mathbf{a}_2 + m_3\mathbf{a}_3$, and sorts
the results by length $|\mathbf{v}|$.

**Example — (001) face of HKUST-1:**
The condition $0\cdot m_1 + 0\cdot m_2 + 1\cdot m_3 = 0$ forces $m_3 = 0$,
so all in-plane vectors are of the form $m_1\mathbf{a}_1 + m_2\mathbf{a}_2$.
The two shortest are $(1,0,0)$ and $(0,1,0)$ with lengths $a = b = 26.34$ Å and
angle $\gamma_{2D} = 90°$, consistent with the cubic face.

**Example — (111) face of HKUST-1 (Fm̄3m):**
The condition $m_1 + m_2 + m_3 = 0$; the shortest solutions with $|m_i|\leq 1$
are $(1,-1,0),\,(0,1,-1),\,(1,0,-1)$ and symmetry equivalents, giving
$|\mathbf{v}| = a/\sqrt{2} \approx 18.63$ Å with $\gamma_{2D} = 60°$
(hexagonal face).

---

## 4  Shortest 2D Basis Selection

### 4.1  The 2D Primitive Cell Problem

The in-plane enumeration yields a finite set of candidate vectors.  We need
the two **shortest linearly independent** ones that form the primitive 2D
unit cell — the smallest repeating surface tile.

This is the 2D version of the **lattice basis reduction** problem.  In 2D
the optimal reduced basis satisfies the **Gauss–Lagrange conditions**:

$$|\mathbf{v}_1| \leq |\mathbf{v}_2| \leq |\mathbf{v}_1 - \mathbf{v}_2|
  \quad\text{and}\quad
  |\mathbf{v}_1| \leq |\mathbf{v}_1 + \mathbf{v}_2|$$

which are equivalent to requiring $|\mathbf{v}_1 \cdot \hat{\mathbf{v}}_2| \leq
\tfrac{1}{2}|\mathbf{v}_1|$ (Lagrange–Gauss reduction).

### 4.2  Implementation Strategy

Rather than running the Gauss–Lagrange algorithm explicitly, the code:

1. Sorts all in-plane candidates by $|\mathbf{v}|$ (ascending).
2. Takes $\mathbf{v}_1$ as the shortest non-zero vector.
3. Scans subsequent candidates for the first $\mathbf{v}_2$ that is
   **linearly independent** of $\mathbf{v}_1$, tested via
   $|\mathbf{v}_1 \times \mathbf{v}_2| > \epsilon\,|\mathbf{v}_1||\mathbf{v}_2|$.

For the modest $N_{\max} = 6$ used in practice (testing 13³ − 1 ≈ 2000
integer triples), this greedy approach reliably returns the Gauss-reduced
basis for all standard crystal systems encountered in MOF epitaxy.  If the
surface unit cell is unusually large, increasing $N_{\max}$ provides more
candidates.

### 4.3  Output

The two selected vectors $\mathbf{v}_1, \mathbf{v}_2$ define the 2D lattice:

$$a_{2D} = |\mathbf{v}_1|, \quad
  b_{2D} = |\mathbf{v}_2|, \quad
  \gamma_{2D} = \arccos\!\left(\frac{\mathbf{v}_1\cdot\mathbf{v}_2}{a_{2D}\,b_{2D}}\right)$$

These are assembled into a `Lattice2D` object and forwarded to the
three-level matching pipeline.

---

## 5  Systematic Absences: Theory and Implementation

### 5.1  Lattice-Centering Extinctions

A non-primitive Bravais lattice has extra lattice points within the
conventional unit cell, creating additional translational symmetry.  In
reciprocal space, this translational symmetry forces the structure factor to
be **zero** (complete destructive interference) for all (hkl) that do not
satisfy the **centering selection rule**.

The selection rules follow directly from the Fourier sum over the centering
translations $\{\boldsymbol{\tau}_j\}$:

$$F_{hkl} \propto \sum_j e^{2\pi i (h\tau_{j1} + k\tau_{j2} + l\tau_{j3})}$$

The sum vanishes unless $h\tau_{j1} + k\tau_{j2} + l\tau_{j3} \in \mathbb{Z}$
for all centering translations simultaneously.

| Bravais type | Centering translations | Reflection condition |
|--------------|----------------------|---------------------|
| **P** | — | All *hkl* present |
| **I** (body) | $(\tfrac{1}{2},\tfrac{1}{2},\tfrac{1}{2})$ | $h+k+l = 2n$ |
| **F** (face) | $(\tfrac{1}{2},\tfrac{1}{2},0),\;(\tfrac{1}{2},0,\tfrac{1}{2}),\;(0,\tfrac{1}{2},\tfrac{1}{2})$ | $h,k,l$ all odd or all even |
| **A** | $(0,\tfrac{1}{2},\tfrac{1}{2})$ | $k+l = 2n$ |
| **B** | $(\tfrac{1}{2},0,\tfrac{1}{2})$ | $h+l = 2n$ |
| **C** | $(\tfrac{1}{2},\tfrac{1}{2},0)$ | $h+k = 2n$ |
| **R** (obverse) | $(\tfrac{2}{3},\tfrac{1}{3},\tfrac{1}{3}),\;(\tfrac{1}{3},\tfrac{2}{3},\tfrac{2}{3})$ | $-h+k+l = 3n$ |

**F case in detail.**  The four centering translations (including the origin)
give

$$\sum_{j=0}^{3} e^{2\pi i(\cdots)} = 1 + e^{i\pi(h+k)} + e^{i\pi(h+l)} + e^{i\pi(k+l)}$$

This equals 4 when $h,k,l$ have mixed parity (e.g., h odd, k even) the
exponentials alternate ±1 and the sum is zero.  Hence the F condition.

### 5.2  Zone-Specific Extinctions: Screw Axes and Glide Planes

**Screw axis** $n_p$ (n-fold rotation + translation p/n along the axis):
the structure factor for reflections along that axis involves a phase factor
$\exp(2\pi i \cdot p \cdot l / n)$.  Constructive interference requires

$$e^{2\pi i \cdot p(l-l') / n} = 1 \quad\Longrightarrow\quad l = \frac{n}{p}\,\mathbb{Z}$$

For example, the **2₁ axis along c** (n=2, p=1) forces $l = 2n$ in the 00l
zone.  The **4₁ axis** (n=4, p=1) forces $l = 4n$.

**Glide plane** with translation $\boldsymbol{\tau}$: affects reflections in
the **zone** perpendicular to the mirror plane.  A c-glide perpendicular to
**a** contributes a phase $\exp(2\pi i \cdot l/2)$ to the h0l reflections;
the condition $l = 2n$ must hold for the reflection to survive.

| Symmetry element | Affected zone | Condition |
|-----------------|--------------|-----------|
| 2₁ along **a** | h00 | $h = 2n$ |
| 2₁ along **b** | 0k0 | $k = 2n$ |
| 2₁ along **c** | 00l | $l = 2n$ |
| 4₁ or 4₃ along **c** | 00l | $l = 4n$ |
| 6₁ or 6₅ along **c** | 00l | $l = 6n$ |
| *a*-glide ⊥ **b** | h0l | $h = 2n$ |
| *b*-glide ⊥ **a** | 0kl | $k = 2n$ |
| *c*-glide ⊥ **a** | 0kl | $l = 2n$ |
| *c*-glide ⊥ **b** | h0l | $l = 2n$ |
| *n*-glide ⊥ **c** | hk0 | $h+k = 2n$ |
| *d*-glide ⊥ **c** | hk0 | $h+k = 4n$ |

### 5.3  Implementation

The centering type is extracted from the CIF via:
1. The Hermann–Mauguin symbol (`_symmetry_space_group_name_H-M`) — first letter
   (P, I, F, A, B, C, R, H).
2. Fallback: the ITA space-group number (`_symmetry_Int_Tables_number`) looked
   up in a precompiled table of all 230 space groups.
3. If neither tag is present: treated as primitive P (conservative — may
   include some absent reflections, but never excludes present ones).

Zone-specific screw/glide conditions are handled through `_SG_ZONE_CONDS`, a
precompiled dictionary mapping ~55 commonly occurring MOF space groups (P2₁/c,
P2₁2₁2₁, Pnma, Pbca, Fd̄3m, …) to their (zone, rule) pairs.

**Why only ~55 space groups?**  There are 230 space groups, but the vast
majority of experimentally observed MOFs belong to a small subset.  Rather
than an exhaustive lookup, the code covers:

- All most common monoclinic SGs (including P2₁/c, No. 14, which accounts for
  ~35% of all organic crystal structures).
- All most common orthorhombic SGs (Pnma, Pbca, P2₁2₁2₁).
- F-centred cubic SGs relevant for isoreticular MOFs (Fm̄3m, Fd̄3m).
- Common hexagonal SGs (P6₃/m, P6₃mc, etc.).

Space groups absent from `_SG_ZONE_CONDS` have no additional zone conditions
beyond centering (all reflections allowed subject to centering rules only).

---

## 6  BFDH Face Pre-selection

### 6.1  Physical Basis

The **Bravais–Friedel–Donnay–Harker (BFDH) rule** ranks crystal faces by
their interplanar spacing $d_{hkl}$: larger $d_{hkl}$ → slower growth rate
→ larger crystal face → surface most likely exposed in experiment.

The interplanar spacing follows directly from the reciprocal lattice vector
magnitude (Section 3.2):

$$d_{hkl} = \frac{2\pi}{|\mathbf{G}_{hkl}|}
          = \frac{2\pi}{\sqrt{h^2|\mathbf{b}_1|^2 + k^2|\mathbf{b}_2|^2 + l^2|\mathbf{b}_3|^2
                              + 2hk\,\mathbf{b}_1\!\cdot\!\mathbf{b}_2
                              + 2hl\,\mathbf{b}_1\!\cdot\!\mathbf{b}_3
                              + 2kl\,\mathbf{b}_2\!\cdot\!\mathbf{b}_3}}$$

### 6.2  Extinction Filtering

Before ranking, all (hkl) faces that are **systematically absent** are
discarded — a face with zero structure factor cannot drive epitaxial
nucleation even if its $d_{hkl}$ is large.  Both centering extinctions
(Section 5.1) and zone-specific extinctions (Section 5.2) are applied.

This filtering is critical for F-centred cubic MOFs such as HKUST-1 (Fm̄3m):
without it, the nominal top face (1,0,0) with d = 26.34 Å would be selected,
but this face is absent (mixed parity) and carries no Bragg intensity.  The
true top face is (1,1,1) with d = 15.21 Å.

### 6.3  Friedel Symmetry Deduplication

Bragg's law in centrosymmetric crystals (Friedel's law) equates the intensity
of (hkl) and (h̄k̄l̄): $I_{hkl} = I_{\bar{h}\bar{k}\bar{l}}$.  To avoid
listing both, only the **canonical representative** of each Friedel pair is
retained — the one whose first non-zero index is positive.

### 6.4  Output

`bfdh_faces(cif_path, n_faces, hkl_max)` returns a list of
$(\text{hkl}, d_{hkl})$ pairs, ranked by $d_{hkl}$ descending.
`best_surface_lattice(cif_path, n_faces)` chains `bfdh_faces()` with
`surface_lattice()` to return both the structural ranking and the
corresponding `Lattice2D` objects in a single call.

---

## 7  Complete Workflow

```
CIF file
   │
   ├─ _parse_cell()          Read (a, b, c, α, β, γ) from CIF text
   │
   ├─ _cartesian_matrix()    Build 3×3 Cartesian matrix A₃
   │     a₁ ∥ x,  a₂ in xy,  a₃ general
   │
   ├─ [BFDH path]            _recip_3d() → G_hkl for each (hkl)
   │   bfdh_faces()          d_hkl = 2π/|G_hkl|
   │                         filter: _is_allowed() + _is_zone_allowed()
   │                         dedup: _canonical_hkl()
   │                         → ranked [(hkl, d_hkl)] list
   │
   ├─ _in_plane_vectors()    Enumerate (m₁,m₂,m₃): h·m₁+k·m₂+l·m₃=0
   │                         Compute Cartesian v = Σ mᵢ aᵢ
   │                         Sort by |v|
   │
   ├─ _pick_two_shortest()   Select first two linearly independent vectors
   │                         |v₁×v₂| > ε |v₁||v₂|
   │
   └─ Lattice2D(a, b, γ)     a₂D=|v₁|, b₂D=|v₂|, γ₂D=∠(v₁,v₂)
         → EpitaxyMatcher
```

---

## 8  Numerical Stability Notes

**c_z guard.**  In pathological cases (extreme cell angles or nearly
co-planar bases) $c_z^2 = c^2 - c_x^2 - c_y^2$ can become slightly negative
due to floating-point rounding.  The code clamps to `max(cz_sq, 0.0)` before
taking the square root.

**cos guard.**  The angle $\gamma_{2D}$ is computed as
$\arccos(\mathbf{v}_1 \cdot \mathbf{v}_2 / (a_{2D} b_{2D}))$.  Cancellation
errors can push the argument fractionally outside $[-1, 1]$; the code clamps
it before calling `arccos`.

**n_max sensitivity.**  For very high-index surfaces (|h|+|k|+|l| ≥ 5), or
for lattices with large unit cells and the shortest in-plane vector involving
large integer coefficients, the default $N_{\max} = 6$ may miss the shortest
vector.  Raising to $N_{\max} = 9$ or 12 is safe but increases enumeration
cost as $O(N_{\max}^3)$.
