# 2D Epitaxial Lattice Matching: Mathematical and Physicochemical Foundations

> This document derives the complete theoretical basis for the three-level
> hierarchical screening algorithm implemented in `miqrophi`.

---

## 1  Physical Problem

### 1.1  Interface Commensurability in Epitaxial Growth

When a MOF thin film grows epitaxially on a single-crystal substrate, two
two-dimensional periodic lattices are superimposed at the interface.  The
binding energy at the interface is governed by the **commensurability** of
these lattices — the existence and quality of a shared periodic supercell.

The classical Hooks–Lattice classification distinguishes three regimes:

| Type | Epitaxy matrix R | Physical meaning |
|------|-----------------|-----------------|
| **Commensurate** | R ∈ Z^{2×2} | Every overlayer atom sits on a substrate lattice point; lowest interface energy |
| **Coincident** | At least one row of R ∈ Z^{1×2} | Overlayer atoms fall on a family of substrate lattice lines |
| **Incommensurate** | R ∉ Q^{2×2} | No periodic match; governed by van der Waals forces |

The central question addressed by this algorithm is: **given the 2D lattice
of a MOF surface and a substrate, does a finite near-commensurate supercell
exist, and how large is the epitaxial strain?**

### 1.2  Why 2D Matching Cannot Be Reduced to 1D

A one-dimensional check compares only |v_MOF| ≈ |v_sub|, ignoring the
angular constraint between the two basis vectors.  For example, the diagonal
of a square lattice (γ = 90°) may coincidentally match the length of a
supercell vector of a hexagonal lattice (γ = 120°), but it is impossible to
**simultaneously** match both basis vectors, because no integer-coefficient
affine map relates 90° to 120°.  This is precisely the HKUST-1 (square) /
Au(111) (hexagonal) case: a 1D check yields a false positive, while the 2D
framework intercepts it algebraically.

### 1.3  Lattice Matrix Representation

The 2D lattice is encoded as a 2×2 real matrix A whose rows are the
Cartesian basis vectors.  By convention the first basis vector lies along x:

$$A_s = \begin{pmatrix} a_1 & 0 \\ a_2\cos\gamma & a_2\sin\gamma \end{pmatrix}$$

The unit-cell area is Ω = |det A|.  A supercell is defined by an integer
matrix M ∈ Z^{2×2}; its Cartesian basis vectors are C = M · A, and
|det M| equals the ratio of supercell to unit-cell area.

---

## 2  Architecture: Three-Level Hierarchical Screening

Each level operates in a different mathematical space and answers a different
physical question.  Only candidates passing one level are forwarded to the next.

| Level | Mathematical space | Physical question | Cost |
|-------|--------------------|-------------------|------|
| **0** | Algebraic number theory | Is a match algebraically possible at all? | O(1) |
| **1** | Reciprocal space / harmonic analysis | At which rotation angle θ* is overlap maximised? | O(N_G² + K log K) |
| **2** | 4-D integer lattice / LLL | What is the optimal integer supercell at θ*? | O(d⁴ log B) |

---

## 3  Level 0: Quadratic-Form Discriminant Check

### 3.1  Physical Motivation

Level 0 answers whether an exact-commensurate solution exists **algebraically**,
not merely numerically.  For the square-on-hexagonal case the answer is a
provable impossibility, independent of numerical tolerances.

### 3.2  Metric Tensor and Automorphism Group

All geometric information of a 2D lattice is encoded in its **metric tensor**:

$$G = \begin{pmatrix} a^2 & ab\cos\gamma \\ ab\cos\gamma & b^2 \end{pmatrix}$$

The point group of the lattice is the **automorphism group**

$$\mathrm{Aut}(G) = \{ M \in \mathbb{Z}^{2\times 2} \mid M^T G M = G \}$$

whose order encodes crystallographic symmetry: square → |Aut| = 4 (C4),
hexagonal → |Aut| = 6 (C6).

### 3.3  Algebraic Condition for Commensurability

Lattices A_s and A_o admit an exact commensurate supercell if and only if

$$\exists M, N \in \mathbb{Z}^{2\times 2} : M A_o = N A_s
\quad\iff\quad
T \equiv A_o A_s^{-1} \in \mathbb{Q}^{2\times 2}$$

**Worked example — square on hexagonal:**

$$T = A_o A_s^{-1} = \frac{a_o}{a_s}
\begin{pmatrix} 1 & 1/\sqrt{3} \\ 0 & 2/\sqrt{3} \end{pmatrix}$$

The factor 1/√3 is irrational for every rational ratio a_o/a_s;
hence T ∉ Q^{2×2} and **no integer solution exists**.

### 3.4  Discriminant Criterion

The impossibility stems from the two lattice types belonging to different
**imaginary quadratic fields**.  Each lattice type has a canonical binary
quadratic form Q(m, n) = am² + bmn + cn² with discriminant Δ = b² − 4ac:

| Lattice type | Canonical form | Δ (normalised) | Field |
|--------------|---------------|----------------|-------|
| Square (γ = 90°, a = b) | m² + n² | −4 | Q(i) |
| Hexagonal (γ = 120°, a = b) | m² − mn + n² | −3 | Q(ω), ω = e^{2πi/3} |
| Rectangular (γ = 90°, a ≠ b) | m² + (b/a)²n² | −(2b/a)² | depends on ratio |

The normalised discriminant is δ = −4 sin²γ (scale-independent).

**Core criterion:**

> Two lattices are algebraically compatible for exact commensurability if and
> only if δ_MOF / δ_sub is a rational perfect square.

$$\boxed{ \frac{\delta_o}{\delta_s} \in \mathbb{Q}^2 \;\Longleftrightarrow\; \text{exact commensurability algebraically feasible} }$$

Verification:

- Hexagonal / hexagonal: (−3)/(−3) = 1 → **feasible** (basis of magic-angle graphene)
- Square / square: (−4)/(−4) = 1 → **feasible**
- Square / hexagonal: (−4)/(−3) = 4/3 → **algebraically forbidden**

### 3.5  Level 0 Output

- **Forbidden**: report symmetry incompatibility; optionally continue to
  Level 1/2 to quantify the *approximate* match upper bound.
- **Feasible**: pass A_s, A_o and the transfer matrix T to Level 1.

---

## 4  Level 1: Reciprocal-Space Coincidence Function Φ(θ)

### 4.1  Physical Motivation

RHEED and LEED experiments observe interface matching directly in reciprocal
space.  The optimal epitaxial orientation angle is the one that maximises
overlap between the two diffraction patterns.  Φ(θ) formalises this
experimental observable and converts it into an optimisation target, avoiding
any enumeration of rotation angles.

### 4.2  Definition

Define the weighted reciprocal-lattice density (discrete Bragg peaks modelled
as a Dirac comb):

$$\rho^*(\mathbf{G}) = \sum_{h,k} w(|\mathbf{G}_{hk}|)\,\delta^{(2)}(\mathbf{G} - h\mathbf{b}_1 - k\mathbf{b}_2)$$

The weight w(|G|) = exp(−|G|²/G_c²) implements a **Debye–Waller-like
suppression**: high-frequency reflections are naturally damped by thermal
vibrations and interface roughness, so low-frequency Bragg peaks dominate the
epitaxial behaviour.

The coincidence function at rotation angle θ is their weighted convolution:

$$\Phi(\theta) = \sum_{\mathbf{G}_s}\sum_{\mathbf{G}_m}
  w_s\,w_m \exp\!\left(-\frac{|\mathbf{G}_s - R_\theta \mathbf{G}_m|^2}{2\sigma^2}\right)$$

Φ(θ) is proportional to the number of Bragg peaks shared by the two
lattices at angle θ — directly connected to the periodic overlap of electron
densities at the interface, and hence to the interface binding energy.

### 4.3  Fourier Series Representation (Jacobi–Anger Expansion)

Expand the squared distance:

$$|\mathbf{G}_s - R_\theta \mathbf{G}_m|^2
= |\mathbf{G}_s|^2 + |\mathbf{G}_m|^2 - 2|\mathbf{G}_s||\mathbf{G}_m|\cos(\alpha_{sm} - \theta)$$

where α_{sm} = ∠(G_s, G_m) is the fixed inter-set angle.  Apply the
**Jacobi–Anger identity** e^{z cos φ} = Σ_n I_n(z) e^{inφ}  (I_n = modified
Bessel function of the first kind):

$$\boxed{
  \Phi(\theta) = \sum_{n=-\infty}^{\infty} c_n\,e^{in\theta},
  \qquad
  c_n = \sum_{\mathbf{G}_s,\mathbf{G}_m}
        A_{sm}\,I_n\!\left(\tfrac{|\mathbf{G}_s||\mathbf{G}_m|}{\sigma^2}\right)
        e^{-in\alpha_{sm}}
}$$

### 4.4  Physical Interpretation of Fourier Coefficients

|c_n| reflects the overlap tendency under the **n-fold rotational symmetry
mode**:

- Square lattice: dominant contributions at n = 4, 8 (C4 symmetry)
- Hexagonal lattice: dominant contributions at n = 6, 3 (C6/C3 symmetry)

If one lattice has no C6 character (e.g. a square lattice), |c_6| ≈ 0 and
Φ(θ) has no hexagonal resonance regardless of θ.  This is the Fourier-space
manifestation of the Level 0 symmetry prohibition.

### 4.5  Peak Finding

Maxima of Φ(θ) satisfy Φ'(θ) = Σ_n in c_n e^{inθ} = 0, a truncated
Fourier series equation solved by FFT evaluation on a fine θ grid (O(K log K),
K ≈ 2000 points) followed by standard peak detection.

### 4.6  Level 1 Output

A ranked list of candidate angles {θ_k*, Φ(θ_k*)} forwarded to Level 2.
The peak height Φ(θ_k*) serves as a prior quality score.

---

## 5  Level 2: Floating-Point LLL Reduction and Strain Tensor

### 5.1  Physical Motivation

Given θ* from Level 1, the transfer matrix T(θ*) = R_{θ*} A_o A_s^{-1} is
known.  The task is to find integer row-vector pairs (m, n) ∈ Z^{1×2} × Z^{1×2}
satisfying m · A_o ≈ n · A_s, or equivalently m · T ≈ n.  This is a
**simultaneous Diophantine approximation** problem — finding the lattice
point in a 4-D integer lattice closest to a given target.

### 5.2  4-D Embedding Lattice

Construct the 4-D real lattice basis (rows are basis vectors):

$$B = \begin{pmatrix}
  T_{11} & T_{12} & \lambda & 0 \\
  T_{21} & T_{22} & 0       & \lambda \\
  -1     & 0      & 0       & 0 \\
  0      & -1     & 0       & 0
\end{pmatrix}$$

Any integer linear combination v = (k_1, k_2, k_3, k_4) · B satisfies

$$|\mathbf{v}|^2 =
\underbrace{(m_1 T_{11} + m_2 T_{21} - n_1)^2 + (m_1 T_{12} + m_2 T_{22} - n_2)^2}_{\text{strain penalty}}
+ \underbrace{\lambda^2(m_1^2 + m_2^2)}_{\text{supercell-size penalty}}$$

The **shortest vector** minimises the trade-off between epitaxial strain and
supercell size — exactly the physical balance in epitaxial growth (large
supercells achieve smaller strain at the cost of more interface dislocations).

The parameter λ controls this trade-off: large λ biases toward small-index
solutions; small λ allows large-index solutions with lower strain.  Scanning
a range of λ values explores the Pareto front.

### 5.3  LLL Reduction

The Lenstra–Lenstra–Lovász (LLL) algorithm produces a basis in which the
first vector satisfies

$$|\mathbf{v}_1| \leq 2^{(d-1)/4}\,\lambda_1(L)$$

where λ_1(L) is the true shortest vector length.  For d = 4: 2^{3/4} ≈ 1.68,
i.e. the result is near-optimal.  Time complexity O(d⁴ log B) is polynomial
in both the lattice dimension and the coefficient bit-length, in sharp
contrast to any enumeration approach.

The MOF supercell matrix M is recovered from the reduced basis vectors via
v[2:4] = λ · m, so m = round(v[2:4] / λ).  Two linearly independent rows
give the 2×2 matrix M, and N = round(M · T).

### 5.4  Green–Lagrange Strain Tensor

Define the **deformation gradient** F that maps the substrate supercell
C_s = N · A_s onto the (rotated) MOF supercell C_o = M · A_{o,rot}:

$$F = C_o\,C_s^{-1}$$

The **Green–Lagrange strain tensor** is

$$\boxed{ \varepsilon = \tfrac{1}{2}(F^T F - I) }$$

The choice of Green–Lagrange over the engineering strain (F − I) is
deliberate: ε is **rotation-invariant** (depends only on the symmetric
stretch factor U in the polar decomposition F = RU).  Equivalent orientation
domains therefore yield the same ε, regardless of their in-plane rotation.

Physical interpretation of tensor components:

| Component | Physical meaning |
|-----------|-----------------|
| ε_{11} | Normal strain along the first supercell basis vector (in-plane lattice constant shift measurable by GIXRD) |
| ε_{22} | Normal strain along the second supercell basis vector |
| ε_{12} | Shear strain (proportional to interfacial misfit dislocation density; related to Burgers vector density) |

The scalar mismatch degree η = ‖ε‖_F (Frobenius norm) is used for ranking.
Empirically η < 0.02 corresponds to physically feasible epitaxial growth.

### 5.5  Connection to DSC Lattice Theory

The short vectors produced by LLL are the basis vectors of Bollmann's
**Displacement Shift Complete (DSC) lattice** — the lattice of smallest
translations that do not increase the interface energy.  Level 2 therefore
yields not only the supercell indices but also the Burgers vector candidates
for the interfacial misfit dislocation network.

---

## 6  Information Flow

```
Input: A_s, A_o, config (sigma, eta_tol, lambda_values, ...)
──────────────────────────────────────────────────────────
Level 0: compute delta_s, delta_o
         │
         ├── delta_o/delta_s not a rational perfect square
         │   -> report symmetry prohibition
         │   -> (optional) continue to bound approximate match
         │
         └── pass -> transfer matrix T

Level 1: generate reciprocal-lattice sets {G_s}, {G_m}
         compute Fourier coefficients c_n       [O(N_G²)]
         evaluate Phi(theta) via FFT            [O(K log K)]
         -> ranked candidate angles {theta_k*, Phi(theta_k*)}

Level 2: for each theta_k*:
           build 4-D embedding lattice B(theta*, lambda)
           LLL reduction -> short vectors -> (M, N)
           compute strain tensor eps, eta = ||eps||_F
           if eta < eta_tol -> record match
──────────────────────────────────────────────────────────
Output: DataFrame sorted by eta; each row contains:
        theta (deg), M, N, eta, eps_11, eps_22, eps_12,
        area (A^2), L0_feasible flag
```

---

## 7  Comparison with the Original Algorithm

| Dimension | Original (phase2-2.py) | Three-level algorithm |
|-----------|------------------------|----------------------|
| Search strategy | Four nested loops, blind enumeration | Hierarchical: each level filters candidates |
| Complexity | O(N_m² · N_s²) ~ 10⁸ | O(1) + O(N_G²) + O(d⁴ log B) ≪ 10⁴ |
| False-positive rejection | Numerical tolerance (may miss) | Level 0 algebraic test (deterministic) |
| Rotation angle | Fixed; no search | Level 1 continuous optimisation |
| Duplicate supercells | Many redundant candidates | LLL produces a reduced basis; minimal redundancy |
| Strain metric | Decoupled percentages (length, angle) | Rotation-invariant Green–Lagrange tensor ε |
| Physical interpretability | Purely numerical filter | Every step has explicit physical/mathematical grounding |
| Tuning parameters | max_index, tol_length, tol_angle | sigma (diffraction resolution), eta_tol (strain tolerance), lambda (supercell-size preference) |

---

## 8  BFDH Face Pre-selection and Centering Extinction Rules

### 8.1  Physical Motivation

Before running the three-level matching pipeline, the crystal surface to
expose must be identified.  The **Bravais–Friedel–Donnay–Harker (BFDH) rule**
provides a rapid, structure-based prediction: faces with larger interplanar
spacings d_hkl grow more slowly (lower attachment energy) and therefore
dominate the equilibrium crystal morphology.  These are also the surfaces
most likely to be exposed in thin-film epitaxy experiments.

### 8.2  Interplanar Spacing

For a 3D lattice with reciprocal basis vectors **b**₁, **b**₂, **b**₃, the
spacing for the (*hkl*) family of planes is

$$d_{hkl} = \frac{2\pi}{|\mathbf{G}_{hkl}|}, \qquad
  \mathbf{G}_{hkl} = h\,\mathbf{b}_1 + k\,\mathbf{b}_2 + l\,\mathbf{b}_3$$

Faces are enumerated up to a maximum index |h|, |k|, |l| ≤ N_max, ranked in
decreasing order of d_hkl, and the top-N are forwarded to surface lattice
extraction and matching.

### 8.3  Centering Extinction Rules

For lattices with a non-primitive Bravais centering (I, F, A, B, C, R) the
unit cell contains additional lattice translations beyond the corner
positions.  By Bragg's law, reflections from planes related by these
translations interfere destructively whenever the centering condition is
violated, giving zero structure factor.  Such faces carry no diffraction
intensity and therefore provide no driving force for epitaxial nucleation.

The general conditions for the reflection (*hkl*) to be **present** are:

| Centering | Condition |
|-----------|-----------|
| P | all *hkl* (primitive — no centering extinction) |
| I (body-centred) | *h* + *k* + *l* = 2*n* |
| F (face-centred) | *h*, *k*, *l* all odd **or** all even |
| A | *k* + *l* = 2*n* |
| B | *h* + *l* = 2*n* |
| C | *h* + *k* = 2*n* |
| R (obverse) | −*h* + *k* + *l* = 3*n* |

**Implementation.**  The centering type is determined from the CIF via the
Hermann–Mauguin space group symbol (`_symmetry_space_group_name_H-M`, first
letter); or from `_symmetry_Int_Tables_number` looked up against a precompiled
table of all 230 space groups.  If neither tag is present, the structure is
treated as primitive (P), giving the conservative over-inclusive result.

**Screw-axis and glide-plane extinctions** (which affect zone-specific
reflection rows, e.g. *0kl* for a *b*-glide) are intentionally excluded:
these zone-specific conditions add significant complexity but improve BFDH
morphology predictions only marginally for typical MOF thin-film experiments.

### 8.4  Case Study: HKUST-1 (Fm̄3m, No. 225)

HKUST-1 (Cu₃(BTC)₂) crystallises in space group **Fm̄3m** — a face-centred
cubic (F) lattice with a ≈ 26.34 Å.

**Without extinction**, the three faces with the largest d_hkl are {100},
{010}, {001} (d = 26.34 Å each, square surface lattice).  This is physically
incorrect: for F-centering the condition is all-odd or all-even:

$$\text{(1,0,0):}\quad h=1\,(\text{odd}),\;k=0\,(\text{even}),\;l=0\,(\text{even})
\;\Rightarrow\; \text{mixed parity} \;\Rightarrow\; \textbf{absent}$$

**After the F-centering filter**, the true top BFDH faces are:

| *hkl* | *d* (Å) | Surface lattice | Allowed because |
|-------|---------|-----------------|----------------|
| (1,1,1) | 15.21 | *a* = *b* = 37.26 Å, γ = 60° | all odd |
| (2,0,0) | 13.17 | *a* = *b* = 26.34 Å, γ = 90° | all even |
| (2,2,0) |  9.31 | *a* = 37.26 Å, *b* = 26.34 Å, γ = 90° | all even |

This is in agreement with experimental observations: HKUST-1 thin films grown
by liquid-phase epitaxy on Au(111) expose the (111) face, presenting a
hexagonal copper-paddle-wheel motif that matches the hexagonal Au(111)
surface periodicity.  The corrected BFDH ranking therefore points the
epitaxy matcher in the physically correct direction from the outset.

---

## 9  Reciprocal-Space Visualization and LEED Spot Prediction

### 9.1  Physical Background

**Low-Energy Electron Diffraction (LEED)** probes the 2D periodicity of a
surface by measuring the momenta of backscattered electrons.  An electron
beam at normal incidence scatters elastically whenever the parallel momentum
transfer equals a 2D reciprocal lattice vector:

$$\Delta\mathbf{k}_{\parallel} = \mathbf{G}_{hk} = h\,\mathbf{b}_1 + k\,\mathbf{b}_2$$

where $\mathbf{b}_1, \mathbf{b}_2$ are the 2D reciprocal basis vectors of
the surface unit cell.  In a LEED experiment one therefore observes a **spot
pattern** that is the 2D reciprocal lattice of the topmost periodic layer.

When a MOF overlayer grows on a substrate, both lattices are present
simultaneously: the substrate contributes its own reciprocal lattice set
$\{\mathbf{G}_s\}$ and the oriented MOF overlayer contributes a rotated set
$\{\mathbf{G}_m\}$.  The resulting LEED screen shows spots from both.

### 9.2  Transformation of Reciprocal Vectors Under Rotation

The 2D reciprocal lattice matrix of a lattice with Cartesian matrix $A$ is

$$B = 2\pi\,(A^{-1})^T \quad\Longleftrightarrow\quad A\,B^T = 2\pi\,I$$

When the overlayer is rotated by angle $\theta$ (real-space matrix
$A' = R_\theta\,A$), the rotated reciprocal matrix becomes

$$B' = 2\pi\,(A'^{-1})^T
     = 2\pi\,\bigl[(R_\theta A)^{-1}\bigr]^T
     = 2\pi\,(A^{-1} R_\theta^T)^T
     = 2\pi\,R_\theta\,(A^{-1})^T
     = R_\theta\,B$$

**Reciprocal lattice vectors therefore transform under the same rotation as
the real-space vectors.**  The overlayer Bragg peaks visible in LEED are at

$$\mathbf{G}_m(h,k;\theta) = R_\theta \cdot (h\,\mathbf{b}_1^{(m)} + k\,\mathbf{b}_2^{(m)})$$

### 9.3  Superstructure Spots

In the commensurate case $M A_m = N A_s$, every overlayer reciprocal vector
can be written as a **rational linear combination** of substrate reciprocal
vectors.  Reciprocal-lattice points that belong to **both** sets — i.e.
$\mathbf{G}_s^{(i)} = \mathbf{G}_m^{(j)}$ for some integer indices $(i,j)$
— are scattered coherently by the interface and appear as **extra bright spots**
not indexable by the substrate lattice alone.  These are the **superstructure
reflections** ubiquitous in surface science literature.

Their positions can be predicted directly: every point in the set

$$\{\mathbf{G}_s : \exists\,j,\; |\mathbf{G}_s - \mathbf{G}_m^{(j)}| < \delta\}$$

is a superstructure spot, where $\delta \approx 0.18\,|\mathbf{b}_1^{(s)}|$
is a near-coincidence tolerance that accommodates the finite strain $\eta$
accepted by the matcher and any residual angular discretisation.

For a **coincident** (single-row commensurate) system the superstructure
spots form along a single family of lines in reciprocal space.  For a
**commensurate** (two-row) system they form a complete 2D sublattice — a
rectangular, oblique, or hexagonal net depending on $M$ and $N$.

### 9.4  Spot Count and Supercell Size

The number of non-equivalent superstructure spots within a reciprocal-space
disk of radius $G_{\max}$ scales as

$$N_{\text{super}} \approx \frac{\pi G_{\max}^2}{\Omega^*_{\text{super}}}$$

where $\Omega^*_{\text{super}}$ is the area of the superstructure reciprocal
unit cell.  For a supercell matrix $M$ with $|\det M| = n_{\text{cells}}$,
the superstructure reciprocal cell is $n_{\text{cells}}$ times smaller than
the substrate reciprocal cell:

$$\Omega^*_{\text{super}} = \frac{\Omega^*_s}{n_{\text{cells}}} = \frac{(2\pi)^2}{\Omega_s \cdot n_{\text{cells}}}$$

Large supercells ($n_{\text{cells}} \gg 1$) therefore produce **many dense
superstructure spots** — a signature of large-periodic moiré-type epitaxy,
visible as a closely spaced grid around every primary LEED spot.

### 9.5  (00) Beam

The (00) Bragg condition $\mathbf{G} = \mathbf{0}$ is always satisfied and
corresponds to the specularly reflected beam.  It appears at the centre of
every LEED pattern regardless of surface structure and is plotted as a filled
black circle.

### 9.6  Match-Card LEED Panel

The **bottom-left panel** of the match card produced by `plot_match_card()`
is a direct rendering of the reciprocal space described above:

| Symbol | Physical meaning |
|--------|-----------------|
| Blue circles | Substrate Bragg peaks ($\mathbf{G}_s$) |
| Red triangles | Rotated overlayer Bragg peaks ($\mathbf{G}_m$ at angle θ) |
| Green stars | Superstructure spots (coincident pairs within tolerance δ) |
| Black circle | (00) specular beam |

**What to tell your experimentalist:**
Count the green stars and measure their positions relative to the nearest
blue circle.  In experiment, superstructure spots are visible as extra LEED
spots with intensity proportional to the overlayer surface coverage and the
degree of long-range order.  Their specific positions pin down the supercell
matrix $M$ more precisely than any other measurement.

---

## 10  References

- **Zur & McGill (1984)**, J. Appl. Phys. 55, 378 — classical supercell
  enumeration framework (the approach replaced by this algorithm at Level 2)
- **Bollmann (1970)**, *Crystal Defects and Crystalline Interfaces*,
  Springer — O-lattice and DSC lattice theory; physical basis for Level 2
  dislocation interpretation
- **Lenstra, Lenstra & Lovász (1982)**, Math. Ann. 261, 515 — original LLL
  paper
- **Cao et al. (2018)**, Nature 556, 80 — magic-angle twisted bilayer
  graphene; motivates Level 0 quadratic-form classification
- **Ewald (1921)** — Ewald summation; theoretical basis for accelerating
  the Level 1 reciprocal-space sum
- **Jacobi–Anger expansion** — standard identity used in Level 1 to convert
  the coincidence function into a Fourier series
