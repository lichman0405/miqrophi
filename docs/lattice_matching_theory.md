# 2D Epitaxial Lattice Matching: Mathematical and Physicochemical Foundations

> This document derives the complete theoretical basis for the three-level
> hierarchical screening algorithm implemented in `miqrocal`.

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

## 8  References

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
