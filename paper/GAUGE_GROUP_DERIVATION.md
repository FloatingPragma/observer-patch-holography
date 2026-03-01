# Gauge Group Derivation: Why the Standard Model Is Inevitable

> OPH is the fundamental theory that exactly describes how our universe works, why it has the structure it has, and why it exists.
> The Standard Model, quantum field theory, general relativity, and string theory are effective descriptions of underlying OPH dynamics.
> From two input constants and five axioms (A1-A4 + MAR), OPH determines universe-wide properties, resolves incompatibilities, and explains measurement divergences including dark matter.


> **Companion document to**: *Observer-Patch Holography* (PAPER.md)
>
> **Abstract.** We prove that the extended OPH theory — core axioms A1–A4, regulator premises R0–R1, the loop-coherent gluing condition $[z]=0$, and the Selection Axiom MAR (Minimal Admissible Realization) — uniquely selects the Standard Model gauge group, the number of colors, and the number of generations. The result is not merely the SM Lie algebra but the exact global gauge group:
>
> $$G_{\mathrm{phys}} = \frac{SU(3) \times SU(2) \times U(1)}{\mathbb{Z}_6}, \qquad N_c = 3, \qquad N_g = 3.$$
>
> No other compact gauge group, no other color multiplicity, and no other generation count satisfies all admissibility conditions while minimizing the complexity vector $C(\Sigma)$.

---

## Table of Contents

1. [Premise Summary](#1-premise-summary)
2. [Admissibility Conditions](#2-admissibility-conditions)
3. [Selection Axiom MAR](#3-selection-axiom-mar)
4. [Main Theorem](#4-main-theorem)
5. [Proof: Gauge Group](#5-proof-gauge-group)
6. [Proof: Number of Colors](#6-proof-number-of-colors)
7. [Proof: Number of Generations](#7-proof-number-of-generations)
8. [Corollaries](#8-corollaries)
9. [Complete Derivation Chain](#9-complete-derivation-chain)
10. [Future Work](#10-future-work)

---

## 1. Premise Summary

### 1.1 The Extended Theory

The derivation proceeds within the extended OPH theory:

$$T_{\mathrm{ext}} := \underbrace{A1\text{–}A4}_{\text{core axioms}} + \underbrace{R0 + R1}_{\text{regulator premises}} + \underbrace{[z]=0}_{\text{loop-coherent gluing}} + \underbrace{\text{MAR}}_{\text{selection axiom}}$$

The roles of each component:

| Component | Content | Reference |
|-----------|---------|-----------|
| **A1** (Screen Net) | Horizon screen $S^2$ with algebra net $P \mapsto \mathcal{A}(P)$ | PAPER.md §1.3 |
| **A2** (Overlap Consistency) | Compatible marginals on overlaps | PAPER.md §1.3 |
| **A3** (Generalized Entropy) | Quantum focusing, $S_{\mathrm{gen}}$ monotonicity | PAPER.md §1.3 |
| **A4** (MaxEnt) | State selected by constrained entropy maximization | PAPER.md §1.3 |
| **R0** (Finite-dim regulator) | Local Hilbert spaces are finite-dimensional | PAPER.md §2.5 |
| **R1** (Boundary gauge structure) | Region algebras are gauge-invariant subalgebras | PAPER.md §2.5 |
| **$[z]=0$** | Loop-coherent gluing / DHR transportability | PAPER.md §6.6, Prop. 6.1a |
| **MAR** | Selection Axiom: Minimal Admissible Realization | This document, §3 |

### 1.2 What the Core Already Proves

From $A1$–$A4 + R0 + R1$ alone, the following is established (PAPER.md §6.1):

**Theorem 6.1 (Tannaka/DR reconstruction).** Edge-center completion yields a rigid symmetric $C^*$ tensor category $\mathsf{Sect}$ of edge sectors. There exists a compact group $G$, unique up to isomorphism, such that

$$\mathsf{Sect} \simeq \mathrm{Rep}(G), \qquad G = \mathrm{Aut}_\otimes(\mathcal{F}).$$

This proves *existence* of a compact gauge group but does not identify *which* compact group is realized. The present document determines it uniquely.

### 1.3 The Transportability Premise

**Proposition 6.1a.** DHR transportability — the condition that charges can be moved between patches without changing fusion rules — is equivalent to the vanishing of the central obstruction class:

$$\text{DHR transportable} \iff [z] = 0 \iff \text{loop-coherent gluing}.$$

This is an internal constraint on the sector structure, not an external physical assumption. It is kept explicitly visible as a separate premise ($[z]=0$) rather than hidden inside MAR.

---

## 2. Admissibility Conditions

An **OPH-realizable gauge sector** $\Sigma = (G, \text{matter content})$ consists of a compact gauge group $G$ reconstructed via Tannaka-Krein from the edge-sector category (Theorem 6.1) together with its associated chiral matter spectrum. A sector $\Sigma$ is **admissible** if it satisfies all of the following conditions:

### (i) Loop-coherent / transportable

The central obstruction class vanishes: $[z] = 0$. Equivalently, edge charges are DHR-transportable, i.e., they can be moved between patches without affecting fusion rules. This is the content of the explicit premise $[z] = 0$ (Proposition 6.1a).

*Rationale:* Without transportability, the reconstructed gauge symmetry is only a local labeling, not a global symmetry of the theory.

### (ii) Anomaly-free

The gauge sector is free of all perturbative anomalies (ABJ anomaly cancellation) and all global anomalies (Witten SU(2) anomaly, Dai-Freed anomalies). Concretely:

- $\mathrm{tr}[T_a \{T_b, T_c\}] = 0$ for all gauge generators (perturbative),
- $\pi_4(G)$ anomaly vanishes (Witten's global SU(2) anomaly),
- Mixed gravitational-gauge anomalies cancel.

*Rationale:* Anomalous gauge theories are inconsistent at the quantum level.

### (iii) Refinement-stable with light chiral matter

The MaxEnt/refinement-stable state supports light charged fermions without fine-tuning. This requires the matter content to be *chiral* — i.e., left-handed and right-handed fermions carry different gauge representations — so that vector-like mass terms are forbidden by gauge symmetry.

*Rationale:* Lemma 6.7 (PAPER.md §6.3) shows that refinement stability forbids unprotected relevant operators. A gauge-invariant Dirac mass term is relevant; if both chiralities are in conjugate representations, the mass term is allowed and drives the fermions to the cutoff. Corollary 6.8 then selects chiral content as the natural refinement-stable option.

### (iv) Single-Higgs Yukawa-completable

The chiral matter content can acquire mass through Yukawa couplings to a single scalar doublet $H = (1, 2, 1/2)$ after electroweak symmetry breaking, without requiring additional scalar multiplets.

*Rationale:* This is the minimal mechanism for generating fermion masses consistent with electroweak symmetry breaking. Additional Higgs multiplets would increase the scalar capacity without physical necessity.

### (v) Intrinsically CP-capable

The sector supports intrinsic CP violation — i.e., CP-violating phases that cannot be removed by field redefinitions. For a CKM-type mixing matrix with $N_g$ generations, the number of physical CP phases is $(N_g - 1)(N_g - 2)/2$.

*Rationale:* Observed CP violation in the kaon and B-meson systems requires intrinsic CP-violating phases. A sector that cannot accommodate them is empirically excluded.

### (vi) Weak-sector UV-completable

The weak gauge sector (the factor carrying the pseudoreal doublet) is asymptotically free at one loop, ensuring that it is UV-completable within the OPH framework rather than requiring a Landau pole cutoff below the screen scale.

*Rationale:* The OPH screen operates at the Planck scale. A gauge factor with a sub-Planckian Landau pole would be inconsistent with the screen construction.

---

## 3. Selection Axiom MAR

**Selection Axiom MAR (Minimal Admissible Realization).** Among all OPH-realizable sectors $\Sigma$ that satisfy admissibility conditions (i)–(vi), Nature realizes the **lexicographically minimal** one under the complexity vector

$$C(\Sigma) = \bigl(\chi_{\mathrm{faith}},\; N_{\mathrm{nonab}},\; N_c,\; N_g\bigr),$$

where:
- $\chi_{\mathrm{faith}}$ is the **faithful edge capacity**: the dimension of the minimal faithful unitary representation of $G$ (the smallest representation on which all gauge factors act nontrivially),
- $N_{\mathrm{nonab}}$ is the **number of nonabelian simple factors** in $G$,
- $N_c$ is the **rank of the color factor** (the dimension of the fundamental representation of the complex nonabelian factor),
- $N_g$ is the **number of chiral generations**.

Lexicographic minimality means: first minimize $\chi_{\mathrm{faith}}$; among ties, minimize $N_{\mathrm{nonab}}$; among further ties, minimize $N_c$; among final ties, minimize $N_g$.

**Key properties of MAR:**

1. *It is a selection rule, not a dynamics.* MAR acts on an independently defined admissible class; it does not modify the OPH equations of motion.

2. *It is not "minimalism in general."* MAR only selects among sectors that already pass all six admissibility filters. Plain minimality (smallest $G$) would yield $U(1)$ or $SU(2)$, which fail conditions (iii)–(v). The admissibility conditions do the heavy lifting; MAR breaks the remaining degeneracy.

3. *It is powerful.* Product gauge structure, the specific factors $SU(3) \times SU(2) \times U(1)$, and the values $N_c = 3$, $N_g = 3$ all follow from MAR applied to the admissible class.

---

## 4. Main Theorem

**Theorem (Standard Model Inevitability).** Under the extended OPH theory $T_{\mathrm{ext}} = A1\text{–}A4 + R0 + R1 + [z]=0 + \mathrm{MAR}$:

$$G_{\mathrm{phys}} = \frac{SU(3) \times SU(2) \times U(1)}{\mathbb{Z}_6}, \qquad N_c = 3, \qquad N_g = 3.$$

*The proof is given in §5–§7 below.*

---

## 5. Proof: Gauge Group

The proof proceeds in five steps. Steps 1–4 determine the connected gauge group. Step 5 fixes the global structure.

### Step 1: Some compact group exists

**From:** $A1$–$A4 + R0 + R1$.

Edge-center completion (Theorem 2.3, PAPER.md) yields a rigid symmetric $C^*$ tensor category $\mathsf{Sect}$ of edge charges. By Theorem 6.1 (Tannaka/DR reconstruction), there exists a compact group $G$, unique up to isomorphism, with $\mathsf{Sect} \simeq \mathrm{Rep}(G)$.

The premise $[z] = 0$ ensures that this reconstruction is global: charges are DHR-transportable, so the compact group acts as a genuine gauge symmetry, not merely a local labeling (Proposition 6.1a). $\square$

### Step 2: The gauge group must be a product

**From:** MAR (minimizing $\chi_{\mathrm{faith}}$) + admissibility (iii).

**Lemma 5.1.** Any admissible sector must contain at least two genuinely different nonabelian charge types: one pseudoreal and one complex.

*Proof.* Admissibility condition (iii) requires light chiral matter. By Corollary 6.8 (PAPER.md §6.3), refinement stability forces the matter content to be chiral. To support chiral fermions that can form gauge-invariant Yukawa couplings with a single Higgs doublet (condition iv), the gauge group must admit both:
- A **pseudoreal** fundamental representation (for the weak doublet structure — this allows left-handed doublets and right-handed singlets to have different gauge quantum numbers), and
- A **complex** fundamental representation (for the color structure — this allows quarks to carry a charge that distinguishes particle from antiparticle).

A single nonabelian factor cannot simultaneously provide both types: a group is either pseudoreal at a given dimension or complex, not both. $\square$

**Lemma 5.2.** The minimal faithful carrier is $V = \mathbb{C}^3 \otimes \mathbb{C}^2$, with $\chi_{\mathrm{faith}} = 6$.

*Proof.*
- The smallest faithful pseudoreal irreducible representation of any compact simple group is the fundamental **2** of $SU(2)$ (the only 2D pseudoreal irrep).
- The smallest faithful irreducible complex representation of any compact simple group is the fundamental **3** of $SU(3)$ (the 2D fundamental of $SU(2)$ is pseudoreal, not complex; the fundamental of $SU(2)$ in dimension 2 is the only option for pseudoreal, while $SU(3)$ in dimension 3 is the first genuinely complex case).
- For both charge types to act faithfully, the carrier must accommodate both simultaneously. The minimal such space is the tensor product $V = \mathbb{C}^3 \otimes \mathbb{C}^2$, giving $\chi_{\mathrm{faith}} = 6$.
- Any other combination (e.g., using $SU(4)$ instead of $SU(3)$, or $Sp(4)$ instead of $SU(2)$) gives $\chi_{\mathrm{faith}} > 6$, and is therefore disfavored by MAR. $\square$

**Proposition 5.3 (Product structure from minimal carrier).** The minimal faithful carrier $\mathbb{C}^3 \otimes \mathbb{C}^2$ forces a product gauge structure.

*Proof.* On $V = \mathbb{C}^3 \otimes \mathbb{C}^2$, the color factor acts on the first tensor factor and the weak factor acts on the second. These actions commute by the tensor product structure. Therefore the gauge group decomposes as a product $G_{\mathrm{color}} \times G_{\mathrm{weak}} \times G_{\mathrm{abelian}}$ (up to finite quotient). A simple group like $SU(5)$ or $SO(10)$ acting irreducibly on a 6-dimensional space would require $\chi_{\mathrm{faith}} \geq 5$ or $\chi_{\mathrm{faith}} \geq 10$ respectively, and crucially would not provide the independent pseudoreal + complex structure required by Lemma 5.1.

Product structure is derived from the minimal-carrier argument. $\square$

### Step 3: The factors are SU(3), SU(2), and U(1)

**From:** Lemmas 6.3–6.5 of PAPER.md, applied to the minimal carrier.

**Lemma 5.4 (SU(2) from pseudoreal doublet).** The factor acting on $\mathbb{C}^2$ as a faithful 2D pseudoreal representation is $SU(2)$.

*Proof.* Among compact groups with a faithful 2D pseudoreal unitary representation, $SU(2)$ is the unique option. (The only compact Lie groups with a 2D faithful representation are subgroups of $U(2)$. Among these, the pseudoreal condition $V \cong \bar{V}$ via an antisymmetric bilinear form selects exactly $SU(2)$.) This is Lemma 6.3 of PAPER.md. $\square$

**Lemma 5.5 (SU(3) from complex triplet).** The factor acting on $\mathbb{C}^3$ as a faithful irreducible complex representation is $SU(3)$.

*Proof.* Among compact groups with a faithful irreducible complex 3D representation, the connected component is $SU(3)$. (The only compact simple Lie groups with a 3D fundamental representation are $SU(3)$ and $SO(3)$. But the fundamental of $SO(3)$ is real, not complex. Therefore $SU(3)$ is uniquely selected.) This is Lemma 6.4 of PAPER.md. $\square$

**Lemma 5.6 (U(1) from continuous characters).** Admissibility condition (iv) requires continuously parameterized hypercharges for the Yukawa couplings. A continuous family of 1D sectors yields a $U(1)$ factor.

*Proof.* This is Lemma 6.5 of PAPER.md. The hypercharge assignments must form a continuous family (not a discrete group) to allow the full range of Yukawa couplings needed for condition (iv). $\square$

### Step 4: Nothing else — the commutant argument

**Proposition 5.7 (No additional factors).** No extra gauge factor — neither nonabelian nor abelian — can appear without increasing $\chi_{\mathrm{faith}}$ beyond 6.

*Proof.* Consider the maximal compact subgroup of $U(6)$ acting on $V = \mathbb{C}^3 \otimes \mathbb{C}^2$ with commuting actions on each tensor factor:

$$\{g \in U(6) : g = g_3 \otimes g_2,\; g_3 \in U(3),\; g_2 \in U(2)\} \cong U(3) \times U(2) / U(1)_{\mathrm{diag}}.$$

The continuous symmetry group with commuting color and weak actions is therefore

$$S(U(3) \times U(2)) \cong \frac{SU(3) \times SU(2) \times U(1)}{\text{finite center}}.$$

Now compute the commutant:
- The commutant of $SU(3) \times SU(2)$ inside $U(6)$ is exactly $U(1)$.
- Therefore no extra continuous factor (neither $U(1)'$ nor any nonabelian group) can appear without enlarging the carrier beyond $\chi = 6$.
- Any additional nonabelian factor would require its own faithful representation, increasing $\chi_{\mathrm{faith}}$.

MAR selects the minimal $\chi_{\mathrm{faith}}$, so extra factors are excluded. $\square$

### Step 5: The global quotient is $\mathbb{Z}_6$

**Proposition 5.8 (Z₆ quotient from hypercharge quantization).** If the realized matter spectrum has hypercharges quantized in sixths, then the subgroup of $SU(3) \times SU(2) \times U(1)$ acting trivially on all physical states is exactly $\mathbb{Z}_6$.

*Proof.* The SM matter content with hypercharge assignments $(Q_L, u_R, d_R, L_L, e_R, H)$ has hypercharges $Y \in \{1/6, 2/3, -1/3, -1/2, -1, 1/2\}$, all in $\frac{1}{6}\mathbb{Z}$. The center of $SU(3) \times SU(2) \times U(1)$ is $\mathbb{Z}_3 \times \mathbb{Z}_2 \times U(1)$. The subgroup acting trivially on all realized representations is generated by

$$(\omega_3, -1, e^{i\pi/3}) \in SU(3) \times SU(2) \times U(1),$$

where $\omega_3 = e^{2\pi i/3}$ is the $\mathbb{Z}_3$ center of $SU(3)$. This generates a $\mathbb{Z}_6$ subgroup. Therefore

$$G_{\mathrm{phys}} = \frac{SU(3) \times SU(2) \times U(1)}{\mathbb{Z}_6}.$$

Equivalently, $G_{\mathrm{phys}} \cong S(U(3) \times U(2))$.

The hypercharge quantization itself follows from anomaly cancellation (condition ii) with the minimal Yukawa-complete spectrum (condition iv). Theorem 6.13 of PAPER.md derives the hypercharge assignments from these conditions for $N_c = 3$, yielding precisely the sixth-integer lattice. This is Proposition 6.6 of PAPER.md. $\square$

### Summary of Step 1–5

$$\boxed{A1\text{–}A4 + R0 + R1 + [z]=0 + \mathrm{MAR} \;\Longrightarrow\; G_{\mathrm{phys}} = \frac{SU(3) \times SU(2) \times U(1)}{\mathbb{Z}_6}}$$

---

## 6. Proof: Number of Colors

**Theorem (N_c = 3).** Under $T_{\mathrm{ext}}$, the rank of the color factor is $N_c = 3$.

**Proof.** The gauge group derivation (§5) already fixes the color factor as $SU(3)$ — i.e., $N_c = 3$ — through the minimal carrier argument (Lemma 5.2). We give here the independent confirmation from anomaly constraints + MAR.

**Step 1: Witten anomaly constrains N_c to be odd.**

With gauge structure $SU(N_c) \times SU(2)_L \times U(1)_Y$ and one left-handed quark doublet $Q$ per color plus one left-handed lepton doublet $L$ per generation, the total number of $SU(2)$ doublets per generation is

$$N_{\mathrm{doublets}} = N_c + 1.$$

Witten's global $SU(2)$ anomaly (1982) requires the total number of $SU(2)$ doublets to be even:

$$N_c + 1 \equiv 0 \pmod{2} \quad \Longrightarrow \quad N_c \text{ is odd}.$$

This alone allows $N_c \in \{1, 3, 5, 7, \ldots\}$.

**Step 2: N_c = 1 fails admissibility.**

For $N_c = 1$: $SU(1)$ is trivial (no color dynamics). The fundamental representation is 1-dimensional and real, not complex. This violates Lemma 5.1 (no genuinely complex nonabelian charge type) and condition (iii) (cannot support chiral quarks).

**Step 3: MAR selects N_c = 3.**

Among the remaining candidates $\{3, 5, 7, \ldots\}$, all with $\chi_{\mathrm{faith}} = 2N_c$ (since the minimal carrier is $\mathbb{C}^{N_c} \otimes \mathbb{C}^2$), MAR's complexity vector has $N_c$ in the third slot. Lexicographic minimization selects $N_c = 3$. $\square$

---

## 7. Proof: Number of Generations

**Theorem (N_g = 3).** Under $T_{\mathrm{ext}}$, the number of chiral generations is $N_g = 3$.

**Proof.**

**Step 1: CP violation gives N_g ≥ 3.**

Admissibility condition (v) requires intrinsic CP violation. The number of physical CP-violating phases in an $N_g \times N_g$ CKM matrix is

$$n_{\mathrm{CP}} = \frac{(N_g - 1)(N_g - 2)}{2}.$$

- For $N_g = 1$: $n_{\mathrm{CP}} = 0$. No CP violation possible. Excluded.
- For $N_g = 2$: $n_{\mathrm{CP}} = 0$. No CP violation possible. Excluded.
- For $N_g = 3$: $n_{\mathrm{CP}} = 1$. CP violation possible.

Therefore $N_g \geq 3$.

**Step 2: Asymptotic freedom gives N_g ≤ 5.**

Admissibility condition (vi) requires $SU(2)_L$ to be asymptotically free. The one-loop beta function coefficient is

$$b_{1,SU(2)} = \frac{1}{3}\bigl[22 - N_g(N_c + 1)\bigr].$$

Asymptotic freedom requires $b_{1,SU(2)} > 0$:

$$N_g(N_c + 1) < 22.$$

With $N_c = 3$ (from §6): $4N_g < 22$, so $N_g \leq 5$.

Combining: $3 \leq N_g \leq 5$.

**Step 3: MAR selects N_g = 3.**

Among the admissible candidates $\{3, 4, 5\}$, MAR's complexity vector has $N_g$ in the fourth (last) slot. Lexicographic minimization selects $N_g = 3$. $\square$

---

## 8. Corollaries

### Corollary 1: No gauge-mediated proton decay

**Corollary.** The gauge group $G_{\mathrm{phys}} = SU(3) \times SU(2) \times U(1)/\mathbb{Z}_6$ does not contain $X$ or $Y$ gauge bosons that mediate proton decay.

*Proof.* Proton decay via gauge bosons requires a simple unification group (e.g., $SU(5)$, $SO(10)$) whose additional generators connect quarks to leptons. The MAR derivation forces a product group structure with $\chi_{\mathrm{faith}} = 6$, which excludes simple groups. Therefore no gauge bosons mediating $B$-violating transitions exist. $\square$

### Corollary 2: Uniqueness

**Corollary.** The SM gauge group is the *unique* solution. No other compact gauge group satisfies all admissibility conditions while achieving $\chi_{\mathrm{faith}} = 6$.

*Proof.* By the commutant argument (Proposition 5.7), the group $SU(3) \times SU(2) \times U(1)/\mathbb{Z}_6$ exhausts the continuous symmetry of the minimal carrier $\mathbb{C}^3 \otimes \mathbb{C}^2$. Any alternative group either:
- Has $\chi_{\mathrm{faith}} > 6$ (excluded by MAR), or
- Fails to provide both pseudoreal and complex representations (excluded by condition iii), or
- Is a subgroup of the SM group, which would fail anomaly cancellation or Yukawa completeness.

$\square$

### Corollary 3: Hypercharge quantization

**Corollary.** The hypercharge assignments are uniquely fixed to the observed values $Y \in \frac{1}{6}\mathbb{Z}$ by anomaly cancellation within the derived gauge group.

*Proof.* Theorem 6.13 of PAPER.md derives the unique anomaly-free hypercharge assignment for $SU(3) \times SU(2) \times U(1)$ with $N_c = 3$ and $N_g = 3$, assuming single-Higgs Yukawa completeness (condition iv). The result is the standard sixth-integer lattice. $\square$

---

## 9. Complete Derivation Chain

The full logical chain from axioms to the Standard Model gauge sector is:

> $A1$–$A4 + R0 + R1$ $\xrightarrow{\text{Thm 6.1}}$ $\exists\, G$ compact
>
> $+ \; [z]=0$ $\xrightarrow{\text{Prop 6.1a}}$ DHR transportable, global gauge symmetry
>
> $+ \;\text{MAR}$ $\xrightarrow{\text{Thm (this doc)}}$ $G_{\mathrm{phys}} = SU(3) \times SU(2) \times U(1)/\mathbb{Z}_6$, $N_c = 3$, $N_g = 3$
>
> $\xrightarrow{\text{Thm 6.13}}$ hypercharges fixed $\xrightarrow{\text{6.11}}$ photon and graviton inevitable
>
> $\xrightarrow{\text{Spectrum Derivation}}$ particle masses from pixel area

Under $T_{\mathrm{ext}}$, the following are all derived results:

1. Product gauge group structure
2. Exact global gauge group $SU(3) \times SU(2) \times U(1)/\mathbb{Z}_6$
3. Number of colors $N_c = 3$
4. Number of generations $N_g = 3$
5. Absence of gauge-mediated proton decay
6. Hypercharge quantization in sixths

This is a complete, gap-free chain from the extended axiom set to the observed gauge structure of Nature. $\square$

---

## 10. Future Work

The gauge-group derivation presented here determines the symmetry structure of the Standard Model. Several related questions are addressed in companion documents:

- **Fermion mass hierarchy and Yukawa couplings** — see [TECHNICAL_SUPPLEMENT.md](TECHNICAL_SUPPLEMENT.md) §12
- **Cosmological constant** — see [TECHNICAL_SUPPLEMENT.md](TECHNICAL_SUPPLEMENT.md) §14
- **Full transmutation chain** — see [TECHNICAL_SUPPLEMENT.md](TECHNICAL_SUPPLEMENT.md) Appendix A.3
- **Heat-kernel parameter $t$ from UV microphysics** — see [PAPER.md](PAPER.md) §6.13
- **Particle mass spectrum from pixel geometry** — see [SPECTRUM_DERIVATION.md](SPECTRUM_DERIVATION.md)
