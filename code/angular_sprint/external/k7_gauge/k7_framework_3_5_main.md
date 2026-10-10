# The K₇ Framework: Conditional Topological Relations for Standard-Model Parameters from G₂ Holonomy and E₈×E₈ Structure

**Brieuc de La Fourniere**

*Independent researcher, Beaune, France*

---

## Abstract

The Standard Model requires 19 experimentally determined parameters lacking theoretical explanation. We present a geometric framework in which physical observables emerge as topological invariants of a seven-dimensional G₂ holonomy manifold K₇ coupled to E₈×E₈ gauge structure, with topological inputs treated as discrete and the metric determinant det(g) = 65/32 imposed as a normalization target (§2.3). The 33 Type I relations employ **no continuously adjustable parameter after the declared structural ledger is frozen**; Types II, III, and IV are *conditional extractions* using explicitly-enumerated inputs (§3.2). Building on a Newton–Kantorovich-certified neck/metric model ([A], [B]), a Joyce–Karigiannis route realizing (b₂, b₃) = (21, 77), and the companion analytic scheme [E] (discharged at a normalised datum $\mathcal D_0$), we derive **95 observables** organized in four types: 33 direct algebraic (Type I, mean deviation 0.73%), 19 one-step physical extractions (Type II, 0.17%), 21 multi-step dynamical chains (Type III, 3.4%), and 22 structural diagnostics (Type IV). Of 66 experimentally comparable observables, 11 are exact matches (deviation < 0.01%) and 53 are within 1%. The existence question is discharged conditionally at the *datum-level analytic layer* by [E]; at the *global-analytic layer* it remains conditional on $\mathsf H_{\rm global}$ and hypothesis $(\mathrm J)$ of [E, §8.3] (two-layer boundary: §9).

Beyond the original 33 predictions, the framework adds: (1) a *representation-branching* chain E₈ ⊃ E₆ × SU(3) ⊃ ⋯ ⊃ SU(3) × SU(2) × U(1) compatible with the SM matter content (a chain of subgroup branchings, not a physical breaking mechanism, whose UV completion remains a load-bearing open question; §5.1, §9); (2) a combined lepton mass hierarchy mechanism achieving sub-percent precision from two independent geometric sources; and (3) a **Sieve reading** under which m_H/m_W (Route A) is the unique survivor at exact rank 1 budget-unique after Type I calibration (§7). Of the 95 observables, 55 are formally verified in Lean 4 (213 certificate conjuncts, **15 classified axioms across the A-F taxonomy**, 0 sorry); Lean certifies the algebraic relations and internal consistency, not the physical interpretation or the existence of the complete compact construction (§8.3, §9).

DUNE will test the topological prediction δ_CP = 197°; a measurement outside [182, 212]° would create serious tension. NuFIT 6.1 (2025) reports δ_CP = $207^{+23}_{-20}$°, so the prediction now sits well within the 1σ band. We present this as an exploratory investigation emphasizing falsifiability, not a claim of correctness.

**Keywords**: G₂ holonomy, exceptional Lie algebras, Standard Model parameters, topological field theory, falsifiability, formal verification, Lean 4

---

**Framework epistemic status at a glance (four layers):**

*Formerly presented as the Geometric Information Field Theory (GIFT), versions 1.0-3.4 of this record.*

- **Formally established.** Conditional algebraic identities among topological invariants of $K_7$ (33 Type I relations), verified in Lean 4 (213 conjuncts, 15 classified axioms in the A–F taxonomy, 0 sorry). The Lean layer certifies internal consistency, not physical interpretation.
- **Established at the datum.** Newton–Kantorovich-certified seam metric [A] and JK17 topological/lattice route producing $(b_2, b_3) = (21, 77)$; datum-level analytic scheme [E] (concept DOI 10.5281/zenodo.21215956) at $\mathcal D_0$ with $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$ by outward-rounded interval arithmetic.
- **Open.** (i) Global compact torsion-free $G_2$ construction on a smooth $K_7$, conditional on $\mathsf H_{\rm global}$ + $(\mathrm J)$ of [E]. (ii) UV bridge to non-abelian chiral gauge sector (Acharya–Witten singularity condition [hep-th/0109152] or equivalent heterotic-bundle/duality realisation); see §5.1 caveat and §9.
- **Phenomenology.** 95 observables organised in four types with explicit dependency ranking: Type I algebraic (no continuously tuned parameter after ledger freeze); Types II, III, IV conditional extractions using experimental anchors / dynamical assumptions / structural diagnostics. Falsifiable predictions listed in §11.1.

---

## 1. Introduction

### 1.1 The Parameter Problem

The Standard Model describes fundamental interactions with remarkable precision, yet requires 19 free parameters determined solely through experiment [1]. These parameters (gauge couplings, Yukawa couplings spanning five orders of magnitude, mixing matrices, and Higgs sector values) lack theoretical explanation.

Several tensions motivate the search for deeper structure:

- **Hierarchy problem**: The Higgs mass requires fine-tuning absent new physics [2].
- **Hubble tension**: CMB and local H₀ measurements differ by >4 sigma [3,4].
- **Flavor puzzle**: No mechanism explains three generations or mass hierarchies [5].
- **Koide mystery**: The charged lepton relation Q = 2/3 holds for 43 years without explanation [6].

These challenges suggest examining whether parameters might emerge from geometric or topological structures.

### 1.2 Framework Overview

The K₇ framework (introduced as the Geometric Information Field Theory, GIFT, in versions 1.0-3.4 of this record) proposes that dimensionless parameters represent topological invariants of an eleven-dimensional spacetime:

```
E₈ x E₈ (496D) --> AdS₄ x K₇ (11D) --> Standard Model (4D)
```

The key elements:

1. **E₈ x E₈ gauge structure** (dimension 496)
2. **Compact 7-manifold K₇** with G₂ holonomy (b₂ = 21, b₃ = 77)
3. **Metric normalization constraint**: det(g) = 65/32, structurally motivated by G₂/E₈ constants (§2.3) but imposed as a normalization target, not derived from topology alone
4. **Cohomological mapping**: Betti numbers constrain field content

We emphasize this represents mathematical exploration, not a claim that nature realizes this structure. The framework's merit lies in falsifiable predictions from topological inputs.

### 1.3 Paper Organization

- **Section 2**: Mathematical framework (E₈×E₈, K₇, G₂ structure)
- **Section 3**: Methodology, epistemic status, and type classification
- **Section 4**: 95 observables across 4 types (33(I) + 19(II) + 21(III) + 22(IV))
- **Section 5**: Gauge sector: E₈ → SM breaking, anomaly cancellation, B-test identity, bundle universality
- **Section 6**: Mass hierarchy: Wilson lines, instantons, combined wilson_line+instanton pipeline
- **Section 7**: Statistical uniqueness under the **Sieve reading** (public 4-null battery + Lean 4 formal-identity flag); traditional sensitivity diagnostics
- **Section 8**: Formal verification (213 Lean conjuncts, 15 classified axioms in the A-F taxonomy, of which 4 external data-package axioms, 0 sorry) and statistical uniqueness
- **Section 9**: **Existence status (two-layer boundary)**: datum-level analytic layer discharged conditionally in the companion paper [E]; global-analytic layer conditional on the two-slot pack $\mathsf H_{\rm global}$ + hypothesis $(\mathrm J)$ of [E]
- **Section 10**: **Selection audit** (why $(21, 77)$?): five investigated routes closed; the selection question itself remains open, with the sole residual $b_2 + b_3 = 98 = \dim(K_7) \cdot \dim(G_2)$ flagged as a pre-registration target
- **Section 11**: Discussion, falsifiability, and conclusion

Four supplements provide technical details: S1 (Mathematical Foundations), S2 (Complete Derivations), S3 (Observable Dataset with full 95-entry table), S4 (Sieve diagnostics: archived coincidence-probability tables from §7.5 of the 3.4 version).

---

## 2. Mathematical Framework

### 2.1 The Octonionic Foundation

The K₇ framework emerges from the algebraic fact that **the octonions are the largest normed division algebra**.

| Algebra | Dim | Physics Role | Extends? |
|---------|-----|--------------|----------|
| R | 1 | Classical mechanics | Yes |
| C | 2 | Quantum mechanics | Yes |
| H | 4 | Spin, Lorentz group | Yes |
| **O** | **8** | **Exceptional structures** | **No** |

The octonions terminate this sequence. Their automorphism group G₂ = Aut(O) has dimension 14 and acts naturally on Im(O) = R^7. The exceptional Lie algebras arise from octonionic constructions through a chain established by Dray and Manogue [17]:

| Algebra | Dimension | Connection to O |
|---------|-----------|-----------------|
| G₂ | 14 | Aut(O) |
| F₄ | 52 | Aut(J₃(O)) |
| E₆ | 78 | Collineations of OP² |
| E₈ | 248 | Contains all lower exceptionals |

This chain is not accidental. It reflects the unique algebraic structure of the octonions: Im(O) has dimension 7, the Fano plane encodes the multiplication table, and G₂ preserves this structure. Octonionic approaches to Standard-Model structure form an independent line of research [7, 8, 10, 37]. A G₂-holonomy manifold is therefore the natural geometric home for octonionic physics, just as U(1) holonomy is the natural setting for complex geometry.

The G₂ structure is concretely encoded in the **standard 3-form φ₀** on ℝ⁷ (Bryant-Joyce convention):
$$\varphi_0 = e_{012} + e_{034} + e_{056} + e_{135} - e_{146} - e_{236} - e_{245}$$
where only 7 of the $\binom{7}{3} = 35$ ordered triples carry nonzero coefficient (all ±1). G₂ is precisely the stabilizer of φ₀ in GL(7, ℝ), and its Lie algebra g₂ is the kernel of the linear map $L_{\varphi_0} : \mathfrak{gl}(7) \to \wedge^3(\mathbb{R}^7)^*$, $X \mapsto \mathcal{L}_X \varphi_0$, giving $\dim(g_2) = 49 - \operatorname{rank}(L_{\varphi_0}) = 49 - 35 = 14$. This is fully formalized in Lean (`G2ThreeForm.lean`): all 7 nonzero coefficients are certified by `native_decide`, and the Lie algebra structure (closure under addition and scalar multiplication) is proven.

### 2.2 E₈ x E₈ Structure

E₈ is the largest exceptional simple Lie group with dimension 248 and rank 8 [18]. The product E₈ x E₈ arises in heterotic string theory for anomaly cancellation [19], with total dimension 496.

The first E₈ contains the Standard Model gauge group at the level of representation branchings:

```
E₈ ⊃ E₆ × SU(3) ⊃ SO(10) × U(1)_1 ⊃ SU(5) × U(1)_2 ⊃ SU(3)_c × SU(2)_L × U(1)_Y
```

We emphasise that this is a chain of **subgroup branchings compatible with the SM representations**, not a physical breaking mechanism: an actual heterotic/M-theory breaking would specify the mechanism (Wilson-line, bundle, singular locus) responsible for each step. The intermediate SU(3) factor after the first branching remains a hidden-sector algebra (a plausible dark-sector candidate; not committed to here). The intermediate step SU(5) × U(1)_2 is **not** Pati–Salam (which is $SU(4) \times SU(2)_L \times SU(2)_R$); an alternative chain via Pati–Salam is possible but is not the one used in this framework. The second E₈ provides a hidden sector whose physical interpretation remains an open question.

Wilson (2024) demonstrates that E₈(-248) encodes three fermion generations (128 degrees of freedom) with GUT structure [9]. The product dimension 496 enters the hierarchy parameter tau = (496 x 21)/(27 x 99) = 3472/891, connecting gauge structure to internal topology.

### 2.3 The K₇ Manifold

We consider a compact, simply connected 7-manifold K₇ with G₂ holonomy and Betti numbers (b₂, b₃) = (21, 77). G₂ holonomy preserves N = 1 SUSY in 4D [11], implies Ric(g) = 0, and is the unique holonomy of G₂ = Aut(𝕆) acting on ℝ⁷.

**Metric construction** ([A]): An explicit 7D G₂ metric g(s, θ, ψ, y) = g_seam(s) ⊕ g_{T²} ⊕ g_{K3}(y) on K₇ is constructed via Chebyshev-Cholesky parametrization (169 parameters), with K3 fiber contribution certified at 0.07% of total torsion [C]. The Newton-Kantorovich certificate (h = 8.95 × 10⁻⁹, margin ×56 million, zero finite differences) establishes existence and local uniqueness (within the NK certificate ball) of a nearby torsion-free metric with δg/g ≤ 4.86 × 10⁻⁶. The spectral analysis ([B]) independently confirms 21 near-zero eigenvalues of Δ₂ (gap ratio 14,635) and 77 near-zero eigenvalues of Δ₃, consistent with the Betti numbers (b₂, b₃) = (21, 77).

**Metric decomposition**: A Chebyshev mode analysis reveals that the metric is dominated by its constant mode (k=0), which carries 99.9998% of the L² energy. This mode is a product metric K3 × T² × I with structural parameters g_ss = 19/6, g_{T²} = 7/6, det(g) = 65/32 (fixed by the structural normalization and topological input data, with det(g) = 65/32 treated as a metric normalization target as discussed in S1 §10.3). The k≥1 modes (0.0002% of energy) constitute the minimal perturbation that breaks the product structure and lifts the holonomy from SU(2)×U(1) to full G₂; these corrections drive the torsion ‖T‖ = 2.949×10⁻⁵ to the NK certification threshold, but contribute nothing to the numerical values of the 95 observables, which depend only on topological integers. The Betti numbers of $K_7$, including $b_1 = 0$, are topological invariants of the construction, not features generated by the metric perturbation: they are fixed by the topological data of $K_7$ upstream of any metric choice.

**Topological classification**: Among the *catalogued* compact G₂ examples (Joyce's 252 orbifold types [20] and the ~100 CHNP TCS examples [21, 22]; see [12, 14, 15] for the surrounding construction literature and census) the pair (b₂, b₃) = (21, 77) does not appear, and orthogonal TCS is excluded by parity (b₂+b₃=98 is even; CHNP Lemma 6.7 [22]). However, a Joyce-Karigiannis (JK) Z₂³ orbifold construction [43] *does* produce (b₂, b₃) = (21, 77): a four-phase computer-assisted audit (V4 symplectic screen on CI(2,2,2), anti-symplectic obstruction, K3 lattice existence via Mukai 1988 [44], Garbagnati-Sarti 2007 [45] and Nikulin's involution classification [51], and the Betti formula b₂ = 0 + 21, b₃ = 22 + 55) closes the topological count exactly (Lean-formalized in `JoyceKarigiannisConstruction.lean`, introduced at core v3.4.14; the module was later migrated to the authors' archive and is available on request). The JK route is verified at the topological/lattice level; an explicit closed-form positive G₂-structure ansatz at the neck level, with hyperkähler rotation and five-layer Wirtinger certificate, is established in [D]. The companion analytic paper [E] discharges the datum-level analytic existence scheme at a normalised datum $\mathcal D_0$ from a $K_7 \to S^3$ construction with $N = 77$ round-unlink branch components: seven theorem-grade results assemble into a conditional branched lifting scheme (a special case of [Donaldson 2017, Conjecture 1]), source-side discharge constants are certified with structural provenance ($C_{\rm src} = 27/16$, $C_{\rm nl} = 2/3$, $C_{\rm link} \le 1$, $\gamma = 2$), and $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$ is enclosed by outward-rounded interval arithmetic. The datum-level layer is therefore *discharged conditionally* on the two-slot external structure pack $\mathsf H_{\rm global}$ of [E] (smooth topology of $K_7$ from classification results, and a hyperkähler K3-fibration with prescribed rank-one Picard-Lefschetz monodromy along the 77-unlink); an *anisotropic three-scale perturbation theorem* $(\mathrm J)$ [E, §8.3] would additionally deform each closed positive $3$-form to a torsion-free $G_2$-structure. On the metric side, the Ricci-flat Kähler (Calabi–Yau) metric on the Z₂³-equivariant K3 fibre now has an explicit finitely-parametrised closed-form approximation whose order-3 residual (667 parameters) satisfies an interval-rigorous, assumption-free, machine-checked bound Var(log R) ≤ ε₃ = 1309/10⁷ ≈ 1.309 × 10⁻⁴ < 10⁻³ on a frozen reproducible 4000-point witness (Krawczyk–Rump-certified K3 points evaluated in forward interval arithmetic, no floating-point assumption), formalized in Lean 4 (`K3ClosedFormWitness`, `native_decide`, zero `sorry`, zero added axiom); promoting this sample bound to a global bound over the entire K3 surface remains an open constraint-aware global-positivity (Positivstellensatz/SOS) problem. See §9 for the two-layer claim boundary between datum-level and global-analytic completion. Among the 65 *catalogued* TCS literature examples, (21, 77) remains an outlier (nearest neighbor at distance 7.6); the JK route therefore extends the known catalogue rather than fitting in it. The certified metric is 99.98% block-diagonal K3×T²×I and is numerically compatible with a U(1)² symmetry to tolerance $2 \times 10^{-5}$ (period integrals $S_\theta = S_\psi = 6.1265$, matched to $2.6 \times 10^{-8}$).

**Further structure**: The certified metric is smooth (no singularities). The framework encodes the SM gauge group SU(3)×SU(2)×U(1) at the algebraic/spectral level: g₂ ⊂ so(7) decomposition, so(8) = g₂ ⊕ L ⊕ R triality giving $N_{\rm gen} = 3$, and spectral gap $\lambda_1 = 6\pi^2/475 = 0.12467$ (Richardson extrapolation: $0.12461 \pm 0.00016$, deviation 0.05%) (§5). *Honest scope caveat.* In the standard M-theory dictionary for compactifications on smooth $G_2$-manifolds, non-abelian gauge symmetry and 4D chiral fermions are known to require appropriate ADE or codimension-seven singularities (Acharya–Witten, [hep-th/0109152]; Acharya–Gukov review [hep-th/0409191]). The framework does **not** yet supply an M-theory-style construction realising the chiral non-abelian SM sector from a smooth K₇; it supplies a conditional topological/algebraic dictionary whose phenomenological relations become physically interpretable only under an additional **UV bridge** (a heterotic gauge bundle, a duality identification, or a singular $G_2$ locus enrichment) that remains a load-bearing open question (§9).

The framework imposes det(g) = 65/32 = (Weyl × α_sum)/2⁵, equivalently 2 + 1/32. Joyce's theorem [20] guarantees existence of torsion-free G₂ metrics on suitable compact 7-manifolds; the Chebyshev-Cholesky construction (companion paper [A]) achieves ‖T‖ < 3 × 10⁻⁵ with Newton-Kantorovich certification h = 8.95 × 10⁻⁹ (Chebyshev-certified, margin ×56 million; §8.5).

### 2.4 Topological Constraints on Field Content

#### 2.4.1 Betti Numbers as Capacity Bounds

The Betti numbers provide upper bounds on field multiplicities:

- **b₂(K₇) = 21**: Bounds the number of gauge field degrees of freedom
- **b₃(K₇) = 77**: Bounds the number of matter field degrees of freedom

**Note on gauge group origin**: In M-theory on smooth G₂ manifolds, dimensional reduction yields b₂ abelian U(1) vector multiplets [11]. The non-abelian SM gauge group in the K₇ framework emerges instead from the algebraic/spectral structure of the G₂ holonomy: the g₂ ⊂ so(7) decomposition and so(8) = g₂ ⊕ L ⊕ R triality (§5). The smooth certified metric is at a generic point in G₂ moduli space; codimension-4 ADE singularities would be a different (singular) limit.

#### 2.4.2 Generation Number

The number of chiral fermion generations follows from a topological constraint:

$$({\rm rank}(E_8) + N_{\rm gen}) \times b_2 = N_{\rm gen} \times b_3$$

Solving: (8 + N_gen) x 21 = N_gen x 77 yields **N_gen = 3**.

This derivation is formal; physically, it reflects index-theoretic constraints on chiral zero modes, which in M-theory on G₂ require singular geometries for chirality [25].

---

## 3. Methodology and Epistemic Status

### 3.1 The Derivation Principle and Type Classification

The K₇ framework derives physical observables through algebraic combinations of 20 topological invariants (Appendix A). The 95 observables are organized by derivation directness:

| Type | Count | Derivation | Example |
|------|-------|-----------|---------|
| **I** | 33 | Direct algebra from topology | sin²θ_W = 3/13 |
| **II** | 19 | One physical identification step | m_u from ratio × VEV |
| **III** | 21 | Multi-step dynamical chains | Combined wilson_line+instanton lepton ratios |
| **IV** | 22 | Structural diagnostics | NK certification, Gram conditioning |

Type I predictions are dimensionless ratios of topological integers: they cannot be "fitted" and are either correct or wrong. Type II adds one scale identification. Type III involves dynamical mechanisms (gauge running, eigenvalue splitting, instanton volumes). Type IV provides internal consistency checks.

**Geometric unit convention.** Throughout this paper, we adopt the convention that the geometric (framework) values are the *reference*, and experimental measurements are expressed as offsets from the geometric prediction. This reversal parallels the 2019 SI redefinition: the kilogram is now defined from Planck's constant (exact), and any physical realization has uncertainty relative to it. In the K₇ framework, the topological formula is exact; the experimental value approximates it.

### 3.2 What the K₇ Framework Claims and Does Not Claim

**Inputs**: Existence of K₇ with G₂ holonomy and (b₂, b₃) = (21, 77); E₈×E₈ gauge structure; det(g) = 65/32.

**Outputs**: 95 observables across 4 types (33 Type I, 19 Type II, 21 Type III, 22 Type IV), 66 experimentally testable.

**Structure of predictions**: All 95 observables are algebraic functions of 6 primitive topological integers (b₂, b₃, dim(G₂), dim(E₈), rank(E₈), dim(K₇)) plus standard transcendentals (π, √2, ln 2, ζ, and the golden ratio φ = (1+√5)/2; see §4.1 for the McKay E₈ ↔ 2I origin of φ and its caveat). No observable depends on the detailed G₂ geometry (the metric certifies that the manifold exists at the seam / datum level, but does not enter the prediction formulas; caveat: 6 observables use det_num/det_den = 65/32, which is a metric normalization target with suggestive but not derivational algebraic expressions in terms of topological integers; see §10.3 of S1).

**Counting caveat.** Type I (33 relations) are algebraic identities in the ledger; they employ *no continuously adjustable parameter after the declared structural ledger is frozen*. Type II (19 relations) are **conditional reconstructions**: they multiply a Type I ratio by an *experimental anchor* (e.g. $m_u = (m_u/m_d) \times m_d^{\rm PDG}$). They should not be counted as 19 independent predictions of new observables; they are 19 traceable extractions of physical scales from Type I ratios given an experimental input. Type III (21 relations) are dynamical / geometric chains conditional on the mechanism they invoke (Wilson-line, instanton, combined pipeline; each explicitly enumerated). Type IV (22) are structural diagnostics. We claim that given the inputs, Type I follows algebraically; Types II, III, and IV are conditional in the senses just enumerated. We do **not** claim uniqueness of the geometry or uniqueness of the formula assignments; the (21, 77) selection question, previously left open in 3.4, is treated as an audit with a single explicit residual (§10).

### 3.3 Three Factors Distinguishing the K₇ Framework from Numerology

**Multiplicity**: 95 observables (66 testable), not cherry-picked coincidences. The **Sieve reading** (§7), a public 4-null battery followed by per-survivor budget-uniqueness ranking, with a Lean 4 R2 **formal-identity flag** for pre-registered relations (the flag certifies the conditional identity that Lean has been given, *not* that the relation was a priori), isolates m_H/m_W (Route A) as the unique survivor at exact rank 1 budget-unique after Type I calibration, and Koide's Q = 2/3 as a further survivor with narrow budget margin. Under the previously reported *coincidence-probability* framing (retained in Supplement S4 as a diagnostic), the joint probabilities were $\sim 10^{-346}$ (uniform null) and $\sim 10^{-133}$ (algebraic null on 4.2M random formulas from the same 20 constants); these are indicators of internal consistency rather than structural claims. The Sieve reading supersedes them as the headline methodology, consistent with the Arithmon methodology paper.

**Exactness**: Several predictions are exactly rational, sin²θ_W = 3/13, Q_Koide = 2/3, m_s/m_d = 20, Ω_DM/Ω_b = 43/8. These cannot be fitted; they are correct or wrong.

**Falsifiability**: DUNE (first beam targeted by 2031, per Fermilab 2026-05 milestone [47]; physics run late 2030s–2040s) will test δ_CP = 197° with expected resolution ranging from a few degrees to ~15° depending on exposure and true parameter values [30, 31]. The SUSY spectrum (m_gravitino = 166 GeV, m_moduli = 3.2 TeV) is model-dependent and constrained in standard simplified searches; viable only in compressed or suppressed-coupling realizations requiring dedicated recasting at HL-LHC/FCC. The proton lifetime τ_p = 4.06 × 10³⁸ years exceeds near-term Hyper-K sensitivity; Hyper-K can strengthen lower bounds and constrain nearby GUT-scale alternatives.

### 3.4 Why These Formulas?

The selection question has two levels of answer.

**Mathematical level (deterministic).** The NK-certified metric has zero free parameters; its 169 Chebyshev coefficients appear constrained to a lattice generated by {π², π, 1, e, χ₉₈} / (b₂·b₃) (a numerical observation, not a derived result; see S1 §10.3). Observables are functionals of this metric and are therefore algebraic in the topological invariants. Within the declared formula grammar and structural-constant set, the observed relations are highly constrained and statistically non-generic; the (21, 77) selection question is treated in §10 (five routes audited and closed; the question itself remains open pending the §10.2 residual).

**Statistical level (empirical).** Under the Sieve reading (§7), the 4-null battery (uniform / algebraic / factorised / permuted) is public and passes through the *same* declared 20 structural constants. Survivors are budget-uniqueness-ranked; those that carry a Lean 4 R2 pre-registered identity receive a **formal-identity flag** that raises their traceability and machine-checkability, but not their a priori probability of being physically correct (the identity remains conditional on Lean's input). Structural redundancy adds a third line of evidence: sin²θ_W = 3/13 admits 14 independent derivations, Q_Koide = 2/3 admits 20 (see S2): an overdetermined web inconsistent with post-hoc cherry-picking.

The deeper question (*why G₂ holonomy?*) has the same epistemic status as "why Lorentz invariance?": G₂ is the unique compact exceptional holonomy group admitting a 7-dimensional representation with the stabilizer chain G₂ ⊃ SU(3) ⊃ SU(2) ⊃ U(1) required by the Standard Model. It is an isolated point in the space of compatible structures, not one choice among a continuum.

### 3.5 Posture: orientation, not ontology

It is useful to distinguish three registers at which a framework of this type operates, because the support the framework offers is uneven across them.

**Predictive register.** The K₇ framework specifies a finite list of inputs and derives 95 observables from them; 66 of these are testable, with mean deviation ~1% and a sensitivity profile that places the joint configuration at $> 4.5\sigma$ against random algebraic null models (§7). This is the register where the framework either holds or fails, and it is the only register at which we make any positive claim.

**Architectural register.** The choice of G₂-holonomy 7-manifolds with $E_8 \times E_8$ structure is *motivated* by mathematical considerations (§2, §10): G₂ is the unique exceptional compact holonomy admitting the SU(3) × SU(2) × U(1) stabilizer chain, $E_8$ is the largest exceptional Lie algebra, and the topological constraints fix $(b_2, b_3) = (21, 77)$ and $N_{\mathrm{gen}} = 3$. We argue this architecture is *natural*; we do not argue it is *necessary*.

**Ontological register.** A reader may wonder whether the K₇ framework carries a thesis about reality: for instance, the speculation that geometry, information, and energy are not three correlated aspects of nature but three views of a single underlying configuration, or the resonance with Wheeler's "It from Bit" programme [39] and the holographic principle. Such readings are compatible with the framework and historically motivated its development, but they are **neither required nor demonstrated** by the predictions. A reader who finds the Wheeler-holographic picture persuasive will see the framework as a natural piece of that puzzle; a reader who prefers a strictly empirical stance will see a falsifiable predictive framework. Both readings are defensible; the framework requires neither.

In Worrall's [40] and Ladyman's [41] terms, the K₇ framework is best read as a *moderate structural-realist orientation*: structure carries predictive weight independently of any further ontological commitment. The framework's success or failure is decided by experiment in the predictive register, not by adjudication in the ontological one. A more extended discussion of this posture, written for a non-technical audience, appears in the companion essay *Orientation, not ontology* [essay].

---

## 4. Observables: 95 Relations from 20 Structural Constants

**Measurement conventions.** Unless a table states otherwise: experimental values are PDG 2024 [1]; quark masses and their ratios are $\overline{\rm MS}$ values at the PDG reference scales (2 GeV for light quarks, $m_Q$ for heavy); gauge couplings and sin²θ_W are quoted at $M_Z$ in $\overline{\rm MS}$ (the topological sin²θ_W = 3/13 is compared at $M_Z$; its GUT-scale reading is treated separately in §5.3); PMNS angles and δ_CP in the tables of §4.6 and S3 §3.4 use the NuFIT 6.0 dataset [27] (frozen 2024-10), with the δ_CP discussion of §4.2 additionally tracking the NuFIT 6.1 (2025) release; cosmological parameters are Planck PR4 (Tristram et al. 2024 [49]; see S3 §3.2 for the source table). Each comparison therefore carries an implicit (scheme, scale, dataset, freeze date) tag; departures from these defaults are flagged in place.

### 4.1 Type I: Direct Algebraic Relations (33)

The 33 Type I predictions derive directly from the 20 structural constants (Appendix A). All are dimensionless closed forms with no continuously adjustable parameter: they cannot be "fitted" and are either correct or wrong. Six of the 33 involve the imposed metric normalization det(g) = 65/32 (§2.3): λ_H, α⁻¹, Ω_c/Ω_Λ, Ω_b/Ω_m, Ω_c/Ω_m, σ₈; these are labelled STRUCTURAL rather than TOPOLOGICAL wherever a status column appears (below and in Supplement S2). All 33 are Lean-certified. Representative highlights:

| Observable | Formula | Prediction | Exp. | Dev. | Status |
|-----------|---------|------|------|------|--------|
| sin²θ_W | b₂/(b₃+dim(G₂)) = 21/91 | 3/13 | 0.23122 | 0.195% | TOPOLOGICAL |
| Q_Koide | dim(G₂)/b₂ = 14/21 | 2/3 | 0.666661 | 0.001% | TOPOLOGICAL |
| α⁻¹ | 128 + 9 + det(g)×κ_T | 137.033 | 137.036 | 0.002% | STRUCTURAL (uses metric normalization) |
| m_τ/m_e | 7 + 2480 + 990 | 3477 | 3477.15 | 0.004% | TOPOLOGICAL |
| m_s/m_d | p₂² × Weyl = 4 × 5 | 20 | 20.0 | 0.00% | TOPOLOGICAL |
| n_s | ζ(11)/ζ(5) | 0.9649 | 0.9649 | 0.004% | DERIVED |
| Ω_DM/Ω_b | (1+42)/8 | 43/8 | 5.375 | 0.00% | TOPOLOGICAL |

Status legend: TOPOLOGICAL = pure topological-integer ratio; DERIVED = closed-form involving standard transcendentals (π, √2, ln 2, ζ, golden ratio φ); STRUCTURAL = depends on the imposed metric normalization (det_num / det_den).

Two relations deserve individual mention, plus one caveat.

**Koide parameter.** The charged-lepton relation Q = (m_e + m_μ + m_τ)/(√m_e + √m_μ + √m_τ)² = 2/3 has resisted explanation since 1982 [6] and holds experimentally to six significant figures (Q_exp = 0.666661 ± 0.000007). The K₇ framework provides Q_Koide = dim(G₂)/b₂ = 14/21 = 2/3: an algebraic identity in two topological invariants, no fitting involved. Deviation **0.0009%**, the most precise agreement in the framework.

**Weinberg angle.** sin²θ_W = b₂/(b₃ + dim(G₂)) = 21/91 = 3/13 = 0.230769, vs 0.23122 ± 0.00004 (PDG 2024 [1]): deviation **0.195%**. The numerator counts gauge moduli; the denominator counts matter plus holonomy degrees of freedom.

**Golden-ratio caveat (m_μ/m_e).** m_μ/m_e = dim(J₃(𝕆))^φ = 27^φ = 207.01 (exp 206.768, deviation 0.118%), with φ = (1+√5)/2 arising from the McKay correspondence E₈ ↔ 2I (binary icosahedral group). φ is the *only non-integer input* among the 33 Type I relations and does not appear in the 20 structural constants of S3 §3.3. Its status is accordingly weaker than the other Type I derivations: the formula is algebraically exact and Lean-certified, but φ's derivation from E₈ structure requires the additional McKay step (S2 §11). We flag it separately so that it does not stand in for, or borrow from, the robustness of the integer-ratio relations.

The complete per-sector relation tables (gauge, lepton, quark, neutrino mixing, Higgs, boson, CKM, cosmology), with their derivational prose, are in **Supplement S3 §3.4**; complete derivations for all 33 in **Supplement S2**.

### 4.2 The δ_CP Prediction

The framework's prediction for the CP-violation phase is:

$$\delta_{CP} = \dim(K_7) \times \dim(G_2) + H^* = 7 \times 14 + 99 = 197°$$

decomposing into a local contribution (7 × 14 = 98, fiber-holonomy coupling) and a global contribution (H* = 99, cohomological dimension). This is the framework's canonical prediction, pure topological, zero corrections.

**Experimental status**: NuFIT 5.2 (2022) gave δ_CP ≈ 197°, an exact match. NuFIT 6.0 (Oct 2024) shifted the central value to 177° ± 20°, and NuFIT 6.1 (2025, NO w/o SK-atm) further shifted the central value to $207^{+23}_{-20}$°. Relative to NuFIT 6.1, the prediction 197° now sits *inside* the 1σ band (deviation 4.8% from the 6.1 central value). The experimental uncertainty is still large (~$\pm 20$°) and the central value may shift further as T2K, NOvA, and DUNE accumulate statistics; the 2025 T2K+NOvA joint analysis [26] is the main new input driving the recent NuFIT shifts.

A post-hoc structural observation involving a possible compactification factor 62/69 is recorded in S2 Appendix F for completeness; it is **not** part of the main framework and is not used to revise the canonical 197° prediction.

**Falsification criterion**: If DUNE measures δ_CP outside [182, 212]° at 3σ, the framework faces serious tension. Against NuFIT 6.1 (207$^{+23}_{-20}$°) the prediction 197° now sits inside the 1σ band; against NuFIT 6.0 (177° ± 20°) it sat at the edge of the band. The prediction has not changed; the experimental central value has drifted, and the two most recent NuFIT releases bracket it on both sides.

### 4.3 Type II: Extended Algebraic Predictions (19)

Type II observables require one physical identification step beyond Type I ratios. Representative results:

**Absolute quark masses** (via m_q = ratio × reference_mass): m_u = 2.16 MeV (0.00%), m_d = 4.67 MeV (0.22%), m_s = 93.4 MeV (0.22%), m_c = 1.27 GeV (0.00%), m_b = 4.18 GeV (0.00%), m_t = 172.7 GeV (0.01%).

**CKM magnitudes** (from Wolfenstein parametrization): |V_us| = 0.2253 (0.13%), |V_cb| = 0.0412 (0.24%), |V_ub| = 0.00365 (0.27%), |V_td| = 0.0087 (0.34%), |V_tb| = 0.9991 (0.00%).

**Extended ratios**: m_c/m_s = 246/21 = 11.714 (exp 11.7, dev 0.12%), m_c/m_d = 234.3 (exp 234.0, dev 0.12%), m_μ/m_τ = 5/84 (exp 0.0595, dev 0.04%).

Type II mean deviation: **0.17%** across 19 observables. These inherit the precision of Type I ratios, with the physical identification step contributing negligible additional error.

### 4.4 Type III: Dynamical Predictions (21)

Type III observables involve multi-step dynamical mechanisms, conditional on explicitly-enumerated mechanism and calibration choices (see the status note at the top of §6 and the counting caveat in §3.2). They are grouped by computation:

**wilson_line Non-adiabatic** (3 obs): Wilson line eigenvalue splitting on K3 fiber at c = 0.452 gives raw lepton mass ratios with 0.5–2.1% deviation. These are improved to < 0.4% by the combined wilson_line+instanton pipeline (§6.4).

**RGE_running RGE running** (4 obs): Two-loop MSSM evolution from M_GUT to M_Z. The topological sin²θ_W = 3/13 at M_GUT runs to 0.2377 at M_Z (exp 0.2312, dev 2.78%). The strong coupling α_s(RGE) = 0.1224 deviates 3.7% using G₂-MSSM split-spectrum matching (§5.3). M_GUT = 2 × 10¹⁶ GeV is an exact match.

**spectral Spectral** (5 obs): effective Weyl law exponent (from the adiabatic seam-sector decomposition, extrapolated to d_eff = 7) α = 3.460 (exp 3.5, dev 1.1%), 22,671 KK states below cutoff, 57,578 fiber channels, Poisson level spacing. See [B] for the direct seam-sector result α = 1.998.

**gauge_bundle Gauge bundle** (4 obs): cond(f_IJ) = 1.047 (near-perfect gauge universality), α_ratio = 1.000002, effective Yukawa rank = 3, and κ(gauge) = 1.047 (4.7% departure from exact universality, the largest Type III bundle deviation).

**instanton Instanton + combined** (5 obs): Associative volume differences give ΔV(e-τ) = 8.633 (dev 5.9%), ΔV(e-μ) = 3.271 (dev 15.9%). Combined wilson_line+instanton pipeline with α = e^K (geometric, §6.4) gives τ/e = 3485 (0.23%), τ/μ = 16.69 (0.75%), μ/e = 208.8 (0.98%).

Type III mean deviation: **3.4%** across the 14 experimentally comparable observables; the remaining 7 are internal structural targets (e.g. M_GUT, N_KK, rank(Y)) with no independent experimental comparison, and are excluded from the statistic. Details in §5 (Gauge), §6 (Mass Hierarchy). (M_res and N_QNM from Pinčák et al. 2026 [42] are classified Type IV, see §4.5.)

### 4.5 Type IV: Structural Diagnostics (22)

Type IV observables are internal consistency checks with no experimental comparison:

- **Topology**: b₂ = 21, b₃ = 77, χ(K₇) = 0, H* = 99
- **Newton-Kantorovich certification**: h = 8.95 × 10⁻⁹ (< 0.5, margin ×56M), δg/g ≤ 4.86 × 10⁻⁶
- **Gram conditioning**: cond(G_K3) = 1.05, cond(G_K7) = 1.05, cond(G_35) = 7.66, G_77 positive definite
- **Spectral counts**: 22,671 KK states, 57,578 fiber channels, Poisson spacing confirmed
- **Torsion**: ‖T‖_C⁰ = 2.949 × 10⁻⁵ (×2995 reduction in 5 Joyce steps)
- **Metric eigenvalues**: g_ss = 19/6, g_{T²} = 7/6, g_{K3} ≈ 64/77
- **Instanton & BH diagnostics** (Pinčák et al. 2026 [42]): N_QNM = 98 QNM mode families, b₃/b₃(S³×S⁴) = 77× instanton suppression, M_res = v_EW²/M_Pl (BH remnant mass: no experimental comparison)

These diagnostics confirm the geometric construction is well-conditioned and internally consistent.

### 4.6 Summary Statistics (All 95 Observables)

**Global performance** (95 observables):

| Metric | Type I | Type II | Type III | All (I+II+III) |
|--------|--------|--------|--------|-------------|
| Count | 33 | 19 | 21 | 73* |
| With exp. comparison | 33 | 19 | 14 | 66 |
| Mean deviation | 0.73% | 0.17% | 3.44% | ~1.15% |
| Median deviation | 0.23% | 0.12% | 1.39% | 0.23% |
| Exact matches (<0.01%) | 5 | 3 | 3 | 11 |
| Within 1% | 28 | 19 | 6 | 53 |
| Maximum deviation | 11.3%† | 0.79% | 15.9% | 15.9% |

†δ_CP = 11.3% dominates Type I. The figure is computed against the *frozen* NuFIT 6.0 central value (177°, per the §4 conventions note); against the NuFIT 6.1 central (207°) the same prediction deviates by 4.8% and sits inside the 1σ band, see §4.2. Excluding δ_CP, Type I mean = 0.40% (28 within-1%, 5 exact matches, max 2.77%).

*66 with experimental comparison out of 73 (I+II+III); 22 Type IV observables are structural.

**Sector breakdown** (11 sectors; 40 observables listed, see §5–6 for remaining 26 comparable):

| Sector | N_obs | Mean dev | Best | Worst |
|--------|-------|----------|------|-------|
| Electroweak | 3 | 0.11% | α⁻¹ 0.002% | sin²θ_W 0.20% |
| Boson | 3 | 0.13% | m_H/m_W 0.02% | m_H/m_t 0.31% |
| Lepton | 3 | 0.04% | Q_Koide 0.001% | m_μ/m_e 0.12% |
| Quark | 4 | 0.24% | m_s/m_d 0.00% | m_b/m_t 0.79% |
| Cosmology | 7 | 0.15% | n_s 0.004% | Y_p 0.37% |
| PMNS | 4 | 0.29% | θ₁₂ 0.03% | sin²θ₁₃ 0.81% |
| CKM | 3 | 0.59% | A_Wolf 0.29% | sin²θ₂₃ 1.13% |
| Gauge (running) | 4 | 2.3% | M_GUT 0.00% | α_s(RGE) 3.7% |
| Instanton | 2 | 10.9% | ΔV(e-τ) 5.9% | ΔV(e-μ) 15.9% |
| Combined | 3 | 0.66% | τ/e 0.23% | μ/e 0.98% |
| Bundle | 4 | 1.6% | α_ratio 0.00% | κ(gauge) 4.7% |

---

## 5. Gauge Sector: From E₈ to the Standard Model

The gauge sector derives the Standard Model gauge group from the E₈×E₈ structure of heterotic M-theory on K₇. This section presents the complete breaking chain, anomaly cancellation, gauge coupling running, and bundle universality, all from the topological data of K₇ and an explicit 7D G₂ metric (169 optimized Chebyshev-Cholesky parameters capturing the dominant seam sector, K3 fiber certified at 0.07% [C]), with no free parameters in the physical predictions.

### 5.1 The E₈ → Standard Model Breaking Chain

![**Figure 1.** Representation-branching chain from $E_8$ (dim 496) to the Standard-Model gauge group $SU(3)_c \times SU(2)_L \times U(1)_Y$ (plus a hidden sector). Each arrow is a subgroup embedding compatible with the SM matter content; this is a chain of representation branchings, not a physical breaking mechanism. Three generations arise from the $\mathbb{Z}_3$ lattice action combined with the $(27, 3)$ branching (§6.1), with the Wilson-line rank-2 SVD as auxiliary data. The UV bridge to the non-abelian chiral gauge sector (Acharya–Witten singularity or heterotic bundle) is load-bearing open (§9.3).](fig_e8_branching_chain.pdf)

The first E₈ factor breaks to the Standard Model through a six-level chain:

| Level | Group | Dimension | Mechanism | Scale |
|-------|-------|-----------|-----------|-------|
| 0 | E₈×E₈ | 496 | Heterotic structure | M_Pl |
| 1 | E₈ | 248 | Second E₈ = hidden sector | M_string |
| 2 | E₆ × SU(3) | 78 + 8 = 86 | Adjoint branching | M_string |
| 3 | SO(10) × U(1)_1 | 45 + 1 = 46 | E₆ ⊃ SO(10) × U(1) branching | M_GUT |
| 4 | SU(5) × U(1)_2 | 24 + 1 = 25 | SO(10) ⊃ SU(5) × U(1) branching (Georgi–Glashow-type; distinct from Pati–Salam $SU(4) \times SU(2)_L \times SU(2)_R$) | M_GUT |
| 5 | SU(3) × SU(2) × U(1) | 8 + 3 + 1 = 12 | Standard Model | M_Z |

The E₈ adjoint decomposes under E₆ × SU(3) as:

$$248 = (78, 1) + (1, 8) + (27, 3) + (\overline{27}, \overline{3})$$

yielding 78 + 8 + 81 + 81 = 248. The fundamental of SU(3) has dimension 3, giving **N_gen = 3** chiral families from the (27, 3) representation. This reproduces the index-theorem result of §2.

**Fundamental group**: $b_1(K_7) = 0$ constrains only the abelianisation of $\pi_1(K_7)$; the framework additionally *postulates* $\pi_1(K_7) = \{1\}$ as a construction-level input (satisfied by the $K_7 \to S^3$ scheme of [E], Supplementary Appendix C, where simple connectivity is established via Van Kampen on the K3-fibration + local Kovalev–Lefschetz models). Traditional Wilson line breaking via $\pi_1$ is therefore trivial. Instead, the breaking proceeds through the Z₃ lattice action on the K3 fiber (§6.1).

**Scales**: M_Pl = 1.22 × 10¹⁹ GeV, M_string = 4 × 10¹⁷ GeV, M_GUT = 2 × 10¹⁶ GeV, M_Z = 91.19 GeV. The scale hierarchy M_GUT/M_Z ~ 2 × 10¹⁴ emerges from the geometry without fine-tuning. (The 1-loop exact-unification identity of §5.4 gives M_GUT = 2.7 × 10¹⁶ GeV; both values sit within the standard MSSM range, and the discrepancy reflects the approximate nature of that identity.)

### 5.2 Anomaly Cancellation

All six Standard Model gauge anomalies vanish:

| Anomaly | Value | Status |
|---------|-------|--------|
| SU(3)³ | 0.0 | Exact |
| SU(2)³ | 0.0 | Exact |
| U(1)³ | 2.3 × 10⁻¹⁶ | Machine zero |
| grav-U(1) | 0.0 | Exact |
| SU(3)²-U(1) | 0.0 | Exact |
| SU(2)²-U(1) | 0.0 | Exact |

The Green-Schwarz mechanism provides additional consistency:
- χ(K₇) = 0 (tadpole cancellation for compact odd-dimensional manifolds)
- Tadpole value = 0 (no net D-brane charge)
- G₄ flux lattice rank = 77 = b₃(K₇) (all flux degrees of freedom realized)

The 10/10 verification certificate confirms all anomaly, tadpole, and lattice checks pass. Lean: `TCSGaugeBreaking.lean` (10 conjuncts, 0 sorry).

### 5.3 Gauge Coupling Running

**Boundary conditions at M_GUT**: α_GUT⁻¹ = 25.3, sin²θ_W = 3/13 = 0.23077. These are topological: α_GUT derives from the gauge kinetic function on K₇, and sin²θ_W from the Betti number ratio b₂/(b₃ + dim(G₂)).

**Two-loop MSSM RGE** (tan β = 2, M_SUSY = m_moduli = 3165 GeV):

| Coupling | Prediction at M_Z | Experimental | Deviation |
|----------|-------------|--------------|-----------|
| sin²θ_W | 0.2377 | 0.23122 | 2.78% |
| α_em⁻¹ | 131.19 | 127.951 | 2.53% |
| α_s (all-MSSM) | 0.1038 | 0.1180 | 12.1% |
| **α_s (split-spectrum)** | **0.1224** | **0.1180** | **3.7%** |

A naive degenerate-spectrum treatment (all superpartners at M_SUSY) yields α_s = 0.1038 (12.1% deviation). The physical G₂-MSSM spectrum is split: squarks and sleptons decouple at M_squark = 3165 GeV, while gauginos remain light. Below M_squark, the effective theory is SM + gauginos with b₃ = −5 (instead of −3 for the full MSSM). This step-function matching gives α_s = 0.1224 (3.7% deviation), which is the value adopted throughout.

The RGE deviations reflect sensitivity to the SUSY spectrum, particularly:

- **tan β dependence**: The 2-loop terms involving top Yukawa are sensitive to tan β
- **Threshold corrections**: The split-spectrum matching captures the dominant effect; residual corrections from the detailed spectrum are subdominant
- **3-loop effects**: Omitted from the current calculation

**Topological vs. infrared**: The topological sin²θ_W = 3/13 = 0.23077 (dev 0.19% from PDG) is more accurate than the RGE-evolved value (dev 2.78%). This suggests the topological value may represent an infrared fixed point rather than a UV boundary condition: an observation also made in the G₂-MSSM literature [33].

**KK threshold corrections**: With 71 KK modes above M_GUT and Σ ln(m/M_GUT) = 89.2, the net correction to α_GUT⁻¹ is 0.0 (smooth K₇ manifolds have vanishing KK threshold corrections due to the cancellation between towers).

**SUSY spectrum implications**:

| Particle | Mass | Source |
|----------|------|--------|
| m_gravitino | 166 GeV | F-term from gaugino condensation |
| m_moduli | 3165 GeV (3.2 TeV) | SU(8) gaugino condensation |
| m_gluino | 442 GeV | cMSSM with m₁/₂ = m_gravitino |
| m_squark | 1094 GeV | cMSSM (running mass at low scale; decoupling matched at M_SUSY = m_moduli = 3165 GeV) |
| m_slepton | 175 GeV | cMSSM |
| LSP (Bino) | 70 GeV | Lightest neutralino (pure Bino; viable: LEP chargino bound 103.5 GeV [46] does not exclude pure Bino LSP with suppressed gauge couplings) |

**Phenomenological caveat**: These masses are computed within the cMSSM approximation. Standard simplified SUSY searches at ATLAS/CMS already exclude portions of this parameter space; the spectrum is viable only in compressed or suppressed-coupling realizations. Definitive testing requires a dedicated recast of current ATLAS/CMS results against the G₂-MSSM split-spectrum scenario.

### 5.4 The B-Test Identity and Holonomy Sequence

The gauge coupling predictions sin²θ_W = 3/13 and α_s = √2/12, combined with the MSSM beta-function structure, yield an algebraic identity connecting the fine structure constant to G₂ representation theory.

**The B parameter**: In any GUT framework, the "B-test" quantifies consistency of gauge coupling unification:

$$B = \frac{\alpha_1^{-1} - \alpha_2^{-1}}{\alpha_2^{-1} - \alpha_3^{-1}}$$

For the MSSM with N_gen = 3 generations and one Higgs doublet pair, the beta-function coefficients are (b₁, b₂, b₃) = (33/5, 1, −3), giving B = (b₁−b₂)/(b₂−b₃) = (28/5)/4 = **7/5**. The number of generations enters through N_gen = b₂/dim(K₇) = 21/7 = 3.

**Theorem** (B-test): *Given sin²θ_W = b₂/(b₃+dim(G₂)) = 3/13 and α_s = √2/(dim(G₂)−p₂) = √2/12, the MSSM relation B = 7/5 holds if and only if*

$$\alpha_{\mathrm{em}}^{-1}(M_Z) = (b_3 + \dim(G_2)) \cdot \sqrt{2} = 91\sqrt{2} = 128.693\ldots$$

*where 91 = dim(Λ²𝔤₂) is the dimension of the exterior square of the G₂ Lie algebra.*

*Proof sketch*. The topological sin²θ_W = 3/13 and the GUT normalization factor 3/5 force α₁⁻¹ = 2α₂⁻¹. Then B = α₂⁻¹/(α₂⁻¹ − α₃⁻¹), and B = 7/5 requires α₂/α₃ = 2/7 = p₂/dim(K₇). Substituting α₃⁻¹ = (dim(G₂) − p₂)/√2 = 6√2 gives α_em⁻¹ = (7 × 13)√2 = 91√2. The factor 91 = b₃ + dim(G₂) = 77 + 14 = dim(Λ²𝔤₂) is a G₂ representation-theoretic invariant. *(Note: the B-test identity gives the GUT-scale α_em⁻¹ ≈ 128.7, valid at M_GUT where B = 7/5; the observed low-energy value α_em⁻¹ ≈ 137 is a separate derivation via §4.1. Both are used consistently in the framework.)*

**The holonomy sequence**: At the B = 7/5 exact scale, all three inverse couplings are integer multiples of √2:

| Coupling | Value | = integer × √2 | Topological origin |
|----------|-------|-----------------|-------------------|
| α₁⁻¹ | 59.40 | 42√2 | 2b₂ |
| α₂⁻¹ | 29.70 | 21√2 | b₂ |
| α₃⁻¹ | 8.49 | 6√2 | dim(K₇)−1 |

Their ratios in lowest terms:

$$\alpha_1^{-1} : \alpha_2^{-1} : \alpha_3^{-1} = \dim(G_2) : \dim(K_7) : p_2 = 14 : 7 : 2$$

This **holonomy sequence** encodes the G₂ structure of K₇ directly in the gauge couplings: the holonomy group dimension, the manifold dimension, and the Pontryagin ratio p₂ = dim(G₂)/dim(K₇).

**Exact unification**: Using the 1-loop MSSM RGE α_i⁻¹(μ) = α_i⁻¹(M_Z) − b_i/(2π)·ln(μ/M_Z), the three couplings unify exactly at:

$$t_{\mathrm{GUT}} = \frac{15\pi\sqrt{2}}{2}, \qquad \alpha_{\mathrm{GUT}}^{-1} = \frac{\sqrt{2}\,(b_3 - \mathrm{rank}(E_8))}{4} = \frac{69\sqrt{2}}{4} \approx 24.4$$

where 69 = b₃ − rank(E₈) = 77 − 8 = |PSL(2,7)| − H* = 168 − 99. This gives M_GUT = M_Z·e^{t} = 2.7 × 10¹⁶ GeV, consistent with the standard MSSM estimate. The GUT coupling α_GUT⁻¹ = 24.4 compares to the gauge kinetic function value 25.3 used in §5.3 (3.6% discrepancy, reflecting the approximate nature of the identity).

**Numerical status**: Using the framework's topological couplings at M_Z, B = 1.4033 (0.23% from 7/5). This is closer to 7/5 than the purely experimental value B = 1.3948 (0.37% off), suggesting the identity has genuine geometric content. Two distinct α_em⁻¹ values appear in this framework and should not be conflated: the B-test identity 91√2 = 128.69 is a GUT-scale quantity (§4.1 note), while the two-loop RGE chain of §5.3 gives α_em⁻¹(M_Z) = 131.19 at the electroweak scale [both distinct from the topological low-energy value α_em⁻¹ = 137.033]. Reconciling the algebraic GUT-scale identity with the end-to-end RGE chain remains an open question: the gap may trace to threshold effects and the SUSY-spectrum sensitivity discussed in §5.3.

### 5.5 Bundle Universality and Gram Conditioning

The gauge kinetic function f_IJ on K₇ determines gauge coupling universality. From the 22 harmonic 2-forms on K₃:

**f_IJ eigenvalue spectrum**: 22 eigenvalues in [0.733, 0.767], with **cond(f_IJ) = 1.047**, near-perfect universality. The gauge coupling ratio α_ratio = α_max/α_min = 1.000002, confirming that all gauge couplings are effectively identical at the compactification scale.

**Gram matrices** quantify orthonormality of the harmonic form basis:

| Matrix | Size | Condition | PD | Off-diag max |
|--------|------|-----------|----|----|
| G_K3(22) | 22×22 | **1.05** | Yes | 0.012 |
| G_K7(22) | 22×22 | **1.05** | Yes | 0.012 |
| G_35 | 35×35 | **7.66** | Yes | n/a |
| G_77 | 77×77 | **7.66** | Yes | 7×10⁻⁵ |

The K3 and K₇ 22-form bases are nearly orthonormal (condition ~1.05). The full 77-form basis has moderate conditioning (7.66), with cross-block coupling between constant and fiber modes bounded by 7 × 10⁻⁵. Gram-Schmidt orthogonalization residuals are < 5 × 10⁻¹⁶ (machine precision).

Gauge-bundle diagnostics (eigenvalue spectrum of $f_{IJ}$, Gram-matrix conditioning) are provided as a repository figure (`fig_gauge_bundle.png`, `canonical/figures/`), not reproduced in this PDF.

**Lean certification**: `TCSGaugeBreaking.lean` (0 axioms, 14 theorems, 10-conjunct master certificate) + `GaugeBundleData.lean` (0 axioms, 12 theorems, 11-conjunct master certificate).

### 5.6 Summary

The gauge sector pipeline covers E₈ → SM with N_gen = 3, all anomalies cancelled, near-perfect bundle universality (cond 1.047), and Lean-certified results. The topological gauge predictions (sin²θ_W, α_s) are more precise than the dynamical RGE running, suggesting the topological values may represent infrared fixed points rather than UV boundary conditions. The B-test identity (§5.4) reveals that these two predictions, combined with the MSSM structure, encode the holonomy sequence dim(G₂):dim(K₇):p₂ = 14:7:2 in the gauge coupling ratios: a direct imprint of G₂ geometry on low-energy physics.

---

## 6. Mass Hierarchy: From Geometry to Generations

The five-order-of-magnitude lepton mass hierarchy (m_e : m_μ : m_τ ~ 1 : 207 : 3477) is one of the deepest puzzles in particle physics. This section presents two independent geometric mechanisms that individually reproduce the hierarchy to ~2–6% and, when combined, achieve sub-percent precision.

**Status: conditional calibration (Type III).** The mechanisms of this section consume explicit calibration choices: the non-adiabatic coupling c = 0.452 and the optimized positions [0.0, 0.693, 1.400] of §6.2, and the generation-to-cycle assignment of §6.3, selected among the 57 associative cycles by minimizing the combined deviation. The search spaces, objective function, and selection protocol are documented with the dataset (S3 §3.8; available from the author on request pending public deposit). The no-continuously-tuned-parameter guarantee of the abstract applies to the Type I ledger only; the results of this section are **conditional reconstructions** given these declared mechanism choices, not independent predictions (see the counting caveat in §3.2). The one cross-mechanism normalization, α = e^K (§6.4), is a geometric quantity, not a fitted one.

### 6.1 Three Generations from the Z₃ Mechanism

Since π₁(K₇) = {1}, traditional Wilson line breaking is trivial. Instead, three generations emerge from a Z₃ lattice action on the K3 fiber of the neck region.

**Wilson line theorem** (numerical statement). The SVD of the fiber-level Wilson-line operator yields singular values $[5.71,\ 0.62,\ 2.4 \times 10^{-15}]$: the numerical rank is **2**, with the third singular value at machine zero. This is the operator whose two nontrivial modes lift the fiber degeneracy along two independent directions inside the K3 fiber (the third mode carries no perturbative signal within numerical tolerance). The counting of three lepton generations does *not* come from a rank-3 Wilson-line operator on the K3 fiber; it comes from the Z₃ lattice action on the fiber (§6.1) together with the (27, 3) branching of §5.1. We flag this as a language correction relative to the 3.4 draft, which conflated the numerical rank of the SVD with the physical generation count.

**K3 metric properties**: The K3 fiber metric is nearly flat: conformal range 0.018% across the fiber, mean anisotropy 1.4%. This near-flatness ensures the eigenvalue splitting is controlled by the fiber geometry rather than by large-scale metric fluctuations.

### 6.2 Lepton Mass Hierarchy: Non-Adiabatic Mechanism (wilson_line)

The non-adiabatic eigenvalue splitting mechanism operates on the K3 fiber at coupling c = 0.452 and optimized positions [0.0, 0.693, 1.400]:

**Eigenvalues**: [0.03383, 0.00205, 9.94 × 10⁻⁶]

| Ratio | Prediction | Experimental | Deviation |
|-------|------|-------------|-----------|
| m_τ/m_μ | 16.54 | 16.82 | 1.7% |
| m_τ/m_e | 3403 | 3477 | 2.1% |
| m_μ/m_e | 205.7 | 206.7 | 0.5% |

The critical coupling c* = 10⁻³/⁴ = 0.1778 marks the transition between adiabatic (small splitting, ~2 generations) and non-adiabatic (large splitting, 3 generations) regimes. The physical coupling c = 0.452 > c* places K₇ firmly in the three-generation regime.

**Adiabatic limit**: At c → 0, only two generations are distinguishable (m₁/m₂ = 77.6, m₁/m₃ → ∞). The three-generation structure requires c > c*, which is satisfied by the neck geometry.

### 6.3 Instanton Volume Differences (instanton)

An independent mechanism generates the mass hierarchy from associative 3-cycle volumes (calibrated submanifolds in the sense of Harvey-Lawson [23]). On K₇, there are **57 associative 3-cycles** with volumes in [0.00075, 11.109]. The mass relation m_i/m_j = exp(ΔV) assigns each generation to a cycle.

**Optimal assignment** (minimizing combined deviation):

| Assignment | Volume | ΔV | Target | Deviation |
|------------|--------|-----|--------|-----------|
| V_e | 11.109 | n/a | n/a | n/a |
| V_μ | 7.838 | ΔV(e-μ) = 3.271 | ln(16.82) = 2.823 | 15.9% |
| V_τ | 2.476 | ΔV(e-τ) = 8.633 | ln(3477) = 8.154 | 5.9% |

The e-τ hierarchy (5 orders of magnitude) is reproduced to 5.9%, while the e-μ hierarchy (2.3 orders) shows 15.9% deviation from the optimal assignment. The total volume range ΔV_range = 8.92 spans the correct order of magnitude.

The instanton mass-hierarchy diagnostics (volume spectrum of the 57 associative cycles with generation assignments) are provided as a repository figure (`fig_instanton_hierarchy.png`).

### 6.4 Combined wilson_line+instanton Pipeline

The key insight is that wilson_line and instanton probe different aspects of K₇ geometry:
- **wilson_line**: Fiber geometry → eigenvalue spacing (relative structure)
- **instanton**: Cycle volumes → exponential hierarchy (absolute scale)

The two mechanisms are connected by α = e^K = exp(K₀) = V̂^{−3}, where K₀ = −5.891 is the Kähler potential of K₇ (§6.5). This is a purely geometric quantity (the instanton action normalization derived from the compactification volume) not a fit parameter:

| Ratio | wilson_line raw | instanton raw | Combined (α = e^K) | Experimental | Dev |
|-------|---------|---------|----------|-------------|------|
| m_τ/m_e | 3403 (2.1%) | exp(8.63) | **3485** | 3477 | **0.23%** |
| m_τ/m_μ | 16.54 (1.7%) | exp(3.27) | **16.69** | 16.82 | **0.75%** |
| m_μ/m_e | 205.7 (0.5%) | n/a | **208.8** | 206.8 | **0.98%** |

All three ratios within 1%: a significant improvement over the individual mechanisms. The combined pipeline works because the mechanisms are **complementary**: wilson_line provides fine structure from the eigenvalue splitting, instanton provides the overall exponential scale from cycle volumes. The key insight is that α = e^K has a natural M-theory interpretation: for M2-branes wrapping associative 3-cycles, the instanton action scales as S_inst = e^K × Vol(Σ) = Vol(K₇)^{−3} × Vol(Σ), giving the correct suppression without any free parameters.

**Lean certification**: `AssociativeVolumes.lean` (0 axioms, 19 theorems, 14-conjunct master certificate) certifies the combined wilson_line+instanton results including all three mass ratios and the geometric α = e^K.

### 6.5 4D Effective Theory

Dimensional reduction of the G₂ metric yields an N = 1 supergravity theory in 4D:

**Particle content**: From Betti numbers directly:
- 1 gravity multiplet
- 21 vector multiplets (= b₂)
- 77 chiral multiplets (= b₃)
- Total: 99 = H* (the effective cohomological dimension)

**Kähler potential**: K₀ = −5.891, with e^K = 0.002763 = V̂^{−3}. This geometric quantity serves as the instanton action normalization α = e^K in the combined wilson_line+instanton pipeline (§6.4). The 7D volume V̂ = 7.126. The metric determinant det(g₇) = 65/32 = 2.03125 is locked to its topological value.

**Moduli metric**: Total dimension 79 (= 77 physical + 2 gauge), rank 77 (2 null directions from gauge redundancy). The condition number cond = 7.66 measures the ratio of largest to smallest moduli mass. The 35 constant-block eigenvalues span [1.66, 12.69]; the 42 fiber-block eigenvalues cluster around 6.1.

**Gaugino condensation**: The SU(8) hidden sector is preferred:

| Condensate | b₀ | Λ (GeV) | m_moduli (GeV) |
|------------|-----|---------|----------------|
| SU(3) | 9 | 4.3 × 10⁸ | 1.3 × 10⁻¹¹ |
| SU(5) | 15 | 5.0 × 10¹¹ | 0.021 |
| **SU(8)** | **24** | **2.66 × 10¹³** | **3165** |

The SU(8) sector gives m_moduli = 3165 GeV (3.2 TeV), m_gravitino = 166 GeV, consistent with the G₂-MSSM spectrum.

**Physical predictions from the effective theory**:
- Proton lifetime: τ_p = 4.06 × 10³⁸ years (well above current bound ~10³⁴ years)
- Yukawa effective rank = 3 (three massive generations)
- λ_physical = 1.358 (physical Yukawa coupling)

The Kaluza-Klein tower structure and 7D Weyl-law diagnostics are provided as a repository figure (`fig_kk_spectrum.png`).

---

## 7. Statistical Uniqueness (Sieve reading)

![**Figure 2.** Sieve reading pipeline (methodology swap 3.4 $\to$ 3.5). The 33 conditional Type I identities of the ledger are pushed through a public 4-null battery (uniform / algebraic / factorised / permuted); certified survivors carry a per-survivor budget-uniqueness rank; the Lean 4 core (v3.4.29, commit `667c8b9`) attaches a **formal-identity flag** to each survivor with a pre-registered R2 identity. The Sieve is a distinguishability test, not a probability statement; the 3.4 coincidence-probability exponents ($10^{-346}$ uniform / $10^{-133}$ algebraic) are archived as diagnostics in Supplement S4.](fig_sieve_pipeline.pdf)

A framework claiming 95 observables from 20 structural constants must address the question: how constrained are these predictions? This section is organised around the **Sieve reading** (the methodology of the Arithmon companion paper on counting coincidences, Zenodo 10.5281/zenodo.20666879); the previously-headlined coincidence-probability figures ($\sim 10^{-346}$ / $\sim 10^{-133}$) are demoted to a diagnostic in Supplement S4 (see the CHANGELOG for the version delta).

**Sieve reading in one paragraph.** (Terminology: R1/R2 are the rebate conditions of the Arithmon methodology paper: R1 = the relation holds as a machine-certified identity; R2 = it was additionally pre-registered as a Lean 4 identity before the sieve was run.) A relation is stated as a *R2 pre-registered identity* (Lean 4 signature `Sieve/GStruct.lean`, budget-uniqueness rank certified via `native_decide`). It is then passed through a public four-null battery: (N1) uniform floats in $[0, 50\%]$; (N2) random algebraic formulas in the same declared 20 structural constants (haystack sizes $|E_{\rm rat}| = 32, 91, 1416, 9782, 105329$ for depth $k = 1, \ldots, 5$, Lean-certified); (N3) factorised nulls (route grammar preserved, coefficients shuffled); (N4) permutation nulls (targets shuffled within sector). Survivors of all four nulls are ranked by *budget-uniqueness rank* (the smallest complexity budget under which the survivor is the unique passing relation). A survivor carrying a R2 pre-registered identity receives a **formal-identity flag**: its sieve-passage is machine-checked as a conditional identity (given the ledger inputs Lean is fed), rather than merely as a coincidence-probability tail. The flag increases traceability, not physical prior. Full construction: Arithmon methodology paper §5–6.

**Headline outcome at 3.5.** Under the sieve applied to the 33 Type I relations (uniform-null passage baseline, then N2–N4), the unique survivor at exact rank 1 budget-unique is $m_H/m_W$ (Route A, R1 audit note $r1\_audit\_ns\_2026\_06\_14$); Koide's $Q = 2/3$ passes the battery with a narrow budget margin; $n_s$ passes without a distinguishing budget signature. The overdetermination structure of the framework (§7.3), the cross-coupling web (§7.4), and the PSLQ residual analysis (§7.6) remain intact as independent diagnostics; they are not the headline any more.

Accordingly, the framework's statistical claim is layered: the 33 Type I relations are **conditional algebraic identities of the frozen ledger** (machine-checked, but not individually sieve-distinguished); $m_H/m_W$ and, with narrow margin, Koide's $Q = 2/3$ are the **sieve-distinguished survivors**; the remaining relations stand as exploratory identities relative to the declared grammar. Statements of the form "33 predictions" elsewhere in the paper are to be read through this layering.

**Pre-registration traceability.** The term *pre-registered* is backed by frozen public artifacts, not by narrative: the structural-constant ledger, formula grammar, and full target list are frozen in the Sieve repository release v1.0 (github.com/arithmon/Sieve, tag `v1.0`, commit `247e165`, 2026-06-12; archived as the frozen-inputs record of the Arithmon methodology paper, Zenodo 10.5281/zenodo.20666879), and the Lean 4 formal-identity ledger is the core release v3.4.29 (github.com/gift-framework/core, commit `667c8b9`). The four-null battery and the budget-uniqueness ranking run against these frozen inputs; relations examined after the freeze are labelled exploratory.

### Traditional sensitivity diagnostics

The remainder of §7 lists the sensitivity analyses inherited from 3.4: SVD structure, effective rank, cross-correlations, PSLQ residual. These continue to serve as internal-consistency checks; they no longer stand as the statistical headline. The coincidence-probability calculations from 3.4 §7.5 are archived in Supplement S4 (Sieve diagnostics).

### 7.1 Formula Structure Analysis

**Question**: Do more complex formulas systematically produce smaller deviations? If so, the framework might be "fitting" by formula complexity.

**Method**: Random Forest regression with leave-one-out cross-validation, predicting |deviation| from formula features (number of constants used, maximum constant value, arithmetic operations, sector membership).

**Result**: R² = −0.518 (LOO-CV), worse than a mean predictor. The top feature importance is max_constant_value (0.41), followed by n_expressions (0.13) and sector_PMNS (0.12). Formula complexity does **not** predict deviation magnitude.

**Verdict**: This is inconsistent with systematic cherry-picking. Complex formulas do not perform better than simple ones; the precision is distributed uniformly across formula types.

### 7.2 Topological Constant Sensitivity

**Method**: Perturb each of the 20 structural constants by ±1 → 20 × 33 sensitivity matrix S_ij = ∂(observable_i)/∂(constant_j).

**SVD analysis**: 19 significant singular values (1 zero, corresponding to a redundant constant combination). The singular value spectrum decays smoothly: [3.50, 2.54, 2.39, 2.14, 2.11, 2.01, 1.93, 1.76, 1.68, 1.63, ...].

**Sensitivity-deviation correlation**: ρ = −0.083, **no systematic pattern**. Observables with high sensitivity to constant perturbations do not have larger deviations. The most sensitive constants are dim(K₇) (0.068), dim(G₂) (0.062), and H* (0.062).

**NK ball rigidity**: The analytical NK ball radius is δg/g ≤ 4.86 × 10⁻⁶ (the certified analytic bound). A tighter numerical eigenvalue-variation estimate gives δg/g ≤ 1.35 × 10⁻⁷ (not used as the certified analytic bound, but consistent with rigidity). The coefficient of variation of metric eigenvalues is ~10⁻⁷, confirming rigidity.

### 7.3 Effective Degrees of Freedom

**Method**: SVD of the 20 × 33 constant-usage Jacobian (binary: 1 if constant appears in formula, 0 otherwise). The effective rank is defined by the singular value decay profile.

**Results**:
- **r_eff = 15.53** effective parameters (out of 20 structural constants)
- 12 singular values capture 90% of variance
- 17 singular values capture 99% of variance
- 1 zero singular value (exact linear dependence)

**Overdetermination ratio**: 33 observables / 15.53 effective parameters = **2.13×**. The system has more than twice as many constraints as degrees of freedom: a hallmark of an overdetermined (not fitted) system.

**Most-used constants**: b₂ appears in 9 observables, dim(G₂) in 7, dim(E₈) in 7, Weyl in 6, H* in 6. No constant is used in more than 27% of observables, confirming broad coverage rather than concentration.

The constant-usage diagnostics (SVD spectrum, effective rank, per-constant usage frequency) are provided as a repository figure (`fig_constant_usage.png`).

### 7.4 Cross-Correlations

**Jaccard similarity**: Of the C(33,2) = 528 observable pairs, **155 share at least one structural constant** (29.4% coupled). All 33 observables belong to a single connected component: no prediction is isolated.

**Strong correlations**: 51 pairs have |ρ| ≥ 0.5. Notable examples:
- m_b/m_t vs m_μ/m_τ: ρ = −0.816 (both involve 2b₂ = 42)
- α⁻¹ vs Ω_DM/Ω_b: ρ = −0.721 (both involve rank(E₈))
- sin²θ_W vs Q_Koide: ρ = +0.673 (both involve b₂ and dim(G₂))

**Mean sector Jaccard** = 0.293: sectors share ~29% of structural constants, creating a web of inter-sector constraints.

The constant-observable sensitivity heatmap (cross-sector coupling matrix) is provided as a repository figure (`fig_sensitivity_heatmap.png`).

The pairwise Pearson correlations between observable deviations are provided as a repository figure (`fig_observable_correlations.png`).

### 7.5 Coincidence Diagnostics (archived to Supplement S4)

The 3.4 headline coincidence-probability figures ($\sim 10^{-346}$ under the uniform null, $\sim 10^{-133}$ under the algebraic null on 4.2M random formulas from the same 20 constants) are preserved as internal-consistency diagnostics in **Supplement S4 (Sieve diagnostics)**. They do not stand as the 3.5 statistical headline: as noted at the top of §7, coincidence-probability tails do not distinguish a genuine survivor from a well-tuned formula grammar, whereas the Sieve reading isolates the *pre-registered rank-1 budget-unique* survivors and separately certifies them via Lean 4 R2. The archived tables and null distributions remain available for readers who want the previous framing.

### 7.6 PSLQ Residual Analysis (PSLQ_residual)

Beyond the statistical null models, we apply PSLQ integer relation detection to the relative residuals r_i = (pred_i - exp_i)/exp_i of the 33 Type I observables with experimental comparison. The goal is to identify whether deviations have structural content, i.e., whether residuals match algebraic expressions built from the same topological constants.

**Method**: For each observable, we compute the residual r_i and test it against:
(1) rational approximations p/q with small denominators,
(2) structural fractions involving the framework's constants (dim(G₂), dim(E₈), b₂, b₃, H*, PSL(2,7), etc.),
(3) PSLQ integer relations with the 20 structural constants,
(4) mpmath.identify() for closed-form recognition.

**Key findings**:

| Observable | Residual r | Best match | Match error |
|------------|-----------|------------|-------------|
| δ_CP | +0.1130 | dim(G₂)/(dim(E₈)/2) = 14/124 | 0.08% |
| m_b/m_t | -0.0079 | -1/(N_gen × 2b₂) = -1/126 | exact |
| Y_p | +0.00368 | 1/(φ × PSL₂₇) = 1/(φ × 168) | 0.04% |
| sin²θ₁₂(CKM) | +0.00372 | 6/(b₂ × b₃) = 6/1617 | 0.2% |

**Documented structural observation**: Only δ_CP. A possible compactification factor (62/69 = dim(E₈) / (dim(E₈) + 4 dim(K₇))) is recorded as a post-hoc observation in S2 Appendix F; it is not adopted as a revision and the canonical prediction remains 197°.

**Not adopted**: The m_b/m_t correction (-1/126) and Y_p correction (1/(φ × 168)) are documented for future work pending structural derivation from the compactification geometry. The CKM correction (6/(b₂ × b₃)) is suggestive but at lower precision.

**Epistemological note**: The framework evolved iteratively from 6 free parameters (v1) through 4, 3, and finally 0 (v3). Each refinement added structural content while removing degrees of freedom. The δ_CP compactification factor is documented because it has geometric meaning independent of the experimental data it happens to match, but is not adopted as a revision: the raw topological 197° stands as the prediction.

---

## 8. Formal Verification and Statistical Analysis

### 8.1 Lean 4 Verification

The K₇ framework is formally verified in Lean 4 [28] with Mathlib [29]:

| Category | Count |
|----------|-------|
| Source files | 146 under `GIFT/` (140 core + 6 generated) |
| Build jobs | 8394 |
| Unproven (sorry) | 0 |
| Classified axioms | 15 in the A-F taxonomy, of which 4 external data packages (see "Axiom accounting" below) |
| Certificate conjuncts | **213** |

The conjuncts (counting top-level conjunctions of each master certificate) cover metric/torsion/topology (39), couplings/masses/mixing (56), KK spectrum (45), metric eigenvalues (15), spectral invariants (10), δ_CP compactification (6), gauge breaking (10), bundle universality (11), instanton hierarchy (14), and 7D Weyl law (7). Full per-file breakdown in Supplement S3, §3.7.

```lean
theorem weinberg_relation :
  b2 * 13 = 3 * (b3 + dim_G2) := by native_decide

theorem koide_relation :
  dim_G2 * 3 = b2 * 2 := by native_decide
```

The E₈ root system is fully proven (12/12 theorems, basis generation). G₂ differential geometry (exterior algebra, Hodge star, torsion-free condition) is axiom-free. G₂ group structure (`g2_mul_closed`, `g2_subset_SO7`, `g2_det_mul_gram`) is proven by `native_decide` (v3.4.5).

**Axiom accounting** (Lean core v3.4.29, `#print axioms` over the ledger). The core exposes **15 classified axioms** across an A–F taxonomy (**A** definitional, **B** standard mathlib-adjacent, **C** geometric, **D** literature, **E** GIFT-claims, **F** numerical) plus the standard Mathlib foundations (`propext`, `Classical.choice`, `Quot.sound`). Of the 15, **4 are the *external data-package* axioms** listed at 3.4 headline: `K7_analysis_data` (HarmonicForms.lean), `K7_spectral_data` (SpectralTheory.lean), `literature_package` (LiteratureAxioms.lean), and `KK_YM_EFT` (KKSpectralBridge.lean). The remaining 11 sit at classes A–D + F (definitional / standard / geometric / literature / numerical). Neither class E (GIFT-claims) nor a sorry is invoked. Both counts are correct and refer to different slices of the same ledger; §8 states both to avoid the discrepancy between the 3.4 abstract ("4 axioms") and the v3.4.29 release audit ("15 classified axioms").

Three certified results anchor the formalization: (1) **G₂ three-form** (Bryant-Joyce φ₀ on ℝ⁷ formalized in `G2ThreeForm.lean`, all 7 nonzero coefficients certified by `native_decide`, dim(g₂) = 14 proven); (2) **ν̄(K₇, g) = 0** CGN invariant certified zero via rectangular TCS (k₊ = k₋ = 1, CGN Main Corollary [16]); (3) **KK spectral bridge**, 4D mass gap formally conditional on KK_YM_EFT alone, all spectral ingredients Lean-certified.

### 8.2 Observable Coverage

Of the 95 observables, **55 are Lean-certified**:

| Type | Certified | Total | Coverage |
|------|-----------|-------|----------|
| **I** | 33 | 33 | 100% |
| **II** | 0 | 19 | 0% |
| **III** | 14 | 21 | 67% |
| **IV** | 8 | 22 | 36% |
| **Total** | **55** | **95** | **58%** |

Type II observables are Type I ratios × experimental VEVs (e.g., m_u = m_u/m_d × m_d(PDG)). The algebraic step (the ratio) is Lean-certified for all 33 core Type I formulas; only the physical scale identification step is uncertified. Axiomatizing VEV inputs would be circular (they are experimental inputs, not predictions). Type III coverage includes the new gauge (10+11 conjuncts) and instanton (14 conjuncts) certificates.

### 8.3 Scope of Verification

**What is proven**: Arithmetic identities relating topological integers. Given b₂ = 21, b₃ = 77, dim(G₂) = 14, etc., the numerical relations are machine-verified.

**What is not proven**: Existence of K₇, physical interpretation of ratios as SM parameters, uniqueness of formula assignments. The verification establishes **internal consistency**, not physical truth.

### 8.4 Statistical Uniqueness (Sieve baseline)

Among 192,349 alternative (gauge group, holonomy, Betti) configurations tested by Monte Carlo, zero outperform the framework: mean deviation 0.73% vs 32.9% for alternatives (P < 5 × 10⁻⁶, > 4.5σ). E₈×E₈ beats the next-best tested gauge product by ~12×; G₂ holonomy beats SU(3) (Calabi-Yau) by ~6×. Only rank 8 gives N_gen = 3 exactly. Under the Sieve reading of §7 this baseline provides the uniform-null passage figure; the algebraic-null coincidence probabilities from the 3.4 abstract ($\sim 10^{-346}$ uniform, $\sim 10^{-133}$ algebraic) are archived as diagnostics in Supplement S4. Full gauge-group and holonomy rankings in Supplement S2, §23.

### 8.5 The G₂ Metric

The predictions in §4 depend only on topological invariants. However, the G₂ metric constrained by det(g) = 65/32 is numerically constructed as a Chebyshev-Cholesky expansion with 169 parameters (companion paper [A]):

| Quantity | Value |
|----------|-------|
| ‖T‖_C⁰ (torsion) | 2.949 × 10⁻⁵ |
| NK parameter h (analytical, [A]) | 8.95 × 10⁻⁹ (β = 0.321 exact, margin ×56M) |
| NK parameter h (numerical, tighter) | 1.43 × 10⁻⁹ (β = 0.0296 numerical, margin ×350M) |
| δg/g (analytical ball radius) | ≤ 4.86 × 10⁻⁶ |
| det(g) | 65/32 (exact, by construction) |

The interval-arithmetic Newton-Kantorovich certification [A] establishes existence and uniqueness of a torsion-free G₂ metric **at the seam / datum level**, within the certified NK ball around the approximant. It is *not* a construction of a smooth compact $K_7$ carrying a globally torsion-free $G_2$-structure; see §9 for the two-layer boundary between the seam-level NK certificate here and the global-analytic layer discharged conditionally in [E]. The analytical certificate uses β = 0.321 derived from det(g) = 65/32 (exact constraint); a tighter numerical estimate β = 0.0296 yields h = 1.43 × 10⁻⁹ (margin ×350M) but rests on numerical Lipschitz estimates rather than the analytical bound. Both values certify h < 0.5 (Kantorovich threshold) by a factor of at least 56 million. Torsion reduction ×2995 in 5 Joyce iterations is well within Joyce's perturbative regime.

**Closed-form Calabi–Yau residual on the K3 fibre (interval-rigorous).** Independently of the seam-level NK certificate, the Ricci-flat Kähler metric on the Z₂³-equivariant K3 fibre admits an explicit finitely-parametrised closed-form approximation whose order-3 residual (667 parameters) satisfies a machine-checked theorem Var(log R) ≤ ε₃ = 1309/10⁷ ≈ 1.309 × 10⁻⁴ < 10⁻³ on a frozen reproducible 4000-point witness: the K3 points are Krawczyk–Rump-certified and the closed-form metric is evaluated in forward interval arithmetic, so the bound carries no floating-point assumption (the earlier conservative per-point bracket is eliminated). This is formalized in Lean 4 (`K3ClosedFormWitness`, `native_decide`, zero `sorry`, zero added axiom) and distributed in the public core package. *Honest scope:* the theorem bounds the residual over the frozen sample only; a bound over the entire K3 surface remains an open constraint-aware global-positivity (Positivstellensatz/SOS) problem. See §9 for the two-layer boundary of the compact G₂ existence question.

**Companion analytic paper [E].** The datum-level analytic scheme for a $K_7 \to S^3$ construction with $N = 77$ round-unlink branch components is discharged conditionally in [E] (concept DOI 10.5281/zenodo.21215956): seven theorem-grade results (unconditional weighted edge Fredholm theory, conditional branched lifting scheme, structural source-side discharge constants) and an outward-rounded interval-arithmetic threshold $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$. Two structural hypotheses remain OPEN at the global-analytic layer, the smooth-topology slot and the hyperkähler-K3-fibration-with-periods slot (jointly the pack $\mathsf H_{\rm global}$), plus the standalone anisotropic three-scale perturbation theorem $(\mathrm J)$ for the torsion-free upgrade. §9 records the two-layer boundary explicitly.

**Analytic invariant.** The Crowley-Goette-Nordström invariant ν̄(K₇, g) ∈ ℤ [16] vanishes for any rectangular TCS with twisting numbers k₊ = k₋ = 1, by CGN Main Corollary (gluing angle θ = π/2 forced). This conditional result is certified in Lean (`TCSConstruction.lean`, `K7_nu_bar_zero`): if K₇ is realized as a rectangular TCS, then ν̄(K₇, g) = 0 without any additional geometric computation. The rectangular-TCS building-block identification is superseded at the datum layer by the [E] $K_7 \to S^3$ scheme; the CGN certificate is retained as an independent structural constraint.

**Results from the NK-certified metric.** Analysis of the NK-certified metric yields several structural results:

- **V_min formula**: The minimum associative cycle volume satisfies V_min = √(Vol(K₇)/11), where 11 = b₃/n = 77/7. NK numerical value 219.90; formula gives 221.24 (0.6% agreement). (Volumes here use the global NK normalization of [A]; the per-cycle associative volumes of §6.3 use the fibre-level normalization and are not directly comparable.)
- **Harmonic decompositions**: b₂ = 7+14 = 3+18 = 11+10 (G₂ reps / hyperkähler triple / TCS blocks); b₃ = 35+42 = (1+7+27)+2×21; spectral gap 10522× between zero and non-zero modes.
- **U(1)² symmetry (numerical)**: Period integrals $S_\theta = S_\psi = 6.1265$, numerically compatible with a U(1)² isometry to tolerance $2.6 \times 10^{-8}$; the pattern propagates from the metric to all period integrals within the same tolerance.
- **Universal law**: λ₁ × H* = 12.3364 holds for all 66 known G₂ manifolds (the 65 catalogued literature examples of §2.3 plus K₇ itself).
- **Lepton hierarchy from periods**: ln(m_τ/m_e) = 8.154 from SD associative volumes → e^8.154 = 3477 ✓ (Lean-certified).

Full details in the companion paper [A].

---

## 9. Existence status (two-layer boundary)

This section gathers in one place every claim about *existence of the compact $K_7$ manifold* (in the 3.4 version this language was distributed across §2.3, §5.5, §8.5, §9.3 and §9.4).

### 9.1 The two layers

- **Datum-level analytic layer.** A *normalised datum* $\mathcal D_0$ is fixed: a base $S^3$, a coassociative fibration $\pi: K_7 \to S^3$ with $N = 77$ round-unlink branch components, a rank-one Picard–Lefschetz monodromy along each meridian, and a collar-affine period section $h$. On $\mathcal D_0$, the companion paper [E] discharges a *conditional branched lifting scheme* (a special case of [Donaldson 2017, Conjecture 1]) via seven theorem-grade results. Source-side discharge constants ($C_{\rm src} = 27/16$, $C_{\rm nl} = 2/3$, $C_{\rm link} \le 1$, $\gamma = 2$) carry structural provenance; the tail-contraction threshold $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$ is enclosed by outward-rounded interval arithmetic. This layer is *discharged conditionally*.
- **Global-analytic layer.** Passing from $\mathcal D_0$ to a smooth compact $K_7$ carrying a genuine torsion-free $G_2$-structure requires closing two external structural slots (jointly the pack $\mathsf H_{\rm global}$ of [E]) and one standalone hypothesis:
  1. **(L4) Hyperkähler K3-fibration with prescribed periods.** A hyperkähler K3-fibration on the topological $K_7$ with rank-one Picard–Lefschetz monodromy along the 77-unlink and a period section $h$ consistent with the datum $(A_2^{\rm collar})$–$(A_5)$. The closest current template is Kovalev's Hopf-link coassociative K3 fibration [13].
  2. **(C.1) Smooth-topology classification.** Standard Voisin–Kovalev classification data [21] (Dehn–Seidel twist + K3-family gluing) producing the smooth 7-manifold on which (L4) can live.
  3. **(J) Anisotropic three-scale perturbation.** An anisotropic refinement of Joyce's perturbation theorem in the three-scale collapsing geometry (base $\sim 1$, fibre $\sim \varepsilon$, EH core $\sim \varepsilon^2$) that deforms each closed positive $3$-form of $G_2$-type to a torsion-free $G_2$-structure. Companion long chantier of [E].

### 9.2 Chain-of-locks and certified constants (headline synthesis from [E])

To make this section standalone (a reader interested only in the framework's headlines should not have to open [E] to gauge the status of the analytic layer), we transport here the four-lock status ledger and the source-side certified constants. Preambles, proofs, and the interval-arithmetic verification scripts remain in [E].

**Datum $\mathcal D_0$.** Base $S^3$; coassociative fibration $\pi: K_7 \to S^3$; branch locus $= N = 77$ round unlink components; rank-one Picard–Lefschetz monodromy along each meridian; collar-affine period section $h$ satisfying gates $(A_2^{\rm collar})$–$(A_5)$; source-side envelopes $K_{F_1} = K_{g_1} = K_v = 2.31$, $D_{\rm max} = 7/5$.

**Chain-of-locks.** The datum-level scheme is discharged conditionally on four locks:

| Lock | Content | Status | Source-side certificate |
|------|---------|--------|-------------------------|
| **L1.6** | $K_{\rm Sch} \le 16/3$ (Schauder recollement) | CLOSED (theorem-grade) | $q_{\rm coeff} \le 0.168$ headline / $0.110$ sharp, slack 33–56% of $1/4$; independent checker gates O0–O5 all pass |
| **L2 assembly** | Product-Banach reconstruction on $X_{\rm AR} = X_\omega \times X_\lambda \times X_\mu \times X_\Theta$ | CLOSED (theorem-grade) | $q_{\rm total} = q_{\rm AR} + q_{\rm comm} + q_{\rm proj} + q_{\rm hodge} + q_{\rm gauge} \le 8.20 \cdot 10^{-3} < 1/2$, contraction margin $60.9\times$; gates A0–A5 pass |
| **(L4)** | Hyperkähler K3-fibration with prescribed periods, rank-1 PL monodromy along the 77-unlink | OPEN | Closest template: Kovalev's Hopf-link coassociative K3 fibration [13] |
| **(C.1)** | Smooth-topology classification (7-manifold on which L4 lives) | OPEN | Voisin–Kovalev classification inputs [21] |
| **(J)** | Anisotropic three-scale perturbation theorem (base $\sim 1$, fibre $\sim \varepsilon$, EH core $\sim \varepsilon^2$) | OPEN | Companion long-term project of [E] |

The two CLOSED locks are certified theorem-grade at $\mathcal D_0$; the three OPEN slots constitute $\mathsf H_{\rm global} + (\mathrm J)$ and are the entire residual analytic content of the global-existence problem.

**Source-side certified constants (all with structural provenance).**

| Constant | Value | Provenance |
|----------|-------|------------|
| $C_{\rm src}$ | $27/16$ | Cubique + parity on Donaldson canonical gauge |
| $C_{\rm nl}$ | $2/3$ | Simons identity + collinéarité; strictly below power-counting envelope $\sim 1$ |
| $\gamma$ | $2$ | Proved (spectral order of the rank-1 branched Jacobi) |
| $C_{\rm link}$ | $\le 1$ | Impedance match ($\Lambda_{\rm in}^{\rm model} = -\Lambda_{\rm out}^{\rm model}$) $+$ parity $\vert c_0\vert^2$ $+$ type-IV sectional $\vert K_{\rm sec}\vert \le 1$ |
| $C_{\rm FS}$ | $\approx 0.29$ | Fredholm–Schauder outer coercivity |
| $\kappa_E$ | $\le 1.02$ | KA1 rigidity ($Df_0 = 0$, $\Vert g'\Vert_\infty = 2$, compensation $(1/2) \cdot 2$) |
| $\Vert\mathfrak M^{-1}\Vert$ | $\le 1.246$ | Schur + sector exhaustion at the gate $\Pi_{\rm obs}^{\rm PDE} = \Pi_{\mathbb R^{2N}}$ |
| $K_{\rm Sch}$ | $\le 16/3$ | L1.6 recollement, $(4/3) \cdot 3 \cdot (4/3)$ |
| $q_{\rm total}$ | $\le 8.20 \cdot 10^{-3}$ | L2 assembly (Banach product-max on $X_{\rm AR}$) |
| $R_0(\mathcal D_0)$ | $\le 4.9 \cdot 10^3$ | Newton–Kantorovich tail-contraction threshold; outward-rounded interval arithmetic |

**Operational meaning of $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$.** In the branched adiabatic scheme, $R_0$ is the smallest scale beyond which the Newton–Kantorovich iteration converges uniformly at $\mathcal D_0$. The certified bound $R_0 \le 4.9 \cdot 10^3$ is enclosed by outward-rounded interval arithmetic (mpmath, 60 decimal digits, 4 verification scripts distributed with [E]) and sits $\sim 9\%$ below the power-counting threshold and $\sim 34\%$ below the Codex Stage B recomputation; i.e., the datum-level scheme is discharged *with headroom*, not marginally.

*What this section does not claim.* All quantities above are proved conditional on $\mathcal D_0$ being an admissible datum; passing to a smooth compact $K_7$ requires (L4), (C.1), (J) discharged (see §9.3). No proof is reproduced here; [E] is the citable source.

### 9.3 What the 3.5 headline claims and what it does not

The framework's headline claims (95 observables, Sieve reading, Lean 4 formal ledger) are **datum-conditional**: they do not consume $\mathsf H_{\rm global}$ or $(\mathrm J)$. They rest only on the topological invariants of $K_7$, the E₈×E₈ gauge structure, and the algebraic constraints of §2, all of which are certified at the datum layer.

The framework does not claim to have constructed a smooth compact $K_7$ carrying a globally torsion-free $G_2$-structure. That claim would require (L4), (C.1), (J) discharged as theorems; all three are currently OPEN, honestly flagged in [E, Remark 1.1′] and in the CHANGELOG.

### 9.4 Relation to earlier drafts and independent constructions

The Chebyshev–Cholesky NK-certified metric of [A] is a *seam-level* result: it certifies the seam that would appear in a TCS-style gluing, subject to that gluing being available. The [D] paper strengthens the neck-level ansatz with an explicit closed-form + five-layer Wirtinger certificate. The [E] paper replaces the TCS-style rectangular hypothesis with a rank-one branched Kovalev–Lefschetz hypothesis; the corresponding CGN invariant statement of §8.5 is retained as an independent structural constraint but no longer serves as the load-bearing existence argument. The (b₂, b₃) = (21, 77) outlier position in the catalogued literature is a documentary observation, not an obstruction, and is consistent with the datum-level scheme of [E].

---

## 10. Selection audit: five investigated routes and one open residual

Earlier versions carried three separate "why (21, 77)?" paragraphs (in the 3.4 numbering: §2.3, §3.4, §9.4), all treating the selection question as OPEN. This section collects them into a single statement of **an audit** across five investigated routes: each route is closed as an investigation, but the underlying selection question itself is *not* declared closed. A single explicit residual (§10.2) is the honest gate.

### 10.1 Five investigated routes

- **Q-B landscape (main).** In the Kondō–Shimada (KS-full) census of Picard-lattice-realisable hyperkähler K3 pairs, $(b_2, b_3) = (21, 77)$ is *sub-median*: convenient, not special. It does not survive as a distinguished point of a natural density.
- **Sub-Q1 (σ_B chamber structure).** The associated $\sigma_B$ Piroddi–style analysis yields one chamber with three gauge orbits off 77; the "77 side" is not a distinguished chamber wall.
- **Sub-Q2 (S₅ propagation).** The propagation $S_5 \to 77$ closes negatively by mechanism: the induced flow on the moduli side does not favour 77 over its neighbours.
- **Sub-Q3 (higher-HK distinguishedness).** Across multiple higher-HK classifications, $(21, 77)$ is *distinguishedness-yes-77-no*: the number 21 has invariants making it stand out (Petersen, S₅-related, several classification schemes), but 77 does not follow the same distinguishedness under any of the *presently specified mechanisms*. Sole exception: Camere Lemma 5.8, which distinguishes 77 in a single, structurally isolated context. The arbiter table is kept in a private working note (`selection_principle_master_2026_07_04`, available on request).
- **Zhou series (independent construction).** The Zhou 2023–2025 series produces $(b_2, b_3) = 7 \cdot (3, 11)$ as a *reparametrisation* rather than an alternative construction; a screening constraint $b_2 \ge 9$ is refuted by the KS-full census (b₂ = 3 is present, so the constraint is not universal), and Zhou himself retrograded the associated 05/2026 preprint.

### 10.2 The honest residual

The sole surviving structural observation is the identity
$$
b_2 + b_3 \;=\; 98 \;=\; 7 \cdot 14 \;=\; \dim(K_7) \cdot \dim(G_2).
$$
This identity is currently held as a *pre-registration target*, not a theorem: until a mechanism is exhibited that forces $b_2 + b_3$ to equal $\dim(K_7) \cdot \dim(G_2)$ on structural grounds (not numerological ones), the identity remains a numerological observation to be flagged as such. Deriving it would dissolve the selection question into the framework's single underdetermined input: *why $(7, 3, 8) = (\dim K_7, N_{\rm gen}, \mathrm{rank}(E_8))$?*, which is of the same status as *why Lorentz invariance?* or *why the Dirac equation?*: the minimal consistent mathematical structure compatible with observed symmetries.

Pending that derivation, the operational status of $(b_2, b_3) = (21, 77)$ in this paper is that of a **construction-level input**, on the same footing as the G₂-holonomy postulate itself, rather than a derived quantity. It is not, however, an input *chosen from a continuum*: it is a discrete pair, fixed before the 33-relation ledger was frozen (§7), consumed rigidly by every relation of §4, and tied to $N_{\rm gen} = 3$ by the generation constraint of §2.4.2. The Sieve reading of §7 is calibrated to exactly this situation: it does not ask whether the input is derivable, but whether, *given* the input, the derived relations are distinguishable from grammar-tuned coincidences.

### 10.3 Consequences for §3–§7

The audit having closed each of the five routes (while the selection question itself remains open pending the §10.2 residual) tightens two earlier passages: (i) the 3.4 language "the formula selection principle is understood" / "remains open" is superseded by pointer to this section; (ii) the Sieve reading of §7 is the correct probabilistic language for a framework whose selection question has been *audited* at the level of *mechanism* (survivors of a public battery + pre-registered rank + Lean formal-identity flag), rather than the earlier coincidence-probability language.

---

## 11. Discussion, Falsifiability, and Conclusion

### 11.1 Falsifiable Predictions

![**Figure 3.** Falsifiability roadmap for $\delta_{CP}$. The K₇ framework predicts $\delta_{CP} = 197°$ (pure geometry, $= 7 \times 14 + 99$). The NuFIT 6.1 $1\sigma$ band ${207^{+23}_{-20}}°$ (NO w/o SK-atm) contains the prediction. The consistency window $[182, 212]°$ (± $\sim 15°$ DUNE-representative resolution) around the prediction is the target measurement zone; a $3\sigma$ measurement outside this window would create serious tension with the framework. Timeline (below): first beam targeted 2031 per Fermilab's 2026 milestone; few-degree $\delta_{CP}$ resolution expected mid-2030s, $\sim 15°$ final resolution by late 2030s / 2040s.](fig_delta_cp_roadmap.pdf)

The framework makes concrete, testable predictions. The most critical:

**δ_CP**: The topological prediction is δ_CP = 197° = 7 × 14 + 99 (pure geometry). The experimental central value has drifted between recent NuFIT releases: 197° (5.2, 2022, exact match), 177° ± 20° (6.0, Oct 2024, 11.3% off centre / edge of 1σ), $207^{+23}_{-20}$° (6.1, 2025, NO w/o SK-atm; 197° now well inside the 1σ band). Further drift with T2K / NOvA / DUNE statistics is expected. A possible compactification factor (62/69) is recorded as a post-hoc observation in S2 Appendix F and is not part of the main prediction. **Falsification**: δ_CP outside [182, 212]° at 3σ creates serious tension.

**N_gen = 3**: No flexibility. A fourth generation immediately falsifies.

**m_s/m_d = 20**: Lattice QCD target precision ±0.5 by 2030. Current 20.0 ± 1.0.

**sin²θ_W = 3/13**: FCC-ee precision ~10⁻⁵ (4× improvement over current).

**New tests from Types II/III**:
- Combined lepton ratios (wilson_line+instanton): τ/e = 3485, τ/μ = 16.69, μ/e = 208.8, all within 1%
- Proton lifetime: τ_p = 4.06 × 10³⁸ years (beyond near-term Hyper-K sensitivity; Hyper-K constrains nearby GUT alternatives)
- SUSY spectrum: m_gravitino = 166 GeV, m_moduli = 3.2 TeV (model-dependent; viable only in compressed/suppressed-coupling realizations, requiring dedicated recast at HL-LHC)

| Experiment | Observable | Timeline | Test Level |
|------------|------------|----------|------------|
| DUNE Phase I | δ_CP (3σ) | after first beam (targeted 2031 [47]); exposure-dependent | Critical |
| DUNE Phase II | δ_CP (5σ) | late 2030s–2040s | Definitive |
| Lattice QCD | m_s/m_d | 2028–2030 | Strong |
| HL-LHC | SUSY spectrum | 2029–2040 | Complementary |
| Hyper-Kamiokande | δ_CP, τ_p | 2034+ | Complementary |
| FCC-ee | sin²θ_W, Q_Koide | 2040s | Definitive |

### 11.2 Relation to M-Theory

E₈×E₈ and G₂ holonomy connect directly to M-theory [24, 33, 34]:
- Heterotic string theory requires E₈×E₈ for anomaly cancellation [19]
- M-theory on G₂ manifolds preserves N = 1 SUSY in 4D [24, 34]

The K₇ framework differs from standard M-theory phenomenology [35] by focusing on topological invariants rather than moduli stabilization. Where M-theory faces the landscape problem (~10⁵⁰⁰ vacua), the framework proposes that topological data alone constrain the physics. The G₂-MSSM spectrum (§5.3, §6.5) is consistent with the phenomenology of Acharya et al. [33]: m_gravitino ~ 100 GeV, m_moduli ~ few TeV, gaugino condensation in a hidden sector.

The second E₈ factor is required by anomaly cancellation but has no direct coupling to Standard Model fields. In heterotic M-theory [32], gaugino condensation in this hidden sector drives SUSY breaking [33]. A natural interpretation is that the hidden E₈ sector provides the dark sector: the predicted cosmological ratios Ω_DM/Ω_b = 43/8 and Ω_DE = ln(2) × 98/99 emerge from the same topological invariants (b₂, b₃) that also determine dim(E₈) = 248. Whether this connection is structural or coincidental remains an open question.

### 11.3 Comparison

| Criterion | K₇ framework | String Landscape | Lisi E₈ |
|-----------|------|------------------|---------|
| Falsifiable predictions | Yes (δ_CP = 197°, N_gen) | Not yet (landscape selection) | Not yet (embedding obstruction) |
| Adjustable parameters | 0 | ~10⁵⁰⁰ | 0 |
| Formal verification | 213 conjuncts, 15 axioms (A-F taxonomy, of which 4 external data packages) | No | No |
| Precise predictions | 95 (66 testable) | Qualitative | ~5 |
| Gauge breaking | E₈→SM (6 levels) | Landscape-dependent | Single E₈ |
| Mass hierarchy | 2 mechanisms + combined | n/a | n/a |
| Statistical distinctiveness | Sieve reading (§7): unique rank-1 survivor m_H/m_W; r_eff = 15.53; archived coincidence-probability diagnostics ~10⁻³⁴⁶ uniform / ~10⁻¹³³ algebraic in Supplement S4 | n/a | n/a |

**Distler-Garibaldi obstruction** [36]: Lisi's E₈ ToE attempted to embed all Standard Model fermions (including chirality) directly as roots of a single E₈. Distler and Garibaldi proved this is mathematically impossible: no E₈ representation decomposes into the correct chiral SM spectrum. The K₇ framework avoids this obstruction entirely: E₈×E₈ is the gauge group of heterotic M-theory (not a particle container), the SM gauge group emerges from a six-level breaking chain (§5.1), fermion generations are fixed by the topological constraint (rank(E₈) + N_gen)b₂ = N_gen·b₃, giving N_gen = 3 (§2), and chirality is a topological property of the compactification. The mathematical relationship between E₈ and particles is cohomological (emergence from geometry), not representational (embedding into algebra).

### 11.4 Limitations and Open Questions

*The main open questions have been promoted to their own sections: existence status (§9) and selection principle (§10). This subsection lists residual limitations that do not fit into either.*


| Issue | Status | Section |
|-------|--------|---------|
| K₇ topological classification | Certified metric [A]; (b₂,b₃)=(21,77) absent from all catalogued compact G₂ constructions; orthogonal TCS excluded by parity; datum-level analytic scheme discharged conditionally in [E]; two structural slots (H_global) + (J) remain OPEN, see §9 (two-layer boundary) | §2, §9 |
| Singularity structure | Not yet supplied: the smooth-K₇ framework encodes the SM gauge group at the algebraic/spectral level (§5), but a physical chiral non-abelian realization requires an explicit UV bridge (Acharya-Witten singular locus, heterotic bundle, or duality identification); load-bearing open question | §2.3, §5.1, §9 |
| Formula selection rules | Consequence of metric uniqueness; five selection routes audited and closed (§10), the selection question itself open, with sole residual b₂+b₃ = 98 = dim(K₇)·dim(G₂) flagged as a pre-registration target | §3, §7, §10 |
| α_s(RGE) = 3.7% deviation | SUSY threshold sensitivity; split-spectrum matching (§5.3) | §5 |
| δ_CP deviation | 4.8% vs NuFIT 6.1 central (207 +23/−20 °), inside the 1σ band; 11.3% vs the earlier NuFIT 6.0 central; see S2 Appendix F | §4.2, §11.1 |
| ΔV(e-μ) = 15.9% deviation | Reduced to 0.75% by combined pipeline (α = e^K) | §6 |
| Hidden E₈ sector | Candidate for dark sector (dark matter, dark energy); no direct observable coupling | §11.2 |
| Quantum gravity completion | Not addressed | n/a |

**Note on the selection audit.** A natural objection to any framework of this type is: why *these* observables, and why the specific pair $(b_2, b_3) = (21, 77)$? This paper treats this question as a completed *audit* (five routes investigated and closed, the selection question itself remaining open) with a single honest residual, and gives it its own section (§10). In short: (i) at the *Q-B landscape* level, $(21, 77)$ is sub-median across the KS-full census: convenient, not special; (ii) at the *sub-Q1* level, $\sigma_B$ has one chamber with three gauge orbits off 77; (iii) at the *sub-Q2* level, propagation $S_5 \to 77$ closes negatively by mechanism; (iv) at the *sub-Q3* level, higher-HK distinguishedness is *yes-77-no* across multiple classifications, with a single exception at Camere Lemma 5.8; (v) the Zhou series produces $(b_2, b_3) = 7\cdot(3,11)$ as a reparametrisation, and the screening constraint $b_2 \ge 9$ is refuted by census. The residual is the derivation of the identity $b_2 + b_3 = 98 = \dim(K_7) \cdot \dim(G_2)$, which is currently held as a pre-registered numerology target rather than a theorem. Beyond this residual, the selection question dissolves into "why $(7, 3, 8)$?", the framework's single underdetermined input. Extended grammar analysis is consistent with the audit reading: adding TCS-level atoms (χ(K3)=24, b₂(M₁)=11, b₂(M₂)=10) to the search grammar discovers simpler and more precise formulas for several observables (m_c/m_s exact at 11+7/10; Ω_DE = 53/77 rational, 5× more precise against Planck 2018 though 2.5× less precise against the frozen PR4 dataset (S2 §18); λ_H = 11/87, 7× more precise), suggesting that observables encode TCS structure more directly than the higher-level algebra currently used.

At the philosophical level, the residual question is: why G₂ holonomy? This is the framework's single underdetermined input: the assertion that the compactification manifold is a compact 7-manifold of G₂ holonomy with b₁=0. This question has the same epistemic status as "why Lorentz invariance?" in special relativity or "why the Dirac equation?" in relativistic quantum mechanics: it is recognized as the minimal consistent mathematical structure compatible with observed symmetries, not derived from something more fundamental. G₂ is the unique compact exceptional holonomy group admitting a 7-dimensional representation with the stabilizer chain G₂ ⊃ SU(3) ⊃ SU(2) ⊃ U(1) that naturally contains the subgroups needed for the Standard Model. It is not selected among a continuum of options; it is an isolated point in the space of compatible mathematical structures. Every other unification framework carries equivalent or greater philosophical inputs alongside substantially more free parameters. The K₇ framework makes its single input explicit.

**Honest assessment of outliers**: Two observables deviate by >5%: ΔV(e-μ) (15.9%, reduced to 0.75% by the combined wilson_line+instanton pipeline) and κ(gauge) (4.7%, reflecting genuine K₇ geometry at the percent level). The α_s(RGE) deviation is 3.7% with split-spectrum matching (§5.3); a naive degenerate-spectrum treatment would give 12.1%. The δ_CP prediction (197°) sits inside the NuFIT 6.1 1σ band ($207^{+23}_{-20}$°, deviation 4.8% from centre); relative to the earlier NuFIT 6.0 release (177° ± 20°) it sat at the edge of the band. A post-hoc structural observation involving a 62/69 factor is recorded in S2 Appendix F (not part of the main framework).

### 11.5 Conclusion

We have explored a framework of **conditional topological relations** for Standard-Model parameters, derived from the topology of a compact $G_2$-holonomy manifold $K_7$ and an E₈×E₈ gauge structure. Following the four-layer status statement at the start of the paper:

- **Formally established.** 33 Type I algebraic identities (no continuously adjustable parameter after ledger freeze); Lean 4 layer with 213 conjuncts, 15 classified axioms in the A–F taxonomy (of which 4 external data-package axioms), 0 sorry.
- **Established at the datum.** NK-certified seam metric [A], JK topological/lattice route to $(21, 77)$, datum-level analytic scheme [E] with $R_0(\mathcal D_0) \le 4.9 \cdot 10^3$ interval-enclosed.
- **Open.** (i) A full compact torsion-free $G_2$ realisation on a smooth $K_7$: conditional on the two-slot external pack $\mathsf H_{\rm global}$ of [E] plus the anisotropic three-scale perturbation hypothesis $(\mathrm J)$; see §9. (ii) A UV bridge to non-abelian chiral gauge physics, in the sense of Acharya–Witten: not supplied here. The framework provides an algebraic/spectral dictionary; a physical mechanism (heterotic bundle, duality, or singular $G_2$ locus) remains to be given.
- **Phenomenology.** 66 testable observables; unique rank-1 budget-unique survivor $m_H/m_W$ under the Sieve reading of §7; falsifiable predictions δ_CP = 197°, $N_{\rm gen} = 3$, $\sin^2\theta_W = 3/13$, testable against NuFIT 6.1 / DUNE / FCC-ee.

**We do not claim this framework is correct.** It may represent genuine geometric insight, effective approximation, or elaborate coincidence. Only experiment can discriminate. The deeper questions (why octonionic geometry would determine particle physics parameters, and what UV mechanism bridges the smooth $K_7$ dictionary to non-abelian chiral gauge physics) remain open.

But the empirical pattern of 66 conditional topological relations that agree with experiment at sub-percent to a-few-percent level (using no continuously tuned parameter in the 33 Type I relations after ledger freeze), together with the Sieve reading isolating a certified rank-1 budget-unique survivor, and the overdetermination ratio of 2.13× on the sensitivity baseline, suggests that topology and physics may be more intimately connected than currently understood, though the mechanism of that connection is precisely what the open UV bridge above asks for.

**The ultimate arbiter is experiment.**

---

## Author's note

This framework was developed through sustained collaboration between the author and several AI systems, primarily Claude (Anthropic), with substantial contributions from GPT and Codex (OpenAI), whose independent audit passes, certificate engineering, and adversarial review rounds shaped this version and the companion paper [E], and with contributions from Gemini (Google), Grok (xAI), Kimi (Moonshot AI), DeepSeek, and GLM (Zhipu AI) for specific mathematical insights and review passes. The formal verification in Lean 4, architectural decisions, and many key derivations emerged from iterative dialogue sessions over several months; AI-derived mathematical content is either machine-verified (Lean 4, interval arithmetic) or cross-checked against primary sources. This collaboration follows a transparent crediting approach for AI-assisted research, consistent with publisher AI policies [38]. Mathematical constants underlying these relationships represent timeless logical structures that preceded human discovery. The value of any theoretical proposal depends on mathematical coherence and empirical accuracy, not origin. Mathematics is evaluated on results, not résumés.

---

## Competing Interests

The author declares no competing interests.

---

## References

*The bibliography is shared across the main paper and Supplements S1-S4; entries cited only in a supplement are marked accordingly.*

[1] Particle Data Group, Phys. Rev. D 110, 030001 (2024)
[2] S. Weinberg, Phys. Rev. D 13, 974 (1976)
[3] Planck Collaboration, A&A 641, A6 (2020)
[4] A.G. Riess et al., ApJL 934, L7 (2022)
[5] C.D. Froggatt, H.B. Nielsen, Nucl. Phys. B 147, 277 (1979)
[6] Y. Koide, Lett. Nuovo Cim. 34, 201 (1982)
[7] C. Furey, *Standard Model Physics from an Algebra?*, PhD thesis, U. Waterloo (2015)
[8] C. Furey, M.J. Hughes, Phys. Lett. B 831, 137186 (2022)
[9] R. Wilson, arXiv:2404.18938 (2024)
[10] T.P. Singh et al., arXiv:2206.06911v3 (2024)
[11] B.S. Acharya, S. Gukov, Phys. Rep. 392, 121 (2004)
[12] L. Foscolo et al., Duke Math. J. 170, 3 (2021)
[13] A. Kovalev, "Coassociative K3 fibrations of compact G₂-manifolds," arXiv:math/0511150 (2005)
[14] M. Haskins, J. Nordstrom, arXiv:1809.09083 (2022)
[15] G. Bera, "Associative submanifolds in twisted connected sum G₂-manifolds," arXiv:2209.00156 (2022)
[16] D. Crowley, S. Goette, J. Nordstrom, Invent. Math. (2025)
[17] T. Dray, C.A. Manogue, *The Geometry of the Octonions*, World Scientific (2015)
[18] J.F. Adams, *Lectures on Exceptional Lie Groups*, U. Chicago Press (1996)
[19] D.J. Gross et al., Nucl. Phys. B 256, 253 (1985)
[20] D.D. Joyce, *Compact Manifolds with Special Holonomy*, Oxford U. Press (2000)
[21] A. Kovalev, J. Reine Angew. Math. 565, 125 (2003)
[22] A. Corti et al., Duke Math. J. 164, 1971 (2015)
[23] R. Harvey, H.B. Lawson, Acta Math. 148, 47 (1982)
[24] B.S. Acharya, Class. Quant. Grav. 19, 5619 (2002)
[25] B.S. Acharya, E. Witten, arXiv:hep-th/0109152 (2001)
[26] The NOvA Collaboration, The T2K Collaboration, "Joint neutrino oscillation analysis from the T2K and NOvA experiments," Nature 646, 818-824 (2025), doi:10.1038/s41586-025-09599-3, arXiv:2510.19888
[27] I. Esteban, M.C. Gonzalez-Garcia, M. Maltoni, I. Martinez-Soler, J.P. Pinheiro, T. Schwetz, "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," JHEP 12 (2024) 216, arXiv:2410.05380; NuFIT 6.1 (2025), www.nu-fit.org
[28] L. de Moura, S. Ullrich, CADE 28, 625 (2021)
[29] mathlib Community, github.com/leanprover-community/mathlib4
[30] DUNE Collaboration, FERMILAB-TM-2696 (2020)
[31] DUNE Collaboration, arXiv:2103.04797 (2021)
[32] E. Witten, "Strong coupling expansion of Calabi-Yau compactification," Nucl. Phys. B 471, 135 (1996)
[33] B.S. Acharya et al., Phys. Rev. D 76, 126010 (2007)
[34] M. Atiyah, E. Witten, Adv. Theor. Math. Phys. 6, 1 (2002)
[35] G. Kane, *String Theory and the Real World* (2017)
[36] J. Distler, S. Garibaldi, Commun. Math. Phys. 298, 419 (2010)
[37] J.C. Baez, "Octonions and the Standard Model," math.ucr.edu/home/baez/standard/ (2020-2025)
[38] Springer Nature, "Artificial intelligence (AI) policy," www.springernature.com/gp/policies (2024)

[39] J.A. Wheeler, "Information, physics, quantum: the search for links," in *Complexity, Entropy, and the Physics of Information* (W.H. Zurek, ed.), Addison-Wesley (1990), pp. 3–28.

[40] J. Worrall, "Structural realism: the best of both worlds?," *Dialectica* **43**, 99–124 (1989).

[41] J. Ladyman, D. Ross, *Every Thing Must Go: Metaphysics Naturalized*, Oxford University Press (2007).

[42] R. Pinčák, A. Pigazzini, M. Pudlák, E. Bartoš, "Geometric origin of a stable black hole remnant from torsion in G₂-manifold geometry," Gen. Rel. Grav. **58**, 29 (2026), doi:10.1007/s10714-026-03528-z

[43] D. Joyce, S. Karigiannis, "A new construction of compact torsion-free G₂-manifolds by gluing families of Eguchi-Hanson spaces," *J. Differential Geom.* (2021), arXiv:1707.09325 (2017).

[44] S. Mukai, "Finite groups of automorphisms of K3 surfaces and the Mathieu group," *Invent. Math.* **94**, 183–221 (1988).

[45] A. Garbagnati, A. Sarti, "Symplectic automorphisms of prime order on K3 surfaces," *J. Algebra* **318**, 323–350 (2007), arXiv:math/0603742.

[46] LEP2 SUSY Working Group (ALEPH, DELPHI, L3, OPAL), "Combined LEP chargino results, up to 208 GeV for large m0," note LEPSUSYWG/01-03.1 (2001), lepsusy.web.cern.ch.

[47] Fermilab News, "Fermilab marks major milestone for world-leading DUNE experiment," 7 May 2026, news.fnal.gov.

[48] J. Chen, H. Hong, "Intermediate curvature and splitting theorem," arXiv:2604.26529 (2026). (Cited in Supplement S1, §8.4.)

[49] M. Tristram et al., "Cosmological parameters derived from the final Planck data release (PR4)," A&A 682, A37 (2024). (Frozen cosmology dataset; cited in Supplement S3.)

[50] E. Tiesinga, P.J. Mohr, D.B. Newell, B.N. Taylor, "CODATA recommended values of the fundamental physical constants: 2022," Rev. Mod. Phys. 97, 025002 (2025). (Cited in Supplement S3.)

[51] V. V. Nikulin, "Integer symmetric bilinear forms and some of their geometric applications," Izv. Akad. Nauk SSSR Ser. Mat. 43 (1979); English transl. Math. USSR Izv. 14 (1980).

---

## Author's Related Works

*Code and Lean proofs:* [github.com/Arithmon/K7-Lean](https://github.com/Arithmon/K7-Lean) (Lean module under `core/Lean/`).

- **[A]** B. de La Fournière, "An Explicit Approximate G₂ Metric on a Compact 7-Manifold with Certified Torsion-Free Completion," Zenodo [10.5281/zenodo.19892350](https://doi.org/10.5281/zenodo.19892350) (2026).
- **[B]** B. de La Fournière, "Spectral Geometry of the G₂-GIFT Manifold: Betti Numbers, KK Spectrum, and Spectral Invariants," Zenodo [10.5281/zenodo.19893371](https://doi.org/10.5281/zenodo.19893371) (2026).
- **[C]** B. de La Fournière, "Newton-Kantorovich Certificate for the K3 Donaldson Embedding in the G₂-GIFT Metric," Zenodo [10.5281/zenodo.19708916](https://doi.org/10.5281/zenodo.19708916) (2026).
- **[D]** B. de La Fournière, "An Explicit Closed-Form G₂ Ansatz on a K3-Coassociative Neck with Hyperkähler Rotation and Picard-Lefschetz Wirtinger Certificate," Zenodo [10.5281/zenodo.20039066](https://doi.org/10.5281/zenodo.20039066) (2026).
- **[E]** B. de La Fournière, "A conditional analytic framework for rank-one branched Kovalev–Lefschetz adiabatic limits," Zenodo concept DOI [10.5281/zenodo.21215956](https://doi.org/10.5281/zenodo.21215956) (2026), main paper 68 pp + supplement 25 pp, dataset and four `mpmath.iv` verification scripts (CC BY 4.0). The datum-level analytic scheme for a $K_7 \to S^3$ construction with $N = 77$ round-unlink branch components; the load-bearing external reference of the K₇ framework v3.5.
- **[essay]** B. de La Fournière, *Orientation, not ontology* (companion essay), <https://arithmon.substack.com/p/orientation-not-ontology> (2026).

*This is the founding-framework paper of the [Arithmon program](https://arithmon.com); the master methodology paper on counting coincidences and the Sieve construction is at Zenodo [10.5281/zenodo.20666879](https://doi.org/10.5281/zenodo.20666879).*

## Appendix A: Topological Input Constants

| Symbol | Definition | Value |
|--------|------------|-------|
| dim(E₈) | Lie algebra dimension | 248 |
| rank(E₈) | Cartan subalgebra dimension | 8 |
| dim(K₇) | Manifold dimension | 7 |
| b₂(K₇) | Second Betti number | 21 |
| b₃(K₇) | Third Betti number | 77 |
| dim(G₂) | Holonomy group dimension | 14 |
| dim(J₃(O)) | Jordan algebra dimension | 27 |

The complete set of 20 structural constants (including dim(E₆), dim(E₇), dim(F₄), fund(E₇), |PSL(2,7)|, D_bulk, α_sum, det(g)_den, det(g)_num) is tabulated in Supplement S3, §3.3.

## Appendix B: Derived Structural Constants

| Symbol | Formula | Value |
|--------|---------|-------|
| p₂ | dim(G₂)/dim(K₇) | 2 |
| Weyl | From W(E₈) factorization | 5 |
| N_gen | Index theorem | 3 |
| H* | b₂ + b₃ + 1 | 99 |
| tau | (496 x 21)/(27 x 99) | 3472/891 |
| kappa_T | 1/(b₃ - dim(G₂) - p₂) | 1/61 |
| det(g) | p₂ + 1/(b₂ + dim(G₂) - N_gen) | 65/32 |

## Appendix C: Supplement Reference

| Supplement | Content | Pages |
|------------|---------|-------|
| S1: Foundations | E₈, G₂, K₇ construction details | 27 |
| S2: Derivations | Complete proofs of 33 Type I relations | 38 |
| S3: Observable Dataset | Full 95-entry table, type classification, sensitivity | 15 |
| S4: Sieve diagnostics | Archived coincidence-probability tables and null distributions (3.4 §7.5) | 6 |

---

*The K₇ Framework (formerly GIFT)*
