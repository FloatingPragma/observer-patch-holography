import ObserverPatchHolography.CollarLayer

/-!
# #544 state-side layer: T0 — the identity-channel no-go with real states

This file crosses the rail that `CollarLayer.lean` deliberately does not:
it puts *states* (positive semidefinite, trace-one matrices), the matrix
logarithm (continuous functional calculus), Umegaki relative entropy, and
the closure defect of the paper's `def:closure-defect` on the ℂ-lift of
the model collar layer.  The HARD RAIL of `CollarLayer.lean` ("no
C*-analysis, no GNS/spectra/Hilbert spaces") remains true OF THAT FILE;
the analytic content lives here and only here.

Scope (T0 only, per the scoping report `oph-544-state-side-scoping-report-
2026-07-17.md` §1.1 and §4, bricks S1–S4):

* **S1** — ℂ-lift of the model layer: `CollarC`, `uuC`/`XXC`/`pPlusC`,
  the lifted layer `modelLayerC`, and the cross-cut/flux/invariance facts,
  transferred from the kernel-`decide`d ℤ facts through the entrywise
  cast ring hom (no re-proving over ℂ).
* **S2** — `DensityMatrix`, the spectral projections of the commuting
  pair `{uuC, XXC}`, and the Gibbs states of the retained family via
  `NormedSpace.exp`.
* **S3** — matrix log via the continuous functional calculus
  (`CFC.log = cfc Real.log`; Mathlib has no `Matrix.log`), Umegaki
  relative entropy, `relEntropy_self`, and the family Klein inequality
  (relative entropy is nonnegative *between Gibbs states of the retained
  family*) — the commuting family reduces it to the classical
  finite-probability Gibbs inequality; no Kubo–Mori machinery.
* **S4** — admissible channels (linear, positive, trace-preserving — the
  requirements named in `rem:msascope`, nothing more), the closure defect
  of `def:closure-defect`, and the T0 theorem
  `stateSide_currentAxioms_cannot_force`.

**Claim discipline (non-negotiable).**  T0 *sharpens* #544; it does NOT
close it.  It shows the paper's current admissible-channel requirements
cannot force the central-interface collar clause, because the identity
channel is admissible and discharges refinement closure with closure
defect exactly zero on a non-central family.  Forcing remains impossible
until the admissible class is strengthened at the axiom level (T2 —
explicitly out of scope here).

Anti-smuggling discipline (report §3): `AdmissibleChannel` is
witness-free — its definition mentions positivity and trace preservation
only, never `uuC`, `XXC`, the flux sector, or the clause.
-/

namespace OPH

open Kronecker Matrix

/-! ## S1 — the ℂ-lift of the model layer

Everything is transferred from the kernel-`decide`d ℤ-matrix facts of
`CollarLayer.lean` through the entrywise cast ring hom `liftC`; nothing
is re-proved entrywise over ℂ. -/

/-- The state-side collar algebra: `M₂(ℂ) ⊗ M₂(ℂ)`. -/
abbrev CollarC : Type := Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ

/-- Entrywise scalar lift `ℤ → ℂ` of the model collar ring, as a ring hom. -/
noncomputable def liftC : CollarM →+* CollarC :=
  (Int.castRingHom ℂ).mapMatrix

@[simp] theorem liftC_apply (m : CollarM) (i j : Fin 2 × Fin 2) :
    liftC m i j = (m i j : ℂ) := rfl

/-- The boundary charge `u ⊗ u` over ℂ. -/
noncomputable def uuC : CollarC := liftC uu

/-- The invariant-but-non-central cross coupling `X ⊗ X` over ℂ. -/
noncomputable def XXC : CollarC := liftC XX

/-- The sector projector `(1 + u⊗u)/2` — the element whose image was
    ℤ-impossible for any integral flux retraction
    (`no_integral_flux_retraction`); over ℂ it is an honest projection. -/
noncomputable def pPlusC : CollarC := (2 : ℂ)⁻¹ • (1 + uuC)

/-! ### ℤ-side kernel-`decide` facts (new ones needed by the state layer) -/

theorem uu_sq : uu * uu = 1 := by decide

theorem XX_sq : XX * XX = 1 := by decide

theorem uu_symm : ∀ i j, uu i j = uu j i := by decide

theorem XX_symm : ∀ i j, XX i j = XX j i := by decide

theorem uu_trace_eq : uu.trace = 0 := by decide

theorem XX_trace_eq : XX.trace = 0 := by decide

theorem uu_mul_XX_trace_eq : (uu * XX).trace = 0 := by decide

/-! ### Transferred ℂ facts -/

theorem uuC_sq : uuC * uuC = 1 := by
  rw [uuC, ← map_mul, uu_sq, map_one]

theorem XXC_sq : XXC * XXC = 1 := by
  rw [XXC, ← map_mul, XX_sq, map_one]

theorem uuC_comm_XXC : uuC * XXC = XXC * uuC := by
  rw [uuC, XXC, ← map_mul, ← map_mul, uu_comm_XX]

theorem liftC_trace (m : CollarM) : (liftC m).trace = (m.trace : ℂ) := by
  simp [Matrix.trace, Matrix.diag]

theorem uuC_trace : uuC.trace = 0 := by
  rw [uuC, liftC_trace, uu_trace_eq, Int.cast_zero]

theorem XXC_trace : XXC.trace = 0 := by
  rw [XXC, liftC_trace, XX_trace_eq, Int.cast_zero]

theorem uuC_mul_XXC_trace : (uuC * XXC).trace = 0 := by
  rw [uuC, XXC, ← map_mul, liftC_trace, uu_mul_XX_trace_eq, Int.cast_zero]

theorem liftC_isHermitian_of_symm {m : CollarM} (h : ∀ i j, m i j = m j i) :
    (liftC m).IsHermitian := by
  ext i j
  simp only [Matrix.conjTranspose_apply, liftC_apply]
  rw [h j i]
  exact star_intCast _

theorem uuC_isHermitian : uuC.IsHermitian := liftC_isHermitian_of_symm uu_symm

theorem XXC_isHermitian : XXC.IsHermitian := liftC_isHermitian_of_symm XX_symm

/-- The sector projector is idempotent over ℂ — the direct witness that
    the ℤ-integrality obstruction was a feature of the ring of scalars,
    not of the operator content. -/
theorem pPlusC_idem : pPlusC * pPlusC = pPlusC := by
  rw [pPlusC, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [add_mul, one_mul, mul_add, mul_one, uuC_sq]
  rw [show (1 : CollarC) + uuC + (uuC + 1) = (2 : ℂ) • (1 + uuC) by
    rw [two_smul]
    abel]
  rw [smul_smul]
  norm_num

/-! ### The lifted collar layer -/

/-- Left half-collar embedding over ℂ. -/
noncomputable def kronLeftC : Matrix (Fin 2) (Fin 2) ℂ →+* CollarC where
  toFun m := m ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℂ)
  map_one' := by rw [Matrix.one_kronecker_one]
  map_mul' m n := by rw [← Matrix.mul_kronecker_mul, one_mul]
  map_zero' := by rw [Matrix.zero_kronecker]
  map_add' m n := by rw [Matrix.add_kronecker]

/-- Right half-collar embedding over ℂ. -/
noncomputable def kronRightC : Matrix (Fin 2) (Fin 2) ℂ →+* CollarC where
  toFun n := (1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ n
  map_one' := by rw [Matrix.one_kronecker_one]
  map_mul' m n := by rw [← Matrix.mul_kronecker_mul, one_mul]
  map_zero' := by rw [Matrix.kronecker_zero]
  map_add' m n := by rw [Matrix.kronecker_add]

/-- Entrywise lift commutes with the Kronecker product. -/
theorem liftC_kronecker (a b : Matrix (Fin 2) (Fin 2) ℤ) :
    liftC (a ⊗ₖ b) = (a.map (Int.cast : ℤ → ℂ)) ⊗ₖ (b.map (Int.cast : ℤ → ℂ)) := by
  ext ⟨i, i'⟩ ⟨j, j'⟩
  simp [Matrix.kroneckerMap_apply, Matrix.map_apply]

theorem uuC_eq_kron :
    uuC = (uMat.map (Int.cast : ℤ → ℂ)) ⊗ₖ (uMat.map (Int.cast : ℤ → ℂ)) := by
  rw [uuC, uu, liftC_kronecker]

theorem XXC_eq_kron :
    XXC = (xMat.map (Int.cast : ℤ → ℂ)) ⊗ₖ (xMat.map (Int.cast : ℤ → ℂ)) := by
  rw [XXC, XX, liftC_kronecker]

/-- The model collar layer over ℂ: same interface data as `modelLayer`,
    lifted scalars. -/
noncomputable def modelLayerC : CollarLayer where
  A := CollarC
  K := Subring.closure {uuC}
  ML := kronLeftC.range
  MR := kronRightC.range
  sides_commute := by
    rintro l ⟨m, rfl⟩ r ⟨n, rfl⟩
    show (m ⊗ₖ 1) * ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ n)
      = ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ n) * (m ⊗ₖ 1)
    rw [← Matrix.mul_kronecker_mul, ← Matrix.mul_kronecker_mul,
      mul_one, one_mul, one_mul, mul_one]
  K_le := by
    intro a ha
    refine Subring.closure_le.mpr ?_ ha
    intro g hg
    rw [Set.mem_singleton_iff] at hg
    subst hg
    have huu : uuC = kronLeftC (uMat.map (Int.cast : ℤ → ℂ))
        * kronRightC (uMat.map (Int.cast : ℤ → ℂ)) := by
      show uuC = (uMat.map (Int.cast : ℤ → ℂ) ⊗ₖ 1)
        * ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ uMat.map (Int.cast : ℤ → ℂ))
      rw [← Matrix.mul_kronecker_mul, mul_one, one_mul, uuC_eq_kron]
    rw [huu]
    exact mul_mem
      (Subring.subset_closure (Or.inl ⟨_, rfl⟩))
      (Subring.subset_closure (Or.inr ⟨_, rfl⟩))

/-- Gauge invariance over ℂ (`K̂_Σ`-invariance on the state-side layer). -/
noncomputable abbrev GaugeInvariantC : CollarC → Prop := modelLayerC.Invariant

/-- Cross-cut over ℂ. -/
noncomputable abbrev CrossCutC : CollarC → Prop := modelLayerC.CrossCut

/-- The flux sector over ℂ. -/
noncomputable abbrev FluxC : Set CollarC := modelLayerC.Flux

/-- Anything commuting with the generator `uuC` is invariant (mirror of
    `invariant_of_comm_uu`; the proof is scalar-independent). -/
theorem invariant_of_comm_uuC {a : CollarC} (h : uuC * a = a * uuC) :
    GaugeInvariantC a := by
  intro k hk
  have hle : Subring.closure {uuC} ≤ Subring.centralizer {a} := by
    rw [Subring.closure_le]
    intro g hg
    rw [Set.mem_singleton_iff] at hg
    subst hg
    show uuC ∈ Subring.centralizer {a}
    rw [Subring.mem_centralizer_iff]
    intro m hm
    rw [Set.mem_singleton_iff] at hm
    subst hm
    exact h.symm
  have hk' := hle hk
  rw [Subring.mem_centralizer_iff] at hk'
  exact (hk' a rfl).symm

theorem uuC_invariant : GaugeInvariantC uuC := invariant_of_comm_uuC rfl

theorem XXC_invariant : GaugeInvariantC XXC :=
  invariant_of_comm_uuC uuC_comm_XXC

/-! ### Cross-cut and (non-)centrality over ℂ (mirrors of the ℤ proofs,
same entry choices; the closing numeric contradictions move from `omega`
and `decide` to `norm_num` / `one_ne_zero`). -/

theorem uuC_notMem_MLC : uuC ∉ modelLayerC.ML := by
  rintro ⟨m, hm⟩
  have h1 : m 0 0 * (1 : Matrix (Fin 2) (Fin 2) ℂ) 0 0 = uuC (0, 0) (0, 0) :=
    congrFun (congrFun hm ((0 : Fin 2), (0 : Fin 2))) ((0 : Fin 2), (0 : Fin 2))
  have h2 : m 0 0 * (1 : Matrix (Fin 2) (Fin 2) ℂ) 1 1 = uuC (0, 1) (0, 1) :=
    congrFun (congrFun hm ((0 : Fin 2), (1 : Fin 2))) ((0 : Fin 2), (1 : Fin 2))
  have e1 : uuC (0, 0) (0, 0) = (1 : ℂ) := by
    have h : uu (0, 0) (0, 0) = 1 := by decide
    simp [uuC, h]
  have e2 : uuC (0, 1) (0, 1) = (-1 : ℂ) := by
    have h : uu (0, 1) (0, 1) = -1 := by decide
    simp [uuC, h]
  rw [e1, Matrix.one_apply_eq, mul_one] at h1
  rw [e2, Matrix.one_apply_eq, mul_one] at h2
  rw [h1] at h2
  norm_num at h2

theorem uuC_notMem_MRC : uuC ∉ modelLayerC.MR := by
  rintro ⟨n, hn⟩
  have h1 : (1 : Matrix (Fin 2) (Fin 2) ℂ) 0 0 * n 0 0 = uuC (0, 0) (0, 0) :=
    congrFun (congrFun hn ((0 : Fin 2), (0 : Fin 2))) ((0 : Fin 2), (0 : Fin 2))
  have h2 : (1 : Matrix (Fin 2) (Fin 2) ℂ) 1 1 * n 0 0 = uuC (1, 0) (1, 0) :=
    congrFun (congrFun hn ((1 : Fin 2), (0 : Fin 2))) ((1 : Fin 2), (0 : Fin 2))
  have e1 : uuC (0, 0) (0, 0) = (1 : ℂ) := by
    have h : uu (0, 0) (0, 0) = 1 := by decide
    simp [uuC, h]
  have e2 : uuC (1, 0) (1, 0) = (-1 : ℂ) := by
    have h : uu (1, 0) (1, 0) = -1 := by decide
    simp [uuC, h]
  rw [e1, Matrix.one_apply_eq, one_mul] at h1
  rw [e2, Matrix.one_apply_eq, one_mul] at h2
  rw [h1] at h2
  norm_num at h2

theorem XXC_notMem_MLC : XXC ∉ modelLayerC.ML := by
  rintro ⟨m, hm⟩
  have h1 : m 0 1 * (1 : Matrix (Fin 2) (Fin 2) ℂ) 1 0 = XXC (0, 1) (1, 0) :=
    congrFun (congrFun hm ((0 : Fin 2), (1 : Fin 2))) ((1 : Fin 2), (0 : Fin 2))
  have e1 : XXC (0, 1) (1, 0) = (1 : ℂ) := by
    have h : XX (0, 1) (1, 0) = 1 := by decide
    simp [XXC, h]
  have e0 : (1 : Matrix (Fin 2) (Fin 2) ℂ) 1 0 = 0 :=
    Matrix.one_apply_ne (by decide)
  rw [e1, e0, mul_zero] at h1
  exact one_ne_zero h1.symm

theorem XXC_notMem_MRC : XXC ∉ modelLayerC.MR := by
  rintro ⟨n, hn⟩
  have h1 : (1 : Matrix (Fin 2) (Fin 2) ℂ) 0 1 * n 1 0 = XXC (0, 1) (1, 0) :=
    congrFun (congrFun hn ((0 : Fin 2), (1 : Fin 2))) ((1 : Fin 2), (0 : Fin 2))
  have e1 : XXC (0, 1) (1, 0) = (1 : ℂ) := by
    have h : XX (0, 1) (1, 0) = 1 := by decide
    simp [XXC, h]
  have e0 : (1 : Matrix (Fin 2) (Fin 2) ℂ) 0 1 = 0 :=
    Matrix.one_apply_ne (by decide)
  rw [e1, e0, zero_mul] at h1
  exact one_ne_zero h1.symm

/-- The boundary charge is cross-cut over ℂ. -/
theorem uuC_crossCut : CrossCutC uuC := by
  rintro (h | h)
  · exact uuC_notMem_MLC h
  · exact uuC_notMem_MRC h

/-- The coupling is cross-cut over ℂ. -/
theorem XXC_crossCut : CrossCutC XXC := by
  rintro (h | h)
  · exact XXC_notMem_MLC h
  · exact XXC_notMem_MRC h

/-- The boundary charge is a flux term over ℂ. -/
theorem uuC_mem_fluxC : uuC ∈ FluxC :=
  ⟨Subring.subset_closure rfl, uuC_invariant⟩

/-- Diagonal matrices form a subring of the ℂ collar ring. -/
def diagonalSubringC : Subring CollarC where
  carrier := {m | ∀ i j, i ≠ j → m i j = 0}
  zero_mem' := fun _ _ _ => Matrix.zero_apply _ _
  one_mem' := fun _ _ hij => Matrix.one_apply_ne hij
  add_mem' := fun ha hb i j hij => by
    rw [Matrix.add_apply, ha i j hij, hb i j hij, add_zero]
  neg_mem' := fun ha i j hij => by
    rw [Matrix.neg_apply, ha i j hij, neg_zero]
  mul_mem' := fun {a b} ha hb i j hij => by
    rw [Matrix.mul_apply]
    apply Finset.sum_eq_zero
    intro k _
    rcases eq_or_ne k i with hk | hk
    · rw [hk, hb i j hij, mul_zero]
    · rw [ha i k hk.symm, zero_mul]

theorem uuC_diagonal : ∀ i j : Fin 2 × Fin 2, i ≠ j → uuC i j = 0 := by
  intro i j hij
  show ((uu i j : ℤ) : ℂ) = 0
  rw [uu_diagonal i j hij, Int.cast_zero]

/-- The coupling is NOT a flux term over ℂ: `K = ⟨uuC⟩` is diagonal and
    `XXC` is not (mirror of `XX_notMem_flux`). -/
theorem XXC_notMem_fluxC : XXC ∉ FluxC := by
  rintro ⟨hK, -⟩
  have hle : Subring.closure {uuC} ≤ diagonalSubringC := by
    rw [Subring.closure_le]
    intro g hg
    rw [Set.mem_singleton_iff] at hg
    subst hg
    exact uuC_diagonal
  have hdiag := hle hK
  have h0 := hdiag ((0 : Fin 2), (1 : Fin 2)) ((1 : Fin 2), (0 : Fin 2)) (by decide)
  have h1 : XXC (0, 1) (1, 0) = (1 : ℂ) := by
    have h : XX (0, 1) (1, 0) = 1 := by decide
    simp [XXC, h]
  rw [h1] at h0
  exact one_ne_zero h0

/-! ## S2 — density matrices and the Gibbs states of the retained family

The retained family `{uuC, XXC}` commutes (`uuC_comm_XXC`), so its joint
spectral decomposition is carried by the four projections
`specProjP (s,t) = ¼ (1 + ε_s uuC)(1 + ε_t XXC)` — pure algebra from
`uuC² = XXC² = 1`; no eigenbasis, no `√2`, no entry computations. -/

open scoped ComplexOrder

/-- A state on the collar algebra: positive semidefinite, trace one. -/
structure DensityMatrix where
  /-- The underlying matrix. -/
  ρ : CollarC
  /-- Positivity. -/
  posSemidef : ρ.PosSemidef
  /-- Normalization. -/
  trace_one : ρ.trace = 1

/-- Real scalars act on the collar algebra through their complex cast. -/
theorem ofReal_smul_collarC (r : ℝ) (M : CollarC) :
    ((r : ℝ) : ℂ) • M = r • M := by
  ext i j
  simp [Matrix.smul_apply, Complex.real_smul]

/-- Products of nonnegative complex numbers are nonnegative (the two
    factors are real, by `Complex.nonneg_iff`). -/
theorem complex_mul_nonneg {a b : ℂ} (ha : 0 ≤ a) (hb : 0 ≤ b) : 0 ≤ a * b := by
  rw [Complex.nonneg_iff] at ha hb ⊢
  refine ⟨?_, ?_⟩
  · rw [Complex.mul_re, ← ha.2, ← hb.2]
    simpa using mul_nonneg ha.1 hb.1
  · rw [Complex.mul_im, ← ha.2, ← hb.2]
    simp

/-- Positive semidefiniteness is preserved by nonnegative real scalars
    (stated for the complex cast of the scalar, which is how the spectral
    coefficients appear). -/
theorem posSemidef_ofReal_smul {r : ℝ} (hr : 0 ≤ r) {M : CollarC}
    (hM : M.PosSemidef) : (((r : ℝ) : ℂ) • M).PosSemidef := by
  refine ⟨?_, fun x => ?_⟩
  · show _ᴴ = _
    rw [Matrix.conjTranspose_smul, hM.1.eq]
    congr 1
    rw [Complex.star_def, Complex.conj_ofReal]
  · have h := hM.2 x
    have key : (x.sum fun i xi => x.sum fun j xj =>
          star xi * (((r : ℝ) : ℂ) • M) i j * xj)
        = ((r : ℝ) : ℂ)
          * x.sum fun i xi => x.sum fun j xj => star xi * M i j * xj := by
      rw [Finsupp.mul_sum]
      refine Finsupp.sum_congr fun i _ => ?_
      rw [Finsupp.mul_sum]
      refine Finsupp.sum_congr fun j _ => ?_
      rw [Matrix.smul_apply, smul_eq_mul]
      ring
    rw [key]
    exact complex_mul_nonneg (Complex.zero_le_real.mpr hr) h

/-- Sign labels for the two-point spectra of `uuC` and `XXC`. -/
noncomputable def sgnR : Fin 2 → ℝ := ![1, -1]

@[simp] theorem sgnR_zero : sgnR 0 = 1 := rfl

@[simp] theorem sgnR_one : sgnR 1 = -1 := rfl

theorem sgnR_sq (s : Fin 2) : sgnR s * sgnR s = 1 := by
  fin_cases s <;> norm_num [sgnR]

theorem sgnR_mul_ne {s s' : Fin 2} (h : s ≠ s') : sgnR s * sgnR s' = -1 := by
  fin_cases s <;> fin_cases s' <;> simp_all [sgnR]

/-- The `uuC` half-projections. -/
noncomputable def quP (s : Fin 2) : CollarC :=
  (2 : ℂ)⁻¹ • (1 + ((sgnR s : ℝ) : ℂ) • uuC)

/-- The `XXC` half-projections. -/
noncomputable def qxP (t : Fin 2) : CollarC :=
  (2 : ℂ)⁻¹ • (1 + ((sgnR t : ℝ) : ℂ) • XXC)

/-- The joint spectral projections of the commuting pair `{uuC, XXC}`. -/
noncomputable def specProjP (p : Fin 2 × Fin 2) : CollarC := quP p.1 * qxP p.2

/-- Casts of the sign squares. -/
theorem sgnC_sq (s : Fin 2) : ((sgnR s : ℝ) : ℂ) * ((sgnR s : ℝ) : ℂ) = 1 := by
  rw [← Complex.ofReal_mul, sgnR_sq, Complex.ofReal_one]

theorem sgnC_mul_ne {s s' : Fin 2} (h : s ≠ s') :
    ((sgnR s : ℝ) : ℂ) * ((sgnR s' : ℝ) : ℂ) = -1 := by
  rw [← Complex.ofReal_mul, sgnR_mul_ne h]
  norm_num

theorem quP_mul (s s' : Fin 2) :
    quP s * quP s' = if s = s' then quP s else 0 := by
  rcases eq_or_ne s s' with rfl | h
  · rw [if_pos rfl, quP, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [add_mul, one_mul, mul_add, mul_one, smul_mul_smul_comm, uuC_sq,
      sgnC_sq, one_smul]
    rw [show ((1 : CollarC) + ((sgnR s : ℝ) : ℂ) • uuC)
          + (((sgnR s : ℝ) : ℂ) • uuC + 1)
        = (2 : ℂ) • (1 + ((sgnR s : ℝ) : ℂ) • uuC) by
      rw [two_smul]
      abel]
    rw [smul_smul]
    norm_num
  · rw [if_neg h, quP, quP, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [add_mul, one_mul, mul_add, mul_one, smul_mul_smul_comm, uuC_sq,
      sgnC_mul_ne h]
    rw [show ((1 : CollarC) + ((sgnR s' : ℝ) : ℂ) • uuC)
          + (((sgnR s : ℝ) : ℂ) • uuC + (-1 : ℂ) • 1) = ((((sgnR s : ℝ) : ℂ))
            + (((sgnR s' : ℝ) : ℂ))) • uuC by
      rw [neg_smul, one_smul, add_smul]
      abel]
    rw [show (((sgnR s : ℝ) : ℂ)) + (((sgnR s' : ℝ) : ℂ)) = 0 by
      rw [← Complex.ofReal_add,
        show sgnR s + sgnR s' = 0 by
          fin_cases s <;> fin_cases s' <;> simp_all [sgnR],
        Complex.ofReal_zero]]
    rw [zero_smul, smul_zero]

theorem qxP_mul (t t' : Fin 2) :
    qxP t * qxP t' = if t = t' then qxP t else 0 := by
  rcases eq_or_ne t t' with rfl | h
  · rw [if_pos rfl, qxP, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [add_mul, one_mul, mul_add, mul_one, smul_mul_smul_comm, XXC_sq,
      sgnC_sq, one_smul]
    rw [show ((1 : CollarC) + ((sgnR t : ℝ) : ℂ) • XXC)
          + (((sgnR t : ℝ) : ℂ) • XXC + 1)
        = (2 : ℂ) • (1 + ((sgnR t : ℝ) : ℂ) • XXC) by
      rw [two_smul]
      abel]
    rw [smul_smul]
    norm_num
  · rw [if_neg h, qxP, qxP, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [add_mul, one_mul, mul_add, mul_one, smul_mul_smul_comm, XXC_sq,
      sgnC_mul_ne h]
    rw [show ((1 : CollarC) + ((sgnR t' : ℝ) : ℂ) • XXC)
          + (((sgnR t : ℝ) : ℂ) • XXC + (-1 : ℂ) • 1) = ((((sgnR t : ℝ) : ℂ))
            + (((sgnR t' : ℝ) : ℂ))) • XXC by
      rw [neg_smul, one_smul, add_smul]
      abel]
    rw [show (((sgnR t : ℝ) : ℂ)) + (((sgnR t' : ℝ) : ℂ)) = 0 by
      rw [← Complex.ofReal_add,
        show sgnR t + sgnR t' = 0 by
          fin_cases t <;> fin_cases t' <;> simp_all [sgnR],
        Complex.ofReal_zero]]
    rw [zero_smul, smul_zero]

theorem quP_comm_qxP (s t : Fin 2) : Commute (quP s) (qxP t) := by
  have huu : Commute uuC XXC := uuC_comm_XXC
  have h : Commute ((1 : CollarC) + ((sgnR s : ℝ) : ℂ) • uuC)
      ((1 : CollarC) + ((sgnR t : ℝ) : ℂ) • XXC) := by
    apply Commute.add_left
    · apply Commute.add_right
      · exact Commute.one_left _
      · exact (Commute.one_left _).smul_right _
    · apply Commute.add_right
      · exact (Commute.one_right _).smul_left _
      · exact (huu.smul_left _).smul_right _
  exact (h.smul_left _).smul_right _

theorem specProjP_mul (p q : Fin 2 × Fin 2) :
    specProjP p * specProjP q = if p = q then specProjP p else 0 := by
  have key : specProjP p * specProjP q
      = (quP p.1 * quP q.1) * (qxP p.2 * qxP q.2) := by
    rw [specProjP, specProjP]
    rw [mul_assoc (quP p.1) (qxP p.2) (quP q.1 * qxP q.2)]
    rw [← mul_assoc (qxP p.2) (quP q.1) (qxP q.2)]
    rw [← (quP_comm_qxP q.1 p.2).eq]
    rw [mul_assoc (quP q.1) (qxP p.2) (qxP q.2)]
    rw [← mul_assoc (quP p.1) (quP q.1) (qxP p.2 * qxP q.2)]
  rw [key, quP_mul, qxP_mul]
  rcases eq_or_ne p q with rfl | h
  · rw [if_pos rfl, if_pos rfl, if_pos rfl, specProjP]
  · rw [if_neg h]
    have h' : p.1 ≠ q.1 ∨ p.2 ≠ q.2 := by
      rcases eq_or_ne p.1 q.1 with h1 | h1
      · exact Or.inr fun h2 => h (Prod.ext h1 h2)
      · exact Or.inl h1
    rcases h' with h1 | h2
    · rw [if_neg h1, zero_mul]
    · rw [if_neg h2, mul_zero]

theorem specProjP_idem (p : Fin 2 × Fin 2) :
    specProjP p * specProjP p = specProjP p := by
  rw [specProjP_mul, if_pos rfl]

theorem sum_quP : ∑ s : Fin 2, quP s = 1 := by
  rw [Fin.sum_univ_two]
  simp only [quP, sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg,
    one_smul, neg_smul]
  rw [← smul_add,
    show ((1 : CollarC) + uuC) + (1 + -uuC) = (2 : ℂ) • 1 by
      rw [two_smul]
      abel,
    smul_smul]
  norm_num

theorem sum_qxP : ∑ t : Fin 2, qxP t = 1 := by
  rw [Fin.sum_univ_two]
  simp only [qxP, sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg,
    one_smul, neg_smul]
  rw [← smul_add,
    show ((1 : CollarC) + XXC) + (1 + -XXC) = (2 : ℂ) • 1 by
      rw [two_smul]
      abel,
    smul_smul]
  norm_num

theorem sum_specProjP : ∑ p : Fin 2 × Fin 2, specProjP p = 1 := by
  rw [Fintype.sum_prod_type]
  simp only [specProjP]
  have inner : ∀ s : Fin 2, ∑ t : Fin 2, quP s * qxP t = quP s := by
    intro s
    rw [← Finset.mul_sum, sum_qxP, mul_one]
  simp_rw [inner]
  exact sum_quP

theorem quP_isHermitian (s : Fin 2) : (quP s).IsHermitian := by
  show _ᴴ = _
  rw [quP, Matrix.conjTranspose_smul, Matrix.conjTranspose_add,
    Matrix.conjTranspose_one, Matrix.conjTranspose_smul, uuC_isHermitian.eq]
  congr 1
  · norm_num [Complex.star_def]
  · congr 1
    rw [Complex.star_def, Complex.conj_ofReal]

theorem qxP_isHermitian (t : Fin 2) : (qxP t).IsHermitian := by
  show _ᴴ = _
  rw [qxP, Matrix.conjTranspose_smul, Matrix.conjTranspose_add,
    Matrix.conjTranspose_one, Matrix.conjTranspose_smul, XXC_isHermitian.eq]
  congr 1
  · norm_num [Complex.star_def]
  · congr 1
    rw [Complex.star_def, Complex.conj_ofReal]

theorem specProjP_isHermitian (p : Fin 2 × Fin 2) :
    (specProjP p).IsHermitian := by
  show _ᴴ = _
  rw [specProjP, Matrix.conjTranspose_mul, (quP_isHermitian p.1).eq,
    (qxP_isHermitian p.2).eq, ← (quP_comm_qxP p.1 p.2).eq]

theorem specProjP_posSemidef (p : Fin 2 × Fin 2) : (specProjP p).PosSemidef := by
  have h : specProjP p = (specProjP p)ᴴ * specProjP p := by
    rw [(specProjP_isHermitian p).eq, specProjP_idem]
  rw [h]
  exact Matrix.posSemidef_conjTranspose_mul_self _

theorem trace_one_collarC : (1 : CollarC).trace = 4 := by
  rw [Matrix.trace_one]
  norm_num

theorem specProjP_trace (p : Fin 2 × Fin 2) : (specProjP p).trace = 1 := by
  rw [specProjP, quP, qxP, smul_mul_smul_comm]
  rw [add_mul, one_mul, mul_add, mul_one, smul_mul_smul_comm]
  simp only [Matrix.trace_smul, Matrix.trace_add, trace_one_collarC,
    uuC_trace, XXC_trace, uuC_mul_XXC_trace, smul_eq_mul]
  norm_num

theorem uuC_mul_quP (s : Fin 2) :
    uuC * quP s = ((sgnR s : ℝ) : ℂ) • quP s := by
  rw [quP, Matrix.mul_smul, mul_add, mul_one, Matrix.mul_smul, uuC_sq]
  conv_rhs => rw [smul_smul, mul_comm, ← smul_smul]
  congr 1
  rw [smul_add, smul_smul, sgnC_sq, one_smul, add_comm]

theorem XXC_mul_qxP (t : Fin 2) :
    XXC * qxP t = ((sgnR t : ℝ) : ℂ) • qxP t := by
  rw [qxP, Matrix.mul_smul, mul_add, mul_one, Matrix.mul_smul, XXC_sq]
  conv_rhs => rw [smul_smul, mul_comm, ← smul_smul]
  congr 1
  rw [smul_add, smul_smul, sgnC_sq, one_smul, add_comm]

/-! ### The retained family, its Hamiltonians, and the spectral data -/

/-- The retained non-central family of the T0 witness: the boundary
    charge and the invariant-but-non-central coupling. -/
noncomputable def SFam : Fin 2 → CollarC := ![uuC, XXC]

@[simp] theorem SFam_zero : SFam 0 = uuC := rfl

@[simp] theorem SFam_one : SFam 1 = XXC := rfl

/-- The joint eigenvalue of `λ₀ uuC + λ₁ XXC` on the `(s,t)` sector. -/
noncomputable def eigval (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) : ℝ :=
  lam 0 * sgnR p.1 + lam 1 * sgnR p.2

theorem sum_sgn_smul_quP :
    ∑ s : Fin 2, ((sgnR s : ℝ) : ℂ) • quP s = uuC := by
  rw [Fin.sum_univ_two]
  simp only [quP, sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg,
    one_smul, neg_smul]
  rw [← smul_neg, ← smul_add,
    show ((1 : CollarC) + uuC) + -(1 + -uuC) = (2 : ℂ) • uuC by
      rw [two_smul]
      abel,
    smul_smul]
  norm_num

theorem sum_sgn_smul_qxP :
    ∑ t : Fin 2, ((sgnR t : ℝ) : ℂ) • qxP t = XXC := by
  rw [Fin.sum_univ_two]
  simp only [qxP, sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg,
    one_smul, neg_smul]
  rw [← smul_neg, ← smul_add,
    show ((1 : CollarC) + XXC) + -(1 + -XXC) = (2 : ℂ) • XXC by
      rw [two_smul]
      abel,
    smul_smul]
  norm_num

/-- Spectral decomposition of the Hamiltonian of the retained family. -/
theorem Ham_eq (lam : Fin 2 → ℝ) :
    ∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a
      = ∑ p : Fin 2 × Fin 2, ((eigval lam p : ℝ) : ℂ) • specProjP p := by
  have expand : ∀ p : Fin 2 × Fin 2, ((eigval lam p : ℝ) : ℂ) • specProjP p
      = ((lam 0 * sgnR p.1 : ℝ) : ℂ) • specProjP p
        + ((lam 1 * sgnR p.2 : ℝ) : ℂ) • specProjP p := by
    intro p
    rw [← add_smul, ← Complex.ofReal_add, eigval]
  rw [Finset.sum_congr rfl fun p _ => expand p, Finset.sum_add_distrib]
  have h1 : ∑ p : Fin 2 × Fin 2, ((lam 0 * sgnR p.1 : ℝ) : ℂ) • specProjP p
      = ((lam 0 : ℝ) : ℂ) • uuC := by
    rw [Fintype.sum_prod_type]
    simp only [specProjP]
    have inner : ∀ s : Fin 2,
        ∑ t : Fin 2, ((lam 0 * sgnR s : ℝ) : ℂ) • (quP s * qxP t)
          = ((lam 0 * sgnR s : ℝ) : ℂ) • quP s := by
      intro s
      rw [← Finset.smul_sum, ← Finset.mul_sum, sum_qxP, mul_one]
    simp_rw [inner]
    have split : ∀ s : Fin 2, ((lam 0 * sgnR s : ℝ) : ℂ) • quP s
        = ((lam 0 : ℝ) : ℂ) • (((sgnR s : ℝ) : ℂ) • quP s) := by
      intro s
      rw [smul_smul, ← Complex.ofReal_mul]
    simp_rw [split]
    rw [← Finset.smul_sum, sum_sgn_smul_quP]
  have h2 : ∑ p : Fin 2 × Fin 2, ((lam 1 * sgnR p.2 : ℝ) : ℂ) • specProjP p
      = ((lam 1 : ℝ) : ℂ) • XXC := by
    rw [Fintype.sum_prod_type]
    simp only [specProjP]
    rw [Finset.sum_comm]
    have inner : ∀ t : Fin 2,
        ∑ s : Fin 2, ((lam 1 * sgnR t : ℝ) : ℂ) • (quP s * qxP t)
          = ((lam 1 * sgnR t : ℝ) : ℂ) • qxP t := by
      intro t
      rw [← Finset.smul_sum, ← Finset.sum_mul, sum_quP, one_mul]
    simp_rw [inner]
    have split : ∀ t : Fin 2, ((lam 1 * sgnR t : ℝ) : ℂ) • qxP t
        = ((lam 1 : ℝ) : ℂ) • (((sgnR t : ℝ) : ℂ) • qxP t) := by
      intro t
      rw [smul_smul, ← Complex.ofReal_mul]
    simp_rw [split]
    rw [← Finset.smul_sum, sum_sgn_smul_qxP]
  rw [h1, h2, Fin.sum_univ_two, SFam_zero, SFam_one]

/-- Spectral decomposition of the negated Hamiltonian. -/
theorem neg_Ham_eq (lam : Fin 2 → ℝ) :
    -(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a)
      = ∑ p : Fin 2 × Fin 2, ((-(eigval lam p) : ℝ) : ℂ) • specProjP p := by
  rw [Ham_eq, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Complex.ofReal_neg, neg_smul]

/-! ### The exponential of the Hamiltonian, spectrally -/

open NormedSpace in
/-- `exp` of a linear combination of the joint spectral projections is the
    corresponding combination of scalar exponentials: the tsum argument,
    using only orthogonality and completeness of the projection family. -/
theorem exp_sum_smul_specProjP (c : Fin 2 × Fin 2 → ℂ) :
    exp (∑ p : Fin 2 × Fin 2, c p • specProjP p)
      = ∑ p : Fin 2 × Fin 2, Complex.exp (c p) • specProjP p := by
  have hpow : ∀ n : ℕ,
      (∑ p : Fin 2 × Fin 2, c p • specProjP p) ^ n
        = ∑ p : Fin 2 × Fin 2, (c p) ^ n • specProjP p := by
    intro n
    induction n with
    | zero => simpa using sum_specProjP.symm
    | succ n ih =>
      rw [pow_succ, ih, Finset.sum_mul]
      refine Finset.sum_congr rfl fun p _ => ?_
      rw [Finset.mul_sum]
      rw [Finset.sum_eq_single p
        (fun q _ hq => by
          rw [smul_mul_smul_comm, specProjP_mul, if_neg (Ne.symm hq), smul_zero])
        (fun h => absurd (Finset.mem_univ p) h)]
      rw [smul_mul_smul_comm, specProjP_mul, if_pos rfl, ← pow_succ]
  rw [exp_eq_tsum ℂ]
  simp only [hpow]
  have hsummand : ∀ n : ℕ,
      ((n.factorial : ℂ)⁻¹ • ∑ p : Fin 2 × Fin 2, (c p) ^ n • specProjP p)
        = ∑ p : Fin 2 × Fin 2,
            ((n.factorial : ℂ)⁻¹ * (c p) ^ n) • specProjP p := by
    intro n
    rw [Finset.smul_sum]
    refine Finset.sum_congr rfl fun _ _ => ?_
    rw [smul_smul]
  simp only [hsummand]
  have hscalar : ∀ p : Fin 2 × Fin 2,
      Summable (fun n : ℕ => (n.factorial : ℂ)⁻¹ * (c p) ^ n) := by
    intro p
    have h := expSeries_div_summable (𝔸 := ℂ) (c p)
    refine h.congr fun n => ?_
    rw [div_eq_mul_inv, mul_comm]
  rw [Summable.tsum_finsetSum
    (fun p _ => (hscalar p).smul_const (specProjP p))]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [(hscalar p).tsum_smul_const]
  congr 1
  have h := congrFun (exp_eq_tsum (𝕂 := ℂ) (𝔸 := ℂ)) (c p)
  rw [← Complex.exp_eq_exp_ℂ] at h
  rw [h]
  refine tsum_congr fun n => ?_
  rw [smul_eq_mul]

/-! ### The Gibbs states of the retained family -/

open NormedSpace in
/-- The (normalized) Gibbs state of a two-density family at multipliers
    `lam`, via `NormedSpace.exp`.  Total definition; its state properties
    are proved for the commuting witness family `SFam`. -/
noncomputable def gibbsM (S : Fin 2 → CollarC) (lam : Fin 2 → ℝ) : CollarC :=
  ((exp (-(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • S a))).trace)⁻¹
    • exp (-(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • S a))

/-- The partition function of the witness family. -/
noncomputable def partitionZ (lam : Fin 2 → ℝ) : ℝ :=
  ∑ p : Fin 2 × Fin 2, Real.exp (-(eigval lam p))

theorem partitionZ_pos (lam : Fin 2 → ℝ) : 0 < partitionZ lam :=
  Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty

theorem partitionZ_ne_zero (lam : Fin 2 → ℝ) : partitionZ lam ≠ 0 :=
  ne_of_gt (partitionZ_pos lam)

open NormedSpace in
/-- Spectral form of the Boltzmann factor of the witness family. -/
theorem exp_neg_Ham (lam : Fin 2 → ℝ) :
    exp (-(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a))
      = ∑ p : Fin 2 × Fin 2,
          ((Real.exp (-(eigval lam p)) : ℝ) : ℂ) • specProjP p := by
  rw [neg_Ham_eq, exp_sum_smul_specProjP]
  refine Finset.sum_congr rfl fun _ _ => ?_
  rw [Complex.ofReal_exp]

open NormedSpace in
theorem exp_neg_Ham_trace (lam : Fin 2 → ℝ) :
    (exp (-(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a))).trace
      = ((partitionZ lam : ℝ) : ℂ) := by
  rw [exp_neg_Ham, Matrix.trace_sum]
  simp_rw [Matrix.trace_smul, specProjP_trace, smul_eq_mul, mul_one]
  rw [partitionZ, Complex.ofReal_sum]

open NormedSpace in
theorem exp_neg_Ham_posSemidef (lam : Fin 2 → ℝ) :
    (exp (-(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a))).PosSemidef := by
  rw [exp_neg_Ham]
  refine Finset.sum_induction _ _ (fun A B hA hB => hA.add hB)
    Matrix.PosSemidef.zero fun p _ => ?_
  exact posSemidef_ofReal_smul (Real.exp_pos _).le (specProjP_posSemidef p)

/-- The Gibbs state of the witness family is a density matrix. -/
noncomputable def gibbsDM (lam : Fin 2 → ℝ) : DensityMatrix where
  ρ := gibbsM SFam lam
  posSemidef := by
    rw [gibbsM, exp_neg_Ham_trace, ← Complex.ofReal_inv]
    exact posSemidef_ofReal_smul (inv_nonneg.mpr (partitionZ_pos lam).le)
      (exp_neg_Ham_posSemidef lam)
  trace_one := by
    rw [gibbsM, Matrix.trace_smul, exp_neg_Ham_trace, smul_eq_mul,
      inv_mul_cancel₀ (Complex.ofReal_ne_zero.mpr (partitionZ_ne_zero lam))]

@[simp] theorem gibbsDM_rho (lam : Fin 2 → ℝ) :
    (gibbsDM lam).ρ = gibbsM SFam lam := rfl

/-! ## S3 — matrix log, Umegaki relative entropy, and the family Klein
inequality

The matrix logarithm is `CFC.log = cfc Real.log` — the continuous
functional calculus instance on matrices is norm-free (topological), and
the `NormedRing` structure needed by the `CFC.log`/`NormedSpace.exp`
interplay is supplied by file-local `l∞`-operator-norm instances (the
resulting values are norm-independent).

The Klein inequality proved here is the *family* Klein inequality:
relative entropy is nonnegative between Gibbs states of the retained
family.  The commuting family reduces it to the classical
finite-probability Gibbs inequality (`Real.log_le_sub_one_of_pos`) — no
Kubo–Mori/Duhamel machinery, exactly as the scoping report predicted.
The fully general Klein inequality (arbitrary density-matrix pairs) is
NOT needed for T0 and is not claimed. -/

section StateSide

attribute [local instance] Matrix.linftyOpNormedAddCommGroup
  Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- The matrix logarithm on the collar algebra, via the continuous
    functional calculus (Mathlib has no `Matrix.log`). -/
noncomputable def matLog (M : CollarC) : CollarC := CFC.log M

/-- Umegaki relative entropy on raw matrices (total definition; its
    entropy meaning is carried by the density-matrix arguments used in
    the theorems). -/
noncomputable def relEntropyM (A B : CollarC) : ℝ :=
  ((A * (matLog A - matLog B)).trace).re

/-- Umegaki relative entropy of two states. -/
noncomputable def relEntropy (ρ σ : DensityMatrix) : ℝ :=
  relEntropyM ρ.ρ σ.ρ

/-- Relative entropy of any matrix with itself vanishes — the minimizer
    witness for the identity-channel closure defect. -/
theorem relEntropyM_self (A : CollarC) : relEntropyM A A = 0 := by
  rw [relEntropyM, sub_self, mul_zero, Matrix.trace_zero]
  simp

theorem relEntropy_self (ρ : DensityMatrix) : relEntropy ρ ρ = 0 :=
  relEntropyM_self ρ.ρ

/-! ### The log of a Gibbs state -/

/-- The self-adjoint exponent of the Gibbs state:
    `-(λ·S) - log Z • 1`. -/
noncomputable def gibbsExponent (lam : Fin 2 → ℝ) : CollarC :=
  -(∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a)
    - ((Real.log (partitionZ lam) : ℝ) : ℂ) • 1

/-- The scalar spectral data of the exponent. -/
noncomputable def gLog (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) : ℝ :=
  -(eigval lam p) - Real.log (partitionZ lam)

theorem smul_one_eq_sum_specProjP (c : ℂ) :
    c • (1 : CollarC) = ∑ p : Fin 2 × Fin 2, c • specProjP p := by
  rw [← Finset.smul_sum, sum_specProjP]

theorem gibbsExponent_spectral (lam : Fin 2 → ℝ) :
    gibbsExponent lam
      = ∑ p : Fin 2 × Fin 2, ((gLog lam p : ℝ) : ℂ) • specProjP p := by
  rw [gibbsExponent, neg_Ham_eq, smul_one_eq_sum_specProjP,
    ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [← sub_smul, ← Complex.ofReal_sub, gLog]

theorem gibbsExponent_isSelfAdjoint (lam : Fin 2 → ℝ) :
    IsSelfAdjoint (gibbsExponent lam) := by
  rw [gibbsExponent_spectral]
  show _ᴴ = _
  rw [Matrix.conjTranspose_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Matrix.conjTranspose_smul, (specProjP_isHermitian p).eq,
    Complex.star_def, Complex.conj_ofReal]

open NormedSpace in
/-- The Gibbs state is the exponential of its self-adjoint exponent. -/
theorem gibbsM_eq_exp_gibbsExponent (lam : Fin 2 → ℝ) :
    gibbsM SFam lam = exp (gibbsExponent lam) := by
  rw [gibbsExponent_spectral, exp_sum_smul_specProjP,
    gibbsM, exp_neg_Ham_trace, exp_neg_Ham, Finset.smul_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [smul_smul]
  congr 1
  rw [← Complex.ofReal_exp, gLog, Real.exp_sub,
    Real.exp_log (partitionZ_pos lam), Complex.ofReal_div, inv_mul_eq_div]

/-- **Log of a Gibbs state**: `log ω(λ) = -λ·S - log Z • 1` — the
    manuscript's local-Gibbs form, on machine. -/
theorem matLog_gibbsM (lam : Fin 2 → ℝ) :
    matLog (gibbsM SFam lam) = gibbsExponent lam := by
  rw [gibbsM_eq_exp_gibbsExponent, matLog]
  exact CFC.log_exp _ (gibbsExponent_isSelfAdjoint lam)

/-! ### The relative entropy between Gibbs states of the family -/

/-- The spectral probability weights of the Gibbs state. -/
noncomputable def gibbsProb (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) : ℝ :=
  Real.exp (-(eigval lam p)) / partitionZ lam

theorem gibbsProb_pos (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) :
    0 < gibbsProb lam p :=
  div_pos (Real.exp_pos _) (partitionZ_pos lam)

theorem sum_gibbsProb (lam : Fin 2 → ℝ) :
    ∑ p : Fin 2 × Fin 2, gibbsProb lam p = 1 := by
  rw [show (∑ p : Fin 2 × Fin 2, gibbsProb lam p)
      = (∑ p : Fin 2 × Fin 2, Real.exp (-(eigval lam p))) / partitionZ lam from
    (Finset.sum_div _ _ _).symm]
  rw [← partitionZ, div_self (partitionZ_ne_zero lam)]

/-- The spectral scalars of the exponent are the logs of the weights. -/
theorem gLog_eq_log_gibbsProb (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) :
    gLog lam p = Real.log (gibbsProb lam p) := by
  rw [gibbsProb, Real.log_div (Real.exp_ne_zero _) (partitionZ_ne_zero lam),
    Real.log_exp, gLog]

/-- Spectral form of the Gibbs state itself. -/
theorem gibbsM_spectral (lam : Fin 2 → ℝ) :
    gibbsM SFam lam
      = ∑ p : Fin 2 × Fin 2, ((gibbsProb lam p : ℝ) : ℂ) • specProjP p := by
  rw [gibbsM, exp_neg_Ham_trace, exp_neg_Ham, Finset.smul_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [smul_smul]
  congr 1
  rw [gibbsProb, Complex.ofReal_div, inv_mul_eq_div]

/-- Products of spectral combinations multiply their scalars. -/
theorem sum_smul_specProjP_mul (c d : Fin 2 × Fin 2 → ℂ) :
    (∑ p : Fin 2 × Fin 2, c p • specProjP p)
        * (∑ p : Fin 2 × Fin 2, d p • specProjP p)
      = ∑ p : Fin 2 × Fin 2, (c p * d p) • specProjP p := by
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Finset.mul_sum]
  rw [Finset.sum_eq_single p
    (fun q _ hq => by
      rw [smul_mul_smul_comm, specProjP_mul, if_neg (Ne.symm hq), smul_zero])
    (fun h => absurd (Finset.mem_univ p) h)]
  rw [smul_mul_smul_comm, specProjP_mul, if_pos rfl]

/-- The relative entropy between two Gibbs states of the retained family
    is the classical Kullback–Leibler divergence of their spectral
    weights. -/
theorem relEntropyM_gibbs_eq (lam lam' : Fin 2 → ℝ) :
    relEntropyM (gibbsM SFam lam) (gibbsM SFam lam')
      = ∑ p : Fin 2 × Fin 2,
          gibbsProb lam p
            * (Real.log (gibbsProb lam p) - Real.log (gibbsProb lam' p)) := by
  rw [relEntropyM, matLog_gibbsM, matLog_gibbsM,
    gibbsExponent_spectral, gibbsExponent_spectral,
    ← Finset.sum_sub_distrib]
  have hsub : ∀ p : Fin 2 × Fin 2,
      ((gLog lam p : ℝ) : ℂ) • specProjP p - ((gLog lam' p : ℝ) : ℂ) • specProjP p
        = (((gLog lam p - gLog lam' p : ℝ) : ℂ)) • specProjP p := by
    intro p
    rw [← sub_smul, ← Complex.ofReal_sub]
  simp_rw [hsub]
  rw [gibbsM_spectral, sum_smul_specProjP_mul, Matrix.trace_sum]
  simp_rw [Matrix.trace_smul, specProjP_trace, smul_eq_mul, mul_one,
    ← Complex.ofReal_mul]
  rw [← Complex.ofReal_sum, Complex.ofReal_re]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [gLog_eq_log_gibbsProb, gLog_eq_log_gibbsProb]

/-! ### The classical finite Gibbs inequality (Klein, commuting case) -/

/-- Classical finite-probability Klein/Gibbs inequality: the discrete
    KL divergence of two positive weight systems with equal mass is
    nonnegative.  Proof: `log x ≤ x - 1`. -/
theorem finset_klDiv_nonneg {ι : Type*} [Fintype ι] (p q : ι → ℝ)
    (hp : ∀ i, 0 < p i) (hq : ∀ i, 0 < q i)
    (hpq : ∑ i, p i = ∑ i, q i) :
    0 ≤ ∑ i, p i * (Real.log (p i) - Real.log (q i)) := by
  have key : ∀ i : ι, p i * (Real.log (q i) - Real.log (p i)) ≤ q i - p i := by
    intro i
    have hlog : Real.log (q i / p i) ≤ q i / p i - 1 :=
      Real.log_le_sub_one_of_pos (div_pos (hq i) (hp i))
    have hrw : Real.log (q i / p i) = Real.log (q i) - Real.log (p i) :=
      Real.log_div (ne_of_gt (hq i)) (ne_of_gt (hp i))
    rw [hrw] at hlog
    have := mul_le_mul_of_nonneg_left hlog (hp i).le
    calc p i * (Real.log (q i) - Real.log (p i))
        ≤ p i * (q i / p i - 1) := this
      _ = q i - p i := by
          rw [mul_sub, mul_one, mul_comm (p i), div_mul_cancel₀ _ (ne_of_gt (hp i))]
  have hsum : ∑ i, p i * (Real.log (q i) - Real.log (p i))
      ≤ ∑ i, (q i - p i) :=
    Finset.sum_le_sum fun i _ => key i
  rw [Finset.sum_sub_distrib, ← hpq, sub_self] at hsum
  have flip : ∑ i, p i * (Real.log (p i) - Real.log (q i))
      = -∑ i, p i * (Real.log (q i) - Real.log (p i)) := by
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring
  rw [flip]
  linarith

/-- **Family Klein inequality**: relative entropy is nonnegative between
    Gibbs states of the retained family.  The commuting family reduces it
    to the classical Gibbs inequality; no operator convexity, no
    Kubo–Mori. -/
theorem relEntropyM_gibbs_nonneg (lam lam' : Fin 2 → ℝ) :
    0 ≤ relEntropyM (gibbsM SFam lam) (gibbsM SFam lam') := by
  rw [relEntropyM_gibbs_eq]
  exact finset_klDiv_nonneg _ _ (gibbsProb_pos lam) (gibbsProb_pos lam')
    (by rw [sum_gibbsProb, sum_gibbsProb])

end StateSide

/-! ## S4 — admissible channels, the closure defect, and the T0 theorem

`AdmissibleChannel` is **witness-free** (anti-smuggling, report §3 trap
i): linearity, positivity, trace preservation — the requirements
`rem:msascope` names for the admissible coarse-graining class — and
nothing else.  It never mentions `uuC`, `XXC`, the flux sector, or the
collar clause. -/

open scoped ComplexOrder

/-- An admissible coarse-graining channel in the sense of the paper's
    named requirements: linear, positive, trace-preserving. -/
structure AdmissibleChannel where
  /-- The underlying linear map. -/
  Φ : CollarC →ₗ[ℂ] CollarC
  /-- Positivity. -/
  pos : ∀ m : CollarC, m.PosSemidef → (Φ m).PosSemidef
  /-- Trace preservation. -/
  tracePreserving : ∀ m : CollarC, (Φ m).trace = m.trace

/-- Admissible channels map states to states. -/
noncomputable def AdmissibleChannel.applyD (C : AdmissibleChannel)
    (ρ : DensityMatrix) : DensityMatrix where
  ρ := C.Φ ρ.ρ
  posSemidef := C.pos _ ρ.posSemidef
  trace_one := by rw [C.tracePreserving, ρ.trace_one]

/-- The identity channel is admissible under the paper's current
    requirements — this is the "choose any channels" free-choice
    stipulation of `ax:maxent`'s internal refinement notion, at the state
    level. -/
noncomputable def idChannel : AdmissibleChannel where
  Φ := LinearMap.id
  pos := fun _ h => h
  tracePreserving := fun _ => rfl

@[simp] theorem idChannel_apply (m : CollarC) : idChannel.Φ m = m := rfl

/-- The closure defect of `def:closure-defect` (same-scale form): the
    infimum of relative entropy from the coarse-grained realized state to
    the exponential family of the retained list. -/
noncomputable def closureDefect (S : Fin 2 → CollarC) (C : AdmissibleChannel)
    (lam : Fin 2 → ℝ) : ℝ :=
  ⨅ lam' : Fin 2 → ℝ, relEntropyM (C.Φ (gibbsM S lam)) (gibbsM S lam')

/-- The defect is bounded by the relative entropy to any family member
    (the I-projection is an infimum).  The boundedness hypothesis is the
    honest guard against `Real.iInf` junk values; for the witness family
    it is discharged by the family Klein inequality. -/
theorem closureDefect_le_relEntropy_apply {S : Fin 2 → CollarC}
    {C : AdmissibleChannel} {lam : Fin 2 → ℝ}
    (hbdd : BddBelow (Set.range fun lam' =>
      relEntropyM (C.Φ (gibbsM S lam)) (gibbsM S lam')))
    (lam' : Fin 2 → ℝ) :
    closureDefect S C lam
      ≤ relEntropyM (C.Φ (gibbsM S lam)) (gibbsM S lam') :=
  ciInf_le hbdd lam'

/-- **Identity-channel discharge**: on the non-central witness family,
    the closure defect of the identity channel vanishes identically along
    the whole realized branch — with the genuine state-side objects, not
    their algebraic shadows.  Lower bound: family Klein inequality; upper
    bound: `relEntropyM_self` at the minimizer `lam' = lam`. -/
theorem identityChannel_closureDefect_eq_zero (lam : Fin 2 → ℝ) :
    closureDefect SFam idChannel lam = 0 := by
  have hbdd : BddBelow (Set.range fun lam' =>
      relEntropyM (idChannel.Φ (gibbsM SFam lam)) (gibbsM SFam lam')) := by
    refine ⟨0, ?_⟩
    rintro x ⟨lam', rfl⟩
    simp only [idChannel_apply]
    exact relEntropyM_gibbs_nonneg lam lam'
  refine le_antisymm ?_ ?_
  · have h := closureDefect_le_relEntropy_apply hbdd lam
    simp only [idChannel_apply] at h
    rwa [relEntropyM_self] at h
  · refine le_ciInf fun lam' => ?_
    simp only [idChannel_apply]
    exact relEntropyM_gibbs_nonneg lam lam'

/-- **T0 — the state-side identity-channel no-go.**

Under the paper's current admissible-channel requirements (positivity,
trace preservation, free channel choice), there is a retained family of
gauge-invariant densities containing a genuinely cross-cut, non-central
coupling, together with an admissible channel under which the family's
realized MaxEnt branch is refinement-stable with closure defect exactly
zero at every point of the branch.

So the current axioms cannot *force* the central-interface collar clause:
the non-central family `{uuC, XXC}` admits a refinement-stable realized
MaxEnt branch.  This sharpens #544 — forcing requires strengthening the
admissible class at the axiom level (T2; out of scope here) — and it does
NOT close it. -/
theorem stateSide_currentAxioms_cannot_force :
    ∃ S : Fin 2 → CollarC,
      (∀ a, GaugeInvariantC (S a)) ∧
      (∃ a, CrossCutC (S a) ∧ S a ∉ FluxC) ∧
      ∃ C : AdmissibleChannel, ∀ lam : Fin 2 → ℝ,
        closureDefect S C lam = 0 := by
  refine ⟨SFam, ?_, ⟨1, ?_, ?_⟩, idChannel,
    identityChannel_closureDefect_eq_zero⟩
  · intro a
    fin_cases a
    · exact uuC_invariant
    · exact XXC_invariant
  · exact XXC_crossCut
  · exact XXC_notMem_fluxC

/-! ## Axiom audit

Expected footprint for every entry: `[propext, Classical.choice,
Quot.sound]`.  No `sorry`, no `native_decide`, no extra axioms. -/

#print axioms uuC_sq
#print axioms XXC_crossCut
#print axioms XXC_notMem_fluxC
#print axioms pPlusC_idem
#print axioms specProjP_mul
#print axioms sum_specProjP
#print axioms exp_sum_smul_specProjP
#print axioms gibbsDM
#print axioms matLog_gibbsM
#print axioms relEntropyM_self
#print axioms relEntropyM_gibbs_eq
#print axioms finset_klDiv_nonneg
#print axioms relEntropyM_gibbs_nonneg
#print axioms identityChannel_closureDefect_eq_zero
#print axioms stateSide_currentAxioms_cannot_force

end OPH
