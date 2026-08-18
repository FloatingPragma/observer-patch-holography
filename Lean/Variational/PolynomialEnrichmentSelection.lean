import Variational.EnrichmentCharacterization

set_option autoImplicit false

namespace OPH.Variational

/-!
# Boost-covariance selection in the full polynomial enrichment grammar
(V3, issue #739, register row PR-06)

`EnrichmentCharacterization` classifies corner-invisible enrichments
inside the quadratic two-slot grammar and proves that the quartic
continuation `lam * y^2 * (y-1)^2` preserves corner invisibility, the
committed stationarity, the committed single-site minimum, a positive
velocity Hessian, and momentum invertibility while changing the
Lagrangian, so the committed selector is not proved minimal in any
grammar larger than the quadratics.  This module extends the ambient
family from quadratics to two-slot polynomials of every degree bound
`n`, states one further named clause, and tests it against that
grammar.  All statements are exact real algebra plus the one polynomial
identity package named below; no premise is discharged.

The committed null notion.  A two-point sentence is a discrete null
Lagrangian when every single-site interior replacement leaves the sum
of the two adjacent terms unchanged (`NullLagrangian`).  By the
committed difference identity `localAction_update_diff`, such a
sentence contributes only telescoping and boundary terms to every
committed fixed-endpoint local action
(`nullLagrangian_localAction_invariant`).  Telescoping sentences
`F y - F x + k` are null (`telescoping_nullLagrangian`), and a sentence
of the velocity slot alone is null exactly when it is constant
(`nullLagrangian_velocitySentence_iff`).

The named clause.  The boost difference of an enrichment `q` by a
shift pair `(u, v)` is `q (x+u) (y+v) - q x y` (`boostDiff`).  The
literal first-order clause, every boost difference is null
(`BoostNull`), is unavailable as a selector for the committed
structure: the committed increment `(a*/2) y (y-1)` fails it at the
diagonal shift `(1, 1)` (`committed_increment_not_boostNull`), because
a lattice-anchored quadratic has a linear, non-telescoping boost
residue.  The named candidate clause is therefore second order: every
boost difference of a boost difference is null (`BoostNullCovariant`);
the boost residue of the enrichment is itself boost-invisible up to
telescoping and boundary terms.

What IS proved, all machine-checked:

* **Quadratic sufficiency (S1).**  Every sentence of the committed
  quadratic grammar satisfies the clause
  (`quadPoly_boostNullCovariant`), in particular the full committed
  two-parameter counterfamily (`twoSlotIncrement_boostNullCovariant`).
  The committed stationary line and the committed point pass.

* **Exact characterization on the velocity-only slice (S2).**  For
  every degree bound `n` and coefficient family `c`, the velocity-only
  polynomial enrichment `y ↦ ∑ c j * y^j` satisfies the clause exactly
  when every coefficient of degree three or higher vanishes
  (`velocityPoly_boostNull_iff`): among velocity-only polynomial
  enrichments of arbitrary degree, the clause holds if and only if the
  enrichment is quadratic-affine.

* **The audit quartic fails the clause (S3).**  The committed quartic
  increment fails the clause at the concrete shift pair
  `(0,1), (0,1)` and replacement pair `x = 1, x' = 0`
  (`quarticIncrement_not_boostNullCovariant`), and so does the full
  audit continuation enrichment
  (`auditQuarticEnrichment_not_boostNullCovariant`).  The bundled
  receipt `audit_quartic_passes_committed_clauses_fails_boost`
  restates the committed five preserved properties of the continuation
  next to the boost failure.

* **Composition (S4).**  In the full two-slot polynomial grammar at
  every degree bound, corner invisibility, the committed velocity-only
  clause, the committed stationarity calibration, and the boost clause
  select exactly the committed point: the enriched Lagrangian equals
  the committed mode-extremal rule pointwise, with the velocity-slice
  coefficients pinned (`polyGrammar_selection_eq_committed`).  The
  committed point itself satisfies all four clauses inside the grammar
  (`committed_point_satisfies_all_clauses`), so the axiom system is
  nonvacuous.

* **Load-bearing receipts (S5).**  Dropping the boost clause readmits
  the audit quartic (S3 next to the committed continuation receipts).
  The boost clause alone does not fix the committed point: the whole
  committed two-parameter family passes it, and the corner-invisible
  antisymmetric cubic `(x-y)^3 - (x-y)` passes it together with corner
  invisibility and the committed stationarity calibration while being
  velocity-dependent on the first slot and outside the quadratic
  grammar (`boost_clause_does_not_force_quadraticity`,
  `velocityOnly_load_bearing_in_polynomial_grammar`).  Stationarity
  and corner invisibility are each load-bearing as well
  (`stationarity_clause_load_bearing`, `corner_clause_load_bearing`).

Delimitation of the sketched two-slot characterization.  In the full
two-slot grammar the clause does NOT force quadraticity: the
antisymmetric cubic is a machine-checked counterexample
(`antisymCubic_boostNullCovariant`, `antisymCubic_ne_quadGrammar`).
The exact quadratic-affine characterization holds on the velocity-only
slice (S2), which is where the committed selection lives; the
committed velocity-only clause is exactly what removes the cubic
direction, so `VelocityOnly` stays load-bearing in the enlarged
grammar.

What is NOT proved.  The grammar is polynomial only; smooth
non-polynomial enrichments are outside every theorem here.  The boost
clause is a named candidate clause for register row PR-06, not a
discharge of PR-06 and not a physically derived symmetry: no Galilei
group, no boost action on realized histories, and no invariance of the
committed source law is constructed, and the physical status of the
clause is a register-side question.  The committed velocity-only
clause results of `EnrichmentCharacterization` are unchanged and still
required.  The realized-history no-go of `RealizedHistoryLegendreNoGo`
stays in force.  The full two-slot kernel of the clause is not
classified here beyond the exhibited members: the quadratics and the
antisymmetric cubic direction pass, every velocity-only member is
quadratic-affine, and the audit quartic fails.  Derivative data is
register row PR-45 through the committed packets.
-/

/-! ## The committed null-Lagrangian notion -/

/-- A discrete null Lagrangian in the committed fixed-endpoint sense:
replacing the interior record of two adjacent segments never changes
the sum of the two adjacent terms.  By `localAction_update_diff` this
is exactly the class of sentences that contribute only telescoping and
boundary terms to every committed fixed-endpoint local action. -/
def NullLagrangian (N : ℝ → ℝ → ℝ) : Prop :=
  ∀ a b x x' : ℝ, N a x + N x b = N a x' + N x' b

/-- **Consistency with the committed action.**  A null sentence
contributes the same amount to the committed local action on every
single-site interior variation: its contribution is fixed by the
boundary data. -/
theorem nullLagrangian_localAction_invariant (N : ℝ → ℝ → ℝ)
    (h : NullLagrangian N) {M : ℕ} (γ : Fin (M + 1) → ℝ)
    {k m : Fin M} (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction N (Function.update γ k.succ x) = localAction N γ := by
  have hdiff := localAction_update_diff N γ hkm x
  have hmc : γ m.castSucc = γ k.succ := by rw [← hkm]
  rw [hmc] at hdiff
  have hkey := h (γ k.castSucc) (γ m.succ) x (γ k.succ)
  linarith

/-- Telescoping sentences plus per-step constants are null. -/
theorem telescoping_nullLagrangian (F : ℝ → ℝ) (k : ℝ) :
    NullLagrangian (fun x y => F y - F x + k) := by
  intro a b x x'
  dsimp only
  ring

/-- Constant sentences are null. -/
theorem constant_nullLagrangian (k : ℝ) :
    NullLagrangian (fun _ _ => k) :=
  fun _ _ _ _ => rfl

/-- Null sentences are closed under addition. -/
theorem NullLagrangian.add {N N' : ℝ → ℝ → ℝ} (h : NullLagrangian N)
    (h' : NullLagrangian N') :
    NullLagrangian (fun x y => N x y + N' x y) := by
  intro a b x x'
  have h1 := h a b x x'
  have h2 := h' a b x x'
  dsimp only
  linarith

/-- Null sentences are closed under subtraction. -/
theorem NullLagrangian.sub {N N' : ℝ → ℝ → ℝ} (h : NullLagrangian N)
    (h' : NullLagrangian N') :
    NullLagrangian (fun x y => N x y - N' x y) := by
  intro a b x x'
  have h1 := h a b x x'
  have h2 := h' a b x x'
  dsimp only
  linarith

/-- Nullness transfers along pointwise equality. -/
theorem NullLagrangian.congr {N N' : ℝ → ℝ → ℝ}
    (hEq : ∀ x y, N x y = N' x y) (h : NullLagrangian N) :
    NullLagrangian N' := by
  intro a b x x'
  rw [← hEq a x, ← hEq x b, ← hEq a x', ← hEq x' b]
  exact h a b x x'

/-- A sentence of the velocity slot alone is null exactly when it is
constant: pure velocity content never telescopes. -/
theorem nullLagrangian_velocitySentence_iff (g : ℝ → ℝ) :
    NullLagrangian (fun _ y => g y) ↔ ∀ y : ℝ, g y = g 0 := by
  constructor
  · intro h y
    have hkey := h 0 0 y 0
    dsimp only at hkey
    linarith
  · intro h a b x x'
    dsimp only
    rw [h x, h x']

/-! ## Boost differences and the named clause -/

/-- The boost difference of an enrichment by the shift pair `(u, v)`:
the residue of the enrichment under a rigid shift of the two record
slots. -/
noncomputable def boostDiff (u v : ℝ) (q : ℝ → ℝ → ℝ) : ℝ → ℝ → ℝ :=
  fun x y => q (x + u) (y + v) - q x y

/-- Mutation guard: the first-slot shift component of `boostDiff` is
load-bearing.  For the record-slot square `x ^ 2`, the shift pair
`(1, 0)` produces residue one at the origin; a drifted `boostDiff`
that ignores the first shift component evaluates to zero here. -/
theorem boostDiff_first_slot_guard :
    boostDiff 1 0 (fun x _ => x ^ 2) 0 0 = 1 := by
  norm_num [boostDiff]

/-- The literal first-order boost clause: every boost difference is
null.  Recorded for the delimitation receipt
`committed_increment_not_boostNull`; it is not the named candidate
clause. -/
def BoostNull (q : ℝ → ℝ → ℝ) : Prop :=
  ∀ u v : ℝ, NullLagrangian (boostDiff u v q)

/-- **The named clause: boost-null covariance.**  Every boost
difference of a boost difference of the enrichment is a discrete null
Lagrangian: the boost residue of the enrichment is itself
boost-invisible up to telescoping and boundary terms.  This is the
candidate selection clause proposed for register row PR-06; it is a
declared axiom of a named axiomatization, not a source-produced fact
and not a constructed symmetry group. -/
def BoostNullCovariant (q : ℝ → ℝ → ℝ) : Prop :=
  ∀ u₁ v₁ u₂ v₂ : ℝ,
    NullLagrangian (boostDiff u₁ v₁ (boostDiff u₂ v₂ q))

/-- Boost-null covariance transfers along pointwise equality. -/
theorem BoostNullCovariant.congr {q q' : ℝ → ℝ → ℝ}
    (hEq : ∀ x y, q x y = q' x y) (h : BoostNullCovariant q) :
    BoostNullCovariant q' := by
  intro u₁ v₁ u₂ v₂
  refine NullLagrangian.congr (fun x y => ?_) (h u₁ v₁ u₂ v₂)
  simp only [boostDiff]
  rw [hEq, hEq, hEq, hEq]

/-- Boost-null covariance is closed under subtraction. -/
theorem BoostNullCovariant.sub {q r : ℝ → ℝ → ℝ}
    (hq : BoostNullCovariant q) (hr : BoostNullCovariant r) :
    BoostNullCovariant (fun x y => q x y - r x y) := by
  intro u₁ v₁ u₂ v₂
  refine NullLagrangian.congr (fun x y => ?_)
    ((hq u₁ v₁ u₂ v₂).sub (hr u₁ v₁ u₂ v₂))
  simp only [boostDiff]
  ring

/-- **The first-order boost clause is unavailable (delimitation).**
The committed increment `(a*/2) y (y-1)` fails the literal clause
`BoostNull` at the diagonal shift `(1, 1)`: its boost residue is
`a* y` plus a constant, a pure velocity sentence that is not constant,
hence not null.  A lattice-anchored quadratic can satisfy at most the
second-order clause, which is why the named candidate clause is
`BoostNullCovariant`. -/
theorem committed_increment_not_boostNull :
    ¬ BoostNull (twoSlotIncrement modeExtremalCurvature 0) := by
  intro h
  have hkey := h 1 1 0 0 1 0
  simp only [boostDiff, twoSlotIncrement] at hkey
  norm_num at hkey
  exact absurd hkey (ne_of_gt modeExtremalCurvature_pos)

/-! ## The full polynomial grammar -/

/-- A velocity-only polynomial enrichment sentence with degree bound
`n`: the value `∑ j ≤ n, c j * y ^ j` read off the velocity slot. -/
noncomputable def velocityPoly (n : ℕ) (c : ℕ → ℝ) (y : ℝ) : ℝ :=
  ∑ j ∈ Finset.range (n + 1), c j * y ^ j

/-- The velocity derivative of `velocityPoly`, term by term. -/
noncomputable def velocityPolyDeriv (n : ℕ) (c : ℕ → ℝ) (y : ℝ) : ℝ :=
  ∑ j ∈ Finset.range (n + 1), c j * ((j : ℝ) * y ^ (j - 1))

/-- The derivative certificate of the velocity-only polynomial
grammar: `velocityPolyDeriv` is the exact derivative of
`velocityPoly` at every point. -/
theorem velocityPoly_hasDerivAt (n : ℕ) (c : ℕ → ℝ) (y : ℝ) :
    HasDerivAt (velocityPoly n c) (velocityPolyDeriv n c y) y := by
  have h := HasDerivAt.sum
    (fun j (_ : j ∈ Finset.range (n + 1)) =>
      (hasDerivAt_pow j y).const_mul (c j))
  have hfun : (∑ j ∈ Finset.range (n + 1), fun t : ℝ => c j * t ^ j)
      = velocityPoly n c := by
    funext t
    simp [velocityPoly]
  unfold velocityPolyDeriv
  rw [← hfun]
  exact h

/-- A two-slot polynomial enrichment sentence with degree bound `n` in
each slot: the ambient family of this module.  The committed quadratic
grammar is the `n = 2` fragment with the mixed monomials `x^2 y`,
`x y^2`, `x^2 y^2` removed; every `quadPoly` sentence is a member. -/
noncomputable def polyEnrichment (n : ℕ) (c : ℕ → ℕ → ℝ)
    (x y : ℝ) : ℝ :=
  ∑ i ∈ Finset.range (n + 1), ∑ j ∈ Finset.range (n + 1),
    c i j * (x ^ i * y ^ j)

/-- On the `x = 0` slice, a two-slot polynomial enrichment reads its
first coefficient row as a velocity-only polynomial. -/
theorem polyEnrichment_velocity_slice (n : ℕ) (c : ℕ → ℕ → ℝ)
    (y : ℝ) :
    polyEnrichment n c 0 y = velocityPoly n (fun j => c 0 j) y := by
  unfold polyEnrichment velocityPoly
  rw [Finset.sum_eq_single 0]
  · simp
  · intro i _ hi
    simp [zero_pow hi]
  · intro hmem
    exact absurd (Finset.mem_range.mpr (Nat.succ_pos n)) hmem

/-- Quadratic normal form of a velocity-only polynomial whose
coefficients of degree three and higher vanish, with the matching
derivative value at the committed junction. -/
theorem velocityPoly_quadratic_normal_form (n : ℕ) (c : ℕ → ℝ)
    (h : ∀ j, 3 ≤ j → j ≤ n → c j = 0) :
    ∃ B D : ℝ,
      (∀ y : ℝ, velocityPoly n c y = c 0 + B * y + D * (y * y))
        ∧ velocityPolyDeriv n c 1 = B + 2 * D
        ∧ (2 ≤ n → B = c 1 ∧ D = c 2) := by
  cases n with
  | zero =>
    refine ⟨0, 0, ?_, ?_, by omega⟩
    · intro y
      simp [velocityPoly]
    · simp [velocityPolyDeriv]
  | succ n1 =>
    cases n1 with
    | zero =>
      refine ⟨c 1, 0, ?_, ?_, by omega⟩
      · intro y
        simp only [velocityPoly, Finset.sum_range_succ,
          Finset.sum_range_zero]
        ring
      · simp only [velocityPolyDeriv, Finset.sum_range_succ,
          Finset.sum_range_zero]
        norm_num
    | succ m =>
      have hsub : ∀ (f : ℕ → ℝ),
          (∀ j, 3 ≤ j → j < m + 3 → f j = 0) →
          ∑ j ∈ Finset.range (m + 3), f j = ∑ j ∈ Finset.range 3, f j := by
        intro f hf
        symm
        refine Finset.sum_subset ?_ ?_
        · intro j hj
          simp only [Finset.mem_range] at hj ⊢
          omega
        · intro j hj hj3
          simp only [Finset.mem_range] at hj hj3
          exact hf j (by omega) hj
      refine ⟨c 1, c 2, ?_, ?_, fun _ => ⟨rfl, rfl⟩⟩
      · intro y
        have hvan : ∀ j, 3 ≤ j → j < m + 3 → c j * y ^ j = 0 :=
          fun j h3 hj => by rw [h j h3 (by omega), zero_mul]
        have hred := hsub (fun j => c j * y ^ j) hvan
        unfold velocityPoly
        rw [show m + 1 + 1 + 1 = m + 3 from rfl, hred]
        simp only [Finset.sum_range_succ, Finset.sum_range_zero]
        ring
      · have hvan : ∀ j, 3 ≤ j → j < m + 3 →
            c j * ((j : ℝ) * (1 : ℝ) ^ (j - 1)) = 0 :=
          fun j h3 hj => by rw [h j h3 (by omega), zero_mul]
        have hred := hsub (fun j => c j * ((j : ℝ) * (1 : ℝ) ^ (j - 1)))
          hvan
        unfold velocityPolyDeriv
        rw [show m + 1 + 1 + 1 = m + 3 from rfl, hred]
        simp only [Finset.sum_range_succ, Finset.sum_range_zero]
        norm_num
        ring

/-! ## Quadratic sufficiency: the committed grammar passes -/

/-- Every sentence with a quadratic normal form satisfies the boost
clause: its second boost difference is a constant, and constants are
null. -/
theorem boostNullCovariant_of_quadratic (q : ℝ → ℝ → ℝ)
    (A B G E F D : ℝ)
    (h : ∀ x y, q x y
      = A + B * x + G * y + E * (x * y) + F * (x * x) + D * (y * y)) :
    BoostNullCovariant q := by
  intro u₁ v₁ u₂ v₂ a b x x'
  simp only [boostDiff, h]
  ring

/-- **Quadratic sufficiency (S1).**  Every sentence of the committed
quadratic grammar satisfies the boost clause. -/
theorem quadPoly_boostNullCovariant (α β γ δ ε ζ : ℝ) :
    BoostNullCovariant (quadPoly α β γ δ ε ζ) :=
  boostNullCovariant_of_quadratic _ α β γ δ ε ζ (fun _ _ => rfl)

/-- The full committed two-parameter counterfamily satisfies the boost
clause; in particular the committed stationary line and the committed
point pass. -/
theorem twoSlotIncrement_boostNullCovariant (a c : ℝ) :
    BoostNullCovariant (twoSlotIncrement a c) := by
  refine boostNullCovariant_of_quadratic _ 0 (-(c / 2)) (-(a / 2)) 0
    (c / 2) (a / 2) ?_
  intro x y
  unfold twoSlotIncrement
  ring

/-! ## Exact characterization on the velocity-only slice -/

/-- Coefficient extraction: the boost clause forces every velocity
coefficient of degree three or higher to vanish.  The proof transfers
the constancy of the second boost difference to the coefficient
polynomial and reads off a Hasse-derivative coefficient. -/
theorem velocityPoly_boostNull_highCoeff_vanish (n : ℕ) (c : ℕ → ℝ)
    (h : BoostNullCovariant (fun _ y => velocityPoly n c y)) :
    ∀ j, 3 ≤ j → j ≤ n → c j = 0 := by
  intro j hj3 hjn
  by_contra hcj
  -- the coefficient polynomial
  set P : Polynomial ℝ :=
    ∑ i ∈ Finset.range (n + 1), Polynomial.C (c i) * Polynomial.X ^ i
    with hP
  have heval : ∀ y : ℝ, P.eval y = velocityPoly n c y := by
    intro y
    simp [hP, Polynomial.eval_finset_sum, velocityPoly]
  have hcoeff : ∀ i, i ≤ n → P.coeff i = c i := by
    intro i hi
    rw [hP, Polynomial.finset_sum_coeff]
    simp only [Polynomial.coeff_C_mul, Polynomial.coeff_X_pow,
      mul_ite, mul_one, mul_zero]
    rw [Finset.sum_ite_eq (Finset.range (n + 1)) i c]
    simp [Nat.lt_succ_of_le hi]
  -- constancy of the second boost difference
  have hgconst : ∀ y : ℝ,
      velocityPoly n c (y + 2) - 2 * velocityPoly n c (y + 1)
          + velocityPoly n c y
        = velocityPoly n c 2 - 2 * velocityPoly n c 1
          + velocityPoly n c 0 := by
    intro y
    have hkey := h 0 1 0 1 0 0 y 0
    simp only [boostDiff] at hkey
    ring_nf at hkey ⊢
    linarith
  -- the difference polynomial
  set R : Polynomial ℝ :=
    Polynomial.taylor 2 P - Polynomial.C 2 * Polynomial.taylor 1 P + P
    with hR
  have hReval : ∀ y : ℝ,
      R.eval y = velocityPoly n c (y + 2) - 2 * velocityPoly n c (y + 1)
        + velocityPoly n c y := by
    intro y
    simp [hR, Polynomial.taylor_eval, heval]
  have hRconst : R = Polynomial.C (R.eval 0) := by
    apply Polynomial.funext
    intro y
    rw [Polynomial.eval_C, hReval y, hReval 0, hgconst y]
    norm_num
  -- degree bookkeeping
  have hjd : j ≤ P.natDegree := by
    apply Polynomial.le_natDegree_of_ne_zero
    rw [hcoeff j hjn]
    exact hcj
  set d : ℕ := P.natDegree with hd
  have hd3 : 3 ≤ d := le_trans hj3 hjd
  set k : ℕ := d - 2 with hk
  have hk1 : 1 ≤ k := by omega
  have h2k : 2 + k = d := by omega
  -- the top surviving coefficient of R
  set H : Polynomial ℝ := Polynomial.hasseDeriv k P with hH
  have hHdeg : H.natDegree < 3 := by
    have hle := Polynomial.natDegree_hasseDeriv_le P k
    rw [← hH] at hle
    omega
  have hHeval : ∀ t : ℝ,
      H.eval t = H.coeff 0 + H.coeff 1 * t + H.coeff 2 * (t * t) := by
    intro t
    rw [Polynomial.eval_eq_sum_range' hHdeg t]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero]
    ring
  have hRk : R.coeff k = 2 * H.coeff 2 := by
    rw [hR]
    simp only [Polynomial.coeff_add, Polynomial.coeff_sub,
      Polynomial.coeff_C_mul, Polynomial.taylor_coeff, ← hH]
    have hP0 : P.coeff k = H.eval 0 := by
      rw [hH, ← Polynomial.coeff_zero_eq_eval_zero,
        Polynomial.hasseDeriv_coeff]
      simp
    rw [hP0, hHeval 2, hHeval 1, hHeval 0]
    ring
  have hHtop : H.coeff 2 = (d.choose k : ℝ) * P.coeff d := by
    rw [hH, Polynomial.hasseDeriv_coeff, h2k]
  -- contradiction: R is constant but has a nonzero positive-degree coefficient
  have hPd : P.coeff d ≠ 0 := by
    rw [hd]
    rw [← Polynomial.leadingCoeff]
    rw [Polynomial.leadingCoeff_ne_zero]
    intro hP0
    apply hcj
    rw [← hcoeff j hjn, hP0]
    simp
  have hchoose : (d.choose k : ℝ) ≠ 0 :=
    Nat.cast_ne_zero.mpr (Nat.choose_pos (show k ≤ d by omega)).ne'
  have hzero : R.coeff k = 0 := by
    rw [hRconst, Polynomial.coeff_C]
    simp [show k ≠ 0 by omega]
  rw [hRk, hHtop] at hzero
  exact mul_ne_zero two_ne_zero (mul_ne_zero hchoose hPd) hzero

/-- The quadratic-affine members pass: sufficiency on the velocity-only
slice at every degree bound. -/
theorem velocityPoly_boostNull_of_lowDegree (n : ℕ) (c : ℕ → ℝ)
    (h : ∀ j, 3 ≤ j → j ≤ n → c j = 0) :
    BoostNullCovariant (fun _ y => velocityPoly n c y) := by
  obtain ⟨B, D, hBD, -, -⟩ := velocityPoly_quadratic_normal_form n c h
  refine boostNullCovariant_of_quadratic _ (c 0) 0 B 0 0 D ?_
  intro x y
  rw [hBD y]
  ring

/-- **Exact characterization on the velocity-only slice (S2).**  Among
velocity-only polynomial enrichments of every degree bound `n`, the
boost clause holds if and only if the enrichment is quadratic-affine:
every coefficient of degree three or higher vanishes. -/
theorem velocityPoly_boostNull_iff (n : ℕ) (c : ℕ → ℝ) :
    BoostNullCovariant (fun _ y => velocityPoly n c y)
      ↔ ∀ j, 3 ≤ j → j ≤ n → c j = 0 :=
  ⟨velocityPoly_boostNull_highCoeff_vanish n c,
    velocityPoly_boostNull_of_lowDegree n c⟩

/-! ## The audit quartic fails the clause -/

/-- **The committed quartic increment fails the boost clause (S3).**
The exact witness: at the shift pairs `(0,1), (0,1)`, endpoints
`a = b = 0`, and replacement pair `x = 1, x' = 0`, the null equation
fails by `24 * lam`. -/
theorem quarticIncrement_not_boostNullCovariant (lam : ℝ)
    (hlam : lam ≠ 0) :
    ¬ BoostNullCovariant (fun _ y => quarticIncrement lam y) := by
  intro h
  have hkey := h 0 1 0 1 0 0 1 0
  simp only [boostDiff, quarticIncrement] at hkey
  norm_num at hkey
  exact hlam (by linarith)

/-- The audit continuation as an enrichment over the committed
bilinear base: the committed increment plus the committed quartic
increment. -/
noncomputable def auditQuarticEnrichment (lam : ℝ) : ℝ → ℝ → ℝ :=
  fun x y => twoSlotIncrement modeExtremalCurvature 0 x y
    + quarticIncrement lam y

/-- The audit enrichment reproduces the committed quartic continuation
Lagrangian exactly. -/
theorem auditQuartic_enrichment_eq (lam x y : ℝ) :
    chainLogLagrangian x y + auditQuarticEnrichment lam x y
      = chainQuarticLagrangian lam x y := by
  unfold auditQuarticEnrichment chainQuarticLagrangian
    chainCurvedLagrangian twoSlotIncrement
  ring

/-- The velocity-slice coefficients of the audit enrichment: a degree
bound `4` member of the grammar. -/
noncomputable def auditQuarticCoeffs (lam : ℝ) : ℕ → ℝ := fun j =>
  if j = 1 then -(modeExtremalCurvature / 2)
  else if j = 2 then modeExtremalCurvature / 2 + lam
  else if j = 3 then -(2 * lam)
  else if j = 4 then lam
  else 0

/-- The audit enrichment is a member of the polynomial grammar. -/
theorem auditQuartic_in_grammar (lam x y : ℝ) :
    auditQuarticEnrichment lam x y
      = velocityPoly 4 (auditQuarticCoeffs lam) y := by
  unfold auditQuarticEnrichment twoSlotIncrement quarticIncrement
    velocityPoly auditQuarticCoeffs
  simp [Finset.sum_range_succ]
  ring

/-- The full audit enrichment fails the boost clause: the committed
increment passes, so the failure of the quartic part is inherited. -/
theorem auditQuarticEnrichment_not_boostNullCovariant (lam : ℝ)
    (hlam : lam ≠ 0) :
    ¬ BoostNullCovariant (auditQuarticEnrichment lam) := by
  intro h
  have hq : BoostNullCovariant
      (fun x y => auditQuarticEnrichment lam x y
        - twoSlotIncrement modeExtremalCurvature 0 x y) :=
    h.sub (twoSlotIncrement_boostNullCovariant _ _)
  have hquartic : BoostNullCovariant
      (fun _ y => quarticIncrement lam y) := by
    refine BoostNullCovariant.congr (fun x y => ?_) hq
    unfold auditQuarticEnrichment
    ring
  exact quarticIncrement_not_boostNullCovariant lam hlam hquartic

/-- **The audit continuation passes every committed clause and fails
the boost clause (S3, bundled).**  At the audit parameter `a*/2`, the
continuation keeps the committed corner table, the velocity-only
clause, the committed stationarity, the committed global real
single-site minimum, an everywhere positive velocity Hessian, and a
velocity solver for every momentum, all by the committed receipts of
`EnrichmentCharacterization`; and its enrichment fails the boost
clause while the committed increment satisfies it.  Dropping the boost
clause therefore readmits the quartic: the clause is load-bearing. -/
theorem audit_quartic_passes_committed_clauses_fails_boost :
    ∃ lam : ℝ,
      0 < lam ∧ lam < modeExtremalCurvature
        ∧ CornerInvisible (fun _ y => quarticIncrement lam y)
        ∧ VelocityOnly (fun _ y => quarticIncrement lam y)
        ∧ chainQuarticD2 lam 1 1 + chainCurvedD1 1 1 = 0
        ∧ (∀ M : ℕ, ∀ {k m : Fin M}, k.succ = m.castSucc →
            ∀ x : ℝ,
              localAction (chainQuarticLagrangian lam)
                  (chainEmb (constOneHistory M))
                ≤ localAction (chainQuarticLagrangian lam)
                    (Function.update (chainEmb (constOneHistory M))
                      k.succ x))
        ∧ (∀ y : ℝ,
            0 < modeExtremalCurvature
              + quarticIncrementSecondDeriv lam y)
        ∧ (∃ vel : ℝ → ℝ → ℝ,
            SolvesMomentum (chainQuarticLagrangian lam) vel)
        ∧ ¬ BoostNullCovariant (auditQuarticEnrichment lam)
        ∧ BoostNullCovariant
            (twoSlotIncrement modeExtremalCurvature 0) := by
  have hpos := modeExtremalCurvature_pos
  have hlam0 : 0 < modeExtremalCurvature / 2 := half_pos hpos
  have hlam : modeExtremalCurvature / 2 < modeExtremalCurvature :=
    half_lt_self hpos
  refine ⟨modeExtremalCurvature / 2, hlam0, hlam,
    quarticIncrement_cornerInvisible _,
    quarticIncrement_velocityOnly _,
    chainQuartic_stationary _,
    ?_, ?_,
    ⟨chainQuarticVelocitySolver _ hlam0 hlam,
      chainQuartic_velocitySolver_solves _ hlam0 hlam⟩,
    auditQuarticEnrichment_not_boostNullCovariant _ (ne_of_gt hlam0),
    twoSlotIncrement_boostNullCovariant _ _⟩
  · intro M k m hkm x
    exact chainQuartic_constOne_realMin _ hlam0.le M hkm x
  · intro y
    exact chainQuartic_hessian_pos _ y hlam0 hlam

/-! ## The antisymmetric cubic: the clause does not force
quadraticity in the full two-slot grammar -/

/-- The antisymmetric corner-invisible cubic `(x-y)^3 - (x-y)`. -/
noncomputable def antisymCubic (x y : ℝ) : ℝ :=
  (x - y) ^ 3 - (x - y)

/-- The antisymmetric cubic vanishes at all four binary corners. -/
theorem antisymCubic_cornerInvisible : CornerInvisible antisymCubic := by
  unfold CornerInvisible antisymCubic
  norm_num

/-- The antisymmetric cubic satisfies the boost clause: its second
boost difference is affine with opposite slot slopes, which
telescopes. -/
theorem antisymCubic_boostNullCovariant :
    BoostNullCovariant antisymCubic := by
  intro u₁ v₁ u₂ v₂ a b x x'
  simp only [boostDiff, antisymCubic]
  ring

/-- The antisymmetric cubic is not any sentence of the committed
quadratic grammar: four first-slot evaluations separate it from every
quadratic. -/
theorem antisymCubic_ne_quadGrammar (α β γ δ ε ζ : ℝ) :
    antisymCubic ≠ quadPoly α β γ δ ε ζ := by
  intro hEq
  have h0 := congrFun (congrFun hEq 0) 0
  have h1 := congrFun (congrFun hEq 1) 0
  have h2 := congrFun (congrFun hEq 2) 0
  have h3 := congrFun (congrFun hEq 3) 0
  simp only [antisymCubic, quadPoly] at h0 h1 h2 h3
  norm_num at h0 h1 h2 h3
  linarith

/-- The two-slot coefficient family of the antisymmetric cubic. -/
noncomputable def antisymCoeffs : ℕ → ℕ → ℝ := fun i j =>
  if i = 3 ∧ j = 0 then 1
  else if i = 2 ∧ j = 1 then -3
  else if i = 1 ∧ j = 2 then 3
  else if i = 0 ∧ j = 3 then -1
  else if i = 1 ∧ j = 0 then -1
  else if i = 0 ∧ j = 1 then 1
  else 0

/-- The antisymmetric cubic is a member of the polynomial grammar at
degree bound three: the delimitation below concerns the ambient
grammar itself. -/
theorem antisymCoeffs_realizes_cubic (x y : ℝ) :
    polyEnrichment 3 antisymCoeffs x y = antisymCubic x y := by
  unfold polyEnrichment antisymCoeffs antisymCubic
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  norm_num
  ring

/-- The antisymmetric cubic is velocity-dependent on the first slot:
it fails the committed velocity-only clause. -/
theorem antisymCubic_not_velocityOnly : ¬ VelocityOnly antisymCubic := by
  intro h
  have hkey := h 2 0 0
  unfold antisymCubic at hkey
  norm_num at hkey

/-- The first-slot derivative value of the antisymmetric cubic. -/
noncomputable def antisymCubicD1 (x y : ℝ) : ℝ :=
  3 * (x - y) ^ 2 - 1

/-- The second-slot derivative value of the antisymmetric cubic. -/
noncomputable def antisymCubicD2 (x y : ℝ) : ℝ :=
  -(3 * (x - y) ^ 2) + 1

/-- First-slot derivative certificate of the antisymmetric cubic. -/
theorem antisymCubic_hasDerivAt_slot1 (x y : ℝ) :
    HasDerivAt (fun t => antisymCubic t y) (antisymCubicD1 x y) x := by
  unfold antisymCubic antisymCubicD1
  have hb : HasDerivAt (fun t : ℝ => t - y) 1 x :=
    (hasDerivAt_id x).sub_const y
  have h := (hb.pow 3).sub hb
  convert h using 1
  push_cast
  ring

/-- Second-slot derivative certificate of the antisymmetric cubic. -/
theorem antisymCubic_hasDerivAt_slot2 (x y : ℝ) :
    HasDerivAt (fun t => antisymCubic x t) (antisymCubicD2 x y) y := by
  unfold antisymCubic antisymCubicD2
  have hb : HasDerivAt (fun t : ℝ => x - t) (-1) y :=
    (hasDerivAt_id y).const_sub x
  have h := (hb.pow 3).sub hb
  convert h using 1
  push_cast
  ring

/-- The antisymmetric cubic is neutral for the committed stationarity
calibration: its discrete Euler-Lagrange sum at the constant-one
junction vanishes at every parameter. -/
theorem antisymCubic_stationarity_neutral :
    antisymCubicD2 1 1 + antisymCubicD1 1 1 = 0 := by
  unfold antisymCubicD1 antisymCubicD2
  norm_num

/-- **Delimitation of the sketched two-slot characterization.**  In
the full two-slot grammar, the boost clause does not force
quadraticity: the antisymmetric cubic is corner-invisible, satisfies
the boost clause, is stationarity-neutral with certified slot
derivatives, and equals no sentence of the committed quadratic
grammar.  The exact quadratic-affine characterization of the boost
clause therefore holds on the velocity-only slice and fails in the
full two-slot grammar; the committed velocity-only clause is what
removes this direction. -/
theorem boost_clause_does_not_force_quadraticity :
    CornerInvisible antisymCubic
      ∧ BoostNullCovariant antisymCubic
      ∧ (antisymCubicD2 1 1 + antisymCubicD1 1 1 = 0)
      ∧ (∀ x y : ℝ,
          HasDerivAt (fun t => antisymCubic t y) (antisymCubicD1 x y) x
            ∧ HasDerivAt (fun t => antisymCubic x t)
                (antisymCubicD2 x y) y)
      ∧ ¬ VelocityOnly antisymCubic
      ∧ ∀ α β γ δ ε ζ : ℝ, antisymCubic ≠ quadPoly α β γ δ ε ζ :=
  ⟨antisymCubic_cornerInvisible, antisymCubic_boostNullCovariant,
    antisymCubic_stationarity_neutral,
    fun x y => ⟨antisymCubic_hasDerivAt_slot1 x y,
      antisymCubic_hasDerivAt_slot2 x y⟩,
    antisymCubic_not_velocityOnly, antisymCubic_ne_quadGrammar⟩

/-- The cubic survivor: the committed increment plus the antisymmetric
cubic.  It survives corner invisibility, the boost clause, and the
committed stationarity calibration while failing the velocity-only
clause. -/
noncomputable def cubicSurvivor : ℝ → ℝ → ℝ :=
  fun x y => twoSlotIncrement modeExtremalCurvature 0 x y
    + antisymCubic x y

/-- The cubic survivor fails the velocity-only clause. -/
theorem cubicSurvivor_not_velocityOnly : ¬ VelocityOnly cubicSurvivor := by
  intro h
  have hkey := h 2 0 0
  unfold cubicSurvivor twoSlotIncrement antisymCubic at hkey
  norm_num at hkey

/-- The cubic survivor differs from the committed increment. -/
theorem cubicSurvivor_ne_committed :
    cubicSurvivor ≠ twoSlotIncrement modeExtremalCurvature 0 := by
  intro h
  have hkey := congrFun (congrFun h 2) 0
  unfold cubicSurvivor twoSlotIncrement antisymCubic at hkey
  norm_num at hkey

/-- **The velocity-only clause stays load-bearing in the polynomial
grammar (S5).**  The cubic survivor is corner-invisible, satisfies the
boost clause, satisfies the committed stationarity calibration through
the committed forcing identity plus the neutral cubic, fails the
velocity-only clause, and differs from the committed increment.
Corner invisibility, stationarity, and the boost clause together do
not select the committed point; the velocity-only clause removes both
the committed first-slot direction and this cubic direction. -/
theorem velocityOnly_load_bearing_in_polynomial_grammar :
    CornerInvisible cubicSurvivor
      ∧ BoostNullCovariant cubicSurvivor
      ∧ ((chainCurvedD2 modeExtremalCurvature 1 1 + antisymCubicD2 1 1)
          + (chainCurvedD1 1 1 + antisymCubicD1 1 1) = 0)
      ∧ ¬ VelocityOnly cubicSurvivor
      ∧ cubicSurvivor ≠ twoSlotIncrement modeExtremalCurvature 0 := by
  refine ⟨?_, ?_, ?_, cubicSurvivor_not_velocityOnly,
    cubicSurvivor_ne_committed⟩
  · obtain ⟨h00, h10, h01, h11⟩ :=
      twoSlotIncrement_cornerInvisible modeExtremalCurvature 0
    obtain ⟨g00, g10, g01, g11⟩ := antisymCubic_cornerInvisible
    refine ⟨?_, ?_, ?_, ?_⟩
    · show twoSlotIncrement modeExtremalCurvature 0 0 0
        + antisymCubic 0 0 = 0
      rw [h00, g00]; ring
    · show twoSlotIncrement modeExtremalCurvature 0 1 0
        + antisymCubic 1 0 = 0
      rw [h10, g10]; ring
    · show twoSlotIncrement modeExtremalCurvature 0 0 1
        + antisymCubic 0 1 = 0
      rw [h01, g01]; ring
    · show twoSlotIncrement modeExtremalCurvature 0 1 1
        + antisymCubic 1 1 = 0
      rw [h11, g11]; ring
  · intro u₁ v₁ u₂ v₂
    refine NullLagrangian.congr (fun x y => ?_)
      (((twoSlotIncrement_boostNullCovariant modeExtremalCurvature 0)
          u₁ v₁ u₂ v₂).add
        (antisymCubic_boostNullCovariant u₁ v₁ u₂ v₂))
    simp only [boostDiff, cubicSurvivor]
    ring
  · have hforced :=
      (modeExtremal_forced modeExtremalCurvature).mpr rfl
    have hneutral := antisymCubic_stationarity_neutral
    linarith

/-! ## Composition: selection of the committed point in the full
polynomial grammar -/

/-- **Selection in the full polynomial grammar (S4).**  A two-slot
polynomial enrichment of any degree bound that satisfies the committed
velocity-only clause, corner invisibility, the boost clause, and the
committed stationarity calibration (stated through the certified
derivative `velocityPolyDeriv` of its velocity slice) is exactly the
committed point: the enriched Lagrangian equals the committed
mode-extremal rule pointwise, and the velocity-slice coefficients are
pinned. -/
theorem polyGrammar_selection_eq_committed (n : ℕ) (c : ℕ → ℕ → ℝ)
    (hvel : VelocityOnly (polyEnrichment n c))
    (hcorner : CornerInvisible (polyEnrichment n c))
    (hboost : BoostNullCovariant (polyEnrichment n c))
    (hstat : (chainFiberSlope 1
        + velocityPolyDeriv n (fun j => c 0 j) 1)
        + chainCurvedD1 1 1 = 0) :
    (∀ j, 3 ≤ j → j ≤ n → c 0 j = 0)
      ∧ c 0 0 = 0
      ∧ (2 ≤ n → c 0 1 = -(modeExtremalCurvature / 2)
          ∧ c 0 2 = modeExtremalCurvature / 2)
      ∧ (∀ x y : ℝ, polyEnrichment n c x y
          = twoSlotIncrement modeExtremalCurvature 0 x y)
      ∧ ∀ x y : ℝ,
          chainLogLagrangian x y + polyEnrichment n c x y
            = chainCurvedLagrangian modeExtremalCurvature x y := by
  -- the enrichment is its velocity slice
  have hslice : ∀ x y : ℝ,
      polyEnrichment n c x y = velocityPoly n (fun j => c 0 j) y := by
    intro x y
    rw [hvel x 0 y, polyEnrichment_velocity_slice]
  -- the boost clause pins the slice to a quadratic
  have hboost' : BoostNullCovariant
      (fun _ y => velocityPoly n (fun j => c 0 j) y) :=
    BoostNullCovariant.congr hslice hboost
  have hhigh := velocityPoly_boostNull_highCoeff_vanish n
    (fun j => c 0 j) hboost'
  obtain ⟨B, D, hBD, hderiv, hpin⟩ :=
    velocityPoly_quadratic_normal_form n (fun j => c 0 j) hhigh
  -- corner invisibility pins the constant and links the slopes
  obtain ⟨h00, -, h01, -⟩ := hcorner
  rw [hslice 0 0, hBD 0] at h00
  rw [hslice 0 1, hBD 1] at h01
  have hc00 : c 0 0 = 0 := by linarith
  have hBmD : B = -D := by linarith
  -- the committed stationarity calibration forces the curvature
  have hstat' : chainCurvedD2 (2 * D) 1 1 + chainCurvedD1 1 1 = 0 := by
    rw [hderiv] at hstat
    simp only [chainCurvedD2]
    linarith
  have hDval : 2 * D = modeExtremalCurvature :=
    (modeExtremal_forced (2 * D)).mp hstat'
  -- assemble
  have hpoint : ∀ x y : ℝ, polyEnrichment n c x y
      = twoSlotIncrement modeExtremalCurvature 0 x y := by
    intro x y
    rw [hslice x y, hBD y]
    unfold twoSlotIncrement
    have hB : B = -(modeExtremalCurvature / 2) := by linarith
    have hD : D = modeExtremalCurvature / 2 := by linarith
    rw [hc00, hB, hD]
    ring
  refine ⟨hhigh, hc00, ?_, hpoint, ?_⟩
  · intro hn
    obtain ⟨hB1, hD2⟩ := hpin hn
    constructor
    · rw [← hB1]; linarith
    · rw [← hD2]; linarith
  · intro x y
    rw [hpoint x y]
    unfold twoSlotIncrement chainCurvedLagrangian
    ring

/-- The velocity-slice coefficient family of the committed point. -/
noncomputable def committedCoeffs : ℕ → ℕ → ℝ := fun i j =>
  if i = 0 ∧ j = 1 then -(modeExtremalCurvature / 2)
  else if i = 0 ∧ j = 2 then modeExtremalCurvature / 2
  else 0

/-- The committed increment is a member of the polynomial grammar at
degree bound two. -/
theorem committedCoeffs_realizes_increment (x y : ℝ) :
    polyEnrichment 2 committedCoeffs x y
      = twoSlotIncrement modeExtremalCurvature 0 x y := by
  unfold polyEnrichment committedCoeffs twoSlotIncrement
  simp [Finset.sum_range_succ]
  ring

/-- **Nonvacuity of the axiom system (S4).**  The committed point,
realized in the polynomial grammar, satisfies the velocity-only
clause, corner invisibility, the boost clause, and the committed
stationarity calibration. -/
theorem committed_point_satisfies_all_clauses :
    VelocityOnly (polyEnrichment 2 committedCoeffs)
      ∧ CornerInvisible (polyEnrichment 2 committedCoeffs)
      ∧ BoostNullCovariant (polyEnrichment 2 committedCoeffs)
      ∧ (chainFiberSlope 1
          + velocityPolyDeriv 2 (fun j => committedCoeffs 0 j) 1)
          + chainCurvedD1 1 1 = 0 := by
  have hvelcommitted : VelocityOnly
      (twoSlotIncrement modeExtremalCurvature 0) :=
    (twoSlotIncrement_velocityOnly_iff modeExtremalCurvature 0).mpr rfl
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x x' y
    rw [committedCoeffs_realizes_increment,
      committedCoeffs_realizes_increment]
    exact hvelcommitted x x' y
  · obtain ⟨h00, h10, h01, h11⟩ :=
      twoSlotIncrement_cornerInvisible modeExtremalCurvature 0
    exact ⟨by rw [committedCoeffs_realizes_increment]; exact h00,
      by rw [committedCoeffs_realizes_increment]; exact h10,
      by rw [committedCoeffs_realizes_increment]; exact h01,
      by rw [committedCoeffs_realizes_increment]; exact h11⟩
  · exact BoostNullCovariant.congr
      (fun x y => (committedCoeffs_realizes_increment x y).symm)
      (twoSlotIncrement_boostNullCovariant modeExtremalCurvature 0)
  · have hforced :=
      (modeExtremal_forced modeExtremalCurvature).mpr rfl
    simp only [chainCurvedD2] at hforced
    have hderiv : velocityPolyDeriv 2 (fun j => committedCoeffs 0 j) 1
        = modeExtremalCurvature / 2 := by
      unfold velocityPolyDeriv committedCoeffs
      simp [Finset.sum_range_succ]
      ring
    rw [hderiv]
    linarith

/-! ## Remaining load-bearing receipts -/

/-- **Stationarity is load-bearing (S5).**  The committed counterfamily
point `(a*/2, 0)` is corner-invisible, velocity-only, and satisfies the
boost clause, yet fails the committed stationarity calibration and
differs from the committed rule. -/
theorem stationarity_clause_load_bearing :
    ∃ a c : ℝ,
      0 < a
        ∧ CornerInvisible (twoSlotIncrement a c)
        ∧ VelocityOnly (twoSlotIncrement a c)
        ∧ BoostNullCovariant (twoSlotIncrement a c)
        ∧ chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 ≠ 0
        ∧ chainTwoSlotCurvedLagrangian a c
            ≠ chainCurvedLagrangian modeExtremalCurvature := by
  have hpos := modeExtremalCurvature_pos
  refine ⟨modeExtremalCurvature / 2, 0, half_pos hpos,
    twoSlotIncrement_cornerInvisible _ _,
    (twoSlotIncrement_velocityOnly_iff _ _).mpr rfl,
    twoSlotIncrement_boostNullCovariant _ _, ?_, ?_⟩
  · rw [Ne, chainTwoSlot_modeExtremal_forced]
    intro h
    linarith
  · intro h
    have h2 : chainTwoSlotCurvedLagrangian (modeExtremalCurvature / 2) 0
        = chainTwoSlotCurvedLagrangian modeExtremalCurvature 0 := by
      rw [h, ← chainTwoSlot_zero_eq_curved]
    have hp : ((modeExtremalCurvature / 2 : ℝ), (0 : ℝ))
        = ((modeExtremalCurvature : ℝ), (0 : ℝ)) :=
      chainTwoSlot_parameters_injective h2
    have hfst := congrArg Prod.fst hp
    simp only at hfst
    linarith

/-- **Corner invisibility is load-bearing (S5).**  Shifting the
committed increment by a nonzero constant keeps the velocity-only
clause, the boost clause, and the committed stationarity calibration
(the constant has derivative zero), yet reads a shifted corner table
and yields a Lagrangian different from the committed rule. -/
theorem corner_clause_load_bearing :
    ∃ s : ℝ, s ≠ 0
      ∧ VelocityOnly
          (fun x y => twoSlotIncrement modeExtremalCurvature 0 x y + s)
      ∧ BoostNullCovariant
          (fun x y => twoSlotIncrement modeExtremalCurvature 0 x y + s)
      ∧ ¬ CornerInvisible
          (fun x y => twoSlotIncrement modeExtremalCurvature 0 x y + s)
      ∧ HasDerivAt
          (fun t => twoSlotIncrement modeExtremalCurvature 0 1 t + s)
          (modeExtremalCurvature / 2) 1
      ∧ ((chainFiberSlope 1 + modeExtremalCurvature / 2)
          + chainCurvedD1 1 1 = 0)
      ∧ (fun x y => chainLogLagrangian x y
            + (twoSlotIncrement modeExtremalCurvature 0 x y + s))
          ≠ chainCurvedLagrangian modeExtremalCurvature := by
  have hvelcommitted : VelocityOnly
      (twoSlotIncrement modeExtremalCurvature 0) :=
    (twoSlotIncrement_velocityOnly_iff modeExtremalCurvature 0).mpr rfl
  refine ⟨1, one_ne_zero, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x x' y
    dsimp only
    rw [hvelcommitted x x' y]
  · refine boostNullCovariant_of_quadratic _ 1 0
      (-(modeExtremalCurvature / 2)) 0 0 (modeExtremalCurvature / 2) ?_
    intro x y
    unfold twoSlotIncrement
    ring
  · intro h
    obtain ⟨h00, -, -, -⟩ := h
    unfold twoSlotIncrement at h00
    norm_num at h00
  · have hfun : (fun t : ℝ =>
        twoSlotIncrement modeExtremalCurvature 0 1 t + 1)
        = fun t : ℝ =>
            (modeExtremalCurvature / 2) * (t * t)
              - (modeExtremalCurvature / 2) * t + 1 := by
      funext t
      unfold twoSlotIncrement
      ring
    rw [hfun]
    have h := ((((hasDerivAt_id (1 : ℝ)).mul (hasDerivAt_id 1)).const_mul
      (modeExtremalCurvature / 2)).sub
        ((hasDerivAt_id (1 : ℝ)).const_mul
          (modeExtremalCurvature / 2))).add_const 1
    convert h using 1
    simp only [id_eq]
    ring
  · have hforced :=
      (modeExtremal_forced modeExtremalCurvature).mpr rfl
    simp only [chainCurvedD2] at hforced
    linarith
  · intro h
    have hkey := congrFun (congrFun h 0) 0
    simp only [chainCurvedLagrangian, twoSlotIncrement] at hkey
    norm_num at hkey

end OPH.Variational

#print axioms OPH.Variational.nullLagrangian_localAction_invariant
#print axioms OPH.Variational.telescoping_nullLagrangian
#print axioms OPH.Variational.constant_nullLagrangian
#print axioms OPH.Variational.nullLagrangian_velocitySentence_iff
#print axioms OPH.Variational.committed_increment_not_boostNull
#print axioms OPH.Variational.velocityPoly_hasDerivAt
#print axioms OPH.Variational.polyEnrichment_velocity_slice
#print axioms OPH.Variational.velocityPoly_quadratic_normal_form
#print axioms OPH.Variational.boostNullCovariant_of_quadratic
#print axioms OPH.Variational.quadPoly_boostNullCovariant
#print axioms OPH.Variational.twoSlotIncrement_boostNullCovariant
#print axioms OPH.Variational.velocityPoly_boostNull_highCoeff_vanish
#print axioms OPH.Variational.velocityPoly_boostNull_of_lowDegree
#print axioms OPH.Variational.velocityPoly_boostNull_iff
#print axioms OPH.Variational.quarticIncrement_not_boostNullCovariant
#print axioms OPH.Variational.auditQuartic_enrichment_eq
#print axioms OPH.Variational.auditQuartic_in_grammar
#print axioms OPH.Variational.auditQuarticEnrichment_not_boostNullCovariant
#print axioms OPH.Variational.audit_quartic_passes_committed_clauses_fails_boost
#print axioms OPH.Variational.antisymCubic_cornerInvisible
#print axioms OPH.Variational.antisymCubic_boostNullCovariant
#print axioms OPH.Variational.antisymCubic_ne_quadGrammar
#print axioms OPH.Variational.antisymCoeffs_realizes_cubic
#print axioms OPH.Variational.antisymCubic_not_velocityOnly
#print axioms OPH.Variational.antisymCubic_hasDerivAt_slot1
#print axioms OPH.Variational.antisymCubic_hasDerivAt_slot2
#print axioms OPH.Variational.antisymCubic_stationarity_neutral
#print axioms OPH.Variational.boost_clause_does_not_force_quadraticity
#print axioms OPH.Variational.cubicSurvivor_not_velocityOnly
#print axioms OPH.Variational.cubicSurvivor_ne_committed
#print axioms OPH.Variational.velocityOnly_load_bearing_in_polynomial_grammar
#print axioms OPH.Variational.polyGrammar_selection_eq_committed
#print axioms OPH.Variational.committedCoeffs_realizes_increment
#print axioms OPH.Variational.committed_point_satisfies_all_clauses
#print axioms OPH.Variational.stationarity_clause_load_bearing
#print axioms OPH.Variational.corner_clause_load_bearing

-- Expected axioms for every theorem above: propext, Classical.choice,
-- Quot.sound (real analysis and polynomial algebra via Mathlib).
-- No native_decide, no decide, no new axioms.
