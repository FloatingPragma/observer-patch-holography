import Mathlib

/-!
# Yang–Mills finite repair-gap — Proposition 8.1 (commuting-color spectral gap)

Formalisation of the **finite, real** §8 payoff of B. Müller, *Explaining the
Yang–Mills Mass Gap with Observer-Patch Repair Dynamics* (r1515):

> **Proposition 8.1.** Let `{Eₐ}_{a∈s}` be a finite family of **mutually
> commuting** orthogonal (star) projections on a real Hilbert space (the
> paper's case is finite-dimensional — formalized in the natural generality,
> only completeness is consumed), and let `P₀ = ∏ₐ Eₐ` be their
> (non-commutative) product — the
> orthogonal projection onto the joint fixed space `⋂ₐ Ran Eₐ`.
> Then for any constant rate `c_* > 0`,
> `      c_* · (I − P₀)  ≤  ∑ₐ c_* · (I − Eₐ),`
> a strictly positive finite-stage spectral gap `Δ ≥ c_*`.

The mathematical heart is the commuting-projection operator inequality
`one_sub_noncommProd_le_sum`:
`      I − ∏ₐ Eₐ  ≤  ∑ₐ (I − Eₐ)`   for a commuting family of star projections,
proved by induction on the finite index set with the two-projection step
`two_proj_le`. Everything downstream (`prop_8_1`, `prop_8_1_gap`) is scaling.

## Honest scope (read this first)

This is an **implication**, nothing more.  The star-projection, commutation, and
product-identification clauses are hypotheses here (`hE`, `hc`, `hprod`).  A
further physical identification of their joint range with a constants sector is
not encoded in the signature and is discharged nowhere in Lean.  These are the
paper's modelling obligations for Yang–Mills' actual relaxation.  This gap is
`Δ_rep`, the finite-stage
**representation** gap: it says **nothing** about `Δ_YM`. The continuum
certificate that would bridge them (Müller's **Assumption 9.2**: Schwinger
convergence, reflection positivity, Osterwalder–Schrader, non-triviality) is the
genuine open problem and is **untouched, unassumed, unclaimed**.

The core inequality `I − ∏ Eₐ ≤ ∑ (I − Eₐ)` needs **zero** physical input — it is
a clean finite-dimensional operator fact on any commuting family of orthogonal
projections. That fact is the real residue; the "mass gap" is rhetoric around it.

## Integration

The sibling `ObserverPatchHolography.YangMillsGap` assembly module re-exports
`prop_8_1` with the identical signature carried here; that obligation is
discharged by direct import of
`ObserverPatchHolography.YangMillsProp81.prop_8_1` below:
```
exact ObserverPatchHolography.YangMillsProp81.prop_8_1 s Ec P0 hE hc hprod hcpos
```

`thm_7_3_finite_gap` then combines this inequality with a uniform floor on
directly assumed positive collar rates.  It does not consume Lemma 7.2; no Lean
theorem currently bridges that uniform-fiber matrix coefficient to the collar
rate argument.

SCOPE: machine-checked here is the commuting-projection inequality and its
constant-rate finite representation-gap form.  The variable-rate assembly is
conditional on supplied positive rates, and the continuum reduction remains
conditional on Assumption 9.2, which is not touched here.
-/

namespace ObserverPatchHolography.YangMillsProp81

open scoped RealInnerProductSpace

/-! **Real Hilbert space** carrier of the repair generator. In the paper's
application `E` is finite-dimensional (whence complete), but only completeness
is actually consumed by the positivity/Loewner machinery — so we hypothesize
exactly `[CompleteSpace E]` and the finite-dimensional case is an instance. -/
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
  [CompleteSpace E]

/-! ## Loewner-order plumbing on `E →L[ℝ] E`

`E →L[ℝ] E` carries the Loewner partial order (`instLoewnerPartialOrder`,
`f ≤ g ↔ (g − f).IsPositive`). We only need two monotonicity facts; both are
one-liners from `le_def` + the additive/scalar positivity API. -/

omit [CompleteSpace E] in
/-- Left-add monotonicity for the Loewner order: `a ≤ b → c + a ≤ c + b`. -/
private lemma add_le_add_left_clm {a b : E →L[ℝ] E} (h : a ≤ b) (c : E →L[ℝ] E) :
    c + a ≤ c + b := by
  rw [ContinuousLinearMap.le_def] at h ⊢
  have e : (c + b) - (c + a) = b - a := by abel
  rwa [e]

omit [CompleteSpace E] in
/-- Nonnegative-scalar monotonicity for the Loewner order: `a ≤ b → 0 ≤ c →
`c • a ≤ c • b`. -/
private lemma smul_le_smul_clm {a b : E →L[ℝ] E} (h : a ≤ b) {c : ℝ} (hc : 0 ≤ c) :
    c • a ≤ c • b := by
  rw [ContinuousLinearMap.le_def] at h ⊢
  have e : c • b - c • a = c • (b - a) := by rw [smul_sub]
  rw [e]
  exact ContinuousLinearMap.IsPositive.smul_of_nonneg h hc

/-! ## The two-projection step -/

/-- **Two commuting projections.** For commuting orthogonal (star) projections
`p, q`, `I − p·q ≤ (I − p) + (I − q)`.

Proof: `((I−p)+(I−q)) − (I−p·q) = (I−p)·(I−q)`, and the product of the two
commuting star projections `I−p`, `I−q` is again a star projection, hence a
positive operator. -/
theorem two_proj_le {p q : E →L[ℝ] E}
    (hp : IsStarProjection p) (hq : IsStarProjection q) (h : Commute p q) :
    (1 - p * q) ≤ (1 - p) + (1 - q) := by
  rw [ContinuousLinearMap.le_def]
  have key : ((1 - p) + (1 - q)) - (1 - p * q) = (1 - p) * (1 - q) := by
    noncomm_ring
  rw [key]
  have hcomm : Commute (1 - p) (1 - q) :=
    (Commute.one_left (1 - q)).sub_left ((Commute.one_right p).sub_right h)
  exact ContinuousLinearMap.IsPositive.of_isStarProjection
    (hp.one_sub.mul hq.one_sub hcomm)

/-! ## The commuting product is a star projection -/

/-- The non-commutative product `∏_{a∈s} p a` of a **commuting** family of star
projections is itself a star projection (hence an orthogonal projection onto the
joint fixed space `⋂ Ran(p a)`). Induction on `s` via `Finset.cons_induction`;
each step is `IsStarProjection.mul`, whose commutativity side-condition comes from
`noncommProd_commute`. -/
theorem noncommProd_isStarProjection {ι : Type*} (p : ι → E →L[ℝ] E) :
    ∀ (s : Finset ι) (_hp : ∀ a ∈ s, IsStarProjection (p a))
      (hc : (↑s : Set ι).Pairwise (Function.onFun Commute p)),
      IsStarProjection (s.noncommProd p hc) := by
  intro s
  induction s using Finset.cons_induction with
  | empty =>
      intro _hp hc
      rw [Finset.noncommProd_empty]
      exact IsStarProjection.one _
  | cons a t hat ih =>
      intro hp hc
      have hcsub : (↑t : Set ι).Pairwise (Function.onFun Commute p) :=
        hc.mono fun _ => Finset.mem_cons.2 ∘ .inr
      have hpsub : ∀ x ∈ t, IsStarProjection (p x) :=
        fun x hx => hp x (Finset.mem_cons_of_mem hx)
      have hpa : IsStarProjection (p a) := hp a (Finset.mem_cons_self a t)
      have hcomm_aQ : Commute (p a) (t.noncommProd p hcsub) :=
        Finset.noncommProd_commute t p hcsub (p a)
          (fun x hx => hc (Finset.mem_cons_self a t) (Finset.mem_cons_of_mem hx)
            (by rintro rfl; exact hat hx))
      rw [Finset.noncommProd_cons]
      exact hpa.mul (ih hpsub hcsub) hcomm_aQ

/-! ## The heart inequality:  `I − ∏ Eₐ ≤ ∑ (I − Eₐ)` -/

/-- **Heart of Proposition 8.1.** For a finite **commuting** family of star
projections `p`, `I − ∏_{a∈s} p a ≤ ∑_{a∈s} (I − p a)`.

Induction on `s`: the empty product is `I` (both sides `0`); the `cons a t` step
factors `∏ = p a · ∏_t`, applies the two-projection bound `two_proj_le` with the
commuting pair `(p a, ∏_t)` (the tail product is a star projection by
`noncommProd_isStarProjection`, and commutes with `p a` by `noncommProd_commute`),
then the inductive hypothesis on the tail via left-add monotonicity. -/
theorem one_sub_noncommProd_le_sum {ι : Type*} (p : ι → E →L[ℝ] E) :
    ∀ (s : Finset ι) (_hp : ∀ a ∈ s, IsStarProjection (p a))
      (hc : (↑s : Set ι).Pairwise (Function.onFun Commute p)),
      (1 - s.noncommProd p hc) ≤ ∑ a ∈ s, (1 - p a) := by
  intro s
  induction s using Finset.cons_induction with
  | empty =>
      intro _hp _hc
      simp only [Finset.noncommProd_empty, sub_self, Finset.sum_empty, le_refl]
  | cons a t hat ih =>
      intro hp hc
      have hcsub : (↑t : Set ι).Pairwise (Function.onFun Commute p) :=
        hc.mono fun _ => Finset.mem_cons.2 ∘ .inr
      have hpsub : ∀ x ∈ t, IsStarProjection (p x) :=
        fun x hx => hp x (Finset.mem_cons_of_mem hx)
      have hpa : IsStarProjection (p a) := hp a (Finset.mem_cons_self a t)
      have hQ : IsStarProjection (t.noncommProd p hcsub) :=
        noncommProd_isStarProjection p t hpsub hcsub
      have hcomm_aQ : Commute (p a) (t.noncommProd p hcsub) :=
        Finset.noncommProd_commute t p hcsub (p a)
          (fun x hx => hc (Finset.mem_cons_self a t) (Finset.mem_cons_of_mem hx)
            (by rintro rfl; exact hat hx))
      rw [Finset.noncommProd_cons, Finset.sum_cons]
      calc (1 : E →L[ℝ] E) - p a * t.noncommProd p hcsub
          ≤ (1 - p a) + (1 - t.noncommProd p hcsub) := two_proj_le hpa hQ hcomm_aQ
        _ ≤ (1 - p a) + ∑ x ∈ t, (1 - p x) :=
            add_le_add_left_clm (ih hpsub hcsub) (1 - p a)

/-! ## Proposition 8.1 and corollaries -/

/-- **The joint projector is a projection.** `P₀ = ∏ₐ Eₐ` is a star projection,
so `I − P₀` is a genuine complementary orthogonal projection.  Formally this is
the complement of the joint fixed range; identifying that range with a physical
constants sector is a separate premise not represented here. -/
theorem prod_isStarProjection {ι : Type*} (s : Finset ι) (Ec : ι → (E →L[ℝ] E))
    (P0 : E →L[ℝ] E) (hE : ∀ a ∈ s, IsStarProjection (Ec a))
    (hc : (↑s : Set ι).Pairwise (Function.onFun Commute Ec))
    (hprod : s.noncommProd Ec hc = P0) :
    IsStarProjection P0 := by
  rw [← hprod]
  exact noncommProd_isStarProjection Ec s hE hc

/-- **Proposition 8.1 (commuting-color finite-stage gap).** For a finite family
of mutually commuting star projections `Ec` whose non-commutative product is the
joint projector `P₀`, and any constant rate `c_* > 0`,
`      c_* · (I − P₀)  ≤  ∑ₐ c_* · (I − Eₐ).`

Signature-compatible with `ObserverPatchHolography.YangMillsGap.prop_8_1`; this
theorem is re-exported there by direct import. -/
theorem prop_8_1 {ι : Type*} (s : Finset ι) (Ec : ι → (E →L[ℝ] E)) (P0 : E →L[ℝ] E)
    (hE : ∀ a ∈ s, IsStarProjection (Ec a))
    (hc : (↑s : Set ι).Pairwise (Function.onFun Commute Ec))
    (hprod : s.noncommProd Ec hc = P0)
    {cstar : ℝ} (hcpos : 0 < cstar) :
    cstar • ((1 : E →L[ℝ] E) - P0) ≤ ∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a) := by
  rw [← hprod, ← Finset.smul_sum]
  exact smul_le_smul_clm (one_sub_noncommProd_le_sum Ec s hE hc) hcpos.le

/-- **Rayleigh form of Proposition 8.1 (the "gap ≥ c_*" as usually quoted).**
On the complement of the joint fixed range (`P₀ x = 0`), the constant-rate
generator `L̃ = ∑ₐ c_* · (I − Eₐ)` satisfies `⟪L̃ x, x⟫ ≥ c_* · ‖x‖²`.
No physical constants-space identification is encoded. -/
theorem prop_8_1_gap {ι : Type*} (s : Finset ι) (Ec : ι → (E →L[ℝ] E)) (P0 : E →L[ℝ] E)
    (hE : ∀ a ∈ s, IsStarProjection (Ec a))
    (hc : (↑s : Set ι).Pairwise (Function.onFun Commute Ec))
    (hprod : s.noncommProd Ec hc = P0)
    {cstar : ℝ} (hcpos : 0 < cstar) (x : E) (hx : P0 x = 0) :
    cstar * ‖x‖ ^ 2 ≤ inner ℝ ((∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a)) x) x := by
  have hle : cstar • ((1 : E →L[ℝ] E) - P0) ≤ ∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a) :=
    prop_8_1 s Ec P0 hE hc hprod hcpos
  rw [ContinuousLinearMap.le_def] at hle
  have hnn := hle.inner_nonneg_left x
  have happ :
      ((∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a)) - cstar • ((1 : E →L[ℝ] E) - P0)) x
        = (∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a)) x - cstar • x := by
    rw [ContinuousLinearMap.sub_apply, ContinuousLinearMap.smul_apply,
      ContinuousLinearMap.sub_apply, ContinuousLinearMap.one_apply, hx, sub_zero]
  rw [happ, inner_sub_left, real_inner_smul_left, real_inner_self_eq_norm_sq] at hnn
  linarith

/-! ## Axiom self-audit (build-log visible)

Expected report for every theorem below: exactly
`[propext, Classical.choice, Quot.sound]` — no `sorryAx`, no project axiom. -/

#print axioms two_proj_le
#print axioms noncommProd_isStarProjection
#print axioms one_sub_noncommProd_le_sum
#print axioms prod_isStarProjection
#print axioms prop_8_1
#print axioms prop_8_1_gap

end ObserverPatchHolography.YangMillsProp81
