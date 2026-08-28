import Mathlib

/-!
# Signed-incidence Gram matrix: kinetic identity, spectral bridge, balance-to-injectivity

Model (from the OPH source-gap module, as pinned by the scout report):
for each edge `e` with endpoints `src e`, `tgt e` and sign `σ e`,
  `(D x) e = σ e * x (tgt e) - x (src e)`,
and the "signed Laplacian" is the Gram matrix `L = Dᵀ * D`.
The gap is the least eigenvalue of `L`; it is positive iff `L` is positive definite
iff `D` has trivial kernel.
-/

open Matrix

namespace GapScout

/-! ## Scout statements, verbatim -/

def signedDerivative {V E : Type} [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) : Matrix E V ℝ :=
  fun e v => if v = tgt e then σ e else if v = src e then -1 else 0

/-- One row of `D x`: the signed difference across the edge. Needs `src e ≠ tgt e`
so the two `if` branches never overlap. -/
theorem signedDerivative_mulVec_apply {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e) (x : V → ℝ) (e : E) :
    (signedDerivative src tgt σ *ᵥ x) e = σ e * x (tgt e) - x (src e) := by
  have h : ∀ v, signedDerivative src tgt σ e v * x v =
      (if v = tgt e then σ e * x (tgt e) else 0) +
        (if v = src e then -x (src e) else 0) := by
    intro v
    unfold signedDerivative
    by_cases h1 : v = tgt e
    · subst h1; simp [(hloop e).symm]
    · by_cases h2 : v = src e
      · subst h2; simp [h1]
      · simp [h1, h2]
  simp only [mulVec, dotProduct, h, Finset.sum_add_distrib, Finset.sum_ite_eq',
    Finset.mem_univ, if_true]
  ring

/-- `x ⬝ (Dᵀ D x) = (D x) ⬝ (D x)` for any rectangular real matrix. -/
theorem dotProduct_gram_mulVec {V E : Type} [Fintype V] [Fintype E]
    (D : Matrix E V ℝ) (x : V → ℝ) :
    x ⬝ᵥ ((Dᵀ * D) *ᵥ x) = (D *ᵥ x) ⬝ᵥ (D *ᵥ x) := by
  rw [← mulVec_mulVec, dotProduct_mulVec, vecMul_transpose]

/-- STEP 1 (scout statement, verbatim). -/
theorem signed_kinetic_identity {V E : Type} [Fintype V] [Fintype E] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e) (x : V → ℝ) :
    x ⬝ᵥ (((signedDerivative src tgt σ)ᵀ * signedDerivative src tgt σ) *ᵥ x) =
      ∑ e, (σ e * x (tgt e) - x (src e)) ^ 2 := by
  rw [dotProduct_gram_mulVec]
  simp only [dotProduct, signedDerivative_mulVec_apply src tgt σ hloop x, sq]

/-- STEP 2 (scout statement, verbatim). -/
theorem gram_posDef_of_injective {V E : Type} [Fintype V] [Fintype E] [DecidableEq V]
    (D : Matrix E V ℝ) (hD : Function.Injective fun x : V → ℝ => D *ᵥ x) :
    (Dᵀ * D).PosDef := by
  rw [posDef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · have := isHermitian_conjTranspose_mul_self D
    simpa [conjTranspose_eq_transpose_of_trivial] using this
  · intro x hx
    rw [star_trivial, dotProduct_gram_mulVec]
    have hDx : D *ᵥ x ≠ 0 := by
      intro h0
      apply hx
      apply hD
      simpa [mulVec_zero] using h0
    have hnn : 0 ≤ (D *ᵥ x) ⬝ᵥ (D *ᵥ x) := by
      simpa using dotProduct_star_self_nonneg (D *ᵥ x)
    rcases hnn.lt_or_eq with hlt | heq
    · exact hlt
    · exact absurd (dotProduct_self_eq_zero.mp heq.symm) hDx

/-- STEP 3: the all-negative triangle, L = D + A (signless Laplacian on an odd cycle).
This is a CONTROL computation, not evidence about the load-bearing domain. -/
def triangleL : Matrix (Fin 3) (Fin 3) ℝ :=
  ![![2, 1, 1], ![1, 2, 1], ![1, 1, 2]]

theorem triangle_quadratic_form (x : Fin 3 → ℝ) :
    x ⬝ᵥ (triangleL *ᵥ x) = (x 0 + x 1) ^ 2 + (x 0 + x 2) ^ 2 + (x 1 + x 2) ^ 2 := by
  simp [triangleL, dotProduct, mulVec, Fin.sum_univ_three]
  ring

/-- STEP 3 (scout statement, verbatim). -/
theorem triangle_posDef : triangleL.PosDef := by
  rw [posDef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [triangleL, conjTranspose]
  · intro x hx
    rw [star_trivial, triangle_quadratic_form]
    by_contra hle
    push_neg at hle
    apply hx
    have h01 := sq_nonneg (x 0 + x 1)
    have h02 := sq_nonneg (x 0 + x 2)
    have h12 := sq_nonneg (x 1 + x 2)
    have e01 : x 0 + x 1 = 0 := by nlinarith
    have e02 : x 0 + x 2 = 0 := by nlinarith
    have e12 : x 1 + x 2 = 0 := by nlinarith
    ext i
    fin_cases i <;> simp <;> linarith

/-! ## STEP 4: signed balance to injectivity (the unformalized layer)

A *signed walk* from `u` to `v` traverses edges forwards (`src → tgt`) or backwards
(`tgt → src`), multiplying the sign `σ e` of each traversed edge into a running product.
`SignedReach src tgt σ u v s` says: there is a walk from `u` to `v` whose sign product is `s`.
A *negative cycle at `u`* is a closed walk `SignedReach u u (-1)`.
-/

inductive SignedReach {V E : Type} (src tgt : E → V) (σ : E → ℝ) : V → V → ℝ → Prop
  | refl (u : V) : SignedReach src tgt σ u u 1
  | fwd {u : V} {s : ℝ} (e : E) (h : SignedReach src tgt σ u (src e) s) :
      SignedReach src tgt σ u (tgt e) (s * σ e)
  | bwd {u : V} {s : ℝ} (e : E) (h : SignedReach src tgt σ u (tgt e) s) :
      SignedReach src tgt σ u (src e) (s * σ e)

/-- The kernel equations, edge by edge. -/
theorem kernel_edge_eq {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e) (x : V → ℝ)
    (hx : signedDerivative src tgt σ *ᵥ x = 0) (e : E) :
    σ e * x (tgt e) = x (src e) := by
  have := congrFun hx e
  rw [signedDerivative_mulVec_apply src tgt σ hloop x e] at this
  simp only [Pi.zero_apply] at this
  linarith

/-- Transport: a kernel section propagates along signed walks, picking up the sign product.
Needs `σ e = ±1` so that the forward equation `σ x(tgt) = x(src)` can be inverted. -/
theorem kernel_transport {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1) (x : V → ℝ)
    (hx : signedDerivative src tgt σ *ᵥ x = 0)
    {u v : V} {s : ℝ} (hr : SignedReach src tgt σ u v s) : x v = s * x u := by
  induction hr with
  | refl => simp
  | fwd e h ih =>
    have he := kernel_edge_eq src tgt σ hloop x hx e
    rcases hσ e with h1 | h1 <;> rw [h1] at he ⊢ <;> linarith
  | bwd e h ih =>
    have he := kernel_edge_eq src tgt σ hloop x hx e
    rw [ih] at he
    rcases hσ e with h1 | h1 <;> rw [h1] at he ⊢ <;> linarith

/-- A negative closed walk at `u` kills the section at `u`. -/
theorem kernel_zero_at_negative_cycle {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1) (x : V → ℝ)
    (hx : signedDerivative src tgt σ *ᵥ x = 0)
    {u : V} (hneg : SignedReach src tgt σ u u (-1)) : x u = 0 := by
  have := kernel_transport src tgt σ hloop hσ x hx hneg
  linarith

/-- STEP 4, main theorem. If some vertex `u` carries a negative cycle and every vertex is
reachable from `u` by a signed walk (i.e. the graph is connected), the kernel of `D` is trivial. -/
theorem kernel_trivial_of_negative_cycle {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1)
    (u : V) (hneg : SignedReach src tgt σ u u (-1))
    (hconn : ∀ v, ∃ s, SignedReach src tgt σ u v s)
    (x : V → ℝ) (hx : signedDerivative src tgt σ *ᵥ x = 0) : x = 0 := by
  have hu := kernel_zero_at_negative_cycle src tgt σ hloop hσ x hx hneg
  funext v
  obtain ⟨s, hs⟩ := hconn v
  rw [kernel_transport src tgt σ hloop hσ x hx hs, hu, mul_zero]
  rfl

/-- STEP 4, injectivity form. -/
theorem injective_of_negative_cycle {V E : Type} [Fintype V] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1)
    (u : V) (hneg : SignedReach src tgt σ u u (-1))
    (hconn : ∀ v, ∃ s, SignedReach src tgt σ u v s) :
    Function.Injective fun x : V → ℝ => signedDerivative src tgt σ *ᵥ x := by
  intro x y hxy
  have h : signedDerivative src tgt σ *ᵥ (x - y) = 0 := by
    rw [mulVec_sub]
    simpa [sub_eq_zero] using hxy
  have := kernel_trivial_of_negative_cycle src tgt σ hloop hσ u hneg hconn (x - y) h
  exact sub_eq_zero.mp this

/-- STEP 4 ∘ STEP 2: the signed Laplacian of a connected graph with a negative cycle is
positive definite. -/
theorem gram_posDef_of_negative_cycle {V E : Type} [Fintype V] [Fintype E] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1)
    (u : V) (hneg : SignedReach src tgt σ u u (-1))
    (hconn : ∀ v, ∃ s, SignedReach src tgt σ u v s) :
    ((signedDerivative src tgt σ)ᵀ * signedDerivative src tgt σ).PosDef :=
  gram_posDef_of_injective _ (injective_of_negative_cycle src tgt σ hloop hσ u hneg hconn)

/-- Spectral reading: every eigenvalue of that Laplacian is strictly positive, i.e. the gap
(least eigenvalue) is positive. -/
theorem eigenvalues_pos_of_negative_cycle {V E : Type} [Fintype V] [Fintype E] [DecidableEq V]
    (src tgt : E → V) (σ : E → ℝ) (hloop : ∀ e, src e ≠ tgt e)
    (hσ : ∀ e, σ e = 1 ∨ σ e = -1)
    (u : V) (hneg : SignedReach src tgt σ u u (-1))
    (hconn : ∀ v, ∃ s, SignedReach src tgt σ u v s) (i : V) :
    0 < (gram_posDef_of_negative_cycle src tgt σ hloop hσ u hneg hconn).1.eigenvalues i :=
  (gram_posDef_of_negative_cycle src tgt σ hloop hσ u hneg hconn).eigenvalues_pos i

/-! ## The triangle, re-derived through STEP 4 (so step 3 is not the only evidence) -/

/-- The all-negative 3-cycle: edges `i → i+1`, every sign `-1`. -/
def triSrc : Fin 3 → Fin 3 := id
def triTgt : Fin 3 → Fin 3 := fun i => i + 1
def triσ : Fin 3 → ℝ := fun _ => -1

theorem triangleL_eq_gram :
    triangleL = (signedDerivative triSrc triTgt triσ)ᵀ * signedDerivative triSrc triTgt triσ := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [triangleL, signedDerivative, triSrc, triTgt, triσ, Matrix.mul_apply, Fin.sum_univ_three] <;>
    norm_num

theorem tri_negative_cycle : SignedReach triSrc triTgt triσ 0 0 (-1) := by
  have h1 : SignedReach triSrc triTgt triσ 0 (triTgt 0) (1 * triσ 0) :=
    SignedReach.fwd 0 (SignedReach.refl 0)
  have h2 : SignedReach triSrc triTgt triσ 0 (triTgt 1) (1 * triσ 0 * triσ 1) :=
    SignedReach.fwd 1 h1
  have h3 : SignedReach triSrc triTgt triσ 0 (triTgt 2) (1 * triσ 0 * triσ 1 * triσ 2) :=
    SignedReach.fwd 2 h2
  have e1 : triTgt 2 = (0 : Fin 3) := by decide
  have e2 : (1 * triσ 0 * triσ 1 * triσ 2 : ℝ) = -1 := by simp [triσ]
  rw [e1, e2] at h3
  exact h3

theorem tri_connected : ∀ v : Fin 3, ∃ s, SignedReach triSrc triTgt triσ 0 v s := by
  intro v
  fin_cases v
  · exact ⟨1, SignedReach.refl 0⟩
  · exact ⟨1 * triσ 0, SignedReach.fwd 0 (SignedReach.refl 0)⟩
  · exact ⟨1 * triσ 0 * triσ 1, SignedReach.fwd 1 (SignedReach.fwd 0 (SignedReach.refl 0))⟩

/-- Step 3 again, now as a corollary of the general balance theorem. -/
theorem triangle_posDef_via_balance : triangleL.PosDef := by
  rw [triangleL_eq_gram]
  exact gram_posDef_of_negative_cycle triSrc triTgt triσ (by decide)
    (fun _ => Or.inr rfl) 0 tri_negative_cycle tri_connected

/-! ## NEGATIVE CONTROL: the all-negative square

Edges `i → i+1` on `Fin 4`, every sign `-1`. The 4-cycle has sign product `(-1)^4 = +1`,
so the graph is *balanced*; the alternating section `(1,-1,1,-1)` is a zero mode.
The Gram matrix must NOT be positive definite, and our hypothesis `SignedReach 0 0 (-1)`
must be UNSATISFIABLE. Both are checked below. -/

def sqSrc : Fin 4 → Fin 4 := id
def sqTgt : Fin 4 → Fin 4 := fun i => i + 1
def sqσ : Fin 4 → ℝ := fun _ => -1

/-- The alternating zero mode. -/
def sqMode : Fin 4 → ℝ := ![1, -1, 1, -1]

theorem sqMode_ne_zero : sqMode ≠ 0 := by
  intro h
  have := congrFun h 0
  simp [sqMode] at this

theorem sqMode_in_kernel : signedDerivative sqSrc sqTgt sqσ *ᵥ sqMode = 0 := by
  ext e
  fin_cases e <;>
    simp [signedDerivative, sqSrc, sqTgt, sqσ, sqMode, mulVec, dotProduct, Fin.sum_univ_four]

theorem square_not_injective :
    ¬ Function.Injective fun x : Fin 4 → ℝ => signedDerivative sqSrc sqTgt sqσ *ᵥ x := by
  intro hinj
  apply sqMode_ne_zero
  apply hinj
  simp only [sqMode_in_kernel, mulVec_zero]

theorem square_not_posDef :
    ¬ ((signedDerivative sqSrc sqTgt sqσ)ᵀ * signedDerivative sqSrc sqTgt sqσ).PosDef := by
  intro hpd
  rw [posDef_iff_dotProduct_mulVec] at hpd
  have := hpd.2 sqMode_ne_zero
  rw [star_trivial, dotProduct_gram_mulVec, sqMode_in_kernel] at this
  simp at this

/-- Switching potential on the square: every signed walk from `u` to `v` has sign
`p u * p v`, where `p` is the alternating sign. So no closed walk is negative: the
hypothesis of the balance theorem genuinely fails here, not merely its conclusion. -/
def sqPot : Fin 4 → ℝ := ![1, -1, 1, -1]

theorem sqPot_succ (e : Fin 4) : sqPot (e + 1) = - sqPot e := by
  fin_cases e <;> simp [sqPot]

theorem sqPot_sq (v : Fin 4) : sqPot v * sqPot v = 1 := by
  fin_cases v <;> simp [sqPot]

theorem square_walk_sign {u v : Fin 4} {s : ℝ}
    (h : SignedReach sqSrc sqTgt sqσ u v s) : s = sqPot u * sqPot v := by
  induction h with
  | refl => rw [sqPot_sq]
  | fwd e _ ih =>
    rw [ih]
    simp only [sqSrc, sqTgt, sqσ, id, sqPot_succ]
    ring
  | bwd e _ ih =>
    rw [ih]
    simp only [sqSrc, sqTgt, sqσ, id, sqPot_succ]
    ring

theorem square_no_negative_cycle (u : Fin 4) : ¬ SignedReach sqSrc sqTgt sqσ u u (-1) := by
  intro h
  have := square_walk_sign h
  rw [sqPot_sq] at this
  norm_num at this

end GapScout
