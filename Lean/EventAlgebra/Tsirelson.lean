import EventAlgebra.Basic

/-!
# The Tsirelson bound

The Tsirelson bound in operator-norm form, proved in a nontrivial unital
C*-ring in Mathlib's sense: for self-adjoint involutions `a₀, a₁, b₀, b₁` with each `aᵢ`
commuting with each `bⱼ`, the CHSH combination

  `S = a₀ b₀ + a₀ b₁ + a₁ b₀ − a₁ b₁`

satisfies `‖S‖ ≤ 2 √2`.

The proof is the classical square trick:

1. `chsh_mul_self` — the ring identity
   `S² = 4·1 − (a₀a₁ − a₁a₀)(b₀b₁ − b₁b₀)`, proved with no norm and no
   star, in a bare ring (sign convention: commutators `[x, y] = xy − yx`,
   and the minus sign in front of the product of commutators);
2. self-adjoint involutions have norm one (C*-identity), so each
   commutator has norm at most `2` and `‖S²‖ ≤ 4 + 4 = 8`;
3. `S` is self-adjoint, so `‖S‖² = ‖S* S‖ = ‖S²‖ ≤ 8`, i.e.
   `‖S‖ ≤ √8 = 2√2`.

No attainment (tightness) claim is made anywhere in this module: the bound
is an inequality, and nothing here asserts that `2√2` is achieved.

Mathlib's `Mathlib/Algebra/Star/CHSH.lean` proves the *order-form*
Tsirelson inequality `S ≤ √2^3 • 1` in an ordered star-algebra
(`tsirelson_inequality`).  The present module is the complementary
*norm-form* statement in a C*-ring; it assumes neither completeness, a
scalar-algebra structure, nor an order structure, only the C*-identity for
the norm.

The abstract statement instantiates in particular to
`Matrix (Fin n) (Fin n) ℂ` with the `L2` operator norm
(`matrix_tsirelson_bound` below), which is the ambient algebra of the
event-algebra modules; `EventAlgebra.IsEvent.one_sub_two_smul_involution`
produces self-adjoint involutions from events, connecting the two layers.

## Tagging convention

As in `EventAlgebra.Basic`.  Everything in this module is **algebra-only**:
no state, no trace pairing — only ring, star, and norm structure.
-/

namespace EventAlgebra

/-! ## The ring identity -/

section RingIdentity

variable {R : Type*} [Ring R]

/-- **Algebra-only.** The CHSH square identity in a bare ring: if
`a₀, a₁, b₀, b₁` are involutions and each `aᵢ` commutes with each `bⱼ`,
then

  `S² = 4·1 − (a₀a₁ − a₁a₀)(b₀b₁ − b₁b₀)`

for `S = a₀b₀ + a₀b₁ + a₁b₀ − a₁b₁`.  No star and no norm are involved;
this is the exact algebraic content of the Tsirelson bound. -/
theorem chsh_mul_self (a₀ a₁ b₀ b₁ : R)
    (ha₀ : a₀ * a₀ = 1) (ha₁ : a₁ * a₁ = 1)
    (hb₀ : b₀ * b₀ = 1) (hb₁ : b₁ * b₁ = 1)
    (h₀₀ : a₀ * b₀ = b₀ * a₀) (h₀₁ : a₀ * b₁ = b₁ * a₀)
    (h₁₀ : a₁ * b₀ = b₀ * a₁) (h₁₁ : a₁ * b₁ = b₁ * a₁) :
    (a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁) *
        (a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁) =
      4 • (1 : R) - (a₀ * a₁ - a₁ * a₀) * (b₀ * b₁ - b₁ * b₀) := by
  -- Rewrite systems: move every `b` to the right of every `a` in
  -- right-associated words, and cancel adjacent equal letters.
  have swap₀₀ : ∀ x : R, b₀ * (a₀ * x) = a₀ * (b₀ * x) := fun x => by
    rw [← mul_assoc, ← h₀₀, mul_assoc]
  have swap₀₁ : ∀ x : R, b₁ * (a₀ * x) = a₀ * (b₁ * x) := fun x => by
    rw [← mul_assoc, ← h₀₁, mul_assoc]
  have swap₁₀ : ∀ x : R, b₀ * (a₁ * x) = a₁ * (b₀ * x) := fun x => by
    rw [← mul_assoc, ← h₁₀, mul_assoc]
  have swap₁₁ : ∀ x : R, b₁ * (a₁ * x) = a₁ * (b₁ * x) := fun x => by
    rw [← mul_assoc, ← h₁₁, mul_assoc]
  have inv₀ : ∀ x : R, a₀ * (a₀ * x) = x := fun x => by
    rw [← mul_assoc, ha₀, one_mul]
  have inv₁ : ∀ x : R, a₁ * (a₁ * x) = x := fun x => by
    rw [← mul_assoc, ha₁, one_mul]
  simp only [mul_add, add_mul, mul_sub, sub_mul, mul_assoc,
    swap₀₀, swap₀₁, swap₁₀, swap₁₁, inv₀, inv₁,
    ha₀, ha₁, hb₀, hb₁, mul_one]
  abel

end RingIdentity

section EventCommutation

variable {R : Type*} [Ring R]

/-- **Algebra-only.** Commuting elements produce commuting dichotomic
expressions. This is the adapter used by the event-level CHSH client. -/
theorem one_sub_two_smul_commute {P Q : R} (h : P * Q = Q * P) :
    (1 - 2 • P) * (1 - 2 • Q) = (1 - 2 • Q) * (1 - 2 • P) := by
  calc
    (1 - 2 • P) * (1 - 2 • Q) =
        1 - 2 • P - 2 • Q + 4 • (P * Q) := by noncomm_ring
    _ = 1 - 2 • P - 2 • Q + 4 • (Q * P) := by rw [h]
    _ = (1 - 2 • Q) * (1 - 2 • P) := by noncomm_ring

end EventCommutation

/-! ## The norm bound -/

section CStarNorm

variable {A : Type*} [NormedRing A] [StarRing A] [CStarRing A] [Nontrivial A]

/-- **Algebra-only.** In a nontrivial unital C*-ring, a self-adjoint involution has
norm one: `‖a‖² = ‖a* a‖ = ‖1‖ = 1`. -/
theorem norm_eq_one_of_selfAdjoint_involution {a : A}
    (ha : IsSelfAdjoint a) (h : a * a = 1) : ‖a‖ = 1 := by
  have hsq : ‖a‖ * ‖a‖ = 1 := by
    rw [← CStarRing.norm_star_mul_self, ha.star_eq, h, norm_one]
  rcases mul_self_eq_one_iff.mp hsq with h1 | hneg
  · exact h1
  · exfalso
    have := norm_nonneg a
    rw [hneg] at this
    norm_num at this

omit [StarRing A] [CStarRing A] [Nontrivial A] in
/-- **Algebra-only.** The commutator of two norm-one elements has norm at
most `2`. -/
theorem norm_commutator_le_two {a₀ a₁ : A} (h₀ : ‖a₀‖ = 1) (h₁ : ‖a₁‖ = 1) :
    ‖a₀ * a₁ - a₁ * a₀‖ ≤ 2 :=
  calc ‖a₀ * a₁ - a₁ * a₀‖ ≤ ‖a₀ * a₁‖ + ‖a₁ * a₀‖ := norm_sub_le _ _
    _ ≤ ‖a₀‖ * ‖a₁‖ + ‖a₁‖ * ‖a₀‖ :=
        add_le_add (norm_mul_le _ _) (norm_mul_le _ _)
    _ = 2 := by rw [h₀, h₁]; norm_num

/-- **Algebra-only.** The **Tsirelson bound**, norm form: in a nontrivial
unital C*-ring, for self-adjoint involutions `a₀, a₁, b₀, b₁` with each `aᵢ`
commuting with each `bⱼ`,

  `‖a₀ b₀ + a₀ b₁ + a₁ b₀ − a₁ b₁‖ ≤ 2 √2`.

No attainment claim is made: this is an upper bound only. -/
theorem tsirelson_bound (a₀ a₁ b₀ b₁ : A)
    (sa₀ : IsSelfAdjoint a₀) (sa₁ : IsSelfAdjoint a₁)
    (sb₀ : IsSelfAdjoint b₀) (sb₁ : IsSelfAdjoint b₁)
    (ha₀ : a₀ * a₀ = 1) (ha₁ : a₁ * a₁ = 1)
    (hb₀ : b₀ * b₀ = 1) (hb₁ : b₁ * b₁ = 1)
    (h₀₀ : a₀ * b₀ = b₀ * a₀) (h₀₁ : a₀ * b₁ = b₁ * a₀)
    (h₁₀ : a₁ * b₀ = b₀ * a₁) (h₁₁ : a₁ * b₁ = b₁ * a₁) :
    ‖a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁‖ ≤ 2 * Real.sqrt 2 := by
  -- `S` is self-adjoint.
  have hSsa : IsSelfAdjoint (a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁) := by
    show star _ = _
    simp only [star_sub, star_add, star_mul,
      sa₀.star_eq, sa₁.star_eq, sb₀.star_eq, sb₁.star_eq]
    rw [← h₀₀, ← h₀₁, ← h₁₀, ← h₁₁]
  -- `‖S‖² = ‖S²‖ ≤ 4 + 4 = 8` via the square identity.
  have hnorm_sq :
      ‖a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁‖ ^ 2 ≤ 8 := by
    rw [← hSsa.norm_mul_self,
      chsh_mul_self a₀ a₁ b₀ b₁ ha₀ ha₁ hb₀ hb₁ h₀₀ h₀₁ h₁₀ h₁₁]
    have hcomm_a : ‖a₀ * a₁ - a₁ * a₀‖ ≤ 2 :=
      norm_commutator_le_two
        (norm_eq_one_of_selfAdjoint_involution sa₀ ha₀)
        (norm_eq_one_of_selfAdjoint_involution sa₁ ha₁)
    have hcomm_b : ‖b₀ * b₁ - b₁ * b₀‖ ≤ 2 :=
      norm_commutator_le_two
        (norm_eq_one_of_selfAdjoint_involution sb₀ hb₀)
        (norm_eq_one_of_selfAdjoint_involution sb₁ hb₁)
    calc ‖4 • (1 : A) - (a₀ * a₁ - a₁ * a₀) * (b₀ * b₁ - b₁ * b₀)‖
        ≤ ‖4 • (1 : A)‖ + ‖(a₀ * a₁ - a₁ * a₀) * (b₀ * b₁ - b₁ * b₀)‖ :=
          norm_sub_le _ _
      _ ≤ 4 * ‖(1 : A)‖ +
            ‖a₀ * a₁ - a₁ * a₀‖ * ‖b₀ * b₁ - b₁ * b₀‖ := by
          gcongr
          · simpa using norm_nsmul_le (n := 4) (a := (1 : A))
          · exact norm_mul_le _ _
      _ ≤ 4 * 1 + 2 * 2 := by
          rw [norm_one]
          have hmul : ‖a₀ * a₁ - a₁ * a₀‖ * ‖b₀ * b₁ - b₁ * b₀‖ ≤ 2 * 2 :=
            mul_le_mul hcomm_a hcomm_b (norm_nonneg _) (by norm_num)
          linarith
      _ = 8 := by norm_num
  -- Extract the square root.
  have hle : ‖a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁‖ ≤ Real.sqrt 8 :=
    (Real.le_sqrt (norm_nonneg _) (by norm_num)).mpr hnorm_sq
  have h8 : Real.sqrt 8 = 2 * Real.sqrt 2 := by
    rw [show (8 : ℝ) = 2 ^ 2 * 2 by norm_num,
      Real.sqrt_mul (by positivity) 2, Real.sqrt_sq (by norm_num)]
  rwa [h8] at hle

/-- **Algebra-only.** The Tsirelson bound restated against Mathlib's
`IsCHSHTuple` hypothesis bundle (`Mathlib.Algebra.Star.CHSH`), for direct
interoperability with Mathlib's order-form `tsirelson_inequality`. -/
theorem tsirelson_bound_of_isCHSHTuple {a₀ a₁ b₀ b₁ : A}
    (T : IsCHSHTuple a₀ a₁ b₀ b₁) :
    ‖a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁‖ ≤ 2 * Real.sqrt 2 := by
  have ha₀ : a₀ * a₀ = 1 := by have := T.A₀_inv; rwa [pow_two] at this
  have ha₁ : a₁ * a₁ = 1 := by have := T.A₁_inv; rwa [pow_two] at this
  have hb₀ : b₀ * b₀ = 1 := by have := T.B₀_inv; rwa [pow_two] at this
  have hb₁ : b₁ * b₁ = 1 := by have := T.B₁_inv; rwa [pow_two] at this
  exact tsirelson_bound a₀ a₁ b₀ b₁ T.A₀_sa T.A₁_sa T.B₀_sa T.B₁_sa
    ha₀ ha₁ hb₀ hb₁ T.A₀B₀_commutes T.A₀B₁_commutes
    T.A₁B₀_commutes T.A₁B₁_commutes

end CStarNorm

/-! ## Instantiation to finite matrix algebras -/

section MatrixInstance

open scoped Matrix.Norms.L2Operator

/-- **Algebra-only.** The Tsirelson bound for finite matrix algebras: the
abstract C*-statement instantiated at `Matrix (Fin n) (Fin n) ℂ` with the
`L2` operator norm (Mathlib's scoped `Matrix.Norms.L2Operator` C*-structure).
This covers, in particular, dichotomic observables `1 - 2P` built from
events via `IsEvent.one_sub_two_smul_involution`. -/
theorem matrix_tsirelson_bound {n : ℕ} [NeZero n]
    (a₀ a₁ b₀ b₁ : Matrix (Fin n) (Fin n) ℂ)
    (sa₀ : IsSelfAdjoint a₀) (sa₁ : IsSelfAdjoint a₁)
    (sb₀ : IsSelfAdjoint b₀) (sb₁ : IsSelfAdjoint b₁)
    (ha₀ : a₀ * a₀ = 1) (ha₁ : a₁ * a₁ = 1)
    (hb₀ : b₀ * b₀ = 1) (hb₁ : b₁ * b₁ = 1)
    (h₀₀ : a₀ * b₀ = b₀ * a₀) (h₀₁ : a₀ * b₁ = b₁ * a₀)
    (h₁₀ : a₁ * b₀ = b₀ * a₁) (h₁₁ : a₁ * b₁ = b₁ * a₁) :
    ‖a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁‖ ≤ 2 * Real.sqrt 2 := by
  haveI : Nontrivial (Matrix (Fin n) (Fin n) ℂ) := by
    refine ⟨0, 1, fun h => ?_⟩
    have h00 := Matrix.ext_iff.mpr h ⟨0, Nat.pos_of_ne_zero (NeZero.ne n)⟩
      ⟨0, Nat.pos_of_ne_zero (NeZero.ne n)⟩
    rw [Matrix.zero_apply, Matrix.one_apply_eq] at h00
    exact zero_ne_one h00
  exact tsirelson_bound a₀ a₁ b₀ b₁ sa₀ sa₁ sb₀ sb₁
    ha₀ ha₁ hb₀ hb₁ h₀₀ h₀₁ h₁₀ h₁₁

/-- **Algebra-only.** End-to-end event client: four projection events with
the cross-party commutation pattern give dichotomic observables satisfying
the operator-norm Tsirelson bound. -/
theorem matrix_tsirelson_bound_of_events {n : ℕ} [NeZero n]
    (A₀ A₁ B₀ B₁ : Matrix (Fin n) (Fin n) ℂ)
    (hA₀ : IsEvent A₀) (hA₁ : IsEvent A₁)
    (hB₀ : IsEvent B₀) (hB₁ : IsEvent B₁)
    (h₀₀ : A₀ * B₀ = B₀ * A₀) (h₀₁ : A₀ * B₁ = B₁ * A₀)
    (h₁₀ : A₁ * B₀ = B₀ * A₁) (h₁₁ : A₁ * B₁ = B₁ * A₁) :
    ‖(1 - 2 • A₀) * (1 - 2 • B₀) +
        (1 - 2 • A₀) * (1 - 2 • B₁) +
        (1 - 2 • A₁) * (1 - 2 • B₀) -
        (1 - 2 • A₁) * (1 - 2 • B₁)‖ ≤ 2 * Real.sqrt 2 := by
  obtain ⟨sa₀, ha₀⟩ := hA₀.one_sub_two_smul_involution
  obtain ⟨sa₁, ha₁⟩ := hA₁.one_sub_two_smul_involution
  obtain ⟨sb₀, hb₀⟩ := hB₀.one_sub_two_smul_involution
  obtain ⟨sb₁, hb₁⟩ := hB₁.one_sub_two_smul_involution
  exact matrix_tsirelson_bound
    (1 - 2 • A₀) (1 - 2 • A₁) (1 - 2 • B₀) (1 - 2 • B₁)
    sa₀ sa₁ sb₀ sb₁ ha₀ ha₁ hb₀ hb₁
    (one_sub_two_smul_commute h₀₀)
    (one_sub_two_smul_commute h₀₁)
    (one_sub_two_smul_commute h₁₀)
    (one_sub_two_smul_commute h₁₁)

end MatrixInstance

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms chsh_mul_self
#print axioms one_sub_two_smul_commute
#print axioms norm_eq_one_of_selfAdjoint_involution
#print axioms norm_commutator_le_two
#print axioms tsirelson_bound
#print axioms tsirelson_bound_of_isCHSHTuple
#print axioms matrix_tsirelson_bound
#print axioms matrix_tsirelson_bound_of_events

end EventAlgebra
