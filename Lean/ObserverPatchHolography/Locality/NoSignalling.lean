import Mathlib

namespace OPH.Locality

open Matrix

/-!
# Algebraic no-signalling for local operations

A local operation on one factor of a bipartite system cannot change
the marginal seen by the other factor. This module proves the exact
finite statement in both descriptions used by the corpus. In the
classical description, pushing a joint law through a kernel that acts
on the first factor and fixes the second leaves the second marginal
unchanged, with row normalization as the only hypothesis. In the
finite quantum description, conjugating a joint matrix by a local
Kraus family lifted with the identity leaves the partial trace over
the acted factor unchanged, with Kraus completeness as the only
hypothesis. A computable negative control shows the normalization
hypothesis carries the content: a non-normalized local kernel rescales
the remote marginal.

No physical spacelike separation, tensor-factor identification, or
laboratory operation enters; the bipartite split is a named input.
-/

section Classical

variable {α β : Type*} [Fintype α] [Fintype β]

/-- The second marginal of a joint law. -/
noncomputable def sndMarginal (p : α × β → ℝ) (b : β) : ℝ :=
  ∑ a, p (a, b)

/-- A kernel on the first factor, lifted to act locally on the joint
space: it moves the first component and fixes the second. -/
noncomputable def liftFst [DecidableEq β] (K : α → α → ℝ)
    (z z' : α × β) : ℝ :=
  if z'.2 = z.2 then K z.1 z'.1 else 0

/-- Pushforward of a joint law through a joint kernel. -/
noncomputable def pushJoint (p : α × β → ℝ)
    (K : α × β → α × β → ℝ) (z' : α × β) : ℝ :=
  ∑ z, p z * K z z'

/-- **Classical no-signalling.** A row-normalized kernel acting on the
first factor leaves the second marginal exactly unchanged. -/
theorem sndMarginal_pushJoint_liftFst [DecidableEq β]
    (p : α × β → ℝ) (K : α → α → ℝ)
    (hK1 : ∀ a, ∑ a', K a a' = 1) (b : β) :
    sndMarginal (pushJoint p (liftFst K)) b = sndMarginal p b := by
  unfold sndMarginal pushJoint liftFst
  rw [Finset.sum_comm]
  have hz : ∀ z : α × β,
      ∑ a', p z * (if b = z.2 then K z.1 a' else 0)
        = if b = z.2 then p z else 0 := by
    intro z
    by_cases hb : b = z.2
    · rw [if_pos hb]
      calc ∑ a', p z * (if b = z.2 then K z.1 a' else 0)
          = ∑ a', p z * K z.1 a' :=
            Finset.sum_congr rfl fun a' _ => by rw [if_pos hb]
        _ = p z := by rw [← Finset.mul_sum, hK1 z.1, mul_one]
    · rw [if_neg hb]
      apply Finset.sum_eq_zero
      intro a' _
      rw [if_neg hb, mul_zero]
  calc ∑ z : α × β, ∑ a', p z * (if b = z.2 then K z.1 a' else 0)
      = ∑ z : α × β, (if b = z.2 then p z else 0) :=
        Finset.sum_congr rfl fun z _ => hz z
    _ = ∑ a, p (a, b) := by
        rw [Fintype.sum_prod_type]
        apply Finset.sum_congr rfl
        intro a _
        simp

/-- **Negative control.** Without row normalization the remote
marginal moves: the doubling kernel on a one-point first factor
doubles the second marginal. The normalization hypothesis carries the
entire content of the no-signalling theorem. -/
theorem signalling_without_row_normalization :
    ∃ (K : Fin 1 → Fin 1 → ℝ) (p : Fin 1 × Fin 1 → ℝ) (b : Fin 1),
      sndMarginal (pushJoint p (liftFst K)) b ≠ sndMarginal p b := by
  refine ⟨fun _ _ => 2, fun _ => 1, 0, ?_⟩
  unfold sndMarginal pushJoint liftFst
  simp only [Fintype.sum_prod_type, Fin.sum_univ_one, Fin.isValue]
  norm_num

end Classical

section Quantum

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α]
  [DecidableEq β]

/-- Partial trace over the first factor of a joint matrix. -/
noncomputable def ptraceFst (M : Matrix (α × β) (α × β) ℂ) :
    Matrix β β ℂ :=
  fun b b' => ∑ a, M (a, b) (a, b')

/-- A matrix on the first factor, lifted with the identity on the
second factor. -/
noncomputable def liftLeft (A : Matrix α α ℂ) :
    Matrix (α × β) (α × β) ℂ :=
  fun z z' => A z.1 z'.1 * (if z.2 = z'.2 then 1 else 0)

omit [Fintype α] [Fintype β] [DecidableEq α] in
theorem liftLeft_conjTranspose (A : Matrix α α ℂ) :
    (liftLeft (β := β) A)ᴴ = liftLeft Aᴴ := by
  ext z z'
  simp only [Matrix.conjTranspose_apply, liftLeft]
  by_cases h : z'.2 = z.2
  · simp [h]
  · have h' : ¬ z.2 = z'.2 := fun hc => h hc.symm
    simp [h, h']

/-- The lifted-conjugation summand in explicit entry form. -/
theorem liftLeft_conj_apply (Kk : Matrix α α ℂ)
    (M : Matrix (α × β) (α × β) ℂ) (z z' : α × β) :
    (liftLeft (β := β) Kk * M * (liftLeft (β := β) Kk)ᴴ) z z'
      = ∑ c, ∑ c', Kk z.1 c * M (c, z.2) (c', z'.2)
          * star (Kk z'.1 c') := by
  rw [liftLeft_conjTranspose]
  rw [Matrix.mul_apply]
  have hmid : ∀ w : α × β,
      (liftLeft (β := β) Kk * M) z w
        = ∑ c, Kk z.1 c * M (c, z.2) w := by
    intro w
    rw [Matrix.mul_apply, Fintype.sum_prod_type]
    apply Finset.sum_congr rfl
    intro c _
    rw [Finset.sum_eq_single z.2
      (fun v _ hv => by simp [liftLeft, Ne.symm hv])
      (fun habs => absurd (Finset.mem_univ z.2) habs)]
    simp [liftLeft]
  calc ∑ w, (liftLeft (β := β) Kk * M) z w * liftLeft Kkᴴ w z'
      = ∑ w : α × β, (∑ c, Kk z.1 c * M (c, z.2) w)
          * liftLeft Kkᴴ w z' := by
        apply Finset.sum_congr rfl
        intro w _
        rw [hmid w]
    _ = ∑ c', (∑ c, Kk z.1 c * M (c, z.2) (c', z'.2))
          * Kkᴴ c' z'.1 := by
        rw [Fintype.sum_prod_type]
        apply Finset.sum_congr rfl
        intro c' _
        rw [Finset.sum_eq_single z'.2
          (fun v _ hv => by simp [liftLeft, hv])
          (fun habs => absurd (Finset.mem_univ z'.2) habs)]
        simp [liftLeft]
    _ = ∑ c', ∑ c, Kk z.1 c * M (c, z.2) (c', z'.2)
          * star (Kk z'.1 c') := by
        apply Finset.sum_congr rfl
        intro c' _
        rw [Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro c _
        rw [Matrix.conjTranspose_apply]
    _ = ∑ c, ∑ c', Kk z.1 c * M (c, z.2) (c', z'.2)
          * star (Kk z'.1 c') := Finset.sum_comm

/-- **Quantum no-signalling.** For a Kraus-complete local family on
the first factor, conjugating the joint matrix by the lifted family
leaves the partial trace over that factor exactly unchanged. -/
theorem ptraceFst_local_kraus {ι : Type*} (s : Finset ι)
    (Kr : ι → Matrix α α ℂ)
    (hK : ∑ i ∈ s, (Kr i)ᴴ * Kr i = 1)
    (M : Matrix (α × β) (α × β) ℂ) :
    ptraceFst
      (∑ i ∈ s, liftLeft (β := β) (Kr i) * M
        * (liftLeft (β := β) (Kr i))ᴴ)
      = ptraceFst M := by
  ext b b'
  unfold ptraceFst
  have hgram : ∀ c c' : α,
      ((∑ i ∈ s, (Kr i)ᴴ * Kr i) c' c)
        = ∑ i ∈ s, ∑ a, star (Kr i a c') * Kr i a c := by
    intro c c'
    rw [Matrix.sum_apply]
    apply Finset.sum_congr rfl
    intro i _
    rw [Matrix.mul_apply]
    apply Finset.sum_congr rfl
    intro a _
    rw [Matrix.conjTranspose_apply]
  calc ∑ a, (∑ i ∈ s, liftLeft (β := β) (Kr i) * M
        * (liftLeft (β := β) (Kr i))ᴴ) (a, b) (a, b')
      = ∑ a, ∑ i ∈ s, ∑ c, ∑ c',
          Kr i a c * M (c, b) (c', b') * star (Kr i a c') := by
        apply Finset.sum_congr rfl
        intro a _
        rw [Matrix.sum_apply]
        apply Finset.sum_congr rfl
        intro i _
        exact liftLeft_conj_apply (Kr i) M (a, b) (a, b')
    _ = ∑ i ∈ s, ∑ c, ∑ c', ∑ a,
          Kr i a c * M (c, b) (c', b') * star (Kr i a c') := by
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro i _
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro c _
        rw [Finset.sum_comm]
    _ = ∑ c, ∑ c', M (c, b) (c', b')
          * ∑ i ∈ s, ∑ a, star (Kr i a c') * Kr i a c := by
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro c _
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro c' _
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro i _
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro a _
        ring
    _ = ∑ c, ∑ c', M (c, b) (c', b')
          * ((1 : Matrix α α ℂ) c' c) := by
        apply Finset.sum_congr rfl
        intro c _
        apply Finset.sum_congr rfl
        intro c' _
        rw [← hgram c c', hK]
    _ = ∑ a, M (a, b) (a, b') := by
        apply Finset.sum_congr rfl
        intro c _
        rw [Finset.sum_eq_single c
          (fun c' _ hc' => by rw [Matrix.one_apply_ne hc', mul_zero])
          (fun habs => absurd (Finset.mem_univ c) habs)]
        rw [Matrix.one_apply_eq, mul_one]

end Quantum

end OPH.Locality

#print axioms OPH.Locality.sndMarginal_pushJoint_liftFst
#print axioms OPH.Locality.signalling_without_row_normalization
#print axioms OPH.Locality.ptraceFst_local_kraus
