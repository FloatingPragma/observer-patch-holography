import Mathlib

/-!
# Finite projection-event calculus: events, states, and Born weights

This module develops the elementary theory of finite-dimensional matrix
algebras equipped with projection events over `Matrix (Fin n) (Fin n) ℂ`:

* `EventAlgebra.IsEvent` — an **event** is a Hermitian idempotent
  (an orthogonal projection);
* `EventAlgebra.IsState` — a **state** is a positive-semidefinite matrix of
  trace one (a density matrix);
* `EventAlgebra.bornWeight` — the **Born weight** `Tr(ρ P)` of an event `P`
  under a state `ρ`, with reality, nonnegativity, normalisation, additivity
  on orthogonal events, the complement bound `≤ 1`, and monotonicity under
  subevents.

Nonnegativity and the upper bound are stated in the partial order of `ℂ`
(`0 ≤ z ↔ 0 ≤ z.re ∧ z.im = 0`), which is strictly stronger than the
real-part statements; real-part corollaries are provided.

## Tagging convention

Every lemma in this development is tagged in its doc comment as either
**algebra-only** (a statement about the `*`-algebra of events, consuming no
state or trace functional) or **trace-dependent** (a statement
whose content passes through the trace pairing `(ρ, P) ↦ Tr(ρ P)`).

This module is deliberately free of any interpretational vocabulary: it is
plain finite-dimensional operator algebra.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {n : ℕ}

/-- An **event** of the finite algebra `Matrix (Fin n) (Fin n) ℂ`:
a Hermitian idempotent, i.e. an orthogonal projection. -/
def IsEvent (P : Matrix (Fin n) (Fin n) ℂ) : Prop :=
  P.IsHermitian ∧ P * P = P

/-- A **state** on the finite algebra `Matrix (Fin n) (Fin n) ℂ`:
a positive-semidefinite matrix of trace one (a density matrix). -/
def IsState (ρ : Matrix (Fin n) (Fin n) ℂ) : Prop :=
  ρ.PosSemidef ∧ ρ.trace = 1

/-- A state-facing subtype. APIs returning this type certify normalization
and positivity in their result rather than returning an arbitrary matrix. -/
abbrev StateMatrix (n : ℕ) :=
  {ρ : Matrix (Fin n) (Fin n) ℂ // IsState ρ}

/-- A projection-event subtype used by semantically typed measurement APIs. -/
abbrev ProjectionEvent (n : ℕ) :=
  {P : Matrix (Fin n) (Fin n) ℂ // IsEvent P}

/-- The **Born weight** of the event `P` under the state `ρ`: the trace
pairing `Tr(ρ P)`. -/
noncomputable def bornWeight (ρ P : Matrix (Fin n) (Fin n) ℂ) : ℂ :=
  (ρ * P).trace

/-! ## Events: closure properties (algebra-only) -/

/-- **Algebra-only.** The zero matrix is an event (the impossible event). -/
theorem isEvent_zero : IsEvent (0 : Matrix (Fin n) (Fin n) ℂ) :=
  ⟨isHermitian_zero, by simp⟩

/-- **Algebra-only.** The identity matrix is an event (the sure event). -/
theorem isEvent_one : IsEvent (1 : Matrix (Fin n) (Fin n) ℂ) :=
  ⟨isHermitian_one, one_mul 1⟩

/-- **Algebra-only.** The complement `1 - P` of an event is an event. -/
theorem IsEvent.compl {P : Matrix (Fin n) (Fin n) ℂ} (hP : IsEvent P) :
    IsEvent (1 - P) := by
  refine ⟨isHermitian_one.sub hP.1, ?_⟩
  rw [sub_mul, one_mul, mul_sub, mul_one, hP.2]
  abel

/-- **Algebra-only.** Orthogonality of events is symmetric: if `P * Q = 0`
for Hermitian idempotents `P`, `Q`, then `Q * P = 0` (take the conjugate
transpose of `P * Q = 0`). -/
theorem IsEvent.orthogonal_symm {P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (h : P * Q = 0) : Q * P = 0 := by
  have hstar : (P * Q)ᴴ = Q * P := by
    rw [conjTranspose_mul, hP.1.eq, hQ.1.eq]
  rw [← hstar, h, conjTranspose_zero]

/-- **Algebra-only.** The sum of two orthogonal events is an event. -/
theorem IsEvent.add {P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (h : P * Q = 0) : IsEvent (P + Q) := by
  refine ⟨hP.1.add hQ.1, ?_⟩
  rw [add_mul, mul_add, mul_add, hP.2, hQ.2, h, hP.orthogonal_symm hQ h]
  abel

/-- **Algebra-only.** The product of two commuting events is an event
(the conjunction of compatible events). -/
theorem IsEvent.mul_of_commute {P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (h : P * Q = Q * P) :
    IsEvent (P * Q) := by
  constructor
  · show (P * Q)ᴴ = P * Q
    rw [conjTranspose_mul, hP.1.eq, hQ.1.eq, ← h]
  · calc P * Q * (P * Q) = P * (Q * P) * Q := by
          simp only [mul_assoc]
      _ = P * P * (Q * Q) := by rw [← h]; simp only [mul_assoc]
      _ = P * Q := by rw [hP.2, hQ.2]

/-- **Algebra-only.** If `P` is a subevent of `Q` (in the projection order,
witnessed by `P * Q = P`), then `Q` absorbs `P` on the other side as well:
`Q * P = P`. -/
theorem IsEvent.absorb_of_le {P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (hle : P * Q = P) : Q * P = P := by
  have hstar : (P * Q)ᴴ = Q * P := by
    rw [conjTranspose_mul, hP.1.eq, hQ.1.eq]
  rw [← hstar, hle, hP.1.eq]

/-- **Algebra-only.** The difference `Q - P` of an event `Q` and a subevent
`P` (witnessed by `P * Q = P`) is an event. -/
theorem IsEvent.sub_of_le {P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (hle : P * Q = P) :
    IsEvent (Q - P) := by
  refine ⟨hQ.1.sub hP.1, ?_⟩
  rw [sub_mul, mul_sub, mul_sub, hQ.2, hP.2, hle, hP.absorb_of_le hQ hle]
  abel

/-- **Algebra-only.** Every event is positive semidefinite:
`P = Pᴴ * P` exhibits it as a Gram matrix. -/
theorem IsEvent.posSemidef {P : Matrix (Fin n) (Fin n) ℂ} (hP : IsEvent P) :
    P.PosSemidef := by
  have : P = Pᴴ * P := by rw [hP.1.eq, hP.2]
  rw [this]
  exact posSemidef_conjTranspose_mul_self P

/-- **Algebra-only.** The dichotomic observable of an event: `1 - 2P` is a
selfadjoint involution (eigenvalues `±1`, answering the yes/no question
`P`). This is the bridge from events to the CHSH-type hypotheses of
`EventAlgebra.Tsirelson`. -/
theorem IsEvent.one_sub_two_smul_involution {P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) :
    IsSelfAdjoint (1 - 2 • P) ∧ (1 - 2 • P) * (1 - 2 • P) = 1 := by
  constructor
  · have hPsa : IsSelfAdjoint P := isHermitian_iff_isSelfAdjoint.mp hP.1
    rw [two_smul]
    exact (IsSelfAdjoint.one _).sub (hPsa.add hPsa)
  · rw [two_smul]
    simp only [mul_sub, sub_mul, add_mul, mul_add, mul_one, one_mul, hP.2]
    abel

/-! ## Born weights (trace-dependent) -/

/-- **Trace-dependent.** The Born weight is linear in the event
argument: additivity under sums (no orthogonality needed for the bare trace
identity). -/
theorem bornWeight_add (ρ P Q : Matrix (Fin n) (Fin n) ℂ) :
    bornWeight ρ (P + Q) = bornWeight ρ P + bornWeight ρ Q := by
  simp only [bornWeight, mul_add, trace_add]

/-- **Trace-dependent.** The Born weight is subtractive in the
event argument. -/
theorem bornWeight_sub (ρ P Q : Matrix (Fin n) (Fin n) ℂ) :
    bornWeight ρ (P - Q) = bornWeight ρ P - bornWeight ρ Q := by
  simp only [bornWeight, mul_sub, trace_sub]

/-- **Trace-dependent.** The Born weight commutes with finite sums
of events. -/
theorem bornWeight_sum {ι : Type*} (ρ : Matrix (Fin n) (Fin n) ℂ)
    (s : Finset ι) (P : ι → Matrix (Fin n) (Fin n) ℂ) :
    bornWeight ρ (∑ i ∈ s, P i) = ∑ i ∈ s, bornWeight ρ (P i) := by
  simp only [bornWeight, Finset.mul_sum, trace_sum]

/-- **Trace-dependent.** Scaling the state scales the Born weight. -/
theorem bornWeight_smul (c : ℂ) (ρ P : Matrix (Fin n) (Fin n) ℂ) :
    bornWeight (c • ρ) P = c * bornWeight ρ P := by
  simp only [bornWeight, smul_mul_assoc, trace_smul, smul_eq_mul]

/-- **Trace-dependent.** Sandwich identity: for an idempotent `P`,
`Tr(P ρ P) = Tr(ρ P)`; the compressed and the paired forms of the weight
agree. -/
theorem trace_sandwich {P : Matrix (Fin n) (Fin n) ℂ} (hP : P * P = P)
    (ρ : Matrix (Fin n) (Fin n) ℂ) :
    (P * ρ * P).trace = bornWeight ρ P := by
  rw [bornWeight, trace_mul_cycle, hP, trace_mul_comm]

/-- **Trace-dependent.** Reality of the Born weight: for a
Hermitian state matrix and a Hermitian event matrix, the weight is fixed by
complex conjugation. (Only hermiticity is used; neither positivity nor
normalisation.) -/
theorem star_bornWeight {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.IsHermitian) (hP : P.IsHermitian) :
    star (bornWeight ρ P) = bornWeight ρ P := by
  rw [bornWeight, ← trace_conjTranspose, conjTranspose_mul, hρ.eq, hP.eq,
    trace_mul_comm]

/-- **Trace-dependent.** The Born weight equals its own real part;
the packaged form of reality. -/
theorem bornWeight_eq_re {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.IsHermitian) (hP : P.IsHermitian) :
    (((bornWeight ρ P).re : ℝ) : ℂ) = bornWeight ρ P :=
  Complex.conj_eq_iff_re.mp (star_bornWeight hρ hP)

/-- **Trace-dependent.** Nonnegativity of the Born weight, in the
partial order of `ℂ`: for `ρ` positive semidefinite and `P` an event,
`0 ≤ Tr(ρ P)`. The proof compresses the state: `Tr(ρ P) = Tr(P ρ P)` and
`P ρ P = P ρ Pᴴ` is positive semidefinite. -/
theorem bornWeight_nonneg {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hP : IsEvent P) : 0 ≤ bornWeight ρ P := by
  rw [← trace_sandwich hP.2]
  have hpsd : (P * ρ * Pᴴ).PosSemidef := hρ.mul_mul_conjTranspose_same P
  rw [hP.1.eq] at hpsd
  exact hpsd.trace_nonneg

/-- **Trace-dependent.** Real-part form of nonnegativity. -/
theorem bornWeight_re_nonneg {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hP : IsEvent P) : 0 ≤ (bornWeight ρ P).re := by
  simpa using (Complex.le_def.mp (bornWeight_nonneg hρ hP)).1

/-- **Trace-dependent.** Normalisation: the sure event has Born
weight one under every state. -/
theorem bornWeight_one {ρ : Matrix (Fin n) (Fin n) ℂ} (hρ : IsState ρ) :
    bornWeight ρ 1 = 1 := by
  rw [bornWeight, mul_one, hρ.2]

/-- **Trace-dependent.** Orthogonal events have an event-valued sum and
additive Born weight. Unlike bare trace linearity, both conclusions usefully
package the projection-event hypotheses. -/
theorem isEvent_add_and_bornWeight_add_of_orthogonal
    {ρ P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (hPQ : P * Q = 0) :
    IsEvent (P + Q) ∧
      bornWeight ρ (P + Q) = bornWeight ρ P + bornWeight ρ Q :=
  ⟨hP.add hQ hPQ, bornWeight_add ρ P Q⟩

/-- **Trace-dependent.** Upper bound: the Born weight of any event
under a state is at most `1`, in the partial order of `ℂ` (via the
complement event `1 - P`). -/
theorem bornWeight_le_one {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hP : IsEvent P) : bornWeight ρ P ≤ 1 := by
  have hcompl : 0 ≤ bornWeight ρ (1 - P) := bornWeight_nonneg hρ.1 hP.compl
  rw [bornWeight_sub, bornWeight_one hρ] at hcompl
  exact sub_nonneg.mp hcompl

/-- **Trace-dependent.** Real-part form of the upper bound. -/
theorem bornWeight_re_le_one {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hP : IsEvent P) : (bornWeight ρ P).re ≤ 1 := by
  simpa using (Complex.le_def.mp (bornWeight_le_one hρ hP)).1

/-- **Trace-dependent.** Monotonicity under subevents: if
`P * Q = P` (i.e. `P ≤ Q` in the projection order), then
`Tr(ρ P) ≤ Tr(ρ Q)` in the partial order of `ℂ`. -/
theorem bornWeight_mono {ρ P Q : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hP : IsEvent P) (hQ : IsEvent Q)
    (hle : P * Q = P) : bornWeight ρ P ≤ bornWeight ρ Q := by
  have hdiff : 0 ≤ bornWeight ρ (Q - P) :=
    bornWeight_nonneg hρ (hP.sub_of_le hQ hle)
  rw [bornWeight_sub] at hdiff
  exact sub_nonneg.mp hdiff

/-- **Trace-dependent.** A Born weight that is nonzero has strictly
positive real part (given a positive-semidefinite state matrix and an
event); conversely a weight with positive real part is nonzero. -/
theorem bornWeight_ne_zero_iff_re_pos {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hP : IsEvent P) :
    bornWeight ρ P ≠ 0 ↔ 0 < (bornWeight ρ P).re := by
  constructor
  · intro hw
    have hpos : 0 < bornWeight ρ P :=
      lt_of_le_of_ne (bornWeight_nonneg hρ hP) (Ne.symm hw)
    simpa using (Complex.lt_def.mp hpos).1
  · intro hre h0
    rw [h0] at hre
    simp at hre

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms isEvent_zero
#print axioms isEvent_one
#print axioms IsEvent.compl
#print axioms IsEvent.orthogonal_symm
#print axioms IsEvent.add
#print axioms IsEvent.mul_of_commute
#print axioms IsEvent.absorb_of_le
#print axioms IsEvent.sub_of_le
#print axioms IsEvent.posSemidef
#print axioms IsEvent.one_sub_two_smul_involution
#print axioms bornWeight_add
#print axioms bornWeight_sub
#print axioms bornWeight_sum
#print axioms bornWeight_smul
#print axioms trace_sandwich
#print axioms star_bornWeight
#print axioms bornWeight_eq_re
#print axioms bornWeight_nonneg
#print axioms bornWeight_re_nonneg
#print axioms bornWeight_one
#print axioms isEvent_add_and_bornWeight_add_of_orthogonal
#print axioms bornWeight_le_one
#print axioms bornWeight_re_le_one
#print axioms bornWeight_mono
#print axioms bornWeight_ne_zero_iff_re_pos

end EventAlgebra
