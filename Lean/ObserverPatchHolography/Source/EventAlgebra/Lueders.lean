import EventAlgebra.Basic

/-!
# Lüders conditioning

The Lüders update rule for conditioning a state on an event:

* `EventAlgebra.luedersUpdate` — `ρ ↦ (Tr(ρ P))⁻¹ • (P ρ P)`;
* the update is again a state whenever the Born weight is nonzero
  (equivalently: has strictly positive real part);
* **repeatability** — the conditioned state assigns weight `1` to the
  conditioning event;
* **idempotence** — conditioning twice on the same event equals
  conditioning once;
* **compatibility** — for commuting events, sequential conditioning
  composes to conditioning on the product event, hence is
  order-exchangeable;
* the **classical restriction** — when the state commutes with the event,
  the update is the normalised restriction `(Tr(ρ P))⁻¹ • (ρ P)`.

## Tagging convention

As in `EventAlgebra.Basic`: each lemma is tagged **algebra-only** (pure
`*`-algebra content) or **trace-dependent** (content passing
through the trace pairing). Everything about `luedersUpdate` consumes the
trace through its normalising Born weight.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {n : ℕ}

/-- The **Lüders update** of a state `ρ` by an event `P`: the compressed
state `P ρ P`, renormalised by the Born weight. When the Born weight is
zero (conditioning on an almost-surely-false event) the inverse is zero by
convention and the update degenerates to the zero matrix; all substantive
lemmas below guard against this with `bornWeight ρ P ≠ 0`. -/
noncomputable def luedersUpdate (ρ P : Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ :=
  (bornWeight ρ P)⁻¹ • (P * ρ * P)

/-- **Trace-dependent.** The compressed matrix `P ρ P` is positive
semidefinite, and so is its renormalisation by the (nonnegative real) Born
weight; no nonvanishing guard is needed because the degenerate update is
the zero matrix. -/
theorem luedersUpdate_posSemidef {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hP : IsEvent P) :
    (luedersUpdate ρ P).PosSemidef := by
  have hsand : (P * ρ * P).PosSemidef := by
    have := hρ.mul_mul_conjTranspose_same P
    rwa [hP.1.eq] at this
  rcases eq_or_ne (bornWeight ρ P) 0 with hw | hw
  · rw [luedersUpdate, hw, _root_.inv_zero, zero_smul]
    exact Matrix.PosSemidef.zero
  · have hpos : 0 < bornWeight ρ P :=
      lt_of_le_of_ne (bornWeight_nonneg hρ hP) (Ne.symm hw)
    have hinv : 0 ≤ (bornWeight ρ P)⁻¹ := (RCLike.inv_pos.mpr hpos).le
    exact hsand.smul hinv

/-- **Trace-dependent.** The Lüders update has trace one whenever
the Born weight does not vanish. -/
theorem trace_luedersUpdate {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    (luedersUpdate ρ P).trace = 1 := by
  rw [luedersUpdate, trace_smul, trace_sandwich hP.2, smul_eq_mul,
    inv_mul_cancel₀ hw]

/-- **Trace-dependent.** The Lüders update of a state by an event
of nonzero Born weight is a state. Via
`bornWeight_ne_zero_iff_re_pos`, the guard is equivalent to the weight
having strictly positive real part. -/
theorem luedersUpdate_isState {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    IsState (luedersUpdate ρ P) :=
  ⟨luedersUpdate_posSemidef hρ.1 hP, trace_luedersUpdate hP hw⟩

/-- **Trace-dependent.** **Repeatability**: after conditioning on
`P`, the event `P` holds with Born weight one. -/
theorem bornWeight_luedersUpdate_self {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    bornWeight (luedersUpdate ρ P) P = 1 := by
  have habsorb : P * ρ * P * P = P * ρ * P := by
    rw [mul_assoc (P * ρ), hP.2]
  have hkey : bornWeight (P * ρ * P) P = bornWeight ρ P := by
    rw [bornWeight, habsorb, trace_sandwich hP.2]
  rw [luedersUpdate, bornWeight_smul, hkey, inv_mul_cancel₀ hw]

/-- **Trace-dependent.** **Idempotence**: conditioning twice on the
same event gives the same state as conditioning once. -/
theorem luedersUpdate_idem {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    luedersUpdate (luedersUpdate ρ P) P = luedersUpdate ρ P := by
  have hcompress : P * (P * ρ * P) * P = P * ρ * P := by
    rw [← mul_assoc, ← mul_assoc, hP.2, mul_assoc (P * ρ), hP.2]
  rw [luedersUpdate, bornWeight_luedersUpdate_self hP hw, inv_one, one_smul,
    luedersUpdate, mul_smul_comm, smul_mul_assoc, hcompress]

/-- **Trace-dependent.** **Compatibility**: for commuting events,
sequential conditioning composes to conditioning on the product event
`P * Q` (which is an event by `IsEvent.mul_of_commute`). Only the weight of
the first conditioning needs a nonvanishing guard: if the joint weight
vanishes, both sides degenerate to zero together. -/
theorem luedersUpdate_luedersUpdate_of_commute
    {ρ P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (_hQ : IsEvent Q) (hc : P * Q = Q * P)
    (hw : bornWeight ρ P ≠ 0) :
    luedersUpdate (luedersUpdate ρ P) Q = luedersUpdate ρ (P * Q) := by
  -- The compressed matrix: `Q (P ρ P) Q = (P Q) ρ (P Q)`.
  have hmat : Q * (P * ρ * P) * Q = (P * Q) * ρ * (P * Q) := by
    calc Q * (P * ρ * P) * Q = ((Q * P) * ρ) * (P * Q) := by
          simp only [mul_assoc]
      _ = (P * Q) * ρ * (P * Q) := by rw [← hc]
  -- The intermediate weight: `Tr(σ Q) = w_P⁻¹ · Tr(ρ (P Q))`.
  have hweight : bornWeight (luedersUpdate ρ P) Q =
      (bornWeight ρ P)⁻¹ * bornWeight ρ (P * Q) := by
    have htr : (P * ρ * P) * Q = (P * ρ) * (P * Q) := by
      simp only [mul_assoc]
    have hPQP : (P * Q) * (P * ρ) = (P * Q) * P * ρ := by
      simp only [mul_assoc]
    have hcollapse : (P * Q) * P = P * Q := by
      rw [hc, mul_assoc, hP.2]
    have hkey : bornWeight (P * ρ * P) Q = bornWeight ρ (P * Q) := by
      rw [bornWeight, bornWeight, htr, trace_mul_comm (P * ρ) (P * Q),
        hPQP, hcollapse, trace_mul_comm (P * Q) ρ]
    rw [luedersUpdate, bornWeight_smul, hkey]
  rw [luedersUpdate, hweight, luedersUpdate, mul_smul_comm, smul_mul_assoc,
    hmat, smul_smul, luedersUpdate]
  congr 1
  rw [mul_inv, inv_inv, mul_comm (bornWeight ρ P), mul_assoc,
    mul_inv_cancel₀ hw, mul_one]

/-- **Trace-dependent.** **Order exchange**: commuting events may
be conditioned on in either order. -/
theorem luedersUpdate_comm {ρ P Q : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (hc : P * Q = Q * P)
    (hwP : bornWeight ρ P ≠ 0) (hwQ : bornWeight ρ Q ≠ 0) :
    luedersUpdate (luedersUpdate ρ P) Q =
      luedersUpdate (luedersUpdate ρ Q) P := by
  rw [luedersUpdate_luedersUpdate_of_commute hP hQ hc hwP,
    luedersUpdate_luedersUpdate_of_commute hQ hP hc.symm hwQ, hc]

/-- **Trace-dependent.** The **classical restriction**: when the
state commutes with the event, the Lüders update is the normalised
restriction `(Tr(ρ P))⁻¹ • (ρ P)` — conditioning collapses to classical
conditional probability. -/
theorem luedersUpdate_of_commute {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) (hc : ρ * P = P * ρ) :
    luedersUpdate ρ P = (bornWeight ρ P)⁻¹ • (ρ * P) := by
  rw [luedersUpdate, ← hc, mul_assoc, hP.2]

/-!
## Conditioning as one-step convergence to a fixed point

Repeatability and idempotence assemble into a fixed-point description of
measurement.  For an event `P`, call a state **certain of `P`** when it
assigns `P` Born weight `1` (`certainStates`).  Then:

* conditioning any state of nonvanishing weight lands in `certainStates P`
  in a single step (`luedersUpdate_mem_certainStates`);
* on `certainStates P` conditioning is the identity
  (`luedersUpdate_eq_self_of_mem_certainStates`), so every element is a
  fixed point of the update map; and
* conversely, among states of nonvanishing weight the fixed points of
  conditioning on `P` are *exactly* the states certain of `P`
  (`luedersUpdate_eq_self_iff`).

Conditioning on `P` is thus a retraction of the states of nonvanishing
weight onto its own fixed-point set `certainStates P`: measurement update
is convergence to a fixed point, reached after one step.
-/

/-- The **certainty set** of an event `P`: the states that assign `P` Born
weight `1`. -/
def certainStates (P : Matrix (Fin n) (Fin n) ℂ) :
    Set (Matrix (Fin n) (Fin n) ℂ) :=
  {σ | IsState σ ∧ bornWeight σ P = 1}

/-- **Trace-dependent.** One-step convergence: conditioning any
state of nonvanishing Born weight on `P` lands in the certainty set of `P`.
This is `luedersUpdate_isState` and repeatability, packaged. -/
theorem luedersUpdate_mem_certainStates {ρ P : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    luedersUpdate ρ P ∈ certainStates P :=
  ⟨luedersUpdate_isState hρ hP hw, bornWeight_luedersUpdate_self hP hw⟩

/-- **Trace-dependent.** A state that is certain of `P` is
supported on `P`: the state matrix absorbs the event on either side.  The
proof compresses the state by the complement `1 - P`, whose Born weight
vanishes; positive semidefiniteness of the zero-trace compression
`(1-P) σ (1-P)` forces it to vanish, and Cauchy–Schwarz for the state's
quadratic form then kills `σ (1-P)` itself. -/
theorem mul_eq_self_of_bornWeight_one {σ P : Matrix (Fin n) (Fin n) ℂ}
    (hσ : IsState σ) (hP : IsEvent P) (h1 : bornWeight σ P = 1) :
    σ * P = σ ∧ P * σ = σ := by
  set Q : Matrix (Fin n) (Fin n) ℂ := 1 - P with hQdef
  have hQ : IsEvent Q := hP.compl
  -- The complement has vanishing Born weight.
  have hcw : bornWeight σ Q = 0 := by
    rw [hQdef, bornWeight_sub, bornWeight_one hσ, h1, sub_self]
  -- Hence the compression of `σ` by the complement vanishes.
  have hpsd : (Q * σ * Q).PosSemidef := by
    have := hσ.1.mul_mul_conjTranspose_same Q
    rwa [hQ.1.eq] at this
  have hzero : Q * σ * Q = 0 := by
    rw [← hpsd.trace_eq_zero_iff, trace_sandwich hQ.2, hcw]
  -- Cauchy–Schwarz for the quadratic form of `σ`: `σ (Q x) = 0` for all `x`.
  have hvec : ∀ x : Fin n → ℂ, σ *ᵥ (Q *ᵥ x) = 0 := by
    intro x
    apply (hσ.1.dotProduct_mulVec_zero_iff (Q *ᵥ x)).mp
    calc star (Q *ᵥ x) ⬝ᵥ σ *ᵥ (Q *ᵥ x)
        = (star x ᵥ* Q) ⬝ᵥ (σ *ᵥ (Q *ᵥ x)) := by rw [star_mulVec, hQ.1.eq]
      _ = (star x ᵥ* Q) ⬝ᵥ ((σ * Q) *ᵥ x) := by rw [mulVec_mulVec]
      _ = star x ⬝ᵥ (Q *ᵥ ((σ * Q) *ᵥ x)) := (dotProduct_mulVec _ _ _).symm
      _ = star x ⬝ᵥ ((Q * σ * Q) *ᵥ x) := by rw [mulVec_mulVec, ← mul_assoc]
      _ = 0 := by rw [hzero, zero_mulVec, dotProduct_zero]
  have hmul : σ * Q = 0 := by
    apply ext_of_mulVec_single (fun i => ?_)
    rw [← mulVec_mulVec, hvec, zero_mulVec]
  -- Unfold the complement.
  have hright : σ * P = σ := by
    rw [hQdef, mul_sub, mul_one] at hmul
    exact (sub_eq_zero.mp hmul).symm
  refine ⟨hright, ?_⟩
  have h' := congrArg conjTranspose hright
  rwa [conjTranspose_mul, hσ.1.isHermitian.eq, hP.1.eq] at h'

/-- **Trace-dependent.** On its certainty set, conditioning on `P`
acts as the identity: every state certain of `P` is a fixed point of the
Lüders update. -/
theorem luedersUpdate_eq_self_of_mem_certainStates
    {σ P : Matrix (Fin n) (Fin n) ℂ} (hP : IsEvent P)
    (hσ : σ ∈ certainStates P) : luedersUpdate σ P = σ := by
  obtain ⟨hstate, h1⟩ := hσ
  obtain ⟨hright, hleft⟩ := mul_eq_self_of_bornWeight_one hstate hP h1
  rw [luedersUpdate, h1, inv_one, one_smul, hleft, hright]

/-- **Trace-dependent.** The fixed-point characterisation of
measurement: among states assigning `P` nonzero weight, the fixed points of
conditioning on `P` are exactly the states certain of `P`.  Together with
`luedersUpdate_mem_certainStates`, conditioning is a one-step retraction of
the states of nonvanishing weight onto its fixed-point set. -/
theorem luedersUpdate_eq_self_iff {σ P : Matrix (Fin n) (Fin n) ℂ}
    (hσ : IsState σ) (hP : IsEvent P) (hw : bornWeight σ P ≠ 0) :
    luedersUpdate σ P = σ ↔ σ ∈ certainStates P := by
  constructor
  · intro hfix
    refine ⟨hσ, ?_⟩
    have := bornWeight_luedersUpdate_self hP hw
    rwa [hfix] at this
  · exact luedersUpdate_eq_self_of_mem_certainStates hP

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms luedersUpdate_posSemidef
#print axioms trace_luedersUpdate
#print axioms luedersUpdate_isState
#print axioms bornWeight_luedersUpdate_self
#print axioms luedersUpdate_idem
#print axioms luedersUpdate_luedersUpdate_of_commute
#print axioms luedersUpdate_comm
#print axioms luedersUpdate_of_commute
#print axioms luedersUpdate_mem_certainStates
#print axioms mul_eq_self_of_bornWeight_one
#print axioms luedersUpdate_eq_self_of_mem_certainStates
#print axioms luedersUpdate_eq_self_iff

end EventAlgebra
