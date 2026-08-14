import EventAlgebra.FiniteBuschGleason
import EventAlgebra.Lueders
import EventAlgebra.Tsirelson
import EventAlgebra.PublicRecordAlgebra
import EventAlgebra.SlotLocalitySurface
import Dynamics.StoneConverse
import QFT.SourcePhaseLiftBridge

set_option autoImplicit false

namespace EventAlgebra

/-!
# The composed quantum adequacy surface

One parameterized package for the observation-ledger rows OL-C1 (Born
probabilities and Lueders conditioning), OL-C2 (Schroedinger dynamics as
the unique continuous symmetry flow), OL-C3 (entanglement to the
Tsirelson bound), and OL-C5 (phase completion), V3 issue #730. Every
theorem below is a committed result of the event-algebra, dynamics, and
phase-lift modules, re-exported through named premise bundles so that a
ledger row cites one surface instead of scattered conditionals. Each
proof is a direct application of the cited module; this file adds no
derivation.

`RepresentedPublicFrame` bundles the register rows PR-02 (the
algebra-state representation of public records) and PR-03 (operational
effect additivity). `PhaseInstrument` bundles the register row PR-04
(the phase-operation instrument); the committed algebraic lift inhabits
it (`committedPhaseInstrument`), and the real/Kraus blindness boundary
is carried alongside as the committed limit.

The companion modules `EventAlgebra.SlotLocalitySurface` (the OL-C3
and OL-C4 receipts composed under the PR-44 slot split) and
`EventAlgebra.SchroedingerFrameFlow` (the OL-C1 single-statement
composition with conditioning closure, and the OL-C2 duality under
PR-43) extend this surface inside the same `QuantumSurface` namespace.

Composition boundaries, stated exactly:

* the Tsirelson bound (`tsirelson_of_events`) is algebra-only: it
  consumes four events of the represented algebra and no frame
  valuation, so it does not factor through `RepresentedPublicFrame`;
  it is re-exported side by side inside the namespace, on the same
  represented carrier, and the bound is an inequality with no
  attainment claim; the composed bipartite form, with the cross-party
  commutation derived from the registered PR-44 slot split, is
  `slot_locality_receipts` on the slot-locality surface;
* the unique continuous symmetry flow (`unique_continuous_flow`)
  consumes a pointwise continuous star-automorphism group of the
  represented algebra and no frame valuation; the represented carrier
  of PR-02 is the composition point, the flow concerns one full matrix
  block, and no physical clock or source-selected Hamiltonian is
  claimed, the parameter is a real number; the composed form on the
  represented state is `schroedinger_frame_duality` in
  `EventAlgebra.SchroedingerFrameFlow`: the unique Busch-Gleason state
  of the Heisenberg-shifted frame follows the propagator conjugation
  of the unique generator, with the flow staying the supplied PR-43
  premise and the clock attachment staying open on PR-15;
* the phase-completion theorems consume the committed dimension-two
  payload and no frame valuation; the register row PR-04 records that
  the payload supplies the algebraic lift and no instrument receipts,
  and `phase_boundary_summary` carries the committed limit: realized
  outcome statistics inhabit only the diagonal native context;
* OL-C4 (no-signalling and locality receipts) is composed on the
  slot-locality surface: `slot_locality_receipts` consumes the PR-44
  interface, derives the Tsirelson commutation from the slot split,
  and threads the lifted local channel into the PR-02 trace pairing;
  the physical spacelike reading stays open on register row PR-52;
  OL-C6 (structural field theory) is owed to the regional-net
  modules.

The row dispositions are the register's: PR-02 is an explicit
representation assumption, while PR-03 and PR-04 identify assumptions
that a source-complete account would have to replace by constructions.
Nothing here claims the axioms force the quantum representation of
public records or supply the declared phase instrument.
-/

/-- **Register rows PR-02 and PR-03: the represented public frame.**
PR-02, the algebra-state representation of public records, enters as
the carrier of this structure: record effects are effects of the finite
complex matrix star algebra `Matrix (Fin d) (Fin d) ℂ`, and outcome
weights are read on that declared surface
(`Lean/EventAlgebra/PublicRecordAlgebra.lean` packages the record span;
`recordSurfaceClassical` below re-exports its exact finite classical
form). PR-03, operational effect additivity, enters as the field
`additive`: the frame valuation is additive on every coexistent effect
sum across measurement contexts. The range and normalisation fields
state the probability reading of the frame. The committed
countermodels, the transverse cubic response of `FiniteWebBornNoGo` on
the retained finite unsharp web and the trine battery of
`InterlockingContexts`, are the record of why `additive` is a register
row rather than a theorem: range, per-context normalisation, and
continuity without cross-context additivity admit non-Born valuations.
No theorem in this module asserts a source selection of the
representation or of the additivity. -/
structure RepresentedPublicFrame (d : ℕ) where
  /-- The frame valuation: the outcome-weight readout on the represented
  record effects. -/
  value : Matrix (Fin d) (Fin d) ℂ → ℝ
  /-- Effects receive nonnegative weights. -/
  nonneg : ∀ {E : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → 0 ≤ value E
  /-- Effects receive weights at most one. -/
  le_one : ∀ {E : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → value E ≤ 1
  /-- The sure effect receives weight one. -/
  map_one : value 1 = 1
  /-- **Register row PR-03.** Additivity across contexts: whenever `E`,
  `F`, and `E + F` are all effects, the weight of the sum is the sum of
  the weights. -/
  additive : ∀ {E F : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → IsEffect F →
    IsEffect (E + F) → value (E + F) = value E + value F

/-- **Register row PR-04: the phase-operation instrument.** The declared
phase operation on the dimension-two represented algebra: one projection
event that separates the two Pauli-Y states which every phase-free real
record context weighs equally, and that closes operator tomography
together with the committed record projector and the algebraically
rotated projector. The committed payload supplies the exact algebraic
lift `I/2 - (2*sqrt(3)/3) i (QP - PQ)` and no instrument receipts: no
rotated or phase outcome counts and no common-preparation validation.
The lift inhabits every field of this structure
(`committedPhaseInstrument`); the instrument character of the row is the
declared part, and no theorem in this module asserts it. -/
structure PhaseInstrument where
  /-- The declared phase operation. -/
  lift : Matrix (Fin 2) (Fin 2) ℂ
  /-- The declared operation is a projection event of the represented
  algebra. -/
  isEvent : IsEvent lift
  /-- The declared operation separates the two Pauli-Y states that the
  phase-free real closure cannot distinguish. -/
  separates : bornWeight OPH.QFT.rhoYPlus lift = 1 ∧
    bornWeight OPH.QFT.rhoYMinus lift = 0
  /-- Tomographic completion: together with the committed record
  projector and the committed rotated projector, the declared operation
  determines every matrix on each fixed-trace slice. -/
  completes : ∀ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
    Matrix.trace rho = Matrix.trace sigma →
    bornWeight rho (OPH.QFT.complexifyRealMatrix OPH.QFT.recordProjector) =
      bornWeight sigma
        (OPH.QFT.complexifyRealMatrix OPH.QFT.recordProjector) →
    bornWeight rho
        (OPH.QFT.complexifyRealMatrix (OPH.QFT.conjProjector 3)) =
      bornWeight sigma
        (OPH.QFT.complexifyRealMatrix (OPH.QFT.conjProjector 3)) →
    bornWeight rho lift = bornWeight sigma lift →
    rho = sigma

namespace QuantumSurface

variable {d : ℕ} (F : RepresentedPublicFrame d)

/-! ## The frame is an additive effect valuation

The premise bundle repackages the hypotheses of the finite
Busch-Gleason theorem; the repackaging is definitional. -/

/-- The represented public frame is an additive effect valuation in the
sense of `FiniteBuschGleason`. -/
theorem frame_valuation : IsEffectValuation F.value :=
  ⟨F.nonneg, F.le_one, F.map_one, F.additive⟩

/-! ## OL-C1: Born probabilities on the represented surface -/

/-- The state represented by the frame: the unique density matrix of the
finite Busch-Gleason theorem. -/
noncomputable def frameState : Matrix (Fin d) (Fin d) ℂ :=
  gleasonState F.value

/-- The represented matrix is a state: positive semidefinite with unit
trace. -/
theorem frameState_isState : IsState (frameState F) :=
  gleasonState_isState (frame_valuation F)

/-- **Born probabilities.** Under PR-02 and PR-03 the frame weight of
every represented effect is the Born weight of the frame state,
`v(E) = Re Tr(rho E)`. -/
theorem born_represents {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    F.value E = (bornWeight (frameState F) E).re :=
  gleasonState_represents (frame_valuation F) hE

/-- **The finite Busch-Gleason representation on the frame.** The frame
determines its density matrix uniquely: exactly one state reproduces
every frame weight in Born form. -/
theorem busch_gleason :
    ∃! ρ : Matrix (Fin d) (Fin d) ℂ,
      IsState ρ ∧ ∀ E, IsEffect E → F.value E = (bornWeight ρ E).re :=
  finite_busch_gleason (frame_valuation F)

/-- **Dimension-two coverage.** The representation theorem instantiated
on a dimension-two frame. The effect-level additivity row covers the
dimension excluded by the projective Gleason hypothesis: with sharp
projectors alone, dimension two admits non-Born frame functions, and
the committed finite-web countermodels (`FiniteWebBornNoGo`,
`InterlockingContexts`) exhibit such a valuation on the retained unsharp
web; with PR-03 on all coexistent effect sums the Born form is forced in
every finite dimension including two. -/
theorem busch_gleason_dimension_two (F₂ : RepresentedPublicFrame 2) :
    ∃! ρ : Matrix (Fin 2) (Fin 2) ℂ,
      IsState ρ ∧ ∀ E, IsEffect E → F₂.value E = (bornWeight ρ E).re :=
  busch_gleason F₂

/-- **The exact boundary of PR-03.** Side-by-side re-export, consuming
no frame: the additive effect valuations on the represented algebra are
precisely the Born functionals of density matrices. The forward
direction is the representation theorem; the converse is the reason the
row is a premise about the valuation and nothing stronger. -/
theorem additive_frame_iff_born {v : Matrix (Fin d) (Fin d) ℂ → ℝ} :
    IsEffectValuation v ↔
      ∃ ρ, IsState ρ ∧ ∀ E, IsEffect E → v E = (bornWeight ρ E).re :=
  isEffectValuation_iff_born

/-! ## PR-02 receipt: the record surface is finite classical

The representation row declares public records as the span of the
orthogonal projectors of a projective partition. The re-exports below
consume a partition of the represented algebra; the first consumes no
frame. -/

/-- The public record algebra of a projective partition is star-algebra
equivalent to complex functions on the active record labels: the
declared record surface of PR-02 is an exact finite classical outcome
algebra. Algebra-only; consumes no frame valuation. -/
noncomputable def recordSurfaceClassical {n k : ℕ}
    (part : ProjectivePartition n k) :
    part.publicSubalgebra ≃⋆ₐ[ℂ] (part.ActiveIndex → ℂ) :=
  part.publicRecordFunctionEquiv

/-- The frame weight of each record of a projective partition is of Born
form under the frame state. -/
theorem born_on_records {k : ℕ} (part : ProjectivePartition d k)
    (i : Fin k) :
    F.value (part.proj i) = (bornWeight (frameState F) (part.proj i)).re :=
  born_represents F (part.isEvent i).isEffect

/-- The frame weights across the records of a projective partition sum
to one. -/
theorem records_total_weight {k : ℕ} (part : ProjectivePartition d k) :
    ∑ i, F.value (part.proj i) = 1 :=
  sum_valuation_projective_decomposition (frame_valuation F) Finset.univ
    part.proj (fun i _ => part.isEvent i) part.complete

/-! ## OL-C1: Lueders conditioning on the represented surface

Conditioning acts on the frame state produced by the representation
theorem; every statement carries the nonvanishing-weight guard of the
Lueders module. -/

/-- Conditioning the frame state on an event of nonzero Born weight
yields a state. -/
theorem lueders_state {P : Matrix (Fin d) (Fin d) ℂ} (hP : IsEvent P)
    (hw : bornWeight (frameState F) P ≠ 0) :
    IsState (luedersUpdate (frameState F) P) :=
  luedersUpdate_isState (frameState_isState F) hP hw

/-- **Repeatability.** After conditioning on `P`, the event `P` holds
with Born weight one. -/
theorem lueders_repeatable {P : Matrix (Fin d) (Fin d) ℂ}
    (hP : IsEvent P) (hw : bornWeight (frameState F) P ≠ 0) :
    bornWeight (luedersUpdate (frameState F) P) P = 1 :=
  bornWeight_luedersUpdate_self hP hw

/-- **Idempotence.** Conditioning the frame state twice on the same
event equals conditioning once. -/
theorem lueders_idempotent {P : Matrix (Fin d) (Fin d) ℂ}
    (hP : IsEvent P) (hw : bornWeight (frameState F) P ≠ 0) :
    luedersUpdate (luedersUpdate (frameState F) P) P =
      luedersUpdate (frameState F) P :=
  luedersUpdate_idem hP hw

/-- **Compatibility.** For commuting events, sequential conditioning of
the frame state composes to conditioning on the product event. -/
theorem lueders_compose_of_commute {P Q : Matrix (Fin d) (Fin d) ℂ}
    (hP : IsEvent P) (hc : P * Q = Q * P)
    (hw : bornWeight (frameState F) P ≠ 0) :
    luedersUpdate (luedersUpdate (frameState F) P) Q =
      luedersUpdate (frameState F) (P * Q) :=
  luedersUpdate_luedersUpdate_of_commute hP hc hw

/-- **Order exchange.** Commuting events may be conditioned on in either
order. -/
theorem lueders_order_exchange {P Q : Matrix (Fin d) (Fin d) ℂ}
    (hP : IsEvent P) (hQ : IsEvent Q) (hc : P * Q = Q * P)
    (hwP : bornWeight (frameState F) P ≠ 0)
    (hwQ : bornWeight (frameState F) Q ≠ 0) :
    luedersUpdate (luedersUpdate (frameState F) P) Q =
      luedersUpdate (luedersUpdate (frameState F) Q) P :=
  luedersUpdate_comm hP hQ hc hwP hwQ

/-! ## OL-C3: entanglement to the Tsirelson bound

Algebra-only; consumes four events of the represented algebra and no
frame valuation, so it is re-exported side by side rather than through
`RepresentedPublicFrame`. The bound is an inequality; no attainment
claim is made. -/

section Tsirelson

open scoped Matrix.Norms.L2Operator

/-- **The Tsirelson bound on the represented algebra.** Four projection
events with the cross-party commutation pattern give dichotomic
observables whose CHSH combination has operator norm at most
`2 * sqrt 2`. -/
theorem tsirelson_of_events [NeZero d]
    (A₀ A₁ B₀ B₁ : Matrix (Fin d) (Fin d) ℂ)
    (hA₀ : IsEvent A₀) (hA₁ : IsEvent A₁)
    (hB₀ : IsEvent B₀) (hB₁ : IsEvent B₁)
    (h₀₀ : A₀ * B₀ = B₀ * A₀) (h₀₁ : A₀ * B₁ = B₁ * A₀)
    (h₁₀ : A₁ * B₀ = B₀ * A₁) (h₁₁ : A₁ * B₁ = B₁ * A₁) :
    ‖(1 - 2 • A₀) * (1 - 2 • B₀) +
        (1 - 2 • A₀) * (1 - 2 • B₁) +
        (1 - 2 • A₁) * (1 - 2 • B₀) -
        (1 - 2 • A₁) * (1 - 2 • B₁)‖ ≤ 2 * Real.sqrt 2 :=
  matrix_tsirelson_bound_of_events A₀ A₁ B₀ B₁ hA₀ hA₁ hB₀ hB₁
    h₀₀ h₀₁ h₁₀ h₁₁

end Tsirelson

/-! ## OL-C2: the unique continuous symmetry flow

Consumes a pointwise continuous star-automorphism group of the
represented algebra and no frame valuation; re-exported side by side.
The represented carrier of PR-02 is the composition point. The flow
concerns one full matrix block; no physical clock or source-selected
Hamiltonian is claimed, and the parameter is a real number. -/

section ContinuousFlow

open OPH.Dynamics

/-- **The unique continuous symmetry flow.** A pointwise continuous
real-parameter star-automorphism group of the represented algebra is
conjugation by the coherent unitary family
`stonePropagator Ham t = exp (t * (-i) * Ham)` of one time-independent
self-adjoint generator: the Schroedinger-von Neumann form of the
dynamics on the represented carrier. -/
theorem unique_continuous_flow [NeZero d]
    (α : ℝ → (Matrix (Fin d) (Fin d) ℂ ≃⋆ₐ[ℂ] Matrix (Fin d) (Fin d) ℂ))
    (hzero : ∀ x, α 0 x = x)
    (hgroup : ∀ s t : ℝ, ∀ x, α (s + t) x = α s (α t x))
    (hcont : ∀ x, Continuous fun t : ℝ => α t x) :
    ∃ Ham : Matrix (Fin d) (Fin d) ℂ, IsSelfAdjoint Ham ∧
      ∀ (t : ℝ) (x : Matrix (Fin d) (Fin d) ℂ),
        α t x = stonePropagator Ham t * x * stonePropagator Ham (-t) := by
  haveI : Nonempty (Fin d) := ⟨⟨0, Nat.pos_of_ne_zero (NeZero.ne d)⟩⟩
  exact finiteStoneConverse_of_continuous α hzero hgroup hcont

/-- **Uniqueness of the generator.** Two self-adjoint generators with
the same propagator conjugation differ by one additive real scalar
multiple of the identity. -/
theorem flow_generator_unique [NeZero d]
    {H₁ H₂ : Matrix (Fin d) (Fin d) ℂ}
    (h₁ : IsSelfAdjoint H₁) (h₂ : IsSelfAdjoint H₂)
    (hflow : ∀ (t : ℝ) (x : Matrix (Fin d) (Fin d) ℂ),
      stonePropagator H₁ t * x * stonePropagator H₁ (-t) =
        stonePropagator H₂ t * x * stonePropagator H₂ (-t)) :
    ∃ r : ℝ, H₁ = H₂ + (r : ℂ) • 1 := by
  haveI : Nonempty (Fin d) := ⟨⟨0, Nat.pos_of_ne_zero (NeZero.ne d)⟩⟩
  exact stoneFlow_generator_unique_up_to_real_scalar h₁ h₂ hflow

end ContinuousFlow

/-! ## OL-C5: phase completion and the real/Kraus blindness boundary

The committed dimension-two payload inhabits the PR-04 bundle, and the
blindness boundary is carried alongside as the committed limit. These
statements consume no frame valuation; they are exact operator facts
about the committed payload, re-exported side by side. -/

section PhaseCompletion

open OPH.QFT

/-- The committed algebraic lift inhabits the PR-04 bundle: it is a
projection event, it separates the two Pauli-Y states, and it closes
operator tomography on every fixed-trace slice. The instrument character
of the row (source-produced outcome receipts) is the declared part; this
inhabitant witnesses the algebraic fields only. -/
noncomputable def committedPhaseInstrument : PhaseInstrument where
  lift := sourcePhaseLift
  isEvent := sourcePhaseLift_isEvent
  separates := sourcePhaseLift_distinguishes_Y_states
  completes rho sigma htr hrec hrot hlift :=
    sourcePhaseTomography_injective_on_equalTrace htr (by
      funext i
      fin_cases i
      · simpa [sourcePhaseTomography] using hrec
      · simpa [sourcePhaseTomography] using hrot
      · simpa [sourcePhaseTomography] using hlift)

/-- **Full operator tomography on states.** Under any inhabitant of the
PR-04 bundle, two states agreeing on the record projector, the rotated
projector, and the declared phase operation are equal. -/
theorem phase_states_determined (I : PhaseInstrument)
    {rho sigma : Matrix (Fin 2) (Fin 2) ℂ}
    (hrho : IsState rho) (hsigma : IsState sigma)
    (hrec : bornWeight rho (complexifyRealMatrix recordProjector) =
      bornWeight sigma (complexifyRealMatrix recordProjector))
    (hrot : bornWeight rho (complexifyRealMatrix (conjProjector 3)) =
      bornWeight sigma (complexifyRealMatrix (conjProjector 3)))
    (hlift : bornWeight rho I.lift = bornWeight sigma I.lift) :
    rho = sigma :=
  I.completes rho sigma (hrho.2.trans hsigma.2.symm) hrec hrot hlift

/-- **Real/Kraus blindness.** The generous phase-free closure of the
committed source effects, permitting complements, real linear coarse
graining, real scaling, and pullback by arbitrary real Kraus matrices,
cannot distinguish the two Pauli-Y states: it is tomographically
incomplete. This is the reason PR-04 is a register row rather than a
consequence of the committed payload. -/
theorem real_closure_blind :
    ∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      IsState rho ∧ IsState sigma ∧ rho ≠ sigma ∧
        ∀ E : Matrix (Fin 2) (Fin 2) ℝ,
          RealSourceEffectClosure E →
            bornWeight rho (complexifyRealMatrix E) =
              bornWeight sigma (complexifyRealMatrix E) :=
  realSourceEffectClosure_not_tomographically_complete

/-- **The declared operation is outside the real closure.** No amount of
the phase-free real processing produces the committed lift; the missing
capability is exactly a phase-sensitive operation. -/
theorem phase_lift_outside_real_closure :
    ¬ ∃ E : Matrix (Fin 2) (Fin 2) ℝ,
      RealSourceEffectClosure E ∧
        complexifyRealMatrix E = sourcePhaseLift :=
  sourcePhaseLift_not_in_realSourceEffectClosure

/-- **The committed limit in one statement.** The lift is a projection
in the complex algebra generated by the source-attached pair, it makes
the fixed-trace frame injective, no phase-free real closure produces it,
and the committed run's outcome statistics inhabit only the diagonal
native context. The last conjunct is the exact content of the PR-04
disposition: the algebraic completion is certified while the instrument
receipts are absent from the payload. -/
theorem phase_boundary_summary :
    ComplexSourceAlgebraClosure sourcePhaseLift ∧
      IsEvent sourcePhaseLift ∧
      (∀ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
        Matrix.trace rho = Matrix.trace sigma →
        sourcePhaseTomography rho = sourcePhaseTomography sigma →
        rho = sigma) ∧
      (¬ ∃ E : Matrix (Fin 2) (Fin 2) ℝ,
        RealSourceEffectClosure E ∧
          complexifyRealMatrix E = sourcePhaseLift) ∧
      (∀ c : WebContext,
        (realizedOutcomeCounts c).isSome ↔ c = WebContext.diagonal) :=
  sourcePhaseLift_boundary_summary

end PhaseCompletion

end QuantumSurface

end EventAlgebra

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms EventAlgebra.QuantumSurface.frame_valuation
#print axioms EventAlgebra.QuantumSurface.frameState
#print axioms EventAlgebra.QuantumSurface.frameState_isState
#print axioms EventAlgebra.QuantumSurface.born_represents
#print axioms EventAlgebra.QuantumSurface.busch_gleason
#print axioms EventAlgebra.QuantumSurface.busch_gleason_dimension_two
#print axioms EventAlgebra.QuantumSurface.additive_frame_iff_born
#print axioms EventAlgebra.QuantumSurface.recordSurfaceClassical
#print axioms EventAlgebra.QuantumSurface.born_on_records
#print axioms EventAlgebra.QuantumSurface.records_total_weight
#print axioms EventAlgebra.QuantumSurface.lueders_state
#print axioms EventAlgebra.QuantumSurface.lueders_repeatable
#print axioms EventAlgebra.QuantumSurface.lueders_idempotent
#print axioms EventAlgebra.QuantumSurface.lueders_compose_of_commute
#print axioms EventAlgebra.QuantumSurface.lueders_order_exchange
#print axioms EventAlgebra.QuantumSurface.tsirelson_of_events
#print axioms EventAlgebra.QuantumSurface.unique_continuous_flow
#print axioms EventAlgebra.QuantumSurface.flow_generator_unique
#print axioms EventAlgebra.QuantumSurface.committedPhaseInstrument
#print axioms EventAlgebra.QuantumSurface.phase_states_determined
#print axioms EventAlgebra.QuantumSurface.real_closure_blind
#print axioms EventAlgebra.QuantumSurface.phase_lift_outside_real_closure
#print axioms EventAlgebra.QuantumSurface.phase_boundary_summary
