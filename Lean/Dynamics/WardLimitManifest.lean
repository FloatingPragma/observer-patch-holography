import Dynamics.ProtectedCharge
import RegionalContinuity

namespace OPH.Dynamics

/-!
# Typed manifest for promoting finite continuity to a Ward identity

The B5 incidence and protected-charge theorems stop at a finite regulator.
This structure bundles the additional receipts required before a gravity
argument may consume them as a continuum Ward identity.  Its mathematical
core rules out the constant-regulator, zero-field, and identically-true Ward
shortcuts.  Source derivation, refinement transport, and tower carriage remain
named downstream evidence relations: inhabiting the structure is not by
itself a physical identification.
-/

/-- Exact interface between a family of finite conservation objects and one
candidate continuum stress tensor.  Later realizations choose the meanings of
source derivation, refinement transport, chart change, and tower carriage, but
weak convergence and the Ward identity have fixed real-valued semantics. -/
structure WardLimitManifest
    (Regulator Chart StressTensor TestField Tower : Type*) where
  /-- Nonvacuity: the regulator family has a named member. -/
  regulatorWitness : Regulator
  /-- The actual finite B5 continuity object at each regulator. -/
  finiteContinuity : Regulator → OPH.RegionalContinuity.FiniteContinuityWitness
  /-- At least one finite regulator carries nonzero load, source, or current.
  Thus the zero continuity witness cannot discharge the handoff. -/
  finiteWitnessNonzero : ∃ r,
    (finiteContinuity r).q ≠ 0 ∨
    (finiteContinuity r).qNext ≠ 0 ∨
    (finiteContinuity r).source ≠ 0 ∨
    (finiteContinuity r).current ≠ 0
  /-- Receipt 1: the complete finite witness comes from the declared source. -/
  SourceDerived :
    Regulator → OPH.RegionalContinuity.FiniteContinuityWitness → Prop
  sourceDerived : ∀ r, SourceDerived r (finiteContinuity r)
  /-- Receipt 2: charge and current transport together along refinements. -/
  Refines : Regulator → Regulator → Prop
  refinesRefl : ∀ r, Refines r r
  refinesTrans : ∀ {r s t}, Refines r s → Refines s t → Refines r t
  refinesDirected : ∀ r s, ∃ t, Refines r t ∧ Refines s t
  cofinalSequence : ℕ → Regulator
  sequenceInjective : Function.Injective cofinalSequence
  sequenceMonotone : ∀ n, Refines (cofinalSequence n) (cofinalSequence (n + 1))
  sequenceCofinal : ∀ r, ∃ n, Refines r (cofinalSequence n)
  TransportCompatible : Regulator → Regulator →
    OPH.RegionalContinuity.FiniteContinuityWitness →
    OPH.RegionalContinuity.FiniteContinuityWitness → Prop
  refinementNatural : ∀ {r s}, Refines r s →
    TransportCompatible r s (finiteContinuity r) (finiteContinuity s)
  /-- Receipt 3: one chart realization supplies the temporal and spatial
  components of a single stress tensor. -/
  chart : Regulator → Chart
  realize : Chart → OPH.RegionalContinuity.FiniteContinuityWitness → StressTensor
  realizedStress : Regulator → StressTensor
  stressComponentsRealized : ∀ r,
    realizedStress r = realize (chart r) (finiteContinuity r)
  /-- Receipt 4: the realized tensors converge weakly through an explicit
  real-valued pairing against a nonempty declared test-field class. -/
  testFieldWitness : TestField
  limitStress : StressTensor
  stressPairing : StressTensor → TestField → ℝ
  stressResponseNonzero : ∃ r test,
    stressPairing (realizedStress r) test ≠ 0
  weakConverges : ∀ test,
    Filter.Tendsto
      (fun n => stressPairing (realizedStress (cofinalSequence n)) test)
      Filter.atTop
      (nhds (stressPairing limitStress test))
  /-- The Ward residual has fixed equality-to-zero semantics, and is required
  to distinguish at least one stress/test pair so that the identity is not an
  identically true predicate. -/
  wardResidual : StressTensor → TestField → ℝ
  wardResidualNontrivial : ∃ stress test, wardResidual stress test ≠ 0
  wardIdentity : ∀ test, wardResidual limitStress test = 0
  /-- Receipt 5: the finite stress realization commutes with chart changes. -/
  Recharts : Chart → Chart → Prop
  RechartStress : Chart → Chart → StressTensor → StressTensor
  rechartsRefl : ∀ chart, Recharts chart chart
  rechartsTrans : ∀ {a b c}, Recharts a b → Recharts b c → Recharts a c
  rechartIdentity : ∀ chart stress,
    RechartStress chart chart stress = stress
  rechartCompose : ∀ {a b c} (_hab : Recharts a b) (_hbc : Recharts b c)
      (stress : StressTensor),
    RechartStress b c (RechartStress a b stress) =
      RechartStress a c stress
  chartCovariant : ∀ {r s}, Recharts (chart r) (chart s) →
    RechartStress (chart r) (chart s) (realizedStress r) = realizedStress s
  /-- Receipt 6: one inhabited tower carries the same source objects and
  their realized stress tensor through every regulator. -/
  commonTower : Tower
  Carries : Tower → Regulator →
    OPH.RegionalContinuity.FiniteContinuityWitness → StressTensor → Prop
  commonTowerCarries : ∀ r,
    Carries commonTower r (finiteContinuity r) (realizedStress r)

end OPH.Dynamics
