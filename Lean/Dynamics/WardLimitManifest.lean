import Dynamics.ProtectedCharge
import RegionalContinuity

namespace OPH.Dynamics

/-!
# Typed manifest for promoting finite continuity to a Ward identity

The B5 incidence and protected-charge theorems stop at a finite regulator.
This structure bundles the additional receipts required before a gravity
argument may consume finite continuity as a continuum Ward identity.  It is a
guarded premise manifest, not a source construction.  Source derivation,
refinement transport, residual semantics, chart realization, and tower
carriage remain downstream evidence: inhabiting the structure is not by
itself a physical identification.
-/

/-- Exact interface between a family of finite conservation objects and one
candidate continuum stress tensor.  Later realizations must define the
source, transport, scale, test, residual, chart, and tower semantics.  Once
the finite residual is proved zero and converges to the continuum residual,
the Ward equality is derived rather than stored as an assumption. -/
structure WardLimitManifest
    (Regulator Chart StressTensor TestField Tower : Type*) where
  /-- Nonvacuity: the regulator family has a named member. -/
  regulatorWitness : Regulator
  /-- The actual finite B5 continuity object at each regulator. -/
  finiteContinuity : Regulator → OPH.RegionalContinuity.FiniteContinuityWitness
  /-- At least one regulator has a datum that is not completely all-zero.
  This does not require a nonzero current or a nonzero candidate limit. -/
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
  /-- A concrete positive regulator scale shrinks along the cofinal sequence.
  Distinct regulator labels alone are therefore insufficient. -/
  scale : Regulator → ℝ
  scalePositive : ∀ r, 0 < scale r
  scaleTendsToZero : Filter.Tendsto
    (fun n => scale (cofinalSequence n)) Filter.atTop (nhds 0)
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
  real-valued pairing against an admissible separating test-field class. -/
  AdmissibleTest : TestField → Prop
  testFieldWitness : TestField
  testFieldWitnessAdmissible : AdmissibleTest testFieldWitness
  limitStress : StressTensor
  stressPairing : StressTensor → TestField → ℝ
  testFieldsSeparate : ∀ {s t},
    (∀ test, AdmissibleTest test →
      stressPairing s test = stressPairing t test) → s = t
  stressResponseNonzero : ∃ r test,
    AdmissibleTest test ∧ stressPairing (realizedStress r) test ≠ 0
  weakConverges : ∀ test, AdmissibleTest test →
    Filter.Tendsto
      (fun n => stressPairing (realizedStress (cofinalSequence n)) test)
      Filter.atTop
      (nhds (stressPairing limitStress test))
  /-- F1 must define this as the smeared finite divergence produced from the
  actual continuity witness.  Its vanishing proof is the finite input to the
  limiting Ward theorem. -/
  finiteWardResidual :
    OPH.RegionalContinuity.FiniteContinuityWitness → TestField → ℝ
  finiteResidualVanishes : ∀ r test, AdmissibleTest test →
    finiteWardResidual (finiteContinuity r) test = 0
  /-- The continuum residual is not allowed to be the globally zero function.
  This guard is necessary but does not by itself identify the residual with a
  physical distributional divergence. -/
  wardResidual : StressTensor → TestField → ℝ
  wardResidualNontrivial : ∃ stress test,
    AdmissibleTest test ∧ wardResidual stress test ≠ 0
  /-- The finite residuals converge to the continuum residual on every
  admissible test.  Together with `finiteResidualVanishes`, this derives the
  Ward equality by uniqueness of real limits. -/
  residualConverges : ∀ test, AdmissibleTest test →
    Filter.Tendsto
      (fun n => finiteWardResidual (finiteContinuity (cofinalSequence n)) test)
      Filter.atTop
      (nhds (wardResidual limitStress test))
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

/-- The continuum Ward equality follows from exact finite residual vanishing
and convergence to the declared continuum residual.  Its physical reading
depends on the concrete F1 definitions used to inhabit the manifest. -/
theorem WardLimitManifest.wardIdentity
    {Regulator Chart StressTensor TestField Tower : Type*}
    (manifest : WardLimitManifest Regulator Chart StressTensor TestField Tower)
    (test : TestField) (htest : manifest.AdmissibleTest test) :
    manifest.wardResidual manifest.limitStress test = 0 := by
  have hzero : Filter.Tendsto
      (fun _ : ℕ => (0 : ℝ)) Filter.atTop (nhds 0) :=
    tendsto_const_nhds
  have hfun :
      (fun n => manifest.finiteWardResidual
        (manifest.finiteContinuity (manifest.cofinalSequence n)) test) =
      (fun _ : ℕ => (0 : ℝ)) := by
    funext n
    exact manifest.finiteResidualVanishes
      (manifest.cofinalSequence n) test htest
  have hresidualZero : Filter.Tendsto
      (fun n => manifest.finiteWardResidual
        (manifest.finiteContinuity (manifest.cofinalSequence n)) test)
      Filter.atTop (nhds 0) := by
    rw [hfun]
    exact hzero
  exact tendsto_nhds_unique
    (manifest.residualConverges test htest) hresidualZero

#print axioms WardLimitManifest.wardIdentity

end OPH.Dynamics
