import ObserverPatchHolography.EinsteinBranch.Consensus
import ObserverPatchHolography.EinsteinBranch.Entropy
import ObserverPatchHolography.EinsteinBranch.SmallBall
import ObserverPatchHolography.EinsteinBranch.Tensor
import ObserverPatchHolography.EinsteinBranch.Manifest

/-!
# Corrected typed Einstein-branch composition

This file states the physical and continuum inputs without claiming their
truth.  In particular, `ContinuumEinsteinPremises` contains the null-balance,
Ward, contracted-Bianchi, connected-chart, vacuum-reference, universal-
coupling, and independent-scale premises.  They are theorem arguments, not
global Lean axioms.  No instance of this structure is constructed here.

`composedEinsteinBranch` returns a certificate containing the finite algebra,
boundary/common-domain composition, exact small-ball arithmetic, strict
dependency audit, and the componentwise Einstein equation.  Thus the theorem
proves the implication on the admissible branch; it does not prove that an
Einstein-admissible source-derived tower exists.
-/

namespace OPH.EinsteinBranch

universe u v

noncomputable section

/-- Explicit continuum/physical premises on one discrete chart model of a
connected component.  `ddiv` is the exact finite-difference composition used
to machine-check Ward+Bianchi constancy; the truth of the continuum Ward and
Bianchi identities remains an input represented by `ward` and `bianchi`. -/
structure ContinuumEinsteinPremises (P : Type u) where
  step : Fin 4 → P → P
  base : P
  geometry : P → Mat 3
  stress : P → Mat 3
  entropyStress : P → Mat 3
  coupling : ℝ
  newton : ℝ
  referenceLambda : ℝ
  geometry_symmetric : ∀ p i j, geometry p i j = geometry p j i
  stress_symmetric : ∀ p i j, stress p i j = stress p j i
  /-- Scaling-limit/null-horizon balance, consumed as a premise. -/
  nullBalance : ∀ p k, quadOf (eta 3) k = 0 →
    quadOf (geometry p) k = coupling * quadOf (stress p) k
  /-- Weak Ward identity for the reconstructed stress. -/
  ward : ∀ p j, ddiv step stress p j = 0
  /-- Contracted Bianchi identity for the smooth geometric side. -/
  bianchi : ∀ p j, ddiv step geometry p j = 0
  /-- Connectivity of this component, in the symmetric step closure. -/
  connected : ∀ q, SymmReachable step base q
  /-- Universal coupling: entropy and null tomography use the same source. -/
  universalCoupling : entropyStress = stress
  /-- Vacuum-reference/base condition fixing the pointwise metric residue. -/
  vacuumReference :
    geometry base 0 0 =
      coupling * stress base 0 0 + referenceLambda * eta 3 0 0
  /-- Independent scale identification of the coupling. -/
  physicalScale : coupling = 8 * Real.pi * newton

/-- Ward+Bianchi composition promotes the pointwise metric ambiguity to the
vacuum-referenced constant on the connected component. -/
theorem continuumEinstein_from_explicit_premises {P : Type u}
    (A : ContinuumEinsteinPremises P) :
    ∀ p i j,
      A.geometry p i j =
        8 * Real.pi * A.newton * A.stress p i j +
          A.referenceLambda * eta 3 i j := by
  obtain ⟨Λ, hΛ⟩ := einstein_equation_with_constant_symm
    A.step A.base A.connected A.geometry A.stress A.coupling
    A.geometry_symmetric A.stress_symmetric A.nullBalance A.bianchi A.ward
  have hbase := hΛ A.base (0 : Fin 4) (0 : Fin 4)
  have hη : eta 3 (0 : Fin 4) (0 : Fin 4) = -1 := eta_zero_zero
  have hconst : Λ = A.referenceLambda := by
    have href := A.vacuumReference
    rw [hη] at hbase href
    linarith
  intro p i j
  rw [hΛ p i j, hconst, A.physicalScale]

/-- Outputs of the full algebraic/compositional audit. -/
structure EinsteinCompositionConclusion
    {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T)
    (r : ι) (q q' : T.Quot r)
    {V : Type v} [AddCommGroup V] [Module ℝ V]
    (D : FiniteFirstLawData V) (δρ : V)
    (M : MaxEntTangent) (S : SmallBallPremises)
    (tail : ScalingTailData)
    {P : Type u} (A : ContinuumEinsteinPremises P)
    (candidate : BranchCandidate) : Prop where
  dependencies : ∀ d, DependencyHolds candidate d
  commonReadouts :
    ReadoutAgreement R r (T.normalForm r q) (T.normalForm r q')
  typedArrows : TypedArrowCommutation R
  refinementNaturality : ReadoutNaturality R
  finiteFirstLaw :
    D.entropyDifferential δρ =
      2 * Real.pi * D.bulkPairing δρ + D.edgeDifferential δρ
  maxEntEnvelope : M.entropyRate = M.multiplier
  smallBallArithmetic :
    S.geometryContraction = 8 * Real.pi * S.G * S.stressContraction
  scalingTail : ScalingTailReceipt tail
  universalSource : A.entropyStress = A.stress
  einstein : ∀ p i j,
    A.geometry p i j =
      8 * Real.pi * A.newton * A.stress p i j +
        A.referenceLambda * eta 3 i j

/-- The corrected composed branch-entry theorem.

Every finite, asymptotic, and physical premise is visible in the signature.
`hManifest` is strict and deletion-sensitive.  The theorem constructs no
tower and establishes no premise by fiat. -/
theorem composedEinsteinBranch
    {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T)
    (candidate : BranchCandidate) (manifest : Finset Dependency)
    (hManifest : StrictManifest candidate manifest)
    (hFiber : BoundaryFiber T)
    (hArrows : TypedArrowCommutation R)
    (hNatural : ReadoutNaturality R)
    (r : ι) (q q' : T.Quot r)
    (hq : T.consistent r q) (hq' : T.consistent r q')
    (hBoundary : T.boundary r q = T.boundary r q')
    {V : Type v} [AddCommGroup V] [Module ℝ V]
    (D : FiniteFirstLawData V) (hFirstLaw : FiniteFirstLawPremises D)
    (δρ : V)
    (M : MaxEntTangent) (hEnvelope : MaxEntEnvelopePremises M)
    (S : SmallBallPremises)
    (tail : ScalingTailData) (hTail : ScalingTailReceipt tail)
    {P : Type u} (A : ContinuumEinsteinPremises P) :
    EinsteinCompositionConclusion R r q q' D δρ M S tail A candidate := by
  have hdeps : ∀ d, DependencyHolds candidate d :=
    fun d => strictManifest_dependency hManifest d
  have hread := boundary_fiber_readout_composition R hFiber r q q' hq hq' hBoundary
  have hfirst := finite_bulk_edge_central_first_law D hFirstLaw δρ
  have hmaxent := maxEnt_envelope_identity M hEnvelope
  have hsmall := restFrameEinsteinRelation S
  have heinstein := continuumEinstein_from_explicit_premises A
  exact ⟨hdeps, hread, hArrows, hNatural, hfirst, hmaxent, hsmall, hTail,
    A.universalCoupling, heinstein⟩

/-- Bundled form of `composedEinsteinBranch`.  This is the entry point for an
explicit `EinsteinAdmissibleTower`; it still takes every continuum,
asymptotic, entropy, MaxEnt, manifest, vacuum-reference, coupling, and scale
premise as data rather than constructing any of them. -/
theorem composedEinsteinAdmissibleTower
    {ι : Type u} [Preorder ι] (E : EinsteinAdmissibleTower ι)
    (candidate : BranchCandidate) (manifest : Finset Dependency)
    (hManifest : StrictManifest candidate manifest)
    (r : ι) (q q' : E.consensus.Quot r)
    (hq : E.consensus.consistent r q)
    (hq' : E.consensus.consistent r q')
    (hBoundary : E.consensus.boundary r q = E.consensus.boundary r q')
    {V : Type v} [AddCommGroup V] [Module ℝ V]
    (D : FiniteFirstLawData V) (hFirstLaw : FiniteFirstLawPremises D)
    (δρ : V)
    (M : MaxEntTangent) (hEnvelope : MaxEntEnvelopePremises M)
    (S : SmallBallPremises)
    (tail : ScalingTailData) (hTail : ScalingTailReceipt tail)
    {P : Type u} (A : ContinuumEinsteinPremises P) :
    EinsteinCompositionConclusion E.readouts r q q' D δρ M S tail A candidate := by
  exact composedEinsteinBranch E.readouts candidate manifest hManifest
    E.boundaryFiber E.typedArrows E.refinementNaturality r q q' hq hq'
    hBoundary D hFirstLaw δρ M hEnvelope S tail hTail A

/-! ## Per-theorem axiom audit -/

#print axioms continuumEinstein_from_explicit_premises
#print axioms composedEinsteinBranch
#print axioms composedEinsteinAdmissibleTower

end


end OPH.EinsteinBranch
