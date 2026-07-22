import Mathlib

/-!
# Strict Einstein-branch dependency manifests

The manifest validator below is intentionally syntactic.  It proves that a
declared dependency cannot be deleted while retaining strict validation.  It
does **not** claim semantic minimality of the physical receipt list; that would
require isolated full-tower countermodels.

Every dependency predicate inspects numerical, functional, asymptotic, or set
data.  None is the constant proposition `True`.
-/

namespace OPH.EinsteinBranch

open Filter Asymptotics
open scoped Topology

/-- Required finite, analytic, conservation, and physical dependencies of the
corrected composed theorem. -/
inductive Dependency where
  | finiteRepair
  | boundaryFiber
  | typedCommonDomain
  | sourceFactorization
  | lorentzGeometry
  | nineDirectionTomography
  | entropyFirstLaw
  | maxEntEnvelope
  | diamondKernel
  | fixedVolumeArea
  | uniformTail
  | timelikeCoverage
  | wardIdentity
  | bianchiIdentity
  | metricCompatibility
  | connectedComponent
  | universalCoupling
  | vacuumReference
  | physicalScale
  deriving DecidableEq, Fintype, Repr

/-- Raw evaluator data.  These fields are deliberately not proofs. -/
structure BranchCandidate where
  mismatchBefore : ℕ
  mismatchAfter : ℕ
  boundaryFiberCardinality : ℕ
  arrowSourceDigest : ℕ
  arrowTargetDigest : ℕ
  targetLeakCount : ℕ
  eventDimension : ℕ
  tomographyRank : ℕ
  entropyDefect : ℝ
  envelopeDefect : ℝ
  kernelObserved : ℝ
  kernelExpected : ℝ
  areaObserved : ℝ
  areaExpected : ℝ
  ell : ℕ → ℝ
  totalError : ℕ → ℝ
  uncoveredTimelikeDirections : ℕ
  wardDefect : Fin 4 → ℝ
  bianchiDefect : Fin 4 → ℝ
  metricCompatibilityDefect : ℝ
  componentCount : ℕ
  stressSource : ℝ
  geometrySource : ℝ
  vacuumResidue : ℝ
  scaleStabilizer : Set ℝ

/-- Mathematical interpretation of each receipt key. -/
def DependencyHolds (C : BranchCandidate) : Dependency → Prop
  | .finiteRepair => C.mismatchAfter < C.mismatchBefore
  | .boundaryFiber => C.boundaryFiberCardinality = 1
  | .typedCommonDomain => C.arrowSourceDigest = C.arrowTargetDigest
  | .sourceFactorization => C.targetLeakCount = 0
  | .lorentzGeometry => C.eventDimension = 4
  | .nineDirectionTomography => C.tomographyRank = 9
  | .entropyFirstLaw => C.entropyDefect = 0
  | .maxEntEnvelope => C.envelopeDefect = 0
  | .diamondKernel => C.kernelObserved = C.kernelExpected
  | .fixedVolumeArea => C.areaObserved = C.areaExpected
  | .uniformTail =>
      Tendsto C.ell atTop (𝓝 0) ∧
        C.totalError =o[atTop] fun r => C.ell r ^ 4
  | .timelikeCoverage => C.uncoveredTimelikeDirections = 0
  | .wardIdentity => ∀ i, C.wardDefect i = 0
  | .bianchiIdentity => ∀ i, C.bianchiDefect i = 0
  | .metricCompatibility => C.metricCompatibilityDefect = 0
  | .connectedComponent => C.componentCount = 1
  | .universalCoupling => C.stressSource = C.geometrySource
  | .vacuumReference => C.vacuumResidue = 0
  | .physicalScale => C.scaleStabilizer = ({1} : Set ℝ)

/-- A strict manifest contains every required key and evidence for every key.
Because `Dependency` is a closed enum and the manifest equals `univ`, it also
contains no unexpected keys.  `Finset` supplies duplicate-freedom. -/
def StrictManifest (C : BranchCandidate) (manifest : Finset Dependency) : Prop :=
  manifest = Finset.univ ∧ ∀ d ∈ manifest, DependencyHolds C d

/-- Strict validation exposes every individual dependency. -/
theorem strictManifest_dependency {C : BranchCandidate}
    {manifest : Finset Dependency} (h : StrictManifest C manifest)
    (d : Dependency) : DependencyHolds C d := by
  exact h.2 d (by rw [h.1]; simp)

/-- Deleting any required key makes strict validation fail. -/
theorem strictManifest_erase_required_fails {C : BranchCandidate}
    {manifest : Finset Dependency} (_h : StrictManifest C manifest)
    (d : Dependency) : ¬ StrictManifest C (manifest.erase d) := by
  intro herased
  have hmem : d ∈ manifest.erase d := by
    rw [herased.1]
    simp
  exact (Finset.notMem_erase d manifest) hmem

/-- Replacing one passing scalar equality by a nonzero defect fails its
receipt.  This is a generic negative-control lemma used by the manifest
audit. -/
theorem nonzero_defect_fails_entropy (C : BranchCandidate)
    (h : C.entropyDefect ≠ 0) :
    ¬ DependencyHolds C .entropyFirstLaw := by
  simpa [DependencyHolds] using h

/-- A nontrivial rescaling stabilizer fails the independent-scale receipt. -/
theorem extra_scale_stabilizer_fails (C : BranchCandidate) (a : ℝ)
    (ha : a ∈ C.scaleStabilizer) (ha1 : a ≠ 1) :
    ¬ DependencyHolds C .physicalScale := by
  intro h
  have : a ∈ ({1} : Set ℝ) := by
    rw [← h]
    exact ha
  simpa using ha1 (Set.mem_singleton_iff.mp this)

/-! ## Per-theorem axiom audit -/

#print axioms strictManifest_dependency
#print axioms strictManifest_erase_required_fails
#print axioms nonzero_defect_fails_entropy
#print axioms extra_scale_stabilizer_fails

end OPH.EinsteinBranch
