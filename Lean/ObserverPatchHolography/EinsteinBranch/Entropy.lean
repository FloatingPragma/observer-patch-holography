import Mathlib

/-!
# Exact finite entropy and MaxEnt algebra for the Einstein branch

The central-sector convention here is the corrected one used by the compact
paper:

* `S_bulk = H(p) + sum p_alpha S(rho_bulk,alpha)`;
* `S_edge = sum p_alpha log d_alpha`;
* `Z` acts sectorwise by `log d_alpha`; the `-log p_alpha` term belongs to
  the bulk modular generator.

Theorems in this file are finite algebra.  No continuum limit, area law, KMS
normalisation, physical stress identification, or MaxEnt existence claim is
made.  Those enter the composed theorem as explicit premises.
-/

namespace OPH.EinsteinBranch

open scoped BigOperators

universe u

noncomputable section

/-! ## Direct-sum entropy bookkeeping -/

/-- Scalar data of a finite central direct sum.  `blockEntropy` is the von
Neumann entropy of the normalised bulk block; `edgeDim` is the dimension of
the maximally mixed edge factor. -/
structure CentralSectorData (α : Type u) where
  probability : α → ℝ
  blockEntropy : α → ℝ
  edgeDim : α → ℕ

/-- Shannon contribution carried by the bulk convention. -/
def shannonTerm {α : Type u} [Fintype α] (D : CentralSectorData α) : ℝ :=
  -∑ a, D.probability a * Real.log (D.probability a)

/-- Bulk entropy, including the central Shannon term. -/
def bulkEntropy {α : Type u} [Fintype α] (D : CentralSectorData α) : ℝ :=
  shannonTerm D + ∑ a, D.probability a * D.blockEntropy a

/-- Edge entropy only; it does not contain the Shannon term. -/
def edgeEntropy {α : Type u} [Fintype α] (D : CentralSectorData α) : ℝ :=
  ∑ a, D.probability a * Real.log (D.edgeDim a : ℝ)

/-- Entropy formula obtained by applying the direct-sum and tensor-product
rules to `oplus_a p_a (rho_bulk,a tensor I_edge,a / d_a)`.  It is written
independently of `bulkEntropy + edgeEntropy` so the split below is a theorem,
not definitional truth. -/
def directSumEntropy {α : Type u} [Fintype α] (D : CentralSectorData α) : ℝ :=
  -∑ a, D.probability a * Real.log (D.probability a) +
    ∑ a, D.probability a *
      (D.blockEntropy a + Real.log (D.edgeDim a : ℝ))

/-- Exact direct-sum bulk/edge split in the corrected convention. -/
theorem directSumEntropy_eq_bulk_add_edge {α : Type u} [Fintype α]
    (D : CentralSectorData α) :
    directSumEntropy D = bulkEntropy D + edgeEntropy D := by
  simp only [directSumEntropy, bulkEntropy, edgeEntropy, shannonTerm]
  simp_rw [mul_add]
  rw [Finset.sum_add_distrib]
  ring

/-! ## Edge/central variation -/

/-- First-order variation of the central probabilities with a declared
central modular weight. -/
structure CentralWeightVariation (α : Type u) where
  deltaProbability : α → ℝ
  edgeDim : α → ℕ
  centralWeight : α → ℝ

def centralVariation {α : Type u} [Fintype α]
    (V : CentralWeightVariation α) : ℝ :=
  ∑ a, V.centralWeight a * V.deltaProbability a

def edgeEntropyVariation {α : Type u} [Fintype α]
    (V : CentralWeightVariation α) : ℝ :=
  ∑ a, Real.log (V.edgeDim a : ℝ) * V.deltaProbability a

/-- The edge-normalisation receipt has mathematical content: every central
weight must equal the logarithm of the corresponding edge dimension. -/
def EdgeNormalization {α : Type u} (V : CentralWeightVariation α) : Prop :=
  ∀ a, V.centralWeight a = Real.log (V.edgeDim a : ℝ)

theorem centralVariation_eq_edgeEntropyVariation {α : Type u} [Fintype α]
    (V : CentralWeightVariation α) (h : EdgeNormalization V) :
    centralVariation V = edgeEntropyVariation V := by
  apply Finset.sum_congr rfl
  intro a ha
  rw [h a]

/-- Exact defect when an arbitrary central normalisation is used. -/
theorem central_edge_normalization_defect {α : Type u} [Fintype α]
    (V : CentralWeightVariation α) :
    centralVariation V - edgeEntropyVariation V =
      ∑ a, (V.centralWeight a - Real.log (V.edgeDim a : ℝ)) *
        V.deltaProbability a := by
  simp only [centralVariation, edgeEntropyVariation]
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro a ha
  ring

/-! ## Finite first law as linear-functional composition -/

/-- Raw finite/type-I first-law data on a variation space.  The physical
content is kept in the failable predicate `FiniteFirstLawPremises`: the
entropy differential must equal the modular pairing, the modular generator
must split as `2*pi B + Z`, and `Z` must equal the edge differential. -/
structure FiniteFirstLawData (V : Type u) [AddCommGroup V] [Module ℝ V] where
  entropyDifferential : V →ₗ[ℝ] ℝ
  modularPairing : V →ₗ[ℝ] ℝ
  bulkPairing : V →ₗ[ℝ] ℝ
  centralPairing : V →ₗ[ℝ] ℝ
  edgeDifferential : V →ₗ[ℝ] ℝ

def FiniteFirstLawPremises {V : Type u} [AddCommGroup V] [Module ℝ V]
    (D : FiniteFirstLawData V) : Prop :=
  D.entropyDifferential = D.modularPairing ∧
  D.modularPairing = (2 * Real.pi) • D.bulkPairing + D.centralPairing ∧
  D.centralPairing = D.edgeDifferential

/-- Exact finite bulk/edge/central first law at fixed operators. -/
theorem finite_bulk_edge_central_first_law {V : Type u}
    [AddCommGroup V] [Module ℝ V] (D : FiniteFirstLawData V)
    (h : FiniteFirstLawPremises D) (δρ : V) :
    D.entropyDifferential δρ =
      2 * Real.pi * D.bulkPairing δρ + D.edgeDifferential δρ := by
  rcases h with ⟨hEntropy, hSplit, hEdge⟩
  rw [hEntropy, hSplit, hEdge]
  simp

/-- A concrete passing first-law datum, showing the premise predicate is
satisfiable. -/
def passingFirstLawData : FiniteFirstLawData ℝ where
  entropyDifferential := (3 * Real.pi + 2) • LinearMap.id
  modularPairing := (3 * Real.pi + 2) • LinearMap.id
  bulkPairing := (3 / 2 : ℝ) • LinearMap.id
  centralPairing := (2 : ℝ) • LinearMap.id
  edgeDifferential := (2 : ℝ) • LinearMap.id

theorem passingFirstLawData_passes :
    FiniteFirstLawPremises passingFirstLawData := by
  constructor
  · rfl
  constructor
  · apply LinearMap.ext
    intro x
    simp [passingFirstLawData]
    ring
  · rfl

/-- A concrete failing datum, formally ruling out a constant-`True`
first-law receipt. -/
def failingFirstLawData : FiniteFirstLawData ℝ where
  entropyDifferential := 0
  modularPairing := 0
  bulkPairing := 0
  centralPairing := 0
  edgeDifferential := LinearMap.id

theorem failingFirstLawData_fails :
    ¬ FiniteFirstLawPremises failingFirstLawData := by
  rintro ⟨hEntropy, hSplit, hEdge⟩
  have h := LinearMap.congr_fun hEdge 1
  norm_num [failingFirstLawData] at h

/-! ## MaxEnt envelope identity -/

/-- First-order scalar data along a declared MaxEnt family.  The Gibbs
differential equation below is the finite-dimensional input; normalisation
and unit constraint speed make the envelope identity an exact consequence. -/
structure MaxEntTangent where
  entropyRate : ℝ
  constraintRate : ℝ
  normalizationRate : ℝ
  multiplier : ℝ
  normalizationMultiplier : ℝ

def MaxEntEnvelopePremises (M : MaxEntTangent) : Prop :=
  M.normalizationRate = 0 ∧
  M.constraintRate = 1 ∧
  M.entropyRate =
    M.multiplier * M.constraintRate +
      M.normalizationMultiplier * M.normalizationRate

/-- Exact finite MaxEnt envelope identity `dS/dt = lambda`. -/
theorem maxEnt_envelope_identity (M : MaxEntTangent)
    (h : MaxEntEnvelopePremises M) :
    M.entropyRate = M.multiplier := by
  rcases h with ⟨hnorm, hconstraint, hgibbs⟩
  rw [hnorm, hconstraint] at hgibbs
  simpa using hgibbs

def passingMaxEntTangent : MaxEntTangent where
  entropyRate := 2 * Real.pi
  constraintRate := 1
  normalizationRate := 0
  multiplier := 2 * Real.pi
  normalizationMultiplier := 7

theorem passingMaxEntTangent_passes :
    MaxEntEnvelopePremises passingMaxEntTangent := by
  simp [MaxEntEnvelopePremises, passingMaxEntTangent]

def failingMaxEntTangent : MaxEntTangent where
  entropyRate := 0
  constraintRate := 1
  normalizationRate := 0
  multiplier := 1
  normalizationMultiplier := 0

theorem failingMaxEntTangent_fails :
    ¬ MaxEntEnvelopePremises failingMaxEntTangent := by
  simp [MaxEntEnvelopePremises, failingMaxEntTangent]

/-! ## Per-theorem axiom audit -/

#print axioms directSumEntropy_eq_bulk_add_edge
#print axioms centralVariation_eq_edgeEntropyVariation
#print axioms central_edge_normalization_defect
#print axioms finite_bulk_edge_central_first_law
#print axioms passingFirstLawData_passes
#print axioms failingFirstLawData_fails
#print axioms maxEnt_envelope_identity
#print axioms passingMaxEntTangent_passes
#print axioms failingMaxEntTangent_fails

end

end OPH.EinsteinBranch
