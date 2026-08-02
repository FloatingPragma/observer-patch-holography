import SeamCurrentDirichletGenerator

open scoped BigOperators

namespace OPH.SeamCurrentFreePhotonLift

open OPH.PrimitivePortTranslationBridge

/-- The unique Cartesian carrier used by the transverse-oscillator packet.  The imported
modules expose definitionally equal `Vec3` abbreviations, so this file fixes
one qualified spelling and never relies on namespace-open resolution. -/
abbrev FreePhotonVec3 := OPH.PrimitivePortTranslationBridge.Vec3

/-!
# Basis-free transverse oscillator lift of the seam-current symbol

This file constructs a mathematical transverse oscillator, or free-vector,
candidate carried by the FZ-12 seam-current symbol.  At every nonzero momentum
the declared transverse space is the kernel of the momentum covector.  Its
dimension is two.  The construction uses the global transverse projector
rather than a global choice of two polarization vectors, which cannot be made
continuously over the momentum sphere.

The scalar seam-current symbol acts identically on that rank-two space.  The
declared Hamiltonian mode generator then has second component
`A'' = -Lambda A`.  The declared algebraic first variation of its quadratic
energy vanishes, and `omega^2 = Lambda` on the nonnegative square-root branch.

These are conditional oscillator mathematics.  The file does not identify
the carrier completion with physical position, the Hamiltonian parameter with
a physical clock, or the transverse fiber with any laboratory field.  It does
not establish a gauge quotient, Maxwell dynamics, Lorentz covariance,
locality, causality, reality pairing, quantization, trajectory existence, a
carrier rest frame, a physical length, an electron or positron action, an
interaction, a boost law, a source, or a readout.
-/

/-- The momentum covector associated with the declared Cartesian response
metric. -/
noncomputable def momentumCovector (k : FreePhotonVec3) :
    FreePhotonVec3 →ₗ[ℝ] ℝ where
  toFun v := dot k v
  map_add' x y := dot_add_right k x y
  map_smul' a x := by
    change dot k (a • x) = a * dot k x
    exact dot_smul_right a k x

/-- The basis-free transverse candidate space at momentum `k`. -/
noncomputable def TransversePolarization (k : FreePhotonVec3) :
    Submodule ℝ FreePhotonVec3 :=
  LinearMap.ker (momentumCovector k)

theorem dot_comm (x y : FreePhotonVec3) : dot x y = dot y x := by
  unfold dot
  apply Finset.sum_congr rfl
  intro d _
  ring

theorem dot_sub_right (x y z : FreePhotonVec3) :
    dot x (y - z) = dot x y - dot x z := by
  simp only [sub_eq_add_neg, dot_add_right, dot_neg_right]

theorem dot_self_pos {k : FreePhotonVec3} (hk : k ≠ 0) : 0 < dot k k := by
  unfold dot
  have hexists : ∃ d : Fin 3, k d ≠ 0 := by
    by_contra hall
    apply hk
    funext d
    exact not_ne_iff.mp (fun hd ↦ hall ⟨d, hd⟩)
  obtain ⟨d, hd⟩ := hexists
  exact Finset.sum_pos' (fun e _ ↦ mul_self_nonneg (k e))
    ⟨d, Finset.mem_univ d, mul_self_pos.mpr hd⟩

/-- A nonzero momentum covector is onto the scalar field. -/
theorem momentumCovector_surjective {k : FreePhotonVec3} (hk : k ≠ 0) :
    Function.Surjective (momentumCovector k) := by
  intro r
  have hkk : dot k k ≠ 0 := ne_of_gt (dot_self_pos hk)
  refine ⟨(r / dot k k) • k, ?_⟩
  change dot k ((r / dot k k) • k) = r
  rw [dot_smul_right]
  field_simp [hkk]

/-- The basis-free transverse candidate space has exactly two real
dimensions at every nonzero momentum. -/
theorem transversePolarization_finrank {k : FreePhotonVec3} (hk : k ≠ 0) :
    Module.finrank ℝ (TransversePolarization k) = 2 := by
  have hrange : LinearMap.range (momentumCovector k) = ⊤ :=
    LinearMap.range_eq_top.mpr (momentumCovector_surjective hk)
  have hdim :=
    LinearMap.finrank_range_add_finrank_ker (momentumCovector k)
  rw [hrange] at hdim
  change Module.finrank ℝ (LinearMap.ker (momentumCovector k)) = 2
  have hdim' :
      1 + Module.finrank ℝ (LinearMap.ker (momentumCovector k)) = 3 := by
    simpa [hrange, FreePhotonVec3,
      OPH.PrimitivePortTranslationBridge.Vec3] using hdim
  omega

/-- Orthogonal projection onto the transverse candidate space.  Its
formula is basis-free and is used only for nonzero `k`. -/
noncomputable def transverseProjector
    (k v : FreePhotonVec3) : FreePhotonVec3 :=
  v - (dot k v / dot k k) • k

/-- The transverse projector lands in the momentum-orthogonal plane. -/
theorem transverseProjector_dot_zero {k : FreePhotonVec3} (hk : k ≠ 0)
    (v : FreePhotonVec3) :
    dot k (transverseProjector k v) = 0 := by
  have hkk : dot k k ≠ 0 := ne_of_gt (dot_self_pos hk)
  unfold transverseProjector
  rw [dot_sub_right, dot_smul_right]
  field_simp [hkk]
  ring

theorem transverseProjector_mem {k : FreePhotonVec3} (hk : k ≠ 0)
    (v : FreePhotonVec3) :
    transverseProjector k v ∈ TransversePolarization k := by
  change dot k (transverseProjector k v) = 0
  exact transverseProjector_dot_zero hk v

/-- A transverse vector is fixed by the global projector. -/
theorem transverseProjector_fixed {k : FreePhotonVec3} {v : FreePhotonVec3}
    (hv : v ∈ TransversePolarization k) :
    transverseProjector k v = v := by
  change dot k v = 0 at hv
  simp [transverseProjector, hv]

/-- Projection is idempotent without selecting a polarization basis. -/
theorem transverseProjector_idempotent {k : FreePhotonVec3} (hk : k ≠ 0)
    (v : FreePhotonVec3) :
    transverseProjector k (transverseProjector k v) =
      transverseProjector k v := by
  exact transverseProjector_fixed (transverseProjector_mem hk v)

/-- Scalar mode actions commute with the transverse projector. -/
theorem transverseProjector_smul {k : FreePhotonVec3} (hk : k ≠ 0)
    (c : ℝ) (v : FreePhotonVec3) :
    transverseProjector k (c • v) = c • transverseProjector k v := by
  have hkk : dot k k ≠ 0 := ne_of_gt (dot_self_pos hk)
  unfold transverseProjector
  rw [dot_smul_right]
  ext d
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  field_simp [hkk]

/-! ## The exact scalar symbol on the transverse fiber -/

/-- The FZ-12 spatial action acts as the same scalar on every transverse
amplitude. -/
noncomputable def photonSpatialAction (a : ℝ)
    (k v : FreePhotonVec3) : FreePhotonVec3 :=
  OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k • v

theorem photonSpatialAction_preserves_transverse
    {a : ℝ} {k : FreePhotonVec3} {v : FreePhotonVec3}
    (hv : v ∈ TransversePolarization k) :
    photonSpatialAction a k v ∈ TransversePolarization k := by
  change dot k
    (OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k • v) = 0
  rw [dot_smul_right]
  change dot k v = 0 at hv
  rw [hv]
  ring

theorem photonSpatialAction_commutes_projector
    (a : ℝ) {k : FreePhotonVec3} (hk : k ≠ 0) (v : FreePhotonVec3) :
    transverseProjector k (photonSpatialAction a k v) =
      photonSpatialAction a k (transverseProjector k v) := by
  exact transverseProjector_smul hk
    (OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k) v

/-- The coordinate-dilated symbol is nonnegative.  This is a spatial
operator statement. -/
theorem dilatedCompletionFourierSymbol_nonnegative
    (a : ℝ) (k : FreePhotonVec3) :
    0 ≤ OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k := by
  unfold OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol
  exact mul_nonneg (div_nonneg (by norm_num) (sq_nonneg a))
    (OPH.SeamCurrentDirichletGenerator.completionFourierSymbol_nonnegative _)

/-- Nonnegative square-root frequency of the declared Hamiltonian mode. -/
noncomputable def photonModeFrequency (a : ℝ) (k : FreePhotonVec3) : ℝ :=
  Real.sqrt
    (OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k)

theorem photonModeFrequency_nonnegative (a : ℝ) (k : FreePhotonVec3) :
    0 ≤ photonModeFrequency a k := Real.sqrt_nonneg _

/-- The Hamiltonian frequency squares to the exact FZ-12 spatial symbol. -/
theorem photonModeFrequency_sq (a : ℝ) (k : FreePhotonVec3) :
    photonModeFrequency a k ^ 2 =
      OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k := by
  exact Real.sq_sqrt (dilatedCompletionFourierSymbol_nonnegative a k)

/-- The exact symbol has no constant mode. -/
theorem dilatedCompletionFourierSymbol_zero (a : ℝ) :
    OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a 0 = 0 := by
  unfold OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol
    OPH.SeamCurrentDirichletGenerator.completionFourierSymbol
  simp [dot]

/-- The declared positive-frequency branch is gapless at zero momentum.  A
physical massless-photon statement additionally needs the open sector and
clock attachments. -/
theorem photonModeFrequency_zero (a : ℝ) :
    photonModeFrequency a 0 = 0 := by
  simp [photonModeFrequency, dilatedCompletionFourierSymbol_zero]

/-- The exact scalar action is degenerate on the full rank-two transverse
fiber. -/
theorem photonSpatialAction_eq_frequency_sq
    (a : ℝ) (k v : FreePhotonVec3) :
    photonSpatialAction a k v = photonModeFrequency a k ^ 2 • v := by
  rw [photonModeFrequency_sq]
  rfl

/-! ## Conditional canonical Hamiltonian mode -/

/-- Displacement and velocity coordinates for one transverse oscillator mode. -/
abbrev PhotonModeState := FreePhotonVec3 × FreePhotonVec3

/-- Canonical first-order mode generator `(A', Pi')=(Pi,-Lambda A)`.
The prime notation is an auxiliary Hamiltonian parameter until a physical
clock is attached. -/
noncomputable def photonModeGenerator
    (a : ℝ) (k : FreePhotonVec3)
    (state : PhotonModeState) : PhotonModeState :=
  (state.2, -photonSpatialAction a k state.1)

/-- Algebraic second displacement component after applying the mode generator
twice. -/
noncomputable def photonModeSecond
    (a : ℝ) (k : FreePhotonVec3)
    (state : PhotonModeState) : FreePhotonVec3 :=
  (photonModeGenerator a k (photonModeGenerator a k state)).1

/-- Exact conditional mode equation `A'' + Lambda A = 0`. -/
theorem photonMode_second_order
    (a : ℝ) (k : FreePhotonVec3) (state : PhotonModeState) :
    photonModeSecond a k state + photonSpatialAction a k state.1 = 0 := by
  simp [photonModeSecond, photonModeGenerator]

/-- Both displacement and velocity lie in the basis-free transverse plane. -/
def ModeIsTransverse (k : FreePhotonVec3) (state : PhotonModeState) : Prop :=
  state.1 ∈ TransversePolarization k ∧
    state.2 ∈ TransversePolarization k

/-- The conditional Hamiltonian generator preserves the transverse
constraint. -/
theorem photonModeGenerator_preserves_transverse
    (a : ℝ) {k : FreePhotonVec3} {state : PhotonModeState}
    (hstate : ModeIsTransverse k state) :
    ModeIsTransverse k (photonModeGenerator a k state) := by
  constructor
  · exact hstate.2
  · change -photonSpatialAction a k state.1 ∈ TransversePolarization k
    exact Submodule.neg_mem _
      (photonSpatialAction_preserves_transverse hstate.1)

/-- Canonical quadratic mode energy of the conditional free-field lift. -/
noncomputable def photonModeEnergy
    (a : ℝ) (k : FreePhotonVec3) (state : PhotonModeState) : ℝ :=
  (dot state.2 state.2 +
    OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k *
      dot state.1 state.1) / 2

/-- Algebraic first variation of the canonical quadratic energy evaluated on
the declared mode generator.  This is not a differentiable-trajectory or
flow-existence theorem. -/
noncomputable def photonModeEnergyFirstVariation
    (a : ℝ) (k : FreePhotonVec3) (state : PhotonModeState) : ℝ :=
  dot state.2 (-photonSpatialAction a k state.1) +
    OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k *
      dot state.1 state.2

/-- The algebraic first variation vanishes on the declared generator.  Energy
conservation along a trajectory requires a separately constructed flow. -/
theorem photonModeEnergy_firstVariation_generator_zero
    (a : ℝ) (k : FreePhotonVec3) (state : PhotonModeState) :
    photonModeEnergyFirstVariation a k state = 0 := by
  unfold photonModeEnergyFirstVariation photonSpatialAction
  rw [show
      -(OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k •
        state.1) =
        (-OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k) •
          state.1 by simp,
    dot_smul_right, dot_comm state.2 state.1]
  ring

/-! ## Axiom audit

The construction uses the exact FZ-12 completion symbol and ordinary
finite-dimensional real linear algebra.  It adds no project axiom.  Its
physical interpretation consumes the open position, clock, field-sector,
frame, scale, gluing, interaction, source-selection, and readout premises.
-/

#print axioms transversePolarization_finrank
#print axioms transverseProjector_dot_zero
#print axioms transverseProjector_idempotent
#print axioms photonSpatialAction_commutes_projector
#print axioms photonModeFrequency_sq
#print axioms photonModeFrequency_zero
#print axioms photonMode_second_order
#print axioms photonModeGenerator_preserves_transverse
#print axioms photonModeEnergy_firstVariation_generator_zero

end OPH.SeamCurrentFreePhotonLift
