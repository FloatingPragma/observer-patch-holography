import LightSignalAdequacySurface
import ModalMaxwellFactorizationBoundary

open scoped BigOperators Matrix

namespace OPH.LightSignalMaxwellComposition

open OPH.LightSignalAdequacySurface
open OPH.SeamCurrentHomogeneousAction
open OPH.SeamCurrentDirichletGenerator
open OPH.SeamCurrentCarrierQuotient
open OPH.SeamCurrentFreePhotonLift
open OPH.ModalMaxwellFactorizationBoundary
open OPH.PrimitivePortTranslationBridge
open OPH.CarrierFrequencySpeed

/-!
# Maxwell and Gauss composition of the light-signal surface (issue #733)

Two composed receipts at the one premise bundle
`LightSignalPremiseData` of the committed light-signal adequacy
surface, and one final conjunction carrying the surface receipt, the
Maxwell-shaped receipt, and the Gauss-threaded receipt at the same
bundle.  The bundle fields carry the registered premise rows exactly as
the surface declares them: the equal source-counting measure (register
row PR-20, field `sourceCountingSelection`), the auxiliary oscillator
lift (register row PR-21, fields `modeDynamics` and
`modeDynamics_is_lift`), and the physical-frequency identification at
the supplied chart scale (register row PR-22, fields `chartScale`,
`chartScale_pos`, `physicalFrequency`, `physicalFrequency_nonneg`,
`physicalFrequency_sq_eq_symbol`).  No premise structure is added or
duplicated here.

The Gauss-threaded receipt (`GaussThreadedClauses`,
`lightSignalGaussReceipt`).  The PR-20-forced measure is
orientation-balanced across the committed double cover of the thirty
seam-table rows, and its expected signed port-boundary readout
vanishes: the sum over all sixty directed seams of the selected weight
times the signed incidence boundary is zero.  The boundary object in
that sum, `directedBoundary`, is the same incidence object whose
rational scalar extension satisfies the committed Gauss receipts, and
two exhibition clauses state the identification: the seam-table atom
boundary equals the directed boundary of the canonical orientation
(`seamBoundary_atom`), and the canonical record step equals the
thirty-current readout composite (`directedSeamStep_canonical_control`).
The four rational Gauss receipts (solvable exactly on neutral loads,
declared spanning-tree section, fibre equals one cycle translate, cycle
rank nineteen) enter the same receipt as premise-free exact clauses
about that shared incidence object.  Only the balance and neutrality
clauses consume PR-20; the exhibition clauses and the Gauss clauses
consume no bundle field, and the receipt conjunction carries both at
one bundle.  The Gauss clauses stay over `ℚ`; the neutrality sum is
over `ℝ` with an integer cast, and no rational-to-real identification
of the Gauss receipts is asserted.  The threaded object is the seam
incidence alphabet and its boundary readouts.  This module commits no
Coulomb readout, no discrete Green function, no source-produced gauge
field, no conserved physical current, and no continuum composition
(PR-54 open).

The Maxwell-shaped receipt (`MaxwellCompositionClauses`,
`lightSignalMaxwellReceipt`).  At the same bundle: the selected
generator is the dimensionless Dirichlet generator (PR-20); the
supplied mode dynamics is carried by the explicit boundary map
`modalMaxwellBoundaryMap` onto the committed opposite-sign curl pairing
`maxwellShapedModalGenerator`, exactly, on every transverse real state
at every nonzero chart momentum (PR-21, PR-22); the pairing's second
iterate obeys the second-order wave law with the square of the supplied
readout frequency on both amplitudes (PR-22); the pairing preserves
both modal divergence constraints; the polarization rank is two, the
identified frequency has an exact zero mode and a global one-Lipschitz
bound (PR-22); the pairing commutes with the global complex phase
action on modal amplitudes; and the same-sign curl mutation fails the
wave law at the same identified frequency.

Boundary.  Every statement in this module is finite and modal: fixed
chart momentum in the seam-built symbol, complex modal amplitudes,
no position-space field.  The carrying statement is an intertwining
through the named boundary map, and no equivalence of dynamics is
claimed: the boundary map has nontrivial kernel at nonzero chart
momenta where the identified frequency vanishes (reciprocal-lattice
zeros of the dilated symbol).  No real-coefficient curl intertwiner
with the committed sign pattern exists; the complexification is part
of the committed pairing, and no real-field clause is stated.  The
phase clause is global complex-scalar equivariance of the pairing and
is not a gauge-symmetry statement.  This module commits no U(1)
potential, no gauge quotient, no action, no opposite-momentum reality
pairing, no source coupling, no boost covariance, no continuum
control, and no identification of either amplitude with a laboratory
field (PR-53 and PR-54 open).  The dispersion rows keep their custody
with the frozen instruments; this module restates no dispersion
coefficient and no decision rule.
-/

/-! ## Gauss threading under the register (PR-20)

The clauses of this section bind the PR-20-forced measure to the same
signed incidence object that carries the committed rational Gauss
receipts.  The balance and neutrality clauses consume the bundled
selection data; the exhibition and Gauss clauses are premise-free. -/

/-- The canonical orientation of a seam-table row is the `false` label
of the committed sixty-label presentation.  Definitional unfolding of
`indexedDirectedSeam`. -/
theorem indexedDirectedSeam_false (e : Fin 30) :
    indexedDirectedSeam (e, false) = canonicalDirectedSeam e := rfl

/-- The reversed orientation of a seam-table row is the `true` label of
the committed sixty-label presentation.  Definitional unfolding of
`indexedDirectedSeam`. -/
theorem indexedDirectedSeam_true (e : Fin 30) :
    indexedDirectedSeam (e, true) =
      reverseDirectedSeam (canonicalDirectedSeam e) := rfl

/-- **Orientation balance of the forced measure (PR-20).**  The
selected weight of every seam-table row equals the selected weight of
its reversed orientation.  Both weights are the forced `1/60` of the
committed A2/A3 selection theorem at the bundled data. -/
theorem selected_weight_orientation_balanced (D : LightSignalPremiseData) :
    ∀ e : Fin 30,
      D.sourceCountingSelection.selected (indexedDirectedSeam (e, false)) =
        D.sourceCountingSelection.selected (indexedDirectedSeam (e, true)) := by
  intro e
  rw [a2a3_directed_seam_weight_eq_one_sixtieth D.sourceCountingSelection,
    a2a3_directed_seam_weight_eq_one_sixtieth D.sourceCountingSelection]

/-- **Expected signed boundary of the forced measure vanishes
(PR-20).**  The sum over all sixty directed seams of the selected
weight times the signed incidence boundary is the zero port load.  The
boundary object is the same `directedBoundary` whose rational
extension carries the committed Gauss receipts; reversed orientation
pairs cancel under the forced equal weights.  This is a finite
statement about the seam alphabet, not about a physical current. -/
theorem selected_expected_boundary_zero (D : LightSignalPremiseData) :
    (∑ e : DirectedSeam,
      D.sourceCountingSelection.selected e •
        fun p : Fin 12 ↦ ((directedBoundary e.1.1 e.1.2 1 p : ℤ) : ℝ)) = 0 := by
  funext p
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply]
  have hreindex :
      (∑ i : Fin 30 × Bool,
        D.sourceCountingSelection.selected (indexedDirectedSeam i) *
          ((directedBoundary (indexedDirectedSeam i).1.1
            (indexedDirectedSeam i).1.2 1 p : ℤ) : ℝ)) =
      ∑ e : DirectedSeam,
        D.sourceCountingSelection.selected e *
          ((directedBoundary e.1.1 e.1.2 1 p : ℤ) : ℝ) :=
    indexedDirectedSeam_bijective.sum_comp
      (fun e : DirectedSeam ↦ D.sourceCountingSelection.selected e *
        ((directedBoundary e.1.1 e.1.2 1 p : ℤ) : ℝ))
  rw [← hreindex, Fintype.sum_prod_type]
  apply Finset.sum_eq_zero
  intro e _
  rw [Fintype.sum_bool]
  simp only [indexedDirectedSeam_true, indexedDirectedSeam_false,
    a2a3_directed_seam_weight_eq_one_sixtieth D.sourceCountingSelection]
  show (1 / 60 : ℝ) *
      ((directedBoundary (seamRight e) (seamLeft e) 1 p : ℤ) : ℝ) +
      (1 / 60 : ℝ) *
        ((directedBoundary (seamLeft e) (seamRight e) 1 p : ℤ) : ℝ) = 0
  have hcancel : directedBoundary (seamRight e) (seamLeft e) 1 p +
      directedBoundary (seamLeft e) (seamRight e) 1 p = 0 := by
    simp only [directedBoundary]
    ring
  have hcast : ((directedBoundary (seamRight e) (seamLeft e) 1 p : ℤ) : ℝ) +
      ((directedBoundary (seamLeft e) (seamRight e) 1 p : ℤ) : ℝ) = 0 := by
    rw [← Int.cast_add, hcancel, Int.cast_zero]
  linarith

/-- **Shared incidence object.**  The port boundary of one canonical
seam atom equals the directed boundary of the canonical orientation of
the same seam-table row.  The PR-20 alphabet and the Gauss incidence
operator are built from the same `seamLeft`/`seamRight` table. -/
theorem canonicalDirectedSeam_boundary_shared (e : Fin 30) :
    seamBoundary (seamAtom e 1) =
      directedBoundary (canonicalDirectedSeam e).1.1
        (canonicalDirectedSeam e).1.2 1 :=
  seamBoundary_atom e 1

/-- **The Gauss-threaded clauses at one premise bundle.**  Clause list,
with the consuming register row named per clause: orientation balance
of the selected measure (PR-20); vanishing expected signed boundary of
the selected measure (PR-20); the seam-atom/directed-boundary
exhibition (premise-free); the record-step/thirty-current exhibition
(premise-free); solvability exactly on neutral loads, the declared
spanning-tree section, the solution fibre as one cycle translate, and
cycle rank nineteen (premise-free rational Gauss receipts).  The
conjunction carries the PR-20 clauses and the Gauss clauses at one
bundle, threaded through the exhibited shared incidence objects. -/
def GaussThreadedClauses (D : LightSignalPremiseData) : Prop :=
  (∀ e : Fin 30,
    D.sourceCountingSelection.selected (indexedDirectedSeam (e, false)) =
      D.sourceCountingSelection.selected (indexedDirectedSeam (e, true)))
  ∧ ((∑ e : DirectedSeam,
      D.sourceCountingSelection.selected e •
        fun p : Fin 12 ↦ ((directedBoundary e.1.1 e.1.2 1 p : ℤ) : ℝ)) = 0)
  ∧ (∀ e : Fin 30,
      seamBoundary (seamAtom e 1) =
        directedBoundary (canonicalDirectedSeam e).1.1
          (canonicalDirectedSeam e).1.2 1)
  ∧ (∀ e : Fin 30,
      (directedSeamStep (canonicalDirectedSeam e)).control =
        seamAxisCurrent (seamAtom e 1))
  ∧ (∀ charge : RationalPortLoad,
      (∃ field : RationalSeamCurrent,
          rationalSeamBoundary field = charge) ↔
        rationalPortTotal charge = 0)
  ∧ (∀ charge : RationalPortLoad,
      rationalPortTotal charge = 0 →
        rationalSeamBoundary (rationalBoundarySection charge) = charge)
  ∧ (∀ (field field' : RationalSeamCurrent) (charge : RationalPortLoad),
      rationalSeamBoundary field = charge →
        (rationalSeamBoundary field' = charge ↔
          field' - field ∈ LinearMap.ker rationalSeamBoundary))
  ∧ Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) = 19

/-- **The Gauss-threaded receipt (issue #733).**  Every clause of
`GaussThreadedClauses` holds at every premise bundle.  The receipt
scope is the seam incidence alphabet and its boundary readouts; no
Coulomb readout, discrete Green function, source-produced gauge field,
physical current, or continuum composition is committed (PR-54
open). -/
theorem lightSignalGaussReceipt (D : LightSignalPremiseData) :
    GaussThreadedClauses D :=
  ⟨selected_weight_orientation_balanced D,
    selected_expected_boundary_zero D,
    canonicalDirectedSeam_boundary_shared,
    fun e ↦ directedSeamStep_canonical_control e,
    OPH.DiscreteGauss.gauss_solution_exists_iff_total_zero,
    fun charge hneutral ↦
      OPH.DiscreteGauss.rationalBoundarySection_is_gauss_solution
        charge hneutral,
    fun _ _ _ hfield ↦
      OPH.DiscreteGauss.gauss_solution_iff_difference_is_cycle hfield,
    OPH.DiscreteGauss.gauss_cycle_space_finrank⟩

/-! ## The composed modal Maxwell receipt (PR-20, PR-21, PR-22)

The supplied mode dynamics acts on real displacement/velocity pairs;
the committed Maxwell-shaped pairing acts on complex amplitude pairs.
The bridge is the explicit boundary map below.  The intertwining is an
exact equality on transverse states at nonzero chart momentum; it is
not an equivalence of dynamics, because the boundary map has
nontrivial kernel at nonzero momenta where the identified frequency
vanishes. -/

/-- Scalar extension of one real modal vector preserves the transverse
constraint: the complex modal divergence of the complexified vector is
the cast of the real momentum pairing. -/
theorem complexifyModalVector_transverse {k v : FreePhotonVec3}
    (hv : v ∈ TransversePolarization k) :
    ComplexModeIsTransverse k (complexifyModalVector v) := by
  have hdot : dot k v = 0 := hv
  unfold ComplexModeIsTransverse complexMomentumDot
  have hcast : complexMomentum k ⬝ᵥ complexifyModalVector v =
      ((dot k v : ℝ) : ℂ) := by
    simp only [complexMomentum, complexifyModalVector, dot, dotProduct]
    push_cast
    rfl
  rw [hcast, hdot, Complex.ofReal_zero]

/-- Scalar extension commutes with negation. -/
theorem complexifyModalVector_neg (v : FreePhotonVec3) :
    complexifyModalVector (-v) = -complexifyModalVector v := by
  funext d
  simp [complexifyModalVector]

/-- The boundary map carrying one real oscillator state onto the
complex Maxwell-shaped carrier: the velocity component is complexified
into the first amplitude, and minus the Fourier curl of the
complexified displacement is the second amplitude.  The map is modal
and finite; it identifies neither amplitude with a laboratory field. -/
noncomputable def modalMaxwellBoundaryMap (a : ℝ) (k : FreePhotonVec3)
    (state : PhotonModeState) : ComplexPairedModalState :=
  (complexifyModalVector state.2,
    -fourierCurl a k (complexifyModalVector state.1))

/-- **The supplied dynamics is carried onto the Maxwell-shaped pairing
(PR-21, PR-22).**  At every nonzero chart momentum and every transverse
real state, the boundary map intertwines the bundled mode dynamics with
the committed opposite-sign curl pairing, exactly.  This is an
intertwining through `modalMaxwellBoundaryMap`, not an equivalence of
dynamics: the boundary map has nontrivial kernel at nonzero momenta
where the identified frequency vanishes. -/
theorem modeDynamics_intertwines_maxwellShaped (D : LightSignalPremiseData)
    {k : FreePhotonVec3} (hk : k ≠ 0) {state : PhotonModeState}
    (hstate : ModeIsTransverse k state) :
    modalMaxwellBoundaryMap D.chartScale k (D.modeDynamics k state) =
      maxwellShapedModalGenerator D.chartScale k
        (modalMaxwellBoundaryMap D.chartScale k state) := by
  rw [D.modeDynamics_is_lift k state]
  apply Prod.ext
  · show complexifyModalVector
        (-photonSpatialAction D.chartScale k state.1) =
      fourierCurl D.chartScale k
        (-fourierCurl D.chartScale k (complexifyModalVector state.1))
    rw [complexifyModalVector_neg, fourierCurl_neg,
      fourierCurl_sq_on_transverse hk
        (complexifyModalVector_transverse hstate.1),
      complexPhotonSpatialAction_complexifies]
  · rfl

/-- **Second-order wave law of the pairing at the identified frequency
(PR-22).**  The second iterate of the Maxwell-shaped pairing at the
bundled chart scale obeys the wave equation with the square of the
supplied readout frequency on both amplitudes, on every complex
transverse pair at every nonzero chart momentum. -/
theorem maxwellShaped_second_order_identified (D : LightSignalPremiseData)
    {k : FreePhotonVec3} (hk : k ≠ 0) (state : ComplexPairedModalState)
    (h1 : ComplexModeIsTransverse k state.1)
    (h2 : ComplexModeIsTransverse k state.2) :
    (maxwellShapedModalGenerator D.chartScale k
        (maxwellShapedModalGenerator D.chartScale k state)).1 +
        ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 = 0 ∧
      (maxwellShapedModalGenerator D.chartScale k
          (maxwellShapedModalGenerator D.chartScale k state)).2 +
        ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.2 = 0 := by
  have h := maxwellShapedModalGenerator_sq_wave (a := D.chartScale)
    hk state ⟨h1, h2⟩
  rw [D.physicalFrequency_eq_mode_branch k]
  exact ⟨h.1, h.2⟩

/-- **Global phase equivariance of the pairing.**  The Maxwell-shaped
pairing commutes with every complex scalar action on amplitude pairs;
the unit-circle restriction of the scalar is the global phase reading.
This is equivariance of one modal map and is not a gauge-symmetry
statement; no U(1) potential or gauge quotient is stated in this
module. -/
theorem maxwellShaped_phase_equivariant (a : ℝ) (k : FreePhotonVec3)
    (c : ℂ) (state : ComplexPairedModalState) :
    maxwellShapedModalGenerator a k (c • state) =
      c • maxwellShapedModalGenerator a k state := by
  apply Prod.ext
  · show fourierCurl a k (c • state).2 = c • fourierCurl a k state.2
    rw [Prod.smul_snd, fourierCurl_smul]
  · show -fourierCurl a k (c • state).1 = c • (-fourierCurl a k state.1)
    rw [Prod.smul_fst, fourierCurl_smul, smul_neg]

/-- **Sign forcing at the identified frequency (PR-22).**  The
same-sign curl mutation fails the second-order wave law with the
supplied readout frequency on every transverse state whose
frequency-scaled first amplitude is nonzero. -/
theorem sameSignMutation_fails_identified_wave (D : LightSignalPremiseData)
    {k : FreePhotonVec3} (hk : k ≠ 0) (state : ComplexPairedModalState)
    (h1 : ComplexModeIsTransverse k state.1)
    (h2 : ComplexModeIsTransverse k state.2)
    (haction : ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 ≠ 0) :
    (sameSignCurlMutation D.chartScale k
        (sameSignCurlMutation D.chartScale k state)).1 +
      ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 ≠ 0 := by
  rw [D.physicalFrequency_eq_mode_branch k] at haction ⊢
  exact sameSignCurlMutation_fails_wave hk state ⟨h1, h2⟩ haction

/-- **The Maxwell-shaped clauses at one premise bundle.**  Clause list,
with the consuming register rows named per clause: the selected
generator is the dimensionless Dirichlet generator (PR-20); the
supplied mode dynamics is carried by the boundary map onto the
opposite-sign curl pairing on transverse states at nonzero momentum
(PR-21, PR-22); the pairing's second iterate obeys the wave law with
the square of the supplied readout frequency (PR-22); the pairing
preserves both modal divergence constraints (PR-22, chart scale only);
the polarization rank is two (premise-free), the identified frequency
has an exact zero mode and a global one-Lipschitz bound (PR-22); the
pairing commutes with the global complex phase action (PR-22, chart
scale only); and the same-sign mutation fails the wave law at the
identified frequency (PR-22).  Every clause is finite and modal. -/
def MaxwellCompositionClauses (D : LightSignalPremiseData) : Prop :=
  (selectedCompletionGenerator D.sourceCountingSelection =
      completionDirichletGenerator)
  ∧ (∀ k : FreePhotonVec3, k ≠ 0 →
      ∀ state : PhotonModeState, ModeIsTransverse k state →
        modalMaxwellBoundaryMap D.chartScale k (D.modeDynamics k state) =
          maxwellShapedModalGenerator D.chartScale k
            (modalMaxwellBoundaryMap D.chartScale k state))
  ∧ (∀ k : FreePhotonVec3, k ≠ 0 →
      ∀ state : ComplexPairedModalState,
        ComplexModeIsTransverse k state.1 →
          ComplexModeIsTransverse k state.2 →
            (maxwellShapedModalGenerator D.chartScale k
                (maxwellShapedModalGenerator D.chartScale k state)).1 +
                ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 = 0 ∧
              (maxwellShapedModalGenerator D.chartScale k
                  (maxwellShapedModalGenerator D.chartScale k state)).2 +
                ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.2 = 0)
  ∧ (∀ (k : FreePhotonVec3) (state : ComplexPairedModalState),
      ComplexModeIsTransverse k
          (maxwellShapedModalGenerator D.chartScale k state).1 ∧
        ComplexModeIsTransverse k
          (maxwellShapedModalGenerator D.chartScale k state).2)
  ∧ ((∀ k : FreePhotonVec3, k ≠ 0 →
        Module.finrank ℝ (TransversePolarization k) = 2)
      ∧ D.physicalFrequency 0 = 0
      ∧ (∀ k p : FreePhotonVec3,
          |D.physicalFrequency k - D.physicalFrequency p| ≤
            euclideanMagnitude (k - p)))
  ∧ (∀ (k : FreePhotonVec3) (c : ℂ) (state : ComplexPairedModalState),
      maxwellShapedModalGenerator D.chartScale k (c • state) =
        c • maxwellShapedModalGenerator D.chartScale k state)
  ∧ (∀ k : FreePhotonVec3, k ≠ 0 →
      ∀ state : ComplexPairedModalState,
        ComplexModeIsTransverse k state.1 →
          ComplexModeIsTransverse k state.2 →
            ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 ≠ 0 →
              (sameSignCurlMutation D.chartScale k
                  (sameSignCurlMutation D.chartScale k state)).1 +
                ((D.physicalFrequency k ^ 2 : ℝ) : ℂ) • state.1 ≠ 0)

/-- **The composed modal Maxwell receipt (issue #733).**  Every clause
of `MaxwellCompositionClauses` holds at every premise bundle.  All
clauses are stated about the bundled objects themselves: the supplied
selection data, the supplied mode dynamics, the supplied chart scale,
and the supplied readout frequency.  The receipt is finite and modal
and claims no position-space field, no reality pairing, no potential,
no action, no source coupling, no covariance, and no continuum
control (PR-53 and PR-54 open). -/
theorem lightSignalMaxwellReceipt (D : LightSignalPremiseData) :
    MaxwellCompositionClauses D :=
  ⟨D.selected_generator_eq_dirichlet,
    fun _ hk _ hstate ↦
      modeDynamics_intertwines_maxwellShaped D hk hstate,
    fun _ hk state h1 h2 ↦
      maxwellShaped_second_order_identified D hk state h1 h2,
    fun k state ↦
      maxwellShapedModalGenerator_transverse D.chartScale k state,
    ⟨fun _ hk ↦ transversePolarization_finrank hk,
      D.physicalFrequency_zero_mode,
      fun k p ↦ D.physicalFrequency_global_one_lipschitz k p⟩,
    fun k c state ↦
      maxwellShaped_phase_equivariant D.chartScale k c state,
    fun _ hk state h1 h2 haction ↦
      sameSignMutation_fails_identified_wave D hk state h1 h2 haction⟩

/-! ## The one citable composed statement -/

/-- **The composed light-signal receipt with Maxwell and Gauss
threading (issue #733).**  For every premise bundle, one conjunction:
the committed light-signal adequacy surface receipt, the composed
modal Maxwell receipt, and the Gauss-threaded receipt, all at the same
bundle.  The first conjunct restates the surface receipt clause for
clause; the second and third are `MaxwellCompositionClauses` and
`GaussThreadedClauses`.  The conjunction composes committed statements
at the bundled premises and claims nothing beyond them. -/
theorem lightSignalMaxwellComposition_receipt (D : LightSignalPremiseData) :
    ((selectedCompletionGenerator D.sourceCountingSelection =
        completionDirichletGenerator)
      ∧ (∀ (k : FreePhotonVec3) (state : PhotonModeState),
          (D.modeDynamics k (D.modeDynamics k state)).1 +
              D.physicalFrequency k ^ 2 • state.1 = 0)
      ∧ (∀ (k : FreePhotonVec3) (state : PhotonModeState),
          ModeIsTransverse k state →
            ModeIsTransverse k (D.modeDynamics k state))
      ∧ (∀ k : FreePhotonVec3, k ≠ 0 →
          Module.finrank ℝ (TransversePolarization k) = 2)
      ∧ D.physicalFrequency 0 = 0
      ∧ (∀ k p : FreePhotonVec3,
          |D.physicalFrequency k - D.physicalFrequency p| ≤
            euclideanMagnitude (k - p))
      ∧ (∀ charge : RationalPortLoad,
          (∃ field : RationalSeamCurrent,
              rationalSeamBoundary field = charge) ↔
            rationalPortTotal charge = 0)
      ∧ Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) = 19)
    ∧ MaxwellCompositionClauses D
    ∧ GaussThreadedClauses D :=
  ⟨lightSignalAdequacySurface_receipt D, lightSignalMaxwellReceipt D,
    lightSignalGaussReceipt D⟩

end OPH.LightSignalMaxwellComposition

/-! ## Axiom audit

Statements that touch the sixty-seam alphabet or consume the bundled
committed chain (the orientation balance, the expected-boundary sum,
the selected-generator clause, the intertwining, second-order, and
sign-forcing clauses through the supplied dynamics and frequency, and
the receipts containing them) inherit the disclosed native-decision
certificates of the committed `SeamCurrentHomogeneousAction`; their
audit lines show the corresponding native-decide axioms.  The
definitional exhibition lemmas, the scalar-extension lemmas, the
seam-atom boundary identity, and the phase clause show only `propext`,
`Classical.choice`, and `Quot.sound`. -/

#print axioms OPH.LightSignalMaxwellComposition.indexedDirectedSeam_false
#print axioms OPH.LightSignalMaxwellComposition.indexedDirectedSeam_true
#print axioms OPH.LightSignalMaxwellComposition.selected_weight_orientation_balanced
#print axioms OPH.LightSignalMaxwellComposition.selected_expected_boundary_zero
#print axioms OPH.LightSignalMaxwellComposition.canonicalDirectedSeam_boundary_shared
#print axioms OPH.LightSignalMaxwellComposition.lightSignalGaussReceipt
#print axioms OPH.LightSignalMaxwellComposition.complexifyModalVector_transverse
#print axioms OPH.LightSignalMaxwellComposition.complexifyModalVector_neg
#print axioms OPH.LightSignalMaxwellComposition.modeDynamics_intertwines_maxwellShaped
#print axioms OPH.LightSignalMaxwellComposition.maxwellShaped_second_order_identified
#print axioms OPH.LightSignalMaxwellComposition.maxwellShaped_phase_equivariant
#print axioms OPH.LightSignalMaxwellComposition.sameSignMutation_fails_identified_wave
#print axioms OPH.LightSignalMaxwellComposition.lightSignalMaxwellReceipt
#print axioms OPH.LightSignalMaxwellComposition.lightSignalMaxwellComposition_receipt
