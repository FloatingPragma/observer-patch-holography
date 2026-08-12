import SeamCurrentAuxiliaryOscillatorLift
import SeamCurrentFreePhotonLift
import CarrierFrequencySpeed
import DiscreteGauss

namespace OPH.LightSignalAdequacySurface

open OPH.SeamCurrentHomogeneousAction
open OPH.SeamCurrentDirichletGenerator
open OPH.SeamCurrentAuxiliaryOscillatorLift
open OPH.SeamCurrentFreePhotonLift
open OPH.CarrierFrequencySpeed
open OPH.SeamCurrentCarrierQuotient

/-!
# The composed light-signal adequacy surface (V3, issue #733)

One citable package for the observation-ledger rows OL-F1, OL-F2, and
OL-F3 of the electromagnetism lane: the source-selected Dirichlet
generator, the basis-free transverse completion with its exact zero
mode, the declared oscillator dynamics with the identified physical
frequency and its exact unit speed bound, and the exact rational Gauss
receipts.  The package composes the committed seam-current modules with
every supplied object named.  Apart from the four small composition
theorems listed below, every theorem re-exports a committed result,
either through the premise structure `LightSignalPremiseData` or
premise-free in the side-by-side Gauss and generator sections.

`LightSignalPremiseData` bundles the three declared premises of the V3
program that the committed chain consumes.

The equal source-counting measure.  The seam move law is selected by an
A2-natural, A3-unique projection on the sixty directed source seams,
supplied as the committed structure `A2A3DirectedSeamProjection` (the
feasible move simplex, its natural objective, and its unique selected
minimizer).  Given that data, the committed theorems force the equal
weight `1/60` and the one dimensionless Dirichlet generator.  No source
theorem produces the projection data itself.

The auxiliary oscillator lift.  The dynamics of one transverse mode is
declared to be the first-order lift `(A', Pi') = (Pi, -Lambda A)` of the
spatial symbol.  The bundle carries the supplied dynamics as a field
together with the identification stating that it equals the committed
lift `photonModeGenerator`.  The committed record derives the spatial
generator; the second-order evolution shape enters as this declaration,
and no source theorem produces it.

The physical-frequency identification.  A nonnegative readout frequency
is supplied together with the identification of its square with the
coordinate-dilated carrier symbol at a supplied strictly positive chart
scale.  The chart scale is a declared coordinate; reading it as a
physical length and the readout as a laboratory frequency consumes the
separate calibration anchor rows.

Contents.  On the generator side: the forced equal weight, the equality
of the selected average with the source-counting Markov average, the
equality of the selected generator with the dimensionless Dirichlet
generator, the exact carre-du-champ identity, the plane-wave eigenvalue,
and the normalization of the symbol to the edge-current character.  On
the oscillator side: the scalar lift equation `q'' + L q = 0` with its
plane-wave dispersion, the second-order mode equation of the supplied
dynamics, transversality preservation, and the vanishing first variation
of the quadratic mode energy.  On the transverse side: rank two of the
polarization space at every nonzero momentum, the exact zero mode, and
the projector compatibility of the spatial action.  On the Gauss side:
solvability exactly on neutral loads, the declared spanning-tree
solution, the solution fibre as one affine translate of the cycle
kernel, and the exact rank nineteen of that kernel, which is the
committed source-free solution ambiguity of the discrete receipts (the
discrete analogue of a gauge ambiguity; no gauge potential or quotient
is committed).

Composition theorems proved here, each a short consequence of the
committed statements under the bundled premises: the supplied readout
frequency equals the committed nonnegative branch
(`physicalFrequency_eq_mode_branch`); it equals the exact edge-30
feature norm and is therefore globally one-Lipschitz in the chart
momentum (`physicalFrequency_eq_edge30_frequency`,
`physicalFrequency_global_one_lipschitz`); and the supplied mode
dynamics obeys the second-order wave-shaped evolution with the
identified frequency on every mode (`transverse_wave_second_order`).

The OL-F1 boundary.  The committed construction is a conditional
free-vector model.  The definition `attainedCommittedOrderLaw` states
the complete equation-of-motion content attained at the committed
order: on the rank-two transverse sector the supplied dynamics preserves
transversality and its second iterate obeys
`A'' + omega^2 A = 0` with the identified frequency, mode by mode in
the chart.  Maxwell-shaped field equations contain, in addition: a
coupled first-order evolution of two transverse fields through a curl
pairing produced by the source; the position-space assembly of the mode
fibres into one field on the carrier completion, with a gauge potential
and its quotient; inhomogeneous equations coupling the field to a
conserved source beyond the exact rational Gauss receipts; boost
covariance of the coupled system; and continuum control of the finite
symbol.  None of these statements has a committed producer, and this
module states none of them as a theorem.

The Gauss receipts consume none of the three bundled premises.  They
are exact rational statements about the seam incidence operator, and
they sit side by side inside this module for the composed OL-F3 row.
The same-metric physical attachment record of
`SeamCurrentPhysicalMetricAttachment` is a separate conditional route
and is not consumed here.  The dispersion rows of OL-F4 keep their
custody with the frozen instruments; this module restates no dispersion
coefficient and no decision rule.

Axiom audit.  Every theorem re-exports or composes committed results;
the module adds no project axiom.  Statements that touch the seam
alphabet or the seam-built symbol inherit the disclosed native-decision
certificates of the committed `SeamCurrentHomogeneousAction` (the exact
sixty-seam count, the port-map adjacency, and the orbit transitivity),
as the sibling modules `SeamCurrentDirichletGenerator` and
`SeamCurrentAuxiliaryOscillatorLift` record; those audit lines show the
corresponding native-decide axioms.  The scalar lift statements, the
normalized symbol identity, the transverse rank, and the Gauss receipts
show only `propext`, `Classical.choice`, and `Quot.sound`.
-/

/-- The supplied objects of the composed light-signal surface: the
three declared premises of the V3 electromagnetism lane.  Each field
carries one declared object or one identification; no field is produced
by a source theorem in the imports. -/
structure LightSignalPremiseData where
  /-- The equal source-counting measure premise: A2-natural, A3-unique
  selection data on the sixty directed source seams.  The committed
  theorems force the equal weight `1/60` from this data. -/
  sourceCountingSelection : A2A3DirectedSeamProjection
  /-- The supplied chart scale of the coordinate dilation.  It is a
  declared coordinate; reading it as a physical length consumes the
  separate calibration anchors. -/
  chartScale : ℝ
  /-- The supplied chart scale is strictly positive. -/
  chartScale_pos : 0 < chartScale
  /-- The supplied dynamics of one transverse mode. -/
  modeDynamics : FreePhotonVec3 → PhotonModeState → PhotonModeState
  /-- The auxiliary oscillator lift premise: the supplied mode dynamics
  is the declared first-order lift `(A', Pi') = (Pi, -Lambda A)` of the
  spatial symbol.  The second-order evolution shape enters through this
  declaration. -/
  modeDynamics_is_lift :
    ∀ (k : FreePhotonVec3) (state : PhotonModeState),
      modeDynamics k state = photonModeGenerator chartScale k state
  /-- The supplied nonnegative readout frequency. -/
  physicalFrequency : FreePhotonVec3 → ℝ
  /-- The supplied readout frequency is nonnegative. -/
  physicalFrequency_nonneg : ∀ k : FreePhotonVec3,
    0 ≤ physicalFrequency k
  /-- The physical-frequency identification premise: the square of the
  supplied readout frequency is the coordinate-dilated carrier symbol
  at the supplied chart scale. -/
  physicalFrequency_sq_eq_symbol : ∀ k : FreePhotonVec3,
    physicalFrequency k ^ 2 = dilatedCompletionFourierSymbol chartScale k

namespace LightSignalPremiseData

variable (D : LightSignalPremiseData)

/-! ## The source-selected generator (OL-F1 input) -/

/-- **The forced equal weight.**  A2 covariance and A3 uniqueness on the
exact sixty-seam orbit force one sixtieth on every directed seam.
Re-export of `a2a3_directed_seam_weight_eq_one_sixtieth` at the bundled
selection data. -/
theorem source_counting_weight_forced (e : DirectedSeam) :
    D.sourceCountingSelection.selected e = 1 / 60 :=
  a2a3_directed_seam_weight_eq_one_sixtieth D.sourceCountingSelection e

/-- The average selected by the bundled data is the equal
source-counting Markov average.  Re-export of
`a2a3_selected_completion_average_eq_markov`. -/
theorem selected_average_eq_source_counting :
    weightedCompletionAverage D.sourceCountingSelection.selected =
      completionMarkovAverage :=
  a2a3_selected_completion_average_eq_markov D.sourceCountingSelection

/-- **The exact Dirichlet generator.**  The generator built from the
bundled selection data is the one dimensionless completion Dirichlet
generator.  Re-export of
`a2a3_selected_completion_generator_eq_dirichlet`. -/
theorem selected_generator_eq_dirichlet :
    selectedCompletionGenerator D.sourceCountingSelection =
      completionDirichletGenerator :=
  a2a3_selected_completion_generator_eq_dirichlet D.sourceCountingSelection

/-! ## The declared mode dynamics (OL-F1 and OL-F2 inputs) -/

/-- The supplied mode dynamics applied twice gives the exact
second-order equation `A'' + Lambda A = 0` against the spatial symbol.
Composition of the lift declaration with the committed
`photonMode_second_order`. -/
theorem modeDynamics_second_order (k : FreePhotonVec3)
    (state : PhotonModeState) :
    (D.modeDynamics k (D.modeDynamics k state)).1 +
        photonSpatialAction D.chartScale k state.1 = 0 := by
  simp only [D.modeDynamics_is_lift]
  exact photonMode_second_order D.chartScale k state

/-- The supplied mode dynamics preserves the basis-free transverse
constraint.  Composition of the lift declaration with the committed
`photonModeGenerator_preserves_transverse`. -/
theorem modeDynamics_preserves_transverse {k : FreePhotonVec3}
    {state : PhotonModeState} (hstate : ModeIsTransverse k state) :
    ModeIsTransverse k (D.modeDynamics k state) := by
  rw [D.modeDynamics_is_lift k state]
  exact photonModeGenerator_preserves_transverse D.chartScale hstate

/-- The algebraic first variation of the quadratic mode energy vanishes
on the declared lift.  Re-export of
`photonModeEnergy_firstVariation_generator_zero` at the bundled chart
scale. -/
theorem mode_energy_first_variation_zero (k : FreePhotonVec3)
    (state : PhotonModeState) :
    photonModeEnergyFirstVariation D.chartScale k state = 0 :=
  photonModeEnergy_firstVariation_generator_zero D.chartScale k state

/-- The spatial action commutes with the basis-free transverse
projector at the bundled chart scale.  Re-export of
`photonSpatialAction_commutes_projector`. -/
theorem spatialAction_commutes_transverse_projector
    {k : FreePhotonVec3} (hk : k ≠ 0) (v : FreePhotonVec3) :
    transverseProjector k (photonSpatialAction D.chartScale k v) =
      photonSpatialAction D.chartScale k (transverseProjector k v) :=
  photonSpatialAction_commutes_projector D.chartScale hk v

/-! ## The identified physical frequency (OL-F2 input) -/

/-- **Composition.**  The supplied readout frequency is the committed
nonnegative branch: both are nonnegative and square to the same dilated
carrier symbol. -/
theorem physicalFrequency_eq_mode_branch (k : FreePhotonVec3) :
    D.physicalFrequency k = photonModeFrequency D.chartScale k := by
  have hsq : D.physicalFrequency k ^ 2 =
      photonModeFrequency D.chartScale k ^ 2 := by
    rw [D.physicalFrequency_sq_eq_symbol k, photonModeFrequency_sq]
  calc D.physicalFrequency k
      = Real.sqrt (D.physicalFrequency k ^ 2) :=
        (Real.sqrt_sq (D.physicalFrequency_nonneg k)).symm
    _ = Real.sqrt (photonModeFrequency D.chartScale k ^ 2) := by
        rw [hsq]
    _ = photonModeFrequency D.chartScale k :=
        Real.sqrt_sq (photonModeFrequency_nonnegative D.chartScale k)

/-- **The exact zero mode.**  The identified frequency vanishes at zero
chart momentum.  A physical massless-photon statement additionally
consumes the open sector and clock attachments. -/
theorem physicalFrequency_zero_mode : D.physicalFrequency 0 = 0 := by
  rw [D.physicalFrequency_eq_mode_branch 0, photonModeFrequency_zero]

/-- The dilated carrier symbol has no constant mode at the bundled
chart scale.  Re-export of `dilatedCompletionFourierSymbol_zero`. -/
theorem symbol_zero_mode :
    dilatedCompletionFourierSymbol D.chartScale 0 = 0 :=
  dilatedCompletionFourierSymbol_zero D.chartScale

/-- **Composition.**  The supplied readout frequency is the exact
edge-30 feature norm: both are nonnegative and square to the same
edge-current character at the bundled chart scale. -/
theorem physicalFrequency_eq_edge30_frequency (k : FreePhotonVec3) :
    D.physicalFrequency k = frequency edge30Frame D.chartScale k := by
  have ha : D.chartScale ≠ 0 := ne_of_gt D.chartScale_pos
  have hsq : D.physicalFrequency k ^ 2 =
      frequency edge30Frame D.chartScale k ^ 2 := by
    rw [D.physicalFrequency_sq_eq_symbol k,
      dilatedCompletionFourierSymbol_eq_edgeCurrentCharacterSymbol
        D.chartScale ha k,
      frequency_sq_eq_cosineSymbol edge30Frame ha k,
      edge30_cosineSymbol_eq_fz12 ha k]
  calc D.physicalFrequency k
      = Real.sqrt (D.physicalFrequency k ^ 2) :=
        (Real.sqrt_sq (D.physicalFrequency_nonneg k)).symm
    _ = Real.sqrt (frequency edge30Frame D.chartScale k ^ 2) := by
        rw [hsq]
    _ = frequency edge30Frame D.chartScale k :=
        Real.sqrt_sq (frequency_nonnegative edge30Frame D.chartScale k)

/-- **Composition: the exact unit speed bound.**  The identified
frequency is globally one-Lipschitz in the chart momentum.  Calling the
bound a signal-speed bound consumes the separate position, clock,
wave-packet, frame, and readout attachments. -/
theorem physicalFrequency_global_one_lipschitz (k p : FreePhotonVec3) :
    |D.physicalFrequency k - D.physicalFrequency p| ≤
      euclideanMagnitude (k - p) := by
  rw [D.physicalFrequency_eq_edge30_frequency k,
    D.physicalFrequency_eq_edge30_frequency p]
  exact fz12_frequency_global_one_lipschitz
    (ne_of_gt D.chartScale_pos) k p

/-- **Composition: the wave-shaped evolution.**  Under the bundled
premises the supplied mode dynamics obeys `A'' + omega^2 A = 0` with
the identified physical frequency, mode by mode in the chart. -/
theorem transverse_wave_second_order (k : FreePhotonVec3)
    (state : PhotonModeState) :
    (D.modeDynamics k (D.modeDynamics k state)).1 +
        D.physicalFrequency k ^ 2 • state.1 = 0 := by
  have h := D.modeDynamics_second_order k state
  rw [photonSpatialAction_eq_frequency_sq] at h
  rw [D.physicalFrequency_eq_mode_branch k]
  exact h

end LightSignalPremiseData

/-! ## Side-by-side committed statements

The statements below consume none of the bundled premises.  They are
committed exact results of the imported modules, re-exported here so
the composed surface is one citable package.  The scalar lift equation
and its dispersion sit on the response-Gram completion; the generator
identities are premise-free carrier mathematics; the Gauss receipts are
exact rational statements about the seam incidence operator. -/

/-- The scalar auxiliary lift satisfies `q'' + L q = 0` exactly, with
`q''` the defined algebraic component.  Re-export of
`auxiliary_q_second_order`. -/
theorem scalar_lift_second_order (state : AuxiliaryPhaseState)
    (x : CarrierPoint) :
    auxiliaryQSecond state x +
        completionComplexDirichletGenerator state.1 x = 0 :=
  auxiliary_q_second_order state x

/-- The scalar lift plane-wave dispersion is exactly the canonical
completion symbol.  Re-export of `auxiliaryModeFrequency_sq`. -/
theorem scalar_lift_planeWave_dispersion (k : CarrierPoint) :
    auxiliaryModeFrequency k ^ 2 =
      OPH.SeamCurrentDirichletGenerator.completionFourierSymbol k :=
  auxiliaryModeFrequency_sq k

/-- The exact carre-du-champ identity of the dimensionless Dirichlet
generator.  Re-export of `completionDirichletGenerator_carre_du_champ`. -/
theorem generator_carre_du_champ (f : CarrierPoint → ℝ)
    (x : CarrierPoint) :
    2 * f x *
          OPH.SeamCurrentDirichletGenerator.completionDirichletGenerator
            f x -
        OPH.SeamCurrentDirichletGenerator.completionDirichletGenerator
          (fun y ↦ f y ^ 2) x =
      2 * OPH.SeamCurrentDirichletGenerator.completionCarreDuChamp f x :=
  OPH.SeamCurrentDirichletGenerator.completionDirichletGenerator_carre_du_champ
    f x

/-- Plane waves diagonalize the complex completion generator with the
exact nonnegative symbol.  Re-export of
`completionComplexDirichletGenerator_planeWave`. -/
theorem generator_planeWave (k x : CarrierPoint) :
    OPH.SeamCurrentDirichletGenerator.completionComplexDirichletGenerator
        (OPH.PrimitivePortTranslationBridge.planeWave k) x =
      (OPH.SeamCurrentDirichletGenerator.completionFourierSymbol k : ℂ) *
        OPH.PrimitivePortTranslationBridge.planeWave k x :=
  OPH.SeamCurrentDirichletGenerator.completionComplexDirichletGenerator_planeWave
    k x

/-- The normalized generator symbol is the edge-current character.
Re-export of
`completionFourierSymbol_normalizes_to_edgeCurrentCharacterSymbol`. -/
theorem generator_symbol_eq_edge_current_character (k : CarrierPoint) :
    (6 / OPH.SeamCurrentDirichletGenerator.completionStepScale ^ 2) *
        OPH.SeamCurrentDirichletGenerator.completionFourierSymbol k =
      OPH.SeamCurrentHomogeneousAction.edgeCurrentCharacterSymbol
        OPH.SeamCurrentDirichletGenerator.completionStepScale k :=
  OPH.SeamCurrentDirichletGenerator.completionFourierSymbol_normalizes_to_edgeCurrentCharacterSymbol
    k

/-- **Rank two at nonzero momentum.**  The basis-free transverse space
has exactly two real dimensions at every nonzero chart momentum.
Re-export of `transversePolarization_finrank`. -/
theorem transverse_rank_two
    {k : OPH.SeamCurrentFreePhotonLift.FreePhotonVec3} (hk : k ≠ 0) :
    Module.finrank ℝ
      (OPH.SeamCurrentFreePhotonLift.TransversePolarization k) = 2 :=
  OPH.SeamCurrentFreePhotonLift.transversePolarization_finrank hk

/-- **Discrete Gauss solvability.**  A rational port load admits a
seam-flux solution exactly when its total load vanishes.  Re-export of
`gauss_solution_exists_iff_total_zero`. -/
theorem gauss_solvability_iff_neutral
    (charge : OPH.SeamCurrentCarrierQuotient.RationalPortLoad) :
    (∃ field : OPH.SeamCurrentCarrierQuotient.RationalSeamCurrent,
        OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary field =
          charge) ↔
      OPH.SeamCurrentCarrierQuotient.rationalPortTotal charge = 0 :=
  OPH.DiscreteGauss.gauss_solution_exists_iff_total_zero charge

/-- The declared spanning-tree section solves the discrete Gauss
equation on every neutral load.  Re-export of
`rationalBoundarySection_is_gauss_solution`. -/
theorem gauss_declared_section_solves
    (charge : OPH.SeamCurrentCarrierQuotient.RationalPortLoad)
    (hneutral :
      OPH.SeamCurrentCarrierQuotient.rationalPortTotal charge = 0) :
    OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary
        (OPH.SeamCurrentCarrierQuotient.rationalBoundarySection charge) =
      charge :=
  OPH.DiscreteGauss.rationalBoundarySection_is_gauss_solution
    charge hneutral

/-- The full Gauss solution fibre over a fixed load is one affine
translate of the cycle kernel.  Re-export of
`gauss_solution_iff_difference_is_cycle`. -/
theorem gauss_fibre_is_cycle_translate
    {field field' : OPH.SeamCurrentCarrierQuotient.RationalSeamCurrent}
    {charge : OPH.SeamCurrentCarrierQuotient.RationalPortLoad}
    (hfield :
      OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary field =
        charge) :
    OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary field' =
        charge ↔
      field' - field ∈
        LinearMap.ker
          OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary :=
  OPH.DiscreteGauss.gauss_solution_iff_difference_is_cycle hfield

/-- **The committed source-free solution ambiguity.**  The cycle kernel of
the discrete Gauss receipts has exact rational dimension nineteen.
Re-export of `gauss_cycle_space_finrank`. -/
theorem gauss_cycle_rank_nineteen :
    Module.finrank ℚ
      (LinearMap.ker
        OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary) = 19 :=
  OPH.DiscreteGauss.gauss_cycle_space_finrank

/-! ## The OL-F1 boundary as a definition -/

/-- **The equations of motion attained at the committed order.**  On
the transverse sector the supplied dynamics preserves transversality
and its second iterate obeys the wave equation with the identified
physical frequency, mode by mode in the chart.  This definition is the
complete equation-of-motion content of the composed surface.
Maxwell-shaped field equations contain, in addition, the coupled
two-field curl pairing, the position-space field assembly with a gauge
potential and its quotient, source-coupled inhomogeneous equations,
boost covariance, and continuum control; none of those statements has a
committed producer, and this module states none of them as a theorem. -/
def attainedCommittedOrderLaw (D : LightSignalPremiseData) : Prop :=
  ∀ (k : OPH.SeamCurrentFreePhotonLift.FreePhotonVec3)
    (state : OPH.SeamCurrentFreePhotonLift.PhotonModeState),
    OPH.SeamCurrentFreePhotonLift.ModeIsTransverse k state →
      OPH.SeamCurrentFreePhotonLift.ModeIsTransverse k
          (D.modeDynamics k state) ∧
        (D.modeDynamics k (D.modeDynamics k state)).1 +
          D.physicalFrequency k ^ 2 • state.1 = 0

/-- The attained law holds for every premise bundle. -/
theorem attainedCommittedOrderLaw_holds (D : LightSignalPremiseData) :
    attainedCommittedOrderLaw D :=
  fun k state hstate ↦
    ⟨D.modeDynamics_preserves_transverse hstate,
      D.transverse_wave_second_order k state⟩

/-! ## Composite receipt -/

/-- **The light-signal adequacy surface receipt (issue #733).**  For
every premise bundle: the selected generator is the dimensionless
Dirichlet generator; the supplied mode dynamics obeys the second-order
wave equation with the identified frequency and preserves
transversality; the transverse space has rank two at every nonzero
momentum; the identified frequency has an exact zero mode and is
globally one-Lipschitz in the chart momentum; the discrete Gauss
equation is solvable exactly on neutral loads; and the source-free
solution ambiguity has exact rank nineteen.  The conjunction composes the
committed statements at the bundled premises and claims nothing beyond
them. -/
theorem lightSignalAdequacySurface_receipt (D : LightSignalPremiseData) :
    (selectedCompletionGenerator D.sourceCountingSelection =
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
      ∧ Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) = 19 :=
  ⟨D.selected_generator_eq_dirichlet,
    fun k state ↦ D.transverse_wave_second_order k state,
    fun _ _ hstate ↦ D.modeDynamics_preserves_transverse hstate,
    fun _ hk ↦ transversePolarization_finrank hk,
    D.physicalFrequency_zero_mode,
    fun k p ↦ D.physicalFrequency_global_one_lipschitz k p,
    OPH.DiscreteGauss.gauss_solution_exists_iff_total_zero,
    OPH.DiscreteGauss.gauss_cycle_space_finrank⟩

end OPH.LightSignalAdequacySurface

/-! ## Axiom audit

Statements that touch the seam alphabet or the seam-built symbol
inherit the disclosed native-decision certificates of
`SeamCurrentHomogeneousAction`; their audit lines show the
corresponding native-decide axioms.  The scalar lift statements, the
normalized symbol identity, the transverse rank, and the Gauss receipts
show only `propext`, `Classical.choice`, and `Quot.sound`. -/

#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.source_counting_weight_forced
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.selected_average_eq_source_counting
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.selected_generator_eq_dirichlet
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.modeDynamics_second_order
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.modeDynamics_preserves_transverse
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.mode_energy_first_variation_zero
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.spatialAction_commutes_transverse_projector
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.physicalFrequency_eq_mode_branch
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.physicalFrequency_zero_mode
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.symbol_zero_mode
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.physicalFrequency_eq_edge30_frequency
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.physicalFrequency_global_one_lipschitz
#print axioms OPH.LightSignalAdequacySurface.LightSignalPremiseData.transverse_wave_second_order
#print axioms OPH.LightSignalAdequacySurface.scalar_lift_second_order
#print axioms OPH.LightSignalAdequacySurface.scalar_lift_planeWave_dispersion
#print axioms OPH.LightSignalAdequacySurface.generator_carre_du_champ
#print axioms OPH.LightSignalAdequacySurface.generator_planeWave
#print axioms OPH.LightSignalAdequacySurface.generator_symbol_eq_edge_current_character
#print axioms OPH.LightSignalAdequacySurface.transverse_rank_two
#print axioms OPH.LightSignalAdequacySurface.gauss_solvability_iff_neutral
#print axioms OPH.LightSignalAdequacySurface.gauss_declared_section_solves
#print axioms OPH.LightSignalAdequacySurface.gauss_fibre_is_cycle_translate
#print axioms OPH.LightSignalAdequacySurface.gauss_cycle_rank_nineteen
#print axioms OPH.LightSignalAdequacySurface.attainedCommittedOrderLaw_holds
#print axioms OPH.LightSignalAdequacySurface.lightSignalAdequacySurface_receipt
