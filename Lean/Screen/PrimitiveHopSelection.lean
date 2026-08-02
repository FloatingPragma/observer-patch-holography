import PrimitivePortOperatorSelectionBoundary
import ObserverPatchHolography.EqualSeamSelection

namespace OPH.PrimitiveHopSelection

open OPH.PrimitivePortTranslationBridge
open OPH.PrimitivePortOperatorSelectionBoundary
open ObserverPatchHolography.EqualSeamSelection

/-!
# Conditional selection of the primitive-hop generator

This module states the strongest finite selection implication supported by the
current #655/#663 theorem stack.  The primitive move type is assumed to be the
complete quotient-deduplicated alphabet, with exactly one event identity for
each of the twelve ports.  A natural unique A3 move projection fixes its
weights.  A separate exact realization premise says that one candidate uses
all and only those atomic direct hops, with no additional range or operator
term.  Antipodal realization derives reciprocity, and a supplied unit-quadratic
normalization fixes the remaining overall coefficient.  The candidate is then
the frozen equal-weight primitive-port operator.

The proposed repair clauses would discharge only part of this packet:

* A1-R clauses 8, 9, and 11 can construct the complete atomic event type,
  remove presentation duplicates, and certify its one-to-one port census.
* A2-R presentation naturality together with A3 uniqueness can construct the
  `NaturalUniqueMoveProjection` on a single transitive move orbit.  The
  source-counting alternative in `RepairWordSchedule` gives the same one-event
  law when its stronger full-simplex premises hold.
* A1-R grammar completeness and A2-R completed reconciliation can exclude
  undeclared primitive repair maps on their common repair workspace.  They do
  not by themselves identify that repair generator with the declared spatial
  translation operator; the direct-only equality below remains an explicit
  source-to-translation attachment.
* Antipodal spatial realization comes from the finite carrier geometry plus
  the conditional translation adapter, not from A1-R or A2-R alone.
* Unit quadratic normalization is a separate long-wavelength normalization;
  neither proposed repair clause currently fixes it.

No theorem here adopts A1-R or A2-R, proves their premises from A1--A3,
identifies a physical field, chooses a time equation or scale, or excludes
operators outside the explicitly supplied direct-only realization.
-/

universe u v

/-- A complete quotient-visible primitive move census with exactly one move
    identity per port.  The equivalence is data supplied by the proposed A1-R
    atomicity, completeness, and deduplication audit; this structure does not
    derive those clauses. -/
structure CompletePrimitivePortGrammar (Move : Type u) where
  moveToPort : Move ≃ Fin 12

/-- Transport a move-law weight to its unique port identity. -/
def portWeightFromMoveLaw {Move : Type u}
    (grammar : CompletePrimitivePortGrammar Move)
    (moveWeight : Move → ℝ) : Fin 12 → ℝ :=
  fun i => moveWeight (grammar.moveToPort.symm i)

/-- Direct-only finite-range generator before its overall normalization is
    fixed.  Its exact form is the load-bearing `no additional ranges or terms`
    premise in the composed theorem. -/
noncomputable def directOnlyGenerator
    (normalization a : ℝ) (u : Fin 12 → Vec3) (w : Fin 12 → ℝ)
    (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x => (normalization : ℂ) *
    ∑ i : Fin 12, (w i : ℂ) * directedDifference (a • u i) f x

/-- The supplied frozen quadratic normalization.  Nonzero scale prevents the
    field convention `0⁻¹ = 0` from masquerading as a continuum
    normalization.  The factor `6/a²` is the state-weighted form of the
    primitive operator's `1/(2a²)` prefactor after the twelve weights become
    `1/12`. -/
def FixedQuadraticNormalization (normalization a : ℝ) : Prop :=
  a ≠ 0 ∧ normalization = 6 / a ^ 2

/-- The unit-quadratic normalization used by the frozen operator converts the
    generic direct-only generator to the existing state-weighted direct-hop
    operator.  The equality is supplied, not derived from A1-R or A2-R. -/
theorem directOnlyGenerator_of_fixed_quadratic_normalization
    (normalization a : ℝ) (u : Fin 12 → Vec3) (w : Fin 12 → ℝ)
    (hquadratic : FixedQuadraticNormalization normalization a) :
    directOnlyGenerator normalization a u w =
      stateWeightedDirectedFirstHopOperator a u w := by
  rw [hquadratic.2]
  funext f x
  unfold directOnlyGenerator stateWeightedDirectedFirstHopOperator
  push_cast
  rfl

/-- A natural unique A3 law on a complete twelve-move grammar gives weight
    `1/12` to the unique move attached to every port. -/
theorem selected_port_weight_is_one_twelfth
    {Move : Type u} {Presentation : Type v}
    [Fintype Move] [DecidableEq Move] [Nonempty Move]
    (grammar : CompletePrimitivePortGrammar Move)
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    ∀ i : Fin 12,
      portWeightFromMoveLaw grammar selection.selected i = 1 / 12 := by
  intro i
  rw [portWeightFromMoveLaw, selected_eq_inverse_card selection]
  have hcard : Fintype.card Move = 12 := by
    calc
      Fintype.card Move = Fintype.card (Fin 12) :=
        Fintype.card_congr grammar.moveToPort
      _ = 12 := Fintype.card_fin 12
  rw [hcard]
  norm_num

/-- Conditional uniqueness theorem for the primitive-hop lane.  The exact
    candidate equality `hdirectOnly` simultaneously says that the complete
    quotient-deduplicated move census is realized by one atomic direct hop per
    port and that no second range, potential, or other operator term is
    present.  Under that premise, A3 uniformity, antipodal inverse pairing,
    and fixed quadratic normalization leave the frozen primitive cosine
    operator as the unique candidate. -/
theorem complete_direct_only_a3_generator_eq_frozen
    {Move : Type u} {Presentation : Type v}
    [Fintype Move] [DecidableEq Move] [Nonempty Move]
    (grammar : CompletePrimitivePortGrammar Move)
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation))
    (normalization a : ℝ) (u : Fin 12 → Vec3)
    (hanti : AntipodalSpatialFrame u)
    (hquadratic : FixedQuadraticNormalization normalization a)
    (candidate : (Vec3 → ℂ) → Vec3 → ℂ)
    (hdirectOnly : candidate =
      directOnlyGenerator normalization a u
        (portWeightFromMoveLaw grammar selection.selected)) :
    candidate = primitivePortOperator a u := by
  have huniform : ∀ i : Fin 12,
      portWeightFromMoveLaw grammar selection.selected i = 1 / 12 :=
    selected_port_weight_is_one_twelfth grammar selection
  calc
    candidate = directOnlyGenerator normalization a u
        (portWeightFromMoveLaw grammar selection.selected) := hdirectOnly
    _ = stateWeightedDirectedFirstHopOperator a u
        (portWeightFromMoveLaw grammar selection.selected) :=
      directOnlyGenerator_of_fixed_quadratic_normalization
        normalization a u
          (portWeightFromMoveLaw grammar selection.selected) hquadratic
    _ = primitivePortOperator a u :=
      stateWeightedDirectedFirstHopOperator_of_uniform
        a u (portWeightFromMoveLaw grammar selection.selected) hanti huniform

/-- Any two candidates satisfying the same complete direct-only realization
    packet are equal.  This is uniqueness inside the explicitly stated
    grammar, not an exclusion theorem for unregistered physical operators. -/
theorem complete_direct_only_a3_generator_unique
    {Move : Type u} {Presentation : Type v}
    [Fintype Move] [DecidableEq Move] [Nonempty Move]
    (grammar : CompletePrimitivePortGrammar Move)
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation))
    (normalization a : ℝ) (u : Fin 12 → Vec3)
    (hanti : AntipodalSpatialFrame u)
    (hquadratic : FixedQuadraticNormalization normalization a)
    (candidate₁ candidate₂ : (Vec3 → ℂ) → Vec3 → ℂ)
    (hdirectOnly₁ : candidate₁ =
      directOnlyGenerator normalization a u
        (portWeightFromMoveLaw grammar selection.selected))
    (hdirectOnly₂ : candidate₂ =
      directOnlyGenerator normalization a u
        (portWeightFromMoveLaw grammar selection.selected)) :
    candidate₁ = candidate₂ := by
  rw [complete_direct_only_a3_generator_eq_frozen grammar selection
    normalization a u hanti hquadratic candidate₁ hdirectOnly₁]
  rw [complete_direct_only_a3_generator_eq_frozen grammar selection
    normalization a u hanti hquadratic candidate₂ hdirectOnly₂]

#print axioms directOnlyGenerator_of_fixed_quadratic_normalization
#print axioms selected_port_weight_is_one_twelfth
#print axioms complete_direct_only_a3_generator_eq_frozen
#print axioms complete_direct_only_a3_generator_unique

end OPH.PrimitiveHopSelection
