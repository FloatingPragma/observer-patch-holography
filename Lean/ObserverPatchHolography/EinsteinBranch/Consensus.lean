import Mathlib

/-!
# Finite consensus and the typed Einstein-branch boundary

This module formalises the model-theoretic and compositional boundary of the
compact paper's corrected Einstein branch.  A bare tower contains only finite
presentations, physical quotients, mismatch/repair data, normal forms, protected
boundary data, and coarse maps.  Geometry and stress are deliberately absent.

The geometric extension counterexample proves that no predicate of that bare
language decides the Einstein equation.  `CommonReadoutTower` then records the
six typed readouts used by the conditional branch.  The types share the same
repaired quotient domain; arrow commutation and refinement naturality are
separate predicates, not fields that are true by construction.

No continuum or physical premise is asserted in this file.
-/

namespace OPH.EinsteinBranch

universe u

/-- A regulator tower in the language of bare finite consensus.

`Quot r` is represented by an explicit surjective quotient map rather than by
Lean's `Quotient`, so a concrete implementation may use its own canonical
finite quotient representation.  The kernel of `quotient` is the redundancy
relation `Gamma_r` of the paper. -/
structure BareConsensusTower (ι : Type u) [Preorder ι] where
  /-- Finite presentation states `Sigma_r`. -/
  State : ι → Type
  stateFinite : ∀ r, Fintype (State r)
  /-- Physical quotient states `Q_r`. -/
  Quot : ι → Type
  quotFinite : ∀ r, Fintype (Quot r)
  /-- Protected boundary/sector values. -/
  Boundary : ι → Type
  quotient : ∀ r, State r → Quot r
  quotient_surjective : ∀ r, Function.Surjective (quotient r)
  mismatch : ∀ r, Quot r → ℕ
  step : ∀ r, Quot r → Quot r → Prop
  normalForm : ∀ r, Quot r → Quot r
  consistent : ∀ r, Quot r → Prop
  boundary : ∀ r, Quot r → Boundary r
  consistent_iff : ∀ r x, consistent r x ↔ mismatch r x = 0
  step_descends : ∀ r x y, step r x y → mismatch r y < mismatch r x
  normalForm_consistent : ∀ r x, consistent r (normalForm r x)
  normalForm_idempotent : ∀ r x, normalForm r (normalForm r x) = normalForm r x
  boundary_normalForm : ∀ r x, boundary r (normalForm r x) = boundary r x
  /-- Physical coarse map `c_sr : Q_s -> Q_r`, for `r <= s`. -/
  coarse : ∀ {r s}, r ≤ s → Quot s → Quot r
  coarse_refl : ∀ r x, coarse (show r ≤ r from le_rfl) x = x
  coarse_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) x,
    coarse hrs (coarse hst x) = coarse (hrs.trans hst) x
  normalForm_natural : ∀ {r s} (hrs : r ≤ s) x,
    coarse hrs (normalForm s x) = normalForm r (coarse hrs x)

attribute [instance] BareConsensusTower.stateFinite BareConsensusTower.quotFinite

/-- Redundancy/gauge equivalence is equality after the physical quotient map. -/
def GaugeEquivalent {ι : Type u} [Preorder ι] (T : BareConsensusTower ι)
    (r : ι) (x y : T.State r) : Prop :=
  T.quotient r x = T.quotient r y

/-- The protected-boundary fibre is a single physical quotient class on
consistent states.  This is a genuine predicate on the tower data. -/
def BoundaryFiber {ι : Type u} [Preorder ι] (T : BareConsensusTower ι) : Prop :=
  ∀ r (x y : T.Quot r), T.consistent r x → T.consistent r y →
    T.boundary r x = T.boundary r y → x = y

/-- A boundary-fibre receipt is not tautological: it fails whenever two
distinct consistent quotient states share a protected boundary value. -/
theorem not_boundaryFiber_of_witness {ι : Type u} [Preorder ι]
    (T : BareConsensusTower ι) (r : ι) (x y : T.Quot r)
    (hx : T.consistent r x) (hy : T.consistent r y)
    (hB : T.boundary r x = T.boundary r y) (hne : x ≠ y) :
    ¬ BoundaryFiber T := by
  intro h
  exact hne (h r x y hx hy hB)

/-! ## A concrete non-degenerate bare tower -/

/-- Broken-overlap count on two Boolean patches. -/
def demoMismatch (x : Bool × Bool) : ℕ := if x.1 = x.2 then 0 else 1

/-- A constant regulator tower with a real broken overlap and copy repair.
It witnesses that the non-entailment theorem is not driven by an empty or
one-point consensus model. -/
def demoTower : BareConsensusTower ℕ where
  State := fun _ => Bool × Bool
  stateFinite := fun _ => inferInstance
  Quot := fun _ => Bool × Bool
  quotFinite := fun _ => inferInstance
  Boundary := fun _ => Bool
  quotient := fun _ => id
  quotient_surjective := fun _ => Function.surjective_id
  mismatch := fun _ => demoMismatch
  step := fun _ x y => demoMismatch x = 1 ∧ demoMismatch y = 0
  normalForm := fun _ x => (x.1, x.1)
  consistent := fun _ x => x.1 = x.2
  boundary := fun _ x => x.1
  consistent_iff := by
    intro r x
    simp only [demoMismatch]
    by_cases h : x.1 = x.2 <;> simp [h]
  step_descends := by
    intro r x y h
    omega
  normalForm_consistent := by
    intro r x
    rfl
  normalForm_idempotent := by
    intro r x
    rfl
  boundary_normalForm := by
    intro r x
    rfl
  coarse := fun _ => id
  coarse_refl := by
    intro r x
    rfl
  coarse_trans := by
    intro r s t hrs hst x
    rfl
  normalForm_natural := by
    intro r s hrs x
    rfl

theorem demoTower_nondegenerate :
    ∃ r, ∃ x y : demoTower.Quot r,
      demoTower.consistent r x ∧ ¬ demoTower.consistent r y := by
  refine ⟨0, (true, true), (true, false), rfl, ?_⟩
  intro h
  exact Bool.noConfusion h

theorem demoTower_boundaryFiber : BoundaryFiber demoTower := by
  intro r x y hx hy hB
  rcases x with ⟨x1, x2⟩
  rcases y with ⟨y1, y2⟩
  simp only [demoTower] at hx hy hB ⊢
  subst x2
  subst y2
  subst y1
  rfl

/-! ## Model-extension non-entailment -/

/-- Geometry/stress decoration of a bare tower.  None of these fields occurs
in `BareConsensusTower`.  Integer scalar components suffice for the logical
separation; the result is about definability, not continuum geometry. -/
structure GeometricExtension where
  reduct : BareConsensusTower ℕ
  Point : Type
  point_nonempty : Nonempty Point
  metric : Point → ℤ
  curvature : Point → ℤ
  stress : Point → ℤ
  cosmological : ℤ
  coupling : ℤ

/-- Componentwise Einstein-equation shape on an extension. -/
def EinsteinEq (E : GeometricExtension) : Prop :=
  ∀ p : E.Point,
    E.curvature p + E.cosmological * E.metric p = E.coupling * E.stress p

def demoEinsteinExtension : GeometricExtension where
  reduct := demoTower
  Point := Unit
  point_nonempty := ⟨()⟩
  metric := fun _ => 1
  curvature := fun _ => 0
  stress := fun _ => 0
  cosmological := 0
  coupling := 1

def demoNonEinsteinExtension : GeometricExtension where
  reduct := demoTower
  Point := Unit
  point_nonempty := ⟨()⟩
  metric := fun _ => 1
  curvature := fun _ => 1
  stress := fun _ => 0
  cosmological := 0
  coupling := 1

theorem counterextensions_share_reduct :
    demoEinsteinExtension.reduct = demoNonEinsteinExtension.reduct := rfl

theorem einsteinEq_demoEinsteinExtension : EinsteinEq demoEinsteinExtension := by
  intro p
  simp [demoEinsteinExtension]

theorem not_einsteinEq_demoNonEinsteinExtension :
    ¬ EinsteinEq demoNonEinsteinExtension := by
  intro h
  have hp := h ()
  change (1 : ℤ) = 0 at hp
  omega

/-- Bare finite consensus is not Einstein-complete: no predicate of the bare
tower agrees with `EinsteinEq` on every geometric extension. -/
theorem bare_consensus_not_einstein_complete :
    ¬ ∃ decide : BareConsensusTower ℕ → Prop,
        ∀ E : GeometricExtension, (decide E.reduct ↔ EinsteinEq E) := by
  rintro ⟨decide, hdecide⟩
  have htrue : decide demoTower :=
    (hdecide demoEinsteinExtension).2 einsteinEq_demoEinsteinExtension
  exact not_einsteinEq_demoNonEinsteinExtension
    ((hdecide demoNonEinsteinExtension).1 htrue)

/-! ## Typed common-domain readouts -/

/-- Raw typed readouts from one repaired quotient domain.  The arrows and
refinement maps are data; their commuting laws are separate receipt
predicates below. -/
structure CommonReadoutTower {ι : Type u} [Preorder ι]
    (T : BareConsensusTower ι) where
  Geometry : ι → Type
  Modular : ι → Type
  Event : ι → Type
  Stress : ι → Type
  Entropy : ι → Type
  Scale : ι → Type
  geometryRead : ∀ r, T.Quot r → Geometry r
  modularRead : ∀ r, T.Quot r → Modular r
  eventRead : ∀ r, T.Quot r → Event r
  stressRead : ∀ r, T.Quot r → Stress r
  entropyRead : ∀ r, T.Quot r → Entropy r
  scaleRead : ∀ r, T.Quot r → Scale r
  geometryToModular : ∀ r, Geometry r → Modular r
  modularToEvent : ∀ r, Modular r → Event r
  eventToStress : ∀ r, Event r → Stress r
  modularToEntropy : ∀ r, Modular r → Entropy r
  entropyToGeometry : ∀ r, Entropy r → Geometry r
  geometryToScale : ∀ r, Geometry r → Scale r
  geometryCoarse : ∀ {r s}, r ≤ s → Geometry s → Geometry r
  modularCoarse : ∀ {r s}, r ≤ s → Modular s → Modular r
  eventCoarse : ∀ {r s}, r ≤ s → Event s → Event r
  stressCoarse : ∀ {r s}, r ≤ s → Stress s → Stress r
  entropyCoarse : ∀ {r s}, r ≤ s → Entropy s → Entropy r
  scaleCoarse : ∀ {r s}, r ≤ s → Scale s → Scale r

/-- Every downstream object is connected by a typed arrow on the same
quotient domain.  This is EB2's common-domain condition, as a failable
predicate rather than a constructor field. -/
def TypedArrowCommutation {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T) : Prop :=
  (∀ r q, R.geometryToModular r (R.geometryRead r q) = R.modularRead r q) ∧
  (∀ r q, R.modularToEvent r (R.modularRead r q) = R.eventRead r q) ∧
  (∀ r q, R.eventToStress r (R.eventRead r q) = R.stressRead r q) ∧
  (∀ r q, R.modularToEntropy r (R.modularRead r q) = R.entropyRead r q) ∧
  (∀ r q, R.entropyToGeometry r (R.entropyRead r q) = R.geometryRead r q) ∧
  (∀ r q, R.geometryToScale r (R.geometryRead r q) = R.scaleRead r q)

/-- All six readouts commute with physical refinement on the declared tower. -/
def ReadoutNaturality {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T) : Prop :=
  (∀ {r s} (h : r ≤ s) q,
    R.geometryCoarse h (R.geometryRead s q) = R.geometryRead r (T.coarse h q)) ∧
  (∀ {r s} (h : r ≤ s) q,
    R.modularCoarse h (R.modularRead s q) = R.modularRead r (T.coarse h q)) ∧
  (∀ {r s} (h : r ≤ s) q,
    R.eventCoarse h (R.eventRead s q) = R.eventRead r (T.coarse h q)) ∧
  (∀ {r s} (h : r ≤ s) q,
    R.stressCoarse h (R.stressRead s q) = R.stressRead r (T.coarse h q)) ∧
  (∀ {r s} (h : r ≤ s) q,
    R.entropyCoarse h (R.entropyRead s q) = R.entropyRead r (T.coarse h q)) ∧
  (∀ {r s} (h : r ≤ s) q,
    R.scaleCoarse h (R.scaleRead s q) = R.scaleRead r (T.coarse h q))

/-- The finite/common-domain part of an Einstein-admissible tower.  Its
existence is deliberately not proved: constructing a value requires an
actual bare consensus tower plus boundary-fibre, typed-arrow, and refinement
receipts.  Analytic and physical premises remain arguments of the composed
theorem in `Composition.lean`. -/
structure EinsteinAdmissibleTower (ι : Type u) [Preorder ι] where
  consensus : BareConsensusTower ι
  readouts : CommonReadoutTower consensus
  boundaryFiber : BoundaryFiber consensus
  typedArrows : TypedArrowCommutation readouts
  refinementNaturality : ReadoutNaturality readouts

/-- Equality of every readout emitted from a common quotient state. -/
structure ReadoutAgreement {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T)
    (r : ι) (q q' : T.Quot r) : Prop where
  geometry : R.geometryRead r q = R.geometryRead r q'
  modular : R.modularRead r q = R.modularRead r q'
  event : R.eventRead r q = R.eventRead r q'
  stress : R.stressRead r q = R.stressRead r q'
  entropy : R.entropyRead r q = R.entropyRead r q'
  scale : R.scaleRead r q = R.scaleRead r q'

/-- Boundary-fibre uniqueness composes with the common-domain readouts: two
consistent states with the same protected observation emit identical
geometry, modular, event, stress, entropy, and scale data. -/
theorem boundary_fiber_readout_composition {ι : Type u} [Preorder ι]
    {T : BareConsensusTower ι} (R : CommonReadoutTower T)
    (hfiber : BoundaryFiber T) (r : ι) (q q' : T.Quot r)
    (hq : T.consistent r q) (hq' : T.consistent r q')
    (hB : T.boundary r q = T.boundary r q') :
    ReadoutAgreement R r (T.normalForm r q) (T.normalForm r q') := by
  have hqq' : q = q' := hfiber r q q' hq hq' hB
  subst q'
  exact ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩

/-! ## Per-theorem axiom audit -/

#print axioms not_boundaryFiber_of_witness
#print axioms demoTower_nondegenerate
#print axioms demoTower_boundaryFiber
#print axioms counterextensions_share_reduct
#print axioms einsteinEq_demoEinsteinExtension
#print axioms not_einsteinEq_demoNonEinsteinExtension
#print axioms bare_consensus_not_einstein_complete
#print axioms boundary_fiber_readout_composition

end OPH.EinsteinBranch
