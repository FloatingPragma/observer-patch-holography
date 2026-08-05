import Geometry.EventGermDisplacement

/-!
# C2: celestial soldering of future-null event displacements

The C1 celestial sphere is a quotient of intrinsic future-null vectors by
positive scale.  Here an oriented Lorentz overlap map acts on that quotient,
and a future-null event displacement is soldered to a celestial direction.
The soldered direction is exactly overlap covariant.

Future-nullness of a record displacement remains an input.  No theorem here
turns arbitrary records into null signals or identifies the celestial value
with a physical sky measurement.
-/

namespace OPH.C2Soldering

noncomputable section

open OPH.C1Lorentz

theorem OrientedLorentzEquiv.mapFutureNullVector_rel
    (L : OrientedLorentzEquiv) {u v : FutureNullVector}
    (h : FutureNullRayRel u v) :
    FutureNullRayRel (L.mapFutureNullVector u) (L.mapFutureNullVector v) := by
  rcases h with ⟨a, ha, hva⟩
  refine ⟨a, ha, ?_⟩
  change L (v : Herm2) = a • L (u : Herm2)
  rw [hva, L.map_smul]

/-- Induced action of an oriented Lorentz map on future-null rays. -/
def OrientedLorentzEquiv.mapFutureNullRay (L : OrientedLorentzEquiv) :
    FutureNullRay → FutureNullRay :=
  Quotient.map L.mapFutureNullVector
    (fun _ _ h => L.mapFutureNullVector_rel h)

@[simp] theorem OrientedLorentzEquiv.mapFutureNullRay_refl
    (r : FutureNullRay) :
    OrientedLorentzEquiv.refl.mapFutureNullRay r = r := by
  refine Quotient.inductionOn r ?_
  intro v
  rfl

theorem OrientedLorentzEquiv.mapFutureNullRay_trans
    (L M : OrientedLorentzEquiv) (r : FutureNullRay) :
    (L.trans M).mapFutureNullRay r =
      M.mapFutureNullRay (L.mapFutureNullRay r) := by
  refine Quotient.inductionOn r ?_
  intro v
  rfl

/-- Induced conformal action on the set-level C1 celestial sphere. -/
def OrientedLorentzEquiv.celestialAction (L : OrientedLorentzEquiv)
    (n : CelestialSphere) : CelestialSphere :=
  rayToCelestial (L.mapFutureNullRay (celestialToRay n))

@[simp] theorem OrientedLorentzEquiv.celestialAction_refl
    (n : CelestialSphere) :
    OrientedLorentzEquiv.refl.celestialAction n = n := by
  simp [OrientedLorentzEquiv.celestialAction,
    rayToCelestial_celestialToRay]

theorem OrientedLorentzEquiv.celestialAction_trans
    (L M : OrientedLorentzEquiv) (n : CelestialSphere) :
    (L.trans M).celestialAction n =
      M.celestialAction (L.celestialAction n) := by
  unfold OrientedLorentzEquiv.celestialAction
  rw [L.mapFutureNullRay_trans M]
  rw [celestialToRay_rayToCelestial]

/-- A future-null coordinate displacement as a positive ray. -/
def EventGermAtlas.futureNullDisplacementRay {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p q : Event)
    (hnull : IsFutureNull (A.displacement i p q)) : FutureNullRay :=
  Quotient.mk futureNullRaySetoid ⟨A.displacement i p q, hnull⟩

/-- Celestial direction soldered to a future-null event displacement. -/
def EventGermAtlas.celestialSolder {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p q : Event)
    (hnull : IsFutureNull (A.displacement i p q)) : CelestialSphere :=
  rayToCelestial (A.futureNullDisplacementRay i p q hnull)

theorem EventGermAtlas.futureNullDisplacementRay_overlap
    {Event Chart : Type*} (A : EventGermAtlas Event Chart)
    {i j : Chart} {p q : Event}
    (hip : A.visible i p) (hiq : A.visible i q)
    (hjp : A.visible j p) (hjq : A.visible j q)
    (hi : IsFutureNull (A.displacement i p q))
    (hj : IsFutureNull (A.displacement j p q)) :
    (A.overlap.lorentz i j).mapFutureNullRay
        (A.futureNullDisplacementRay i p q hi) =
      A.futureNullDisplacementRay j p q hj := by
  apply Quotient.sound
  refine ⟨1, by norm_num, ?_⟩
  change A.displacement j p q =
    (1 : ℝ) • A.overlap.lorentz i j (A.displacement i p q)
  simpa using A.displacement_overlap hip hiq hjp hjq

/-- Exact celestial overlap covariance. -/
theorem EventGermAtlas.celestialSolder_overlap
    {Event Chart : Type*} (A : EventGermAtlas Event Chart)
    {i j : Chart} {p q : Event}
    (hip : A.visible i p) (hiq : A.visible i q)
    (hjp : A.visible j p) (hjq : A.visible j q)
    (hi : IsFutureNull (A.displacement i p q))
    (hj : IsFutureNull (A.displacement j p q)) :
    A.celestialSolder j p q hj =
      (A.overlap.lorentz i j).celestialAction
        (A.celestialSolder i p q hi) := by
  unfold EventGermAtlas.celestialSolder
    OrientedLorentzEquiv.celestialAction
  rw [celestialToRay_rayToCelestial]
  rw [A.futureNullDisplacementRay_overlap hip hiq hjp hjq hi hj]

#print axioms EventGermAtlas.futureNullDisplacementRay_overlap
#print axioms EventGermAtlas.celestialSolder_overlap
#print axioms OrientedLorentzEquiv.celestialAction_trans

end

end OPH.C2Soldering
