import Geometry.CanonicalLorentzModule

/-!
# C1: the celestial projective future null cone

This file constructs the actual quotient of the nonzero future null cone by
positive rescaling and proves a set-level equivalence with the unit
two-sphere.  It does not assert a homeomorphism or diffeomorphism and makes no
physical spacetime, clock, or continuum attachment.
-/

namespace OPH.C1Lorentz

noncomputable section

theorem spatialNormSq_smul (a : ℝ) (x : Spatial) :
    spatialNormSq (a • x) = a ^ 2 * spatialNormSq x := by
  simp only [spatialNormSq, Pi.smul_apply, smul_eq_mul, mul_pow]
  rw [Finset.mul_sum]

theorem spatialNormSq_pos {x : Spatial} (hx : x ≠ 0) : 0 < spatialNormSq x := by
  have h := spatial_axis_negative hx
  simpa [lorentzQ] using h

/-- Future nonzero null vectors in the declared `(+---)` Lorentz module. -/
def IsFutureNull (v : Herm2) : Prop := lorentzQ v = 0 ∧ 0 < v.1

abbrev FutureNullVector := {v : Herm2 // IsFutureNull v}

/-- Positive-ray equivalence on the future null cone. -/
def FutureNullRayRel (u v : FutureNullVector) : Prop :=
  ∃ a : ℝ, 0 < a ∧ (v : Herm2) = a • (u : Herm2)

theorem futureNullRayRel_refl (u : FutureNullVector) : FutureNullRayRel u u := by
  exact ⟨1, by norm_num, by simp⟩

theorem futureNullRayRel_symm {u v : FutureNullVector}
    (h : FutureNullRayRel u v) : FutureNullRayRel v u := by
  rcases h with ⟨a, ha, hva⟩
  refine ⟨a⁻¹, inv_pos.mpr ha, ?_⟩
  calc
    (u : Herm2) = a⁻¹ • (a • (u : Herm2)) := by simp [ne_of_gt ha]
    _ = a⁻¹ • (v : Herm2) := by rw [hva]

theorem futureNullRayRel_trans {u v w : FutureNullVector}
    (huv : FutureNullRayRel u v) (hvw : FutureNullRayRel v w) :
    FutureNullRayRel u w := by
  rcases huv with ⟨a, ha, hva⟩
  rcases hvw with ⟨b, hb, hwb⟩
  refine ⟨b * a, mul_pos hb ha, ?_⟩
  rw [hwb, hva, smul_smul]

/-- The actual positive-scaling quotient of the future null cone. -/
def futureNullRaySetoid : Setoid FutureNullVector where
  r := FutureNullRayRel
  iseqv := ⟨futureNullRayRel_refl, futureNullRayRel_symm,
    futureNullRayRel_trans⟩

abbrev FutureNullRay := Quotient futureNullRaySetoid

/-- The unit two-sphere in the three Pauli coordinates. -/
abbrev CelestialSphere := {x : Spatial // spatialNormSq x = 1}

/-- Normalized spatial direction of a future null vector. -/
def normalizedDirection (v : FutureNullVector) : CelestialSphere := by
  refine ⟨(v.1.1)⁻¹ • v.1.2, ?_⟩
  rw [spatialNormSq_smul]
  have hq := v.2.1
  have ht : v.1.1 ≠ 0 := ne_of_gt v.2.2
  simp only [lorentzQ] at hq
  field_simp
  nlinarith

theorem normalizedDirection_eq_of_ray {u v : FutureNullVector}
    (h : FutureNullRayRel u v) : normalizedDirection u = normalizedDirection v := by
  rcases h with ⟨a, ha, hva⟩
  apply Subtype.ext
  funext i
  have ha0 : a ≠ 0 := ne_of_gt ha
  have hut : u.1.1 ≠ 0 := ne_of_gt u.2.2
  change (u.1.1)⁻¹ * u.1.2 i = (v.1.1)⁻¹ * v.1.2 i
  rw [hva]
  simp only [Prod.smul_fst, Prod.smul_snd, Pi.smul_apply, smul_eq_mul]
  field_simp

/-- Canonical map from positive future-null rays to the celestial sphere. -/
def rayToCelestial : FutureNullRay → CelestialSphere :=
  Quotient.lift normalizedDirection
    (fun _ _ h ↦ normalizedDirection_eq_of_ray h)

/-- The future-null representative `(1,n)` of a celestial direction. -/
def celestialRepresentative (n : CelestialSphere) : FutureNullVector := by
  refine ⟨(1, n.1), ?_⟩
  exact ⟨by simp [lorentzQ, n.2], by norm_num⟩

/-- A celestial direction regarded as a positive future-null ray. -/
def celestialToRay (n : CelestialSphere) : FutureNullRay :=
  Quotient.mk futureNullRaySetoid (celestialRepresentative n)

theorem rayToCelestial_celestialToRay (n : CelestialSphere) :
    rayToCelestial (celestialToRay n) = n := by
  apply Subtype.ext
  funext i
  simp [rayToCelestial, celestialToRay, celestialRepresentative,
    normalizedDirection]

theorem celestialToRay_rayToCelestial (r : FutureNullRay) :
    celestialToRay (rayToCelestial r) = r := by
  refine Quotient.inductionOn r ?_
  intro v
  apply Quotient.sound
  refine ⟨v.1.1, v.2.2, ?_⟩
  apply Prod.ext
  · simp [rayToCelestial, celestialRepresentative, normalizedDirection]
  · funext i
    simp only [rayToCelestial, celestialRepresentative, normalizedDirection,
      Quotient.lift_mk, Subtype.coe_mk, Prod.smul_snd, Pi.smul_apply,
      smul_eq_mul]
    field_simp [ne_of_gt v.2.2]

/-- Exact set-level identification of the projectivized future null cone with
`S²`.  No topology or smooth structure is bundled in this equivalence. -/
def futureNullRayEquivCelestial : FutureNullRay ≃ CelestialSphere where
  toFun := rayToCelestial
  invFun := celestialToRay
  left_inv := celestialToRay_rayToCelestial
  right_inv := rayToCelestial_celestialToRay

#print axioms futureNullRayRel_symm
#print axioms normalizedDirection
#print axioms futureNullRayEquivCelestial

end


end OPH.C1Lorentz
