import Geometry.EventFrameSoldering
import PrimitivePortFrameQuotient

/-!
# C2: local spatial readback in an observer rest space

For a future-unit frame `u`, every intrinsic displacement has a canonical
Lorentz-orthogonal component

`v - lorentzB u v • u`.

This file proves that the component lies in the three-dimensional C1 rest
space, gives the exact temporal-plus-rest decomposition, and proves overlap
covariance under the C2 frame/displacement soldering contract.

The readback is frame-local.  It is not absolute space, a rod calibration, a
global spatial slice, or an attachment of the separate rank-three screen
carrier to physical position.
-/

namespace OPH.C2Soldering

noncomputable section

open OPH.C1Lorentz

/-- Lorentz-orthogonal component of an ambient displacement relative to one
future-unit observer frame. -/
def restProjection (u : FrameHyperboloid) (v : Herm2) : RestSpace u := by
  refine ⟨v - lorentzB u.1 v • u.1, ?_⟩
  change lorentzB u.1 (v - lorentzB u.1 v • u.1) = 0
  rw [sub_eq_add_neg, lorentzB_add_right]
  have hneg : lorentzB u.1 (-(lorentzB u.1 v • u.1)) =
      -lorentzB u.1 (lorentzB u.1 v • u.1) := by
    rw [show -(lorentzB u.1 v • u.1) =
      (-1 : ℝ) • (lorentzB u.1 v • u.1) by simp,
      lorentzB_smul_right]
    ring
  rw [hneg, lorentzB_smul_right, lorentzB_self, u.2.1]
  ring

@[simp] theorem restProjection_coe (u : FrameHyperboloid) (v : Herm2) :
    ((restProjection u v : RestSpace u) : Herm2) =
      v - lorentzB u.1 v • u.1 := rfl

/-- Exact temporal-plus-rest decomposition at one frame. -/
theorem restProjection_decomposition (u : FrameHyperboloid) (v : Herm2) :
    lorentzB u.1 v • u.1 + (restProjection u v : Herm2) = v := by
  simp [restProjection]

/-- Rest projection is natural under every oriented Lorentz overlap map. -/
theorem OrientedLorentzEquiv.restProjection_covariant
    (L : OrientedLorentzEquiv) (u : FrameHyperboloid) (v : Herm2) :
    L (restProjection u v : Herm2) =
      (restProjection (L.mapFrame u) (L v) : Herm2) := by
  change L (v - lorentzB u.1 v • u.1) =
    L v - lorentzB (L u.1) (L v) • L u.1
  rw [L.map_sub, L.map_smul, L.map_lorentzB]

/-- Linear transport between the rest fibers of Lorentz-related frames. -/
def OrientedLorentzEquiv.restLinearMap (L : OrientedLorentzEquiv)
    (u : FrameHyperboloid) : RestSpace u →ₗ[ℝ] RestSpace (L.mapFrame u) where
  toFun v := by
    refine ⟨L v.1, ?_⟩
    change lorentzB (L u.1) (L v.1) = 0
    rw [L.map_lorentzB]
    exact v.2
  map_add' v w := by
    apply Subtype.ext
    exact L.map_add v.1 w.1
  map_smul' a v := by
    apply Subtype.ext
    exact L.map_smul a v.1

theorem OrientedLorentzEquiv.restLinearMap_injective
    (L : OrientedLorentzEquiv) (u : FrameHyperboloid) :
    Function.Injective (L.restLinearMap u) := by
  intro v w h
  apply Subtype.ext
  apply L.toLinearEquiv.injective
  exact congrArg Subtype.val h

theorem OrientedLorentzEquiv.restLinearMap_surjective
    (L : OrientedLorentzEquiv) (u : FrameHyperboloid) :
    Function.Surjective (L.restLinearMap u) := by
  intro w
  let v : Herm2 := L.toLinearEquiv.symm w.1
  have hw := w.2
  change lorentzB (L u.1) w.1 = 0 at hw
  have hv : lorentzB u.1 v = 0 := by
    calc
      lorentzB u.1 v = lorentzB (L u.1) (L v) :=
        (L.map_lorentzB u.1 v).symm
      _ = lorentzB (L u.1) w.1 := by simp [v]
      _ = 0 := hw
  refine ⟨⟨v, hv⟩, ?_⟩
  apply Subtype.ext
  simp [OrientedLorentzEquiv.restLinearMap, v]

/-- Rest-space transport is a linear equivalence, not merely a coordinate
map. -/
def OrientedLorentzEquiv.restEquiv (L : OrientedLorentzEquiv)
    (u : FrameHyperboloid) : RestSpace u ≃ₗ[ℝ] RestSpace (L.mapFrame u) :=
  LinearEquiv.ofBijective (L.restLinearMap u)
    ⟨L.restLinearMap_injective u, L.restLinearMap_surjective u⟩

/-- Oriented Lorentz transport preserves the positive rest metric exactly. -/
theorem OrientedLorentzEquiv.restEquiv_preserves_metric
    (L : OrientedLorentzEquiv) (u : FrameHyperboloid)
    (v w : RestSpace u) :
    restMetric (L.mapFrame u) (L.restEquiv u v) (L.restEquiv u w) =
      restMetric u v w := by
  change -lorentzB (L v.1) (L w.1) = -lorentzB v.1 w.1
  rw [L.map_lorentzB]

/-! ## Candidate bridge from the source response quotient

The separate screen theorem gives an abstract rank-three Gram quotient.  The
following equivalence places it in the rest fiber of the *standard internal
frame*.  This is a type-safe candidate readback only; no source theorem
selects the standard frame or identifies the quotient with physical space.
-/

open OPH.PrimitivePortFrameQuotient

/-- Cartesian source-frame coordinates as the rest fiber of the standard C1
observer frame. -/
noncomputable def vec3EquivStandardRest :
    OPH.PrimitivePortTranslationBridge.Vec3 ≃ₗ[ℝ] RestSpace standardFrame where
  toFun x := by
    refine ⟨(0, x), ?_⟩
    change lorentzB (standardFrame : Herm2) (0, x) = 0
    simp [standardFrame, lorentzB, spatialDot]
  invFun v := v.1.2
  left_inv x := rfl
  right_inv v := by
    apply Subtype.ext
    apply Prod.ext
    · have hv := v.2
      change lorentzB (standardFrame : Herm2) v.1 = 0 at hv
      simpa [standardFrame, lorentzB, spatialDot] using hv.symm
    · rfl
  map_add' x y := by
    apply Subtype.ext
    simp
  map_smul' a x := by
    apply Subtype.ext
    simp

/-- Exact linear equivalence from the source Gram quotient to the candidate
standard-frame rest readout. -/
noncomputable def frameQuotientEquivStandardRest :
    FrameQuotient ≃ₗ[ℝ] RestSpace standardFrame :=
  quotientEquivVec3.trans vec3EquivStandardRest

/-- The source quotient Gram is exactly the standard-frame rest metric under
the candidate bridge. -/
theorem frameQuotientEquivStandardRest_preserves_metric
    (q r : FrameQuotient) :
    restMetric standardFrame (frameQuotientEquivStandardRest q)
        (frameQuotientEquivStandardRest r) =
      quotientGram q r := by
  simp [frameQuotientEquivStandardRest, vec3EquivStandardRest,
    restMetric, lorentzB, spatialDot, quotientGram,
    OPH.PrimitivePortTranslationBridge.dot]

/-- Frame-local spatial component of an event displacement. -/
def EventFrameSoldering.spatialReadback {Event Chart : Type*}
    (S : EventFrameSoldering Event Chart) (i : Chart) (p q : Event) :
    RestSpace (S.frame i p) :=
  restProjection (S.frame i p) (S.displacement i p q)

/-- Exact split of a chart displacement into one frame-time coefficient and
its local three-dimensional rest-space readback. -/
theorem EventFrameSoldering.displacement_time_add_spatial
    {Event Chart : Type*} (S : EventFrameSoldering Event Chart)
    (i : Chart) (p q : Event) :
    lorentzB (S.frame i p).1 (S.displacement i p q) • (S.frame i p).1 +
        (S.spatialReadback i p q : Herm2) =
      S.displacement i p q :=
  restProjection_decomposition (S.frame i p) (S.displacement i p q)

/-- Spatial readback is overlap covariant as an ambient vector.  The target
and source live in different dependent rest-space types, so the theorem
compares their canonical `Herm2` representatives. -/
theorem EventFrameSoldering.spatialReadback_overlap
    {Event Chart : Type*} (S : EventFrameSoldering Event Chart)
    {i j : Chart} {p q : Event}
    (hip : S.visible i p) (hiq : S.visible i q)
    (hjp : S.visible j p) (hjq : S.visible j q) :
    (S.spatialReadback j p q : Herm2) =
      S.overlap.lorentz i j (S.spatialReadback i p q : Herm2) := by
  have hf := S.frame_overlap hip hjp
  have hd := S.toEventGermAtlas.displacement_overlap hip hiq hjp hjq
  change (restProjection (S.frame j p) (S.displacement j p q) : Herm2) =
    S.overlap.lorentz i j
      (restProjection (S.frame i p) (S.displacement i p q) : Herm2)
  rw [hf, hd]
  exact (S.overlap.lorentz i j).restProjection_covariant
    (S.frame i p) (S.displacement i p q) |>.symm

#print axioms restProjection_decomposition
#print axioms OrientedLorentzEquiv.restProjection_covariant
#print axioms OrientedLorentzEquiv.restEquiv_preserves_metric
#print axioms frameQuotientEquivStandardRest_preserves_metric
#print axioms EventFrameSoldering.spatialReadback_overlap

end

end OPH.C2Soldering
