import Geometry.SpatialReadbackSoldering
import Geometry.InverseSquareShellLaw

/-!
# The inverse-square shell law transported to declared rest fibers

V3 gravitation lane (issue #729), Newtonian-limit row, rest-fiber transport
step.  The committed shell law lives on the completion carrier
`EuclideanVec3` (`Geometry/InverseSquareShellLaw.lean`), and the committed
rest-space theory gives every future unit-timelike frame of the declared
Hermitian Lorentz module a three-dimensional positive-definite rest space
(`Geometry/ObserverRestSpace.lean`), with the source Gram quotient
identified with the standard-frame rest fiber as a candidate readback
(`Geometry/SpatialReadbackSoldering.lean`).  This module joins the two
committed packages.

For every hyperboloid frame `u = (γ, w)` the map `restFiberIsometry u` is an
explicit linear equivalence from the completion carrier onto `RestSpace u`
whose pullback of the committed rest metric is exactly the Euclidean
coordinate form.  The formula is the boost shear
`x ↦ (⟨w, x⟩, x + (⟨w, x⟩ / (γ + 1)) • w)`; at the standard frame it reduces
to the committed candidate bridge `vec3EquivStandardRest`, so the committed
screen chain (finite source Gram table, integer record controls, completion
carrier) lands in the rest fiber through committed pieces plus this one
typed equivalence.  Under the declared radial readout pulled back through
the equivalence the committed shell-law theorems transport with no additional
proof: the strength at a nonzero rest vector is the charge over
`4 π r ^ 2` with `r` the rest radius, the pulled-back normalization and
scale-free identities hold, and the exponent reads on the module's own
dimension because `finrank ℝ (RestSpace u) = 3` holds through the transport
equivalence itself (agreeing with the committed `finrank_restSpace`).  The
committed oriented Lorentz rest transport (`OrientedLorentzEquiv.restEquiv`)
carries the whole package to the transported frame's rest space exactly.
One composed identity carries the screen chain into the fiber: the finite
source Gram of two integer record controls equals the rest metric of their
transported difference at every declared frame.

The radial readout premise (PR-18) and the shell-flux normalization premise
(PR-19) enter exactly as in the committed claim, as the two fields of
`RadialFluxData`; the frame is a declared choice on the hyperboloid; the
receipt `rest_fiber_shell_receipt` consumes one typed antecedent bundle
naming both.  Nothing here attaches the rest fiber to physical position
(PR-52 open), produces a gravitational field, potential, or mass density,
joins the composed Einstein branch, selects a frame or the radial readout
from source data, or discharges PR-18/PR-19.  OL-B2 stays partial.
-/

namespace OPH.RestFiberShellTransport

open OPH.C1Lorentz OPH.C2Soldering OPH.InverseSquareShellLaw
open OPH.PrimitivePortFrameQuotient
open Real Module

noncomputable section

/-! ## Spatial dot-product algebra

Small exact expansion lemmas for the three-coordinate dot product; each is a
finite ring identity. -/

private theorem spatialDot_comm' (f g : Spatial) : spatialDot f g = spatialDot g f := by
  simp only [spatialDot, Fin.sum_univ_three]
  ring

private theorem spatialDot_shear_right (w f : Spatial) (a : ℝ) :
    spatialDot w (fun i => f i + a * w i) = spatialDot w f + a * spatialDot w w := by
  simp only [spatialDot, Fin.sum_univ_three]
  ring

private theorem spatialDot_section_right (w f : Spatial) (a : ℝ) :
    spatialDot w (fun i => f i - a * w i) = spatialDot w f - a * spatialDot w w := by
  simp only [spatialDot, Fin.sum_univ_three]
  ring

private theorem spatialDot_shear_pair (f g w : Spatial) (a b : ℝ) :
    spatialDot (fun i => f i + a * w i) (fun i => g i + b * w i) =
      spatialDot f g + b * spatialDot f w + a * spatialDot w g +
        a * b * spatialDot w w := by
  simp only [spatialDot, Fin.sum_univ_three]
  ring

/-! ## Frame scalars -/

/-- The frame scalar coordinate shifted by one is nonzero: frames are future
unit-timelike, so the scalar coordinate is at least one. -/
theorem frame_succ_ne_zero (u : FrameHyperboloid) : u.1.1 + 1 ≠ 0 :=
  ne_of_gt (by linarith [frame_time_ge_one u])

/-- The spatial self dot of a frame is the squared scalar coordinate minus
one: the hyperboloid normalization in dot-product form. -/
theorem frame_spatialDot_self (u : FrameHyperboloid) :
    spatialDot u.1.2 u.1.2 = u.1.1 ^ 2 - 1 := by
  have hsq := frame_time_sq_eq_one_add_spatial u
  have hself : spatialDot u.1.2 u.1.2 = spatialNormSq u.1.2 := by
    simp only [spatialDot, spatialNormSq, Fin.sum_univ_three]
    ring
  linarith

/-! ## The boost-shear map

For a frame `(γ, w)` the shear sends a spatial coordinate vector `X` to the
ambient vector `(⟨w, X⟩, X + (⟨w, X⟩ / (γ + 1)) • w)`.  This is the
standard-frame rest embedding `X ↦ (0, X)` in boost-shear form at the frame
`(γ, w)`; every property needed below is an exact finite identity. -/

/-- Boost-shear of a spatial coordinate vector into the ambient module. -/
def shearVec (w : Spatial) (γ : ℝ) (X : Spatial) : Herm2 :=
  (spatialDot w X, fun i => X i + spatialDot w X / (γ + 1) * w i)

private theorem shearVec_lorentzB (w : Spatial) (γ : ℝ) (hγ : γ + 1 ≠ 0)
    (hw : spatialDot w w = γ ^ 2 - 1) (X : Spatial) :
    lorentzB (γ, w) (shearVec w γ X) = 0 := by
  show γ * spatialDot w X -
    spatialDot w (fun i => X i + spatialDot w X / (γ + 1) * w i) = 0
  rw [spatialDot_shear_right w X (spatialDot w X / (γ + 1)), hw]
  field_simp
  ring

private theorem shearVec_metric (w : Spatial) (γ : ℝ) (hγ : γ + 1 ≠ 0)
    (hw : spatialDot w w = γ ^ 2 - 1) (X Y : Spatial) :
    -lorentzB (shearVec w γ X) (shearVec w γ Y) = spatialDot X Y := by
  show -(spatialDot w X * spatialDot w Y -
      spatialDot (fun i => X i + spatialDot w X / (γ + 1) * w i)
        (fun i => Y i + spatialDot w Y / (γ + 1) * w i)) = spatialDot X Y
  rw [spatialDot_shear_pair X Y w (spatialDot w X / (γ + 1))
    (spatialDot w Y / (γ + 1)), hw, spatialDot_comm' X w]
  field_simp
  ring

private theorem shearVec_section (w : Spatial) (γ : ℝ) (hγ : γ + 1 ≠ 0)
    (hw : spatialDot w w = γ ^ 2 - 1) (v : Herm2)
    (hv : lorentzB (γ, w) v = 0) :
    shearVec w γ (fun i => v.2 i - v.1 / (γ + 1) * w i) = v := by
  have hv' : γ * v.1 - spatialDot w v.2 = 0 := hv
  have hs : spatialDot w v.2 = γ * v.1 := by linarith
  have hd : spatialDot w (fun i => v.2 i - v.1 / (γ + 1) * w i) = v.1 := by
    rw [spatialDot_section_right w v.2 (v.1 / (γ + 1)), hs, hw]
    field_simp
    ring
  unfold shearVec
  rw [hd]
  apply Prod.ext
  · rfl
  · funext i
    dsimp only
    ring

private theorem section_shearVec (w : Spatial) (γ : ℝ) (X : Spatial) (i : Fin 3) :
    (shearVec w γ X).2 i - (shearVec w γ X).1 / (γ + 1) * w i = X i := by
  unfold shearVec
  dsimp only
  ring

private theorem shearVec_add (w : Spatial) (γ : ℝ) (X Y : Spatial) :
    shearVec w γ (fun i => X i + Y i) = shearVec w γ X + shearVec w γ Y := by
  have hd : spatialDot w (fun i => X i + Y i) =
      spatialDot w X + spatialDot w Y := by
    simp only [spatialDot, Fin.sum_univ_three]
    ring
  unfold shearVec
  rw [hd]
  apply Prod.ext
  · rfl
  · funext i
    show X i + Y i + (spatialDot w X + spatialDot w Y) / (γ + 1) * w i =
      (X i + spatialDot w X / (γ + 1) * w i) +
        (Y i + spatialDot w Y / (γ + 1) * w i)
    ring

private theorem shearVec_smul (w : Spatial) (γ : ℝ) (a : ℝ) (X : Spatial) :
    shearVec w γ (fun i => a * X i) = a • shearVec w γ X := by
  have hd : spatialDot w (fun i => a * X i) = a * spatialDot w X := by
    simp only [spatialDot, Fin.sum_univ_three]
    ring
  unfold shearVec
  rw [hd]
  apply Prod.ext
  · rfl
  · funext i
    show a * X i + a * spatialDot w X / (γ + 1) * w i =
      a * (X i + spatialDot w X / (γ + 1) * w i)
    ring

/-! ## The typed isometry from the completion carrier onto every rest fiber -/

/-- Boost-shear embedding of the committed completion carrier into the
ambient module at one declared frame. -/
def restShear (u : FrameHyperboloid) (x : EuclideanVec3) : Herm2 :=
  shearVec u.1.2 u.1.1 (fun i => x i)

theorem restShear_mem (u : FrameHyperboloid) (x : EuclideanVec3) :
    restShear u x ∈ restSubspace u := by
  change lorentzB u.1 (restShear u x) = 0
  exact shearVec_lorentzB u.1.2 u.1.1 (frame_succ_ne_zero u)
    (frame_spatialDot_self u) (fun i => x i)

/-- Explicit section of the boost shear at one declared frame. -/
def restSection (u : FrameHyperboloid) (v : RestSpace u) : EuclideanVec3 :=
  WithLp.toLp 2 (fun i => v.1.2 i - v.1.1 / (u.1.1 + 1) * u.1.2 i)

/-- The explicit typed isometry: a linear equivalence from the shell-law
carrier's ambient space onto the three-dimensional rest space of the
declared frame.  The metric identity is `restFiberIsometry_metric`. -/
def restFiberIsometry (u : FrameHyperboloid) :
    EuclideanVec3 ≃ₗ[ℝ] RestSpace u where
  toFun x := ⟨restShear u x, restShear_mem u x⟩
  invFun v := restSection u v
  left_inv x := by
    apply PiLp.ext
    intro i
    exact section_shearVec u.1.2 u.1.1 (fun j => x j) i
  right_inv v := by
    apply Subtype.ext
    have hv := v.2
    change lorentzB u.1 v.1 = 0 at hv
    show restShear u (restSection u v) = v.1
    exact shearVec_section u.1.2 u.1.1 (frame_succ_ne_zero u)
      (frame_spatialDot_self u) v.1 hv
  map_add' x y := by
    apply Subtype.ext
    show restShear u (x + y) = restShear u x + restShear u y
    exact shearVec_add u.1.2 u.1.1 (fun i => x i) (fun i => y i)
  map_smul' a x := by
    apply Subtype.ext
    show restShear u (a • x) = a • restShear u x
    exact shearVec_smul u.1.2 u.1.1 a (fun i => x i)

@[simp] theorem restFiberIsometry_coe (u : FrameHyperboloid) (x : EuclideanVec3) :
    ((restFiberIsometry u x : RestSpace u) : Herm2) = restShear u x := rfl

/-- Exact metric identity: the pullback of the committed positive-definite
rest metric along the equivalence is the Euclidean coordinate form of the
completion carrier.  This is the sense in which `restFiberIsometry` is an
isometry rather than a bare linear equivalence. -/
theorem restFiberIsometry_metric (u : FrameHyperboloid) (x y : EuclideanVec3) :
    restMetric u (restFiberIsometry u x) (restFiberIsometry u y) =
      ∑ i : Fin 3, x i * y i := by
  show -lorentzB (restShear u x) (restShear u y) = ∑ i : Fin 3, x i * y i
  exact shearVec_metric u.1.2 u.1.1 (frame_succ_ne_zero u)
    (frame_spatialDot_self u) (fun i => x i) (fun i => y i)

/-- At the standard frame the boost shear degenerates to the committed
candidate bridge `vec3EquivStandardRest`: the boost-shear equivalence extends the
committed screen-to-rest-fiber identification from the standard frame to
every declared frame. -/
theorem restFiberIsometry_standardFrame (x : EuclideanVec3) :
    restFiberIsometry standardFrame x = vec3EquivStandardRest (fun i => x i) := by
  apply Subtype.ext
  show restShear standardFrame x = ((0 : ℝ), fun i => x i)
  unfold restShear shearVec
  apply Prod.ext
  · show spatialDot standardFrame.1.2 (fun i => x i) = 0
    simp [standardFrame, spatialDot]
  · funext i
    show x i + spatialDot standardFrame.1.2 (fun j => x j) /
        (standardFrame.1.1 + 1) * standardFrame.1.2 i = x i
    simp [standardFrame, spatialDot]

/-! ## Exponent provenance on the module's own dimension -/

/-- The rest-space dimension equals the carrier dimension, through the
transport equivalence itself. -/
theorem restSpace_finrank_eq_carrier (u : FrameHyperboloid) :
    finrank ℝ (RestSpace u) = finrank ℝ EuclideanVec3 :=
  ((restFiberIsometry u).finrank_eq).symm

/-- The rest space is three-dimensional, re-derived through the transport
equivalence and the committed carrier dimension theorem; it agrees with the
committed `finrank_restSpace`. -/
theorem restSpace_finrank_three (u : FrameHyperboloid) :
    finrank ℝ (RestSpace u) = 3 := by
  rw [restSpace_finrank_eq_carrier, carrier_finrank]

/-- The committed ball-volume scaling with the exponent read on the module's
own rest-space dimension. -/
theorem ballVolume_restSpace_scaling (u : FrameHyperboloid) {r : ℝ}
    (hr : 0 ≤ r) :
    ballVolume r = (π * 4 / 3) * r ^ (finrank ℝ (RestSpace u)) := by
  rw [restSpace_finrank_three u]
  exact ballVolume_eq hr

/-- The committed shell content with coefficient and exponent read on the
module's own rest-space dimension: the exponent is the rest-space dimension
minus one. -/
theorem shellContent_restSpace_exponent (u : FrameHyperboloid) {r : ℝ}
    (hr : 0 < r) :
    shellContent r = (finrank ℝ (RestSpace u) : ℝ) * (π * 4 / 3) *
      r ^ (finrank ℝ (RestSpace u) - 1) := by
  rw [restSpace_finrank_three u, shellContent_eq hr]
  have h2 : (3 : ℕ) - 1 = 2 := rfl
  rw [h2]
  push_cast
  ring

/-! ## The transported shell law -/

/-- Declared radial readout on a rest fiber: the square root of the
committed positive-definite rest metric. -/
def restRadius (u : FrameHyperboloid) (v : RestSpace u) : ℝ :=
  Real.sqrt (restMetric u v v)

theorem restMetric_self_nonneg (u : FrameHyperboloid) (v : RestSpace u) :
    0 ≤ restMetric u v v := by
  rcases eq_or_ne v 0 with rfl | hv
  · exact le_of_eq ((restMetric_self_eq_zero_iff u 0).mpr rfl).symm
  · exact le_of_lt (restMetric_pos u hv)

theorem restRadius_sq (u : FrameHyperboloid) (v : RestSpace u) :
    restRadius u v ^ 2 = restMetric u v v :=
  Real.sq_sqrt (restMetric_self_nonneg u v)

theorem restRadius_pos (u : FrameHyperboloid) {v : RestSpace u} (hv : v ≠ 0) :
    0 < restRadius u v :=
  Real.sqrt_pos.mpr (restMetric_pos u hv)

/-- The radial readout pulled back through the isometry is the carrier
norm: the declared radial structure on the fiber restates the committed
Euclidean one. -/
theorem restFiberIsometry_radius (u : FrameHyperboloid) (x : EuclideanVec3) :
    restRadius u (restFiberIsometry u x) = ‖x‖ := by
  unfold restRadius
  rw [restFiberIsometry_metric, EuclideanSpace.norm_eq]
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  rw [Real.norm_eq_abs, sq_abs]
  ring

/-- The transported strength field on a rest fiber: the declared radial
strength profile (PR-18) read through the declared radial readout. -/
def restField (u : FrameHyperboloid) (D : RadialFluxData)
    (v : RestSpace u) : ℝ :=
  D.strength (restRadius u v)

/-- Transported inverse-square law: on every declared rest fiber the
strength at a nonzero vector is the charge over `4 π r ^ 2`, exactly, by
the committed shell-law theorem. -/
theorem restField_inverse_square (u : FrameHyperboloid) (D : RadialFluxData)
    {v : RestSpace u} (hv : v ≠ 0) :
    restField u D v = D.charge / (4 * π * restRadius u v ^ 2) :=
  D.inverse_square (restRadius_pos u hv)

/-- Metric form of the transported law: the strength is the charge over
`4 π` times the rest metric of the vector with itself. -/
theorem restField_inverse_square_metric (u : FrameHyperboloid)
    (D : RadialFluxData) {v : RestSpace u} (hv : v ≠ 0) :
    restField u D v = D.charge / (4 * π * restMetric u v v) := by
  rw [restField_inverse_square u D hv, restRadius_sq]

/-- Transported shell-flux normalization (PR-19 pulled back): strength times
transported shell content is the constant charge on every nonzero fiber
vector. -/
theorem restField_normalized (u : FrameHyperboloid) (D : RadialFluxData)
    {v : RestSpace u} (hv : v ≠ 0) :
    restField u D v * shellContent (restRadius u v) = D.charge :=
  D.normalized (restRadius u v) (restRadius_pos u hv)

/-- Transported scale-free form: strength times squared rest radius is the
same on every pair of nonzero fiber vectors. -/
theorem restField_scale_free (u : FrameHyperboloid) (D : RadialFluxData)
    {v w : RestSpace u} (hv : v ≠ 0) (hw : w ≠ 0) :
    restField u D v * restRadius u v ^ 2 =
      restField u D w * restRadius u w ^ 2 :=
  D.scale_free (restRadius_pos u hv) (restRadius_pos u hw)

/-- The transported exponent statement on the module's own dimension: the
strength times the rest radius to the power rest-space dimension minus one
is the charge over `4 π`. -/
theorem restField_restSpace_exponent (u : FrameHyperboloid)
    (D : RadialFluxData) {v : RestSpace u} (hv : v ≠ 0) :
    restField u D v * restRadius u v ^ (finrank ℝ (RestSpace u) - 1) =
      D.charge / (4 * π) := by
  have hr := restRadius_pos u hv
  have hπ : (0 : ℝ) < π := Real.pi_pos
  rw [restSpace_finrank_three u, restField_inverse_square u D hv]
  have h2 : (3 : ℕ) - 1 = 2 := rfl
  rw [h2]
  field_simp

/-- The transported flux law pulled back through the isometry, in carrier
coordinates: the same exponent, with the carrier norm as radius. -/
theorem restField_isometry_pullback (u : FrameHyperboloid)
    (D : RadialFluxData) (x : EuclideanVec3) :
    restField u D (restFiberIsometry u x) = D.strength ‖x‖ := by
  unfold restField
  rw [restFiberIsometry_radius]

theorem restField_isometry_inverse_square (u : FrameHyperboloid)
    (D : RadialFluxData) {x : EuclideanVec3} (hx : x ≠ 0) :
    restField u D (restFiberIsometry u x) = D.charge / (4 * π * ‖x‖ ^ 2) := by
  have hx' : restFiberIsometry u x ≠ 0 := fun h =>
    hx ((restFiberIsometry u).map_eq_zero_iff.mp h)
  rw [restField_inverse_square u D hx', restFiberIsometry_radius]

/-! ## Chart invariance under every oriented Lorentz map

The committed transport `OrientedLorentzEquiv.restEquiv` carries rest
fibers to rest fibers preserving the rest metric; therefore it carries the
declared radial readout, the transported strength field, and every
transported identity to the transported frame's rest space, exactly. -/

theorem restEquiv_radius (L : OrientedLorentzEquiv) (u : FrameHyperboloid)
    (v : RestSpace u) :
    restRadius (L.mapFrame u) (L.restEquiv u v) = restRadius u v := by
  unfold restRadius
  rw [L.restEquiv_preserves_metric]

theorem restEquiv_field (L : OrientedLorentzEquiv) (u : FrameHyperboloid)
    (D : RadialFluxData) (v : RestSpace u) :
    restField (L.mapFrame u) D (L.restEquiv u v) = restField u D v := by
  unfold restField
  rw [restEquiv_radius]

/-- The transported inverse-square law on the transported frame's rest
space, read through the committed chart transport. -/
theorem restEquiv_inverse_square (L : OrientedLorentzEquiv)
    (u : FrameHyperboloid) (D : RadialFluxData) {v : RestSpace u}
    (hv : v ≠ 0) :
    restField (L.mapFrame u) D (L.restEquiv u v) =
      D.charge / (4 * π * restRadius (L.mapFrame u) (L.restEquiv u v) ^ 2) := by
  have hv' : L.restEquiv u v ≠ 0 := fun h =>
    hv ((L.restEquiv u).map_eq_zero_iff.mp h)
  exact restField_inverse_square (L.mapFrame u) D hv'

/-- Composing the typed carrier isometry with any committed oriented chart
transport is again a metric-exact identification of the carrier with the
transported frame's rest space: the whole package is chart invariant. -/
theorem restEquiv_isometry_composition (L : OrientedLorentzEquiv)
    (u : FrameHyperboloid) (x y : EuclideanVec3) :
    restMetric (L.mapFrame u)
        (L.restEquiv u (restFiberIsometry u x))
        (L.restEquiv u (restFiberIsometry u y)) =
      ∑ i : Fin 3, x i * y i := by
  rw [L.restEquiv_preserves_metric, restFiberIsometry_metric]

/-! ## The screen chain carried into the fiber -/

/-- The committed screen chain composed with the typed fiber isometry: an
integer record control lands in the rest space of the declared frame. -/
def sourceChain (u : FrameHyperboloid) (p : IntegerFramePoint) :
    RestSpace u :=
  restFiberIsometry u (pointEuclideanFrame p)

/-- Composed chain identity: the finite source Gram of two integer record
controls equals the rest metric of their transported difference.  The left
side is the committed screen metric (`point_dist_sq_eq_sourceGram`); the
right side is the declared frame's rest fiber.  No physical-position
attachment is claimed. -/
theorem sourceChain_metric_eq_sourceGram (u : FrameHyperboloid)
    (p q : IntegerFramePoint) :
    restMetric u (sourceChain u p - sourceChain u q)
        (sourceChain u p - sourceChain u q) =
      sourceGram (pointControlDifference p q) (pointControlDifference p q) := by
  have hsub : sourceChain u p - sourceChain u q =
      restFiberIsometry u (pointEuclideanFrame p - pointEuclideanFrame q) := by
    unfold sourceChain
    rw [map_sub]
  rw [hsub, restFiberIsometry_metric, ← point_dist_sq_eq_sourceGram]
  change _ = dist (pointEuclideanFrame p) (pointEuclideanFrame q) ^ 2
  rw [EuclideanSpace.dist_sq_eq]
  apply Finset.sum_congr rfl
  intro i _
  rw [Real.dist_eq, sq_abs]
  show (pointEuclideanFrame p i - pointEuclideanFrame q i) *
      (pointEuclideanFrame p i - pointEuclideanFrame q i) =
    (pointEuclideanFrame p i - pointEuclideanFrame q i) ^ 2
  ring

/-! ## Composed receipt -/

/-- The typed antecedent bundle of the rest-fiber shell transport.  The
`flux` field is the committed premise packet of the shell law: its
`strength` field is the declared spherically symmetric radial readout
(premise register row PR-18) and its `normalized` field is the declared
shell-flux normalization (premise register row PR-19), consumed exactly as
in the committed claim `OPH-GEOMETRY-INVERSE-SQUARE-SHELL`.  The `frame`
field is the declared future unit-timelike frame whose rest space receives
the law; no source theorem selects it.  The transport isometry is the boost
shear of the frame, a function of the frame alone; no uniqueness of the
identification is claimed. -/
structure RestShellAntecedents where
  /-- PR-18 (radial readout) and PR-19 (shell-flux normalization), as the
  committed `RadialFluxData` premise packet. -/
  flux : RadialFluxData
  /-- The declared future unit-timelike frame. -/
  frame : FrameHyperboloid

/-- Composed receipt from the single typed antecedent bundle: the rest
space of the declared frame is three-dimensional; the typed carrier
isometry is metric exact onto it; the declared radial readout pulled back
through the isometry is the carrier norm; the transported inverse-square
law holds on it with the committed exponent; the exponent reads on the
module's own rest-space dimension; every committed oriented chart
transport carries the strength field exactly; and the committed screen
Gram chain lands in the fiber.  Each conjunct reuses a committed theorem;
the two flux premises stay declared. -/
theorem rest_fiber_shell_receipt (A : RestShellAntecedents) :
    finrank ℝ (RestSpace A.frame) = 3 ∧
    (∀ x y : EuclideanVec3,
      restMetric A.frame (restFiberIsometry A.frame x)
        (restFiberIsometry A.frame y) = ∑ i : Fin 3, x i * y i) ∧
    (∀ x : EuclideanVec3,
      restRadius A.frame (restFiberIsometry A.frame x) = ‖x‖) ∧
    (∀ v : RestSpace A.frame, v ≠ 0 →
      restField A.frame A.flux v =
        A.flux.charge / (4 * π * restRadius A.frame v ^ 2)) ∧
    (∀ v : RestSpace A.frame, v ≠ 0 →
      restField A.frame A.flux v *
          restRadius A.frame v ^ (finrank ℝ (RestSpace A.frame) - 1) =
        A.flux.charge / (4 * π)) ∧
    (∀ (L : OrientedLorentzEquiv) (v : RestSpace A.frame),
      restField (L.mapFrame A.frame) A.flux (L.restEquiv A.frame v) =
        restField A.frame A.flux v) ∧
    (∀ p q : IntegerFramePoint,
      restMetric A.frame (sourceChain A.frame p - sourceChain A.frame q)
          (sourceChain A.frame p - sourceChain A.frame q) =
        sourceGram (pointControlDifference p q)
          (pointControlDifference p q)) :=
  ⟨restSpace_finrank_three A.frame,
   restFiberIsometry_metric A.frame,
   restFiberIsometry_radius A.frame,
   fun _ hv => restField_inverse_square A.frame A.flux hv,
   fun _ hv => restField_restSpace_exponent A.frame A.flux hv,
   fun L v => restEquiv_field L A.frame A.flux v,
   sourceChain_metric_eq_sourceGram A.frame⟩

end

end OPH.RestFiberShellTransport

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.RestFiberShellTransport.restShear_mem
#print axioms OPH.RestFiberShellTransport.restFiberIsometry_metric
#print axioms OPH.RestFiberShellTransport.restFiberIsometry_standardFrame
#print axioms OPH.RestFiberShellTransport.restSpace_finrank_three
#print axioms OPH.RestFiberShellTransport.shellContent_restSpace_exponent
#print axioms OPH.RestFiberShellTransport.restFiberIsometry_radius
#print axioms OPH.RestFiberShellTransport.restField_inverse_square
#print axioms OPH.RestFiberShellTransport.restField_normalized
#print axioms OPH.RestFiberShellTransport.restField_scale_free
#print axioms OPH.RestFiberShellTransport.restField_restSpace_exponent
#print axioms OPH.RestFiberShellTransport.restField_isometry_inverse_square
#print axioms OPH.RestFiberShellTransport.restEquiv_field
#print axioms OPH.RestFiberShellTransport.restEquiv_inverse_square
#print axioms OPH.RestFiberShellTransport.restEquiv_isometry_composition
#print axioms OPH.RestFiberShellTransport.sourceChain_metric_eq_sourceGram
#print axioms OPH.RestFiberShellTransport.rest_fiber_shell_receipt
