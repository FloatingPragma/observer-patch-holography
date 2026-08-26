import Geometry.InternalEnergyInertia
import Geometry.ChargeFixedInteraction

/-!
# Frame covariance selects the slope of internal energy in the composite
four-momentum

STATUS.  Kinematic selection on the declared Lorentz module.  Among the
declared one-parameter slope family of composite momentum maps, frame
covariance under oriented Lorentz maps holds exactly at slope one or at
zero internal energy.  The selection consumes the declared family and the
covariance requirement; it is a statement about the momentum map, not about
dynamics, and it does not derive a mass-energy identity.

WHAT IS PROVED.

(1) Declared slope family.  `slopeMomentum m E lam frame = fourMomentum
(m + lam * E) frame + ((1 - lam) * E) • restVector`, with `restVector =
(1, 0)` the vector of the standard frame (`restVector_eq_standardFrame`,
`restVector_eq_fourMomentum_one`).  Endpoints: slope one is the shape-A
composite momentum `fourMomentum (m + E) frame` (`slopeMomentum_one`);
slope zero is `fourMomentum m frame + E • restVector` (`slopeMomentum_zero`).
At the standard frame every member has scalar coordinate `m + E`
(`slopeMomentum_standardFrame_time`) and in fact equals `(m + E, 0)`
(`slopeMomentum_standardFrame`): the whole family agrees at rest.

(2) Covariance.  `FrameCovariant P` means `P (L.mapFrame frame) = L (P frame)`
for every oriented Lorentz map `L`.  The transport identity
`slopeMomentum_mapFrame` shows the family transports up to the term
`((1 - lam) * E) • (restVector - L restVector)`.  An explicit boost
`boost a b` on the declared module (parameters with `a ^ 2 - b ^ 2 = 1`,
`0 < a`; the instance `standardBoost` has `a = 5/4`, `b = 3/4`) is proved
oriented Lorentz and moves `restVector` (`standardBoost_moves_restVector`).
Main theorem `frameCovariant_slopeMomentum_iff`: `slopeMomentum m E lam`
is frame covariant iff `lam = 1 ∨ E = 0`.

(3) Mass shell.  Exact Lorentz square `lorentzQ_slopeMomentum`:
`lorentzQ (slopeMomentum m E lam frame) = (m + E) ^ 2 +
2 * ((1 - lam) * E) * (m + lam * E) * (frame.1.1 - 1)`; the defect is
displayed against the spatial part through the mass-shell identity
(`slopeMomentum_shell_defect_spatial`).  `slopeMomentum_shell_all_iff`:
the shell identity at the rest energy holds at every frame iff
`(1 - lam) * E * (m + lam * E) = 0`; under `0 < m`, `0 < m + E`,
`0 ≤ lam ≤ 1` this is `lam = 1 ∨ E = 0`
(`slopeMomentum_shell_all_iff_of_pos`).

(4) Linearity.  `slopeMomentum_add` and `slopeMomentum_smul`: at fixed
slope and frame the map is additive and homogeneous in `(m, E)`.  At slope
one the sum of two composites is the shape-A momentum at the summed rest
energy with invariant mass the summed rest energy
(`slopeMomentum_one_add_shapeA`); for `(1 - lam) * E ≠ 0` and a frame with
nonzero spatial part the member is no `fourMomentum` of any parameter
(`slopeMomentum_not_fourMomentum`).

(5) Join.  `compositeFourMomentum m E frame = slopeMomentum m E 1 frame`
(`compositeFourMomentum_eq_slopeMomentum_one`), with no factor;
`discreteMomentum_eq_two_slopeMomentum`: the discrete momentum of the
quadratic worldline action equals `2 • slopeMomentum m 0 1 frame` on a
frame increment, the factor `2` being the inertial coefficient `2 m` of
that action.

(6) Conclusion `covariance_selects_slope_one`: a conjunction of (1), (2),
(5) exactly scoped, together with the open-row register.

ROWS TOUCHED.  Coupled-action row (the composite shape is declared; this
module reads the momentum map of the kinematics island only), physical
spacetime attachment row (`Herm2` and `FrameHyperboloid` are the declared
module and hyperboloid; the frame is not attached to a physical observer),
source clock and duration row (the frame increment of the discrete momentum
join is the declared step), and the laboratory clock and energy calibration
import (the ledger `E` is a declared scalar in tick-energy units; its
identification with a physical internal energy, and `m` with a physical
mass, are the import).  None discharged.

NEGATIVES CITED.  The Legendre non-identifiability at its scope, re-exported
as `legendre_scope_cited`: realized histories select no velocity curvature
or Legendre map, so every shape in the slope family is a declared
enrichment.  The selection here is kinematic, on the momentum map, and is
independent of and consistent with the dynamical shape selection of
`Geometry/InternalEnergyInertia.lean`.

CONVENTIONS.  Signature `(+---)`, `Herm2 = ℝ × (Fin 3 → ℝ)`, invariant-speed
conversion one.  `boostFun a b` acts on the scalar coordinate and spatial
coordinate `0` by the matrix `[[a, b], [b, a]]` and fixes coordinates `1`,
`2`.  The frame `frame.1.1` is the scalar coordinate of the frame vector.

FALSIFIER.  The module fails if the transport identity misses a term, if
`standardBoost` fixes `restVector`, if the exact Lorentz square differs
from the displayed formula, if the shape-A join carries a hidden factor,
or if a member with `lam ≠ 1` and `E ≠ 0` is frame covariant.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.CompositeMomentumCovariance

open OPH.C1Lorentz OPH.C2Soldering OPH.InternalEnergyInertia

noncomputable section

/-! ## An explicit boost on the declared module -/

/-- The declared boost map along spatial coordinate `0`, with parameters
`a`, `b`. -/
def boostFun (a b : ℝ) (v : Herm2) : Herm2 :=
  (a * v.1 + b * v.2 0, ![b * v.1 + a * v.2 0, v.2 1, v.2 2])

theorem boostFun_fst (a b : ℝ) (v : Herm2) :
    (boostFun a b v).1 = a * v.1 + b * v.2 0 := rfl

theorem boostFun_snd_zero (a b : ℝ) (v : Herm2) :
    (boostFun a b v).2 0 = b * v.1 + a * v.2 0 := rfl

theorem boostFun_snd_one (a b : ℝ) (v : Herm2) :
    (boostFun a b v).2 1 = v.2 1 := rfl

theorem boostFun_snd_two (a b : ℝ) (v : Herm2) :
    (boostFun a b v).2 2 = v.2 2 := rfl

/-- The boost preserves the Lorentz pairing when `a ^ 2 - b ^ 2 = 1`. -/
theorem lorentzB_boostFun (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (u v : Herm2) :
    lorentzB (boostFun a b u) (boostFun a b v) = lorentzB u v := by
  simp only [lorentzB, spatialDot, Fin.sum_univ_three, boostFun_fst,
    boostFun_snd_zero, boostFun_snd_one, boostFun_snd_two]
  linear_combination (u.1 * v.1 - u.2 0 * v.2 0) * h

theorem boostFun_boostFun (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (v : Herm2) :
    boostFun a (-b) (boostFun a b v) = v := by
  apply Prod.ext
  · simp only [boostFun_fst, boostFun_snd_zero]
    linear_combination v.1 * h
  · funext i
    fin_cases i
    · show -b * (a * v.1 + b * v.2 0) + a * (b * v.1 + a * v.2 0) = v.2 0
      linear_combination v.2 0 * h
    · rfl
    · rfl

theorem boostFun_add (a b : ℝ) (u v : Herm2) :
    boostFun a b (u + v) = boostFun a b u + boostFun a b v := by
  apply Prod.ext
  · simp only [boostFun_fst, Prod.fst_add, Prod.snd_add, Pi.add_apply]
    ring
  · funext i
    fin_cases i
    · show b * (u + v).1 + a * (u + v).2 0 = (b * u.1 + a * u.2 0) + (b * v.1 + a * v.2 0)
      simp only [Prod.fst_add, Prod.snd_add, Pi.add_apply]
      ring
    · rfl
    · rfl

theorem boostFun_smul (a b c : ℝ) (v : Herm2) :
    boostFun a b (c • v) = c • boostFun a b v := by
  apply Prod.ext
  · simp only [boostFun_fst, Prod.smul_fst, Prod.smul_snd, Pi.smul_apply,
      smul_eq_mul]
    ring
  · funext i
    fin_cases i
    · show b * (c • v).1 + a * (c • v).2 0 = c * (b * v.1 + a * v.2 0)
      simp only [Prod.smul_fst, Prod.smul_snd, Pi.smul_apply, smul_eq_mul]
      ring
    · rfl
    · rfl

/-- The boost as a linear equivalence, for `a ^ 2 - b ^ 2 = 1`. -/
def boostEquiv (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) : Herm2 ≃ₗ[ℝ] Herm2 where
  toFun := boostFun a b
  invFun := boostFun a (-b)
  map_add' := boostFun_add a b
  map_smul' := boostFun_smul a b
  left_inv := boostFun_boostFun a b h
  right_inv := by
    intro v
    have h' : a ^ 2 - (-b) ^ 2 = 1 := by rw [neg_sq]; exact h
    have := boostFun_boostFun a (-b) h' v
    rwa [neg_neg] at this

theorem coord_zero_sq_le_spatialNormSq (x : Spatial) :
    (x 0) ^ 2 ≤ spatialNormSq x :=
  Finset.single_le_sum (fun j (_ : j ∈ Finset.univ) ↦ sq_nonneg (x j))
    (Finset.mem_univ 0)

/-- A boost with `0 < a` keeps the scalar coordinate of a future causal
vector positive. -/
theorem boostFun_time_pos (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a)
    {v : Herm2} (hq : 0 ≤ lorentzQ v) (ht : 0 < v.1) :
    0 < (boostFun a b v).1 := by
  rw [boostFun_fst]
  have hx : (v.2 0) ^ 2 ≤ spatialNormSq v.2 := coord_zero_sq_le_spatialNormSq v.2
  have hvq : spatialNormSq v.2 ≤ v.1 ^ 2 := by
    simp only [lorentzQ] at hq
    linarith
  have hprod : 0 < (a * v.1 + b * v.2 0) * (a * v.1 - b * v.2 0) := by
    have hb : 0 ≤ b ^ 2 := sq_nonneg b
    have h1 : (a * v.1 + b * v.2 0) * (a * v.1 - b * v.2 0) =
        v.1 ^ 2 + b ^ 2 * (v.1 ^ 2 - (v.2 0) ^ 2) := by
      linear_combination (v.1 ^ 2) * h
    rw [h1]
    have h2 : 0 ≤ b ^ 2 * (v.1 ^ 2 - (v.2 0) ^ 2) :=
      mul_nonneg hb (by linarith)
    have h3 : 0 < v.1 ^ 2 := pow_pos ht 2
    linarith
  by_contra hneg
  have hneg : a * v.1 + b * v.2 0 ≤ 0 := le_of_not_gt hneg
  have hsum : 0 < a * v.1 := mul_pos ha ht
  have hminus : 0 < a * v.1 - b * v.2 0 := by linarith
  have : (a * v.1 + b * v.2 0) * (a * v.1 - b * v.2 0) ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg hneg hminus.le
  linarith

theorem lorentzQ_boostFun (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (v : Herm2) :
    lorentzQ (boostFun a b v) = lorentzQ v := by
  rw [← lorentzB_self, lorentzB_boostFun a b h, lorentzB_self]

theorem boostFun_futureNull (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a)
    {v : Herm2} (hv : IsFutureNull v) : IsFutureNull (boostFun a b v) :=
  ⟨by rw [lorentzQ_boostFun a b h]; exact hv.1,
    boostFun_time_pos a b h ha hv.1.ge hv.2⟩

theorem boostFun_futureUnit (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a)
    {v : Herm2} (hv : IsFutureUnitTimelike v) :
    IsFutureUnitTimelike (boostFun a b v) :=
  ⟨by rw [lorentzQ_boostFun a b h]; exact hv.1,
    boostFun_time_pos a b h ha (by rw [hv.1]; exact zero_le_one) hv.2⟩

/-- The boost as an oriented Lorentz map, for `a ^ 2 - b ^ 2 = 1` and
`0 < a`. -/
def boost (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a) :
    OrientedLorentzEquiv where
  toLinearEquiv := boostEquiv a b h
  map_lorentzB := lorentzB_boostFun a b h
  map_futureNull_iff := by
    intro v
    constructor
    · intro hv
      have h' : a ^ 2 - (-b) ^ 2 = 1 := by rw [neg_sq]; exact h
      have := boostFun_futureNull a (-b) h' ha hv
      rwa [show (boostEquiv a b h) v = boostFun a b v from rfl,
        boostFun_boostFun a b h] at this
    · exact boostFun_futureNull a b h ha
  map_futureUnit_iff := by
    intro v
    constructor
    · intro hv
      have h' : a ^ 2 - (-b) ^ 2 = 1 := by rw [neg_sq]; exact h
      have := boostFun_futureUnit a (-b) h' ha hv
      rwa [show (boostEquiv a b h) v = boostFun a b v from rfl,
        boostFun_boostFun a b h] at this
    · exact boostFun_futureUnit a b h ha

theorem boost_apply (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a)
    (v : Herm2) : (boost a b h ha) v = boostFun a b v := rfl

/-- The rest direction of the declared module: the vector of the standard
frame. -/
def restVector : Herm2 := (1, 0)

theorem restVector_eq_standardFrame : restVector = (standardFrame : Herm2) := rfl

theorem restVector_eq_fourMomentum_one :
    restVector = fourMomentum 1 standardFrame := by
  rw [fourMomentum_standardFrame]; rfl

theorem boost_restVector (a b : ℝ) (h : a ^ 2 - b ^ 2 = 1) (ha : 0 < a) :
    (boost a b h ha) restVector = (a, ![b, 0, 0]) := by
  rw [boost_apply]
  apply Prod.ext
  · show a * 1 + b * 0 = a
    ring
  · funext i
    fin_cases i
    · show b * 1 + a * 0 = b
      ring
    · rfl
    · rfl

/-- The explicit instance with `a = 5/4`, `b = 3/4`. -/
def standardBoost : OrientedLorentzEquiv :=
  boost (5 / 4) (3 / 4) (by norm_num) (by norm_num)

theorem standardBoost_restVector :
    standardBoost restVector = (5 / 4, ![3 / 4, 0, 0]) :=
  boost_restVector (5 / 4) (3 / 4) (by norm_num) (by norm_num)

/-- The explicit boost moves the rest direction. -/
theorem standardBoost_moves_restVector : standardBoost restVector ≠ restVector := by
  rw [standardBoost_restVector]
  intro heq
  have h0 := congrArg Prod.fst heq
  change (5 / 4 : ℝ) = 1 at h0
  norm_num at h0

/-! ## (1) The declared slope family -/

/-- **Declared slope family.**  Internal energy `E` enters the composite
four-momentum with slope `lam` on the frame direction and slope `1 - lam`
on the rest direction. -/
def slopeMomentum (m E lam : ℝ) (frame : FrameHyperboloid) : Herm2 :=
  fourMomentum (m + lam * E) frame + ((1 - lam) * E) • restVector

/-- Slope one is the shape-A composite momentum. -/
theorem slopeMomentum_one (m E : ℝ) (frame : FrameHyperboloid) :
    slopeMomentum m E 1 frame = fourMomentum (m + E) frame := by
  simp [slopeMomentum]

/-- Slope zero is the declared-mass momentum plus the internal energy on
the rest direction: the momentum-map analogue of the zero endpoint of the
slope family named in `Geometry/InternalEnergyInertia.lean` (the action-level
shape B there carries no ledger term in its momentum). -/
theorem slopeMomentum_zero (m E : ℝ) (frame : FrameHyperboloid) :
    slopeMomentum m E 0 frame = fourMomentum m frame + E • restVector := by
  simp [slopeMomentum]

theorem slopeMomentum_time (m E lam : ℝ) (frame : FrameHyperboloid) :
    (slopeMomentum m E lam frame).1 = (m + lam * E) * frame.1.1 + (1 - lam) * E := by
  simp [slopeMomentum, restVector]

theorem slopeMomentum_spatial (m E lam : ℝ) (frame : FrameHyperboloid) :
    (slopeMomentum m E lam frame).2 = (m + lam * E) • frame.1.2 := by
  simp [slopeMomentum, restVector]

/-- Every member equals `(m + E, 0)` at the standard frame. -/
theorem slopeMomentum_standardFrame (m E lam : ℝ) :
    slopeMomentum m E lam standardFrame = (m + E, 0) := by
  apply Prod.ext
  · rw [slopeMomentum_time]
    show (m + lam * E) * 1 + (1 - lam) * E = m + E
    ring
  · rw [slopeMomentum_spatial]
    show (m + lam * E) • (0 : Spatial) = 0
    simp

/-- The rest energy of every member is `m + E`. -/
theorem slopeMomentum_standardFrame_time (m E lam : ℝ) :
    (slopeMomentum m E lam standardFrame).1 = m + E := by
  rw [slopeMomentum_standardFrame]

/-! ## (2) Frame covariance -/

/-- A frame-indexed vector is frame covariant when every oriented Lorentz
map transports it along the transported frame. -/
def FrameCovariant (P : FrameHyperboloid → Herm2) : Prop :=
  ∀ (L : OrientedLorentzEquiv) (frame : FrameHyperboloid),
    P (L.mapFrame frame) = L (P frame)

/-- The shape-A momentum is frame covariant (the corpus transport law). -/
theorem frameCovariant_fourMomentum (mass : ℝ) :
    FrameCovariant (fourMomentum mass) :=
  fun L frame => (oriented_map_fourMomentum L mass frame).symm

/-- Exact transport identity for the slope family. -/
theorem slopeMomentum_mapFrame (m E lam : ℝ) (L : OrientedLorentzEquiv)
    (frame : FrameHyperboloid) :
    slopeMomentum m E lam (L.mapFrame frame) =
      L (slopeMomentum m E lam frame) +
        ((1 - lam) * E) • (restVector - L restVector) := by
  unfold slopeMomentum
  rw [OrientedLorentzEquiv.map_add, OrientedLorentzEquiv.map_smul,
    ← oriented_map_fourMomentum, smul_sub]
  abel

theorem frameCovariant_slopeMomentum_of (m E lam : ℝ)
    (h : (1 - lam) * E = 0) : FrameCovariant (slopeMomentum m E lam) := by
  intro L frame
  rw [slopeMomentum_mapFrame, h, zero_smul, add_zero]

/-- **Covariance selects slope one.**  The slope-`lam` member is frame
covariant iff `lam = 1` or the internal energy vanishes. -/
theorem frameCovariant_slopeMomentum_iff (m E lam : ℝ) :
    FrameCovariant (slopeMomentum m E lam) ↔ lam = 1 ∨ E = 0 := by
  constructor
  · intro hcov
    have hkey := hcov standardBoost standardFrame
    rw [slopeMomentum_mapFrame] at hkey
    have hz : ((1 - lam) * E) • (restVector - standardBoost restVector) = 0 := by
      have := congrArg (fun w => w - standardBoost (slopeMomentum m E lam standardFrame)) hkey
      simpa using this
    rcases smul_eq_zero.mp hz with hc | hd
    · rcases mul_eq_zero.mp hc with h1 | h2
      · left; linarith
      · right; exact h2
    · exact absurd (sub_eq_zero.mp hd).symm standardBoost_moves_restVector
  · intro h
    apply frameCovariant_slopeMomentum_of
    rcases h with h | h
    · rw [h]; ring
    · rw [h]; ring

/-! ## (3) Mass-shell consequence -/

/-- Exact Lorentz square of the slope family: the shell value at the rest
energy plus a defect proportional to `(1 - lam) * E`, to `m + lam * E`,
and to the excess of the frame's scalar coordinate over one. -/
theorem lorentzQ_slopeMomentum (m E lam : ℝ) (frame : FrameHyperboloid) :
    lorentzQ (slopeMomentum m E lam frame) =
      (m + E) ^ 2 + 2 * ((1 - lam) * E) * (m + lam * E) * (frame.1.1 - 1) := by
  have hframe := frame_time_sq_eq_one_add_spatial frame
  rw [lorentzQ, slopeMomentum_time, slopeMomentum_spatial, spatialNormSq_smul]
  have hs : spatialNormSq frame.1.2 = frame.1.1 ^ 2 - 1 := by linarith
  rw [hs]
  ring

/-- The shell defect against the spatial part of the shape-A momentum at
parameter `m + lam * E`, through the kinematic mass-shell identity. -/
theorem slopeMomentum_shell_defect_spatial (m E lam : ℝ) (frame : FrameHyperboloid) :
    (lorentzQ (slopeMomentum m E lam frame) - (m + E) ^ 2) *
        ((fourMomentum (m + lam * E) frame).1 + (m + lam * E)) =
      2 * ((1 - lam) * E) * spatialNormSq (fourMomentum (m + lam * E) frame).2 := by
  have hshell := fourMomentum_time_sq_eq_mass_sq_add_spatial (m + lam * E) frame
  rw [fourMomentum_time] at hshell
  rw [lorentzQ_slopeMomentum, fourMomentum_time]
  have hsp : spatialNormSq (fourMomentum (m + lam * E) frame).2 =
      ((m + lam * E) * frame.1.1) ^ 2 - (m + lam * E) ^ 2 := by linarith
  rw [hsp]
  ring

/-- The boosted standard frame has scalar coordinate `5 / 4`. -/
theorem standardBoost_standardFrame_time :
    (standardBoost.mapFrame standardFrame).1.1 = 5 / 4 := by
  rw [OrientedLorentzEquiv.mapFrame_coe, ← restVector_eq_standardFrame,
    standardBoost_restVector]

/-- The shell identity at the rest energy holds at every frame iff the
product `(1 - lam) * E * (m + lam * E)` vanishes. -/
theorem slopeMomentum_shell_all_iff (m E lam : ℝ) :
    (∀ frame : FrameHyperboloid,
        lorentzQ (slopeMomentum m E lam frame) = (m + E) ^ 2) ↔
      (1 - lam) * E * (m + lam * E) = 0 := by
  constructor
  · intro hall
    have h := hall (standardBoost.mapFrame standardFrame)
    rw [lorentzQ_slopeMomentum, standardBoost_standardFrame_time] at h
    linarith
  · intro h frame
    rw [lorentzQ_slopeMomentum]
    have : 2 * ((1 - lam) * E) * (m + lam * E) * (frame.1.1 - 1) =
        ((1 - lam) * E * (m + lam * E)) * (2 * (frame.1.1 - 1)) := by ring
    rw [this, h, zero_mul, add_zero]

/-- For positive rest parameter and rest energy and slope in `[0, 1]`, the
member is on the shell of its rest energy at every frame iff `lam = 1` or
`E = 0`: only the slope-one member has one invariant mass equal to its rest
energy. -/
theorem slopeMomentum_shell_all_iff_of_pos {m E lam : ℝ} (hm : 0 < m)
    (hmE : 0 < m + E) (h0 : 0 ≤ lam) (h1 : lam ≤ 1) :
    (∀ frame : FrameHyperboloid,
        lorentzQ (slopeMomentum m E lam frame) = (m + E) ^ 2) ↔
      lam = 1 ∨ E = 0 := by
  rw [slopeMomentum_shell_all_iff]
  have hmu : 0 < m + lam * E := by
    have : m + lam * E = (1 - lam) * m + lam * (m + E) := by ring
    rw [this]
    rcases eq_or_lt_of_le h0 with hz | hpos
    · rw [← hz]; simpa using hm
    · have := mul_pos hpos hmE
      have := mul_nonneg (by linarith : 0 ≤ 1 - lam) hm.le
      linarith
  constructor
  · intro h
    rcases mul_eq_zero.mp h with hc | hmu0
    · rcases mul_eq_zero.mp hc with ha | hb
      · left; linarith
      · right; exact hb
    · exact absurd hmu0 hmu.ne'
  · intro h
    rcases h with h | h
    · rw [h]; ring
    · rw [h]; ring

/-- The defect is nonzero at a moving frame for `(1 - lam) * E ≠ 0` and
`m + lam * E ≠ 0`. -/
theorem slopeMomentum_off_shell {m E lam : ℝ} (hc : (1 - lam) * E ≠ 0)
    (hmu : m + lam * E ≠ 0) (frame : FrameHyperboloid) (hx : frame.1.2 ≠ 0) :
    lorentzQ (slopeMomentum m E lam frame) ≠ (m + E) ^ 2 := by
  rw [lorentzQ_slopeMomentum]
  intro h
  have ht : frame.1.1 - 1 ≠ 0 := by
    intro h1
    have hu : frame.1.1 = 1 := by linarith
    have hframe := frame_time_sq_eq_one_add_spatial frame
    rw [hu] at hframe
    have hs : spatialNormSq frame.1.2 = 0 := by linarith
    exact (spatialNormSq_pos hx).ne' hs
  have : 2 * ((1 - lam) * E) * (m + lam * E) * (frame.1.1 - 1) = 0 := by linarith
  rcases mul_eq_zero.mp this with h2 | h3
  · rcases mul_eq_zero.mp h2 with h4 | h5
    · rcases mul_eq_zero.mp h4 with h6 | h7
      · norm_num at h6
      · exact hc h7
    · exact hmu h5
  · exact ht h3

/-! ## (4) Linearity in `(m, E)` at fixed slope and frame -/

theorem slopeMomentum_add (m₁ E₁ m₂ E₂ lam : ℝ) (frame : FrameHyperboloid) :
    slopeMomentum (m₁ + m₂) (E₁ + E₂) lam frame =
      slopeMomentum m₁ E₁ lam frame + slopeMomentum m₂ E₂ lam frame := by
  apply Prod.ext
  · simp only [Prod.fst_add, slopeMomentum_time]
    ring
  · funext i
    simp only [Prod.snd_add, Pi.add_apply, slopeMomentum_spatial, Pi.smul_apply,
      smul_eq_mul]
    ring

theorem slopeMomentum_smul (s m E lam : ℝ) (frame : FrameHyperboloid) :
    slopeMomentum (s * m) (s * E) lam frame = s • slopeMomentum m E lam frame := by
  apply Prod.ext
  · simp only [Prod.smul_fst, smul_eq_mul, slopeMomentum_time]
    ring
  · funext i
    simp only [Prod.smul_snd, Pi.smul_apply, smul_eq_mul, slopeMomentum_spatial]
    ring

/-- At slope one the sum of two composites at one frame is the shape-A
momentum at the summed rest energy, with invariant mass the summed rest
energy. -/
theorem slopeMomentum_one_add_shapeA (m₁ E₁ m₂ E₂ : ℝ) (frame : FrameHyperboloid) :
    slopeMomentum m₁ E₁ 1 frame + slopeMomentum m₂ E₂ 1 frame =
        fourMomentum ((m₁ + m₂) + (E₁ + E₂)) frame ∧
      lorentzQ (slopeMomentum m₁ E₁ 1 frame + slopeMomentum m₂ E₂ 1 frame) =
        ((m₁ + m₂) + (E₁ + E₂)) ^ 2 := by
  have h : slopeMomentum m₁ E₁ 1 frame + slopeMomentum m₂ E₂ 1 frame =
      fourMomentum ((m₁ + m₂) + (E₁ + E₂)) frame := by
    rw [← slopeMomentum_add, slopeMomentum_one]
  exact ⟨h, by rw [h, lorentzQ_fourMomentum]⟩

/-- For `(1 - lam) * E ≠ 0` and a frame with nonzero spatial part, the
member is not the shape-A momentum of any parameter: only the slope-one
sums are shape-A momenta at moving frames. -/
theorem slopeMomentum_not_fourMomentum {m E lam : ℝ} (hc : (1 - lam) * E ≠ 0)
    (frame : FrameHyperboloid) (hx : frame.1.2 ≠ 0) (μ : ℝ) :
    slopeMomentum m E lam frame ≠ fourMomentum μ frame := by
  intro heq
  have hsp := congrArg Prod.snd heq
  rw [slopeMomentum_spatial, fourMomentum_spatial] at hsp
  have hmu : m + lam * E = μ := by
    have : (m + lam * E - μ) • frame.1.2 = 0 := by
      rw [sub_smul, hsp, sub_self]
    rcases smul_eq_zero.mp this with h | h
    · linarith
    · exact absurd h hx
  have ht := congrArg Prod.fst heq
  rw [slopeMomentum_time, fourMomentum_time, hmu] at ht
  apply hc
  linarith

/-! ## (5) Join to the inertia precursor and the coupled worldline -/

/-- The shape-A composite momentum of `Geometry/InternalEnergyInertia.lean`
is the slope-one member, with no factor. -/
theorem compositeFourMomentum_eq_slopeMomentum_one (m E : ℝ)
    (frame : FrameHyperboloid) :
    compositeFourMomentum m E frame = slopeMomentum m E 1 frame := by
  rw [slopeMomentum_one]
  rfl

/-- The discrete momentum of the quadratic worldline action on a frame
increment is twice the slope-one member at zero internal energy: the factor
`2` is the inertial coefficient `2 m` of that action and is not a slope. -/
theorem discreteMomentum_eq_two_slopeMomentum (m : ℝ) (frame : FrameHyperboloid)
    (x : ℕ → Herm2) (k : ℕ) (hx : x (k + 1) - x k = (frame : Herm2)) :
    OPH.ChargeFixedInteraction.discreteMomentum m x k =
      (2 : ℝ) • slopeMomentum m 0 1 frame := by
  rw [OPH.ChargeFixedInteraction.discreteMomentum_frame m frame x k hx,
    slopeMomentum_one, add_zero]

/-! ## (6) Conclusion, Legendre scope, row register -/

/-- **Frame covariance selects slope one (kinematic).**  Among the declared
slope family, frame covariance of the composite four-momentum together
with the frame-independent scalar ledger `E` holds exactly at slope one or
at zero ledger; every member has rest energy `m + E`; the slope-one member
is the shape-A composite momentum of the inertia precursor; and the shell
identity at the rest energy holds at every frame iff
`(1 - lam) * E * (m + lam * E) = 0`.  The selection
consumes the declared family and the covariance requirement.  It is a
selection of the momentum map, independent of and consistent with the
dynamical shape selection of `Geometry/InternalEnergyInertia.lean`, and it
does not derive a mass-energy identity: the identification of `E` with a
physical internal energy and of `m` with a physical mass is the laboratory
clock and energy calibration import, the identification of the frame with
a physical observer is the physical spacetime attachment row, and the
composite shape is the coupled-action row; none is discharged here. -/
theorem covariance_selects_slope_one :
    (∀ m E lam : ℝ, FrameCovariant (slopeMomentum m E lam) ↔ lam = 1 ∨ E = 0) ∧
    (∀ m E lam : ℝ, (slopeMomentum m E lam standardFrame).1 = m + E) ∧
    (∀ (m E : ℝ) (frame : FrameHyperboloid),
      compositeFourMomentum m E frame = slopeMomentum m E 1 frame) ∧
    (∀ m E lam : ℝ,
      (∀ frame : FrameHyperboloid,
          lorentzQ (slopeMomentum m E lam frame) = (m + E) ^ 2) ↔
        (1 - lam) * E * (m + lam * E) = 0) :=
  ⟨frameCovariant_slopeMomentum_iff, slopeMomentum_standardFrame_time,
    compositeFourMomentum_eq_slopeMomentum_one, slopeMomentum_shell_all_iff⟩

/-- The Legendre non-identifiability at its scope, re-exported from the
inertia precursor: realized histories select no velocity solver and the
regular enrichments at curvatures `1` and `2` differ, so every member of
the slope family is a declared enrichment and the covariance selection
above is kinematic, on the momentum map only. -/
theorem legendre_scope_cited :
    (¬ ∃ vel : ℝ → ℝ → ℝ,
        OPH.Variational.SolvesMomentum OPH.Variational.chainLogLagrangian vel) ∧
      OPH.Variational.chainCurvedLagrangian 1 ≠
        OPH.Variational.chainCurvedLagrangian 2 :=
  legendre_nonidentifiability_cited

/-- The rows and the import this module touches. -/
inductive TouchedRow : Type
  /-- Source clock and duration row. -/
  | sourceClock
  /-- Physical spacetime attachment row. -/
  | spacetimeAttachment
  /-- Coupled-action row. -/
  | coupledAction
  /-- Laboratory clock and energy calibration import. -/
  | calibrationImport

/-- Register label of the rows touched; a label, not a discharge. -/
def touchedRows : List TouchedRow :=
  [TouchedRow.sourceClock, TouchedRow.spacetimeAttachment,
    TouchedRow.coupledAction, TouchedRow.calibrationImport]

/-- The rows this module discharges: none. -/
def dischargedRows : List TouchedRow := []

theorem dischargedRows_empty : dischargedRows = [] := rfl

end

end OPH.CompositeMomentumCovariance

/- Axiom audit: expected at most `propext`, `Classical.choice`, `Quot.sound`
per line; no native decision procedure. -/

#print axioms OPH.CompositeMomentumCovariance.standardBoost_moves_restVector
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_one
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_zero
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_standardFrame
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_standardFrame_time
#print axioms OPH.CompositeMomentumCovariance.frameCovariant_fourMomentum
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_mapFrame
#print axioms OPH.CompositeMomentumCovariance.frameCovariant_slopeMomentum_iff
#print axioms OPH.CompositeMomentumCovariance.lorentzQ_slopeMomentum
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_shell_defect_spatial
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_shell_all_iff
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_shell_all_iff_of_pos
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_off_shell
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_add
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_smul
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_one_add_shapeA
#print axioms OPH.CompositeMomentumCovariance.slopeMomentum_not_fourMomentum
#print axioms OPH.CompositeMomentumCovariance.compositeFourMomentum_eq_slopeMomentum_one
#print axioms OPH.CompositeMomentumCovariance.discreteMomentum_eq_two_slopeMomentum
#print axioms OPH.CompositeMomentumCovariance.covariance_selects_slope_one
#print axioms OPH.CompositeMomentumCovariance.legendre_scope_cited
#print axioms OPH.CompositeMomentumCovariance.dischargedRows_empty
