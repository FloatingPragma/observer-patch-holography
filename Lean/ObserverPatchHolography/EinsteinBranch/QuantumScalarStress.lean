import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# Directional scalar stress, constraint propagation and diagonal volume focusing

These local derivative identities audit the moment and metric equations used by
the homogeneous quantum scalar benchmark. The expansion rates and stress are
expressed in the same proper time. The gravity equation appears as an explicit
hypothesis of constraint propagation. The volume identity uses geodesic
deviation only in identifying the trace of measured axis accelerations as Ricci.
-/

namespace OPH.EinsteinBranch.QuantumScalarStress

noncomputable section

/-- Cell-averaged energy density of one standing mode. -/
def energyDensity (m b U K : ℝ) : ℝ := (K + (m ^ 2 + b) * U) / 2

/-- Pressure perpendicular to the gradient axis of the standing mode. -/
def transversePressure (m b U K : ℝ) : ℝ := (K - (m ^ 2 + b) * U) / 2

/-- The aligned pressure includes an additional gradient contribution. -/
def alignedPressure (m b U K : ℝ) : ℝ := transversePressure m b U K + b * U

/-- The oscillator moment equations with redshifting squared wave number
imply directional covariant conservation of the mode's mean stress. The
first expansion rate is along its gradient. -/
theorem standingMode_stressContinuity
    (U K b : ℝ → ℝ) (t m C hx hy hz : ℝ)
    (hdU : HasDerivAt U (2 * C) t)
    (hdK : HasDerivAt K (-2 * (hx + hy + hz) * K t - 2 * (m ^ 2 + b t) * C) t)
    (hdb : HasDerivAt b (-2 * hx * b t) t) :
    HasDerivAt (fun s => energyDensity m (b s) (U s) (K s))
      (-(hx * (energyDensity m (b t) (U t) (K t) + alignedPressure m (b t) (U t) (K t)) +
         hy * (energyDensity m (b t) (U t) (K t) + transversePressure m (b t) (U t) (K t)) +
         hz * (energyDensity m (b t) (U t) (K t) + transversePressure m (b t) (U t) (K t)))) t := by
  apply ((hdK.add (((hasDerivAt_const t (m ^ 2)).add hdb).mul hdU)).div_const 2).congr_deriv
  dsimp [energyDensity, alignedPressure, transversePressure]
  ring

/-- The trace-reversed Einstein evolution of one directional expansion rate. -/
def expansionSlope (θ Λ κ H p ρ pressureSum : ℝ) : ℝ :=
  -θ * H + Λ + κ * (p + (ρ - pressureSum) / 2)

/-- The time-time Einstein constraint in proper time, with the source on the right. -/
def hamiltonianDefect (Λ κ hx hy hz ρ : ℝ) : ℝ :=
  hx * hy + hx * hz + hy * hz - Λ - κ * ρ

/-- The directional Einstein equations and the matter continuity equation give
a homogeneous first-order equation for the Hamiltonian-constraint defect. -/
theorem bianchiI_constraintPropagation
    (Hx Hy Hz ρ : ℝ → ℝ) (t Λ κ px py pz : ℝ)
    (hdx : HasDerivAt Hx
      (expansionSlope (Hx t + Hy t + Hz t) Λ κ (Hx t) px (ρ t) (px + py + pz)) t)
    (hdy : HasDerivAt Hy
      (expansionSlope (Hx t + Hy t + Hz t) Λ κ (Hy t) py (ρ t) (px + py + pz)) t)
    (hdz : HasDerivAt Hz
      (expansionSlope (Hx t + Hy t + Hz t) Λ κ (Hz t) pz (ρ t) (px + py + pz)) t)
    (hdρ : HasDerivAt ρ
      (-(Hx t * (ρ t + px) + Hy t * (ρ t + py) + Hz t * (ρ t + pz))) t) :
    HasDerivAt (fun s => hamiltonianDefect Λ κ (Hx s) (Hy s) (Hz s) (ρ s))
      (-2 * (Hx t + Hy t + Hz t) * hamiltonianDefect Λ κ (Hx t) (Hy t) (Hz t) (ρ t)) t := by
  apply (((((hdx.mul hdy).add (hdx.mul hdz)).add (hdy.mul hdz)).sub_const Λ).sub
    (hdρ.const_mul κ)).congr_deriv
  dsimp [expansionSlope, hamiltonianDefect]
  ring


/-- Physical volume in a diagonal orthogonal Jacobi frame with unit initial
coordinate increments. The three functions give the proper lengths. -/
def boxVolume (a b c : ℝ → ℝ) (t : ℝ) : ℝ := a t * b t * c t

/-- Given the diagonal geodesic-deviation interpretation, the timelike Ricci
contraction is the negative trace of fractional Jacobi-axis accelerations. -/
def jacobiRicci (a b c ax ay az : ℝ) : ℝ := -(ax / a + ay / b + az / c)

/-- The acceleration of a diagonal volume contains a kinematic cross term
from its three measured first-order expansion rates. -/
theorem diagonalBox_volumeAcceleration
    (a b c da db dc : ℝ → ℝ) (t ax ay az : ℝ)
    (ha : ∀ s, HasDerivAt a (da s) s) (hb : ∀ s, HasDerivAt b (db s) s)
    (hc : ∀ s, HasDerivAt c (dc s) s)
    (hda : HasDerivAt da ax t) (hdb : HasDerivAt db ay t) (hdc : HasDerivAt dc az t) :
    HasDerivAt (deriv (boxVolume a b c))
      (ax * b t * c t + a t * ay * c t + a t * b t * az +
       2 * (da t * db t * c t + da t * b t * dc t + a t * db t * dc t)) t := by
  have first : deriv (boxVolume a b c) =
      fun s => da s * b s * c s + a s * db s * c s + a s * b s * dc s := by
    funext s
    apply HasDerivAt.deriv
    apply (((ha s).mul (hb s)).mul (hc s)).congr_deriv
    dsimp
    ring
  rw [first]
  apply ((((hda.mul (hb t)).mul (hc t)).add (((ha t).mul hdb).mul (hc t))).add
    (((ha t).mul (hb t)).mul hdc)).congr_deriv
  dsimp
  ring

/-- Recover the timelike Ricci contraction from the observed proper volume
acceleration and directional log-rates. The Jacobi interpretation supplies
the curvature meaning of the left side. The right side uses the volume
acceleration and first derivatives of the individual arms. -/
theorem diagonalBox_timelikeRicci
    (a b c da db dc : ℝ → ℝ) (t ax ay az : ℝ)
    (ha : ∀ s, HasDerivAt a (da s) s) (hb : ∀ s, HasDerivAt b (db s) s)
    (hc : ∀ s, HasDerivAt c (dc s) s)
    (hda : HasDerivAt da ax t) (hdb : HasDerivAt db ay t) (hdc : HasDerivAt dc az t)
    (a_ne : a t ≠ 0) (b_ne : b t ≠ 0) (c_ne : c t ≠ 0) :
    jacobiRicci (a t) (b t) (c t) ax ay az =
      -deriv (deriv (boxVolume a b c)) t / boxVolume a b c t +
      2 * (da t / a t * (db t / b t) + da t / a t * (dc t / c t) +
           db t / b t * (dc t / c t)) := by
  rw [(diagonalBox_volumeAcceleration a b c da db dc t ax ay az ha hb hc hda hdb hdc).deriv]
  dsimp [boxVolume, jacobiRicci]
  field_simp [a_ne, b_ne, c_ne]
  ring

end

end OPH.EinsteinBranch.QuantumScalarStress
