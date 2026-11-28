import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall

/-!
# Enclosed mass and pressureless spherical shells

In the marginally bound comoving metric the orthonormal Einstein entries can
be written in terms of the areal radius R, its proper-time velocity V and
acceleration A, and their first shell-label derivatives. The functions below
give these geometric expressions. The shell evolution and its label
derivatives imply the pressure-free tensor and the matched residual at G=1/4.
The mass profile and falling clock/metric are the inputs to this calculation.
-/

namespace OPH.EinsteinBranch.SphericalDust

open Filter
open scoped Topology

noncomputable section

/-- Data on one proper-time slice and one regular comoving shell. R, V, A
mean areal radius, its first and second time derivatives. The derivatives in
the fields below are with respect to the shell label; laws need hold only
near the chosen label, permitting a separate interior and exterior. -/
structure RegularShell where
  label : ℝ
  gain : ℝ
  R : ℝ → ℝ
  V : ℝ → ℝ
  A : ℝ → ℝ
  M : ℝ → ℝ
  radial : ℝ
  velocitySlope : ℝ
  accelSlope : ℝ
  massSlope : ℝ
  radius_pos : 0 < R label
  radial_pos : 0 < radial
  radius_deriv : HasDerivAt R radial label
  velocity_deriv : HasDerivAt V velocitySlope label
  accel_deriv : HasDerivAt A accelSlope label
  mass_deriv : HasDerivAt M massSlope label
  shellEnergy : (fun q => R q * (V q)^2) =ᶠ[𝓝 label] fun q => (2*gain) * M q
  shellForce : (fun q => (R q)^2 * A q) =ᶠ[𝓝 label] fun q => (-gain) * M q

namespace RegularShell

variable (s : RegularShell)

def density : ℝ := s.massSlope/(4*Real.pi*(s.R s.label)^2*s.radial)

/-- The standard orthonormal curvature entries of
ds² = -dt² + Rχ² dχ² + R² dΩ², with equal angular entries. -/
def energyEntry : ℝ :=
  (s.radial*(s.V s.label)^2 + 2*(s.R s.label)*(s.V s.label)*s.velocitySlope) /
    ((s.R s.label)^2*s.radial)

def radialEntry : ℝ := -(2*(s.R s.label)*(s.A s.label)+(s.V s.label)^2)/(s.R s.label)^2

def angularEntry : ℝ :=
  -(s.radial*(s.A s.label)+(s.R s.label)*s.accelSlope+(s.V s.label)*s.velocitySlope) /
    ((s.R s.label)*s.radial)

theorem energy_slope :
    s.radial*(s.V s.label)^2 + 2*(s.R s.label)*(s.V s.label)*s.velocitySlope =
      2*s.gain*s.massSlope := by
  have h := (s.radius_deriv.mul (s.velocity_deriv.fun_pow 2)).unique
    ((s.mass_deriv.const_mul (2*s.gain)).congr_of_eventuallyEq s.shellEnergy)
  norm_num at h
  linarith

theorem force_slope :
    2*(s.R s.label)*s.radial*(s.A s.label)+(s.R s.label)^2*s.accelSlope =
      -s.gain*s.massSlope := by
  have h := ((s.radius_deriv.fun_pow 2).mul s.accel_deriv).unique
    ((s.mass_deriv.const_mul (-s.gain)).congr_of_eventuallyEq s.shellForce)
  norm_num at h
  linarith

/-- This calculation also applies at zero velocity and zero mass. -/
theorem curvature_entries :
    s.energyEntry = 8*Real.pi*s.gain*s.density ∧ s.radialEntry = 0 ∧ s.angularEntry = 0 := by
  have hE : s.R s.label*(s.V s.label)^2 = (2*s.gain)*s.M s.label :=
    s.shellEnergy.eq_of_nhds
  have hF : (s.R s.label)^2*s.A s.label = (-s.gain)*s.M s.label :=
    s.shellForce.eq_of_nhds
  have hrad : 2*(s.R s.label)*(s.A s.label)+(s.V s.label)^2 = 0 := by
    have hz : (s.R s.label)*(2*(s.R s.label)*(s.A s.label)+(s.V s.label)^2) = 0 := by
      linear_combination hE + 2*hF
    exact (mul_eq_zero.mp hz).resolve_left s.radius_pos.ne'
  have hang : s.radial*(s.A s.label)+(s.R s.label)*s.accelSlope+
      (s.V s.label)*s.velocitySlope = 0 := by
    have hz : (s.R s.label)*(s.radial*(s.A s.label)+(s.R s.label)*s.accelSlope+
        (s.V s.label)*s.velocitySlope) = 0 := by
      linear_combination s.force_slope + s.energy_slope/2 - (s.radial/2)*hrad
    exact (mul_eq_zero.mp hz).resolve_left s.radius_pos.ne'
  constructor
  · rw [energyEntry, s.energy_slope, density]
    field_simp [s.radius_pos.ne', s.radial_pos.ne', Real.pi_ne_zero]
    ring
  · simp [radialEntry, angularEntry, hrad, hang]

def metricTensor : Mat 3 := fun i j => if i = j then
  (if i = 0 then s.energyEntry else if i = 1 then s.radialEntry else s.angularEntry) else 0

def matterTensor : Mat 3 := fun i j => if i = j ∧ i = 0 then s.density else 0

/-- Ricci contraction and trace removal in a (-+++) tetrad. tr, ta, ra
and aa denote R₀₁₀₁, R₀₂₀₂=R₀₃₀₃, R₁₂₁₂=R₁₃₁₃ and R₂₃₂₃. -/
def sectionalEinstein (tr ta ra aa : ℝ) (i : Fin 4) : ℝ :=
  let ricTime := tr+2*ta
  let ricRadial := -tr+2*ra
  let ricAngular := -ta+ra+aa
  let scalar := -ricTime+ricRadial+2*ricAngular
  if i = 0 then ricTime+scalar/2 else
    (if i = 1 then ricRadial else ricAngular)-scalar/2

/-- Evaluating this contraction on the sectional curvatures of the
marginally bound metric gives precisely the displayed metric tensor. -/
theorem metric_diagonal_from_sections (i : Fin 4) :
    s.metricTensor i i = sectionalEinstein
      (-s.accelSlope/s.radial) (-s.A s.label/s.R s.label)
      (s.V s.label*s.velocitySlope/(s.R s.label*s.radial))
      ((s.V s.label)^2/(s.R s.label)^2) i := by
  by_cases h₀ : i = 0 <;> by_cases h₁ : i = 1 <;>
    simp [metricTensor, sectionalEinstein, h₀, h₁, energyEntry, radialEntry, angularEntry] <;>
    field_simp [s.radius_pos.ne',s.radial_pos.ne'] <;> ring

/-- At the ordinary gain the package's matched four-tensor residual has every
entry zero. This is the direct shell-law route, with the same normalization
as the independent fixed-volume area/entropy theorem. -/
theorem matched_residual (gain : s.gain = 1/4) (i j : Fin 4) :
    Quantitative.residual s.metricTensor s.matterTensor 0 1 i j = 0 := by
  obtain ⟨he,hr,ha⟩ := s.curvature_entries
  by_cases h : i = j
  · subst j
    by_cases hz : i = 0
    · simp [Quantitative.residual, metricTensor, matterTensor, hz, he, gain]
      ring
    · by_cases h₁ : i = 1 <;>
        simp [Quantitative.residual, metricTensor, matterTensor, hz, h₁, hr, ha]
  · simp [Quantitative.residual, metricTensor, matterTensor, h]

theorem matched_charges (gain : s.gain = 1/4) (frame : Fin 10) :
    quadOf (Quantitative.residual s.metricTensor s.matterTensor 0 1)
      (Quantitative.probes frame) = 0 := by
  simp [quadOf, bilinOf, s.matched_residual gain]

end RegularShell

/-- Differentiate the interior clock y³ = 1 − 3 h t / 2 with y>0.
This derives both signed proper-time derivatives of the executable R=q*y². -/
theorem interior_time_derivatives (q h t : ℝ) (y : ℝ → ℝ) (hy : DifferentiableAt ℝ y t)
    (pos : 0 < y t) (clock : (fun x => y x ^ 3) =ᶠ[𝓝 t] fun x => 1-(3*h/2)*x) :
    HasDerivAt (fun x => q*y x^2) (-q*h/y t) t ∧
    HasDerivAt (fun x => -q*h/y x) (-q*h^2/(2*(y t)^4)) t := by
  have hclock : HasDerivAt (fun x : ℝ => 1-(3*h/2)*x) (-(3*h/2)) t := by
    simpa only [zero_sub] using (hasDerivAt_const t (1:ℝ)).sub (hasDerivAt_const_mul (x := t) (3*h/2))
  have hcube := (hy.hasDerivAt.fun_pow 3).unique (hclock.congr_of_eventuallyEq clock)
  norm_num at hcube
  have derivative : deriv y t = -h/(2*y t^2) := by
    apply (eq_div_iff (by positivity : 2*y t^2 ≠ 0)).mpr
    linarith
  have dy : HasDerivAt y (-h/(2*y t^2)) t := derivative ▸ hy.hasDerivAt
  constructor
  · convert (dy.fun_pow 2).const_mul q using 1
    field_simp [pos.ne']; ring
  · convert (hasDerivAt_const t (-q*h)).div dy pos.ne' using 1
    field_simp [pos.ne']; ring

/-- Positive real cube root in the executable interior clock, before focusing. -/
def clockRoot (h t : ℝ) : ℝ := (1-(3*h/2)*t)^(1/3 : ℝ)

theorem executable_time_derivatives (q h t : ℝ) (before_focus : 0 < 1-(3*h/2)*t) :
    HasDerivAt (fun x => q*(clockRoot h x)^2) (-q*h/clockRoot h t) t ∧
    HasDerivAt (fun x => -q*h/clockRoot h x) (-q*h^2/(2*(clockRoot h t)^4)) t := by
  have base : DifferentiableAt ℝ (fun x : ℝ => 1-(3*h/2)*x) t := by fun_prop
  have pos : 0 < clockRoot h t := Real.rpow_pos_of_pos before_focus _
  apply interior_time_derivatives q h t (clockRoot h) (base.rpow_const (Or.inl before_focus.ne')) pos
  have near : ∀ᶠ x in 𝓝 t, 0 < 1-(3*h/2)*x :=
    base.continuousAt.eventually (Ioi_mem_nhds before_focus)
  filter_upwards [near] with x hx
  change ((1-(3*h/2)*x)^(1/3 : ℝ))^3 = 1-(3*h/2)*x
  rw [← Real.rpow_natCast _ 3, ← Real.rpow_mul hx.le]
  norm_num

/-- For y=(1-3*h*t/2)^(1/3)>0 the interior R=q*y² has V=-q*h/y and
A=-q*h²/(2*y⁴). The resulting energy and acceleration laws have enclosed
mass M0*q³. The collapsing sign belongs to V. -/
theorem homogeneous_shell_laws (q y h G M₀ : ℝ) (hy : 0 < y) (hh : h^2 = 2*G*M₀) :
    (q*y^2)*(-q*h/y)^2 = 2*G*(M₀*q^3) ∧
    (q*y^2)^2*(-q*h^2/(2*y^4)) = -G*(M₀*q^3) := by
  constructor
  · calc
      _ = q^3*h^2 := by field_simp [hy.ne']
      _ = _ := by rw [hh]; ring
  · calc
      _ = -(q^3*h^2)/2 := by field_simp [hy.ne']
      _ = _ := by rw [hh]; ring

/-- A fixed-time interior slice of the homogeneous source. Its input is the
shell clock and the enclosed mass profile, including their measured gain. -/
def homogeneousShell (q y h G M₀ : ℝ) (hq : 0 < q) (hy : 0 < y)
    (hh : h^2 = 2*G*M₀) : RegularShell where
  label := q
  gain := G
  R := fun x => x*y^2
  V := fun x => x*(-h/y)
  A := fun x => x*(-h^2/(2*y^4))
  M := fun x => M₀*x^3
  radial := y^2
  velocitySlope := -h/y
  accelSlope := -h^2/(2*y^4)
  massSlope := 3*M₀*q^2
  radius_pos := mul_pos hq (pow_pos hy 2)
  radial_pos := pow_pos hy 2
  radius_deriv := by simpa using (hasDerivAt_id q).mul_const (y^2)
  velocity_deriv := by simpa using (hasDerivAt_id q).mul_const (-h/y)
  accel_deriv := by simpa using (hasDerivAt_id q).mul_const (-h^2/(2*y^4))
  mass_deriv := by
    convert ((hasDerivAt_id q).fun_pow 3).const_mul M₀ using 1
    simp; ring
  shellEnergy := Filter.Eventually.of_forall fun x => by
    convert (homogeneous_shell_laws x y h G M₀ hy hh).1 using 1
    ring
  shellForce := Filter.Eventually.of_forall fun x => by
    convert (homogeneous_shell_laws x y h G M₀ hy hh).2 using 1
    ring

/-- With the ordinary G=1/4, this actual interior source inhabits the
matched residual on all tensor components at each regular shell. -/
theorem homogeneous_residual (q y h M₀ : ℝ) (hq : 0 < q) (hy : 0 < y)
    (hh : h^2 = 2*(1/4)*M₀) (i j : Fin 4) :
    let s := homogeneousShell q y h (1/4) M₀ hq hy hh
    Quantitative.residual s.metricTensor s.matterTensor 0 1 i j = 0 :=
  (homogeneousShell q y h (1/4) M₀ hq hy hh).matched_residual rfl i j

/-- Along a spacelike launch, the FLRW connection gives dγ/ds=−bγ with
b=H dt/ds, and da/ds=ba. Their product is conserved by parallel transport. -/
theorem comoving_boost_derivative (a γ : ℝ → ℝ) (s b : ℝ)
    (ha : HasDerivAt a (b*a s) s) (hγ : HasDerivAt γ (-b*γ s) s) :
    HasDerivAt (fun x => a x*γ x) 0 s := by
  convert ha.mul hγ using 1
  ring

theorem transported_dust_stress (a density gamma : ℝ) (ha : 0 < a) :
    (density/a^3)*(gamma/a)^2 = density*gamma^2/a^5 := by
  field_simp [ha.ne']

/-- A constant vacuum offset and a comoving pressureless source have unique
coefficients in two unit timelike contractions with distinct squared boosts.
The contractions and the dust form of the stress are physical premises. -/
theorem two_observer_dust_vacuum (x y m n l k : ℝ) (hxy : x ≠ y)
    (hx : m*x + l = n*x + k) (hy : m*y + l = n*y + k) :
    m = n ∧ l = k := by
  have h : (m-n)*(x-y) = 0 := by linear_combination hx - hy
  have hmn := sub_eq_zero.mp ((mul_eq_zero.mp h).resolve_right (sub_ne_zero.mpr hxy))
  exact ⟨hmn, by simpa [hmn] using hx⟩

/-- In the one-unit coincidence mechanism every extra mass ticket requires
an ordinary traveller ticket. Any common nonnegative aggregation retains
the ceiling, even across regions with different posting frequencies. -/
theorem coincidence_mass_le_ordinary {ι : Type*} (sites : Finset ι)
    (weight left right : ι → ℝ) (hw : ∀ i ∈ sites, 0 ≤ weight i)
    (hl : ∀ i ∈ sites, 0 ≤ left i) (hr : ∀ i ∈ sites, right i ≤ 1) :
    ∑ i ∈ sites, weight i * left i * right i ≤ ∑ i ∈ sites, weight i * left i := by
  apply Finset.sum_le_sum
  intro i hi
  exact mul_le_of_le_one_right (mul_nonneg (hw i hi) (hl i hi)) (hr i hi)

assert_no_sorry coincidence_mass_le_ordinary
assert_no_sorry two_observer_dust_vacuum
assert_no_sorry RegularShell.energy_slope
assert_no_sorry RegularShell.force_slope
assert_no_sorry RegularShell.curvature_entries
assert_no_sorry RegularShell.metric_diagonal_from_sections
assert_no_sorry RegularShell.matched_residual
assert_no_sorry RegularShell.matched_charges
assert_no_sorry interior_time_derivatives
assert_no_sorry executable_time_derivatives
assert_no_sorry homogeneous_shell_laws
assert_no_sorry homogeneous_residual
assert_no_sorry comoving_boost_derivative
assert_no_sorry transported_dust_stress
#print axioms RegularShell.curvature_entries
#print axioms RegularShell.matched_residual
#print axioms RegularShell.matched_charges
#print axioms homogeneous_shell_laws
#print axioms homogeneous_residual

end
end OPH.EinsteinBranch.SphericalDust
