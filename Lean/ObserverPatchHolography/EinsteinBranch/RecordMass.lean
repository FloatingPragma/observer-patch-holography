import ObserverPatchHolography.EinsteinBranch.Consensus

/-! Identification of masses in the chart reconstructed from a declared
record menu. Exact record equality is used for descent. Experimental
confidence regions and differentiability of physical branches are separate. -/

namespace OPH.EinsteinBranch.RecordMass

/-- A scalar has a well-defined value on the physical quotient exactly
when changing presentation leaves it unchanged. -/
theorem descent {P Q V : Type*} (post : P → Q) (onto : Function.Surjective post)
    (read : P → V) :
    (∃ f : Q → V, ∀ p, f (post p) = read p) ↔
    ∀ p p', post p = post p' → read p = read p' := by
  classical
  constructor
  · rintro ⟨f, hf⟩ p p' h
    rw [← hf p, ← hf p', h]
  · intro stable
    refine ⟨fun q => read (Classical.choose (onto q)), ?_⟩
    intro p
    exact stable _ _ (Classical.choose_spec (onto (post p)))

/-- Removing a common sector offset leaves every within-sector energy
difference invariant. A comparison to a separately specified reference
does change if it does not undergo the offset. -/
theorem reference_shift (a b ref shift : ℝ) :
    (a+shift)-(b+shift)=a-b ∧ (a+shift)-ref=(a-ref)+shift := by
  constructor <;> ring

/-- For a local symmetric branch, any nonzero momentum on the mass shell
determines its rest energy from a measured energy increment. The statement
tests consistency with independently referenced rest spectroscopy. -/
theorem mass_shell_increment (e0 de p c : ℝ) (hd : de ≠ 0)
    (shell : (e0+de)^2=e0^2+c^2*p^2) :
    e0=(c^2*p^2-de^2)/(2*de) := by
  apply (eq_div_iff (mul_ne_zero (by norm_num) hd)).2
  nlinarith

/-- Squared dimensionless lapse-transfer kernels have a strictly ordered
ratio at a larger lapse reading v than at u. This is the pointwise step
in the injectivity proof for the three-station inverse. -/
theorem radial_kernel_order (x y u v : ℝ) (hxy : x < y) (huv : u < v) :
    (1-x*u)*(1-y*v) < (1-x*v)*(1-y*u) := by
  nlinarith [mul_pos (sub_pos.mpr hxy) (sub_pos.mpr huv)]

/-- Free fall determines the ratio of passive to inertial mass once the
field has an independent calibration. Rescaling both masses together
preserves the acceleration on that body's entire set of weak-field probes. -/
theorem passive_scale (mp mi k : ℝ) (hi : mi ≠ 0) (hk : k ≠ 0) :
    (k*mp)/(k*mi)=mp/mi := by field_simp

end OPH.EinsteinBranch.RecordMass
