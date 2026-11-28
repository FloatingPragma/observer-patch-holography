import ObserverPatchHolography.EinsteinBranch.TrailAperture

/-! Joint clock and flight comparisons on the finite routed source family.
The equipped endpoint inequalities below are hypotheses supplied by simultaneous
count bounds. The bare endpoint inequalities use verified integer propagation. -/

namespace OPH.EinsteinBranch.JointClockFlight

noncomputable def response (u f g : ℝ) : ℝ :=
  ((4*u-1)/3)*f + (1-(4*u-1)/3)*g

/-- For every retained record a nonnegative linear statistic transfers upper
bounds on the two physical endpoint preparations to the entire source segment. -/
theorem upper (u f g a b : ℝ) (lo : 1/4 ≤ u) (hi : u ≤ 1)
    (hf : f ≤ a) (hg : g ≤ b) : response u f g ≤ response u a b := by
  unfold response
  exact add_le_add (mul_le_mul_of_nonneg_left hf (by linarith))
    (mul_le_mul_of_nonneg_left hg (by linarith))

/-- A measured source agreement forces inward mean readings in this finite band.
Each pair of equipped bounds represents one mirror launch, normalized by six hops. -/
theorem inward (u f g a b c d : ℝ) (lo : 0.67 ≤ u) (hi : u ≤ 0.69)
    (hf : f ≤ -0.134266648) (hg : g ≤ 0.170790329)
    (ha : a ≤ 0.003235) (hb : b ≤ -0.004699)
    (hc : c ≤ -0.000420) (hd : d ≤ -0.001048) :
    response u f g < 0 ∧ response u a b < 0 ∧ response u c d < 0 := by
  have h1 : (1:ℝ)/4 ≤ u := by linarith
  have h2 : u ≤ 1 := by linarith
  refine ⟨(upper u _ _ _ _ h1 h2 hf hg).trans_lt ?_,
    (upper u _ _ _ _ h1 h2 ha hb).trans_lt ?_,
    (upper u _ _ _ _ h1 h2 hc hd).trans_lt ?_⟩ <;> unfold response <;> linarith

/-- Any nonnegative assignment of an effective source strength must have the
sign of the load-standard witness. This includes all positive dilutions and
any dependence of the strength on additional measured source information. -/
theorem excludes_positive_load (u f g standard amplitude : ℝ)
    (lo : 1/4 ≤ u) (hi : u ≤ 1)
    (hf : f ≤ -0.00121) (hg : g ≤ -0.000137)
    (hs : 0 ≤ standard) (ha : 0 ≤ amplitude) :
    response u f g < amplitude*standard := by
  have h := upper u _ _ _ _ lo hi hf hg
  have hn : response u (-0.00121) (-0.000137) < 0 := by
    unfold response
    linarith
  exact lt_of_lt_of_le (h.trans_lt hn) (mul_nonneg ha hs)

/-- The late two-clock score excludes every nonnegative incoherent source
profile, including zero and arbitrary overall strength. The numerical
hypotheses are the outward integer propagation bounds for this apparatus. -/
theorem excludes_any_port_profile (u f g : ℝ) (load weight : Fin 24 → ℝ)
    (lo : 1/4 ≤ u) (hi : u ≤ 1)
    (hf : 0.02453 ≤ f) (hg : 0.02453 ≤ g)
    (hload : ∀ i, load i ≤ -0.00000248) (hw : ∀ i, 0 ≤ weight i) :
    ∑ i, weight i * load i < response u f g := by
  have hm : ∑ i, weight i * load i ≤ 0 := Finset.sum_nonpos
    (fun i _ => mul_nonpos_of_nonneg_of_nonpos (hw i) (by linarith [hload i]))
  have ha : 0 ≤ (4*u-1)/3 := by linarith
  have hb : 0 ≤ 1-(4*u-1)/3 := by linarith
  unfold response
  nlinarith [mul_nonneg ha (sub_nonneg.mpr hf), mul_nonneg hb (sub_nonneg.mpr hg)]

/-- The generous galaxy inventory and the cosmological paid ratio have a
fixed difference for every universal baryon-weighted mixture of these stocks. -/
theorem mixed_inventory (a b c : ℝ) (h : a+b+c=1) :
    (a*476+b*1088+c*2176)/202 = (a*425+b*1037+c*2125)/202 + 51/202 := by
  linarith

end OPH.EinsteinBranch.JointClockFlight
