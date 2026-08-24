import Mathlib

/-!
# The fixed-capacity dark-energy w-law

The de Sitter capacity relation of the corpus reads the dark-energy density
as inversely proportional to the record capacity: `rho_DE(a) = kappa / N(a)`,
with `kappa` the single positive constant packaging
`c^4 * 3 * pi / (8 * pi * G * l_P^2)` (the SI constants are premises and are
not unfolded here).  The effective equation of state is DEFINED through the
FRW continuity relation,

  `w(a) = -1 - (1/3) * (a * rho'(a)) / rho(a)`,

which is `w(a) = -1 - (1/3) d ln rho / d ln a` wherever `rho` is positive.
Per the P4 remark of the DK-01 theorem file
(`proof/epic_wins/dk01_wlaw/DK01_FIXED_N_WLAW.md` in the oph-meta
metarepository), this definition-through-continuity is exactly the object a
`w0waCDM` analyst reconstructs: the fit reads the dark-energy density from
the shape of `H^2(a)` beyond the matter and radiation terms and defines
`w_eff` from that density through continuity, so the theorems below are
stated for exactly the object the comparator measures.

Formalization choice.  The domain is all scale factors `a > 0`, taken as
`Set.Ioi (0 : ℝ)`; this is the cleanest single domain.  Every pointwise
proof below uses only data in a neighborhood of the point, or order data to
the right of it, so hypotheses restricted to any open scored subinterval
give the same conclusions on that subinterval by the same arguments.

Results:

* Fixed-N w-law (`CapacityLaw.w_eq_neg_one_of_constN`): if `N` is constant
  on the domain then `w(a) = -1` at every point.  The constant value is a
  universally quantified variable carrying no hypothesis, so the capacity
  value cancels: the statement takes no input about `N` beyond the
  positivity packaged in the structure.  Zero free parameters.
* Converse drift map, exact (`CapacityLaw.w_eq_drift`):
  `w(a) = -1 + (1/3) * (a * N'(a)) / N(a)` on the domain.  Growing capacity
  gives quintessence-side `w`, shrinking capacity gives phantom-side `w`,
  matching the sign audit of the DK-01 theorem file.
* No-phantom bound (`CapacityLaw.w_ge_neg_one_of_monotone`): if `N` is
  monotone nondecreasing on the domain then `w(a) >= -1` everywhere.
  Capacity does not decrease, so the OPH capacity reading forbids any
  phantom epoch on the drift branch as well; on the fixed branch equality
  holds.  The strict readback
  (`CapacityLaw.deriv_N_neg_of_w_lt_neg_one`): `w(a0) < -1` at a point of
  the domain forces `N'(a0) < 0`, capacity loss at `a0`.  (Every point of
  the open domain is interior.)
* CPL projection (`cpl_at_neg_one_unique`,
  `CapacityLaw.cpl_forced_of_constN`): a CPL form `w0 + wa * (1 - a)` that
  equals `-1` at two distinct scale factors has `(w0, wa) = (-1, 0)`
  exactly; under constant capacity the CPL representation of `w` is forced
  to `(-1, 0)`.  This is the exact bridge to the comparator's parameter
  plane.
* Inhabitants: the constant-capacity instance (`constCapacity`, with
  `w = -1` computed), the power-law family (`powerCapacity` with
  `N(a) = N0 * a ^ eps`, giving `w = -1 + eps / 3` exactly, an exact
  thawing-side line for `eps > 0` and monotone for `eps >= 0`), and the
  decreasing instance `eps = -1` with `w = -4/3 < -1`, which shows the
  monotone hypothesis of the no-phantom bound is load-bearing.

What is not proved here.  P1 (fixed capacity over the scored range) and
monotone capacity are named branch premises of the corpus, not theorems;
this module derives their consequences and nothing selects between the
branches.  The de Sitter capacity relation `Lambda * l_P^2 = 3 * pi / N`
and the constancy of `l_P` are premises (Assumption A-GEOM of
`cosmology/oph_boltzmann_transport_derivation.tex`; principle SL-4 and the
capacity-closure section as cited in the DK-01 theorem file), packaged here
as the single positive constant `kappa`; the FRW continuity relation is the
definition of `w`, not a derivation from field equations.  Nothing in this
module scores data.  DESI DR1 and DR2 are already-examined data, named only
as provenance, and are excluded from any future comparison this surface
feeds; every decision rule lives in the separate registration proposal and
is a proposal pending the owner's freeze.
-/

namespace OPH.EinsteinBranch

open Filter

noncomputable section

/-! ## The capacity dark-energy density and the continuity-defined w -/

/-- Capacity data for the dark-energy sector: a record-capacity function
`N` of the scale factor, positive and differentiable on the domain
`Set.Ioi 0`, together with the single positive constant `kappa` packaging
`c^4 * 3 * pi / (8 * pi * G * l_P^2)`.  The de Sitter capacity relation and
the constancy of `l_P` are the premises that make `kappa` one constant; the
SI constants are not unfolded. -/
structure CapacityLaw where
  /-- Record capacity as a function of the scale factor. -/
  N : ℝ → ℝ
  /-- The bundled constant `c^4 * 3 * pi / (8 * pi * G * l_P^2)`. -/
  kappa : ℝ
  kappa_pos : 0 < kappa
  N_pos : ∀ a ∈ Set.Ioi (0 : ℝ), 0 < N a
  N_diff : ∀ a ∈ Set.Ioi (0 : ℝ), DifferentiableAt ℝ N a

/-- The capacity dark-energy density, `rho_DE(a) = kappa / N(a)`. -/
def CapacityLaw.rho (C : CapacityLaw) (a : ℝ) : ℝ := C.kappa / C.N a

/-- The effective-fluid equation of state, DEFINED through the FRW
continuity relation: `w(a) = -1 - (1/3) * (a * rho'(a)) / rho(a)`.  This is
`-1 - (1/3) d ln rho / d ln a` wherever `rho` is positive, and per the P4
remark of the DK-01 theorem file it is exactly the object a `w0waCDM`
analyst reconstructs from the expansion history. -/
def CapacityLaw.w (C : CapacityLaw) (a : ℝ) : ℝ :=
  -1 - (1 / 3) * (a * deriv C.rho a) / C.rho a

/-- Derivative of the capacity density through the quotient rule:
`rho'(a) = -(kappa * N'(a)) / N(a)^2` on the domain. -/
theorem CapacityLaw.rho_hasDerivAt (C : CapacityLaw) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    HasDerivAt C.rho (-(C.kappa * deriv C.N a) / C.N a ^ 2) a := by
  have hN : HasDerivAt C.N (deriv C.N a) a := (C.N_diff a ha).hasDerivAt
  have hNa : C.N a ≠ 0 := (C.N_pos a ha).ne'
  have h := (hasDerivAt_const a C.kappa).div hN hNa
  have heq : (0 * C.N a - C.kappa * deriv C.N a) / C.N a ^ 2
      = -(C.kappa * deriv C.N a) / C.N a ^ 2 := by ring
  rw [heq] at h
  exact h

/-! ## Theorem: the fixed-N w-law -/

/-- Fixed-N w-law: if the capacity is constant on the domain then
`w(a) = -1` at every point of the domain.  The constant value `c` carries
no hypothesis, so the capacity value cancels: no input about `N` enters
beyond the positivity packaged in the structure.  Zero free parameters. -/
theorem CapacityLaw.w_eq_neg_one_of_constN (C : CapacityLaw) (c : ℝ)
    (hc : ∀ x ∈ Set.Ioi (0 : ℝ), C.N x = c) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    C.w a = -1 := by
  have hd : deriv C.rho a = 0 := by
    have hmem : Set.Ioi (0 : ℝ) ∈ nhds a := isOpen_Ioi.mem_nhds ha
    have heq : C.rho =ᶠ[nhds a] fun _ => C.kappa / c := by
      filter_upwards [hmem] with x hx
      show C.kappa / C.N x = C.kappa / c
      rw [hc x hx]
    rw [heq.deriv_eq]
    exact deriv_const a (C.kappa / c)
  unfold CapacityLaw.w
  rw [hd]
  simp

/-! ## Theorem: the converse drift map, exact -/

/-- Converse drift map, exact: on the domain,
`w(a) = -1 + (1/3) * (a * N'(a)) / N(a)`.  The sign direction matches the
DK-01 sign audit: `rho` proportional to `1 / N` and
`d ln rho / d ln a = -3 (1 + w)` give growing capacity `w > -1`
(quintessence side) and shrinking capacity `w < -1` (phantom side). -/
theorem CapacityLaw.w_eq_drift (C : CapacityLaw) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    C.w a = -1 + (1 / 3) * (a * deriv C.N a) / C.N a := by
  have hd : deriv C.rho a = -(C.kappa * deriv C.N a) / C.N a ^ 2 :=
    (C.rho_hasDerivAt ha).deriv
  have hval : C.rho a = C.kappa / C.N a := rfl
  have hNa : C.N a ≠ 0 := (C.N_pos a ha).ne'
  have hk : C.kappa ≠ 0 := C.kappa_pos.ne'
  unfold CapacityLaw.w
  rw [hd, hval]
  field_simp
  ring

/-! ## Theorem: the no-phantom bound on the monotone-capacity branch -/

/-- A function monotone nondecreasing on `Set.Ioi 0` and differentiable at
an interior point has nonnegative derivative there.  Proved through the
one-sided slope limit; every point of the open domain is interior. -/
theorem deriv_nonneg_of_monotoneOn {f : ℝ → ℝ} {a : ℝ}
    (hf : MonotoneOn f (Set.Ioi (0 : ℝ))) (ha : a ∈ Set.Ioi (0 : ℝ))
    (hdiff : DifferentiableAt ℝ f a) : 0 ≤ deriv f a := by
  have hslope : Tendsto (slope f a) (nhdsWithin a (Set.Ioi a))
      (nhds (deriv f a)) := by
    have h := hasDerivAt_iff_tendsto_slope.mp hdiff.hasDerivAt
    exact h.mono_left (nhdsWithin_mono a fun x hx =>
      Set.mem_compl_singleton_iff.mpr (ne_of_gt (Set.mem_Ioi.mp hx)))
  have hev : ∀ᶠ x in nhdsWithin a (Set.Ioi a), 0 ≤ slope f a x := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    have hax : a < x := Set.mem_Ioi.mp hx
    have hx0 : x ∈ Set.Ioi (0 : ℝ) := lt_trans (Set.mem_Ioi.mp ha) hax
    have hfa : f a ≤ f x := hf ha hx0 (le_of_lt hax)
    rw [slope_def_field]
    exact div_nonneg (sub_nonneg.mpr hfa) (le_of_lt (sub_pos.mpr hax))
  exact ge_of_tendsto hslope hev

/-- No-phantom bound: if the capacity is monotone nondecreasing on the
domain then `w(a) >= -1` at every point.  Capacity does not decrease, so
the OPH capacity reading forbids any phantom epoch on the drift branch as
well; on the fixed branch equality holds
(`CapacityLaw.w_eq_neg_one_of_constN`). -/
theorem CapacityLaw.w_ge_neg_one_of_monotone (C : CapacityLaw)
    (hmono : MonotoneOn C.N (Set.Ioi (0 : ℝ))) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    -1 ≤ C.w a := by
  have hN' : 0 ≤ deriv C.N a :=
    deriv_nonneg_of_monotoneOn hmono ha (C.N_diff a ha)
  have hterm : 0 ≤ (1 / 3) * (a * deriv C.N a) / C.N a :=
    div_nonneg
      (mul_nonneg (by norm_num)
        (mul_nonneg (le_of_lt (Set.mem_Ioi.mp ha)) hN'))
      (le_of_lt (C.N_pos a ha))
  rw [C.w_eq_drift ha]
  linarith

/-- Strict converse readback: `w(a0) < -1` at a point of the domain forces
`N'(a0) < 0`, capacity loss at `a0`.  Every point of the open domain is
interior, so no separate interiority hypothesis is needed. -/
theorem CapacityLaw.deriv_N_neg_of_w_lt_neg_one (C : CapacityLaw) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) (hw : C.w a < -1) : deriv C.N a < 0 := by
  by_contra hcon
  push Not at hcon
  have hterm : 0 ≤ (1 / 3) * (a * deriv C.N a) / C.N a :=
    div_nonneg
      (mul_nonneg (by norm_num)
        (mul_nonneg (le_of_lt (Set.mem_Ioi.mp ha)) hcon))
      (le_of_lt (C.N_pos a ha))
  rw [C.w_eq_drift ha] at hw
  linarith

/-! ## Theorem: the CPL projection -/

/-- A CPL form `w0 + wa * (1 - a)` that equals `-1` at two distinct scale
factors has `(w0, wa) = (-1, 0)` exactly: the functions `1` and `1 - a`
are linearly independent on any two-point set. -/
theorem cpl_at_neg_one_unique {w0 wa a1 a2 : ℝ} (hne : a1 ≠ a2)
    (h1 : w0 + wa * (1 - a1) = -1) (h2 : w0 + wa * (1 - a2) = -1) :
    w0 = -1 ∧ wa = 0 := by
  have hwa : wa * (a2 - a1) = 0 := by linear_combination h1 - h2
  have hwa0 : wa = 0 := by
    rcases mul_eq_zero.mp hwa with h | h
    · exact h
    · exact absurd (sub_eq_zero.mp h) (Ne.symm hne)
  refine ⟨?_, hwa0⟩
  rw [hwa0] at h1
  linarith

/-- Exact bridge to the comparator's parameter plane: if the capacity is
constant on the domain and `w` agrees on the domain with a CPL form
`w0 + wa * (1 - a)`, then `(w0, wa) = (-1, 0)` exactly.  The domain
contains the two distinct points `1` and `2`, which suffice. -/
theorem CapacityLaw.cpl_forced_of_constN (C : CapacityLaw) (c w0 wa : ℝ)
    (hc : ∀ x ∈ Set.Ioi (0 : ℝ), C.N x = c)
    (hcpl : ∀ a ∈ Set.Ioi (0 : ℝ), C.w a = w0 + wa * (1 - a)) :
    w0 = -1 ∧ wa = 0 := by
  have h1mem : (1 : ℝ) ∈ Set.Ioi (0 : ℝ) := Set.mem_Ioi.mpr one_pos
  have h2mem : (2 : ℝ) ∈ Set.Ioi (0 : ℝ) := Set.mem_Ioi.mpr (by norm_num)
  have h1 : w0 + wa * (1 - (1 : ℝ)) = -1 := by
    rw [← hcpl 1 h1mem]
    exact C.w_eq_neg_one_of_constN c hc h1mem
  have h2 : w0 + wa * (1 - (2 : ℝ)) = -1 := by
    rw [← hcpl 2 h2mem]
    exact C.w_eq_neg_one_of_constN c hc h2mem
  exact cpl_at_neg_one_unique (by norm_num) h1 h2

/-! ## Inhabitants

The constant-capacity instance realizes the fixed branch with `w = -1`
computed.  The power-law family realizes the drift branch exactly:
`N(a) = N0 * a ^ eps` gives `w = -1 + eps / 3` at every point, an exact
thawing-side line for `eps > 0`, monotone for `eps >= 0`, and for
`eps = -1` a decreasing capacity with `w = -4/3 < -1`, which shows the
monotone hypothesis of the no-phantom bound is load-bearing. -/

/-- Constant-capacity instance: `N = 1`, `kappa = 1`. -/
def constCapacity : CapacityLaw where
  N := fun _ => 1
  kappa := 1
  kappa_pos := one_pos
  N_pos := fun _ _ => one_pos
  N_diff := fun _ _ => differentiableAt_const 1

instance : Inhabited CapacityLaw := ⟨constCapacity⟩

/-- The constant-capacity instance has `w = -1` at every point of the
domain. -/
theorem constCapacity_w {a : ℝ} (ha : a ∈ Set.Ioi (0 : ℝ)) :
    constCapacity.w a = -1 :=
  constCapacity.w_eq_neg_one_of_constN 1 (fun _ _ => rfl) ha

/-- Power-law capacity instance: `N(a) = N0 * a ^ eps` (real power), with
`N0 > 0` and any real exponent `eps`. -/
def powerCapacity (N0 : ℝ) (hN0 : 0 < N0) (eps : ℝ) : CapacityLaw where
  N := fun a => N0 * a ^ eps
  kappa := 1
  kappa_pos := one_pos
  N_pos := fun _ ha =>
    mul_pos hN0 (Real.rpow_pos_of_pos (Set.mem_Ioi.mp ha) eps)
  N_diff := fun _ ha =>
    ((Real.hasDerivAt_rpow_const (p := eps)
      (Or.inl (ne_of_gt (Set.mem_Ioi.mp ha)))).const_mul N0).differentiableAt

/-- Derivative of the power-law capacity on the domain. -/
theorem powerCapacity_N_deriv (N0 : ℝ) (hN0 : 0 < N0) (eps : ℝ) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    deriv (powerCapacity N0 hN0 eps).N a = N0 * (eps * a ^ (eps - 1)) := by
  have h : HasDerivAt (fun x : ℝ => N0 * x ^ eps)
      (N0 * (eps * a ^ (eps - 1))) a :=
    (Real.hasDerivAt_rpow_const (p := eps)
      (Or.inl (ne_of_gt (Set.mem_Ioi.mp ha)))).const_mul N0
  have hNe : (powerCapacity N0 hN0 eps).N = fun x : ℝ => N0 * x ^ eps := rfl
  rw [hNe]
  exact h.deriv

/-- The power-law capacity gives `w = -1 + eps / 3` exactly, at every point
of the domain: the exact drift line of the DK-01 converse map. -/
theorem powerCapacity_w (N0 : ℝ) (hN0 : 0 < N0) (eps : ℝ) {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    (powerCapacity N0 hN0 eps).w a = -1 + eps / 3 := by
  have ha0 : (0 : ℝ) < a := Set.mem_Ioi.mp ha
  have hNval : (powerCapacity N0 hN0 eps).N a = N0 * a ^ eps := rfl
  rw [(powerCapacity N0 hN0 eps).w_eq_drift ha,
    powerCapacity_N_deriv N0 hN0 eps ha, hNval]
  have hpow : a ^ (eps - 1) * a = a ^ eps := by
    have h1 : a ^ (eps - 1 + 1) = a ^ (eps - 1) * a ^ (1 : ℝ) :=
      Real.rpow_add ha0 (eps - 1) 1
    rw [Real.rpow_one] at h1
    have h2 : eps - 1 + 1 = eps := by ring
    rw [h2] at h1
    exact h1.symm
  have hkey : a * (N0 * (eps * a ^ (eps - 1))) = eps * (N0 * a ^ eps) := by
    calc a * (N0 * (eps * a ^ (eps - 1)))
        = N0 * eps * (a ^ (eps - 1) * a) := by ring
      _ = N0 * eps * a ^ eps := by rw [hpow]
      _ = eps * (N0 * a ^ eps) := by ring
  rw [hkey]
  have hNa : N0 * a ^ eps ≠ 0 :=
    (mul_pos hN0 (Real.rpow_pos_of_pos ha0 eps)).ne'
  field_simp

/-- The power-law capacity is monotone nondecreasing on the domain for
`eps >= 0`, so it inhabits the monotone-branch premise nontrivially. -/
theorem powerCapacity_monotoneOn (N0 : ℝ) (hN0 : 0 < N0) {eps : ℝ}
    (heps : 0 ≤ eps) :
    MonotoneOn (powerCapacity N0 hN0 eps).N (Set.Ioi (0 : ℝ)) := by
  intro x hx _ _ hxy
  have hx0 : (0 : ℝ) < x := Set.mem_Ioi.mp hx
  show N0 * x ^ eps ≤ N0 * _ ^ eps
  exact mul_le_mul_of_nonneg_left
    (Real.rpow_le_rpow (le_of_lt hx0) hxy heps) (le_of_lt hN0)

/-- Growing capacity, thawing side: `eps > 0` places `w` strictly above
`-1`, at the exact value `-1 + eps / 3`. -/
theorem powerCapacity_w_gt_neg_one (N0 : ℝ) (hN0 : 0 < N0) {eps : ℝ}
    (heps : 0 < eps) {a : ℝ} (ha : a ∈ Set.Ioi (0 : ℝ)) :
    -1 < (powerCapacity N0 hN0 eps).w a := by
  rw [powerCapacity_w N0 hN0 eps ha]
  linarith

/-- Decreasing capacity crosses to the phantom side: the instance
`N(a) = a ^ (-1)` has `w = -4/3 < -1` at every point of the domain.  The
monotone hypothesis of the no-phantom bound is therefore load-bearing. -/
theorem decreasingCapacity_w_lt_neg_one {a : ℝ} (ha : a ∈ Set.Ioi (0 : ℝ)) :
    (powerCapacity 1 one_pos (-1)).w a < -1 := by
  rw [powerCapacity_w 1 one_pos (-1) ha]
  norm_num

/-- The strict readback fires on the decreasing instance: its phantom-side
`w` certifies capacity loss, `N' (a) < 0`, at every point of the domain. -/
theorem decreasingCapacity_deriv_N_neg {a : ℝ}
    (ha : a ∈ Set.Ioi (0 : ℝ)) :
    deriv (powerCapacity 1 one_pos (-1)).N a < 0 :=
  (powerCapacity 1 one_pos (-1)).deriv_N_neg_of_w_lt_neg_one ha
    (decreasingCapacity_w_lt_neg_one ha)

/-! ## Per-theorem axiom audit -/

#print axioms CapacityLaw.rho_hasDerivAt
#print axioms CapacityLaw.w_eq_neg_one_of_constN
#print axioms CapacityLaw.w_eq_drift
#print axioms deriv_nonneg_of_monotoneOn
#print axioms CapacityLaw.w_ge_neg_one_of_monotone
#print axioms CapacityLaw.deriv_N_neg_of_w_lt_neg_one
#print axioms cpl_at_neg_one_unique
#print axioms CapacityLaw.cpl_forced_of_constN
#print axioms constCapacity_w
#print axioms powerCapacity_N_deriv
#print axioms powerCapacity_w
#print axioms powerCapacity_monotoneOn
#print axioms powerCapacity_w_gt_neg_one
#print axioms decreasingCapacity_w_lt_neg_one
#print axioms decreasingCapacity_deriv_N_neg

end

end OPH.EinsteinBranch
