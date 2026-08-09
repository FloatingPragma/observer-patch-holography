import InformationProjection.LogTransitionAction
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Real.Sqrt
import Mathlib.Data.Fin.VecNotation

/-!
# The history-action binding (E9 sector 2, issue #716)

`InformationProjection.LogTransitionAction` rewrites the history law
of a supplied Markov kernel as an exponential tilt.  Here, however,
the common object is defined in the opposite direction: for a chosen
step cost `q`, this module **constructs** the translation-invariant
kernel `gibbsKernel q` and then verifies that its log-transition action
agrees with the increment action of `q`, up to a path-independent
constant.  Thus this is a consistency/binding theorem for the kernel
constructed from `q`; it does not independently derive `q` or that
kernel from source dynamics.  Calling the resulting path law a source
law requires a separate identification of the constructed kernel with
a physically supplied source process.

## Layer 1: the general binding theorem

For any function `q` on a finite nonempty additive abelian group, the
Gibbs kernel `gibbsKernel q x y = exp (-q (y - x)) / gibbsNorm q` is
strictly positive (`gibbsKernel_pos`) and row-stochastic with the same
normalizer in every row (`gibbsKernel_row_sum`); the row sum reindexes
by group translation, and this translation invariance carries the
whole binding.  The keystone `logTransition_eq_increment_add_const`
computes the log-transition action of this constructed Gibbs kernel: it
equals the increment action
`incrementAction q n γ = ∑ i, q (γ (i + 1) - γ i)` plus the per-step
constant `n * log (gibbsNorm q)`.  Consequences:

* `increment_action_reproduces_law`: the tilt of the B7 reference by
  the increment action at multiplier one is the Markov path law of the
  Gibbs kernel; additive constants leave the tilt unchanged.
* `binding_unique_up_to_gauge`: an action-multiplier pair reproduces
  that path law over the reference precisely when its
  multiplier-weighted action is the increment action plus one
  constant.
* `same_start_most_probable_iff_least_increment`: for two paths with a
  common start, path-law order is reversed increment-action order.
  The most probable histories from a common start are the least-action
  histories for the same `q` used to construct the kernel.
* `coupling_scale_is_multiplier_gauge`: the path law of the kernel
  with rescaled cost `g * q` is the tilt of the reference by the
  unrescaled increment action at multiplier `g`.  This is one common
  rescaling freedom for the single cost being considered.  It can
  absorb an overall scale, or the scale of an isolated factor, but it
  does not determine independent relative couplings among factors.

## Layer 2: one certified three-dimensional P-factor increment

The concrete witness transcribes literals from the certified kinetic record
`code/e9_kinetic/kinetic_form_v1.json`, block
`families.P.derived_algebra.kinetic_form_matrix_upper_entries`, gauge
block indices `[0, 3)`: exactly one three-dimensional `three_plus`
factor of the P-family derived algebra `so(3) ⊕ so(3)`.  It is not a
witness for the full six-dimensional P derived algebra.  Scalars encode as
`[i, j, n1, d1, n2, d2]` meaning `n1/d1 + (n2/d2) * sqrt 5`
(`conventions.field`, `conventions.matrix_encoding`), and every
recorded kinetic matrix equals minus one quarter of the corresponding
Killing matrix on the gauge directions (`conventions.kinetic_shape`).
The block entries are `K[0][0] = K[1][1] = K[2][2] = 1200 - 400 √5`,
`K[0][1] = 800 - 400 √5`, `K[0][2] = K[1][2] = -800 + 400 √5`,
symmetric, held here as the integer pair matrices `sectorKineticRat`
and `sectorKineticIrr` with
`K = sectorKineticRat + sectorKineticIrr * sqrt 5`.  The step group is
`(ZMod 3)³` with the centered lift `0 ↦ 0`, `1 ↦ 1`, `2 ↦ -1`: the
certified basis increment lattice truncated to single steps `-1, 0,
+1` per direction, an exact finite carrier.  `sectorIncrement v` is
the kinetic increment `(lift v)ᵀ K (lift v)`.  Receipts:
`sectorIncrement_zero` (value `0` at the zero step),
`sectorIncrement_pos` (strictly positive on all 26 nonzero steps,
evaluated exactly on this three-dimensional block), and
the exact values `sectorIncrement_unit_step` (`1200 - 400 √5` at
`(1, 0, 0)`) and `sectorIncrement_double_step` (`4000 - 1600 √5` at
`(1, 1, 0)`).  Every layer-1 theorem is instantiated at
`sectorIncrement` as the `sector_*` theorems.

## Layer 3: negative control

`dependentKernel` on `ZMod 2` is strictly positive and row-stochastic
with state-dependent rows `(1/3, 2/3)` and `(1/4, 3/4)`.  Its derived
log-transition action lies outside the additive-constant gauge orbit
of every increment action (`binding_requires_invariance`, and
`binding_requires_invariance_multiplier` at every multiplier): the
constant paths at `0` and at `1` carry the same increment and
different transition weights.  Translation invariance is the
load-bearing restriction on kernels represented by this
increment-only construction.

## Claim boundary

The binding holds at representation level over a declared finite step
group with the centered lift and a declared initial law.  The
reference measure is the declared object inherited from the B7 packet,
the initial law times the uniform-step counting weight.  The step
group is a truncation of the certified basis increment lattice to
single steps.  No continuum limit, no physical units, and no
laboratory gauge field are claimed.  The concrete kinetic witness
covers only one three-dimensional P factor.  Its one multiplier gives
only a common overall rescaling (or an isolated-factor rescaling).
The companion `QFT.TwoFactorHistoryBinding` binds both P factors and one
complete color-bearing F family for kernels constructed from the displayed
costs, and proves that only common scaling is multiplier gauge.  Identifying
either constructed kernel with an independent OPH source law, selecting the
relative coefficient, extending to the mirror G family, and supplying
physical units or a continuum limit remain open.
-/

namespace OPH.QFT

open OPH.InformationProjection

/-! ## Layer 1: the general binding theorem -/

/-- The Gibbs normalizer of a step cost `q`: the total exponential
weight of one step over the increment group. -/
noncomputable def gibbsNorm {Ω : Type*} [Fintype Ω] (q : Ω → ℝ) : ℝ :=
  ∑ u, Real.exp (-q u)

theorem gibbsNorm_pos {Ω : Type*} [Fintype Ω] [Nonempty Ω]
    (q : Ω → ℝ) : 0 < gibbsNorm q :=
  Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty

/-- The translation-invariant Gibbs step kernel **constructed from** a
chosen step cost `q` on a finite additive abelian group: the transition
weight from `x` to `y` depends on the increment `y - x` alone.  This
definition is not an independent derivation of a source kernel. -/
noncomputable def gibbsKernel {Ω : Type*} [AddCommGroup Ω] [Fintype Ω]
    (q : Ω → ℝ) (x y : Ω) : ℝ :=
  Real.exp (-q (y - x)) / gibbsNorm q

/-- The Gibbs kernel is strictly positive, for every step cost. -/
theorem gibbsKernel_pos {Ω : Type*} [AddCommGroup Ω] [Fintype Ω]
    [Nonempty Ω] (q : Ω → ℝ) (x y : Ω) : 0 < gibbsKernel q x y :=
  div_pos (Real.exp_pos _) (gibbsNorm_pos q)

/-- **Row stochasticity with constant normalizer.**  The row sum of
the Gibbs kernel is `1` at every state: the sum over targets reindexes
by the group translation `y ↦ y - x`, so the normalizer is one and the
same constant in every row.  This translation invariance is the
load-bearing input of the binding theorems below. -/
theorem gibbsKernel_row_sum {Ω : Type*} [AddCommGroup Ω] [Fintype Ω]
    [Nonempty Ω] (q : Ω → ℝ) (x : Ω) : ∑ y, gibbsKernel q x y = 1 := by
  unfold gibbsKernel
  rw [← Finset.sum_div,
    Fintype.sum_equiv (Equiv.subRight x)
      (fun y => Real.exp (-q (y - x))) (fun u => Real.exp (-q u))
      (fun y => by rw [Equiv.subRight_apply])]
  exact div_self (ne_of_gt (gibbsNorm_pos q))

/-- The increment action of a step cost `q`: the sum of the cost of
the successive increments along the path.  This is the discrete
field-configuration kinetic action once `q` is a quadratic form on the
increment group. -/
noncomputable def incrementAction {Ω : Type*} [AddCommGroup Ω]
    (q : Ω → ℝ) (n : ℕ) (γ : PathSpace Ω n) : ℝ :=
  ∑ i : Fin n, q (γ i.succ - γ i.castSucc)

/-- The increment action is linear in the step cost: rescaling the
cost rescales the action. -/
theorem increment_action_mul {Ω : Type*} [AddCommGroup Ω] (g : ℝ)
    (q : Ω → ℝ) (n : ℕ) (γ : PathSpace Ω n) :
    incrementAction (fun v => g * q v) n γ
      = g * incrementAction q n γ := by
  unfold incrementAction
  rw [Finset.mul_sum]

/-- Log of one Gibbs transition weight: minus the increment cost minus
the log normalizer. -/
theorem gibbsKernel_log {Ω : Type*} [AddCommGroup Ω] [Fintype Ω]
    [Nonempty Ω] (q : Ω → ℝ) (x y : Ω) :
    Real.log (gibbsKernel q x y)
      = -q (y - x) - Real.log (gibbsNorm q) := by
  unfold gibbsKernel
  rw [Real.log_div (Real.exp_ne_zero _) (ne_of_gt (gibbsNorm_pos q)),
    Real.log_exp]

/-- **Keystone for the constructed kernel.**  The log-transition
action of the explicitly defined kernel `gibbsKernel q` equals the
increment action of the chosen `q` plus the per-step constant
`n * log (gibbsNorm q)`.  This verifies consistency of the two
representations; it does not infer `q` or the kernel from an
independently supplied source. -/
theorem logTransition_eq_increment_add_const {Ω : Type*}
    [AddCommGroup Ω] [Fintype Ω] [Nonempty Ω] (q : Ω → ℝ) (n : ℕ)
    (γ : PathSpace Ω n) :
    logTransitionAction (gibbsKernel q) n γ
      = incrementAction q n γ + n * Real.log (gibbsNorm q) := by
  unfold logTransitionAction incrementAction
  rw [Finset.sum_congr rfl fun i _ =>
      gibbsKernel_log q (γ i.castSucc) (γ i.succ),
    Finset.sum_sub_distrib, Finset.sum_neg_distrib, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  ring

/-- Closed form of the Gibbs path law: initial mass times the
exponential of minus the increment action, over the `n`-step
normalizer. -/
theorem markovPathLaw_gibbs_closed {Ω : Type*} [AddCommGroup Ω]
    [Fintype Ω] (pi q : Ω → ℝ) (n : ℕ) (γ : PathSpace Ω n) :
    markovPathLaw pi (gibbsKernel q) n γ
      = pi (γ 0) * Real.exp (-incrementAction q n γ)
        / gibbsNorm q ^ n := by
  unfold markovPathLaw incrementAction gibbsKernel
  rw [Finset.prod_div_distrib, Finset.prod_const, Finset.card_univ,
    Fintype.card_fin, ← Real.exp_sum, Finset.sum_neg_distrib]
  ring

/-- **The increment action reproduces the constructed-kernel law.**
The tilt of the declared B7 reference by the increment action at
multiplier one is the Markov path law of `gibbsKernel q`, which was
constructed from the same `q`: the per-step normalizer constant
aggregates to a path constant and leaves the tilt unchanged.  No
independent source-selection statement is made. -/
theorem increment_action_reproduces_law {Ω : Type*} [AddCommGroup Ω]
    [Fintype Ω] [Nonempty Ω] (pi q : Ω → ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction q n) 1
      = markovPathLaw pi (gibbsKernel q) n := by
  rw [markov_path_law_eq_gibbs pi (gibbsKernel q) n (gibbsKernel_pos q)
      hpi1 (gibbsKernel_row_sum q),
    show logTransitionAction (gibbsKernel q) n
        = fun γ => incrementAction q n γ
          + ↑n * Real.log (gibbsNorm q) from
      funext fun γ => logTransition_eq_increment_add_const q n γ,
    tilt_action_add_const (stepUniformRef pi n) (incrementAction q n) 1
      (↑n * Real.log (gibbsNorm q)) (stepUniformRef_pos pi n hpi)]

/-- **The binding is unique up to gauge.**  An action-multiplier pair
reproduces the Gibbs-kernel path law over the declared reference
precisely when its multiplier-weighted action equals the increment
action plus one additive constant.  This composes
`action_unique_up_to_gauge` with the keystone and folds the two
constants. -/
theorem binding_unique_up_to_gauge {Ω : Type*} [AddCommGroup Ω]
    [Fintype Ω] [Nonempty Ω] (pi q : Ω → ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1)
    (S' : PathSpace Ω n → ℝ) (lam' : ℝ) :
    tilt (stepUniformRef pi n) S' lam'
        = markovPathLaw pi (gibbsKernel q) n
      ↔ ∃ c : ℝ, ∀ γ, lam' * S' γ = incrementAction q n γ + c := by
  rw [action_unique_up_to_gauge pi (gibbsKernel q) n hpi
    (gibbsKernel_pos q) hpi1 (gibbsKernel_row_sum q) S' lam']
  constructor
  · rintro ⟨c, hc⟩
    refine ⟨↑n * Real.log (gibbsNorm q) + c, fun γ => ?_⟩
    rw [hc γ, logTransition_eq_increment_add_const q n γ]
    ring
  · rintro ⟨c, hc⟩
    refine ⟨c - ↑n * Real.log (gibbsNorm q), fun γ => ?_⟩
    rw [hc γ, logTransition_eq_increment_add_const q n γ]
    ring

/-- **Least-action readout.**  For two paths with a common start, the
Gibbs path law orders them oppositely to the increment action: the
most probable paths from a common start are the least-action paths of
the same `q` used to construct the kernel.  Equal initial weight and
equal normalizer cancel, and the exponential is monotone. -/
theorem same_start_most_probable_iff_least_increment {Ω : Type*}
    [AddCommGroup Ω] [Fintype Ω] [Nonempty Ω] (pi q : Ω → ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (γ γ' : PathSpace Ω n) (h0 : γ 0 = γ' 0) :
    markovPathLaw pi (gibbsKernel q) n γ'
        ≤ markovPathLaw pi (gibbsKernel q) n γ
      ↔ incrementAction q n γ ≤ incrementAction q n γ' := by
  rw [markovPathLaw_gibbs_closed pi q n γ,
    markovPathLaw_gibbs_closed pi q n γ', ← h0,
    div_le_div_iff_of_pos_right (pow_pos (gibbsNorm_pos q) n),
    mul_le_mul_iff_of_pos_left (hpi (γ 0)), Real.exp_le_exp,
    neg_le_neg_iff]

/-- **One common cost scale is the multiplier gauge.**  The tilt of the
declared reference by the increment action of `q` at multiplier `g`
equals the Markov path law of the Gibbs kernel with rescaled cost
`g * q`.  This single scalar absorbs a common overall scale, or the
scale of an isolated factor.  It does not absorb independent relative
couplings between P/F/G or between several simple factors. -/
theorem coupling_scale_is_multiplier_gauge {Ω : Type*}
    [AddCommGroup Ω] [Fintype Ω] [Nonempty Ω] (pi q : Ω → ℝ) (n : ℕ)
    (g : ℝ) (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction q n) g
      = markovPathLaw pi (gibbsKernel fun v => g * q v) n := by
  rw [markov_path_law_eq_gibbs pi (gibbsKernel fun v => g * q v) n
      (gibbsKernel_pos _) hpi1 (gibbsKernel_row_sum _),
    show logTransitionAction (gibbsKernel fun v => g * q v) n
        = fun γ => g * incrementAction q n γ
          + ↑n * Real.log (gibbsNorm fun v => g * q v) from
      funext fun γ => by
        rw [logTransition_eq_increment_add_const (fun v => g * q v) n γ,
          increment_action_mul g q n γ],
    tilt_action_add_const (stepUniformRef pi n)
      (fun γ => g * incrementAction q n γ) 1
      (↑n * Real.log (gibbsNorm fun v => g * q v))
      (stepUniformRef_pos pi n hpi),
    tilt_action_multiplier_rescale (stepUniformRef pi n)
      (incrementAction q n) g (stepUniformRef_pos pi n hpi)]

/-! ## Layer 2: one certified three-dimensional P-factor increment -/

/-- The centered lift of the step group `ZMod 3` into the integer
increment lattice: `0 ↦ 0`, `1 ↦ 1`, `2 ↦ -1`. -/
def centeredLift : ZMod 3 → ℤ :=
  fun x => if x = 0 then 0 else if x = 1 then 1 else -1

/-- The sector step group: one `ZMod 3` factor per gauge direction of
the single three-dimensional P block `[0, 3)`, with the basis increment
lattice truncated to steps `-1, 0, +1` under the centered lift.  This
is not the whole six-dimensional P derived algebra. -/
abbrev SectorStep : Type := ZMod 3 × ZMod 3 × ZMod 3

/-- The lifted coordinate vector of a sector step. -/
def sectorLift (v : SectorStep) : Fin 3 → ℤ :=
  ![centeredLift v.1, centeredLift v.2.1, centeredLift v.2.2]

/-- Rational parts of the certified kinetic matrix on gauge block
`[0, 3)` of the P family: entry `[i, j, n1, d1, n2, d2]` of
`families.P.derived_algebra.kinetic_form_matrix_upper_entries` in
`code/e9_kinetic/kinetic_form_v1.json` contributes `n1/d1` here; all
recorded denominators are `1`. -/
def sectorKineticRat : Fin 3 → Fin 3 → ℤ :=
  ![![1200, 800, -800], ![800, 1200, -800], ![-800, -800, 1200]]

/-- `sqrt 5` coefficients of the certified kinetic matrix on gauge
block `[0, 3)` of the P family: entry `[i, j, n1, d1, n2, d2]`
contributes `n2/d2` here; all recorded denominators are `1`. -/
def sectorKineticIrr : Fin 3 → Fin 3 → ℤ :=
  ![![-400, -400, 400], ![-400, -400, 400], ![400, 400, -400]]

/-- The certified kinetic matrix over `ℝ`:
`K i j = sectorKineticRat i j + sectorKineticIrr i j * sqrt 5`. -/
noncomputable def sectorKinetic (i j : Fin 3) : ℝ :=
  (sectorKineticRat i j : ℝ) + (sectorKineticIrr i j : ℝ) * Real.sqrt 5

/-- **The one-factor P kinetic increment.**  The quadratic form
`(lift v)ᵀ K (lift v)` of the certified kinetic matrix on the centered
lift of a sector step. -/
noncomputable def sectorIncrement (v : SectorStep) : ℝ :=
  ∑ i, ∑ j, (sectorLift v i : ℝ) * sectorKinetic i j * (sectorLift v j : ℝ)

/-- Rational part of the kinetic increment, over `ℤ`. -/
def sectorIncrementRat (v : SectorStep) : ℤ :=
  ∑ i, ∑ j, sectorLift v i * sectorKineticRat i j * sectorLift v j

/-- `sqrt 5` coefficient of the kinetic increment, over `ℤ`. -/
def sectorIncrementIrr (v : SectorStep) : ℤ :=
  ∑ i, ∑ j, sectorLift v i * sectorKineticIrr i j * sectorLift v j

/-- The kinetic increment splits into its exact `ℤ[√5]` components. -/
theorem sectorIncrement_eq_pair (v : SectorStep) :
    sectorIncrement v
      = (sectorIncrementRat v : ℝ)
        + (sectorIncrementIrr v : ℝ) * Real.sqrt 5 := by
  unfold sectorIncrement sectorIncrementRat sectorIncrementIrr
    sectorKinetic
  push_cast
  rw [Finset.sum_mul, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_mul, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  ring

/-- Positivity transfer for exact `ℤ[√5]` values: `p + q √5` is
positive once `p` is positive and either `q` is nonnegative or
`5 q² < p²`. -/
theorem pair_pos_of_criterion (p q : ℤ) (hp : 0 < p)
    (h : 0 ≤ q ∨ 5 * q ^ 2 < p ^ 2) :
    0 < (p : ℝ) + (q : ℝ) * Real.sqrt 5 := by
  have hpr : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  rcases h with hq | hq
  · have hqr : (0 : ℝ) ≤ (q : ℝ) := by exact_mod_cast hq
    have hprod := mul_nonneg hqr (Real.sqrt_nonneg 5)
    linarith
  · have hltr : 5 * (q : ℝ) ^ 2 < (p : ℝ) ^ 2 := by exact_mod_cast hq
    by_contra hcon
    have hle : (p : ℝ) + (q : ℝ) * Real.sqrt 5 ≤ 0 := not_lt.mp hcon
    have h1 : (p : ℝ) ≤ -(q : ℝ) * Real.sqrt 5 := by linarith
    have h2 : (p : ℝ) * (p : ℝ)
        ≤ (-(q : ℝ) * Real.sqrt 5) * (-(q : ℝ) * Real.sqrt 5) :=
      mul_self_le_mul_self (le_of_lt hpr) h1
    have h3 : (-(q : ℝ) * Real.sqrt 5) * (-(q : ℝ) * Real.sqrt 5)
        = (q : ℝ) ^ 2 * (Real.sqrt 5 * Real.sqrt 5) := by ring
    rw [h3, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 5)] at h2
    nlinarith

set_option maxRecDepth 8000 in
/-- Kernel-decided sign data of the 26 nonzero sector steps: the
rational part is positive and dominates the `sqrt 5` part in the exact
integer sense of `pair_pos_of_criterion`. -/
theorem sector_pairs_criterion :
    ∀ v : SectorStep, v ≠ 0 →
      0 < sectorIncrementRat v
        ∧ (0 ≤ sectorIncrementIrr v
            ∨ 5 * sectorIncrementIrr v ^ 2 < sectorIncrementRat v ^ 2) := by
  decide

/-- **Positive definiteness receipt for the selected block.**  The
kinetic increment is strictly positive on every nonzero step of this
27-element, three-dimensional step group.  This exhaustive receipt is
for block `[0, 3)` only; it is not a six-dimensional P-family or
multi-family inertia theorem. -/
theorem sectorIncrement_pos (v : SectorStep) (hv : v ≠ 0) :
    0 < sectorIncrement v := by
  rw [sectorIncrement_eq_pair v]
  exact pair_pos_of_criterion _ _ (sector_pairs_criterion v hv).1
    (sector_pairs_criterion v hv).2

/-- The kinetic increment vanishes at the zero step. -/
theorem sectorIncrement_zero : sectorIncrement 0 = 0 := by
  rw [sectorIncrement_eq_pair,
    show sectorIncrementRat 0 = 0 from by decide,
    show sectorIncrementIrr 0 = 0 from by decide]
  norm_num

/-- Exact receipt on the first unit step: the diagonal entry
`1200 - 400 √5`. -/
theorem sectorIncrement_unit_step :
    sectorIncrement (1, 0, 0) = 1200 - 400 * Real.sqrt 5 := by
  rw [sectorIncrement_eq_pair,
    show sectorIncrementRat (1, 0, 0) = 1200 from by decide,
    show sectorIncrementIrr (1, 0, 0) = -400 from by decide]
  push_cast
  ring

/-- Exact receipt on the double step `(1, 1, 0)`:
`2 * (1200 - 400 √5) + 2 * (800 - 400 √5) = 4000 - 1600 √5`. -/
theorem sectorIncrement_double_step :
    sectorIncrement (1, 1, 0) = 4000 - 1600 * Real.sqrt 5 := by
  rw [sectorIncrement_eq_pair,
    show sectorIncrementRat (1, 1, 0) = 4000 from by decide,
    show sectorIncrementIrr (1, 1, 0) = -1600 from by decide]
  push_cast
  ring

/-- The unit-step value is positive as a displayed literal, through
the same exact sign analysis. -/
theorem sector_unit_step_value_pos :
    (0 : ℝ) < 1200 - 400 * Real.sqrt 5 := by
  rw [← sectorIncrement_unit_step]
  exact sectorIncrement_pos _ (by decide)

/-! ### The earned binding instance at the certified kinetic form -/

/-- **One-factor P keystone.**  The log-transition action of the Gibbs
kernel constructed from this three-dimensional P increment is the
kinetic action plus the per-step gauge constant. -/
theorem sector_logTransition_eq_kinetic_action (n : ℕ)
    (γ : PathSpace SectorStep n) :
    logTransitionAction (gibbsKernel sectorIncrement) n γ
      = incrementAction sectorIncrement n γ
        + n * Real.log (gibbsNorm sectorIncrement) :=
  logTransition_eq_increment_add_const sectorIncrement n γ

/-- The certified one-factor kinetic action reproduces the path law of
the Gibbs kernel constructed from that same increment, over the
declared reference at multiplier one. -/
theorem sector_kinetic_action_reproduces_law (pi : SectorStep → ℝ)
    (n : ℕ) (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction sectorIncrement n) 1
      = markovPathLaw pi (gibbsKernel sectorIncrement) n :=
  increment_action_reproduces_law pi sectorIncrement n hpi hpi1

/-- The sector binding is unique up to the additive-constant gauge:
any pair reproducing the sector path law carries the kinetic action. -/
theorem sector_binding_unique_up_to_gauge (pi : SectorStep → ℝ)
    (n : ℕ) (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1)
    (S' : PathSpace SectorStep n → ℝ) (lam' : ℝ) :
    tilt (stepUniformRef pi n) S' lam'
        = markovPathLaw pi (gibbsKernel sectorIncrement) n
      ↔ ∃ c : ℝ, ∀ γ,
          lam' * S' γ = incrementAction sectorIncrement n γ + c :=
  binding_unique_up_to_gauge pi sectorIncrement n hpi hpi1 S' lam'

/-- Least-action readout at the certified form: from a common start,
more probable sector histories are histories of smaller kinetic
action. -/
theorem sector_most_probable_iff_least_kinetic (pi : SectorStep → ℝ)
    (n : ℕ) (hpi : ∀ x, 0 < pi x) (γ γ' : PathSpace SectorStep n)
    (h0 : γ 0 = γ' 0) :
    markovPathLaw pi (gibbsKernel sectorIncrement) n γ'
        ≤ markovPathLaw pi (gibbsKernel sectorIncrement) n γ
      ↔ incrementAction sectorIncrement n γ
        ≤ incrementAction sectorIncrement n γ' :=
  same_start_most_probable_iff_least_increment pi sectorIncrement n hpi
    γ γ' h0

/-- The single factor's common scale is the multiplier gauge:
rescaling this one certified three-dimensional P form by `g` produces
the constructed path law read off by the unrescaled action at
multiplier `g`. This one-factor theorem alone does not fix relative couplings.
The companion `TwoFactorHistoryBinding` module binds full P and color-F inside
the constructed Gibbs family and proves relative-coefficient identifiability
without source-selecting its value. -/
theorem sector_coupling_is_multiplier (pi : SectorStep → ℝ) (n : ℕ)
    (g : ℝ) (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction sectorIncrement n) g
      = markovPathLaw pi (gibbsKernel fun v => g * sectorIncrement v) n :=
  coupling_scale_is_multiplier_gauge pi sectorIncrement n g hpi hpi1

/-! ## Layer 3: negative control -/

/-- A strictly positive row-stochastic kernel on `ZMod 2` with
state-dependent rows `(1/3, 2/3)` and `(1/4, 3/4)`. -/
noncomputable def dependentKernel : ZMod 2 → ZMod 2 → ℝ :=
  fun x y =>
    if x = 0 then (if y = 0 then 1 / 3 else 2 / 3)
    else if y = 0 then 1 / 4 else 3 / 4

theorem dependentKernel_apply_00 : dependentKernel 0 0 = 1 / 3 := by
  unfold dependentKernel
  rw [if_pos rfl, if_pos rfl]

theorem dependentKernel_apply_01 : dependentKernel 0 1 = 2 / 3 := by
  unfold dependentKernel
  rw [if_pos rfl, if_neg (by decide : ¬(1 : ZMod 2) = 0)]

theorem dependentKernel_apply_10 : dependentKernel 1 0 = 1 / 4 := by
  unfold dependentKernel
  rw [if_neg (by decide : ¬(1 : ZMod 2) = 0), if_pos rfl]

theorem dependentKernel_apply_11 : dependentKernel 1 1 = 3 / 4 := by
  unfold dependentKernel
  rw [if_neg (by decide : ¬(1 : ZMod 2) = 0),
    if_neg (by decide : ¬(1 : ZMod 2) = 0)]

theorem zmod2_cases (x : ZMod 2) : x = 0 ∨ x = 1 := by
  revert x
  decide

theorem zmod2_univ : (Finset.univ : Finset (ZMod 2)) = {0, 1} := by
  decide

theorem zmod2_sum (f : ZMod 2 → ℝ) : ∑ y, f y = f 0 + f 1 := by
  rw [zmod2_univ, Finset.sum_insert (by decide), Finset.sum_singleton]

/-- The control kernel is strictly positive. -/
theorem dependentKernel_pos (x y : ZMod 2) : 0 < dependentKernel x y := by
  rcases zmod2_cases x with hx | hx <;>
    rcases zmod2_cases y with hy | hy <;> subst hx <;> subst hy <;>
    first
      | (rw [dependentKernel_apply_00]; norm_num)
      | (rw [dependentKernel_apply_01]; norm_num)
      | (rw [dependentKernel_apply_10]; norm_num)
      | (rw [dependentKernel_apply_11]; norm_num)

/-- The control kernel is row-stochastic. -/
theorem dependentKernel_row_sum (x : ZMod 2) :
    ∑ y, dependentKernel x y = 1 := by
  rw [zmod2_sum (dependentKernel x)]
  rcases zmod2_cases x with hx | hx <;> subst hx
  · rw [dependentKernel_apply_00, dependentKernel_apply_01]
    norm_num
  · rw [dependentKernel_apply_10, dependentKernel_apply_11]
    norm_num

/-- The rows of the control kernel differ: the kernel is a function of
the pair, and of no increment. -/
theorem dependentKernel_state_dependent :
    dependentKernel 0 0 ≠ dependentKernel 1 0 := by
  rw [dependentKernel_apply_00, dependentKernel_apply_10]
  norm_num

/-- One-transition evaluation of the log-transition action. -/
theorem logTransitionAction_one {Ω : Type*} [Fintype Ω]
    (P : Ω → Ω → ℝ) (γ : PathSpace Ω 1) :
    logTransitionAction P 1 γ = -Real.log (P (γ 0) (γ 1)) := by
  unfold logTransitionAction
  rw [Fin.sum_univ_one,
    show ((0 : Fin 1).castSucc) = (0 : Fin 2) from by decide,
    show ((0 : Fin 1).succ) = (1 : Fin 2) from by decide]

/-- One-transition evaluation of the increment action. -/
theorem incrementAction_one {Ω : Type*} [AddCommGroup Ω] (q : Ω → ℝ)
    (γ : PathSpace Ω 1) :
    incrementAction q 1 γ = q (γ 1 - γ 0) := by
  unfold incrementAction
  rw [Fin.sum_univ_one,
    show ((0 : Fin 1).castSucc) = (0 : Fin 2) from by decide,
    show ((0 : Fin 1).succ) = (1 : Fin 2) from by decide]

/-- **Keystone negative control: the binding requires translation
invariance.**  The log-transition action of the state-dependent
control kernel lies outside the additive-constant gauge orbit of every
increment action: the constant paths at `0` and at `1` carry the same
increment `0` and the distinct transition weights `1/3` and `3/4`. -/
theorem binding_requires_invariance :
    ¬ ∃ (q : ZMod 2 → ℝ) (c : ℝ), ∀ γ : PathSpace (ZMod 2) 1,
        logTransitionAction dependentKernel 1 γ
          = incrementAction q 1 γ + c := by
  rintro ⟨q, c, hqc⟩
  have h00 := hqc fun _ => 0
  have h11 := hqc fun _ => 1
  rw [logTransitionAction_one, incrementAction_one] at h00 h11
  rw [dependentKernel_apply_00, sub_self] at h00
  rw [dependentKernel_apply_11, sub_self] at h11
  have hlog : Real.log ((1 : ℝ) / 3) = Real.log ((3 : ℝ) / 4) := by
    linarith
  have h13 : ((1 : ℝ) / 3) = 3 / 4 := by
    rw [← Real.exp_log (by norm_num : (0 : ℝ) < 1 / 3), hlog,
      Real.exp_log (by norm_num : (0 : ℝ) < 3 / 4)]
  norm_num at h13

/-- The control escapes at every multiplier as well: no rescaled
increment action plus constant reproduces the log-transition action of
the state-dependent kernel. -/
theorem binding_requires_invariance_multiplier :
    ¬ ∃ (q : ZMod 2 → ℝ) (lam c : ℝ), ∀ γ : PathSpace (ZMod 2) 1,
        logTransitionAction dependentKernel 1 γ
          = lam * incrementAction q 1 γ + c := by
  rintro ⟨q, lam, c, h⟩
  exact binding_requires_invariance
    ⟨fun v => lam * q v, c, fun γ => by
      rw [h γ, increment_action_mul lam q 1 γ]⟩

/-! ## Axiom audit -/

#print axioms gibbsKernel_pos
#print axioms gibbsKernel_row_sum
#print axioms logTransition_eq_increment_add_const
#print axioms increment_action_reproduces_law
#print axioms binding_unique_up_to_gauge
#print axioms same_start_most_probable_iff_least_increment
#print axioms coupling_scale_is_multiplier_gauge
#print axioms sectorIncrement_eq_pair
#print axioms sector_pairs_criterion
#print axioms sectorIncrement_pos
#print axioms sectorIncrement_zero
#print axioms sectorIncrement_unit_step
#print axioms sectorIncrement_double_step
#print axioms sector_unit_step_value_pos
#print axioms sector_logTransition_eq_kinetic_action
#print axioms sector_kinetic_action_reproduces_law
#print axioms sector_binding_unique_up_to_gauge
#print axioms sector_most_probable_iff_least_kinetic
#print axioms sector_coupling_is_multiplier
#print axioms dependentKernel_pos
#print axioms dependentKernel_row_sum
#print axioms binding_requires_invariance
#print axioms binding_requires_invariance_multiplier

end OPH.QFT
