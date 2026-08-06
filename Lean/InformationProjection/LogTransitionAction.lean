import InformationProjection.SourceHistoryPacket

namespace OPH.InformationProjection

open OPH.Thermodynamics

/-!
# The log-transition action: the source dynamics selects its action (B7, #683)

`SourceHistoryPacket` inhabits the `HistoryLaw` interface at literals from
one preregistered run and leaves the action slot and the multiplier slot as
declared inputs.  This module closes the which-action half at the
representation level: for Markov source dynamics the action is derived from
the dynamics, with no free slot.

For a Markov chain with strictly positive transition matrix `P` and initial
law `pi` on a finite state space, the path law
`pi(g0) * P(g0,g1) * ... * P(g_{n-1},g_n)` is itself a Gibbs law: it equals
the exponential tilt of the declared reference (the initial law times the
uniform-step counting weight `card^(-n)`) by the log-transition action
`S(γ) = -∑ i, log P(γ_i, γ_{i+1})` at multiplier one, with the partition
constant computed exactly as `card^(-n)` (`markov_path_law_eq_gibbs`,
`markov_tiltZ_eq`).  Strict positivity of `P` matches the committed chain,
whose integer transition counts are all positive
(`mixingChainCounts_pos`); the zero-transition case is excluded by that
hypothesis rather than handled by support restriction.

The gauge freedom is characterized exactly: two action-multiplier pairs
produce the same tilt at a fixed strictly positive reference precisely when
their multiplier-weighted actions differ by one additive constant
(`tilt_eq_tilt_iff_gauge`).  Additive shifts of the action and the joint
rescaling of the multiplier against the action scale are the two closure
instances (`tilt_action_add_const`, `tilt_action_multiplier_rescale`), and
per-step additive constants aggregate to a path constant, so they fall under
the same characterization.  Combining with the derived Gibbs form: an
action-multiplier pair reproduces the Markov path law over the declared
reference exactly when its multiplier-weighted action equals the
log-transition action plus a constant (`action_unique_up_to_gauge`).  The
source dynamics selects its action up to this gauge, and the normalization
convention that the action is the bare log-transition sum pins the
multiplier to one; the residual freedom is the exact trade of multiplier
against action scale (`markov_multiplier_action_tradeoff`), so the B7
multiplier-selection gap reduces to the same gauge fixing.

The earned instantiation reads the committed two-state mixing chain through
this derivation.  The transition products of the eight length-three paths
are exact rational literals tied to the committed integer counts by
kernel-decided identities, the log-transition action takes the value
`-log(transition product)` on each path, and the tilt of the declared
reference by that action at multiplier one equals the committed
`sourceTauChainR` path law exactly, with partition constant `1/4`
(`sourceLogTilt_eq_chain`, `sourceLogTiltZ_eq`).  A negative control shows
the committed repair-count action reproduces the chain path law at no
multiplier: it takes equal values on two paths that the chain path law
separates (`sourceRepairAction_no_multiplier`).

Scope boundary: representation-level derivation over named finite data and,
in the instantiation, over the committed literals of one preregistered run.
The reference stays a declared object (initial law times uniform-step
counting weight); the theorems state that once the reference is declared,
the source dynamics forces the action up to gauge and the convention pins
the multiplier.
-/

/-! ## Gauge freedom of the exponential tilt -/

section Gauge

variable {Γ : Type*} [Fintype Γ]

/-- Gauge construction half: if the multiplier-weighted actions differ by a
constant, the tilts agree.  The constant absorbs into the normalizer. -/
theorem tilt_gauge_shift [Nonempty Γ] (τ S S' : Γ → ℝ) (lam lam' : ℝ)
    (hτ : ∀ g, 0 < τ g) (c : ℝ)
    (hc : ∀ g, lam' * S' g = lam * S g + c) :
    tilt τ S' lam' = tilt τ S lam := by
  have hZ : tiltZ τ S' lam' = Real.exp (-c) * tiltZ τ S lam := by
    unfold tiltZ
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun g _ => ?_
    rw [hc g, show -(lam * S g + c) = -c + -(lam * S g) from by ring,
      Real.exp_add]
    ring
  funext g
  unfold tilt
  rw [hZ, hc g, show -(lam * S g + c) = -c + -(lam * S g) from by ring,
    Real.exp_add]
  have hZpos := tiltZ_pos τ S lam hτ
  have hepos := Real.exp_pos (-c)
  field_simp

/-- **Exact gauge characterization.**  At a fixed strictly positive
reference, two action-multiplier pairs produce the same exponential tilt
precisely when their multiplier-weighted actions differ by one additive
constant.  These affine reparameterizations are the entire invariance group
of the Gibbs representation at fixed reference. -/
theorem tilt_eq_tilt_iff_gauge [Nonempty Γ] (τ S S' : Γ → ℝ)
    (lam lam' : ℝ) (hτ : ∀ g, 0 < τ g) :
    tilt τ S' lam' = tilt τ S lam
      ↔ ∃ c : ℝ, ∀ g, lam' * S' g = lam * S g + c := by
  constructor
  · intro h
    refine ⟨Real.log (tiltZ τ S lam) - Real.log (tiltZ τ S' lam'),
      fun g => ?_⟩
    have h1 := tilt_log τ S lam hτ g
    have h2 := tilt_log τ S' lam' hτ g
    rw [h] at h2
    linarith
  · rintro ⟨c, hc⟩
    exact tilt_gauge_shift τ S S' lam lam' hτ c hc

/-- Gauge instance one: adding a constant to the action leaves the tilt
invariant.  Per-step additive constants on a path action aggregate to a path
constant, so they fall under this instance. -/
theorem tilt_action_add_const [Nonempty Γ] (τ S : Γ → ℝ) (lam c : ℝ)
    (hτ : ∀ g, 0 < τ g) :
    tilt τ (fun g => S g + c) lam = tilt τ S lam :=
  tilt_gauge_shift τ S _ lam lam hτ (lam * c) fun g => by ring

/-- Gauge instance two: the multiplier trades exactly against the action
scale.  The tilt at multiplier `lam` equals the tilt of the rescaled action
`lam * S` at multiplier one, so multiplier selection and action-scale
normalization are one choice, made once. -/
theorem tilt_action_multiplier_rescale [Nonempty Γ] (τ S : Γ → ℝ)
    (lam : ℝ) (hτ : ∀ g, 0 < τ g) :
    tilt τ (fun g => lam * S g) 1 = tilt τ S lam :=
  tilt_gauge_shift τ S _ lam 1 hτ 0 fun g => by ring

/-- Two-point reading of the gauge: tilts agreeing at a positive reference
force equal multiplier-weighted action differences on every pair of
histories. -/
theorem tilt_gauge_two_point [Nonempty Γ] (τ S S' : Γ → ℝ)
    (lam lam' : ℝ) (hτ : ∀ g, 0 < τ g)
    (h : tilt τ S' lam' = tilt τ S lam) (g g' : Γ) :
    lam' * S' g - lam' * S' g' = lam * S g - lam * S g' := by
  obtain ⟨c, hc⟩ := (tilt_eq_tilt_iff_gauge τ S S' lam lam' hτ).mp h
  rw [hc g, hc g']
  ring

end Gauge

/-! ## The Markov path law is a Gibbs law of the log-transition action -/

section MarkovGibbs

variable {Ω : Type*} [Fintype Ω]

/-- The Markov path law on `n`-transition paths: initial mass at the first
state times the product of the traversed transition probabilities. -/
noncomputable def markovPathLaw (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (γ : PathSpace Ω n) : ℝ :=
  pi (γ 0) * ∏ i : Fin n, P (γ i.castSucc) (γ i.succ)

/-- **The log-transition action.**  Minus the sum of the log transition
probabilities along the path: the action the source dynamics itself
produces, with no declared slot. -/
noncomputable def logTransitionAction (P : Ω → Ω → ℝ) (n : ℕ)
    (γ : PathSpace Ω n) : ℝ :=
  -∑ i : Fin n, Real.log (P (γ i.castSucc) (γ i.succ))

/-- The declared reference: the initial law times the uniform-step counting
weight `card^(-n)`.  This is the path law of pure counting dynamics, with
every step uniform. -/
noncomputable def stepUniformRef (pi : Ω → ℝ) (n : ℕ)
    (γ : PathSpace Ω n) : ℝ :=
  pi (γ 0) * ((Fintype.card Ω : ℝ))⁻¹ ^ n

/-- The reference is itself a Markov path law: the one with the uniform
transition kernel. -/
theorem stepUniformRef_eq_uniform_markov (pi : Ω → ℝ) (n : ℕ) :
    stepUniformRef pi n
      = markovPathLaw pi (fun _ _ => ((Fintype.card Ω : ℝ))⁻¹) n := by
  funext γ
  unfold stepUniformRef markovPathLaw
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]

omit [Fintype Ω] in
theorem markovPathLaw_pos (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (hP : ∀ x y, 0 < P x y) (γ : PathSpace Ω n) :
    0 < markovPathLaw pi P n γ :=
  mul_pos (hpi _) (Finset.prod_pos fun _ _ => hP _ _)

theorem stepUniformRef_pos [Nonempty Ω] (pi : Ω → ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (γ : PathSpace Ω n) :
    0 < stepUniformRef pi n γ := by
  have hcard : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
    exact_mod_cast Fintype.card_pos
  exact mul_pos (hpi _) (pow_pos (inv_pos.mpr hcard) n)

omit [Fintype Ω] in
/-- One-step extension of the Markov path law: appending a state multiplies
the weight by the corresponding transition probability. -/
theorem markovPathLaw_snoc (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (γ : PathSpace Ω n) (x : Ω) :
    markovPathLaw pi P (n + 1) (Fin.snoc γ x)
      = markovPathLaw pi P n γ * P (γ (Fin.last n)) x := by
  unfold markovPathLaw
  have h0 : (Fin.snoc γ x : PathSpace Ω (n + 1)) 0 = γ 0 := by
    rw [show (0 : Fin (n + 2)) = Fin.castSucc 0 from (Fin.castSucc_zero).symm,
      Fin.snoc_castSucc]
  have hstep : ∀ i : Fin n,
      P ((Fin.snoc γ x : PathSpace Ω (n + 1)) (Fin.castSucc i).castSucc)
          ((Fin.snoc γ x : PathSpace Ω (n + 1)) (Fin.castSucc i).succ)
        = P (γ i.castSucc) (γ i.succ) := by
    intro i
    rw [Fin.snoc_castSucc, Fin.succ_castSucc, Fin.snoc_castSucc]
  have hlast :
      P ((Fin.snoc γ x : PathSpace Ω (n + 1)) (Fin.last n).castSucc)
          ((Fin.snoc γ x : PathSpace Ω (n + 1)) (Fin.last n).succ)
        = P (γ (Fin.last n)) x := by
    rw [Fin.snoc_castSucc, Fin.succ_last, Fin.snoc_last]
  rw [Fin.prod_univ_castSucc, h0, hlast,
    Finset.prod_congr rfl fun i _ => hstep i]
  ring

/-- The Markov path law is normalized: row stochasticity of the kernel and
normalization of the initial law give unit total mass on every path
length. -/
theorem markovPathLaw_sum (pi : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hpi1 : ∑ x, pi x = 1) (hP1 : ∀ x, ∑ y, P x y = 1) (n : ℕ) :
    ∑ γ : PathSpace Ω n, markovPathLaw pi P n γ = 1 := by
  induction n with
  | zero =>
    rw [Fintype.sum_equiv (Equiv.funUnique (Fin 1) Ω)
      (fun γ => markovPathLaw pi P 0 γ) pi (fun γ => by
        simp [markovPathLaw])]
    exact hpi1
  | succ n ih =>
    have hbij : Function.Bijective
        (fun p : PathSpace Ω n × Ω =>
          (Fin.snoc p.1 p.2 : PathSpace Ω (n + 1))) :=
      Function.bijective_iff_has_inverse.mpr
        ⟨fun γ => (Fin.init γ, γ (Fin.last (n + 1))),
          fun p => by simp,
          fun γ => by simp [Fin.snoc_init_self]⟩
    rw [← Fintype.sum_bijective _ hbij
      (fun p : PathSpace Ω n × Ω =>
        markovPathLaw pi P (n + 1) (Fin.snoc p.1 p.2))
      (fun γ => markovPathLaw pi P (n + 1) γ) (fun p => rfl),
      Fintype.sum_prod_type]
    have hinner : ∀ γ : PathSpace Ω n,
        ∑ x, markovPathLaw pi P (n + 1) (Fin.snoc γ x)
          = markovPathLaw pi P n γ := by
      intro γ
      rw [Finset.sum_congr rfl fun x _ => markovPathLaw_snoc pi P n γ x,
        ← Finset.mul_sum, hP1 (γ (Fin.last n)), mul_one]
    rw [Finset.sum_congr rfl fun γ _ => hinner γ]
    exact ih

omit [Fintype Ω] in
/-- The Gibbs factor of the log-transition action at multiplier one is the
transition product exactly. -/
theorem exp_neg_logTransitionAction (P : Ω → Ω → ℝ) (n : ℕ)
    (hP : ∀ x y, 0 < P x y) (γ : PathSpace Ω n) :
    Real.exp (-(1 * logTransitionAction P n γ))
      = ∏ i : Fin n, P (γ i.castSucc) (γ i.succ) := by
  rw [one_mul]
  unfold logTransitionAction
  rw [neg_neg, Real.exp_sum]
  exact Finset.prod_congr rfl fun i _ => Real.exp_log (hP _ _)

/-- Pointwise Gibbs identity: reference times Gibbs factor equals the
uniform-step weight times the Markov path law. -/
theorem stepUniformRef_mul_exp (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (hP : ∀ x y, 0 < P x y) (γ : PathSpace Ω n) :
    stepUniformRef pi n γ * Real.exp (-(1 * logTransitionAction P n γ))
      = ((Fintype.card Ω : ℝ))⁻¹ ^ n * markovPathLaw pi P n γ := by
  rw [exp_neg_logTransitionAction P n hP γ]
  unfold stepUniformRef markovPathLaw
  ring

/-- **The partition constant is computed exactly.**  The tilt normalizer of
the declared reference at the log-transition action and multiplier one is
the uniform-step weight `card^(-n)`. -/
theorem markov_tiltZ_eq (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (hP : ∀ x y, 0 < P x y)
    (hpi1 : ∑ x, pi x = 1) (hP1 : ∀ x, ∑ y, P x y = 1) :
    tiltZ (stepUniformRef pi n) (logTransitionAction P n) 1
      = ((Fintype.card Ω : ℝ))⁻¹ ^ n := by
  unfold tiltZ
  rw [Finset.sum_congr rfl fun γ _ => stepUniformRef_mul_exp pi P n hP γ,
    ← Finset.mul_sum, markovPathLaw_sum pi P hpi1 hP1 n, mul_one]

/-- **The Markov path law equals the Gibbs law of the log-transition
action.**  Over any finite state space, for a strictly positive row
stochastic kernel and a normalized initial law, the Markov path law on
`n`-transition paths is the exponential tilt of the declared reference
(initial law times uniform-step counting weight) by the log-transition
action at multiplier one.  The action is derived from the dynamics; the
exponential form carries no declared slot beyond the reference. -/
theorem markov_path_law_eq_gibbs [Nonempty Ω] (pi : Ω → ℝ)
    (P : Ω → Ω → ℝ) (n : ℕ) (hP : ∀ x y, 0 < P x y)
    (hpi1 : ∑ x, pi x = 1) (hP1 : ∀ x, ∑ y, P x y = 1) :
    markovPathLaw pi P n
      = tilt (stepUniformRef pi n) (logTransitionAction P n) 1 := by
  have hcard : ((Fintype.card Ω : ℝ))⁻¹ ^ n ≠ 0 := by
    have h : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
      exact_mod_cast Fintype.card_pos
    positivity
  funext γ
  unfold tilt
  rw [markov_tiltZ_eq pi P n hP hpi1 hP1,
    stepUniformRef_mul_exp pi P n hP γ, mul_comm, mul_div_assoc,
    div_self hcard, mul_one]

/-- **The source dynamics selects the action (derivation corollary).**  An
action-multiplier pair reproduces the Markov path law as its Gibbs law over
the declared reference exactly when its multiplier-weighted action is the
log-transition action plus a constant.  Together with
`tilt_eq_tilt_iff_gauge` this closes the which-action question at the
representation level: the log-transition action is the unique action with
this property up to the additive-constant and multiplier-rescaling gauge,
and both gauge directions are realized (`tilt_action_add_const`,
`tilt_action_multiplier_rescale`). -/
theorem action_unique_up_to_gauge [Nonempty Ω] (pi : Ω → ℝ)
    (P : Ω → Ω → ℝ) (n : ℕ) (hpi : ∀ x, 0 < pi x)
    (hP : ∀ x y, 0 < P x y)
    (hpi1 : ∑ x, pi x = 1) (hP1 : ∀ x, ∑ y, P x y = 1)
    (S' : PathSpace Ω n → ℝ) (lam' : ℝ) :
    tilt (stepUniformRef pi n) S' lam' = markovPathLaw pi P n
      ↔ ∃ c : ℝ, ∀ γ, lam' * S' γ = logTransitionAction P n γ + c := by
  rw [markov_path_law_eq_gibbs pi P n hP hpi1 hP1,
    tilt_eq_tilt_iff_gauge (stepUniformRef pi n)
      (logTransitionAction P n) S' 1 lam'
      (stepUniformRef_pos pi n hpi)]
  exact exists_congr fun c => forall_congr' fun γ => by rw [one_mul]

/-- **The multiplier is pinned by the action-scale convention.**  At every
nonzero multiplier the rescaled action `lam⁻¹ * S_log` reproduces the Markov
path law, and by `action_unique_up_to_gauge` nothing outside this orbit
does: the multiplier trades exactly against the action scale, so fixing the
action as the bare log-transition sum pins the multiplier to one.  The B7
multiplier-selection gap over Markov source dynamics is this gauge
fixing. -/
theorem markov_multiplier_action_tradeoff [Nonempty Ω] (pi : Ω → ℝ)
    (P : Ω → Ω → ℝ) (n : ℕ) (hpi : ∀ x, 0 < pi x)
    (hP : ∀ x y, 0 < P x y)
    (hpi1 : ∑ x, pi x = 1) (hP1 : ∀ x, ∑ y, P x y = 1)
    (lam : ℝ) (hlam : lam ≠ 0) :
    tilt (stepUniformRef pi n)
        (fun γ => lam⁻¹ * logTransitionAction P n γ) lam
      = markovPathLaw pi P n := by
  rw [markov_path_law_eq_gibbs pi P n hP hpi1 hP1]
  exact tilt_gauge_shift (stepUniformRef pi n) (logTransitionAction P n) _
    1 lam (stepUniformRef_pos pi n hpi) 0
    fun γ => by rw [← mul_assoc, mul_inv_cancel₀ hlam]; ring

end MarkovGibbs

/-! ## The earned instantiation: the committed two-state mixing chain

The eight length-three paths of `SourceHistoryPacket` over the recurrent
class of the committed chain, encoded `g = 4*s0 + 2*s1 + s2`.  The integer
layer (transition-product numerators and denominators against the committed
counts, and the cleared-denominator Gibbs identity with the committed path
law) is kernel-decided on the literals.  The rational layer is elaborated
from the same literals by `norm_num`.  The real layer casts the literals
and lands the derived-action identity: the tilt of the declared reference
by the log-transition action at multiplier one is the committed
`sourceTauChainR` path law, with partition constant `1/4`. -/

/-- Numerators of the per-path transition products over the committed
integer counts: `count(s0,s1) * count(s1,s2)`. -/
def sourceTransProdNum : Fin 8 → ℕ :=
  ![1752976, 141668, 535, 53821, 6620, 535, 2515, 253009]

/-- Denominators of the per-path transition products: products of the
committed row sums `rowSum(s0) * rowSum(s1)`. -/
def sourceTransProdDen : Fin 8 → ℕ :=
  ![2047761, 2047761, 726948, 726948, 726948, 726948, 258064, 258064]

/-- The per-path transition products as exact rationals:
`P(s0,s1) * P(s1,s2)` over the committed chain entries. -/
def sourceTransProdQ : Fin 8 → ℚ :=
  ![1752976/2047761, 141668/2047761, 535/726948, 53821/726948,
    6620/726948, 535/726948, 2515/258064, 253009/258064]

/-- The transition-product numerators are the products of the committed
integer transition counts, kernel-decided. -/
theorem sourceTransProdNum_eq_counts :
    ∀ g : Fin 8, sourceTransProdNum g
      = mixingChainCounts (sourceState0 g) (sourceState1 g)
        * mixingChainCounts (sourceState1 g) (sourceState2 g) := by decide

/-- The transition-product denominators are the products of the committed
row sums, kernel-decided. -/
theorem sourceTransProdDen_eq_rowSums :
    ∀ g : Fin 8, sourceTransProdDen g
      = mixingChainRowSums (sourceState0 g)
        * mixingChainRowSums (sourceState1 g) := by decide

/-- **Cleared-denominator Gibbs identity, kernel-decided.**  Over ℕ, every
path satisfies
`tauChainNum g * (61511 * transDen g) = piNum(s0) * transNum g * tauChainDen g`:
the identity `tau_chain(g) = pi(s0) * (transition product)(g)`, the Gibbs
form of the committed path law at the log-transition action, with every
denominator cleared against the committed literals. -/
theorem sourceTauChain_gibbs_integer_identity :
    ∀ g : Fin 8, sourceTauChainNum g * (61511 * sourceTransProdDen g)
      = mixingChainStationaryNum (sourceState0 g) * sourceTransProdNum g
        * sourceTauChainDen g := by decide

theorem sourceTransProdQ_apply_0 : sourceTransProdQ 0 = 1752976/2047761 := rfl
theorem sourceTransProdQ_apply_1 : sourceTransProdQ 1 = 141668/2047761 := rfl
theorem sourceTransProdQ_apply_2 : sourceTransProdQ 2 = 535/726948 := rfl
theorem sourceTransProdQ_apply_3 : sourceTransProdQ 3 = 53821/726948 := rfl
theorem sourceTransProdQ_apply_4 : sourceTransProdQ 4 = 6620/726948 := rfl
theorem sourceTransProdQ_apply_5 : sourceTransProdQ 5 = 535/726948 := rfl
theorem sourceTransProdQ_apply_6 : sourceTransProdQ 6 = 2515/258064 := rfl
theorem sourceTransProdQ_apply_7 : sourceTransProdQ 7 = 253009/258064 := rfl

/-- The rational transition products are the quotients of the kernel-layer
integer literals, path by path. -/
theorem sourceTransProdQ_eq_num_div_den :
    ∀ g, sourceTransProdQ g
      = (sourceTransProdNum g : ℚ) / (sourceTransProdDen g : ℚ) := by
  have h0 : sourceTransProdQ 0
      = (sourceTransProdNum 0 : ℚ) / (sourceTransProdDen 0 : ℚ) := by
    rw [sourceTransProdQ_apply_0, show sourceTransProdNum 0 = 1752976 from rfl,
      show sourceTransProdDen 0 = 2047761 from rfl]
    norm_num
  have h1 : sourceTransProdQ 1
      = (sourceTransProdNum 1 : ℚ) / (sourceTransProdDen 1 : ℚ) := by
    rw [sourceTransProdQ_apply_1, show sourceTransProdNum 1 = 141668 from rfl,
      show sourceTransProdDen 1 = 2047761 from rfl]
    norm_num
  have h2 : sourceTransProdQ 2
      = (sourceTransProdNum 2 : ℚ) / (sourceTransProdDen 2 : ℚ) := by
    rw [sourceTransProdQ_apply_2, show sourceTransProdNum 2 = 535 from rfl,
      show sourceTransProdDen 2 = 726948 from rfl]
    norm_num
  have h3 : sourceTransProdQ 3
      = (sourceTransProdNum 3 : ℚ) / (sourceTransProdDen 3 : ℚ) := by
    rw [sourceTransProdQ_apply_3, show sourceTransProdNum 3 = 53821 from rfl,
      show sourceTransProdDen 3 = 726948 from rfl]
    norm_num
  have h4 : sourceTransProdQ 4
      = (sourceTransProdNum 4 : ℚ) / (sourceTransProdDen 4 : ℚ) := by
    rw [sourceTransProdQ_apply_4, show sourceTransProdNum 4 = 6620 from rfl,
      show sourceTransProdDen 4 = 726948 from rfl]
    norm_num
  have h5 : sourceTransProdQ 5
      = (sourceTransProdNum 5 : ℚ) / (sourceTransProdDen 5 : ℚ) := by
    rw [sourceTransProdQ_apply_5, show sourceTransProdNum 5 = 535 from rfl,
      show sourceTransProdDen 5 = 726948 from rfl]
    norm_num
  have h6 : sourceTransProdQ 6
      = (sourceTransProdNum 6 : ℚ) / (sourceTransProdDen 6 : ℚ) := by
    rw [sourceTransProdQ_apply_6, show sourceTransProdNum 6 = 2515 from rfl,
      show sourceTransProdDen 6 = 258064 from rfl]
    norm_num
  have h7 : sourceTransProdQ 7
      = (sourceTransProdNum 7 : ℚ) / (sourceTransProdDen 7 : ℚ) := by
    rw [sourceTransProdQ_apply_7, show sourceTransProdNum 7 = 253009 from rfl,
      show sourceTransProdDen 7 = 258064 from rfl]
    norm_num
  intro g
  fin_cases g
  · exact h0
  · exact h1
  · exact h2
  · exact h3
  · exact h4
  · exact h5
  · exact h6
  · exact h7

/-- The rational transition products are the products of committed chain
entries, path by path. -/
theorem sourceTransProdQ_eq_chain_product :
    ∀ g, sourceTransProdQ g
      = mixingChain (sourceState0 g) (sourceState1 g)
        * mixingChain (sourceState1 g) (sourceState2 g) := by
  have h0 : sourceTransProdQ 0 = mixingChain 0 0 * mixingChain 0 0 := by
    rw [sourceTransProdQ_apply_0]
    norm_num [mixingChain]
  have h1 : sourceTransProdQ 1 = mixingChain 0 0 * mixingChain 0 1 := by
    rw [sourceTransProdQ_apply_1]
    norm_num [mixingChain]
  have h2 : sourceTransProdQ 2 = mixingChain 0 1 * mixingChain 1 0 := by
    rw [sourceTransProdQ_apply_2]
    norm_num [mixingChain]
  have h3 : sourceTransProdQ 3 = mixingChain 0 1 * mixingChain 1 1 := by
    rw [sourceTransProdQ_apply_3]
    norm_num [mixingChain]
  have h4 : sourceTransProdQ 4 = mixingChain 1 0 * mixingChain 0 0 := by
    rw [sourceTransProdQ_apply_4]
    norm_num [mixingChain]
  have h5 : sourceTransProdQ 5 = mixingChain 1 0 * mixingChain 0 1 := by
    rw [sourceTransProdQ_apply_5]
    norm_num [mixingChain]
  have h6 : sourceTransProdQ 6 = mixingChain 1 1 * mixingChain 1 0 := by
    rw [sourceTransProdQ_apply_6]
    norm_num [mixingChain]
  have h7 : sourceTransProdQ 7 = mixingChain 1 1 * mixingChain 1 1 := by
    rw [sourceTransProdQ_apply_7]
    norm_num [mixingChain]
  intro g
  fin_cases g
  · exact h0
  · exact h1
  · exact h2
  · exact h3
  · exact h4
  · exact h5
  · exact h6
  · exact h7

/-- The transition products are strictly positive: the committed chain has
all entries positive, matching the positivity hypothesis of the general
theorem. -/
theorem sourceTransProdQ_pos : ∀ g, 0 < sourceTransProdQ g := by
  intro g
  rw [sourceTransProdQ_eq_chain_product g]
  exact mul_pos (mixingChain_entries_pos _ _) (mixingChain_entries_pos _ _)

/-- The committed chain path law in Gibbs form over ℚ: stationary mass at
the initial state times the transition product. -/
theorem sourceTauChainQ_eq_pi_mul_transProd :
    ∀ g, sourceTauChainQ g
      = mixingChainStationary (sourceState0 g) * sourceTransProdQ g := by
  intro g
  rw [sourceTransProdQ_eq_chain_product g, ← mul_assoc]
  exact sourceTauChainQ_eq_markov_product g

/-! ### The real layer: reference, action, and the derived-action identity -/

/-- The declared reference over ℝ: stationary mass at the initial state
times the uniform-step counting weight `(1/2)^2 = 1/4` of the two-state,
two-transition path space. -/
noncomputable def sourceLogRefR : Fin 8 → ℝ :=
  fun g => ((mixingChainStationary (sourceState0 g) : ℚ) : ℝ) / 4

/-- The transition product over ℝ. -/
noncomputable def sourceTransProdR : Fin 8 → ℝ :=
  fun g => ((sourceTransProdQ g : ℚ) : ℝ)

/-- **The log-transition action of the committed chain.**  Minus the log of
the transition product of the path: the action the committed source
dynamics selects, evaluated on the eight length-three paths. -/
noncomputable def sourceLogActionR : Fin 8 → ℝ :=
  fun g => -Real.log (sourceTransProdR g)

theorem sourceTransProdR_pos : ∀ g, 0 < sourceTransProdR g := fun g => by
  simp only [sourceTransProdR]
  exact_mod_cast sourceTransProdQ_pos g

theorem sourceLogRefR_pos : ∀ g, 0 < sourceLogRefR g := by
  intro g
  simp only [sourceLogRefR]
  have h : (0 : ℝ) < ((mixingChainStationary (sourceState0 g) : ℚ) : ℝ) := by
    exact_mod_cast mixingChainStationary_pos (sourceState0 g)
  positivity

/-- The declared reference is normalized: each initial state heads four of
the eight paths, and the stationary law is normalized. -/
theorem sourceLogRefR_sum : ∑ g, sourceLogRefR g = 1 := by
  simp only [sourceLogRefR]
  rw [Fin.sum_univ_eight,
    show sourceState0 0 = 0 from rfl, show sourceState0 1 = 0 from rfl,
    show sourceState0 2 = 0 from rfl, show sourceState0 3 = 0 from rfl,
    show sourceState0 4 = 1 from rfl, show sourceState0 5 = 1 from rfl,
    show sourceState0 6 = 1 from rfl, show sourceState0 7 = 1 from rfl,
    show mixingChainStationary 0 = 7155/61511 from rfl,
    show mixingChainStationary 1 = 54356/61511 from rfl]
  push_cast
  norm_num

/-- Per-path value of the log-transition action: path 0, the double stay at
state 18, scores `-log(1752976/2047761)`. -/
theorem sourceLogActionR_apply_0 :
    sourceLogActionR 0 = -Real.log ((1752976 : ℝ) / 2047761) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_0]
  norm_num

/-- Path 1 (stay, then move 18 to 19) scores `-log(141668/2047761)`. -/
theorem sourceLogActionR_apply_1 :
    sourceLogActionR 1 = -Real.log ((141668 : ℝ) / 2047761) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_1]
  norm_num

/-- Path 2 (move, then move back) scores `-log(535/726948)`. -/
theorem sourceLogActionR_apply_2 :
    sourceLogActionR 2 = -Real.log ((535 : ℝ) / 726948) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_2]
  norm_num

/-- Path 3 (move 18 to 19, then stay) scores `-log(53821/726948)`. -/
theorem sourceLogActionR_apply_3 :
    sourceLogActionR 3 = -Real.log ((53821 : ℝ) / 726948) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_3]
  norm_num

/-- Path 4 (move 19 to 18, then stay) scores `-log(6620/726948)`. -/
theorem sourceLogActionR_apply_4 :
    sourceLogActionR 4 = -Real.log ((6620 : ℝ) / 726948) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_4]
  norm_num

/-- Path 5 (move, then move back) scores `-log(535/726948)`. -/
theorem sourceLogActionR_apply_5 :
    sourceLogActionR 5 = -Real.log ((535 : ℝ) / 726948) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_5]
  norm_num

/-- Path 6 (stay, then move 19 to 18) scores `-log(2515/258064)`. -/
theorem sourceLogActionR_apply_6 :
    sourceLogActionR 6 = -Real.log ((2515 : ℝ) / 258064) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_6]
  norm_num

/-- Path 7, the double stay at state 19, scores `-log(253009/258064)`. -/
theorem sourceLogActionR_apply_7 :
    sourceLogActionR 7 = -Real.log ((253009 : ℝ) / 258064) := by
  simp only [sourceLogActionR, sourceTransProdR, sourceTransProdQ_apply_7]
  norm_num

/-- The Gibbs factor of the derived action at multiplier one is the
transition product, path by path. -/
theorem sourceExp_neg_logAction (g : Fin 8) :
    Real.exp (-(1 * sourceLogActionR g)) = sourceTransProdR g := by
  rw [one_mul]
  simp only [sourceLogActionR]
  rw [neg_neg, Real.exp_log (sourceTransProdR_pos g)]

/-- Pointwise Gibbs identity at the packet: reference times Gibbs factor is
one quarter of the committed chain path law. -/
theorem sourceLogRef_mul_exp (g : Fin 8) :
    sourceLogRefR g * Real.exp (-(1 * sourceLogActionR g))
      = sourceTauChainR g / 4 := by
  rw [sourceExp_neg_logAction g]
  simp only [sourceLogRefR, sourceTransProdR, sourceTauChainR]
  rw [div_mul_eq_mul_div, ← Rat.cast_mul,
    ← sourceTauChainQ_eq_pi_mul_transProd g]

/-- **The partition constant of the packet is exactly `1/4`.**  The tilt
normalizer of the declared reference at the derived action and multiplier
one is the uniform-step weight `(1/2)^2` of the committed path space. -/
theorem sourceLogTiltZ_eq :
    tiltZ sourceLogRefR sourceLogActionR 1 = 1 / 4 := by
  unfold tiltZ
  rw [Finset.sum_congr rfl fun g _ => sourceLogRef_mul_exp g,
    ← Finset.sum_div, sourceTauChainR_sum]

/-- **The committed chain path law is the Gibbs law of its log-transition
action.**  The tilt of the declared reference by the derived action at
multiplier one equals the committed `sourceTauChainR` path law exactly:
the earned instantiation of `markov_path_law_eq_gibbs` on the literals of
the preregistered run. -/
theorem sourceLogTilt_eq_chain :
    tilt sourceLogRefR sourceLogActionR 1 = sourceTauChainR := by
  funext g
  unfold tilt
  rw [sourceLogTiltZ_eq, sourceLogRef_mul_exp g]
  ring

/-- **The which-action question closes at the packet.**  An
action-multiplier pair reproduces the committed chain path law as its Gibbs
law over the declared reference exactly when its multiplier-weighted action
is the log-transition action plus a constant.  The committed source
dynamics selects its action up to the additive-constant and
multiplier-rescaling gauge. -/
theorem sourceAction_unique_up_to_gauge (S' : Fin 8 → ℝ) (lam' : ℝ) :
    tilt sourceLogRefR S' lam' = sourceTauChainR
      ↔ ∃ c : ℝ, ∀ g, lam' * S' g = sourceLogActionR g + c := by
  rw [← sourceLogTilt_eq_chain,
    tilt_eq_tilt_iff_gauge sourceLogRefR sourceLogActionR S' 1 lam'
      sourceLogRefR_pos]
  exact exists_congr fun c => forall_congr' fun g => by rw [one_mul]

/-- The residual gauge at the packet: every nonzero multiplier reproduces
the committed chain path law with the correspondingly rescaled action, so
the multiplier slot of the packet is the action-scale convention and
nothing else. -/
theorem sourceMultiplier_action_tradeoff (lam : ℝ) (hlam : lam ≠ 0) :
    tilt sourceLogRefR (fun g => lam⁻¹ * sourceLogActionR g) lam
      = sourceTauChainR := by
  rw [← sourceLogTilt_eq_chain]
  exact tilt_gauge_shift sourceLogRefR sourceLogActionR _ 1 lam
    sourceLogRefR_pos 0 fun g => by rw [← mul_assoc, mul_inv_cancel₀ hlam]; ring

/-! ### Negative control: the repair-count action fails at every multiplier -/

/-- Paths 1 and 3 share their initial state, so the declared reference
agrees on them. -/
theorem sourceLogRefR_one_eq_three : sourceLogRefR 1 = sourceLogRefR 3 := rfl

/-- The committed chain path law separates paths 1 and 3. -/
theorem sourceTauChainR_one_ne_three :
    sourceTauChainR 1 ≠ sourceTauChainR 3 := by
  rw [sourceTauChainR_apply_1, sourceTauChainR_apply_3]
  norm_num

/-- **Negative control.**  The committed repair-count action of
`SourceHistoryPacket` reproduces the chain path law at no multiplier: it
takes equal values on paths 1 and 3 and the reference agrees there, so
every Gibbs law built from it assigns those two paths equal mass, while the
committed chain path law separates them.  The exact two-path discrepancy is
`708340/88022241` against `269105/31247588`. -/
theorem sourceRepairAction_no_multiplier (lam : ℝ) :
    tilt sourceLogRefR sourceActionR lam ≠ sourceTauChainR := by
  intro h
  have h1 := congrFun h 1
  have h3 := congrFun h 3
  have heq : tilt sourceLogRefR sourceActionR lam 1
      = tilt sourceLogRefR sourceActionR lam 3 := by
    unfold tilt
    rw [sourceLogRefR_one_eq_three, sourceActionR_apply_1,
      sourceActionR_apply_3]
  rw [h1, h3] at heq
  exact sourceTauChainR_one_ne_three heq

/-- The repair-count action lies outside the gauge orbit of the derived
action: no multiplier and additive constant reconcile the two. -/
theorem sourceRepairAction_not_gauge :
    ¬ ∃ lam c : ℝ, ∀ g, lam * sourceActionR g = sourceLogActionR g + c := by
  rintro ⟨lam, c, hc⟩
  exact sourceRepairAction_no_multiplier lam
    ((sourceAction_unique_up_to_gauge sourceActionR lam).mpr ⟨c, hc⟩)

/-! ### The committed chain under the general theorem -/

/-- The committed stationary law over ℝ. -/
noncomputable def sourcePiR : Fin 2 → ℝ :=
  fun i => ((mixingChainStationary i : ℚ) : ℝ)

/-- The committed chain over ℝ. -/
noncomputable def sourcePR : Fin 2 → Fin 2 → ℝ :=
  fun i j => ((mixingChain i j : ℚ) : ℝ)

theorem sourcePiR_sum : ∑ i, sourcePiR i = 1 := by
  simp only [sourcePiR]
  rw [← Rat.cast_sum, mixingChainStationary_sum, Rat.cast_one]

theorem sourcePR_pos : ∀ i j, 0 < sourcePR i j := fun i j => by
  simp only [sourcePR]
  exact_mod_cast mixingChain_entries_pos i j

theorem sourcePR_rows : ∀ i, ∑ j, sourcePR i j = 1 := fun i => by
  simp only [sourcePR]
  rw [← Rat.cast_sum, mixingChain_row_stochastic i, Rat.cast_one]

/-- **General-theorem instantiation at the committed chain.**  The
two-transition Markov path law of the cast chain is the Gibbs law of its
log-transition action at multiplier one over the declared reference. -/
theorem sourceMarkov_eq_gibbs :
    markovPathLaw sourcePiR sourcePR 2
      = tilt (stepUniformRef sourcePiR 2) (logTransitionAction sourcePR 2) 1 :=
  markov_path_law_eq_gibbs sourcePiR sourcePR 2 sourcePR_pos
    sourcePiR_sum sourcePR_rows

/-- Two-transition evaluation of the Markov path law. -/
theorem markovPathLaw_two {Ω : Type*} [Fintype Ω] (pi : Ω → ℝ)
    (P : Ω → Ω → ℝ) (γ : PathSpace Ω 2) :
    markovPathLaw pi P 2 γ = pi (γ 0) * (P (γ 0) (γ 1) * P (γ 1) (γ 2)) := by
  unfold markovPathLaw
  rw [Fin.prod_univ_two]
  rfl

/-- The committed eight-path literals are the general Markov path law read
through the committed encoding: `tau_chain(g)` is the two-transition path
law of the cast chain on the decoded path of `g`. -/
theorem sourceTauChainR_eq_markovPathLaw (g : Fin 8) :
    sourceTauChainR g
      = markovPathLaw sourcePiR sourcePR 2
          ![sourceState0 g, sourceState1 g, sourceState2 g] := by
  rw [markovPathLaw_two]
  show sourceTauChainR g
    = sourcePiR (sourceState0 g)
      * (sourcePR (sourceState0 g) (sourceState1 g)
        * sourcePR (sourceState1 g) (sourceState2 g))
  simp only [sourceTauChainR, sourcePiR, sourcePR]
  rw [← Rat.cast_mul, ← Rat.cast_mul]
  have hq : sourceTauChainQ g
      = mixingChainStationary (sourceState0 g)
        * (mixingChain (sourceState0 g) (sourceState1 g)
          * mixingChain (sourceState1 g) (sourceState2 g)) := by
    rw [← mul_assoc]
    exact sourceTauChainQ_eq_markov_product g
  exact_mod_cast hq

/-! ### Packet receipt -/

/-- **Derived-action packet receipt (B7, issue #683).**  On the committed
literals of the preregistered run: the tilt of the declared reference by
the log-transition action at multiplier one is the committed chain path law
with partition constant exactly `1/4`; an action-multiplier pair reproduces
that path law precisely when its multiplier-weighted action is the derived
action plus a constant; and the committed repair-count action reproduces it
at no multiplier.  The which-action slot of the B7 packet is closed at the
representation level by derivation, and the multiplier slot reduces to the
action-scale gauge fixing. -/
theorem logTransitionActionPacket_receipt :
    tilt sourceLogRefR sourceLogActionR 1 = sourceTauChainR
      ∧ tiltZ sourceLogRefR sourceLogActionR 1 = 1 / 4
      ∧ (∀ (S' : Fin 8 → ℝ) (lam' : ℝ),
          tilt sourceLogRefR S' lam' = sourceTauChainR
            ↔ ∃ c : ℝ, ∀ g, lam' * S' g = sourceLogActionR g + c)
      ∧ (∀ lam : ℝ, tilt sourceLogRefR sourceActionR lam ≠ sourceTauChainR) :=
  ⟨sourceLogTilt_eq_chain, sourceLogTiltZ_eq,
    fun S' lam' => sourceAction_unique_up_to_gauge S' lam',
    sourceRepairAction_no_multiplier⟩

end OPH.InformationProjection

#print axioms OPH.InformationProjection.tilt_eq_tilt_iff_gauge
#print axioms OPH.InformationProjection.tilt_action_add_const
#print axioms OPH.InformationProjection.tilt_action_multiplier_rescale
#print axioms OPH.InformationProjection.tilt_gauge_two_point
#print axioms OPH.InformationProjection.markovPathLaw_sum
#print axioms OPH.InformationProjection.markov_tiltZ_eq
#print axioms OPH.InformationProjection.markov_path_law_eq_gibbs
#print axioms OPH.InformationProjection.action_unique_up_to_gauge
#print axioms OPH.InformationProjection.markov_multiplier_action_tradeoff
#print axioms OPH.InformationProjection.sourceTransProdNum_eq_counts
#print axioms OPH.InformationProjection.sourceTransProdDen_eq_rowSums
#print axioms OPH.InformationProjection.sourceTauChain_gibbs_integer_identity
#print axioms OPH.InformationProjection.sourceTauChainQ_eq_pi_mul_transProd
#print axioms OPH.InformationProjection.sourceLogTiltZ_eq
#print axioms OPH.InformationProjection.sourceLogTilt_eq_chain
#print axioms OPH.InformationProjection.sourceAction_unique_up_to_gauge
#print axioms OPH.InformationProjection.sourceMultiplier_action_tradeoff
#print axioms OPH.InformationProjection.sourceRepairAction_no_multiplier
#print axioms OPH.InformationProjection.sourceRepairAction_not_gauge
#print axioms OPH.InformationProjection.sourceMarkov_eq_gibbs
#print axioms OPH.InformationProjection.sourceTauChainR_eq_markovPathLaw
#print axioms OPH.InformationProjection.logTransitionActionPacket_receipt
