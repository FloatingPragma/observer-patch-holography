import WhitneyTimeBridge
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# Full-metric time-prism Maxwell dynamics

The mass and stiffness pairings, cochain gradient, source covectors and time
step are supplied. Exact quadratic action variations give the temporal-gauge
recurrence; the gradient-null stiffness identity gives a weak Noether identity
and propagation, rather than repair, of the initial Gauss residual. The
conserved/work form is that of the exactly integrated prism action, not the
counting-metric or right-endpoint action.

The scalar modal analysis keeps positive stiffness, zero modes, the critical
Jordan mode, and endpoint boundary-value resonance separate. No source-selected
geometry, metric, physical clock, continuum limit, Born rule or quantization is
asserted by these finite algebraic statements.
-/

namespace OPH.WhitneyMaxwellDynamics

noncomputable section
open scoped BigOperators

variable {V P : Type*} [AddCommGroup V] [Module ℝ V]
  [AddCommGroup P] [Module ℝ P]

abbrev Pairing (V : Type*) [AddCommGroup V] [Module ℝ V] :=
  V →ₗ[ℝ] V →ₗ[ℝ] ℝ

/-- The temporal-gauge field-plus-source action of one prism slab. -/
def slabAction (h : ℝ) (M K : Pairing V) (J : V →ₗ[ℝ] ℝ) (x y : V) : ℝ :=
  M (y - x) (y - x) / (2 * h) -
    h / 6 * (K x x + K x y + K y y) + h * J y

/-- Exact middle-node first variation of two adjacent slabs. No stationarity
is assumed. The force at the later endpoint does not couple to this node. -/
theorem two_slab_action_expansion
    (h t : ℝ) (M K : Pairing V)
    (hM : ∀ x y, M x y = M y x) (hK : ∀ x y, K x y = K y x)
    (J₀ J₁ : V →ₗ[ℝ] ℝ) (x y z v : V) :
    slabAction h M K J₀ x (y + t • v) +
      slabAction h M K J₁ (y + t • v) z =
    slabAction h M K J₀ x y + slabAction h M K J₁ y z +
      t * (M ((2 : ℝ) • y - x - z) v / h -
        h / 6 * K (x + (4 : ℝ) • y + z) v + h * J₀ v) +
      t ^ 2 * (M v v / h - h / 3 * K v v) := by
  simp only [slabAction, map_add, map_sub, map_smul, LinearMap.add_apply,
    LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
  rw [hM v x, hM v y, hM v z, hK v y, hK v z]
  ring

/-- The scalar-potential variation gives the full mass-metric Gauss
functional, independently of the magnetic action. -/
theorem scalar_potential_action_expansion
    (h t : ℝ) (M : Pairing V) (hM : ∀ x y, M x y = M y x)
    (D : P →ₗ[ℝ] V) (ρ : P →ₗ[ℝ] ℝ) (E : V) (φ χ : P) :
    h / 2 * M (E - t • D χ) (E - t • D χ) + h * ρ (φ + t • χ) =
      h / 2 * M E E + h * ρ φ +
      t * h * (ρ χ - M E (D χ)) + t ^ 2 * h / 2 * M (D χ) (D χ) := by
  simp only [map_add, map_sub, map_smul, LinearMap.sub_apply,
    LinearMap.smul_apply, smul_eq_mul]
  rw [hM (D χ) E]
  ring

/-- The zero middle-node action derivative is the implicit prism recurrence
against that same test covector. Nonzero `h` is required for this equivalence. -/
theorem action_derivative_iff_recurrence
    (h : ℝ) (hh : h ≠ 0) (M K : Pairing V)
    (J : V →ₗ[ℝ] ℝ) (x y z v : V) :
    (M ((2 : ℝ) • y - x - z) v / h -
      h / 6 * K (x + (4 : ℝ) • y + z) v + h * J v = 0) ↔
    M (z - (2 : ℝ) • y + x) v +
      h ^ 2 / 6 * K (x + (4 : ℝ) • y + z) v = h ^ 2 * J v := by
  simp only [map_add, map_sub, map_smul, LinearMap.add_apply,
    LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
  constructor <;> intro hEq
  · field_simp at hEq
    nlinarith [hEq]
  · field_simp
    nlinarith [hEq]

/-- The left recurrence coefficient is `F=M+h²K/6`; this form is positive
on every nonzero vector if mass is positive and stiffness is nonnegative.
This does not claim invertibility on an infinite-dimensional space. -/
theorem effective_mass_positive
    (h : ℝ) (M K : Pairing V)
    (hM : ∀ v, v ≠ 0 → 0 < M v v) (hK : ∀ v, 0 ≤ K v v)
    (v : V) (hv : v ≠ 0) :
    0 < M v v + h ^ 2 / 6 * K v v := by
  have := hM v hv
  have := mul_nonneg (show 0 ≤ h ^ 2 / 6 by positivity) (hK v)
  linarith

/-- At most one next potential slice solves all the weak recurrence tests.
The positivity of the effective mass is sufficient; no step-size restriction
is needed for this uniqueness statement. -/
theorem recurrence_next_unique
    (h : ℝ) (M K : Pairing V)
    (hM : ∀ v, v ≠ 0 → 0 < M v v) (hK : ∀ v, 0 ≤ K v v)
    (J : V →ₗ[ℝ] ℝ) (x y z₁ z₂ : V)
    (h₁ : ∀ v, M (z₁ - (2 : ℝ) • y + x) v +
      h ^ 2 / 6 * K (x + (4 : ℝ) • y + z₁) v = h ^ 2 * J v)
    (h₂ : ∀ v, M (z₂ - (2 : ℝ) • y + x) v +
      h ^ 2 / 6 * K (x + (4 : ℝ) • y + z₂) v = h ^ 2 * J v) : z₁ = z₂ := by
  by_contra hne
  have hp := effective_mass_positive h M K hM hK (z₁ - z₂) (sub_ne_zero.mpr hne)
  have heq₁ := h₁ (z₁ - z₂)
  have heq₂ := h₂ (z₁ - z₂)
  simp only [map_add, map_sub, map_smul, LinearMap.add_apply,
    LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul] at hp heq₁ heq₂
  nlinarith

/-- A supplied quadratic spectral bound makes the corrected electric form
positive below its sharp sufficient step bound. This theorem does not infer
the bound from a sampled numerical eigensystem. -/
theorem corrected_kinetic_positive
    (h Λ : ℝ) (M K : Pairing V)
    (hM : ∀ v, v ≠ 0 → 0 < M v v)
    (hbound : ∀ v, K v v ≤ Λ * M v v) (hstep : h ^ 2 * Λ < 12)
    (v : V) (hv : v ≠ 0) : 0 < M v v - h ^ 2 / 12 * K v v := by
  have hm := hM v hv
  have hc : 0 < 1 - h ^ 2 * Λ / 12 := by linarith
  have hp := mul_pos hc hm
  have hb := mul_le_mul_of_nonneg_left (hbound v) (show 0 ≤ h ^ 2 / 12 by positivity)
  nlinarith

/-- Full weak Ampere functional, with the fourfold middle magnetic sample. -/
def ampereWeak (h : ℝ) (M K : Pairing V) (J : V →ₗ[ℝ] ℝ)
    (E₀ E₁ x y z v : V) : ℝ :=
  M (E₁ - E₀) v - h / 6 * K (x + (4 : ℝ) • y + z) v + h * J v

/-- With electric cochains obtained from the temporal-gauge potential,
Ampere is exactly the middle-node action coefficient above. -/
theorem temporal_gauge_ampere_is_action_derivative
    (h : ℝ) (M K : Pairing V) (J : V →ₗ[ℝ] ℝ) (x y z v : V) :
    ampereWeak h M K J (-(h⁻¹ • (y - x))) (-(h⁻¹ • (z - y))) x y z v =
      M ((2 : ℝ) • y - x - z) v / h -
        h / 6 * K (x + (4 : ℝ) • y + z) v + h * J v := by
  simp only [ampereWeak, map_add, map_sub, map_neg, map_smul,
    LinearMap.add_apply, LinearMap.sub_apply, LinearMap.neg_apply,
    LinearMap.smul_apply, smul_eq_mul, div_eq_mul_inv]
  ring

/-- This join lets the recurrence feed the weak Noether propagation theorem
using the same potential history, without a separately supplied Ampere law. -/
theorem temporal_gauge_ampere_iff_recurrence
    (h : ℝ) (hh : h ≠ 0) (M K : Pairing V)
    (J : V →ₗ[ℝ] ℝ) (x y z v : V) :
    (ampereWeak h M K J (-(h⁻¹ • (y - x))) (-(h⁻¹ • (z - y))) x y z v = 0) ↔
      M (z - (2 : ℝ) • y + x) v +
        h ^ 2 / 6 * K (x + (4 : ℝ) • y + z) v = h ^ 2 * J v := by
  rw [temporal_gauge_ampere_is_action_derivative]
  exact action_derivative_iff_recurrence h hh M K J x y z v

def gaussWeak (M : Pairing V) (D : P →ₗ[ℝ] V)
    (E : V) (ρ : P →ₗ[ℝ] ℝ) (χ : P) : ℝ := M E (D χ) - ρ χ

def continuityWeak (h : ℝ) (D : P →ₗ[ℝ] V)
    (ρ₀ ρ₁ : P →ₗ[ℝ] ℝ) (J : V →ₗ[ℝ] ℝ) (χ : P) : ℝ :=
  ρ₁ χ - ρ₀ χ + h * J (D χ)

/-- For a stiffness obtained by pulling back a magnetic pairing through the
cochain curl, the chain identity supplies the Noether cancellation. No metric
identity between different cochain degrees is used. -/
theorem curl_stiffness_gradient_null
    {W : Type*} [AddCommGroup W] [Module ℝ W]
    (D : P →ₗ[ℝ] V) (C : V →ₗ[ℝ] W) (G : Pairing W) (K : Pairing V)
    (hK : ∀ v w, K v w = G (C v) (C w))
    (hCD : ∀ χ, C (D χ) = 0) (v : V) (χ : P) : K v (D χ) = 0 := by
  rw [hK, hCD, map_zero]

/-- Discrete weak Noether identity. The stiffness-gradient cancellation is
the only cochain hypothesis; physical divergence is not identified here. -/
theorem noether_gauss_identity
    (h : ℝ) (M K : Pairing V) (D : P →ₗ[ℝ] V)
    (hKD : ∀ v χ, K v (D χ) = 0)
    (ρ₀ ρ₁ : P →ₗ[ℝ] ℝ) (J : V →ₗ[ℝ] ℝ) (E₀ E₁ x y z : V) (χ : P) :
    gaussWeak M D E₁ ρ₁ χ - gaussWeak M D E₀ ρ₀ χ =
      ampereWeak h M K J E₀ E₁ x y z (D χ) - continuityWeak h D ρ₀ ρ₁ J χ := by
  simp only [gaussWeak, ampereWeak, continuityWeak, hKD, map_sub,
    LinearMap.sub_apply]
  ring

/-- Conserved sources and Ampere preserve each initial Gauss residual;
they do not set a nonzero residual to zero. -/
theorem gauss_residual_preserved
    (h : ℝ) (M K : Pairing V) (D : P →ₗ[ℝ] V)
    (hKD : ∀ v χ, K v (D χ) = 0)
    (A E : ℕ → V) (ρ : ℕ → P →ₗ[ℝ] ℝ) (J : ℕ → V →ₗ[ℝ] ℝ)
    (hAmp : ∀ n v, ampereWeak h M K (J n) (E n) (E (n + 1))
      (A n) (A (n + 1)) (A (n + 2)) v = 0)
    (hCon : ∀ n χ, continuityWeak h D (ρ n) (ρ (n + 1)) (J n) χ = 0)
    (n : ℕ) (χ : P) : gaussWeak M D (E n) (ρ n) χ = gaussWeak M D (E 0) (ρ 0) χ := by
  induction n with
  | zero => rfl
  | succ n ih =>
      have hn := noether_gauss_identity h M K D hKD (ρ n) (ρ (n + 1)) (J n)
        (E n) (E (n + 1)) (A n) (A (n + 1)) (A (n + 2)) χ
      rw [hAmp n (D χ), hCon n χ] at hn
      linarith

/-- Prism invariant expressed using adjacent potential slices. -/
def pairEnergy (h : ℝ) (M K : Pairing V) (x y : V) : ℝ :=
  M (y - x) (y - x) / (2 * h ^ 2) +
    (K x x + 4 * K x y + K y y) / 12

/-- Centered magnetic form plus the corrected electric kinetic form.
The coefficient is `1-h² lambda/12` on a mass-normalized positive mode. -/
theorem pairEnergy_centered
    (h : ℝ) (hh : h ≠ 0) (M K : Pairing V)
    (hK : ∀ x y, K x y = K y x) (x y : V) :
    pairEnergy h M K x y =
      (M (h⁻¹ • (y - x)) (h⁻¹ • (y - x)) -
        h ^ 2 / 12 * K (h⁻¹ • (y - x)) (h⁻¹ • (y - x))) / 2 +
      K (x + y) (x + y) / 8 := by
  simp only [pairEnergy, map_add, map_sub, map_smul, LinearMap.add_apply,
    LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
  rw [hK y x]
  field_simp
  ring

/-- Exact polarization of the invariant difference before imposing motion. -/
theorem pairEnergy_difference
    (h : ℝ) (M K : Pairing V)
    (hM : ∀ x y, M x y = M y x) (hK : ∀ x y, K x y = K y x)
    (x y z : V) :
    pairEnergy h M K y z - pairEnergy h M K x y =
      M (z - (2 : ℝ) • y + x) (z - x) / (2 * h ^ 2) +
      K (x + (4 : ℝ) • y + z) (z - x) / 12 := by
  simp only [pairEnergy, map_add, map_sub, map_smul, LinearMap.add_apply,
    LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
  rw [hM z y, hM y x, hM z x, hK y x, hK z x]
  ring

/-- The sourced prism recurrence has the exact work law. -/
theorem pairEnergy_work
    (h : ℝ) (hh : h ≠ 0) (M K : Pairing V)
    (hM : ∀ x y, M x y = M y x) (hK : ∀ x y, K x y = K y x)
    (J : V →ₗ[ℝ] ℝ) (x y z : V)
    (hRec : ∀ v, M (z - (2 : ℝ) • y + x) v +
      h ^ 2 / 6 * K (x + (4 : ℝ) • y + z) v = h ^ 2 * J v) :
    pairEnergy h M K y z - pairEnergy h M K x y = J (z - x) / 2 := by
  rw [pairEnergy_difference h M K hM hK]
  have heq := hRec (z - x)
  field_simp
  nlinarith [heq]

/-- The sign of work in electric variables is minus the average electric
field paired with the supplied current. -/
theorem work_electric_sign
    (h : ℝ) (hh : h ≠ 0) (J : V →ₗ[ℝ] ℝ) (x y z : V) :
    J (z - x) / 2 =
      -(h / 2) * J (-(h⁻¹ • (z - y)) - h⁻¹ • (y - x)) := by
  simp only [map_sub, map_neg, map_smul, smul_eq_mul]
  field_simp
  ring

/-- Every source-free prism solution conserves its adjacent-slice form. -/
theorem pairEnergy_conserved
    (h : ℝ) (hh : h ≠ 0) (M K : Pairing V)
    (hM : ∀ x y, M x y = M y x) (hK : ∀ x y, K x y = K y x)
    (A : ℕ → V)
    (hRec : ∀ n v, M (A (n + 2) - (2 : ℝ) • A (n + 1) + A n) v +
      h ^ 2 / 6 * K (A n + (4 : ℝ) • A (n + 1) + A (n + 2)) v = 0)
    (n : ℕ) : pairEnergy h M K (A n) (A (n + 1)) = pairEnergy h M K (A 0) (A 1) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      have hn := pairEnergy_work h hh M K hM hK 0 (A n) (A (n + 1)) (A (n + 2))
        (by simpa using hRec n)
      simp only [LinearMap.zero_apply, zero_div] at hn
      change pairEnergy h M K (A (n + 1)) (A (n + 2)) = _
      linarith

/-! ## Positive scalar modes, the sharp boundary, and endpoint resonance -/

/-- Dimensionless modal parameter is `z=h² lambda`, with mass normalized. -/
def modalEnergy (z x y : ℝ) : ℝ :=
  (y - x) ^ 2 / 2 + z / 12 * (x ^ 2 + 4 * x * y + y ^ 2)

theorem modalEnergy_centered (z x y : ℝ) :
    modalEnergy z x y = (1 - z / 12) * (y - x) ^ 2 / 2 + z * (x + y) ^ 2 / 8 := by
  unfold modalEnergy
  ring

theorem modalEnergy_conserved_step (z x y w : ℝ)
    (hRec : (1 + z / 6) * (w + x) = (2 - 2 * z / 3) * y) :
    modalEnergy z y w = modalEnergy z x y := by
  unfold modalEnergy
  have heq := congrArg (fun q : ℝ => q * (w - x)) hRec
  nlinarith [heq]

/-- Strict coercivity on the two-dimensional state of a positive mode. -/
theorem modalEnergy_positive_iff (z : ℝ) :
    (∀ x y : ℝ, x ≠ 0 ∨ y ≠ 0 → 0 < modalEnergy z x y) ↔ 0 < z ∧ z < 12 := by
  constructor
  · intro h
    have hdiag := h 1 1 (Or.inl one_ne_zero)
    have hanti := h 1 (-1) (Or.inl one_ne_zero)
    simp only [modalEnergy] at hdiag hanti
    constructor <;> nlinarith
  · rintro ⟨hz, hz12⟩ x y hxy
    rw [modalEnergy_centered]
    have ha : 0 < 1 - z / 12 := by linarith
    have hn₁ := mul_nonneg ha.le (sq_nonneg (y - x))
    have hn₂ := mul_nonneg hz.le (sq_nonneg (x + y))
    by_cases hdiff : y - x = 0
    · have hsum : x + y ≠ 0 := by
        intro heq
        rcases hxy with hx | hy <;> apply_assumption <;> linarith
      have := sq_pos_of_ne_zero hsum
      nlinarith [mul_pos hz this]
    · have := sq_pos_of_ne_zero hdiff
      nlinarith [mul_pos ha this]

/-- An explicit positive coercivity constant for the full two-slice state. -/
theorem modalEnergy_coercive (z x y : ℝ) (hz : 0 < z) (hz12 : z < 12) :
    0 < min (1 - z / 12) (z / 4) ∧
    min (1 - z / 12) (z / 4) * (x ^ 2 + y ^ 2) ≤ modalEnergy z x y := by
  constructor
  · exact lt_min (by linarith) (by positivity)
  · rw [modalEnergy_centered]
    have h₁ := mul_nonneg
      (sub_nonneg.mpr (min_le_left (1 - z / 12) (z / 4))) (sq_nonneg (y - x))
    have h₂ := mul_nonneg
      (sub_nonneg.mpr (min_le_right (1 - z / 12) (z / 4))) (sq_nonneg (x + y))
    nlinarith

/-- Every solution in the positive modal window has a uniform state bound,
with no unstated diagonalizability hypothesis on a larger matrix system. -/
theorem modal_solution_bound (z : ℝ) (hz : 0 < z) (hz12 : z < 12)
    (a : ℕ → ℝ)
    (hRec : ∀ n, (1 + z / 6) * (a (n + 2) + a n) =
      (2 - 2 * z / 3) * a (n + 1)) (n : ℕ) :
    (a n) ^ 2 + (a (n + 1)) ^ 2 ≤
      modalEnergy z (a 0) (a 1) / min (1 - z / 12) (z / 4) := by
  have hCons : ∀ m, modalEnergy z (a m) (a (m + 1)) = modalEnergy z (a 0) (a 1) := by
    intro m
    induction m with
    | zero => rfl
    | succ m ih =>
        calc
          _ = modalEnergy z (a m) (a (m + 1)) :=
            modalEnergy_conserved_step z (a m) (a (m + 1)) (a (m + 2)) (hRec m)
          _ = _ := ih
  obtain ⟨hc, hbound⟩ := modalEnergy_coercive z (a n) (a (n + 1)) hz hz12
  rw [hCons n] at hbound
  exact (le_div_iff₀ hc).2 (by simpa only [mul_comm] using hbound)

/-- The numerator of the characteristic discriminant changes sign exactly
at zero and twelve. This is distinct from endpoint boundary-value resonance. -/
theorem characteristic_discriminant (z : ℝ) :
    (2 - 2 * z / 3) ^ 2 - 4 * (1 + z / 6) ^ 2 = z * (z - 12) / 3 := by
  ring

/-- Every parameter above the critical value has a real expanding root.
The equality case is handled separately by the Jordan solution below. -/
theorem supercritical_expanding_root (z : ℝ) (hz : 12 < z) :
    ∃ r : ℝ, r < -1 ∧
      (1 + z / 6) * r ^ 2 - (2 - 2 * z / 3) * r + (1 + z / 6) = 0 := by
  let d : ℝ := 1 + z / 6
  let b : ℝ := 2 - 2 * z / 3
  let α : ℝ := b / d
  have hd : 0 < d := by dsimp [d]; linarith
  have hα : α < -2 := by
    dsimp [α]
    apply (div_lt_iff₀ hd).2
    dsimp [b, d]
    linarith
  have hdisc : 0 ≤ α ^ 2 - 4 := by nlinarith
  have hs := Real.sq_sqrt hdisc
  have hs0 := Real.sqrt_nonneg (α ^ 2 - 4)
  have hmul : d * α = b := by
    dsimp [α]
    field_simp
  let r : ℝ := (α - Real.sqrt (α ^ 2 - 4)) / 2
  have hr : r < -1 := by dsimp [r]; linarith
  have hroot : r ^ 2 - α * r + 1 = 0 := by dsimp [r]; nlinarith [hs]
  refine ⟨r, hr, ?_⟩
  change d * r ^ 2 - b * r + d = 0
  rw [← hmul]
  linear_combination d * hroot

/-- The expanding root gives an actual recurrence solution whose squared
potential grows without bound. Thus strict modal coercivity is also a sharp
dynamical boundary above zero, not merely a sufficient energy test. -/
theorem supercritical_unbounded_solution (z : ℝ) (hz : 12 < z) :
    ∃ r : ℝ,
      (∀ n : ℕ, (1 + z / 6) * (r ^ (n + 2) + r ^ n) =
        (2 - 2 * z / 3) * r ^ (n + 1)) ∧
      Filter.Tendsto (fun n : ℕ => (r ^ n) ^ 2) Filter.atTop Filter.atTop := by
  obtain ⟨r, hr, hroot⟩ := supercritical_expanding_root z hz
  refine ⟨r, ?_, ?_⟩
  · intro n
    rw [pow_add, pow_add]
    simp only [pow_one]
    linear_combination (r ^ n) * hroot
  · have hr2 : 1 < r ^ 2 := by nlinarith
    convert tendsto_pow_atTop_atTop_of_one_lt hr2 using 1
    ext n
    rw [← pow_mul, ← pow_mul, Nat.mul_comm]

/-- At the critical parameter the generalized `-1` mode solves the recurrence. -/
def criticalMode (n : ℕ) : ℝ := (-1 : ℝ) ^ n * n

theorem criticalMode_recurrence (n : ℕ) :
    criticalMode (n + 2) + 2 * criticalMode (n + 1) + criticalMode n = 0 := by
  simp only [criticalMode, pow_add, Nat.cast_add, Nat.cast_ofNat]
  ring

theorem criticalMode_difference (n : ℕ) :
    (criticalMode (n + 1) - criticalMode n) ^ 2 = (2 * (n : ℝ) + 1) ^ 2 := by
  simp only [criticalMode, pow_add, Nat.cast_add, Nat.cast_one, pow_one]
  have hp : ((-1 : ℝ) ^ n) ^ 2 = 1 := by rw [← pow_mul, mul_comm n 2, pow_mul]; simp
  calc
    _ = ((-1 : ℝ) ^ n) ^ 2 * (2 * (n : ℝ) + 1) ^ 2 := by ring
    _ = _ := by rw [hp, one_mul]

/-- At the boundary the electric difference is itself unbounded; the
critical growth is not merely an unobservable constant-potential shift. -/
theorem criticalMode_unbounded_difference :
    Filter.Tendsto (fun n : ℕ => (criticalMode (n + 1) - criticalMode n) ^ 2)
      Filter.atTop Filter.atTop := by
  apply Filter.tendsto_atTop_mono (f := fun n : ℕ => (n : ℝ))
  · intro n
    rw [criticalMode_difference]
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    nlinarith [sq_nonneg (2 * (n : ℝ))]
  · exact tendsto_natCast_atTop_atTop

/-- A zero stiffness mode may drift in potential while its electric
difference remains constant; it is excluded from positive-mode coercivity. -/
theorem zeroMode_drift (a b : ℝ) (n : ℕ) :
    ((a + (n + 2 : ℕ) * b) - 2 * (a + (n + 1 : ℕ) * b) +
      (a + (n : ℝ) * b) = 0) ∧
    (a + (n + 1 : ℕ) * b) - (a + (n : ℝ) * b) = b := by
  push_cast
  constructor <;> ring

/-- Away from the two-slab modal resonance, fixed endpoints and one impulse
give exactly one middle value. This is separate from the IVP stability test. -/
theorem two_slab_endpoint_solution
    (σ x y z f : ℝ) (hσ : σ ≠ 3) :
    ((1 + σ / 6) * (z + x) = (2 - 2 * σ / 3) * y + f) ↔
      y = ((1 + σ / 6) * (z + x) - f) / (2 - 2 * σ / 3) := by
  have hd : 2 - 2 * σ / 3 ≠ 0 := by intro heq; apply hσ; linarith
  rw [eq_div_iff hd]
  constructor <;> intro heq <;> nlinarith

/-- At resonance only endpoint/source compatibility remains, leaving the
middle value free when compatible. -/
theorem two_slab_resonance_compatibility (x y z f : ℝ) :
    ((1 + (3 : ℝ) / 6) * (z + x) = (2 - 2 * 3 / 3) * y + f) ↔
      f = 3 / 2 * (z + x) := by
  norm_num
  constructor <;> intro heq <;> linarith

/-- At `z=3`, zero endpoint data leave the middle scalar amplitude free
for zero impulse and permit no solution for a nonzero impulse. The IVP
is inside the strict modal coercivity window. -/
theorem two_slab_endpoint_resonance :
    (0 : ℝ) < 3 ∧ (3 : ℝ) < 12 ∧
    (∀ y : ℝ, (1 + (3 : ℝ) / 6) * (0 + 0) = (2 - 2 * 3 / 3) * y) ∧
    (∀ f : ℝ, f ≠ 0 → ¬ ∃ y : ℝ,
      (1 + (3 : ℝ) / 6) * (0 + 0) = (2 - 2 * 3 / 3) * y + f) := by
  norm_num
  exact fun f hf heq => hf heq.symm

end
end OPH.WhitneyMaxwellDynamics

#print axioms OPH.WhitneyMaxwellDynamics.two_slab_action_expansion
#print axioms OPH.WhitneyMaxwellDynamics.action_derivative_iff_recurrence
#print axioms OPH.WhitneyMaxwellDynamics.recurrence_next_unique
#print axioms OPH.WhitneyMaxwellDynamics.corrected_kinetic_positive
#print axioms OPH.WhitneyMaxwellDynamics.temporal_gauge_ampere_iff_recurrence
#print axioms OPH.WhitneyMaxwellDynamics.noether_gauss_identity
#print axioms OPH.WhitneyMaxwellDynamics.gauss_residual_preserved
#print axioms OPH.WhitneyMaxwellDynamics.pairEnergy_work
#print axioms OPH.WhitneyMaxwellDynamics.pairEnergy_conserved
#print axioms OPH.WhitneyMaxwellDynamics.modalEnergy_positive_iff
#print axioms OPH.WhitneyMaxwellDynamics.modal_solution_bound
#print axioms OPH.WhitneyMaxwellDynamics.supercritical_unbounded_solution
#print axioms OPH.WhitneyMaxwellDynamics.criticalMode_difference
#print axioms OPH.WhitneyMaxwellDynamics.criticalMode_unbounded_difference
#print axioms OPH.WhitneyMaxwellDynamics.two_slab_endpoint_solution
#print axioms OPH.WhitneyMaxwellDynamics.two_slab_endpoint_resonance
