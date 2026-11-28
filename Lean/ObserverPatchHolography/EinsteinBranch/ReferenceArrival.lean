import ObserverPatchHolography.EinsteinBranch.ContactCensus
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-! Finite comparison of a named observer record against independent reference
records. The length-k prefix has no match with probability (1-p)^k. -/
namespace OPH.EinsteinBranch.ReferenceArrival

open Finset

noncomputable def cost (p : ℝ) (n : ℕ) :=
  ∑ i ∈ range n, (1-p)^(i+1) / (i+1 : ℝ)

/-- The expected number of new rank minima before the first match sums to
the modular cost of the observed word. Ranks are independent of the matches. -/
theorem recurrence_cost (p : ℝ) (hp : 0 < p) (hp1 : p ≤ 1) :
    HasSum (fun i : ℕ => (1-p)^(i+1)/(i+1 : ℝ)) (-Real.log p) := by
  have h : |1-p| < 1 := by rw [abs_of_nonneg (by linarith)]; linarith
  simpa using Real.hasSum_pow_div_log_of_abs_lt_one h

/-- Uniform finite stopping error, with the denominator from the first
omitted comparison retained. -/
theorem finite_cost (p : ℝ) (hp : 0 < p) (hp1 : p ≤ 1) (n : ℕ) :
    0 ≤ -Real.log p-cost p n ∧
      -Real.log p-cost p n ≤ (1-p)^(n+1)/((n+1 : ℝ)*p) := by
  have hn : 0 ≤ 1-p := by linarith
  have htail := (hasSum_nat_add_iff' n).mpr (recurrence_cost p hp hp1)
  change HasSum (fun i : ℕ => (1-p)^(i+n+1)/((i+n : ℕ)+1 : ℝ))
    (-Real.log p - cost p n) at htail
  constructor
  · exact htail.nonneg (fun j => by positivity)
  · have hg := (hasSum_geometric_of_lt_one hn (by linarith : 1-p<1)).mul_left
        ((1-p)^(n+1)/(n+1 : ℝ))
    have hb := hasSum_le (g := fun i : ℕ => ((1-p)^(n+1)/(n+1 : ℝ))*(1-p)^i)
      (fun i => (show (1-p)^(i+n+1)/((i+n : ℕ)+1 : ℝ) ≤
          (1-p)^(n+1)/(n+1 : ℝ)*(1-p)^i from by
        rw [show i+n+1 = (n+1)+i by omega, pow_add, div_mul_eq_mul_div]
        exact div_le_div_of_nonneg_left (by positivity) (by positivity)
          (by push_cast; linarith))) htail hg
    calc
      _ ≤ _ := hb
      _ = _ := by
        rw [show 1-(1-p)=p by ring, ← div_eq_mul_inv, div_div]

/-- With errors between zero and e on every sharp record, any contrast
of normalized record preparations inherits the same bound. -/
theorem all_preparations {n : ℕ} (u v err : Fin n → ℝ) (e : ℝ)
    (hu : ∀ i, 0 ≤ u i) (hv : ∀ i, 0 ≤ v i)
    (su : ∑ i, u i = 1) (sv : ∑ i, v i = 1)
    (lo : ∀ i, 0 ≤ err i) (hi : ∀ i, err i ≤ e) :
    |(∑ i, u i*err i) - (∑ i, v i*err i)| ≤ e := by
  have avg (w : Fin n → ℝ) (hw : ∀ i, 0 ≤ w i) (hs : ∑ i, w i=1) :
      0 ≤ ∑ i, w i*err i ∧ ∑ i, w i*err i ≤ e := by
    constructor
    · exact Finset.sum_nonneg fun i _ => mul_nonneg (hw i) (lo i)
    · calc
        _ ≤ ∑ i, w i*e := Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (hi i) (hw i)
        _ = e := by rw [← Finset.sum_mul, hs, one_mul]
  obtain ⟨h₁,h₂⟩ := avg u hu su
  obtain ⟨h₃,h₄⟩ := avg v hv sv
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- Every downstream event conditioned only on receipt has the same
transfer ratio, including an event about the complete archived history. -/
theorem receiver_contrast (p q yes no : ℝ) :
   (q*yes+(1-q)*no)-(p*yes+(1-p)*no) = (q-p)*(yes-no) := by ring

/-- Changing only the admission probability bounds the change in any
enclosed ordinary population by the absolute admission change. -/
theorem ordinary_stock (p q yes no : ℝ) (y0 : 0 ≤ yes) (y1 : yes ≤ 1)
    (n0 : 0 ≤ no) (n1 : no ≤ 1) :
    |(q-p)*(yes-no)| ≤ |q-p| := by
  rw [abs_mul]
  calc
    _ ≤ |q-p| * 1 := mul_le_mul_of_nonneg_left (abs_le.mpr ⟨by linarith, by linarith⟩) (abs_nonneg _)
    _ = _ := mul_one _

/-- Any finite orbit satisfying this linear recurrence stays inside the
quiet span. A detector orthogonal to the span gives zero on every delay
whose propagator acts in that span. -/
theorem quiet_orbit {V : Type*} [AddCommGroup V] [Module ℝ V] (K : Submodule ℝ V)
    (h : V →ₗ[ℝ] V) (invariant : ∀ v ∈ K, h v ∈ K) (x : V) (start : x ∈ K)
    (m : V →ₗ[ℝ] ℝ) (zero : ∀ v ∈ K, m v = 0) (n : ℕ) :
    m ((h^[n]) x) = 0 := by
  apply zero
  induction n with
  | zero => exact start
  | succ n ih =>
    rw [Function.iterate_succ_apply']
    exact invariant _ ih

/-- At the central phase flip every incident hopping term changes by twice
its measured real coherence. The first cut reverses and the outer cut stays
fixed. Their equally weighted coherent current loses half the injected energy. -/
theorem phase_energy_flux {n : ℕ} (c : Fin n → ℝ) (outer : ℝ) :
    (∑ i, 2*c i) = -2*(((-(∑ i, c i))+outer)/2-((∑ i, c i)+outer)/2) := by
  rw [← Finset.mul_sum]
  ring

/-- A lower bound on the joined entropy slope and a positive calibrated
area coefficient impose this upper bound on the fixed-volume area slope. -/
theorem required_contraction (η area entropy witness : ℝ) (pos : 0 < η)
    (measure : witness ≤ entropy) (stationary : η*area+entropy=0) :
    area ≤ -witness/η := by
  apply (le_div_iff₀ pos).mpr
  linarith [mul_comm area η]

/-- The public choice of outer radius giving the prescribed mean volume. -/
theorem fixed_volume_coin (v₀ v₁ c target : ℝ) (hc : c ≠ 0) (hv : v₁-v₀ ≠ 0) :
    let b := (target/c-v₀)/(v₁-v₀)
    ((1-b)*v₀+b*v₁)*c = target := by
  dsimp
  field_simp
  grind

/-- For two sharp radius choices, removing the volume residual with the
measured slope gives the observable at the target coin exactly. -/
theorem equal_volume_read (v₀ v₁ y₀ y₁ b target : ℝ) (hv : v₁-v₀ ≠ 0) :
    let v := fun t : ℝ => (1-t)*v₀+t*v₁
    let y := fun t : ℝ => (1-t)*y₀+t*y₁
    y b - (v b-v target)*((y₁-y₀)/(v₁-v₀)) = y target := by
  dsimp
  field_simp
  ring

/-- Two perturbations of the same region satisfying the same entropy
balance have a vanishing area-entropy determinant. -/
theorem common_stationarity (η a b s t : ℝ) (h₁ : η*a+s=0) (h₂ : η*b+t=0) :
    a*t-b*s=0 := by
  linear_combination a*h₂-b*h₁

/-- Positive and negative changes of enclosed ordinary mass admit a genuine
randomization with exactly zero change. Every other mean uses this same coin. -/
theorem mass_neutral_coin (w d : ℝ) (hw : 0 < w) (hd : d < 0) :
    let t := w/(w-d)
    0 < t ∧ t < 1 ∧ (1-t)*w+t*d = 0 := by
  have h : 0 < w-d := by linarith
  dsimp
  refine ⟨div_pos hw h, (div_lt_one h).mpr (by linarith), ?_⟩
  field_simp
  ring

/-- For the marked center and its neighbors, a decrease of shell occupancy
by half the increase at the center adds enclosed mass and increases the
clocked volume. The shell alone changes the center clock. -/
theorem growing_occupied_region (center shell cells r : ℝ)
    (hc : 0 < center) (hV : 0 < cells) (hr : 0 ≤ r) (hr1 : r < 1)
    (mix : shell+center/2=0) :
    0 < center+shell ∧ 0 < cells*(r^3-1)*shell := by
  have hs : shell<0 := by linarith
  have hr3 : r^3 < 1 := by simpa using pow_lt_one₀ hr hr1 (by decide : (3:ℕ)≠0)
  exact ⟨by linarith, mul_pos_of_neg_of_neg (mul_neg_of_pos_of_neg hV (by linarith)) hs⟩

/-- Positive shell entry contracts the same marked volume at the same clock. -/
theorem shrinking_occupied_region (center shell cells r : ℝ)
    (hc : 0 ≤ center) (hs : 0 < shell) (hV : 0 < cells) (hr : 0 ≤ r) (hr1 : r < 1) :
    0 < center+shell ∧ cells*(r^3-1)*shell < 0 := by
  have hr3 : r^3 < 1 := by simpa using pow_lt_one₀ hr hr1 (by decide : (3:ℕ)≠0)
  exact ⟨by linarith, mul_neg_of_neg_of_pos (mul_neg_of_pos_of_neg hV (by linarith)) hs⟩

/-- Both observers use the same two sharp clock readings. Converting the
mean frequency by the secant of the cube equals the mean of cubed rates. -/
theorem two_rate_clock (r p : ℝ) :
    1+(1+r+r^2)*(((1-p)+p*r)-1) = (1-p)+p*r^3 := by ring

/-- A single coincidence of independently encoded local bits gives a
rank-one interaction after quantities confined to either station cancel. -/
theorem single_coincidence_minor (a₀ a₁ a₂ b₀ b₁ b₂ gain : ℝ) :
    let d := fun a b : ℝ => gain*(a*b-a*b₀-a₀*b+a₀*b₀)
    d a₁ b₁ * d a₂ b₂ - d a₁ b₂ * d a₂ b₁ = 0 := by
  dsimp
  ring

end OPH.EinsteinBranch.ReferenceArrival
