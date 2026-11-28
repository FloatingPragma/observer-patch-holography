import ObserverPatchHolography.EinsteinBranch.TrailAperture
import Mathlib.Analysis.Convex.Deriv
import Mathlib.LinearAlgebra.Matrix.Rank

/-! The radial weight in the flat rest-frame ball stress integral. Distances,
nonnegative stress increments, and the common physical conversion are hypotheses. -/
namespace OPH.EinsteinBranch.SourceApertures

/-- An independent positive read of the connected charge forces negative
survey-area response when it exceeds the certified stationarity and volume
error, with the positive calibration fixed by the parked ordinary source. -/
theorem area_from_source (q area gain err lo : ℝ)
    (hgain : 0 < gain) (measured : lo ≤ q) (resolved : err < lo)
    (stationarity : |q+gain*area| ≤ err) : area < 0 := by
  have bound := (abs_le.mp stationarity).2
  by_contra! nonneg
  linarith [mul_nonneg (le_of_lt hgain) nonneg]

/-- Eliminate the area and volume gains using two independent calibration
directions, also when the local record statistic responds to both. -/
theorem calibrated_record_transfer (s z q η v ao vo az vz ad vd : ℝ)
    (ordinary : s + η * ao - v * vo = 0)
    (rotation : z + η * az - v * vz = 0)
    (held : q + η * ad - v * vd = 0) :
    q * (ao * vz - az * vo) = s * (ad * vz - az * vd) + z * (ao * vd - ad * vo) := by
  linear_combination (ao * vz - az * vo) * held - (ad * vz - az * vd) * ordinary +
    (ad * vo - ao * vd) * rotation

/-- A volume-adjusting entropy-preserving rotation and an independent
ordinary intervention determine a necessary held-source response for each
survey, even when every survey has its own area gain and volume multiplier. -/
theorem calibrated_survey_transfer (s q η volumeMultiplier ao vo az vz ad vd : ℝ)
    (ordinary : s + η * ao - volumeMultiplier * vo = 0)
    (rotation : η * az - volumeMultiplier * vz = 0)
    (held : q + η * ad - volumeMultiplier * vd = 0) :
    q * (ao * vz - az * vo) = s * (ad * vz - az * vd) := by
  simpa using calibrated_record_transfer s 0 q η volumeMultiplier ao vo az vz ad vd
    ordinary (by simpa using rotation) held

noncomputable def weight (R r : ℝ) := max 0 ((R^2-r^2)/(2*R))

/-- Doubling the radius at least doubles every positive shell's weight,
including shells crossing the old boundary. -/
theorem weight_double (R r : ℝ) (hR : 0 < R) :
    2 * weight R r ≤ weight (2*R) r := by
  unfold weight
  by_cases h : 0 ≤ (R^2-r^2)/(2*R)
  · rw [max_eq_right h]
    apply le_trans _ (le_max_right _ _)
    field_simp
    nlinarith [sq_nonneg r]
  · rw [max_eq_left (le_of_not_ge h), mul_zero]
    exact le_max_left _ _

/-- No assumption about which shells contain the stress is needed. -/
theorem ball_double {ι : Type*} (s : Finset ι) (R : ℝ) (r m : ι → ℝ)
    (hR : 0 < R) (hm : ∀ i ∈ s, 0 ≤ m i) :
    2 * ∑ i ∈ s, weight R (r i) * m i ≤
      ∑ i ∈ s, weight (2*R) (r i) * m i := by
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro i hi
  rw [← mul_assoc]
  exact mul_le_mul_of_nonneg_right (weight_double R (r i) hR) (hm i hi)

/-- A positive saturated charge forces at least a third of that charge
as error at one radius, in the same modular units. -/
theorem saturation_error (Q a b ε : ℝ) (growth : 2*a ≤ b)
    (small : |Q-a| ≤ ε) (large : |Q-b| ≤ ε) : Q/3 ≤ ε := by
  obtain ⟨_, hu⟩ := abs_le.mp small
  obtain ⟨hl, _⟩ := abs_le.mp large
  linarith

/-- Given genuine diamond-weighted readings from the same mass shells,
two radii outside them recover the mass independently of its second moment. -/
theorem exterior_inversion (R S mass moment bR bS : ℝ)
    (hR : 0 < R) (hS : R < S)
    (inner : 2*R*bR = R^2*mass-moment)
    (outer : 2*S*bS = S^2*mass-moment) :
    (2*S*bS-2*R*bR)/(S^2-R^2) = mass := by
  have hd : 0 < S^2-R^2 := by nlinarith [mul_pos (sub_pos.mpr hS) (show 0 < S+R by linarith)]
  apply (div_eq_iff (ne_of_gt hd)).mpr
  rw [inner, outer]
  ring


open StrictQuantum StrictQuantum.Native

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
/-- The native writer has zero amplitude on either orientation of a lone
double load at every written head. -/
theorem double_outside_writer : ∀ m : Fin 10, ∀ i : Input,
    J (TrailAperture.vac m + 2) i = 0 ∧ J (TrailAperture.vac m + 6) i = 0 := by decide

/-- Each live seam is oriented away from the written set. Distance is zero
precisely on that set, and changes by at most one hop across a seam. -/
theorem cold_edge (r s m n : ℕ) (order : r ≤ s) (step : s ≤ r+1)
    (near : r = 0 ↔ 0 < m) (far : s = 0 ↔ 0 < n) :
    ((r : ℤ) - s) * ((m : ℤ) - n) = if r=0 ∧ s=1 then -(m : ℤ) else 0 := by
  by_cases eq : r = s
  · subst s
    simp
    omega
  · by_cases zero : r=0
    · have one : s=1 := by omega
      have hn : n=0 := by omega
      simp [zero, one, hn]
    · have hm : m=0 := by omega
      have hn : n=0 := by omega
      simp [hm, hn, zero]

/-- For the uniform cold preparation the outward distance acceleration is the
negative of a sum of positive incidences at the boundary of the written set. -/
theorem cold_inward {ι : Type*} (boundary : Finset ι) (m : ι → ℝ) (gain : ℝ)
    (nonempty : boundary.Nonempty) (marked : ∀ e ∈ boundary, 0 < m e)
    (stiffness : 0 < gain) : -gain * ∑ e ∈ boundary, m e < 0 := by
  exact mul_neg_of_neg_of_pos (neg_lt_zero.mpr stiffness) (Finset.sum_pos marked nonempty)

/-- A measured energy below the free band forces occupation of the written
support at every appointment where the conserved energy is E. -/
theorem cold_binding (kinetic saved E depth probability deficit : ℝ)
    (balance : E=kinetic-saved) (positive : 0 ≤ kinetic)
    (localization : saved ≤ depth*probability) (measured : E ≤ -deficit)
    (well : 0 < depth) : deficit/depth ≤ probability := by
  apply (div_le_iff₀ well).mpr
  linarith [mul_comm depth probability]


/-- A native single neighbor and an incident written contact have the same
coefficient in the insertion frequency of a double load. The exchange
diagonal loses half a unit for either shared exchange. -/
theorem insertion_clock (κ degree marks adjacent : ℝ) :
    (degree-marks)/2 - 4*κ*marks - (4*κ+1/2)*adjacent =
    degree/2 - (4*κ+1/2)*(marks+adjacent) := by ring

def bracket {A : Type*} [Ring A] (a b : A) := a*b-b*a

/-- Local change V preserves the total load inside the read surface; when
the cut's exchange current is disjoint from V the full acceleration
operator at release is unchanged, for every prepared medium. -/
theorem exterior_release {A : Type*} [Ring A] (h v p : A)
    (total : v*p=p*v) (cut : v*bracket h p=bracket h p*v) :
    bracket (h+v) (bracket (h+v) p) = bracket h (bracket h p) := by
  have first : bracket (h+v) p = bracket h p := by
    simp only [bracket, add_mul, mul_add, total]
    abel
  rw [first]
  change (h+v)*bracket h p-bracket h p*(h+v)=h*bracket h p-bracket h p*h
  rw [add_mul,mul_add,cut]
  abel

/-- To change a diagonal thermal moment a walk must both reach and return
from the altered seam. Left and right incident moments agree through d. -/
theorem return_moments {A : Type*} [Monoid A] (h g p : A) (d : ℕ)
    (left : ∀ i ≤ d, p*h^i=p*g^i) (right : ∀ i ≤ d, h^i*p=g^i*p)
    (n : ℕ) (hn : n ≤ 2*d) : p*h^n*p=p*g^n*p := by
  obtain ⟨a, b, ha, hb, hab⟩ : ∃ a b : ℕ, a ≤ d ∧ b ≤ d ∧ n=a+b := by
    exact ⟨min n d,n-min n d,by omega,by omega,by omega⟩
  rw [hab, pow_add, pow_add]
  calc
    _ = (p*h^a)*(h^b*p) := by simp only [mul_assoc]
    _ = (p*g^a)*(g^b*p) := by rw [left a ha, right b hb]
    _ = _ := by simp only [mul_assoc]

/-- The scalar-clock cut prediction on each exterior edge has a linear
ordinary-distribution term and a quadratic exclusion correction. -/
theorem clock_cut_split (u v su sv : ℝ) (hu : u≠1) (hv : v≠1) :
    (1-u-v)*(su/(1-u)-sv/(1-v)) =
    su-sv-v*su/(1-u)+u*sv/(1-v) := by
  field_simp [sub_ne_zero.mpr (Ne.symm hu), sub_ne_zero.mpr (Ne.symm hv)]
  ring


/-- Summing spectator and exchange alternatives after a sharp source read
on a live edge. The marginal ordinary probabilities at its ends are u,v;
su,sv count their live neighbors. The result is the physical edge current
derivative times the number of allowed probe addresses. -/
theorem mediated_edge (κ mu mv u v su sv : ℝ) :
    4*κ*((mu-mv)*(1-u-v)+(su-v)-(sv-u))+(u-v)/2 =
    4*κ*(mu-mv)*(1-u-v)+(4*κ*su+(4*κ+1/2)*u)-(4*κ*sv+(4*κ+1/2)*v) := by ring

/-- On an exterior edge, a coherent ordinary preparation contributes the
independently measurable derivative of its imaginary bond coherence.
At stationarity that derivative vanishes. -/
theorem coherent_exterior_edge (κ u v su sv fringeDerivative : ℝ) :
    4*κ*((su-v)-(sv-u))-fringeDerivative =
    4*κ*(su+u)-4*κ*(sv+v)-fringeDerivative := by ring


/-- At fixed injection, deleting an edge subtracts its potential drop times
its two-ended incidence from the grounded Laplacian. The measured change
in potential therefore has a source at one end and a sink at the other. -/
theorem passive_edge_source {V : Type*} [AddCommGroup V] [Module ℝ V]
    (L : V →ₗ[ℝ] V) (old changed incidence : V) (drop : ℝ)
    (balance : L changed - drop • incidence = L old) :
    L (changed-old) = drop • incidence := by
  rw [map_sub]
  rw [← balance]
  abel

/-- With a fixed grounded graph and the same two-ended incidence, any two
supplies induce proportional changes. Opposite drops reverse every reading. -/
theorem passive_supply_reversal {V : Type*} [AddCommGroup V] [Module ℝ V]
    (L : V →ₗ[ℝ] V) (unique : Function.Injective L)
    (x y incidence : V) (a b : ℝ)
    (hx : L x = a • incidence) (hy : L y = b • incidence) :
    b • x = a • y := by
  apply unique
  simp only [map_smul, hx, hy, smul_smul, mul_comm b a]

/-- Summing measured closed-neighborhood probabilities counts each ordinary
port once plus once for every live neighbor. -/
theorem source_sum {ι : Type*} (s : Finset ι) (p : ι → ℝ) (a : ι → ι → ℝ) :
    ∑ i ∈ s, (p i + ∑ j ∈ s, a i j * p j) =
    ∑ j ∈ s, (1 + ∑ i ∈ s, a i j) * p j := by
  simp only [Finset.sum_add_distrib, add_mul, one_mul, Finset.sum_mul]
  congr 1
  exact Finset.sum_comm

/-- Any normalized one-unit source on a degree-six chart has total native
stationary field at most 224 at the measured-model coupling κ = 8. -/
theorem single_source_ceiling {ι : Type*} (s : Finset ι)
    (p : ι → ℝ) (a : ι → ι → ℝ)
    (positive : ∀ j ∈ s, 0 ≤ p j) (normalized : ∑ j ∈ s, p j = 1)
    (degree : ∀ j ∈ s, ∑ i ∈ s, a i j ≤ 6) :
    32 * (∑ i ∈ s, (p i + ∑ j ∈ s, a i j * p j)) ≤ 224 := by
  rw [source_sum]
  have h : ∑ j ∈ s, (1 + ∑ i ∈ s, a i j) * p j ≤ ∑ j ∈ s, 7 * p j := by
    apply Finset.sum_le_sum
    intro j hj
    exact mul_le_mul_of_nonneg_right (by linarith [degree j hj]) (positive j hj)
  rw [← Finset.mul_sum, normalized] at h
  linarith

/-- A weighted combination of orbital work readings gives a source budget
bound for every mixture whose single-offset work is between -L and U. -/
theorem signed_work_ceiling {ι : Type*} (s : Finset ι) (k p n : ι → ℝ)
    (L U A : ℝ) (hL : 0 ≤ L) (hU : 0 ≤ U)
    (hk : ∀ i ∈ s, -L ≤ k i ∧ k i ≤ U)
    (hp : ∀ i ∈ s, 0 ≤ p i) (hn : ∀ i ∈ s, 0 ≤ n i)
    (bp : ∑ i ∈ s, p i ≤ A) (bn : ∑ i ∈ s, n i ≤ A) :
    ∑ i ∈ s, (p i - n i) * k i ≤ A * (U + L) := by
  calc
    _ ≤ ∑ i ∈ s, (U * p i + L * n i) := by
      apply Finset.sum_le_sum
      intro i hi
      nlinarith [mul_nonneg (hp i hi) (sub_nonneg.mpr (hk i hi).2),
                 mul_nonneg (hn i hi) (show 0 ≤ k i + L by linarith [(hk i hi).1])]
    _ = U * (∑ i ∈ s, p i) + L * (∑ i ∈ s, n i) := by
      rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.mul_sum]
    _ ≤ A * (U + L) := by
      nlinarith [mul_le_mul_of_nonneg_left bp hU, mul_le_mul_of_nonneg_left bn hL]

/-- Two controllable source preparations averaging to an invariant comparison
have opposite responses for every linear observer read. -/
theorem complementary_readings {V : Type*} [AddCommGroup V] [Module ℝ V]
    (read : V →ₗ[ℝ] ℝ) (a b uniform : V)
    (partition : a + b = uniform + uniform) (invariant : read uniform = 0) :
    read a = -read b := by
  have h := congrArg read partition
  simp only [map_add, invariant, add_zero] at h
  linarith

/-- Ground-state change of variables on one row of a symmetric kinetic operator.
The two neighbor sums use the same blank hopping coefficients. -/
theorem ground_poisson_row {ι : Type*} (s : Finset ι) (a f g : ι → ℝ)
    (x y d e₀ e source : ℝ)
    (blank : ∑ j ∈ s, a j * f j = (d-e₀)*x)
    (held : ∑ j ∈ s, a j * g j = (d-e)*y + source) :
    ∑ j ∈ s, a j * (f j*y-x*g j) = (e-e₀)*x*y-x*source := by
  calc
    _ = (∑ j ∈ s, a j * f j)*y - x*(∑ j ∈ s, a j * g j) := by
      rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl (fun j _ => by ring)
    _ = _ := by rw [blank, held]; ring

/-- Derivative of the normalized occupation at a flat blank mode:
the subtracted source average enforces conservation of one ordinary unit. -/
theorem weak_occupation_poisson (κ m mean n lapTangent lapOccupation : ℝ)
    (eigenDerivative : lapTangent = 2*κ*(m-mean))
    (bornDerivative : n*lapOccupation = 2*lapTangent) :
    n*lapOccupation=4*κ*(m-mean) := by linarith

/-- The measured probability of a real stationary mode has a Laplacian
equal to twice the amplitude term less a sum of squared differences. -/
theorem occupation_laplacian {ι : Type*} (s : Finset ι) (a f : ι → ℝ) (x : ℝ) :
    ∑ j ∈ s, a j*(x^2-(f j)^2) =
    2*x*(∑ j ∈ s, a j*(x-f j)) - ∑ j ∈ s, a j*(x-f j)^2 := by
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl (fun j _ => by ring)

/-- The degree commutator in the release field is dominated by the
binding gap outside one clean collar. -/
theorem screened_release_sign (prob near lapProb nearLap comm gap spread : ℝ)
    (hp : 0 ≤ prob) (hn : 0 ≤ near) (hg : 0 ≤ gap) (hd : spread ≤ 4*gap)
    (siteBound : lapProb ≤ -4*gap*prob) (neighborBound : nearLap ≤ -4*gap*near)
    (degreeBound : comm ≤ spread*near) :
    lapProb+nearLap+comm ≤ 0 := by
  nlinarith [mul_nonneg hg hp, mul_nonneg (sub_nonneg.mpr hd) hn]

/-- A bound mode loses at least this amplitude factor per hop on any
degree-six chart when its spectral energy is at most -11/2. -/
theorem bound_hop (value previous degree gap : ℝ)
    (hv : 0 ≤ value) (hd : degree ≤ 6) (hg : 11/2 ≤ gap)
    (order : value ≤ previous) (row : (degree+2*gap)*value ≤ degree*previous) :
    value ≤ (6/17)*previous := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hd) (sub_nonneg.mpr order),
             mul_nonneg (sub_nonneg.mpr hg) hv]

/-- Exterior radial work required by two orbital records sets a positive
mass increase even if interior source locations are completely free. -/
theorem exterior_enclosed_floor (inner outer hi lo : ℝ) (a : inner ≤ hi)
    (b : lo ≤ outer) : lo-hi ≤ outer-inner := by linarith

/-- A measured contrast on the active preparation transfers to its precisely
    enrolled fraction of copies at the same sign and scaled resolution. -/
theorem diluted_release_band (q active blank lo hi : ℝ) (hq : 0 ≤ q)
    (lower : lo ≤ active - blank) (upper : active - blank ≤ hi) :
    q * lo ≤ ((1-q)*blank + q*active) - blank ∧
    ((1-q)*blank + q*active) - blank ≤ q * hi := by
  have h : ((1-q)*blank + q*active) - blank = q*(active-blank) := by ring
  rw [h]
  exact ⟨mul_le_mul_of_nonneg_left lower hq, mul_le_mul_of_nonneg_left upper hq⟩

/-- Accumulating a bound on each interval of a unitary comparison bounds
    its final error, even when the signed contact time changes from tick to tick. -/
theorem timed_comparison_bound (err cost : ℕ → ℝ) (start : err 0 = 0)
    (step : ∀ n, err (n+1) ≤ err n + cost n) (n : ℕ) :
    err n ≤ ∑ i ∈ Finset.range n, cost i := by
  induction n with
  | zero => simp [start]
  | succ n ih =>
    rw [Finset.sum_range_succ]
    exact (step n).trans (by linarith [ih])

/-- When both one-record marginals are matched, any receiver change is the
    probability change in the joint positive cell times its mixed response. -/
theorem binary_joint_receiver (d₀₀ d₀₁ d₁₀ d₁₁ y₀₀ y₀₁ y₁₀ y₁₁ : ℝ)
    (left : d₀₀+d₀₁=0) (right : d₀₀+d₁₀=0) (last : d₁₀+d₁₁=0) :
    d₀₀*y₀₀ + d₀₁*y₀₁ + d₁₀*y₁₀ + d₁₁*y₁₁ =
      d₁₁*(y₀₀-y₀₁-y₁₀+y₁₁) := by
  have h₁ : d₁₀ = -d₁₁ := by linarith
  have h₂ : d₀₀ = d₁₁ := by linarith
  have h₃ : d₀₁ = -d₁₁ := by linarith
  rw [h₁, h₂, h₃]
  ring

/-- A correlated two-path ordinary phase with complex visibility x+i y has
    this squared-distance response after matching both single-head marginals. -/
theorem pair_fringe_mixed_response (x y : ℝ) :
    1-(1+x)/2-(1+x)/2+(1+x^2+y^2)/2=((1-x)^2+y^2)/2 := by ring

/-- Eliminate the measured joint-probability difference between a calibrated
    connected logarithm and a material receiver on the same preparations. -/
theorem pair_modular_transfer (q odds delta response mixed : ℝ)
    (faithfulContrast : odds ≠ 0) (reference : q = -delta * odds)
    (receiver : response = delta * mixed) : response = -(q/odds)*mixed := by
  rw [receiver, reference]
  field_simp

/-- Among accepted records sharing a witness, independent re-pairing changes
    the joint positive cell by the signed determinant divided by the stratum mass. -/
theorem witness_repair_difference (a b c d : ℝ) (accepted : a+b+c+d ≠ 0) :
    (c+d)*(b+d)/(a+b+c+d) - d = (b*c-a*d)/(a+b+c+d) := by
  field_simp
  ring

/-- A path through successive observational bands exists locally exactly
    on the intersection with the previous attainable band plus the shell interval. -/
theorem extend_shell_band (l u a b x : ℝ) (prior : l ≤ u) (shell : a ≤ b)
    (lower : l+a ≤ x) (upper : x ≤ u+b) :
    ∃ y, l ≤ y ∧ y ≤ u ∧ a ≤ x-y ∧ x-y ≤ b := by
  refine ⟨max l (x-b), le_max_left _ _, max_le prior (by linarith), ?_, ?_⟩
  · have h : max l (x-b) ≤ x-a := max_le (by linarith) (by linarith)
    linarith
  · have h := le_max_right l (x-b)
    linarith

/-- At fixed positive connected increment the signs of a mixed detector
    oppose one another when the two measured reference log odds have opposite signs. -/
theorem matched_charge_opposite_fringes {q m odds₁ odds₂ : ℝ}
    (charge : 0 < q) (mixed : 0 < m) (one : 0 < odds₁) (two : odds₂ < 0) :
    -(q/odds₁)*m < 0 ∧ 0 < -(q/odds₂)*m :=
  ⟨mul_neg_of_neg_of_pos (neg_neg_of_pos (div_pos charge one)) mixed,
   mul_pos (neg_pos.mpr (div_neg_of_pos_of_neg charge two)) mixed⟩

private theorem positive_walk_mono {ι : Type*} [Fintype ι] [DecidableEq ι]
    (K L : Matrix ι ι ℝ) (pos : ∀ i j, 0 ≤ K i j)
    (order : ∀ i j, K i j ≤ L i j) (n : ℕ) :
    ∀ i j, (K^n) i j ≤ (L^n) i j := by
  induction n with
  | zero => simp
  | succ n ih =>
    intro i j
    simp only [pow_succ, Matrix.mul_apply]
    exact Finset.sum_le_sum fun k _ => mul_le_mul (ih i k) (order k j) (pos k j)
      (le_trans (Matrix.pow_apply_nonneg pos n i k) (ih i k))

/-- The mixed response to two positive insertions in a shifted hopping
    matrix is entrywise nonnegative at every finite polynomial order. -/
theorem positive_walk_mixed {ι : Type*} [Fintype ι] [DecidableEq ι]
    (K A B : Matrix ι ι ℝ) (hk : ∀ i j, 0 ≤ K i j)
    (ha : ∀ i j, 0 ≤ A i j) (hb : ∀ i j, 0 ≤ B i j) (n : ℕ) :
    ∀ i j, 0 ≤ ((K+A+B)^n) i j - ((K+A)^n) i j - ((K+B)^n) i j + (K^n) i j := by
  have ka : ∀ i j, 0 ≤ (K+A) i j := fun i j => add_nonneg (hk i j) (ha i j)
  have kb : ∀ i j, 0 ≤ (K+B) i j := fun i j => add_nonneg (hk i j) (hb i j)
  have a : ∀ i j, (K+A) i j ≤ (K+A+B) i j := fun i j => by
    simp only [Matrix.add_apply]; linarith [hb i j]
  have b : ∀ i j, (K+B) i j ≤ (K+A+B) i j := fun i j => by
    simp only [Matrix.add_apply]; linarith [ha i j]
  induction n with
  | zero => simp
  | succ n ih =>
    intro i j
    simp only [pow_succ, Matrix.mul_apply, Matrix.add_apply, ← Finset.sum_sub_distrib,
      ← Finset.sum_add_distrib]
    apply Finset.sum_nonneg
    intro k _
    nlinarith [mul_nonneg (ih i k) (hk k j),
      mul_nonneg (sub_nonneg.mpr (positive_walk_mono _ _ ka a n i k)) (ha k j),
      mul_nonneg (sub_nonneg.mpr (positive_walk_mono _ _ kb b n i k)) (hb k j)]

/-- Separate equilibration at fixed source can reverse the sign. This is the
    normalized mixed response for additive unnormalized weights at one aperture. -/
theorem normalized_thermal_mix (z x y w u v : ℝ) (hz : 0 < z) (hx : 0 < x) (hy : 0 < y) :
    (w+u+v)/(z+x+y) - (w+u)/(z+x) - (w+v)/(z+y) + w/z =
      x*y/((z+x)*(z+x+y))*(w/z-u/x) + x*y/((z+y)*(z+x+y))*(w/z-v/y) := by
  have a : z+x ≠ 0 := by positivity
  have b : z+y ≠ 0 := by positivity
  have c : z+x+y ≠ 0 := by positivity
  field_simp
  ring

/-- The two binary head marginals are bounded by the supply of loaded patches
    on both arms. This bounds the per-patch energy of any re-pairing. -/
theorem matched_supply_mass_bound (old new loaded p e : ℝ) (ho : 0 ≤ old) (hn : 0 ≤ new)
    (a : old ≤ loaded) (b : new ≤ loaded) (supply : loaded ≤ 4*p) (he : 0 ≤ e) :
    |(new-old)*e| ≤ 4*p*e := by
  rw [abs_mul, abs_of_nonneg he]
  apply mul_le_mul_of_nonneg_right _ he
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- Independent re-pairing of binary readings within a witness class changes
    the doubly positive mass by at most a quarter of the entire stratum. -/
theorem witness_covariance_cap (a b c d n : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hc : 0 ≤ c) (hd : 0 ≤ d) (h : n = a+b+c+d) (hn : 0 < n) :
    |(b*c-a*d)/n| ≤ n/4 := by
  rw [abs_div, abs_of_pos hn, div_le_iff₀ hn]
  apply abs_le.mpr
  have ad : (a+d)^2 ≤ n^2 := by
    apply sq_le_sq' <;> linarith
  have bc : (b+c)^2 ≤ n^2 := by
    apply sq_le_sq' <;> linarith
  constructor <;> nlinarith [sq_nonneg (a-d), sq_nonneg (b-c),
    mul_nonneg ha hd, mul_nonneg hb hc]

/-- Sorting sharp ordinary records has the same acceptance rate in every
head class when each retained cell has the same weight across classes. -/
theorem sorted_head_mass {H X : Type*} [Fintype X] (prior : H → ℝ)
    (law keep : H → X → ℝ) (common : X → ℝ)
    (cells : ∀ h x, law h x * keep h x = common x) (h : H) :
    ∑ x, prior h * (law h x * keep h x) = prior h * ∑ x, common x := by
  simp only [cells, Finset.mul_sum]

/-- Two source preparations with the same normalized ordinary-position
law give the same expectation for every instrument depending only on
those positions. The instrument may itself average over private coins. -/
theorem common_position_readout {H X : Type*} [Fintype H] [Fintype X]
    (p q : H → ℝ) (w read : X → ℝ) (mass : ∑ h, p h = ∑ h, q h) :
    (∑ h, p h * ∑ x, w x * read x) = ∑ h, q h * ∑ x, w x * read x := by
  rw [← Finset.sum_mul, ← Finset.sum_mul, mass]

/-- A single pair of coefficients for two independently counted boundary
classes forces the three measured response rows to have zero determinant. -/
theorem two_cut_stationarity {g k a b c d e f s t u : ℝ}
    (one : s + g*a + k*b = 0) (two : t + g*c + k*d = 0)
    (three : u + g*e + k*f = 0) :
    a*(d*u-t*f) + b*(t*e-c*u) + s*(c*f-d*e) = 0 := by
  linear_combination (c*f-d*e)*one + (b*e-a*f)*two + (a*d-b*c)*three

/-- With fixed individual success frequencies, a change D in joint success
changes the earliest of two conditional receipts and a bypass by the smaller
of the two available time savings, even if one saving is zero. -/
theorem earliest_receipt_mixed (D clock x y : ℝ) :
    D*(clock-(clock-x)-(clock-y)+(clock-max x y)) = D*min x y := by
  rcases le_total x y with h | h
  · rw [max_eq_right h, min_eq_left h]
    ring
  · rw [max_eq_left h, min_eq_right h]
    ring

/-- Even with a separate area gain for each ball, a proposed additional
local record entropy with weight g must explain the held and ordinary
variations at that same weight. -/
theorem additional_record_balance {ao ah so sh fo fh g k : ℝ}
    (ordinary : so + g*fo + k*ao = 0) (held : sh + g*fh + k*ah = 0) :
    ao*sh-ah*so + g*(ao*fh-ah*fo) = 0 := by
  linear_combination ao*held - ah*ordinary

/-- In any common material preparation, a diagonal frozen source can enter
the initial occupation acceleration only linearly. Thus a change of joint
source law at fixed one-contact marginals has zero instantaneous contrast. -/
theorem diagonal_source_acceleration {R : Type*} [Ring R] (K V M : R)
    (diagonal : V*M = M*V) :
    ((K+V)*M - M*(K+V))*(K+V) - (K+V)*((K+V)*M - M*(K+V))
      = (K*M-M*K)*K - K*(K*M-M*K) + (K*M-M*K)*V - V*(K*M-M*K) := by
  have first : (K+V)*M-M*(K+V) = K*M-M*K := by
    rw [add_mul, mul_add, diagonal]
    abel
  rw [first]
  noncomm_ring

/-- Averaging the history-dependent part of an acceleration instead of
independently permuting the frozen source among flights leaves exactly
the covariance of the source reading with the material response. -/
theorem retained_source_covariance {ι : Type*} [Fintype ι]
    (w x y : ι → ℝ) (u v : ℝ) (norm : ∑ i, w i = 1)
    (meanX : ∑ i, w i * x i = u) (meanY : ∑ i, w i * y i = v) :
    ∑ i, w i*(x i-u)*(y i-v) = (∑ i, w i*x i*y i)-u*v := by
  calc
    _ = ∑ i, (w i*x i*y i - w i*x i*v - u*(w i*y i) + w i*(u*v)) := by
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ = _ := by
      rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
        ← Finset.sum_mul, ← Finset.mul_sum, ← Finset.sum_mul, meanX, meanY, norm]
      ring

/-- On each nondegenerate frequency pair the triangular average of the
boundary acceleration integrates to the terminal-minus-launch occupation.
The coincident pair vanishes on both sides by the zero gap. -/
theorem triangular_frequency_transfer (x : ℝ) (hx : x ≠ 0) :
    -(x^2/2) * (Real.sin (x/2)/(x/2))^2 = Real.cos x - 1 := by
  have h := Real.cos_two_mul (x/2)
  have h' := Real.sin_sq_add_cos_sq (x/2)
  rw [show 2*(x/2) = x by ring] at h
  field_simp
  nlinarith

/-- Measured edge drops generated by a scalar potential increment have
zero oriented sum on a triangle, after subtracting the native arm. -/
theorem scalar_increment_circulation (a b c p q r : ℝ) :
    ((a-b)-(p-q)) + ((b-c)-(q-r)) - ((a-c)-(p-r)) = 0 := by ring

/-- On the 192-marker fixture the crossing count is 52 for an empty clamp
and 44 for an occupied clamp. The same clamp marginal fixes the entire
area law under any redistribution of source records. A resolved positive
entropy tangent then excludes stationarity for every finite area gain. -/
theorem clamp_only_area {ι : Type*} [Fintype ι] (w c : ι → ℝ)
    (mass : ∑ i, w i = 0) (fixedClamp : ∑ i, w i*c i = 0)
    {s lo : ℝ} (measured : lo ≤ s) (positive : 0 < lo) (gain : ℝ) :
    0 < s + gain * ∑ i, w i*(52-8*c i) := by
  have expand : ∀ i, w i*(52-8*c i) = w i*52-8*(w i*c i) := by
    intro i
    ring
  simp_rw [expand]
  rw [Finset.sum_sub_distrib, ← Finset.sum_mul, ← Finset.mul_sum, mass, fixedClamp]
  simpa using (lt_of_lt_of_le positive measured)

/-- A resolved balanced-port contrast gives a strictly positive contribution
to the symmetric relative entropy of the two classical output laws. -/
theorem resolved_fringe_entropy {m z : ℝ} (visible : 0 < z) (support : z < m) :
    0 < z * Real.log ((m+z)/(m-z)) := by
  apply mul_pos visible
  apply Real.log_pos
  exact (one_lt_div (by linarith : 0 < m-z)).2 (by linarith)

/-- Diluting each local read by the other port's fraction of their summed
ordinary work gives exactly the same mean heat in both interventions. -/
theorem matched_pinching_work (x y : ℝ) :
    (y/(x+y))*x = (x/(x+y))*y := by ring

/-- Three entropy readings bound the native tangent even when the regional
observer receives no source labels. Concavity holds for a material mixture
conditioned on its classical ball and all three arms have the same support. -/
theorem measured_entropy_chords {f : ℝ → ℝ} {s lo hi : ℝ}
    (concave : ConcaveOn ℝ (Set.Icc (-1) 1) f) (tangent : HasDerivAt f s 0)
    (forward : lo ≤ f 1 - f 0) (backward : f 0 - f (-1) ≤ hi) :
    lo ≤ s ∧ s ≤ hi := by
  have zero : (0 : ℝ) ∈ Set.Icc (-1) 1 := ⟨by norm_num, by norm_num⟩
  have plus : (1 : ℝ) ∈ Set.Icc (-1) 1 := ⟨by norm_num, by norm_num⟩
  have minus : (-1 : ℝ) ∈ Set.Icc (-1) 1 := ⟨by norm_num, by norm_num⟩
  have lower := concave.slope_le_of_hasDerivAt zero plus (by norm_num) tangent
  have upper := concave.le_slope_of_hasDerivAt minus zero (by norm_num) tangent
  norm_num [slope_def_field] at lower upper
  constructor <;> linarith

/-- A common scalar with two independently calibrated ordinary gains transfers
a two-location covariance even if the field fluctuates arbitrarily. -/
theorem clock_covariance_transfer {v e c vv ve : ℝ} (hv : v ≠ 0)
    (same : vv = v*v*c) (mixed : ve = v*e*c) : ve = e/v*vv := by
  rw [same, mixed]
  field_simp

/-- If the same additive line shift is allowed for both types of clock,
equality of the two observed shifts eliminates the ordinary scalar component
whenever the ordinary gains differ. -/
theorem equal_clock_lines {v e phi common : ℝ} (gain : v ≠ e)
    (observed : v*phi+common = e*phi+common) : phi = 0 := by
  have h : (v-e)*phi = 0 := by linear_combination observed
  exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr gain)

/-- The resolvent change when adding one written interaction needs only a
left inverse on the blank chart and a right inverse on the written chart. -/
theorem inverse_mark_drop {R : Type*} [Ring R] (M A X Y : R)
    (blank : X*M = 1) (written : (M+A)*Y = 1) : X-Y = X*A*Y := by
  calc
    X-Y = X*((M+A)*Y) - (X*M)*Y := by rw [blank, written, one_mul, mul_one]
    _ = X*A*Y := by noncomm_ring

/-- The mixed change contains two ordered routes visiting both source
contacts. This identity is exact without a weak-contact expansion. -/
theorem mixed_resolvent_routes {R : Type*} [Ring R] (o a b c A B : R)
    (oa : o-a = o*A*a) (ob : o-b = o*B*b)
    (bc : b-c = b*A*c) (ac : a-c = a*B*c) :
    o-a-b+c = o*B*b*A*a + b*A*a*B*c := by
  calc
    o-a-b+c = (o-a)-(b-c) := by abel
    _ = o*A*a-b*A*c := by rw [oa, bc]
    _ = (o-b)*A*a + b*A*(a-c) := by noncomm_ring
    _ = o*B*b*A*a + b*A*a*B*c := by rw [ob, ac]; noncomm_ring

/-- Entrywise positivity of the below-band legs makes an ordered route
strictly visible between any source and receiver addresses. -/
theorem positive_clock_routes {ι : Type*} (A B : Finset ι) (hA : A.Nonempty)
    (hB : B.Nonempty) (l m r : ι → ι → ℝ) (x y : ι)
    (hl : ∀ i ∈ A, 0 < l x i) (hm : ∀ i ∈ A, ∀ j ∈ B, 0 < m i j)
    (hr : ∀ j ∈ B, 0 < r j y) :
    0 < ∑ i ∈ A, ∑ j ∈ B, l x i*m i j*r j y := by
  apply Finset.sum_pos _ hA
  intro i hi
  apply Finset.sum_pos _ hB
  intro j hj
  exact mul_pos (mul_pos (hl i hi) (hm i hi j hj)) (hr j hj)

/-- The weak-rail work is strictly positive at finite exchange
when the measured leading work dominates the bounded level correction. -/
theorem finite_rail_work {w response error : ℝ} (visible : error < response)
    (level : |w-response| ≤ error) : 0 < w := by
  linarith [(abs_le.mp level).1]

/-- Independent one-head fields reweight rows and columns and cancel from
the measured cross-product odds. -/
theorem local_fields_cancel (a b c d x y z w : ℝ)
    (hx : x ≠ 0) (hy : y ≠ 0) (hz : z ≠ 0) (hw : w ≠ 0) :
    (x*z*a)*(y*w*d)/((x*w*b)*(y*z*c)) = a*d/(b*c) := by
  by_cases hb : b = 0
  · simp [hb]
  by_cases hc : c = 0
  · simp [hc]
  field_simp

/-- After summing the held head lines the odds are a weighted mean
of their sharp odds, with weights set by the two single-contact traces. -/
theorem pooled_pair_odds {ι : Type*} [Fintype ι] (u r : ι → ℝ)
    {z : ℝ} (norm : ∑ i, u i = z) (hz : 0 < z) (hu : ∀ i, 0 ≤ u i)
    {floor : ℝ} (hr : ∀ i, floor ≤ r i) :
    floor ≤ (∑ i, u i*r i)/z := by
  apply (le_div_iff₀ hz).2
  rw [← norm, Finset.mul_sum]
  apply Finset.sum_le_sum
  intro i _
  simpa [mul_comm] using mul_le_mul_of_nonneg_left (hr i) (hu i)

/-- If every pair with an inactive line has unit odds, only the
doubly active portion of the normalized weight can reduce the odds. -/
theorem inactive_line_floor {ι : Type*} [Fintype ι] (u r mask : ι → ℝ)
    (norm : ∑ i, u i = 1) (pos : ∀ i, 0 ≤ u i)
    (inactive : ∀ i, 1-mask i ≤ r i) :
    1-(∑ i, u i*mask i) ≤ ∑ i, u i*r i := by
  rw [← norm, ← Finset.sum_sub_distrib]
  apply Finset.sum_le_sum
  intro i _
  simpa [mul_sub] using mul_le_mul_of_nonneg_left (inactive i) (pos i)

/-- The probability of seeing two active lines under the product of
the two separately tilted head laws is the product of active fractions. -/
theorem active_line_product {ι κ : Type*} [Fintype ι] [Fintype κ]
    (u : ι → ℝ) (v : κ → ℝ) :
    (∑ i, ∑ j, u i*v j) = (∑ i, u i)*(∑ j, v j) := by
  simp_rw [← Finset.mul_sum]
  rw [Finset.sum_mul]

/-- Per-bath comparisons survive arbitrary integer inventories when
both materials interrogate the same pair of sharp records. -/
theorem shared_line_power (a b u v : ℝ) (n m k : ℕ)
    (ha : 0 ≤ a) (hu : 0 ≤ u) (hv : 0 ≤ v)
    (au : u ≤ a^k) (bv : v ≤ b^k) :
    u^n*v^m ≤ (a^n*b^m)^k := by
  rw [mul_pow, ← pow_mul a n k, Nat.mul_comm n k, pow_mul a k n,
    ← pow_mul b m k, Nat.mul_comm m k, pow_mul b k m]
  exact mul_le_mul (pow_le_pow_left₀ hu au n) (pow_le_pow_left₀ hv bv m)
    (pow_nonneg hv m) (pow_nonneg (pow_nonneg ha k) n)

/-- At fixed one-head marginals, either equal-outcome cell reads the
same signed displacement. On the 2200 support the middle nonunit load
is the zero-zero cell. -/
theorem closed_middle_source {a b c d : ℝ} (mass : a+b+c+d = 0)
    (row : c+d = 0) (col : b+d = 0) : a = d := by linarith

/-- The idle port and the occupied middle port share a blank native bond.
The energy is diagonal in expectation on sharp preparations. The receiver
parks an occupied target before the echoed swap. -/
theorem middle_energy_and_echo (z o t : ℝ) (mass : z+o+t = 0) :
    (17/2*o + 65/2*t = 31/4*(z+t) + 65/4*(o+2*t)) ∧
    (o+t = ((o+2*t)-(z+t))/2) := by constructor <;> linarith

/-- Comparing the complete spectator fuel histograms sees at least twice
the joint displacement in total variation across the three bins. -/
theorem spectator_fuel_variation (d : ℝ) :
    (|d|+|-(2*d)|+|d|)/2 = 2*|d| := by
  rw [abs_neg, abs_mul, abs_of_nonneg (by norm_num : (0:ℝ) ≤ 2)]
  ring

/-- Proportional sharp-word weights commute with each closed Hamiltonian
when its nonzero matrix entries connect only words of equal weight. -/
theorem closed_sharp_stationarity (p q e : ℝ) (support : e ≠ 0 → p = q) :
    p*e=e*q := by
  by_cases h : e = 0
  · simp [h]
  · rw [support h, mul_comm]

/-- A possible signed port density must satisfy any weighted combination
of the radius-two records. The error of an approximate coefficient identity
is bounded by its residual against the same density. -/
theorem held_ball_residual {ι κ : Type*} [Fintype ι] [Fintype κ]
    (w : ι → ℝ) (k : ι → κ → ℝ) (d held : κ → ℝ) :
    (∑ x, held x*d x) - ∑ i, w i * ∑ x, k i x*d x =
      ∑ x, (held x - ∑ i, w i*k i x)*d x := by
  simp_rw [Finset.mul_sum, ← mul_assoc, sub_mul, Finset.sum_mul]
  rw [Finset.sum_sub_distrib, Finset.sum_comm (s := Finset.univ) (t := Finset.univ)
    (f := fun i x => w i*k i x*d x)]

/-- Grounded contact update, checked after application of the closed Laplacian
at an arbitrary port. Here b is its incidence, a the old contact drop,
r the unit-current contact resistance, and y the prescribed forcing. -/
theorem passage_contact_update (a b r y : ℝ) (hr : 0 ≤ r) :
    (y - a / (1+r) * b) + b * (a - a / (1+r) * r) = y := by
  have h : 1+r ≠ 0 := by linarith
  field_simp
  nlinarith

/-- Discharging the local blank-frequency premise predicts the change in
a material round trip between every pair with calibrated voltage drop d. -/
theorem passage_mixture (p q t d k : ℝ) :
    (q*t+(1-q)*(t+k*d^2)) - (p*t+(1-p)*(t+k*d^2)) = -(q-p)*k*d^2 := by ring

/-- The exterior source of a squared clock voltage has a definite sign.
The premise is current conservation at a port outside the changed contacts. -/
theorem harmonic_clock_square {ι : Type*} [Fintype ι] (w v : ι → ℝ) (x : ℝ)
    (harmonic : ∑ i, w i*(x-v i) = 0) :
    ∑ i, w i*(x^2-(v i)^2) = -(∑ i, w i*(x-v i)^2) := by
  have term (i : ι) : w i*(x^2-(v i)^2) =
      (2*x)*(w i*(x-v i)) - w i*(x-v i)^2 := by ring
  simp_rw [term, Finset.sum_sub_distrib, ← Finset.mul_sum, harmonic, mul_zero, zero_sub]

/-- Probability recursion for the finite even return ring, starting at zero
half-length. The native crossing probability at one appointment is p. -/
noncomputable def returnTransit (p : ℝ) : ℕ → ℝ
  | 0 => 0
  | n+1 => (returnTransit p n + p)/(1+returnTransit p n)

theorem returnTransit_nonneg (p : ℝ) (hp : 0 ≤ p) (n : ℕ) :
    0 ≤ returnTransit p n := by
  induction n with
  | zero => simp [returnTransit]
  | succ n ih =>
    simp only [returnTransit]
    positivity

/-- The measured finite-rate departure from a square-root response,
in a form without division or a limiting-rate assumption. -/
theorem returnTransit_error (s : ℝ) (n : ℕ) :
    (s-returnTransit (s^2) n)*(1+s)^n =
    (s+returnTransit (s^2) n)*(1-s)^n := by
  induction n with
  | zero => simp [returnTransit]
  | succ n ih =>
    have hn := returnTransit_nonneg (s^2) (sq_nonneg s) n
    have hd : 1+returnTransit (s^2) n ≠ 0 := by linarith
    simp only [returnTransit, pow_succ]
    field_simp [hd]
    linear_combination (1-s^2) * ih

/-- A conserved centered regional reading has variance equal to the
negative cross-covariance with its complement, on each sampled word. -/
theorem complementary_fluctuation (inside outside : ℝ) (conserved : inside+outside=0) :
    inside^2 = -(inside*outside) := by
  linear_combination inside * conserved

/-- If a returned cut sees either all centered home readings or none,
any response vanishing at zero integrates to zero on those cuts. -/
theorem silent_nested_cuts {ι : Type*} [Fintype ι] (weight reading : ι → ℝ)
    (response : ℝ → ℝ) (null : response 0=0) (silent : ∀ i, reading i=0) :
    ∑ i, weight i*response (reading i)=0 := by simp [silent, null]

/-- Subtract the empty contact and the two isolated contrast prices read by
spectroscopy. The remaining signed local work is proportional to test load. -/
theorem blank_contrast_field (k n t : ℝ) :
    k*(n-t)^2-k*n^2-k*t^2=t*(-2*k*n) := by ring

/-- On each eligible pair a nonzero middle deviation has an opposite partner.
This identifies the squared signed field difference through two local counts. -/
theorem paired_field_square (x y : ℝ) (partner : y ≠ 0 → x = -y) :
    (x-y)^2 = x^2+3*y^2 := by
  by_cases h : y=0
  · simp [h]
  · rw [partner h]
    ring

/-- Sharp returned records with the measured middle frequency and mean one
have a strictly positive minimum compression-work increment. -/
theorem compression_work_floor (v work : ℝ) (subset : 5/28 ≤ v)
    (calibration : work = (65/4:ℝ)*(1/14)+(3/4)*(v-11/28)) :
    1 ≤ work := by linarith

/-- In inverse scale factor s, conservation gives s e' = 3p/rho_b.
The product derivative is the active source in units of baryon density. -/
theorem active_density_derivative (e : ℝ → ℝ) (s de pressure : ℝ)
    (slope : HasDerivAt e de s) (conservation : s*de=pressure) :
    HasDerivAt (fun u => u*e u) (e s+pressure) s := by
  convert (hasDerivAt_id s).mul slope using 1
  simp [conservation]

/-- The integrated active reading of any consecutive finite set of bins is
fixed by the first and last density. The work in a bin may have either sign. -/
theorem active_density_bins (s e work : ℕ → ℝ)
    (balance : ∀ i, work i=s (i+1)*e (i+1)-s i*e i) (n : ℕ) :
    ∑ i ∈ Finset.range n, work i = s n*e n-s 0*e 0 := by
  induction n with
  | zero => simp
  | succ k ih => rw [Finset.sum_range_succ, ih, balance k]; ring

theorem active_density_cap (slo shi elo ehi E A : ℝ)
    (lo : 0 ≤ slo) (hi : 0 ≤ shi) (left : -E ≤ elo) (right : ehi ≤ E)
    (balance : (shi-slo)*A=shi*ehi-slo*elo) :
    (shi-slo)*A ≤ (shi+slo)*E := by
  rw [balance]
  nlinarith [mul_nonneg hi (sub_nonneg.mpr right), mul_nonneg lo (sub_nonneg.mpr left)]

/-- The clock waits for its first flag on each retained sharp source.
At mean load one, its two nonunit outcomes have equal weight. -/
theorem retained_first_flag_delay (f g T0 T1 T2 : ℝ) :
    ((1-g)*T1+g*(T0+T2)/2)-((1-f)*T1+f*(T0+T2)/2)
    =(g-f)*((T0+T2)/2-T1) := by ring

/-- Weighting a reversible source by its mean number of slots per event
gives detailed balance when the source advances only on a successful flag. -/
theorem tick_service_balance (a b p q k l : ℝ) (hp : p ≠ 0) (hq : q ≠ 0)
    (event_balance : a*k=b*l) : (a/p)*(p*k)=(b/q)*(q*l) := by
  field_simp [hp,hq]
  exact event_balance

/-- A finite distribution of local clock rates. Its scatter and its
first moment fix the second moment with no small-scatter expansion. -/
theorem clock_scatter_identity {ι : Type*} [Fintype ι] (w r : ι → ℝ) (m : ℝ)
    (norm : ∑ i, w i=1) (mean : ∑ i, w i*r i=m) :
    ∑ i, w i*(r i-m)^2 = (∑ i, w i*(r i)^2)-m^2 := by
  calc
    _ = ∑ i, (w i*(r i)^2-(w i*r i)*(2*m)+w i*m^2) := by
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ = _ := by
      rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.sum_mul,
        ← Finset.sum_mul, norm, mean]
      ring

/-- A finite survey corrected for the longer residence of slow clocks.
The second moment in slot time measures the arithmetic rate per renewal. -/
theorem renewal_second_moment {ι : Type*} [Fintype ι] (w r : ι → ℝ)
    (R : ℝ) (live : ∀ i, r i ≠ 0) :
    (∑ i, (w i/r i*R)*(r i)^2) = (∑ i, w i*r i)*R := by
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro i _
  field_simp [live i]

/-- Two independent uniform deadlines miss a port precisely when both are
below its measured hop distance. The normalized opening rate is its ball weight. -/
theorem two_deadline_weight (R d : ℝ) (hR : 0 < R) (hd : 0 ≤ d) (inside : d ≤ R) :
    R / 2 * (1 - (d / R)^2) = weight R d := by
  unfold weight
  have hprod := mul_nonneg (sub_nonneg.mpr inside) (show 0 ≤ R+d by linarith)
  have hnum : 0 ≤ R^2-d^2 := by nlinarith
  rw [max_eq_right (div_nonneg hnum (by linarith))]
  field_simp

/-- Growth of the finite boost weight at general nested radii. Multiplying
by a nonnegative density and summing gives the lower bound for the next ball. -/
theorem nested_weight_growth (R S d : ℝ) (hR : 0 < R) (hS : R ≤ S) :
    S * weight R d ≤ R * weight S d := by
  have hp : 0 < S := lt_of_lt_of_le hR hS
  unfold weight
  by_cases hw : 0 ≤ (R^2-d^2)/(2*R)
  · rw [max_eq_right hw]
    apply le_trans _ (mul_le_mul_of_nonneg_left (le_max_right _ _) hR.le)
    have hprod := mul_nonneg (sub_nonneg.mpr hS) (show 0 ≤ S+R by linarith)
    have hsquare : 0 ≤ S^2-R^2 := by nlinarith
    field_simp
    nlinarith [mul_nonneg (sq_nonneg d) hsquare]
  · rw [max_eq_left (le_of_not_ge hw), mul_zero]
    exact mul_nonneg hR.le (le_max_left _ _)

/-- At fixed mean one a nonzero middle load gives one binary spectral
parity. The increase in nonunit probability reduces its conditional entropy. -/
theorem spectral_parity_tangent (f g gate bit : ℝ) :
    gate*((1-g+g/2)*bit-(1-f+f/2)*bit) = -gate*(g-f)*bit/2 := by ring

/-- With four load units, m doubles leave 4-2m singles. Every long-run
ordered-pair disagreement is counted by the three population readings. -/
theorem four_load_disagreement (V m : ℝ) :
    V^2-((V-4+m)^2+(4-2*m)^2+m^2) =
      8*V-32+(24-2*V)*m-6*m^2 := by ring

/-- Increasing the record of a binary half exchange predicts a deficit of
full-pulse word returns on a boundary of blank contacts at the same activity. -/
theorem blank_history_return_balance (interior boundary du : ℝ) :
    boundary*(interior*du) + interior*(boundary*(-du)) = 0 := by ring

/-- For rational log coefficients, the program separately verifies exact
calibration on each prime. This bilinear identity is the necessary test
independent of the unknown area coefficient. -/
theorem two_source_history_response (ordinary source ao ac gain : ℝ)
    (o : ordinary = -gain*ao) (s : source = -gain*ac) :
    source*ao-ordinary*ac=0 := by rw [o,s]; ring

/-- All nonnegative prices for the mixed 1--2 current can only lower
the response once the two vacuum-current prices have been measured. -/
theorem calibrated_jump_upper (a b c p q mixed : ℝ) (decline : c ≤ 0)
    (price : 0 ≤ mixed) : a*p+b*q+c*mixed ≤ a*p+b*q := by
  linarith [mul_nonpos_of_nonpos_of_nonneg decline price]

/-- Incidence counts each internal working channel twice, each crossing
channel once. On a written four-way input the history has g additional bits. -/
theorem history_channel_residual (s a v inside g f : ℝ)
    (history : s=inside+g*f) (incidence : v=2*inside+a) :
    s+a/2-v/2=g*f := by rw [history, incidence]; ring

/-- Two ordinary aperture comparisons on the 240-port chart uniquely fix
the coefficient for a crossing channel and that for a station incidence. -/
theorem ordinary_channel_calibration (area volume : ℝ)
    (first : 12*area-32*volume= -10) (second : 15*area-37*volume= -11) :
    area=1/2 ∧ volume=1/2 := by constructor <;> linarith

/-- A controlled admixture of an ordinary donor matches the volume
measured by the channel survey. The donor has zero four-way residual. -/
theorem matched_channel_area (s a v so ao vo mix g f : ℝ)
    (source : s+a/2-v/2=g*f) (ordinary : so+ao/2-vo/2=0)
    (sameVolume : v=mix*vo) :
    a-mix*ao=2*(g*f-(s-mix*so)) := by
  rw [sameVolume] at source
  linear_combination 2*source-2*mix*ordinary

/-- Closed monitored release: blank layers act linearly on the vector of
expected material readings. Each written-contact current adds one incidence
vector before the next such layer. No independence of currents is assumed. -/
theorem two_release_currents {E : Type*} [AddCommGroup E] [Module ℝ E]
    (P Q : E →ₗ[ℝ] E) (x b : E) (j k : ℝ) :
    Q (P (x + j • b) + k • b) =
      Q (P x) + j • Q (P b) + k • Q b := by simp

/-- Opposite signed current contrasts at both visits cancel in every
terminal linear material census, with arbitrary common transport gains. -/
theorem balanced_release_census (g h j₁ j₂ k₁ k₂ : ℝ)
    (first : j₁+k₁=0) (second : j₂+k₂=0) :
    ((g*j₁+h*j₂)+(g*k₁+h*k₂))/2 = 0 := by
  have id : (g*j₁+h*j₂)+(g*k₁+h*k₂) = g*(j₁+k₁)+h*(j₂+k₂) := by ring
  rw [id, first, second]
  ring

/-- A shortcut between stations at the same original distance from the
observer leaves each radial reading unchanged for every nonnegative cost. -/
theorem equal_radius_shortcut (r s a b cost : ℝ) (hc : 0 ≤ cost)
    (viaA : r ≤ s+a) (viaB : r ≤ s+b) :
    min r (min (s+cost+a) (s+cost+b)) = r := by
  apply min_eq_left
  exact le_min (by linarith) (by linarith)

/-- Renew until the terminal observer heralds a decoded bit. A local
transfer and the occupied-target correction then fix the mean range. -/
theorem relay_from_material (p u correction time : ℝ)
    (positive : 0 < p) (localRead : p=2*u+correction)
    (renew : time=1+(1-p)*time) : time=1/(2*u+correction) := by
  rw [← localRead]
  apply (eq_div_iff (ne_of_gt positive)).2
  nlinarith

/-- Once all material placements are uniform at fixed inventory, each
station sees the same increment for a replacement on both preparations. -/
theorem uniform_release_excess (n a b load : ℝ) :
    ((a+load)/n-a/n)-((b+load)/n-b/n)=0 := by ring

/-- A swap changes the coincidence at a spectator by its signed current
times that spectator's read. Its incidence on the two endpoints is opposite. -/
theorem swap_pair_incidence (a b c wi wj : ℝ) :
    (b*c*wi+a*c*wj)-(a*c*wi+b*c*wj)=(a-b)*c*(wj-wi) := by ring

/-- The connected two-probe release contains one additional unordered pair
at fixed background count. -/
theorem connected_pair_count (n : ℝ) :
    (n+2)*(n+1)/2-2*((n+1)*n/2)+n*(n-1)/2=1 := by ring

/-- On an equal-radius circle a central acceleration with measured scalar
coefficient k acts on the chord between the two positions. -/
theorem circular_chord_tide {E : Type*} [AddCommGroup E] [Module ℝ E]
    (x y : E) (k : ℝ) : (-k) • y - (-k) • x = (-k) • (y-x) := by
  rw [smul_sub]

/-- The microscopic full-pulse swap and disagreement frequencies fix the
coherent swap coefficient in ordinary blank-clock units. The enumerated
native columns separately supply the three possible rates. -/
theorem coherent_swap_precursor (p D : ℝ) (pairs : D=2*p*(1-p)) :
    p + 144*D = p+288*p*(1-p) := by rw [pairs]; ring

/-- Two additional measured shares suffice even without resolving the
coincidence with a single target, when the supplied upper bound follows
from the conserved inventory. -/
theorem relay_margin_from_local_shares (a b c d e f cap : ℝ)
    (nonneg : 0 ≤ c) (upper : f ≤ cap*e)
    (margin : d+(1+cap)*e/4 < a+b/4) :
    d+(e+f)/4 < a+(b+c)/4 := by linarith

/-- Substituting the two measured blank-star sums into the chord-coordinate
acceleration leaves exactly the cross-ranges and the common endpoint terms. -/
theorem exterior_chord_formula (K D degU degV crossU crossV fluxU fluxV : ℝ)
    (departureU : 2*D^2*fluxU=degU*(1+D^2)-crossU)
    (departureV : 2*D^2*fluxV=crossV-degV*(1+D^2))
    (accel : K=2*(fluxV-fluxU)) :
    D^2*K=crossU+crossV-(degU+degV)*(1+D^2) := by
  rw [accel]
  linear_combination departureV-departureU

/-- Two exterior launches have common initial chord D and unit local blank
distances on both surveys. After subtracting these common terms, their
projected acceleration change is a sum of changes in squared cross-ranges.
Strict shortening of at least one of the finitely many neighbors suffices. -/
theorem exterior_chord_focus {ι : Type*} (neighbors : Finset ι)
    (old new : ι → ℝ) (D change : ℝ) (positive : 0 < D)
    (rates : change*D^2 = ∑ i ∈ neighbors, ((new i)^2-(old i)^2))
    (nearer : ∀ i ∈ neighbors, 0 ≤ new i ∧ new i ≤ old i)
    (witness : ∃ i ∈ neighbors, new i < old i) : change < 0 := by
  have strict : (∑ i ∈ neighbors, (new i)^2) < ∑ i ∈ neighbors, (old i)^2 := by
    apply Finset.sum_lt_sum
    · intro i hi
      obtain ⟨ho, hz⟩ := nearer i hi
      exact pow_le_pow_left₀ ho hz 2
    · obtain ⟨i, hi, lt⟩ := witness
      exact ⟨i, hi, pow_lt_pow_left₀ lt (nearer i hi).1 (by decide)⟩
  rw [Finset.sum_sub_distrib] at rates
  exact neg_of_mul_neg_left (by linarith : change*D^2 < 0) (le_of_lt (sq_pos_of_pos positive))

/-- At fixed imaged mass a nonnegative distribution over physical lengths
cannot produce more projected acceleration than that mass times the
pointwise upper response. Each source annulus can use its own bound. -/
theorem finite_scale_mixture_bound {ι : Type*} (scales : Finset ι)
    (mass response : ι → ℝ) (cap : ℝ)
    (nonneg : ∀ i ∈ scales, 0 ≤ mass i)
    (bounded : ∀ i ∈ scales, response i ≤ cap) :
    ∑ i ∈ scales, mass i * response i ≤ (∑ i ∈ scales, mass i) * cap := by
  rw [Finset.sum_mul]
  exact Finset.sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (bounded i hi) (nonneg i hi)

/-- Recording the first occupied aperture by differencing consecutive
no-peer probabilities conserves all exposure when the last complement is kept. -/
theorem first_peer_mass_budget (p : ℕ → ℝ) (n : ℕ) :
    (∑ i ∈ Finset.range n, (p i - p (i+1))) + p n = p 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ]
    linarith

/-- A common pool of material positions contributes its tested
partial-head/single-target coincidence to the native message probability.
A uniform extra ordinary single changes that probability by b/(4V). -/
theorem pooled_ordinary_step (a b joint extra V : ℝ) :
    (a+b/4+(joint+b*extra)/(4*V))-(a+b/4+joint/(4*V))=b*extra/(4*V) := by ring

/-- A locally measured repair probability determines an ordinary count
under fresh enrollment, including randomization between integer counts. -/
theorem pooled_ordinary_equivalent (p q b V : ℝ) (hb : b ≠ 0) (hv : V ≠ 0) :
    p+b*(4*V*(q-p)/b)/(4*V)=q := by field_simp; ring

/-- If one entirely repaired route beats a blank detour that every
single repair preserves, its integer blank saving lies between the
two measured paid-route totals. -/
theorem collective_route_window (blank detour joint single : ℝ)
    (wins : blank+joint<detour) (silent : detour≤blank+single) :
    And (joint < detour - blank) (detour - blank ≤ single) := by constructor <;> linarith

/-- Compose the circuit threshold with the exterior release identity.
A local time measurement below the blank detour suffices even when the
laboratory never measures the tiny change of material acceleration. -/
theorem collective_relay_focus {ι : Type*} (star : Finset ι)
    (old new : ι → ℝ) (w : ι) (D change blank toll threshold : ℝ)
    (hD : 0<D) (hw : w∈star) (faster : toll<threshold)
    (before : old w=blank+threshold)
    (after : new w=min (blank+threshold) (blank+toll))
    (monotone : ∀ i ∈ star, 0≤new i ∧ new i≤old i)
    (rates : change*D^2=∑ i ∈ star, ((new i)^2-(old i)^2)) :
    change<0 := by
  apply exterior_chord_focus star old new D change hD rates monotone
  refine ⟨w,hw,?_⟩
  rw [before,after,min_eq_right (by linarith)]
  linarith

/-- The port-volume term and the constant term are paid by an occupied
pool once. Zero occupancy contributes nothing to its ordinary equivalent. -/
theorem pooled_source_expectation {ι : Type*} (s : Finset ι)
    (prob occupied number : ι → ℝ) (a b c : ℝ) :
    ∑ i ∈ s, prob i * ((a+c)*occupied i+b*number i) =
      (a+c)*(∑ i ∈ s, prob i*occupied i)+b*(∑ i ∈ s, prob i*number i) := by
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i hi
  ring

/-- A symmetric renewed relay adds its success probability to the
conductance of its edge. Compare the two grounded inverse updates. -/
theorem renewed_contact_difference (p q R : ℝ) (hp : 0 ≤ p) (hq : 0 ≤ q) (hR : 0 ≤ R) :
    q/(1+q*R)-p/(1+p*R)=(q-p)/((1+p*R)*(1+q*R)) := by
  have p0 : 0<1+p*R := by positivity
  have q0 : 0<1+q*R := by positivity
  field_simp
  ring

/-- The range change to a fixed landmark is minus a squared voltage
drop. Material release at a harmonic exterior star gives the same change
at every landmark, including a landmark on the far side of the launch. -/
theorem unrestricted_landmark_release {ι : Type*} [Fintype ι]
    (w v : ι → ℝ) (x landmark gain : ℝ) (h : ∑ i, w i*(x-v i)=0) :
    2*(∑ i, w i*((-gain)*(v i-landmark)^2-(-gain)*(x-landmark)^2)) =
    -2*gain*(∑ i, w i*(x-v i)^2) := by
  have shifted : ∑ i, w i*((x-landmark)-(v i-landmark))=0 := by
    simpa only [sub_sub_sub_cancel_right] using h
  have eqn := harmonic_clock_square w (fun i => v i-landmark) (x-landmark) shifted
  simp only [sub_sub_sub_cancel_right] at eqn
  have term (i : ι) : w i*((-gain)*(v i-landmark)^2-(-gain)*(x-landmark)^2)
      = gain*(w i*((x-landmark)^2-(v i-landmark)^2)) := by ring
  simp_rw [term]
  rw [← Finset.mul_sum, eqn]
  ring

/-- At two absorbing beacons the half-transmitting link joins two
copies of the bypass. The starting observer has equal weight on both.
The open and partial hitting-probability changes have the same numerator. -/
theorem two_sheet_arrival_fraction (num S transfer : ℝ)
    (hS : 0 ≤ S) (hR : 0 ≤ S-2*transfer) :
    num/(2*(2+S)) =
      ((1+S-2*transfer)/(2*(2+S)))*(num/(1+S-2*transfer)) := by
  have hn : 1+S-2*transfer ≠ 0 := by linarith
  have hd : 2+S ≠ 0 := by linarith
  field_simp

/-- Nonnegative Dirichlet Green cross-voltage bounds the fraction of
the open response contributed by a partial head in an equal layer mixture. -/
theorem two_sheet_fraction_bounds (S transfer : ℝ)
    (h : 0 ≤ transfer) (hR : 0 ≤ S-2*transfer) :
    0 < (1+S-2*transfer)/(2*(2+S)) ∧
    (1+S-2*transfer)/(2*(2+S)) < 1/2 := by
  have hd : 0 < 2*(2+S) := by linarith
  constructor
  · exact div_pos (by linarith) hd
  · apply (div_lt_iff₀ hd).2
    linarith

/-- The measured partial-head shift is minus nine sevenths of the blank
shift on the canonical repair. Its two probe responses share a sign. -/
theorem canonical_probe_ratio (fraction : ℝ) (lo : 0<fraction) (hi : fraction<1/2) :
    5/14 < 1-(9/7)*fraction ∧ 1-(9/7)*fraction < 1 := by
  constructor <;> linarith

/-- Fix the probability of every selective ordinary-probe history. For
any measurable trajectory event, an increase in open histories displaces
exactly the same probability of blocked histories for either load. -/
theorem matched_ordinary_trajectory (a delta b opened closed selective : ℝ) :
    ((a+delta)*opened+b*selective+(1-(a+delta)-b)*closed) -
    (a*opened+b*selective+(1-a-b)*closed) = delta*(opened-closed) := by ring

/-- Squared residuals along a feasible affine path are convex. Applying
this identity coordinate by coordinate certifies the Gibbs chi-squared
bound throughout a closed interval from its two feasible endpoint reads. -/
theorem entropy_segment_square (r left right : ℝ) (hr : 0≤r) (hr' : r≤1) :
    ((1-r)*left+r*right)^2 ≤ (1-r)*left^2+r*right^2 := by
  have h := mul_nonneg (mul_nonneg hr (sub_nonneg.mpr hr')) (sq_nonneg (left-right))
  nlinarith

/-- The exact linear response dominates the quadratic uncertainty at
all sufficiently weak positive displacements, including arbitrarily
small ones that a finite laboratory sample may not resolve. -/
theorem weak_entropy_source (s upper slope loss floor gain : ℝ)
    (positive : 0<s) (small : s≤upper) (nonneg : 0≤loss)
    (margin : floor<slope-loss*upper) (estimate : s*slope-loss*s^2≤gain) :
    s*floor<gain := by
  have h := mul_le_mul_of_nonneg_left small nonneg
  have h' := mul_lt_mul_of_pos_left (lt_of_lt_of_le margin (by linarith : slope-loss*upper ≤ slope-loss*s)) positive
  nlinarith

/-- Given fixed peer reports, a fresh proposal matches their cell with
probability cell. The first matching word has its ordinary conditional
frequency. This uses independence of each renewed proposal. -/
theorem renewed_peer_record (word cell released : ℝ) (positive : 0 < cell)
    (renew : released = word + (1-cell)*released) : released = word / cell := by
  apply (eq_div_iff (ne_of_gt positive)).2
  linarith

/-- If the observed bearing obeys the first-step equation on a stopped
native kernel, every finite number of additional clock ticks preserves
its expectation. This also applies to a fixed ordinary chart throughout
an exterior region on which the two stopped kernels coincide. -/
theorem stopped_harmonic_bearing {ι : Type*} (step : (ι → ℝ) → ι → ℝ)
    (bearing : ι → ℝ) (ordinary : step bearing = bearing) (n : ℕ) :
    step^[n] bearing = bearing := by
  induction n with
  | zero => rfl
  | succ n ih => simp [Function.iterate_succ_apply, ordinary, ih]

/-- A thermal energy line in a sector of retained probability p has this
modular reading. The second term records the reference's sector weights. -/
theorem retained_sector_modular_line (p Z β E : ℝ) (hp : 0 < p) (hz : 0 < Z) :
    -Real.log (p * Real.exp (-β*E) / Z) = β*E + Real.log (Z/p) := by
  rw [Real.log_div (mul_ne_zero (ne_of_gt hp) (Real.exp_ne_zero _)) (ne_of_gt hz),
      Real.log_mul (ne_of_gt hp) (Real.exp_ne_zero _), Real.log_exp,
      Real.log_div (ne_of_gt hz) (ne_of_gt hp)]
  ring

/-- The four changing centered spectral shapes have weights
(-s-t,-s+t,s,s). Their convolution identity cancels s in the log trace.
The recorded variance change is -4225 t/8 in ordinary squared energy. -/
theorem centered_spectral_variance_charge (a b c d s t : ℝ)
    (convolution : a+b=c+d) :
    (-s-t)*a+(-s+t)*b+s*c+s*d = ((8:ℝ)/4225)*(a-b)*(-(4225:ℝ)/8*t) := by
  rw [show d=a+b-c by linarith]
  ring

/-- The ordinary Ramsey calibration of the broad centered spectrum
exceeds the narrow visibility whenever the two positive arms differ. -/
theorem broad_centered_trace_gain (x y broad narrow : ℝ) (neq : x ≠ y)
    (formula : broad-narrow=(x-y)^2/4) : narrow < broad := by
  have pos : 0<(x-y)^2 := sq_pos_of_ne_zero (sub_ne_zero.mpr neq)
  linarith

/-- A three-outcome home whose mean is fixed has only one independent
marginal change. The ordinary reference gives the three modular prices. -/
theorem fixed_home_variance_modular (a b c s k0 k1 k2 : ℝ)
    (mass : a+b+c=0) (mean : b+2*c=0) (squared : s=a+c) :
    a*k0+b*k1+c*k2 = s*(k0-2*k1+k2)/2 := by
  have ha : a=s/2 := by linarith
  have hb : b= -s := by linarith
  have hc : c=s/2 := by linarith
  rw [ha,hb,hc]
  ring

/-- When the next shared record can be updated from the current shared record
and the next ticket, every finite transcript follows the same rule. The
record may include all earlier reported positions and an observer's clock. -/
theorem shared_stopped_record {Carrier Record Ticket : Type*}
    (read : Carrier → Record) (step : Carrier → Ticket → Carrier)
    (recordStep : Record → Ticket → Record)
    (one : ∀ x t, read (step x t) = recordStep (read x) t)
    (tickets : List Ticket) (x : Carrier) :
    read (tickets.foldl step x) = tickets.foldl recordStep (read x) := by
  induction tickets generalizing x with
  | nil => rfl
  | cons t ts ih =>
    change read (ts.foldl step (step x t)) = ts.foldl recordStep (recordStep (read x) t)
    rw [ih, one]

/-- Reading only the events in one regional bin gives the retained-reference
conditional average of the whole likelihood on that bin. -/
theorem regional_likelihood_tower {ι : Type*} (s : Finset ι) (a b : ι → ℝ)
    (ha : 0 < ∑ i ∈ s, a i) (hb : 0 < ∑ i ∈ s, b i)
    (hpos : ∀ i ∈ s, 0 < b i) :
    Real.exp (Real.log (∑ i ∈ s, a i) - Real.log (∑ i ∈ s, b i)) =
      (∑ i ∈ s, b i * (a i / b i)) / (∑ i ∈ s, b i) := by
  rw [Real.exp_sub, Real.exp_log ha, Real.exp_log hb]
  congr 1
  apply Finset.sum_congr rfl
  intro i hi
  field_simp [(hpos i hi).ne']

/-- Two equal-volume surveys whose boundaries are unchanged can both be
stationary only when both corresponding regional entropy tangents vanish. -/
theorem incompatible_fixed_volume_records (x y a b η₁ η₂ : ℝ)
    (ha : a=0) (hb : b=0) (visible : x ≠ 0 ∨ y ≠ 0) :
    ¬ (x+η₁*a=0 ∧ y+η₂*b=0) := by
  rintro ⟨hx,hy⟩
  rw [ha, mul_zero, add_zero] at hx
  rw [hb, mul_zero, add_zero] at hy
  rcases visible with h | h
  · exact h hx
  · exact h hy

/-- A protected binary marker survey. A shift toward empty reports changes
the boundary by exactly the negative of the ordinary fill operation. -/
theorem protected_binary_survey (δ Aempty Afilled Kempty Kfilled η : ℝ)
    (ordinary : (Kfilled-Kempty)+η*(Afilled-Aempty)=0) :
    (δ*(Kempty-Kfilled)+η*(δ*(Aempty-Afilled))=0) ∧
    (δ*(Aempty-Afilled)=(-δ)*(Afilled-Aempty)) := by
  constructor
  · calc
      _ = -δ*((Kfilled-Kempty)+η*(Afilled-Aempty)) := by ring
      _ = 0 := by rw [ordinary]; ring
  · ring

/-- Two calibrated regional modular readings recover the ordinary rest
inventory. Here x and y are changes in the empty and double frequencies.
The inversion is observable when the calibrated product b*c is nonzero. -/
theorem regional_charge_rest_inventory (x y a b c qSix qHome rest : ℝ)
    (six : qSix = c*x) (home : qHome = a*x+b*y)
    (cal : rest = 48*(-x+3*y)) :
    b*c*rest = 48*(3*c*qHome-(b+3*a)*qSix) := by
  rw [six, home, cal]
  ring

/-- The independent clean bond spectrum and protected echo give the
same rest inventory in ordinary clock units. -/
theorem clean_bond_rest_inventory (x y energy echo : ℝ)
    (bond : energy = -(17/2)*x + 24*y) (timing : echo = -2*x) :
    6*energy-(3/2)*echo = 48*(-x+3*y) := by
  rw [bond, timing]
  ring

/-- Any nonnegative mixture of blank-to-occupied ordinary deliveries
has the opposite echo sign to a source with increased empty frequency. -/
theorem ordinary_deposit_echo_sign (x singles doubles : ℝ)
    (source : 0 < x) (hs : 0 ≤ singles) (hd : 0 ≤ doubles) :
    -2*x < 2*(singles+doubles) := by linarith

/-- Any probability law on the returned inventories has mean rest at
most eight single units. The uniform native reference has mean 37/7. -/
theorem returned_inventory_gain_bound {ι : Type*} (s : Finset ι)
    (p rest : ι → ℝ) (hp : ∀ i ∈ s, 0 ≤ p i)
    (norm : ∑ i ∈ s, p i = 1) (cap : ∀ i ∈ s, rest i ≤ 8) :
    ((∑ i ∈ s, p i * rest i) - 37/7) / (37/7) ≤ 19/37 := by
  have h : ∑ i ∈ s, p i * rest i ≤ ∑ i ∈ s, p i * 8 :=
    Finset.sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (cap i hi) (hp i hi)
  rw [← Finset.sum_mul, norm] at h
  linarith

/-- The clock potential of a journey cancels on a closed tour. -/
theorem passage_clock_closed_tour (a b c p q r gain : ℝ) :
    (gain*a+p-q)+(gain*b+q-r)+(gain*c+r-p) = gain*(a+b+c) := by ring

/-- On a block containing two measured binary occupancies the mixed
response is determined by ordinary calibration on the four inputs. -/
theorem occupied_block_expectation (r00 r10 r01 r11 x y joint : ℝ) :
    (1-x-y+joint)*r00+(x-joint)*r10+(y-joint)*r01+joint*r11 =
      r00+(r10-r00)*x+(r01-r00)*y+(r11-r10-r01+r00)*joint := by ring

/-- Pairing two independent copies subtracts the product of the marginals.
The remainder of a two-port modular reading then predicts the remainder
of any ordinarily calibrated two-port resistance. -/
theorem paired_modular_exterior_clock (joint x y ell mixed charge delay : ℝ)
    (q : charge = -ell*(joint-x*y)) (r : delay = mixed*(joint-x*y)) :
    ell*delay = -mixed*charge := by rw [q,r]; ring

/-- If a change of ruler at each accessible pair of endpoints is its
ordinary mixed coefficient times the measured coincidence change, the
same factor passes through both launches' neighbor differences. -/
theorem paired_launch_tidal_response {ι : Type*} (s : Finset ι)
    (w k d : ι → ℝ) (k₀ c : ℝ) (ruler : ∀ i ∈ s, d i = k i * c) :
    (∑ i ∈ s, w i * (d i - k₀*c)) = (∑ i ∈ s, w i * (k i-k₀)) * c := by
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro i hi
  rw [ruler i hi]
  ring

/-- At an unchanged launch the Green difference is harmonic both in the
column of the launch and in the column of the far terminal. The sum of
the ruler differences then equals the sum of incident resistance changes. -/
theorem resistance_launch_incident_sum {ι : Type*} (s : Finset ι)
    (diag home far : ι → ℝ) (h₀ f₀ : ℝ)
    (hh : ∑ i ∈ s, (home i - h₀) = 0) (hf : ∑ i ∈ s, (far i - f₀) = 0) :
    (∑ i ∈ s, ((diag i - h₀) - 2 * (far i - f₀))) =
    ∑ i ∈ s, (diag i + h₀ - 2 * home i) := by
  have e (i : ι) : diag i + h₀ - 2 * home i =
      (diag i - h₀) - 2 * (home i - h₀) := by ring
  have reduce (v : ι → ℝ) (z : ℝ) :
      (∑ i ∈ s, ((diag i - h₀) - 2 * (v i - z))) =
        (∑ i ∈ s, (diag i - h₀)) - 2 * ∑ i ∈ s, (v i - z) := by
    rw [Finset.sum_sub_distrib, Finset.mul_sum]
  simp_rw [e]
  rw [reduce far f₀, reduce home h₀, hf, hh]

/-- A harmonic beacon coordinate has zero initial displacement under any
finite generator whose position acceleration is proportional to this
ordinary nearest-neighbor difference. -/
theorem harmonic_launch_acceleration {ι : Type*} (s : Finset ι)
    (p : ι → ℝ) (f : ι → ℝ) (f0 gain accel : ℝ)
    (weights : ∑ i ∈ s, p i = 1) (harmonic : ∑ i ∈ s, p i*f i = f0)
    (release : accel = gain*(∑ i ∈ s, p i*(f i-f0))) : accel = 0 := by
  simp_rw [mul_sub] at release
  rw [Finset.sum_sub_distrib, ← Finset.sum_mul, weights, harmonic] at release
  simpa using release

/-- On a route visiting both homes, re-pairing at fixed marginals moves
equal weight to the first and last receipts. For every time readout
the connected response is its second difference at those three times. -/
theorem front_time_readout (x y joint f₀ f₂ f₄ : ℝ) :
    ((1-x-y+joint)*f₀+(x+y-2*joint)*f₂+joint*f₄) -
    ((1-x-y+x*y)*f₀+(x+y-2*x*y)*f₂+x*y*f₄) =
      (joint-x*y)*(f₀-2*f₂+f₄) := by ring

/-- Any centering time gives the same change in second moment. The first
moment cancels, so the expression is also the change in time variance. -/
theorem front_delay_variance (t μ c fraction : ℝ) :
    fraction*c*((t-μ)^2-2*(t+2-μ)^2+(t+4-μ)^2) = 8*fraction*c := by ring

/-- Even with arbitrary dependence between registers, two normalized
occupancy totals in [0,1] have joint moment at most their smaller mean.
Their centered covariance is then at most one quarter. -/
theorem capped_joint_covariance (x y joint : ℝ) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (jx : joint ≤ x) (jy : joint ≤ y) : joint - x*y ≤ 1/4 := by
  rcases le_total x y with h | h
  · nlinarith [sq_nonneg (x-1/2), mul_nonneg hx (sub_nonneg.mpr h)]
  · nlinarith [sq_nonneg (y-1/2), mul_nonneg hy (sub_nonneg.mpr h)]

/-- Re-pairing two delay totals with the same marginals replaces their
cross moment by the product of their means. With two extra ticks per
occupied home the mean squared disagreement changes by minus eight
times the occupancy covariance. -/
theorem front_clock_agreement (xx yy xy x y : ℝ) :
    4*(xx+yy-2*xy)-4*(xx+yy-2*x*y) = -8*(xy-x*y) := by ring

/-- An additive two-member modular price is blind to every signed
joint intervention whose row and column marginals vanish. -/
theorem ordered_register_entropy {ι κ : Type*} (s : Finset ι) (t : Finset κ)
    (d : ι → κ → ℝ) (a : ι → ℝ) (b : κ → ℝ)
    (rows : ∀ i ∈ s, ∑ j ∈ t, d i j = 0)
    (cols : ∀ j ∈ t, ∑ i ∈ s, d i j = 0) :
    ∑ i ∈ s, ∑ j ∈ t, d i j * (a i + b j) = 0 := by
  simp_rw [mul_add, Finset.sum_add_distrib]
  have h₁ : ∑ i ∈ s, ∑ j ∈ t, d i j * a i = 0 := by
    apply Finset.sum_eq_zero
    intro i hi
    rw [← Finset.sum_mul, rows i hi, zero_mul]
  rw [h₁, zero_add, Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro j hj
  rw [← Finset.sum_mul, cols j hj, zero_mul]

/-- A common observation keeps likelihood bounds when it pools finitely
many components with nonnegative reference weights. Division is only
needed for outcomes of positive reference probability. -/
theorem pool_likelihood_bounds {ι : Type*} (s : Finset ι) (w r : ι → ℝ)
    (lo hi : ℝ) (hw : ∀ i ∈ s, 0 ≤ w i)
    (hr : ∀ i ∈ s, lo ≤ r i ∧ r i ≤ hi) :
    (∑ i ∈ s, w i) * lo ≤ ∑ i ∈ s, w i*r i ∧
    (∑ i ∈ s, w i*r i) ≤ (∑ i ∈ s, w i) * hi := by
  simp_rw [Finset.sum_mul]
  exact ⟨Finset.sum_le_sum fun i hi' => mul_le_mul_of_nonneg_left (hr i hi').1 (hw i hi'),
    Finset.sum_le_sum fun i hi' => mul_le_mul_of_nonneg_left (hr i hi').2 (hw i hi')⟩

/-- For fixed one-member occupations, a same-home coincidence fixes the
change of the accepted departure frequency, in units of the blank rate. -/
theorem departure_from_same_home_coincidence (d₀ d₁ d₂ : ℝ)
    (prob : d₀ + d₁ + d₂ = 0) (mean : d₁ + 2*d₂ = 0) :
    d₀ + (2/3)*d₁ + (1/2)*d₂ = d₂/6 := by linarith

/-- The constant term in the return resolvent fixes the change in the
variance of a return to a uniformly weighted station. -/
theorem renewal_variance_difference (n h₁ h₂ v₁ v₂ : ℝ)
    (one : v₁ = 2*n^2*h₁+n-n^2) (two : v₂ = 2*n^2*h₂+n-n^2) :
    v₁-v₂ = 2*n^2*(h₁-h₂) := by rw [one,two]; ring

/-- Outside a conductance intervention the changed Green columns are
harmonic. On the commute ruler the remaining Laplacians of the diagonal
are measured by the two local return-variance contrasts. -/
theorem exterior_release_from_return_variance (n a b v w acceleration : ℝ)
    (first : v = 12*n^2*a) (second : w = 12*n^2*b)
    (release : acceleration = -(a+b)/2) :
    24*n^2*acceleration = -(v+w) := by rw [first,second,release]; ring

/-- A change in the same-station ratio requires more than a common
multiplicative change in both clocks' frequencies. -/
theorem common_rate_cross_product (f₀ g₀ f₁ g₁ scale : ℝ)
    (f : f₁ = scale*f₀) (g : g₁ = scale*g₀) : f₁*g₀ = g₁*f₀ := by
  rw [f,g]
  ring

/-- Change in ordinary relative release on the commute ruler from two
renewed neighbor-clock readings, at clean exterior launch sites. -/
theorem neighbor_period_release (n d e t u a b acceleration : ℝ)
    (first : d*t = 6*n*a) (second : e*u = 6*n*b)
    (release : acceleration = (a+b)/2) :
    12*n*acceleration = d*t+e*u := by rw [first, second, release]; ring

/-- Subtract the two Foster traces. Every increase of edge resistance
is paid for by the lost conductance on the changed links. -/
theorem foster_contrast_budget {ι : Type*} (s : Finset ι) (r blank c : ι → ℝ)
    (trace_eq : ∑ i ∈ s, c i*r i = ∑ i ∈ s, blank i) :
    ∑ i ∈ s, (r i - blank i) = ∑ i ∈ s, (1-c i)*r i := by
  simp_rw [sub_mul, one_mul, Finset.sum_sub_distrib]
  rw [trace_eq]

/-- On the straight variation from a correlated reference to the product
of its head marginals, the lost modular score is the drop in the sum of
the two directed relative entropies. Both contract under observation. -/
theorem erasure_reading_bounds (f r f' r' q : ℝ)
    (forward : f' ≤ f) (reverse : r' ≤ r)
    (nf : 0 ≤ f') (nr : 0 ≤ r') (reading : q = (f+r)-(f'+r')) :
    0 ≤ q ∧ q ≤ f+r := by constructor <;> linarith

/-- A proper ball in a three-dimensional metric cone has V = Ωr³/3
and A = Ωr². Eliminating the first-order radial adjustment gives the
fixed-volume area change from the solid-angle change. -/
theorem cone_ball_fixed_volume_area (r omega dr dw da : ℝ) (hr : 0 < r)
    (volume : dw*r^3+3*omega*r^2*dr = 0)
    (area : da = dw*r^2+2*omega*r*dr) : 3*da = r^2*dw := by
  apply mul_left_cancel₀ (ne_of_gt hr)
  linear_combination 3*r*area + 2*volume

/-- A positive ordinary calibration and a nonnegative erased source
require a negative additional entropy response on an expanding ball. -/
theorem expanding_ball_missing_charge (eta area q extra : ℝ)
    (he : 0 < eta) (ha : 0 < area) (hq : 0 ≤ q)
    (balance : q + extra + eta*area = 0) :
    extra ≤ -(eta*area) ∧ extra < 0 := by
  constructor <;> linarith [mul_pos he ha]

/-- Radius changes at fixed cone-ball volume are proportional to the area
slope for every intervention. Calibrating on ordinary matter eliminates
any source-price term that depends differentiably on radius alone. -/
theorem moving_radius_ordinary_calibration (g h a b x y ratio : ℝ)
    (ordinary : x + (g+h)*a = 0) (matching : b = ratio*a) :
    y+(g+h)*b = y-ratio*x := by rw [matching]; linear_combination ratio*ordinary

/-- On a fixed-volume ordinary intervention, independently measured
entropy per parked energy fixes area per energy at a stationary survey. -/
theorem ordinary_energy_area_slope (entropy area energy beta eta : ℝ)
    (response : entropy = beta*energy) (balance : entropy+eta*area = 0) :
    eta*area = -(beta*energy) := by linarith

/-- A cone-ball area measured at fixed volume is quadratic in its
reference radius. Halving that radius multiplies the apparent
fourth-order coefficient by four for the same solid-angle variation. -/
theorem cone_tower_rescaling (ell slope : ℝ) (positive : 0 < ell) :
    (((ell/2)^2*slope)/((ell/2)^4)) = 4*((ell^2*slope)/ell^4) := by
  field_simp
  ring

/-- If the two sharp counterflow outcomes have identical ranging records,
the expected value of every fixed ruler on them is independent of the
probability of exchange. -/
theorem equal_ranging_exchange {R : Type*} (r s : R) (eq : r = s)
    (area : R → ℝ) (p q : ℝ) :
    (1-p)*area r+p*area s = (1-q)*area r+q*area s := by
  rw [eq]
  ring

/-- A measured finite increase in entropy on a mixture, together with a
bound on the reference objective deficit at the same survey volume,
forces a decrease of the independently calibrated area. -/
theorem finite_pairing_area_bound (s₀ s₁ a₀ a₁ eta t u g₁ g₂ gap : ℝ)
    (gain : 0 < eta)
    (concavity : t*g₁+u*g₂ ≤ s₁-s₀)
    (ceiling : s₁+eta*a₁ ≤ s₀+eta*a₀+gap) :
    (t*g₁+u*g₂-gap)/eta ≤ a₀-a₁ := by
  apply (div_le_iff₀ gain).2
  linarith

/-- Two probe readings explained by the same ordinary dose must have
zero two-by-two calibration determinant, irrespective of their units. -/
theorem common_ordinary_force_determinant (j₁ j₂ o₁ o₂ h : ℝ)
    (p₁ : j₁ = h * o₁) (p₂ : j₂ = h * o₂) : j₁ * o₂ - j₂ * o₁ = 0 := by
  rw [p₁, p₂]
  ring

/-- An observed source cannot factor through a collection of surveying
records when two realizable preparations share the records and differ in source. -/
theorem source_varies_on_equal_survey {P R : Type*} (read : P → R)
    (source : P → ℝ) (x y : P) (same : read x = read y) (gain : source x < source y) :
    ¬ ∃ price : R → ℝ, ∀ z, source z = price (read z) := by
  rintro ⟨price, h⟩
  have eq : source x = source y := (h x).trans ((congrArg price same).trans (h y).symm)
  exact (ne_of_lt gain) eq

/-- A lower entropy bound inferred by a source reader implies a finite
area bound after including the independently observed change in volume. -/
theorem finite_source_measured_volume (s₀ s₁ a₀ a₁ v₀ v₁ eta price budget gain : ℝ)
    (ordinary : 0 < eta) (source : gain ≤ s₁ - s₀)
    (ceiling : s₁ + eta*a₁ - price*v₁ ≤ s₀ + eta*a₀ - price*v₀ + budget) :
    (gain - price*(v₁-v₀) - budget)/eta ≤ a₀-a₁ := by
  apply (div_le_iff₀ ordinary).2
  nlinarith

/-- Substituting the finite entropy balance makes explicit the marginal price,
relative-entropy cost, volume change and reference deficit in a charge test. -/
theorem finite_connected_cost_contraction (s₀ s₁ a₀ a₁ v₀ v₁ eta price budget q m cost : ℝ)
    (ordinary : 0 < eta) (balance : s₁ - s₀ = q + m - cost)
    (ceiling : s₁ + eta*a₁ - price*v₁ ≤ s₀ + eta*a₀ - price*v₀ + budget) :
    (q+m-cost - price*(v₁-v₀) - budget)/eta ≤ a₀-a₁ :=
  finite_source_measured_volume _ _ _ _ _ _ _ _ _ _ ordinary balance.ge ceiling

/-- Four source readings at opposite pointer phases determine every affine
two-quadrature response relative to the dephased arm. -/
theorem source_pointer_quadratures (d u v a b : ℝ) :
    (d+a*u+b*v)-d =
      a*((d+u)-(d-u))/2 + b*((d+v)-(d-v))/2 := by ring

/-- A source fringe and its finite remainder give a two-sided area
forecast at the ordinary price, including the measured volume change. -/
theorem source_fringe_area_interval (q lo hi x v eta price eps : ℝ)
    (gain : 0 < eta) (lower : lo ≤ q) (upper : q ≤ hi)
    (left : q-eps ≤ price*v-eta*x) (right : price*v-eta*x ≤ q+eps) :
    (price*v-hi-eps)/eta ≤ x ∧ x ≤ (price*v-lo+eps)/eta := by
  constructor
  · apply (div_le_iff₀ gain).2
    linarith [mul_comm eta x]
  · apply (le_div_iff₀ gain).2
    linarith [mul_comm eta x]

/-- An ordinary calibration of a linear sum of two local readings
requires the same linear sum of their two-probe determinants to vanish. -/
theorem regional_charge_weight_determinant (m o c p a b w k : ℝ)
    (regional : m+w*c=k*a) (ordinary : o+w*p=k*b) :
    (m*b-o*a)+w*(c*b-p*a)=0 := by
  linear_combination b*regional-a*ordinary

/-- Averaging opposite probes removes an arbitrary charge-odd force. If every
ordinary preparation has nonnegative even impulse, a positive mixture must too. -/
theorem charge_even_ordinary_impulse {ι : Type*} (s : Finset ι) (w p n : ι → ℝ)
    (a b odd : ℝ) (hw : ∀ i ∈ s, 0 ≤ w i) (himp : ∀ i ∈ s, 0 ≤ p i + n i)
    (ha : a = (∑ i ∈ s, w i * p i) + odd)
    (hb : b = (∑ i ∈ s, w i * n i) - odd) : 0 ≤ (a + b)/2 := by
  have h : 0 ≤ ∑ i ∈ s, w i * (p i + n i) :=
    Finset.sum_nonneg fun i hi => mul_nonneg (hw i hi) (himp i hi)
  simp only [mul_add, Finset.sum_add_distrib] at h
  linarith

/-- A scalar reading that increases on every member increases by at least
the same amount on each normalized lottery, including empirically calibrated ones. -/
theorem service_lottery_drift {ι : Type*} (s : Finset ι) (w d : ι → ℝ) (gap : ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (unit : (∑ i ∈ s, w i) = 1)
    (drift : ∀ i ∈ s, gap ≤ d i) : gap ≤ ∑ i ∈ s, w i * d i := by
  have h := Finset.sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (drift i hi) (hw i hi)
  calc
    gap = (∑ i ∈ s, w i) * gap := by rw [unit, one_mul]
    _ = ∑ i ∈ s, w i * gap := Finset.sum_mul ..
    _ ≤ _ := h

/-- A stationary mixture of a drifting service and a bounded restoring
operation has a minimum restoring probability. -/
theorem necessary_restoring_probability (p d r gap span : ℝ)
    (p0 : 0 ≤ p) (p1 : p ≤ 1) (hd : gap ≤ d) (hr : -span ≤ r)
    (pos : 0 < gap+span) (steady : (1-p)*d+p*r=0) : gap/(gap+span) ≤ p := by
  apply (div_le_iff₀ pos).2
  have left := mul_le_mul_of_nonneg_left hd (sub_nonneg.mpr p1)
  have right := mul_le_mul_of_nonneg_left hr p0
  linarith

/-- The perpendicular to the two quadrature mismatch coefficients
equalizes the two ordinary-normalized material impulses. The acquisition
separately establishes nonzero denominators, amplitude, and sign. -/
theorem matched_pointer_impulses (r₁ r₂ i₁ i₂ o₁ o₂ : ℝ) :
    ((i₁*o₂-i₂*o₁)*r₁ - (r₁*o₂-r₂*o₁)*i₁)*o₂ =
    ((i₁*o₂-i₂*o₁)*r₂ - (r₁*o₂-r₂*o₁)*i₂)*o₁ := by ring

/-- A positive ordinary response fixes the normalization used in a
simultaneous clock and material comparison. -/
theorem normalized_pointer_response (a b r i o k : ℝ) (pos : 0 < k)
    (response : a*r+b*i=k*o) : (a/k)*r+(b/k)*i=o := by
  field_simp [ne_of_gt pos]
  nlinarith [response]

/-- On a fixed conical chart the fixed-volume area variation at radii
r, r/2 and r/4 has relative sizes 1, 1/4 and 1/16. An additive site
reading and the central flat-ball boost are affine in radius. The
three-radius contrast removes both while retaining the area signal. -/
theorem three_radius_cone_contrast (x site boost : ℝ) :
    (x + site + boost) - 3*(x/4 + site + boost/2) +
      2*(x/16 + site + boost/4) = 3*x/8 := by ring

/-- At a common reference a differentiable law of regional records has a
linear derivative. Two variations agreeing on every measured source
coordinate must therefore agree on its derivative, including for signed
redistributions of ordinary preparations. -/
theorem matched_regional_linear_read {ι κ : Type*} (hist : Finset ι)
    (sites : Finset κ) (w : ι → ℝ) (src : ι → κ → ℝ) (v gain : κ → ℝ)
    (meter : ι → ℝ)
    (matched : ∀ j ∈ sites, v j = ∑ i ∈ hist, w i * src i j)
    (ordinary : ∀ i ∈ hist, meter i = ∑ j ∈ sites, gain j * src i j) :
    (∑ j ∈ sites, gain j * v j) = ∑ i ∈ hist, w i * meter i := by
  calc
    _ = ∑ j ∈ sites, ∑ i ∈ hist, gain j * (w i * src i j) := by
      apply Finset.sum_congr rfl
      intro j hj
      rw [matched j hj, Finset.mul_sum]
    _ = ∑ i ∈ hist, ∑ j ∈ sites, gain j * (w i * src i j) := Finset.sum_comm
    _ = ∑ i ∈ hist, w i * meter i := by
      apply Finset.sum_congr rfl
      intro i hi
      rw [ordinary i hi, Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro j _
      ring

/-- On directions preserving the measured regional source, a response
mediated by k further real readings has rank at most k. -/
theorem additional_readings_rank_bound {m n k r : ℕ}
    (response : Matrix (Fin m) (Fin n) ℝ)
    (gain : Matrix (Fin m) (Fin k) ℝ) (source : Matrix (Fin k) (Fin n) ℝ)
    (factors : response = gain * source) (observed : r ≤ response.rank) : r ≤ k := by
  rw [factors] at observed
  exact observed.trans ((Matrix.rank_mul_le_left gain source).trans gain.rank_le_width)

/-- The two conjugate cross terms at a partial observation give the
native work from the measured real bond and its directed current. -/
theorem pointer_cut_work (g z : ℂ) :
    (g*z + star (g*z)).re = g.re*(2*z.re) - g.im*(2*z.im) := by
  simp [Complex.mul_re]
  ring

/-- An orientation perpendicular to the measured work row has zero
mean work and a signed modular increment fixed by its determinant. -/
theorem zero_work_modular_read (j w qR qI scale : ℝ) :
    (scale*j)*w - (scale*w)*j = 0 ∧
    (scale*j)*qR + (scale*w)*qI = scale*(j*qR+w*qI) := by
  constructor <;> ring

/-- A continuous adjustment bracketed by the volume meter admits an
energy-neutral, equal-volume point. A measured full modular increment above
the stationarity error forces its survey area to decrease. -/
theorem fixed_volume_flight_contraction (a b gain err : ℝ)
    (volume work price area : ℝ → ℝ) (hab : a ≤ b) (hg : 0 < gain)
    (hc : ContinuousOn volume (Set.Icc a b))
    (left : 0 ≤ volume a) (right : volume b ≤ 0)
    (neutral : ∀ u ∈ Set.Icc a b, work u = 0)
    (signal : ∀ u ∈ Set.Icc a b, err < price u)
    (balance : ∀ u ∈ Set.Icc a b, volume u = 0 → |price u + gain * area u| ≤ err) :
    ∃ u ∈ Set.Icc a b, volume u = 0 ∧ work u = 0 ∧ area u < 0 := by
  obtain ⟨u, hu, hv⟩ := intermediate_value_Icc' hab hc ⟨right, left⟩
  exact ⟨u, hu, hv, neutral u hu,
    area_from_source _ _ _ _ _ hg le_rfl (signal u hu) (balance u hu hv)⟩

/-- Shrinking a cap while retaining one central receiver leaves a finite
energy times the radius. A nonzero constant modular term then excludes
bounded inferred stress at all sufficiently small positive radii. -/
theorem one_receiver_stress_lower (q c r : ℝ) (hr : 0 < r)
    (small : |c*r| ≤ |q|/2) : |q|/(2*r^4) ≤ |q+c*r|/r^4 := by
  have h := abs_sub (q+c*r) (c*r)
  have he : q+c*r-c*r = q := by ring
  rw [he] at h
  have bound : |q| / 2 ≤ |q+c*r| := by linarith
  calc
    |q|/(2*r^4) = (|q|/2)/r^4 := by ring
    _ ≤ |q+c*r|/r^4 := div_le_div_of_nonneg_right bound (by positivity)

/-- With squared radii in units of 1/64, a unit shell of any radius
contributes between zero and twelve to the three-radius boost contrast. -/
theorem nested_shell_kernel (x : ℝ) :
    0 ≤ max (16-x) 0 - 5 * max (4-x) 0 + 4 * max (1-x) 0 ∧
    max (16-x) 0 - 5 * max (4-x) 0 + 4 * max (1-x) 0 ≤ 12 := by
  rcases le_total x 1 with h | h
  · rw [max_eq_left (by linarith : 0 ≤ 16-x), max_eq_left (by linarith : 0 ≤ 4-x),
      max_eq_left (by linarith : 0 ≤ 1-x)]
    constructor <;> linarith
  rcases le_total x 4 with h' | h'
  · rw [max_eq_left (by linarith : 0 ≤ 16-x), max_eq_left (by linarith : 0 ≤ 4-x),
      max_eq_right (by linarith : 1-x ≤ 0)]
    constructor <;> linarith
  rcases le_total x 16 with h'' | h''
  · rw [max_eq_left (by linarith : 0 ≤ 16-x), max_eq_right (by linarith : 4-x ≤ 0),
      max_eq_right (by linarith : 1-x ≤ 0)]
    constructor <;> linarith
  rw [max_eq_right (by linarith : 16-x ≤ 0), max_eq_right (by linarith : 4-x ≤ 0),
    max_eq_right (by linarith : 1-x ≤ 0)]
  norm_num

/-- An affine one-receiver cap read has a shell contrast independent
of all finite choices of centered energy. -/
theorem nested_cap_contrast (q c : ℝ) :
    (1/2 : ℝ)*(q+c/2) - 5*(1/4 : ℝ)*(q+c/4) + 4*(1/8 : ℝ)*(q+c/8) = -q/4 := by
  ring

/-- The regional word price and the clock fringe count the same eligible
pairs. After calibrating the conditional price, one reading predicts
the other without knowing the probability of admission. -/
theorem word_price_predicts_fringe (p κ f q : ℝ) (price : q = p*κ)
    (measurable : κ ≠ 0) : p*f = q*f/κ := by
  apply (eq_div_iff measurable).2
  rw [price]
  ring

/-- A time-sensitive balanced read with zero source contrast, together
with a nonzero contrast on the other read, excludes a common time shift. -/
theorem word_reads_exclude_time_shift (balanced other signal : ℝ)
    (sensitive : balanced ≠ 0) (observed : signal ≠ 0) :
    ¬ ∃ t : ℝ, balanced*t = 0 ∧ other*t = signal := by
  rintro ⟨t, ht, he⟩
  have hz := (mul_eq_zero.mp ht).resolve_left sensitive
  apply observed
  calc
    signal = other*t := he.symm
    _ = 0 := by rw [hz, mul_zero]

/-- A single positive mass measure on the surveyed positions must satisfy
the same two-aperture bound as every possible unit mass. The functions
are the kernels summed over all centers before the mass is integrated. -/
theorem all_center_positive_source_bound {ι : Type*} (s : Finset ι)
    (mass small large : ι → ℝ) (c : ℝ)
    (positive : ∀ i ∈ s, 0 ≤ mass i)
    (survey : ∀ i ∈ s, large i ≤ c * small i) :
    (∑ i ∈ s, mass i * large i) ≤ c * ∑ i ∈ s, mass i * small i := by
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun i hi => by
    have h := mul_le_mul_of_nonneg_left (survey i hi) (positive i hi)
    linarith only [h, show mass i * (c * small i) = c * (mass i * small i) by ring]

/-- Changing correlation at fixed one-word frequencies leaves an
additive position-acceleration read unchanged. -/
theorem correlation_annihilates_additive_read {ι κ : Type*}
    (s : Finset ι) (t : Finset κ) (δ : ι → κ → ℝ) (f : ι → ℝ) (g : κ → ℝ)
    (rows : ∀ i ∈ s, ∑ j ∈ t, δ i j = 0)
    (cols : ∀ j ∈ t, ∑ i ∈ s, δ i j = 0) :
    ∑ i ∈ s, ∑ j ∈ t, δ i j * (f i + g j) = 0 := by
  simp_rw [mul_add, Finset.sum_add_distrib]
  have left : (∑ i ∈ s, ∑ j ∈ t, δ i j * f i) = 0 := by
    apply Finset.sum_eq_zero
    intro i hi
    rw [← Finset.sum_mul, rows i hi, zero_mul]
  have right : (∑ i ∈ s, ∑ j ∈ t, δ i j * g j) = 0 := by
    rw [Finset.sum_comm]
    apply Finset.sum_eq_zero
    intro j hj
    rw [← Finset.sum_mul, cols j hj, zero_mul]
  rw [left, right, add_zero]

/-- Two changes on disjoint cells have no mixed initial acceleration
for a read that preserves those cells. The untouched generator is arbitrary. -/
theorem disjoint_word_initial_acceleration {R : Type*} [Ring R] (h a b x : R)
    (ab : a*b = 0) (ba : b*a = 0) (axb : a*x*b = 0) (bxa : b*x*a = 0) :
    let acc := fun k : R => 2*(k*x*k) - k*k*x - x*(k*k)
    acc (h+a+b) - acc (h+a) - acc (h+b) + acc h = 0 := by
  dsimp only
  calc
    _ = 2*(a*x*b + b*x*a) - (a*b+b*a)*x - x*(a*b+b*a) := by noncomm_ring
    _ = 0 := by simp [ab, ba, axb, bxa]

/-- If a fair tag is shared instead of drawn independently at two
sites, the change of the second moment of their summed clock readings
is a product of two local timing contrasts. The marginal second
moments cancel. Equal local responses therefore force a square. -/
theorem shared_tag_clock_variance (a b c d : ℝ) :
    2 * ((a*c+b*d)/2 - ((a+b)/2)*((c+d)/2)) = (a-b)*(c-d)/2 := by ring

/-- Every finite-rate deadline in the symmetric word-clock experiment
has a strictly positive timely-return contrast whenever the two
addresses have positive probability of joint admission. -/
theorem word_deadline_survival_gain (p r : ℝ) (hp : 0 < p) (hr : 0 < r)
    (hlt : r < 1) : 0 < p * (r*(1-r))^2 / 1296 := by
  have : 0 < 1-r := sub_pos.mpr hlt
  positivity

/-- A common change of duration would square the half-rate survival
ratio at the quarter-rate clock. The finite lower bound on its baseline
and the measured gains give a strict inequality on every positive aperture. -/
theorem two_word_deadlines_exclude_common_duration (u v x : ℝ) (hx : 0 < x)
    (hu : 0 < u) (hr : 9*u < 32*v) :
    1 + (9*x/16)/v < (1+x/u)^2 := by
  have hv : 0 < v := by linarith
  have h : (9*x/16)/v < 2*(x/u) := by
    apply (div_lt_iff₀ hv).2
    field_simp [ne_of_gt hu]
    nlinarith [mul_pos hx (show 0 < 32*v-9*u by linarith)]
  nlinarith [sq_nonneg (x/u)]

/-- The squared separation of two unit directions determines the depth
of their inferred midpoint on every individual recorded pair. -/
theorem pair_midpoint_depth {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (a b : E) (ha : ‖a‖ = 1) (hb : ‖b‖ = 1) :
    ‖(1/2 : ℝ) • (a+b)‖^2 = 1 - ‖a-b‖^2/4 := by
  rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/2),
    mul_pow, norm_add_sq_real, norm_sub_sq_real, ha, hb]
  ring

/-- Integrating the angularly averaged squared chord kernel gives the
fifth power of the aperture in the boundary-centered moment. -/
theorem midpoint_shell_integral (r : ℝ) :
    (∫ u in (0:ℝ)..r, (r^2-u^2)^2) = 8*r^5/15 := by
  let f := fun u : ℝ => r^4*u-2*r^2*u^3/3+u^5/5
  have h (u : ℝ) : HasDerivAt f ((r^2-u^2)^2) u := by
    dsimp [f]
    convert (((hasDerivAt_id u).const_mul (r^4)).sub
      ((((hasDerivAt_id u).pow 3).const_mul (2*r^2)).div_const 3)).add
      (((hasDerivAt_id u).pow 5).div_const 5) using 1
    simp
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun u _ => h u)
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  dsimp [f]
  ring

/-- Once the projected midpoint read and the conditional word price are
known, the boundary entropy has the same normalization for every radius
on the stated round, isotropic branch. -/
theorem midpoint_boost_matches_boundary (κ r : ℝ) :
    (15*r/32)*(κ*(2*r^3/15)) = κ*(r^2/4)^2 := by ring

/-- On a crossing with exchange element -1/2, any positive scalar mass
coupling to a resolved potential difference changes the phase-sensitive
initial acceleration. -/
theorem exterior_potential_read (gain near far : ℝ)
    (positive : 0 < gain) (different : near < far) :
    0 < (-1/2 : ℝ)*(gain*near-gain*far) := by
  nlinarith [mul_pos positive (sub_pos.mpr different)]

/-- A positive shell of midpoint sources cut off at squared radius s
contributes min(u,s) to enclosed mass, with decreasing secant slopes.
Positive sums obey the same three-ring test. -/
theorem midpoint_mass_kernel_concave (u v w s : ℝ) (hu : u ≤ v) (hv : v ≤ w) :
    (v-u)*(min w s-min v s) ≤ (w-v)*(min v s-min u s) := by
  rcases le_total s u with h | h
  · simp [min_eq_right h, min_eq_right (h.trans hu), min_eq_right ((h.trans hu).trans hv)]
  rcases le_total s v with h' | h'
  · rw [min_eq_right (h'.trans hv), min_eq_right h', min_eq_left h]
    simpa using mul_nonneg (sub_nonneg.mpr hv) (sub_nonneg.mpr h)
  rcases le_total s w with h'' | h''
  · rw [min_eq_right h'', min_eq_left h', min_eq_left h]
    nlinarith [mul_nonneg (sub_nonneg.mpr hu) (sub_nonneg.mpr h'')]
  rw [min_eq_left h'', min_eq_left h', min_eq_left h, mul_comm]

/-- A pair of material clock modes with unchanged restricted generator
has identical action at every order of the free propagator. The two
matrices can have different numbers of ports and clock levels. -/
theorem clock_subspace_powers {ι κ : Type*} [Fintype ι] [Fintype κ]
    [DecidableEq ι] [DecidableEq κ] (H : Matrix ι ι ℂ) (D : Matrix κ κ ℂ)
    (T : Matrix ι κ ℂ) (fixed : H * T = T * D) (n : ℕ) :
    H ^ n * T = T * D ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
    calc
      H ^ (n+1) * T = H ^ n * (H * T) := by rw [pow_succ, Matrix.mul_assoc]
      _ = H ^ n * (T * D) := by rw [fixed]
      _ = (T * D ^ n) * D := by rw [← Matrix.mul_assoc, ih]
      _ = T * D ^ (n+1) := by rw [Matrix.mul_assoc, pow_succ]

/-- The squared distances used to compare the three records obey this
identity for every potential source position on the declared chart. -/
theorem three_region_distances {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (x n : E) (unit : ‖n‖ = 1) :
    2*((9/16 : ℝ)-‖x-(1/2 : ℝ) • n‖^2) =
      ((11/16 : ℝ)-‖x-n‖^2) + ((15/16 : ℝ)-‖x‖^2) := by
  simp only [norm_sub_sq_real, norm_smul, Real.norm_eq_abs, real_inner_smul_right,
    abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/2), unit]
  ring

/-- Positive parts keep this bound at every point. It applies to
arbitrary positions and nonnegative weights, including off-screen mass. -/
theorem positive_three_regions {ι : Type*} (sites : Finset ι)
    (mass mid near center : ι → ℝ) (positive : ∀ i ∈ sites, 0 ≤ mass i)
    (chart : ∀ i ∈ sites, 2*mid i = near i + center i) :
    ∑ i ∈ sites, mass i * (2 * max 0 (mid i)) ≤
      ∑ i ∈ sites, mass i * (max 0 (near i) + max 0 (center i)) := by
  apply Finset.sum_le_sum
  intro i hi
  apply mul_le_mul_of_nonneg_left _ (positive i hi)
  by_cases hm : 0 ≤ mid i
  · rw [max_eq_right hm, chart i hi]
    exact add_le_add (le_max_right _ _) (le_max_right _ _)
  · rw [max_eq_left (le_of_not_ge hm), mul_zero]
    exact add_nonneg (le_max_left _ _) (le_max_left _ _)

/-- The flat rest-ball weight is one-Lipschitz in positive aperture
radius at every source location. Multiplication by 2π restores the
rest-ball normalization used in the observer comparisons. -/
theorem aperture_lipschitz (a b d : ℝ) (ha : 0 < a) (hab : a ≤ b) (hd : 0 ≤ d) :
    0 ≤ weight b d - weight a d ∧ weight b d - weight a d ≤ b-a := by
  have hb : 0 < b := lt_of_lt_of_le ha hab
  have hd2 := sq_nonneg d
  have monotone : (a^2-d^2)/(2*a) ≤ (b^2-d^2)/(2*b) := by
    apply (div_le_div_iff₀ (by positivity) (by positivity)).2
    nlinarith [mul_nonneg (sub_nonneg.mpr hab) (show 0 ≤ a*b + d^2 by positivity)]
  constructor
  · exact sub_nonneg.mpr (max_le_max_left 0 monotone)
  unfold weight
  by_cases hda : d ≤ a
  · have hdb := hda.trans hab
    have posa : 0 ≤ (a^2-d^2)/(2*a) := div_nonneg (by nlinarith [sq_nonneg (a-d), mul_self_le_mul_self hd hda]) (by positivity)
    have posb := posa.trans monotone
    rw [max_eq_right posa, max_eq_right posb]
    field_simp
    nlinarith [mul_nonneg (sub_nonneg.mpr hab)
      (show 0 ≤ a*b-d^2 by nlinarith [mul_self_le_mul_self hd hda, mul_nonneg (le_of_lt ha) (sub_nonneg.mpr hab)])]
  · have had : a ≤ d := le_of_not_ge hda
    have neg : (a^2-d^2)/(2*a) ≤ 0 := div_nonpos_of_nonpos_of_nonneg
      (by nlinarith [mul_self_le_mul_self (le_of_lt ha) had]) (by positivity)
    rw [max_eq_left neg, sub_zero]
    apply max_le (sub_nonneg.mpr hab)
    apply (div_le_iff₀ (by positivity)).2
    nlinarith [mul_self_le_mul_self (le_of_lt ha) had, sq_nonneg (b-a)]

/-- For two possible surface encodings of each address, the slower
angular delivery cannot beat the inner product of the hemisphere lifts. -/
theorem paired_delivery_bound {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (x y u v : E) :
    min (inner ℝ (x+u) (y+v)) (inner ℝ (x-u) (y-v)) ≤
      inner ℝ x y + ‖u‖ * ‖v‖ := by
  have h := real_inner_le_norm u v
  have left := min_le_left (inner ℝ (x+u) (y+v)) (inner ℝ (x-u) (y-v))
  have right := min_le_right (inner ℝ (x+u) (y+v)) (inner ℝ (x-u) (y-v))
  have sum : inner ℝ (x+u) (y+v) + inner ℝ (x-u) (y-v) =
      2*(inner ℝ x y + inner ℝ u v) := by
    simp only [inner_add_left,inner_add_right,inner_sub_left,inner_sub_right]
    ring
  linarith

/-- A common normal gives both transfers the same angular score. -/
theorem paired_delivery_attained {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (x y n : E) (p q : ℝ) (hx : inner ℝ x n = 0) (hy : inner ℝ n y = 0)
    (hn : ‖n‖ = 1) :
    min (inner ℝ (x+p • n) (y+q • n)) (inner ℝ (x-p • n) (y-q • n)) =
      inner ℝ x y + p*q := by
  simp [inner_add_left,inner_add_right,inner_sub_left,inner_sub_right,
    real_inner_smul_left,real_inner_smul_right,hx,hy,hn,mul_comm]

/-- One shell of proper energy e contributes c_i e to the spherical
mass aspect. Volume compensation gives the ordinary Einstein area factor. -/
theorem curved_shell_balance (G L e ci c r : ℝ) (hc : c ≠ 0) (hr : r ≠ 0) :
    (-2*c/r)*(4*Real.pi*G*L^2*e*(ci/c-1)) =
      -4*G*(2*Real.pi*L^2*(ci-c)*e/r) := by
  field_simp

/-- A matched pair of apertures in an empty radial gap gives the
enclosed gravitational mass even when the proper mass is unknown. -/
theorem paired_aperture_mass (c d y z m proper : ℝ) (hd : d<c)
    (first : y=m-c*proper) (second : z=m-d*proper) :
    (c*z-d*y)/(c-d)=m := by
  apply (div_eq_iff (ne_of_gt (sub_pos.mpr hd))).2
  rw [first,second]
  ring

/-- A public mixture of two deadlines attains exactly the original
mean volume, including at a finite change of writer. -/
theorem pair_volume_coin (lo hi target : ℝ) (different : lo<hi) :
    (1-(target-lo)/(hi-lo))*lo + ((target-lo)/(hi-lo))*hi = target := by
  field_simp [ne_of_gt (sub_pos.mpr different)]
  ring

/-- The measured energy and one internal transit distinguish the complete
nine-word calibration. The calibration values come from native swap clocks. -/
def localEnergyTime : Fin 9 → ℕ × ℕ :=
  ![(0,2), (32,2), (33,2), (32,3), (64,2), (65,3), (33,3), (65,2), (66,2)]

theorem localEnergyTime_decodes :
    ∀ i j : Fin 9, localEnergyTime i = localEnergyTime j → i = j := by decide

/-- The same sharp words are identified by their summed energy and the
terminal position bit after one fixed native seam pulse on ports 3 and 6. -/
def localEnergyArrival : Fin 9 → ℕ × Bool :=
  ![(0,true), (32,true), (33,false), (32,false), (64,true),
    (65,true), (33,true), (65,false), (66,true)]

theorem localEnergyArrival_decodes :
    ∀ i j : Fin 9, localEnergyArrival i = localEnergyArrival j → i = j := by decide

/-- Rewriting the alphabet through a calibrated bijection preserves each
finite sum, including the entropy or the modular expectation on a pair. -/
theorem measured_record_sum {A B : Type*} [Fintype A] [Fintype B]
    (rec : A ≃ B) (score : B → ℝ) :
    (∑ a, score (rec a)) = ∑ b, score b := rec.sum_comp score

/-- Even moving every pair to the most favorable depth cannot avoid the
minimum price growth required by an observed expanding-screen density. -/
theorem expanding_pair_price_bound (growth f₁ f₂ density : ℝ)
    (hg : 0 ≤ growth) (hf : f₂ ≤ 1)
    (observed : 16*density*f₁ = growth*f₂) :
    16*density*f₁ ≤ growth := by
  rw [observed]
  simpa using mul_le_mul_of_nonneg_left hf hg

end OPH.EinsteinBranch.SourceApertures
