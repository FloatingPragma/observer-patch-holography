import ObserverPatchHolography.EinsteinBranch.WrittenPopulation

/-! Energy inside an arbitrary set of written edges, its load-free modular
charge, and the finite range of a released ordinary probe. -/

namespace OPH.EinsteinBranch.TrailAperture
open scoped BigOperators
open NativeBoundary

def vac (i : Fin 10) : Nat := 9*(i.val+1)

def vacuumEntry (k : ℤ) (i j : Fin 10) : ℤ :=
  2*((if i=j then 1+2*k*(contrast (vac j))^2 else 0) -
    (if vac i = reflection (vac j) then 1 else 0))

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
/-- All written-memory matrix entries on empty ordinary endpoints at the
native stiffness and on the zero-stiffness control. Blank memory has zero energy. -/
theorem vacuum_block : ∀ i j : Fin 10,
    h4 8 (vac i) (vac j) = vacuumEntry 8 i j ∧
    h4 0 (vac i) (vac j) = vacuumEntry 0 i j := by decide

/-- The complete-board positive deficit bounds every subset of the recorded
edges, also when matter and memories have become entangled. -/
theorem aperture {V : Type*} [DecidableEq V] (p s : Finset V) (hs : s ⊆ p)
    (loss : V → ℝ) (hl : ∀ i ∈ p, 0 ≤ loss i) (cap : ℝ)
    (ht : ∑ i ∈ p, loss i ≤ cap) (g : ℝ) :
    g*(s.card : ℝ) - cap ≤ ∑ i ∈ s, (g-loss i) ∧
    ∑ i ∈ s, (g-loss i) ≤ g*(s.card : ℝ) := by
  have small := Finset.sum_le_sum_of_subset_of_nonneg hs (fun i hi _ => hl i hi)
  have pos := Finset.sum_nonneg (fun i hi => hl i (hs hi))
  rw [Finset.sum_sub_distrib, Finset.sum_const, nsmul_eq_mul]
  constructor <;> nlinarith

/-- Exact positive deficit interval for the 50.5-unit donor and the
original equal-degree endpoint boards. The matched blank operator is
48 plus one half of a degree-six graph Laplacian. -/
theorem total_deficit (n h loss : ℝ) (low : 48 ≤ h) (high : h ≤ 54)
    (balance : 50.5 + 8.5*n - 12.75 = h+8.5*n-loss) :
    10.25 ≤ loss ∧ loss ≤ 16.25 := by constructor <;> linarith

/-- For two independent equilibrium factors the modular energy of the
product and the sum of the individual logarithms coincide. -/
theorem product_charge (p q : ℝ) (hp : 0 < p) (hq : 0 < q) :
    -Real.log (p*q) + Real.log p + Real.log q = 0 := by
  rw [Real.log_mul hp.ne' hq.ne']
  ring

/-- Summing energy readings with a missing flag in the Gibbs-reference
measurement cancels common log-probability offsets between preparations. -/
theorem flagged (p k a b : ℝ) :
    p*(a+k)+(1-p)*b - (p*a+(1-p)*b)=p*k := by ring

/-- All outcome operators of an arbitrary remote measurement are invariant
under the isolated source unitary. The same identity applies to an entire
remote protocol's effective outcome operator. -/
theorem inaccessible {A : Type*} [Monoid A] (u adj o : A)
    (unitary : adj*u=1) (separate : o*u=u*o) :
    adj*o*u=o := by rw [mul_assoc, separate, ←mul_assoc, unitary, one_mul]

/-- If the source perturbation annihilates the freely reachable states up
to distance d, the first d evolution moments of the full generator agree. -/
theorem flight_moments {A : Type*} [Ring A] (h v x : A) (d : ℕ)
    (far : ∀ k < d, v*(h^k*x)=0) :
    ∀ k ≤ d, (h-v)^k*x=h^k*x := by
  intro k
  induction k with
  | zero => simp
  | succ k ih =>
    intro hk
    rw [pow_succ', pow_succ' h, mul_assoc, ih (by omega), ←mul_assoc,
      sub_mul, sub_mul, mul_assoc v, far k (by omega), sub_zero]

/-- Squared velocity of a compact ordinary mass plus a truncated, constant
linear density, at radius ratio s. The inner and outer enclosed dark
fractions are x and y. The result allows every positive hop length. -/
theorem two_radii (s q x y : ℝ) (hs : 1 ≤ s) (hq : 0 ≤ q)
    (tail : y ≤ 1) (growth : y ≤ s*x) :
    (s+q)*(1+q*y) ≤ s*(1+q)*(1+q*x) := by
  have a := mul_nonneg (mul_nonneg hq (sub_nonneg.mpr hs)) (sub_nonneg.mpr tail)
  have b := mul_nonneg (mul_nonneg hq (show 0 ≤ 1+q by linarith)) (sub_nonneg.mpr growth)
  nlinarith

/-- A change of source basis preserves both the spectrum of a density and
the mean of every invariant logarithm. The spectral entropy is determined
by this characteristic polynomial for a positive Hermitian density. -/
theorem source_spectrum_and_mean {n : Type*} [Fintype n] [DecidableEq n]
    (U : (Matrix n n ℂ)ˣ) (rho K : Matrix n n ℂ)
    (symmetry : K * (U : Matrix n n ℂ) = (U : Matrix n n ℂ) * K) :
    (U.val * rho * U⁻¹.val).charpoly = rho.charpoly ∧
    (U.val * rho * U⁻¹.val * K).trace = (rho*K).trace := by
  constructor
  · exact Matrix.charpoly_units_conj U rho
  · calc
      _ = (U.val * (rho * U⁻¹.val * K)).trace := by simp only [mul_assoc]
      _ = (rho * U⁻¹.val * K * U.val).trace := Matrix.trace_mul_comm _ _
      _ = (rho*K).trace := by simp [mul_assoc, symmetry]

/-- The source symmetry also preserves the filled-in reference at every
probability, provided it fixes the product of the two reference marginals. -/
theorem commutes_thinned {n : Type*} [Fintype n] (U A B : Matrix n n ℂ) (p : ℂ)
    (ha : U*A=A*U) (hb : U*B=B*U) :
    U*((1-p) • A + p • B)=((1-p) • A + p • B)*U := by
  simp [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul, ha, hb]

/-- Propagating the two rigorous material intervals to the routed Bell endpoint
uses the exact affine response on the conserved-charge blocks. -/
theorem routed_displacement (e f g : ℝ)
    (loe : 0.069104635 ≤ e) (hie : e ≤ 0.069104665)
    (lof : -0.134266663 ≤ f) (hif : f ≤ -0.134266648) (relation : g=(3*e-f)/2) :
    0.1707902765 ≤ g ∧ g ≤ 0.170790329 := by constructor <;> linarith

/-- A nonnegative post-BBN photon work history bounds the fraction of cold
stock made later than any cutoff by its mean formation scale factor. -/
theorem late_stock {V : Type*} [DecidableEq V] (late all : Finset V)
    (sub : late ⊆ all) (f a : V → ℝ) (cut : ℝ) (hf : ∀ i ∈ all, 0 ≤ f i)
    (ha : ∀ i ∈ all, 0 ≤ a i) (hc : ∀ i ∈ late, cut ≤ a i) :
    cut*(∑ i ∈ late, f i) ≤ ∑ i ∈ all, a i * f i := by
  rw [Finset.mul_sum]
  exact (Finset.sum_le_sum (fun i hi => mul_le_mul_of_nonneg_right (hc i hi)
    (hf i (sub hi)))).trans (Finset.sum_le_sum_of_subset_of_nonneg sub
    (fun i hi _ => mul_nonneg (ha i hi) (hf i hi)))

end OPH.EinsteinBranch.TrailAperture
