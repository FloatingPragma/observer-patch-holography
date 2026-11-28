import ObserverPatchHolography.EinsteinBranch.NativeBoundary
import Mathlib

/-! Written energy and regional charges at fixed incident energy. The native
integer matrix and finite cell inequalities are verified with the source
calculation. The operator and positive-stress inequalities here state their
premises explicitly. -/

namespace OPH.EinsteinBranch.WrittenPopulation
open Matrix
open scoped Matrix ComplexOrder BigOperators

/-- Splitting the written space by its write image releases only the part of
the memory-energy ruler outside that image. -/
theorem complement_split {A : Type*} [Ring A] (G P : A) (hP : P*P=P)
    (hG : P*G=G*P) : (1-P)*G*(1-P) = G-P*G*P := by
  have h : P*G*P = P*G := by rw [hG,mul_assoc,hP,←hG]
  rw [sub_mul,one_mul,mul_sub,mul_one,sub_mul,h,←hG]
  abel

/-- Bound in arbitrary complex states. Columns of X comprise an ensemble
amplitude; taking this trace includes every mixed density matrix. -/
theorem memory_ceiling {s t u : Type*} [Fintype s] [Fintype t] [Fintype u]
    (G H : Matrix s s ℂ) (W : Matrix s t ℂ) (B : Matrix t t ℂ)
    (X : Matrix s u ℂ) (hB : B.PosSemidef) (identity : G-H = W*B*Wᴴ) :
    ((Xᴴ*H*X).trace).re ≤ ((Xᴴ*G*X).trace).re := by
  have h := ((hB.mul_mul_conjTranspose_same W).conjTranspose_mul_mul_same X).trace_nonneg
  rw [← identity, Matrix.mul_sub, Matrix.sub_mul, Matrix.trace_sub] at h
  exact (Complex.le_def.mp (sub_nonneg.mp h)).1

/-- Positive-cap matching costs sum without independence of local inputs. -/
theorem matching_bound {V : Type*} (s : Finset V) (gain base : V → ℝ) (total cap : ℝ)
    (hc : 0≤cap) (each : ∀ v ∈ s, gain v ≤ cap*base v)
    (incident : ∑ v ∈ s, base v ≤ total) : ∑ v ∈ s, gain v ≤ cap*total := by
  calc
    ∑ v ∈ s, gain v ≤ ∑ v ∈ s, cap*base v := Finset.sum_le_sum each
    _ = cap*∑ v ∈ s, base v := (Finset.mul_sum ..).symm
    _ ≤ cap*total := by gcongr

theorem federation_bound (gain external other gates : ℝ)
    (stock : gain ≤ external) (coupling : external ≤ 9*gates)
    (interior : 12*gates ≤ other) : gain ≤ 3/7*(external+other) := by linarith

theorem sharp_federation_bound (gain external other mismatch gates : ℝ)
    (stock : gain ≤ 17/4*mismatch) (energy : 8*mismatch ≤ external)
    (ends : mismatch ≤ gates) (interior : 12*gates ≤ other) :
    gain ≤ 17/80*(external+other) := by linarith

/-- Channel-resolved equilibrium readouts cannot separate equal channel
populations even when every spectral outcome is retained. -/
theorem indistinguishable {V : Type*} (s : Finset V) (a b f : V → ℝ)
    (same : ∀ i ∈ s, a i=b i) : ∑ i ∈ s, a i*f i = ∑ i ∈ s, b i*f i := by
  apply Finset.sum_congr rfl
  intro i hi
  rw [same i hi]

/-- Mirror-averaged radial coefficient calibrated on an unmatched sharp word,
a matched sharp word and a bright native write. -/
theorem two_source_reads (unmatched matched bright u v : ℝ) :
    (1-u)*unmatched + (2*u-2*v)*matched + (2*v-u)*bright =
    unmatched+(2*matched-bright-unmatched)*u+2*(bright-matched)*v := by ring

/-- The four-dimensional source response allowed by the conserved left load
minus memory coin and mirror exchange decomposes into exactly the native
coincidence and echo-blank observables plus a common scalar. -/
theorem mirror_operator (a b c : ℝ) :
    (!![a,0,0,0; 0,b,c,0; 0,c,b,0; 0,0,0,a] : Matrix (Fin 4) (Fin 4) ℝ) =
    a • (1 : Matrix (Fin 4) (Fin 4) ℝ) +
    (b-c-a) • !![0,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,0] +
    (2*c) • !![0,0,0,0; 0,1/2,1/2,0; 0,1/2,1/2,0; 0,0,0,0] := by
  ext i j
  fin_cases i <;> fin_cases j <;> norm_num [Matrix.one_apply] <;> ring

/-- The two nontrivial pieces of the native intertwining relation on all
eighteen branches of the sparse write, including zero and double loads. -/
theorem native_intertwining :
    ∀ (i : StrictQuantum.Native.Input) (c : Bool),
    NativeBoundary.reflection (StrictQuantum.Native.row i c) =
      StrictQuantum.Native.row (i.2, i.1) c ∧
    NativeBoundary.contrast (StrictQuantum.Native.row i c) =
      (StrictQuantum.Native.n i.1 : ℤ) - StrictQuantum.Native.n i.2 := by decide

/-- Each symmetric regional outcome is an affine mixture in the bright
probability. The two memory-energy parities have equal compressions. -/
theorem charge_law (v b₈ d₈ b₉ d₉ : ℝ) :
    (v/2)*b₈ + ((1-v)/2)*d₈ + (v/2)*b₉ + ((1-v)/2)*d₉ =
    (1-v)*((d₈+d₉)/2) + v*((b₈+b₉)/2) := by ring

/-- A calibrated regional expectation can replace the echo read when the
bright and dark means differ. Coupling to coincidence remains independent. -/
theorem modular_and_coincidence (a b c k kd kb u v : ℝ)
    (cal : k = (1-v)*kd+v*kb) (resolved : kb ≠ kd) :
    a+b*u+c*((k-kd)/(kb-kd)) = a+b*u+c*v := by
  have h : (k-kd)/(kb-kd)=v := by
    apply (div_eq_iff (sub_ne_zero.mpr resolved)).2
    rw [cal]
    ring
  rw [h]

/-- Positive mixtures of unmatched, matched-dark and matched-bright
sources realize every allowed coincidence and unwrite pair. -/
theorem attainable (u v : ℝ) (hv : 0 ≤ v) (hvu : v ≤ u) (hu : u ≤ 1) :
    0 ≤ 1-u ∧ 0 ≤ u-v ∧ 0 ≤ v ∧ (1-u)+(u-v)+v = 1 ∧ (u-v)+v = u := by
  exact ⟨by linarith, by linarith, hv, by ring, by ring⟩

/-- Exact extent of the radial ambiguity at fixed modular distribution. -/
theorem response_diameter (a b c v : ℝ) (hi : v ≤ 1) :
    abs ((a+b+c*v)-(a+b*v+c*v)) = abs b * (1-v) := by
  rw [show (a+b+c*v)-(a+b*v+c*v)=b*(1-v) by ring, abs_mul,
    abs_of_nonneg (sub_nonneg.mpr hi)]

/-- With one mobile load, at most two bonds on a simple path can touch
it. Empty written bonds contribute at least eight per unit head at all times. -/
theorem path_storage {V : Type*} (s : Finset V) (touch cost : V → ℝ)
    (per_edge : ∀ i ∈ s, 8*(1-touch i) ≤ cost i)
    (one_load : ∑ i ∈ s, touch i ≤ 2) (blank : ℝ) (blank_cap : blank ≤ 18) :
    8*((s.card : ℝ)-2)-18 ≤ (∑ i ∈ s, cost i)-blank := by
  have h := Finset.sum_le_sum per_edge
  simp only [← Finset.mul_sum, Finset.sum_sub_distrib, Finset.sum_const,
    nsmul_eq_mul, mul_one] at h
  linarith

/-- Static conservation multiplied by radius squared, in geometric units.
Primes are evaluated at one radius. -/
theorem stress_divergence (r rho pr pt grad pr' : ℝ) (hr : 0 < r)
    (balance : pr' = -(rho+pr)*grad+2*(pt-pr)/r) :
    2*r*pr+r^2*pr'=2*r*pt-r^2*(rho+pr)*grad := by
  rw [balance]
  field_simp
  ring

/-- One shell's contribution to 2πr³p_r is at most (1+ε) times
its positive exterior mass. x=r/s, t=-p_t/ρ, p=p_r/ρ, g=s·Φ'. -/
theorem shell_support (x t p g eps : ℝ) (hx : 0 ≤ x) (hxr : x≤1)
    (ht : t≤1) (hp : -1≤p) (hpu : p≤1) (hg : g≤eps) (he : 0≤eps) :
    x*(t+(1+p)*g/2) ≤ 1+eps := by
  have h : (1+p)*g ≤ (1+p)*eps := mul_le_mul_of_nonneg_left hg (by linarith)
  have h' : (1+p)*eps ≤ 2*eps := mul_le_mul_of_nonneg_right (by linarith) he
  calc
    _ ≤ x*(1+eps) := mul_le_mul_of_nonneg_left (by linarith) hx
    _ ≤ 1*(1+eps) := mul_le_mul_of_nonneg_right hxr (by linarith)
    _ = _ := by ring

/-- Finite integrated exterior contributions include a possible outer
surface traction. Isolated support takes boundary = 0. -/
theorem tail_mass {V : Type*} (s : Finset V) (mass source : V → ℝ) (eps boundary pr : ℝ)
    (balance : pr=boundary + ∑ i ∈ s, source i)
    (shells : ∀ i ∈ s, source i ≤ (1+eps)*mass i) :
    pr ≤ boundary+(1+eps)*(∑ i ∈ s, mass i) := by
  rw [balance, Finset.mul_sum]
  exact add_le_add (le_refl _) (Finset.sum_le_sum shells)

theorem lens_total (lens mass exterior p eps eta : ℝ) (hm : 0 ≤ mass) (he : 0 ≤ eps)
    (compact : 0 ≤ eta) (weak : eta < 1)
    (pseudo : (1-eta)*lens = mass*(1-eta/2)+p) (pressure : p ≤ (1+eps)*exterior) :
    lens ≤ (1+eps)*(mass+exterior)/(1-eta) := by
  apply (le_div_iff₀ (by linarith)).2
  nlinarith [mul_nonneg hm compact, mul_nonneg hm he]

end OPH.EinsteinBranch.WrittenPopulation
