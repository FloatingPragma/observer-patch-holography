import ObserverPatchHolography.ObservablePriority.NativeColumns
import Mathlib.Tactic

/-! Energy on the original native writer, the binary owner boundary, and a
finite budget for two deposited records per write. -/

namespace OPH.EinsteinBranch.NativeBoundary
open StrictQuantum StrictQuantum.Native
open scoped BigOperators

/-- Memory label from one coin branch, entire outgoing load pair from another.
Indices are the global row numbers of the deposited sparse writer. -/
def mixed (i : Input) (c d : Bool) : Nat := (row i c / 9)*9 + row i d % 9

def reflection (s : Nat) : Nat :=
  9 * (1+5*((s/9-1)/5)+(4-(s/9-1)%5)) + s%9

def contrast (s : Nat) : ℤ := (((s/9-1)%5 : ℕ) : ℤ)-2

def inputs : List Input :=
  [(.empty,.empty),(.empty,.single),(.empty,.double),(.single,.empty),
   (.single,.single),(.single,.double),(.double,.empty),(.double,.single),(.double,.double)]

/-- Written block of twice the complement of the native write image. -/
def qentry (r s : Nat) : ℤ :=
  (if r=s then 2 else 0) - (inputs.map (fun i => J r i * J s i)).sum

/-- Written block of four times the local Hamiltonian. Its second factor
acts on a column s as the diagonal contrast and a single reflected label. -/
def h4 (k : ℤ) (r s : Nat) : ℤ :=
  qentry r s*(1+2*k*(contrast s)^2)-qentry r (reflection s)

def mixedDiag (k : ℤ) (i : Input) (c d : Bool) : ℤ :=
  h4 k (mixed i c d) (mixed i c d)

/-- Sixteen times the sharp independent-minus-retained energy. Both marginals
are averaged over the two actual native columns of the input word. -/
def sharpGain16 (k : ℤ) (i : Input) : ℤ :=
  mixedDiag k i false true + mixedDiag k i true false -
  mixedDiag k i false false - mixedDiag k i true true

def oddInput (i : Input) : Bool := (n i.1 + n i.2) % 2 == 1

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
theorem sharp_native : ∀ i : Input,
    sharpGain16 8 i = (if oddInput i then 34 else 0) ∧
    sharpGain16 0 i = (if oddInput i then 2 else 0) := by decide

theorem sharp_affine (k : ℤ) (i : Input) :
    sharpGain16 k i = sharpGain16 0 i + k*(sharpGain16 1 i-sharpGain16 0 i) := by
  unfold sharpGain16 mixedDiag h4
  ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
theorem sharp_unit : ∀ i : Input, sharpGain16 1 i = (if oddInput i then 6 else 0) := by
  decide

/-- All integer stiffnesses, including removal of the contrast cost. -/
theorem native_law (k : ℤ) (i : Input) :
    sharpGain16 k i = (if oddInput i then 4*k+2 else 0) := by
  rw [sharp_affine, (sharp_native i).2, sharp_unit]
  cases oddInput i <;> simp
  ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
/-- Coherent retained written inputs are in the zero-energy block. -/
theorem written_energy : ∀ i : Input,
    h4 8 (row i false) (row i false) + h4 8 (row i false) (row i true) +
    h4 8 (row i true) (row i false) + h4 8 (row i true) (row i true) = 0 := by decide

/-- A common layer threshold converts a number of unequal binary services to
the discrete total variation. Each integer threshold belongs to the levels
strictly above the lower endpoint and at or below the upper endpoint. -/
theorem layer_count (a b : ℕ) :
    (Finset.Ioc (min a b) (max a b)).card = max a b - min a b := by simp

theorem binary_contrast (a b : Bool) :
    oddInput ((if a then .single else .empty),(if b then .single else .empty)) =
      (a != b) := by cases a <;> cases b <;> decide

/-- Compact face descendants have source coefficient 3L for N=L²; equal-area
assemblies with different perimeter have different coefficients. -/
theorem face_scaling (L : ℝ) : ((17/8 : ℝ)*(3*L))^2 = (51/8)^2*L^2 := by ring

/-- At the ordinary degree-completed energy of a binary owner patch, the
one-layer stress fraction is universally bounded even at perfect delivery. -/
theorem source_ceiling (N b : ℝ) (hn : 0<N) (hb : 0≤b) (degree : b≤3*N) :
    (17/8*b)/(72*N+17/2*b) ≤ 17/260 := by
  apply (div_le_iff₀ (by positivity)).mpr
  linarith

/-- Independent courier routes that stop at most once can deposit only a
subprobability of coincident pairs over an entire physical chart. -/
theorem collision_budget {V : Type*} (S : Finset V) (p q : V → ℝ)
    (hp : ∀ x ∈ S, 0 ≤ p x) (hq : ∀ x ∈ S, q x ≤ 1)
    (stock : ∑ x ∈ S, p x ≤ 1) :
    ∑ x ∈ S, p x * q x ≤ 1 := by
  apply le_trans _ stock
  exact Finset.sum_le_sum fun x hx => by
    simpa using mul_le_mul_of_nonneg_left (hq x hx) (hp x hx)

/-- In an independent stream of writes, pairing the memory and loads with
independent routing hits with probability p*q preserves their connected
contrast by exactly that factor. -/
theorem transported_difference (paired separate p q : ℝ) :
    separate - (p*q*paired+(1-p*q)*separate) = p*q*(separate-paired) := by ring

end OPH.EinsteinBranch.NativeBoundary
