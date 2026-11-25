import Mathlib.Tactic

/-! Cardinality, directional occupation and conservation for a finite octant
readout. The index of a voxel in the layered graph is one plus its largest
coordinate. The separate injection vertex has layer zero. -/
namespace OPH.RoutedOctant

/-- A step increasing one coordinate crosses the boundary of the n-cube
exactly when it leaves the last coordinate plane with both transverse
coordinates inside the face. -/
theorem crosses_face (k i j n : ℕ) :
    (max k (max i j) < n ∧ n ≤ max (k+1) (max i j)) ↔
    (k+1 = n ∧ i < n ∧ j < n) := by omega

/-- Choice of exiting axis and the two transverse addresses. -/
theorem outward_addresses (n : ℕ) :
    Fintype.card (Fin 3 × Fin n × Fin n) = 3*n^2 := by
  simp [pow_two]

/-- Number of launch settings whose ray intersects a specified transverse
bin on a face. The same twelve settings propagate to every shell. -/
def rayBin (n b : ℕ) : ℕ :=
  ((List.range 12).filter fun i => n*(2*i+1)/24 == b).length

theorem midpoint_counts :
    (List.range 4).all (fun k => (List.range (k+1)).all
      (fun b => rayBin (k+1) b == 12/(k+1))) = true := by decide

theorem midpoint_count_at {n b : ℕ} (hn : 0 < n) (hd : n ≤ 4) (hb : b < n) :
    rayBin n b = 12/n := by
  have h := List.all_eq_true.mp midpoint_counts (n-1) (List.mem_range.mpr (by omega))
  have he : n - 1 + 1 = n := by omega
  rw [he] at h
  exact beq_iff_eq.mp (List.all_eq_true.mp h b (List.mem_range.mpr hb))

theorem pulse_count_per_edge {n i j : ℕ} (hn : 0 < n) (hd : n ≤ 4)
    (hi : i < n) (hj : j < n) :
    (3*n^2) * (rayBin n i * rayBin n j) = 432 := by
  rw [midpoint_count_at hn hd hi, midpoint_count_at hn hd hj]
  interval_cases n <;> decide

/-- Halving the shell-four transverse bin reproduces the shell-two bin
of the same ray, for every launch rank. -/
theorem coarse_fine_bin (i : ℕ) :
    (4*(2*i+1)/24)/2 = 2*(2*i+1)/24 := by omega

variable {α : Type*} [DecidableEq α]

def indicator (v a : α) : ℤ := if v = a then 1 else 0

def lastVertex (a : α) : List α → α
  | [] => a
  | b :: rest => lastVertex b rest

/-- One signed departure and arrival per observed hop. -/
def outflow (v a : α) : List α → ℤ
  | [] => 0
  | b :: rest => indicator v a - indicator v b + outflow v b rest

theorem path_balance (v a : α) (walk : List α) :
    outflow v a walk = indicator v a - indicator v (lastVertex a walk) := by
  induction walk generalizing a with
  | nil => simp [outflow, lastVertex]
  | cons b rest ih => simp only [outflow, lastVertex, ih]; ring

/-- Any interior station outside injection and collection has zero balance
on each completed pulse history, including paths with repeated visits. -/
theorem interior_balance (v a : α) (walk : List α)
    (hs : v ≠ a) (ht : v ≠ lastVertex a walk) : outflow v a walk = 0 := by
  simp [path_balance, indicator, hs, ht]

theorem rectilinear_census (q : ℕ) :
    Fintype.card (Fin q × Fin q × Fin q × Fin q) = q^4 := by
  simp [pow_succ, Nat.mul_assoc]

#print axioms crosses_face
#print axioms outward_addresses
#print axioms midpoint_count_at
#print axioms pulse_count_per_edge
#print axioms coarse_fine_bin
#print axioms path_balance
#print axioms interior_balance
#print axioms rectilinear_census
end OPH.RoutedOctant
