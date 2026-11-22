import ObserverPatchHolography.ObservablePriority.QuantumBridge
import ObserverPatchHolography.ObservablePriority.AddressedKernel

/-! Same finite screen and same addressed hit mask under two legal schedules.
    Output coordinates are final loads and a field of actual sign/coin collars.
    The proof supplies the classical premises from the physical sharp run. -/
namespace ObserverPatchSharp.FullChannel

open QuantumBridge
open scoped BigOperators

def signBit (z : Int) : StrictQuantum.Sign := decide (z = 2)

private theorem inverse_sign (v : Option StrictQuantum.Sign) :
    (v.map contrast).map signBit = v := by
  cases v with
  | none => rfl
  | some a => cases a <;> rfl

private theorem inverse_load (a : StrictQuantum.Load) : encode (decode a) = a := by
  cases a <;> rfl

theorem addressed_signature {n : Nat} (rank : Edge → Nat) (es : List (BEdge n))
    (h : Ordered rank (board es)) (x : StrictQuantum.Loads (Fin n)) :
    (StrictQuantum.scan (screen n) es x).2 = es.map (fun e =>
      ((sharpRun rank (board es) (source x)).collar e.val).map signBit) := by
  have ht := congrArg (List.map (Option.map signBit)) (signed_scan_sharp rank es x h).2
  simpa only [List.map_map, Function.comp_def, inverse_sign, List.map_id'] using ht

/-- Equality at every actual physical matrix coordinate of the branch Kraus
    operator. Duplicate freedom and local priority order are included in Ordered. -/
theorem actual_element_schedule_agreement {n : Nat} (rank : Edge → Nat)
    {as bs : List (BEdge n)} (p : as.Perm bs)
    (ha : Ordered rank (board as)) (hb : Ordered rank (board bs))
    (M : BEdge n → Bool) (x y : StrictQuantum.Loads (Fin n))
    (coll : BEdge n → StrictQuantum.Addressed.Collar) :
    StrictQuantum.Addressed.element (screen n) as M x y coll =
      StrictQuantum.Addressed.element (screen n) bs M x y coll := by
  have bp : (board as).Perm (board bs) :=
    (List.reverse_perm _).trans ((p.map Subtype.val).trans (List.reverse_perm _).symm)
  have hstate := sharpRun_eq_of_legal_perm rank (source x) bp ha hb
  let record := fun e : BEdge n =>
    ((sharpRun rank (board bs) (source x)).collar e.val).map signBit
  apply StrictQuantum.Addressed.element_schedule_agreement (screen n) p M x y coll record
  · funext v
    have h₁ := (signed_scan_sharp rank as x ha).1 v
    have h₂ := (signed_scan_sharp rank bs x hb).1 v
    rw [hstate] at h₁
    simpa only [inverse_load] using congrArg encode (h₁.symm.trans h₂)
  · simpa only [record, hstate] using addressed_signature rank as ha x
  · exact addressed_signature rank bs hb x

/-- The matrix product K_M ρ K_M† in the common physical output coordinates. -/
noncomputable def densityElement {n : Nat} (es : List (BEdge n)) (M : BEdge n → Bool)
    (ρ : StrictQuantum.Observed.density (Fin n)) (y z : StrictQuantum.Loads (Fin n))
    (cy cz : BEdge n → StrictQuantum.Addressed.Collar) : ℂ :=
  ∑ x : StrictQuantum.Loads (Fin n), ∑ x' : StrictQuantum.Loads (Fin n),
    (StrictQuantum.Addressed.element (screen n) es M x y cy : ℂ) * ρ x x' *
      star (StrictQuantum.Addressed.element (screen n) es M x' z cz : ℂ)

theorem actual_density_schedule_agreement {n : Nat} (rank : Edge → Nat)
    {as bs : List (BEdge n)} (p : as.Perm bs)
    (ha : Ordered rank (board as)) (hb : Ordered rank (board bs))
    (M : BEdge n → Bool) (ρ : StrictQuantum.Observed.density (Fin n))
    (y z : StrictQuantum.Loads (Fin n)) (cy cz : BEdge n → StrictQuantum.Addressed.Collar) :
    densityElement as M ρ y z cy cz = densityElement bs M ρ y z cy cz := by
  unfold densityElement
  apply Finset.sum_congr rfl
  intro x _
  apply Finset.sum_congr rfl
  intro x' _
  rw [actual_element_schedule_agreement rank p ha hb M x y cy,
      actual_element_schedule_agreement rank p ha hb M x' z cz]

#print axioms addressed_signature
#print axioms actual_element_schedule_agreement
#print axioms actual_density_schedule_agreement

end ObserverPatchSharp.FullChannel
