import Mathlib
import A5PortAction

namespace OPH.A5PortModule

/-! # The twelve-port permutation module and the S5 receipt

Issue #604, second target.  Step S5 of the selection skeleton
(`Screen/A5OPH.lean`, Part IV) needs: trivial-isotypic multiplicity one in
the twelve-port module `P12` forces `dim Z <= 1` under the
central-triviality hypothesis.  This file supplies that receipt on the
concrete port module.

CONTENT.

* `P60`: the sixty proper icosahedral rotations as permutations of
  `Fin 12`, row-for-row identical to the vetted incidence-automorphism
  list `OPH.A5PortAction.perms` (`P60_matches_port_rows`), so the closure,
  inverse, adjacency-preservation, antipode-equivariance, and transitivity
  receipts of that file apply to this list verbatim.
* `Fixed`: the subspace of port fields fixed by every listed rotation.
* `fixed_eq_span`, `finrank_fixed`: the fixed subspace is exactly the line
  of constant port fields — trivial-isotypic multiplicity one, forced by
  transitivity alone (`exists_map_zero`).
* `centrally_trivial_finrank_le_one`: the S5 receipt.  Any submodule on
  which every listed rotation acts as the identity — in particular the
  centre of a compact current algebra when the icosahedral action is
  inner, since inner automorphisms fix the centre pointwise — has
  dimension at most one.

BOUNDARY.  The companion half of S5 (`dim Z /= 0`, the `su(2)^4`
exclusion) uses the compact-simple dimension arithmetic and the
triviality of `A5`-actions on at most four objects, proved in
`Screen/A5OPH.lean` (S2, S4); the Lie-theoretic inputs (reductive
decomposition, compact-simple dimension list) remain declared classical
inputs there.  No physical content is at stake in this receipt lane. -/

/-- Build a permutation of `Fin 12` from its value vector and inverse
vector; the two inverse laws are discharged by `decide` at each literal. -/
def perm (f g : Fin 12 → Fin 12)
    (h1 : Function.LeftInverse g f := by decide)
    (h2 : Function.RightInverse g f := by decide) : Equiv.Perm (Fin 12) :=
  ⟨f, g, h1, h2⟩

/-- The sixty proper icosahedral rotations on the twelve ports, as
permutations of `Fin 12`. -/
def P60 : List (Equiv.Perm (Fin 12)) := [
    perm ![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] ![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    perm ![0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11] ![0, 3, 1, 6, 2, 7, 4, 9, 5, 10, 8, 11],
    perm ![0, 3, 1, 6, 2, 7, 4, 9, 5, 10, 8, 11] ![0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11],
    perm ![0, 4, 6, 2, 3, 10, 1, 8, 9, 5, 7, 11] ![0, 6, 3, 4, 1, 9, 2, 10, 7, 8, 5, 11],
    perm ![0, 6, 3, 4, 1, 9, 2, 10, 7, 8, 5, 11] ![0, 4, 6, 2, 3, 10, 1, 8, 9, 5, 7, 11],
    perm ![1, 0, 3, 2, 7, 6, 5, 4, 9, 8, 11, 10] ![1, 0, 3, 2, 7, 6, 5, 4, 9, 8, 11, 10],
    perm ![1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10] ![2, 0, 1, 4, 5, 3, 8, 6, 7, 10, 11, 9],
    perm ![1, 3, 7, 0, 5, 9, 2, 6, 11, 4, 8, 10] ![3, 0, 6, 1, 9, 4, 7, 2, 10, 5, 11, 8],
    perm ![1, 5, 2, 7, 0, 8, 3, 11, 4, 9, 6, 10] ![4, 0, 2, 6, 8, 1, 10, 3, 5, 9, 11, 7],
    perm ![1, 7, 5, 3, 2, 11, 0, 9, 8, 6, 4, 10] ![6, 0, 4, 3, 10, 2, 9, 1, 8, 7, 11, 5],
    perm ![2, 0, 1, 4, 5, 3, 8, 6, 7, 10, 11, 9] ![1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10],
    perm ![2, 1, 5, 0, 8, 7, 4, 3, 11, 6, 10, 9] ![3, 1, 0, 7, 6, 2, 9, 5, 4, 11, 10, 8],
    perm ![2, 4, 0, 8, 1, 6, 5, 10, 3, 11, 7, 9] ![2, 4, 0, 8, 1, 6, 5, 10, 3, 11, 7, 9],
    perm ![2, 5, 8, 1, 4, 11, 0, 7, 10, 3, 6, 9] ![6, 3, 0, 9, 4, 1, 10, 7, 2, 11, 8, 5],
    perm ![2, 8, 4, 5, 0, 10, 1, 11, 6, 7, 3, 9] ![4, 6, 0, 10, 2, 3, 8, 9, 1, 11, 5, 7],
    perm ![3, 0, 6, 1, 9, 4, 7, 2, 10, 5, 11, 8] ![1, 3, 7, 0, 5, 9, 2, 6, 11, 4, 8, 10],
    perm ![3, 1, 0, 7, 6, 2, 9, 5, 4, 11, 10, 8] ![2, 1, 5, 0, 8, 7, 4, 3, 11, 6, 10, 9],
    perm ![3, 6, 9, 0, 7, 10, 1, 4, 11, 2, 5, 8] ![3, 6, 9, 0, 7, 10, 1, 4, 11, 2, 5, 8],
    perm ![3, 7, 1, 9, 0, 5, 6, 11, 2, 10, 4, 8] ![4, 2, 8, 0, 10, 5, 6, 1, 11, 3, 9, 7],
    perm ![3, 9, 7, 6, 1, 11, 0, 10, 5, 4, 2, 8] ![6, 4, 10, 0, 9, 8, 3, 2, 11, 1, 7, 5],
    perm ![4, 0, 2, 6, 8, 1, 10, 3, 5, 9, 11, 7] ![1, 5, 2, 7, 0, 8, 3, 11, 4, 9, 6, 10],
    perm ![4, 2, 8, 0, 10, 5, 6, 1, 11, 3, 9, 7] ![3, 7, 1, 9, 0, 5, 6, 11, 2, 10, 4, 8],
    perm ![4, 6, 0, 10, 2, 3, 8, 9, 1, 11, 5, 7] ![2, 8, 4, 5, 0, 10, 1, 11, 6, 7, 3, 9],
    perm ![4, 8, 10, 2, 6, 11, 0, 5, 9, 1, 3, 7] ![6, 9, 3, 10, 0, 7, 4, 11, 1, 8, 2, 5],
    perm ![4, 10, 6, 8, 0, 9, 2, 11, 3, 5, 1, 7] ![4, 10, 6, 8, 0, 9, 2, 11, 3, 5, 1, 7],
    perm ![5, 1, 7, 2, 11, 3, 8, 0, 9, 4, 10, 6] ![7, 1, 3, 5, 9, 0, 11, 2, 6, 8, 10, 4],
    perm ![5, 2, 1, 8, 7, 0, 11, 4, 3, 10, 9, 6] ![5, 2, 1, 8, 7, 0, 11, 4, 3, 10, 9, 6],
    perm ![5, 7, 11, 1, 8, 9, 2, 3, 10, 0, 4, 6] ![9, 3, 6, 7, 10, 0, 11, 1, 4, 5, 8, 2],
    perm ![5, 8, 2, 11, 1, 4, 7, 10, 0, 9, 3, 6] ![8, 4, 2, 10, 5, 0, 11, 6, 1, 9, 7, 3],
    perm ![5, 11, 8, 7, 2, 10, 1, 9, 4, 3, 0, 6] ![10, 6, 4, 9, 8, 0, 11, 3, 2, 7, 5, 1],
    perm ![6, 0, 4, 3, 10, 2, 9, 1, 8, 7, 11, 5] ![1, 7, 5, 3, 2, 11, 0, 9, 8, 6, 4, 10],
    perm ![6, 3, 0, 9, 4, 1, 10, 7, 2, 11, 8, 5] ![2, 5, 8, 1, 4, 11, 0, 7, 10, 3, 6, 9],
    perm ![6, 4, 10, 0, 9, 8, 3, 2, 11, 1, 7, 5] ![3, 9, 7, 6, 1, 11, 0, 10, 5, 4, 2, 8],
    perm ![6, 9, 3, 10, 0, 7, 4, 11, 1, 8, 2, 5] ![4, 8, 10, 2, 6, 11, 0, 5, 9, 1, 3, 7],
    perm ![6, 10, 9, 4, 3, 11, 0, 8, 7, 2, 1, 5] ![6, 10, 9, 4, 3, 11, 0, 8, 7, 2, 1, 5],
    perm ![7, 1, 3, 5, 9, 0, 11, 2, 6, 8, 10, 4] ![5, 1, 7, 2, 11, 3, 8, 0, 9, 4, 10, 6],
    perm ![7, 3, 9, 1, 11, 6, 5, 0, 10, 2, 8, 4] ![7, 3, 9, 1, 11, 6, 5, 0, 10, 2, 8, 4],
    perm ![7, 5, 1, 11, 3, 2, 9, 8, 0, 10, 6, 4] ![8, 2, 5, 4, 11, 1, 10, 0, 7, 6, 9, 3],
    perm ![7, 9, 11, 3, 5, 10, 1, 6, 8, 0, 2, 4] ![9, 6, 10, 3, 11, 4, 7, 0, 8, 1, 5, 2],
    perm ![7, 11, 5, 9, 1, 8, 3, 10, 2, 6, 0, 4] ![10, 4, 8, 6, 11, 2, 9, 0, 5, 3, 7, 1],
    perm ![8, 2, 5, 4, 11, 1, 10, 0, 7, 6, 9, 3] ![7, 5, 1, 11, 3, 2, 9, 8, 0, 10, 6, 4],
    perm ![8, 4, 2, 10, 5, 0, 11, 6, 1, 9, 7, 3] ![5, 8, 2, 11, 1, 4, 7, 10, 0, 9, 3, 6],
    perm ![8, 5, 11, 2, 10, 7, 4, 1, 9, 0, 6, 3] ![9, 7, 3, 11, 6, 1, 10, 5, 0, 8, 4, 2],
    perm ![8, 10, 4, 11, 2, 6, 5, 9, 0, 7, 1, 3] ![8, 10, 4, 11, 2, 6, 5, 9, 0, 7, 1, 3],
    perm ![8, 11, 10, 5, 4, 9, 2, 7, 6, 1, 0, 3] ![10, 9, 6, 11, 4, 3, 8, 7, 0, 5, 2, 1],
    perm ![9, 3, 6, 7, 10, 0, 11, 1, 4, 5, 8, 2] ![5, 7, 11, 1, 8, 9, 2, 3, 10, 0, 4, 6],
    perm ![9, 6, 10, 3, 11, 4, 7, 0, 8, 1, 5, 2] ![7, 9, 11, 3, 5, 10, 1, 6, 8, 0, 2, 4],
    perm ![9, 7, 3, 11, 6, 1, 10, 5, 0, 8, 4, 2] ![8, 5, 11, 2, 10, 7, 4, 1, 9, 0, 6, 3],
    perm ![9, 10, 11, 6, 7, 8, 3, 4, 5, 0, 1, 2] ![9, 10, 11, 6, 7, 8, 3, 4, 5, 0, 1, 2],
    perm ![9, 11, 7, 10, 3, 5, 6, 8, 1, 4, 0, 2] ![10, 8, 11, 4, 9, 5, 6, 2, 7, 0, 3, 1],
    perm ![10, 4, 8, 6, 11, 2, 9, 0, 5, 3, 7, 1] ![7, 11, 5, 9, 1, 8, 3, 10, 2, 6, 0, 4],
    perm ![10, 6, 4, 9, 8, 0, 11, 3, 2, 7, 5, 1] ![5, 11, 8, 7, 2, 10, 1, 9, 4, 3, 0, 6],
    perm ![10, 8, 11, 4, 9, 5, 6, 2, 7, 0, 3, 1] ![9, 11, 7, 10, 3, 5, 6, 8, 1, 4, 0, 2],
    perm ![10, 9, 6, 11, 4, 3, 8, 7, 0, 5, 2, 1] ![8, 11, 10, 5, 4, 9, 2, 7, 6, 1, 0, 3],
    perm ![10, 11, 9, 8, 6, 7, 4, 5, 3, 2, 0, 1] ![10, 11, 9, 8, 6, 7, 4, 5, 3, 2, 0, 1],
    perm ![11, 5, 7, 8, 9, 1, 10, 2, 3, 4, 6, 0] ![11, 5, 7, 8, 9, 1, 10, 2, 3, 4, 6, 0],
    perm ![11, 7, 9, 5, 10, 3, 8, 1, 6, 2, 4, 0] ![11, 7, 9, 5, 10, 3, 8, 1, 6, 2, 4, 0],
    perm ![11, 8, 5, 10, 7, 2, 9, 4, 1, 6, 3, 0] ![11, 8, 5, 10, 7, 2, 9, 4, 1, 6, 3, 0],
    perm ![11, 9, 10, 7, 8, 6, 5, 3, 4, 1, 2, 0] ![11, 9, 10, 7, 8, 6, 5, 3, 4, 1, 2, 0],
    perm ![11, 10, 8, 9, 5, 4, 7, 6, 2, 3, 1, 0] ![11, 10, 8, 9, 5, 4, 7, 6, 2, 3, 1, 0]
  ]

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Row-for-row agreement with the vetted incidence-automorphism list. -/
theorem P60_matches_port_rows :
    P60.map (fun g => (List.finRange 12).map fun k => (g k : Nat)) =
      OPH.A5PortAction.perms := by decide

set_option maxHeartbeats 2000000 in
/-- Transitivity witness: every port is the image of port `0`. -/
theorem exists_map_zero : ∀ i : Fin 12, ∃ g ∈ P60, g 0 = i := by decide

/-- A permutation acts linearly on port fields by right composition with
its inverse. -/
def pAct (g : Equiv.Perm (Fin 12)) : (Fin 12 → ℚ) →ₗ[ℚ] (Fin 12 → ℚ) where
  toFun v i := v (g⁻¹ i)
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

/-- The subspace of port fields fixed by every listed rotation. -/
def Fixed : Submodule ℚ (Fin 12 → ℚ) where
  carrier := {v | ∀ g ∈ P60, pAct g v = v}
  add_mem' := by
    intro v w hv hw g hg
    rw [(pAct g).map_add v w, hv g hg, hw g hg]
  zero_mem' := by
    intro g hg
    exact (pAct g).map_zero
  smul_mem' := by
    intro c v hv g hg
    rw [(pAct g).map_smul c v, hv g hg]

/-- The constant port field. -/
def ones : Fin 12 → ℚ := fun _ => 1

theorem ones_mem_fixed : ones ∈ Fixed := fun _ _ => rfl

theorem ones_ne_zero : ones ≠ 0 := by
  intro h
  exact one_ne_zero (congrFun h 0)

/-- Trivial-isotypic multiplicity one: the fixed subspace is exactly the
line of constant port fields, by transitivity alone. -/
theorem fixed_eq_span : Fixed = Submodule.span ℚ {ones} := by
  apply le_antisymm
  · intro v hv
    rw [Submodule.mem_span_singleton]
    refine ⟨v 0, funext fun i => ?_⟩
    obtain ⟨g, hg, hgi⟩ := exists_map_zero i
    have hfix : v (g⁻¹ i) = v i := congrFun (hv g hg) i
    have h0 : g⁻¹ i = 0 := by
      rw [← hgi]
      exact g.symm_apply_apply 0
    have : v i = v 0 := by rw [← hfix, h0]
    simp [ones, this]
  · rw [Submodule.span_le, Set.singleton_subset_iff]
    exact ones_mem_fixed

theorem finrank_fixed : Module.finrank ℚ Fixed = 1 := by
  rw [fixed_eq_span]
  exact finrank_span_singleton ones_ne_zero

/-- S5 receipt: any submodule of the port module on which every listed
rotation acts as the identity has dimension at most one. -/
theorem centrally_trivial_finrank_le_one (Z : Submodule ℚ (Fin 12 → ℚ))
    (hZ : ∀ g ∈ P60, ∀ z ∈ Z, pAct g z = z) :
    Module.finrank ℚ Z ≤ 1 := by
  have hle : Z ≤ Fixed := fun z hz g hg => hZ g hg z hz
  have h := Submodule.finrank_mono hle
  rwa [finrank_fixed] at h

end OPH.A5PortModule

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5PortModule.P60_matches_port_rows
#print axioms OPH.A5PortModule.exists_map_zero
#print axioms OPH.A5PortModule.fixed_eq_span
#print axioms OPH.A5PortModule.finrank_fixed
#print axioms OPH.A5PortModule.centrally_trivial_finrank_le_one
