import Mathlib

open scoped BigOperators

namespace OPH.Z6Exact

abbrev Lattice := Fin 6 → ℤ

def total : Lattice →+ ℤ where
  toFun x := ∑ i, x i
  map_zero' := by simp
  map_add' x y := by simp [Finset.sum_add_distrib]

def ones : Lattice := fun _ ↦ 1

def diagonal : AddSubgroup Lattice := AddSubgroup.zmultiples ones

def balanced : AddSubgroup Lattice := total.ker

def gauge : AddSubgroup Lattice := diagonal ⊔ balanced

def residue : Lattice →+ ZMod 6 :=
  (Int.castAddHom (ZMod 6)).comp total

@[simp] theorem total_ones : total ones = 6 := by
  simp [total, ones]

@[simp] theorem total_zsmul_ones (k : ℤ) : total (k • ones) = 6 * k := by
  simp [total, ones]

theorem gauge_eq_kernel : gauge = residue.ker := by
  ext x
  constructor
  · intro hx
    rw [gauge, AddSubgroup.mem_sup] at hx
    obtain ⟨d, hd, b, hb, rfl⟩ := hx
    rw [AddMonoidHom.mem_ker]
    rw [diagonal, AddSubgroup.mem_zmultiples_iff] at hd
    obtain ⟨k, rfl⟩ := hd
    rw [balanced, AddMonoidHom.mem_ker] at hb
    change ((total (k • ones + b) : ℤ) : ZMod 6) = 0
    rw [map_add, total_zsmul_ones, hb]
    push_cast
    have h6 : (6 : ZMod 6) = 0 := ZMod.natCast_self 6
    simp [h6]
  · intro hx
    rw [AddMonoidHom.mem_ker] at hx
    change ((total x : ℤ) : ZMod 6) = 0 at hx
    have hdiv : (6 : ℤ) ∣ total x :=
      (ZMod.intCast_zmod_eq_zero_iff_dvd (total x) 6).mp hx
    obtain ⟨k, hk⟩ := hdiv
    rw [gauge, AddSubgroup.mem_sup]
    refine ⟨k • ones, ?_, x - k • ones, ?_, by abel⟩
    · exact (AddSubgroup.mem_zmultiples_iff).2 ⟨k, rfl⟩
    · rw [balanced, AddMonoidHom.mem_ker]
      rw [map_sub, total_zsmul_ones, hk]
      ring

theorem residue_surjective : Function.Surjective residue := by
  intro z
  obtain ⟨k, rfl⟩ := ZMod.intCast_surjective z
  refine ⟨fun i ↦ if i = 0 then k else 0, ?_⟩
  simp [residue, total]

noncomputable def quotientEquivZMod6 : Lattice ⧸ gauge ≃+ ZMod 6 :=
  (QuotientAddGroup.quotientAddEquivOfEq gauge_eq_kernel).trans
    (QuotientAddGroup.quotientKerEquivOfSurjective residue residue_surjective)

theorem representative_formula (x : Lattice) :
    quotientEquivZMod6 (QuotientAddGroup.mk x) = residue x := by
  rfl

#print axioms gauge_eq_kernel
#print axioms residue_surjective
#print axioms quotientEquivZMod6
#print axioms representative_formula

end OPH.Z6Exact
