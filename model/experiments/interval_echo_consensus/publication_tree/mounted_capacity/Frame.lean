import Mathlib.GroupTheory.GroupAction.Basic

namespace OPH.MountedCapacity

variable {G X M : Type*} [Group G] [MulAction G X] [MulAction G M]

/-- Without a common mount label, two launches are confusable when they can give
the same payload and reference readings. -/
def Shared (G : Type*) [Group G] [MulAction G X] [MulAction G M] (m : M) (x y : X) : Prop :=
  ∃ g h : G, g • x = h • y ∧ g • m = h • m

theorem shared_iff (m : M) (x y : X) :
    Shared G m x y ↔ ∃ k : G, k • m = m ∧ k • x = y := by
  constructor
  · rintro ⟨g,h,hx,hm⟩
    exact ⟨h⁻¹*g, by rw [mul_smul,hm,inv_smul_smul], by rw [mul_smul,hx,inv_smul_smul]⟩
  · rintro ⟨k,hm,hx⟩
    exact ⟨k,1,by simpa using hx,by simpa using hm⟩

/-- A reference with trivial stabilizer resolves the entire input. -/
theorem faithful_reference (m : M) (free : ∀ k : G, k • m = m → k = 1) (x y : X) :
    Shared G m x y ↔ x = y := by
  rw [shared_iff]
  constructor
  · rintro ⟨k,hm,hx⟩
    simpa [free k hm] using hx
  · rintro rfl
    exact ⟨1,one_smul _ _,one_smul _ _⟩

/-- If the reference can be remounted independently for each launch, it
cannot remove a payload ambiguity. -/
theorem independent_reference (m : M) (x y : X) :
    (∃ g h u v : G, g • x = h • y ∧ u • m = v • m) ↔ ∃ g h : G, g • x = h • y := by
  constructor
  · rintro ⟨g,h,_,_,hx,_⟩
    exact ⟨g,h,hx⟩
  · rintro ⟨g,h,hx⟩
    exact ⟨g,h,1,1,hx,rfl⟩

theorem duplicate_readout {T : Type*} (a b : T) : (a,a) = (b,b) ↔ a = b :=
  ⟨fun h => congrArg Prod.fst h,fun h => h ▸ rfl⟩

/-- A veto has a public implementation exactly when equal observed records force equal decisions. -/
theorem decision_factors {Run R : Type*} (read : Run → R) (veto : Run → Bool) :
    (∃ f : R → Bool, ∀ x, f (read x) = veto x) ↔
      ∀ x y, read x = read y → veto x = veto y := by
  classical
  constructor
  · rintro ⟨f,hf⟩ x y h
    rw [← hf x, ← hf y, h]
  · intro h
    refine ⟨fun r => if w : ∃ x, read x = r then veto w.choose else false, ?_⟩
    intro x
    have available : ∃ y, read y = read x := ⟨x,rfl⟩
    simp only [dif_pos available]
    exact h _ _ available.choose_spec

end OPH.MountedCapacity
