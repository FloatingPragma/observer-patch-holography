import Geometry.ObserverRestSpace

/-!
# C2: oriented Lorentz and affine-overlap cocycles

This file gives the algebraic transition object used by the event-frame
soldering lane.  An `OrientedLorentzEquiv` is deliberately stronger than a
bare linear equivalence: preservation of the Lorentz pairing and of the two
future-oriented loci is explicit data.  A `LorentzOverlapCocycle` then adds
affine translations and the exact identity and triple-overlap laws.

The transition data are a contract.  Nothing here proves that observer-patch
records produce charts, transitions, an event population, a topology, or a
physical spacetime.
-/

namespace OPH.C2Soldering

noncomputable section

open OPH.C1Lorentz

/-- A time-orientation-preserving Lorentz linear equivalence on the intrinsic
`Herm2` module.  The future-null and future-unit clauses are stated as
equivalences so that no hidden surjectivity/orientation argument is needed
when transitions are inverted or composed. -/
structure OrientedLorentzEquiv where
  toLinearEquiv : Herm2 ≃ₗ[ℝ] Herm2
  map_lorentzB : ∀ u v, lorentzB (toLinearEquiv u) (toLinearEquiv v) = lorentzB u v
  map_futureNull_iff : ∀ v, IsFutureNull (toLinearEquiv v) ↔ IsFutureNull v
  map_futureUnit_iff : ∀ v, IsFutureUnitTimelike (toLinearEquiv v) ↔
    IsFutureUnitTimelike v

instance : CoeFun OrientedLorentzEquiv (fun _ => Herm2 → Herm2) :=
  ⟨fun L => L.toLinearEquiv⟩

@[simp] theorem OrientedLorentzEquiv.map_add (L : OrientedLorentzEquiv)
    (u v : Herm2) : L (u + v) = L u + L v :=
  L.toLinearEquiv.map_add u v

@[simp] theorem OrientedLorentzEquiv.map_sub (L : OrientedLorentzEquiv)
    (u v : Herm2) : L (u - v) = L u - L v :=
  L.toLinearEquiv.map_sub u v

@[simp] theorem OrientedLorentzEquiv.map_smul (L : OrientedLorentzEquiv)
    (a : ℝ) (v : Herm2) : L (a • v) = a • L v :=
  L.toLinearEquiv.map_smul a v

/-- Identity oriented Lorentz map. -/
def OrientedLorentzEquiv.refl : OrientedLorentzEquiv where
  toLinearEquiv := LinearEquiv.refl ℝ Herm2
  map_lorentzB := by simp
  map_futureNull_iff := by simp
  map_futureUnit_iff := by simp

/-- Composition, in the order "first `L`, then `M`". -/
def OrientedLorentzEquiv.trans (L M : OrientedLorentzEquiv) :
    OrientedLorentzEquiv where
  toLinearEquiv := L.toLinearEquiv.trans M.toLinearEquiv
  map_lorentzB := by
    intro u v
    change lorentzB (M (L u)) (M (L v)) = lorentzB u v
    rw [M.map_lorentzB, L.map_lorentzB]
  map_futureNull_iff := by
    intro v
    rw [LinearEquiv.trans_apply, M.map_futureNull_iff, L.map_futureNull_iff]
  map_futureUnit_iff := by
    intro v
    rw [LinearEquiv.trans_apply, M.map_futureUnit_iff, L.map_futureUnit_iff]

/-- Transport a future-null vector through an oriented Lorentz map. -/
def OrientedLorentzEquiv.mapFutureNullVector (L : OrientedLorentzEquiv)
    (v : FutureNullVector) : FutureNullVector :=
  ⟨L v.1, (L.map_futureNull_iff v.1).2 v.2⟩

/-- Transport an observer frame through an oriented Lorentz map. -/
def OrientedLorentzEquiv.mapFrame (L : OrientedLorentzEquiv)
    (u : FrameHyperboloid) : FrameHyperboloid :=
  ⟨L u.1, (L.map_futureUnit_iff u.1).2 u.2⟩

@[simp] theorem OrientedLorentzEquiv.mapFrame_coe
    (L : OrientedLorentzEquiv) (u : FrameHyperboloid) :
    (L.mapFrame u : Herm2) = L (u : Herm2) := rfl

/-- An affine Lorentz transition between two chart labels. -/
structure LorentzOverlapCocycle (Chart : Type*) where
  lorentz : Chart → Chart → OrientedLorentzEquiv
  translation : Chart → Chart → Herm2
  lorentz_self : ∀ i x, lorentz i i x = x
  lorentz_cocycle : ∀ i j k x,
    lorentz i k x = lorentz j k (lorentz i j x)
  translation_self : ∀ i, translation i i = 0
  translation_cocycle : ∀ i j k,
    translation i k =
      lorentz j k (translation i j) + translation j k

/-- Affine action of the transition from chart `i` to chart `j`. -/
def LorentzOverlapCocycle.act {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i j : Chart) (x : Herm2) : Herm2 :=
  T.lorentz i j x + T.translation i j

@[simp] theorem LorentzOverlapCocycle.act_self {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i : Chart) (x : Herm2) :
    T.act i i x = x := by
  simp [LorentzOverlapCocycle.act, T.lorentz_self, T.translation_self]

/-- Exact affine triple-overlap cocycle. -/
theorem LorentzOverlapCocycle.act_cocycle {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i j k : Chart) (x : Herm2) :
    T.act i k x = T.act j k (T.act i j x) := by
  rw [LorentzOverlapCocycle.act, T.lorentz_cocycle,
    T.translation_cocycle]
  simp only [LorentzOverlapCocycle.act, OrientedLorentzEquiv.map_add]
  abel

/-- The reverse chart transition is a left inverse, derived from the cocycle
and identity laws rather than assumed separately. -/
theorem LorentzOverlapCocycle.act_reverse_left {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i j : Chart) (x : Herm2) :
    T.act j i (T.act i j x) = x := by
  rw [← T.act_cocycle i j i x, T.act_self]

/-- The reverse chart transition is also a right inverse. -/
theorem LorentzOverlapCocycle.act_reverse_right {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i j : Chart) (x : Herm2) :
    T.act i j (T.act j i x) = x := by
  rw [← T.act_cocycle j i j x, T.act_self]

/-- The affine cocycle transports a coordinate difference only through its
Lorentz-linear part; translations cancel exactly. -/
theorem LorentzOverlapCocycle.act_sub_act {Chart : Type*}
    (T : LorentzOverlapCocycle Chart) (i j : Chart) (x y : Herm2) :
    T.act i j y - T.act i j x = T.lorentz i j (y - x) := by
  simp only [LorentzOverlapCocycle.act, OrientedLorentzEquiv.map_sub]
  abel

#print axioms LorentzOverlapCocycle.act_cocycle
#print axioms LorentzOverlapCocycle.act_reverse_left
#print axioms LorentzOverlapCocycle.act_sub_act

end

end OPH.C2Soldering
