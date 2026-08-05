import Geometry.LorentzOverlapCocycle

/-!
# C2: event-germ coordinates and displacement covariance

An `EventGermAtlas` is the exact handoff required before the intrinsic C1
Lorentz module can be soldered to event-like records.  It contains chart
visibility, one coordinate readback, affine transitions, and the overlap
compatibility law.  The principal theorem proves that coordinate differences
are translation-blind and transform covariantly.

The atlas is supplied data.  This module does not construct it from the three
OPH axioms or assert population, separation, topology, causal reachability, or
physical distance.
-/

namespace OPH.C2Soldering

noncomputable section

open OPH.C1Lorentz

/-! ## Coincidence descent before charting -/

/-- A raw readback descends through a declared coincidence setoid exactly
when it is invariant on coincidence classes. -/
def CoincidenceInvariant {Raw : Type*} (s : Setoid Raw)
    (readback : Raw → Herm2) : Prop :=
  ∀ a b, s.r a b → readback a = readback b

/-- Coordinate readback on coincidence classes. -/
def descendedReadback {Raw : Type*} (s : Setoid Raw)
    (readback : Raw → Herm2) (h : CoincidenceInvariant s readback) :
    Quotient s → Herm2 :=
  Quotient.lift readback h

@[simp] theorem descendedReadback_mk {Raw : Type*} (s : Setoid Raw)
    (readback : Raw → Herm2) (h : CoincidenceInvariant s readback)
    (x : Raw) :
    descendedReadback s readback h (Quotient.mk s x) = readback x := rfl

/-- Exact descent criterion, including uniqueness of the chart readback on
the quotient.  The theorem requires an actual setoid; a merely reflexive and
symmetric pairwise-overlap test is not silently closed transitively. -/
theorem coincidenceInvariant_iff_existsUnique_descendedReadback
    {Raw : Type*} (s : Setoid Raw) (readback : Raw → Herm2) :
    CoincidenceInvariant s readback ↔
      ∃! f : Quotient s → Herm2,
        ∀ x, f (Quotient.mk s x) = readback x := by
  constructor
  · intro h
    refine ⟨descendedReadback s readback h, by simp, ?_⟩
    intro f hf
    funext q
    refine Quotient.inductionOn q ?_
    intro x
    rw [hf]
    rfl
  · rintro ⟨f, hf, _⟩
    intro a b hab
    rw [← hf a, ← hf b]
    exact congrArg f (Quotient.sound hab)

/-! A three-germ negative control: pairwise/cofinal-style overlap can be
reflexive and symmetric without being transitive. -/

def overlapControl (a b : Fin 3) : Prop :=
  a = b ∨
    (a = 0 ∧ b = 1) ∨ (a = 1 ∧ b = 0) ∨
    (a = 1 ∧ b = 2) ∨ (a = 2 ∧ b = 1)

theorem overlapControl_refl (a : Fin 3) : overlapControl a a :=
  Or.inl rfl

theorem overlapControl_symm {a b : Fin 3} (h : overlapControl a b) :
    overlapControl b a := by
  rcases h with h | h | h | h | h
  · exact Or.inl h.symm
  · exact Or.inr (Or.inr (Or.inl ⟨h.2, h.1⟩))
  · exact Or.inr (Or.inl ⟨h.2, h.1⟩)
  · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨h.2, h.1⟩)))
  · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨h.2, h.1⟩)))

theorem overlapControl_not_transitive :
    overlapControl 0 1 ∧ overlapControl 1 2 ∧ ¬overlapControl 0 2 := by
  refine ⟨Or.inr (Or.inl ⟨rfl, rfl⟩),
    Or.inr (Or.inr (Or.inr (Or.inl ⟨rfl, rfl⟩))), ?_⟩
  intro h
  rcases h with h | h | h | h | h <;> omega

/-- A charted record-germ interface over the intrinsic Hermitian Lorentz
module. -/
structure EventGermAtlas (Event Chart : Type*) where
  visible : Chart → Event → Prop
  coordinate : Chart → Event → Herm2
  overlap : LorentzOverlapCocycle Chart
  coordinate_overlap : ∀ {i j e}, visible i e → visible j e →
    coordinate j e = overlap.act i j (coordinate i e)

/-- Coordinate displacement from `p` to `q` in one visible chart. -/
def EventGermAtlas.displacement {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p q : Event) : Herm2 :=
  A.coordinate i q - A.coordinate i p

@[simp] theorem EventGermAtlas.displacement_self {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p : Event) :
    A.displacement i p p = 0 := by
  simp [EventGermAtlas.displacement]

theorem EventGermAtlas.displacement_reverse {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p q : Event) :
    A.displacement i q p = -A.displacement i p q := by
  simp [EventGermAtlas.displacement]

theorem EventGermAtlas.displacement_chain {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (i : Chart) (p q r : Event) :
    A.displacement i p q + A.displacement i q r =
      A.displacement i p r := by
  simp [EventGermAtlas.displacement]

/-- Exact overlap covariance of event-germ displacement.  This is where the
affine translation cancels. -/
theorem EventGermAtlas.displacement_overlap {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) {i j : Chart} {p q : Event}
    (hip : A.visible i p) (hiq : A.visible i q)
    (hjp : A.visible j p) (hjq : A.visible j q) :
    A.displacement j p q =
      A.overlap.lorentz i j (A.displacement i p q) := by
  have hq := A.coordinate_overlap (i := i) (j := j) (e := q) hiq hjq
  have hp := A.coordinate_overlap (i := i) (j := j) (e := p) hip hjp
  rw [EventGermAtlas.displacement, hq, hp]
  exact A.overlap.act_sub_act i j (A.coordinate i p) (A.coordinate i q)

/-- Lorentz intervals of visible displacements are chart independent. -/
theorem EventGermAtlas.interval_overlap {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) {i j : Chart} {p q : Event}
    (hip : A.visible i p) (hiq : A.visible i q)
    (hjp : A.visible j p) (hjq : A.visible j q) :
    lorentzQ (A.displacement j p q) =
      lorentzQ (A.displacement i p q) := by
  rw [A.displacement_overlap hip hiq hjp hjq]
  simpa only [← lorentzB_self] using
    A.overlap.lorentz i j |>.map_lorentzB
      (A.displacement i p q) (A.displacement i p q)

#print axioms EventGermAtlas.displacement_overlap
#print axioms EventGermAtlas.interval_overlap
#print axioms coincidenceInvariant_iff_existsUnique_descendedReadback
#print axioms overlapControl_not_transitive

end

end OPH.C2Soldering
