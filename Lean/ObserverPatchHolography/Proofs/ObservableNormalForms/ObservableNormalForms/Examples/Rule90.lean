import Mathlib

/-!
# Width-three Rule 90 regression example

Self-contained proofs of the exact finite claims used in the manuscript's
Rule-90 example: kernel, image, good and deficient coordinate readouts, and
the obstruction to a total reverse repair that preserves an arbitrary target
row while changing only the seed.
-/

namespace ObservableNormalForms.Rule90

abbrev Row := Bool × Bool × Bool

/-- One Rule-90 step on a width-three row with zero exterior boundary. -/
def step (s : Row) : Row :=
  (s.2.1, Bool.xor s.1 s.2.2, s.2.1)

def zeroRow : Row := (false, false, false)

def kernelGenerator : Row := (true, false, true)

def read01 (t : Row) : Bool × Bool := (t.1, t.2.1)

def read02 (t : Row) : Bool × Bool := (t.1, t.2.2)

/-- Exact two-element kernel, equivalently the Boolean-linear span of
`(1,0,1)`. -/
theorem kernel_exact (s : Row) :
    step s = zeroRow ↔ s = zeroRow ∨ s = kernelGenerator := by
  rcases s with ⟨a, b, c⟩
  cases a <;> cases b <;> cases c <;> decide

theorem outer_coordinates_equal (s : Row) :
    (step s).1 = (step s).2.2 := rfl

/-- Exact image characterization: precisely the rows `(u,v,u)`. -/
theorem mem_range_iff_outer_coordinates_equal (t : Row) :
    t ∈ Set.range step ↔ t.1 = t.2.2 := by
  constructor
  · rintro ⟨s, rfl⟩
    exact outer_coordinates_equal s
  · intro h
    rcases t with ⟨u, v, w⟩
    simp only at h
    subst w
    exact ⟨(false, u, v), by simp [step]⟩

theorem image_exact :
    Set.range step = {t : Row | t.1 = t.2.2} := by
  ext t
  exact mem_range_iff_outer_coordinates_equal t

/-- Reading coordinates 0 and 1 determines an image row. -/
theorem read01_injective_on_image : Set.InjOn read01 (Set.range step) := by
  intro x hx y hy hread
  have hxOuter : x.1 = x.2.2 :=
    (mem_range_iff_outer_coordinates_equal x).mp hx
  have hyOuter : y.1 = y.2.2 :=
    (mem_range_iff_outer_coordinates_equal y).mp hy
  change (x.1, x.2.1) = (y.1, y.2.1) at hread
  have h0 : x.1 = y.1 := congrArg (fun p : Bool × Bool => p.1) hread
  have h1 : x.2.1 = y.2.1 := congrArg (fun p : Bool × Bool => p.2) hread
  have h2 : x.2.2 = y.2.2 := hxOuter.symm.trans (h0.trans hyOuter)
  exact Prod.ext h0 (Prod.ext h1 h2)

/-- Reading the two repeated outer coordinates does not determine an image
row. -/
theorem read02_not_injective_on_image :
    ¬ Set.InjOn read02 (Set.range step) := by
  intro hinj
  let x : Row := (false, false, false)
  let y : Row := (false, true, false)
  have hx : x ∈ Set.range step := ⟨(false, false, false), rfl⟩
  have hy : y ∈ Set.range step := ⟨(false, false, true), rfl⟩
  have hxy : x = y := hinj hx hy rfl
  exact absurd hxy (by decide)

/-- A total reverse repair preserving every target row would be a right
inverse of `step`; the out-of-image row `(0,0,1)` rules it out. -/
theorem no_total_reverse_repair :
    ¬ ∃ repair : Row → Row, ∀ target : Row, step (repair target) = target := by
  rintro ⟨repair, hrepair⟩
  let bad : Row := (false, false, true)
  have houter := outer_coordinates_equal (repair bad)
  rw [hrepair bad] at houter
  exact absurd houter (by decide)

end ObservableNormalForms.Rule90
