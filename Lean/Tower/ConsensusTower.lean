import EventAlgebra.PublicRecordAlgebra
import Time.TimeOrderLedger

/-!
# Timeless finite observer-consensus towers

This module supplies the structural root requested by completion-plan issue
`#698`.  A `ConsensusTower` is a covariant refinement system of finite matrix
algebras.  Each regulator carries finite observer and record types,
observer-indexed record orders, a commutative public star-subalgebra, selected
density states, and linear generators.  Refinement maps and every compatibility
law are fields of the structure; none is inferred from a limit or a physical
interpretation.

The word *timeless* is literal here.  The tower carries the A1
`ObserverRecordOrder` type but no worldline, clock readout, proper time,
Lorentzian event, causal region, continuum limit, or global time function.
Those are downstream decorations with named realization maps.

`ConsensusTower.constantConsensusTower` adapts the existing finite event
algebra partition machinery to a constant one-regulator tower.  It reuses the
partition's public algebra and active record labels without changing either
definition.  The caller supplies one certified density state; the
constant tower uses the zero generator.
-/

namespace OPH.Tower

open EventAlgebra
open OPH.TimeOrderLedger

universe u

/-- A timeless refinement system carrying finite observer-facing record and
matrix-algebra data.

The maps point from a coarser regulator `r` to a refinement `s` when `r ≤ s`.
The selected states and generators are structural data, not a claim that a
physical source has produced them.  State compatibility is contravariant:
the refined density has the same expectation on every mapped coarse
observable.  It is not the generally false claim that a star-algebra
homomorphism pushes a normalized density matrix to another one. -/
structure ConsensusTower (ι : Type u) [Preorder ι] where
  /-- The regulator system is inhabited. -/
  indexNonempty : Nonempty ι
  /-- Any two regulators have a declared common refinement. -/
  commonRefinement : ∀ r s : ι, ∃ t, r ≤ t ∧ s ≤ t
  /-- Matrix size of the private finite algebra at each regulator. -/
  dim : ι → ℕ
  /-- Finite observer labels at each regulator. -/
  Observer : ι → Type u
  observerFinite : ∀ r, Fintype (Observer r)
  /-- Finite public record labels at each regulator. -/
  Record : ι → Type u
  recordFinite : ∀ r, Fintype (Record r)
  /-- Each observer has an explicit strict order on the common record type. -/
  recordOrder : ∀ r, Observer r → ObserverRecordOrder (Record r)
  /-- Observer-indexed commutative public records inside the private matrix
  algebra. -/
  publicAlgebra : ∀ r, Observer r →
    StarSubalgebra ℂ (Matrix (Fin (dim r)) (Fin (dim r)) ℂ)
  public_mul_comm : ∀ r o (X Y : publicAlgebra r o), X * Y = Y * X
  /-- A record label is represented by a named element of the corresponding
  public algebra.  Injectivity or geometric separation is not assumed. -/
  recordElement : ∀ r o, Record r → publicAlgebra r o
  /-- One selected finite density state per observer. -/
  state : ∀ r, Observer r → Matrix (Fin (dim r)) (Fin (dim r)) ℂ
  state_isState : ∀ r o, IsState (state r o)
  /-- One selected algebraic generator per observer.  No semigroup, locality,
  Lindblad, modular, or physical-time interpretation is included. -/
  generator : ∀ r, Observer r →
    Matrix (Fin (dim r)) (Fin (dim r)) ℂ →ₗ[ℂ]
      Matrix (Fin (dim r)) (Fin (dim r)) ℂ

  /- Refinement data. -/
  observerRefine : ∀ {r s}, r ≤ s → Observer r → Observer s
  recordRefine : ∀ {r s}, r ≤ s → Record r → Record s
  algebraRefine : ∀ {r s}, r ≤ s →
    Matrix (Fin (dim r)) (Fin (dim r)) ℂ →⋆ₐ[ℂ]
      Matrix (Fin (dim s)) (Fin (dim s)) ℂ

  /- Identity and composition receipts: the three refinement families are
  functorial rather than merely a collection of unrelated maps. -/
  observer_refine_refl : ∀ r o,
    observerRefine (show r ≤ r from le_rfl) o = o
  observer_refine_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) o,
    observerRefine hst (observerRefine hrs o) =
      observerRefine (hrs.trans hst) o
  record_refine_refl : ∀ r x,
    recordRefine (show r ≤ r from le_rfl) x = x
  record_refine_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) x,
    recordRefine hst (recordRefine hrs x) = recordRefine (hrs.trans hst) x
  algebra_refine_refl : ∀ r X,
    algebraRefine (show r ≤ r from le_rfl) X = X
  algebra_refine_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) X,
    algebraRefine hst (algebraRefine hrs X) =
      algebraRefine (hrs.trans hst) X

  /- Cross-layer compatibility receipts. -/
  record_order_natural : ∀ {r s} (hrs : r ≤ s) o {x y},
    (recordOrder r o).precedes x y →
      (recordOrder s (observerRefine hrs o)).precedes
        (recordRefine hrs x) (recordRefine hrs y)
  public_natural : ∀ {r s} (hrs : r ≤ s) o X,
    X ∈ publicAlgebra r o →
      algebraRefine hrs X ∈ publicAlgebra s (observerRefine hrs o)
  record_element_natural : ∀ {r s} (hrs : r ≤ s) o x,
    algebraRefine hrs (recordElement r o x).1 =
      (recordElement s (observerRefine hrs o) (recordRefine hrs x)).1
  state_natural : ∀ {r s} (hrs : r ≤ s) o X,
    (state s (observerRefine hrs o) * algebraRefine hrs X).trace =
      (state r o * X).trace
  generator_natural : ∀ {r s} (hrs : r ≤ s) o X,
    algebraRefine hrs (generator r o X) =
      generator s (observerRefine hrs o) (algebraRefine hrs X)

attribute [instance] ConsensusTower.observerFinite ConsensusTower.recordFinite

namespace ConsensusTower

variable {ι : Type u} [Preorder ι] (T : ConsensusTower ι)

/-- The private finite matrix algebra at a regulator. -/
abbrev PrivateAlgebra (r : ι) :=
  Matrix (Fin (T.dim r)) (Fin (T.dim r)) ℂ

/-- Public elements remain public under any declared refinement. -/
theorem public_mem_refine {r s : ι} (hrs : r ≤ s) (o : T.Observer r)
    (X : T.publicAlgebra r o) :
    T.algebraRefine hrs X.1 ∈
      T.publicAlgebra s (T.observerRefine hrs o) :=
  T.public_natural hrs o X.1 X.2

/-- Refinement commutes with the public representative of a record label. -/
theorem refine_recordElement {r s : ι} (hrs : r ≤ s)
    (o : T.Observer r) (x : T.Record r) :
    T.algebraRefine hrs (T.recordElement r o x).1 =
      (T.recordElement s (T.observerRefine hrs o)
        (T.recordRefine hrs x)).1 :=
  T.record_element_natural hrs o x

/-- Observer record precedence is preserved by every declared refinement. -/
theorem refine_precedes {r s : ι} (hrs : r ≤ s)
    (o : T.Observer r) {x y : T.Record r}
    (hxy : (T.recordOrder r o).precedes x y) :
    (T.recordOrder s (T.observerRefine hrs o)).precedes
      (T.recordRefine hrs x) (T.recordRefine hrs y) :=
  T.record_order_natural hrs o hxy

/-- Generator naturality is an intertwining equation on private algebras. -/
theorem refine_generator {r s : ι} (hrs : r ≤ s)
    (o : T.Observer r) (X : T.PrivateAlgebra r) :
    T.algebraRefine hrs (T.generator r o X) =
      T.generator s (T.observerRefine hrs o) (T.algebraRefine hrs X) :=
  T.generator_natural hrs o X

/-- Refined states restrict to the selected coarse state on mapped coarse
observables. This is dual to the covariant star-algebra refinement map;
injectivity is not part of the interface. -/
theorem refine_state_pairing {r s : ι} (hrs : r ≤ s)
    (o : T.Observer r) (X : T.PrivateAlgebra r) :
    (T.state s (T.observerRefine hrs o) * T.algebraRefine hrs X).trace =
      (T.state r o * X).trace :=
  T.state_natural hrs o X

variable {n k : ℕ}

/-- Adapt an existing projective partition and certified density matrix to a
constant, one-regulator consensus tower.

The private algebra is `M_n(ℂ)`, the sole observer's public algebra is exactly
`part.publicSubalgebra`, and the record type is exactly `part.ActiveIndex`.
The discrete order and zero generator prevent this adaptor from smuggling in
history, dynamics, or a physical clock. -/
noncomputable def constantConsensusTower
    (part : EventAlgebra.ProjectivePartition n k)
    (ρ : EventAlgebra.StateMatrix n) : ConsensusTower Unit where
  indexNonempty := ⟨()⟩
  commonRefinement := by
    intro r s
    exact ⟨(), by trivial, by trivial⟩
  dim := fun _ => n
  Observer := fun _ => Unit
  observerFinite := fun _ => inferInstance
  Record := fun _ => part.ActiveIndex
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ _ => ObserverRecordOrder.discrete part.ActiveIndex
  publicAlgebra := fun _ _ => part.publicSubalgebra
  public_mul_comm := fun _ _ X Y => part.publicSubalgebra_mul_comm X Y
  recordElement := fun _ _ i =>
    ⟨part.proj i.1, part.proj_mem_span i.1⟩
  state := fun _ _ => ρ.1
  state_isState := fun _ _ => ρ.2
  generator := fun _ _ => 0
  observerRefine := fun _ o => o
  recordRefine := fun _ x => x
  algebraRefine := fun _ => StarAlgHom.id ℂ _
  observer_refine_refl := by intros; rfl
  observer_refine_trans := by intros; rfl
  record_refine_refl := by intros; rfl
  record_refine_trans := by intros; rfl
  algebra_refine_refl := by intros; rfl
  algebra_refine_trans := by intros; rfl
  record_order_natural := by simp [ObserverRecordOrder.discrete]
  public_natural := by
    intro r s hrs o X hX
    exact hX
  record_element_natural := by intros; rfl
  state_natural := by intros; rfl
  generator_natural := by intros; simp

@[simp]
theorem constantConsensusTower_public
    (part : EventAlgebra.ProjectivePartition n k)
    (ρ : EventAlgebra.StateMatrix n) (o : Unit) :
    (constantConsensusTower part ρ).publicAlgebra () o =
      part.publicSubalgebra :=
  rfl

@[simp]
theorem constantConsensusTower_recordElement
    (part : EventAlgebra.ProjectivePartition n k)
    (ρ : EventAlgebra.StateMatrix n) (o : Unit) (i : part.ActiveIndex) :
    ((constantConsensusTower part ρ).recordElement () o i).1 =
      part.proj i.1 :=
  rfl

end ConsensusTower

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.Tower.ConsensusTower.public_mem_refine
#print axioms OPH.Tower.ConsensusTower.refine_recordElement
#print axioms OPH.Tower.ConsensusTower.refine_precedes
#print axioms OPH.Tower.ConsensusTower.refine_generator
#print axioms OPH.Tower.ConsensusTower.refine_state_pairing
#print axioms OPH.Tower.ConsensusTower.constantConsensusTower_public
#print axioms OPH.Tower.ConsensusTower.constantConsensusTower_recordElement

end OPH.Tower
