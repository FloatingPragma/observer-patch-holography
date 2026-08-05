import Tower.ConsensusTower

/-!
# Finite public-world quotients

This module implements the quotient half of completion-plan issue `#699`.
At one finite regulator, a raw configuration is publicly equivalent to another
exactly when their complete observer-indexed public-record signatures agree.
The resulting `PublicWorld` is therefore a literal kernel quotient, not a
metaphor and not an assertion that a physical world has been constructed.

`atRegulator` is the typed handoff from `ConsensusTower`: it fixes the public
signature to `∀ o : T.Observer r, T.Record r` while leaving the raw
configuration carrier and its readback map explicit.  In particular, the A3
tower does not manufacture a source-selected configuration space or readback.

The hidden-bit adaptor is a nonvacuity receipt.  It adds one private bit to a
public signature and proves that the two raw configurations are distinct but
become equal in the public quotient.  It is a finite mathematical example,
not a model of physical gauge redundancy.
-/

namespace OPH.Tower

universe u v

/-- A finite raw configuration carrier equipped with a public readback.

Only the raw carrier is required to be finite.  The public signature can be
any type; the quotient is finite because it is a quotient of `Config`. -/
structure PublicWorldPresentation (Signature : Type v) where
  Config : Type u
  configFintype : Fintype Config
  configNonempty : Nonempty Config
  readback : Config → Signature

attribute [instance] PublicWorldPresentation.configFintype
attribute [instance] PublicWorldPresentation.configNonempty

namespace PublicWorldPresentation

variable {Signature : Type v} (P : PublicWorldPresentation Signature)

/-- Public-record equivalence is equality of the declared readback. -/
def publicSetoid : Setoid P.Config := Setoid.ker P.readback

/-- The finite public world represented by the declared public records. -/
abbrev PublicWorld := Quotient P.publicSetoid

/-- Send a raw configuration to its public equivalence class. -/
def toPublicWorld (x : P.Config) : P.PublicWorld :=
  Quotient.mk P.publicSetoid x

/-- Two raw configurations define the same public world exactly when their
declared public signatures agree. -/
theorem toPublicWorld_eq_iff (x y : P.Config) :
    P.toPublicWorld x = P.toPublicWorld y ↔ P.readback x = P.readback y := by
  constructor
  · intro h
    exact Quotient.exact h
  · intro h
    exact Quotient.sound h

/-- The public signature descends to the kernel quotient. -/
def publicSignature : P.PublicWorld → Signature :=
  Quotient.lift P.readback (fun _ _ h => h)

@[simp]
theorem publicSignature_mk (x : P.Config) :
    P.publicSignature (P.toPublicWorld x) = P.readback x :=
  rfl

/-- The descended public signature is a complete invariant of the quotient.
This is the exact fixed-regulator meaning of quotienting by public records. -/
theorem publicSignature_injective : Function.Injective P.publicSignature := by
  intro q₁ q₂ h
  refine Quotient.inductionOn₂ q₁ q₂ ?_ h
  intro x y hxy
  exact (P.toPublicWorld_eq_iff x y).2 hxy

/-- The public quotient is exactly the type of signatures actually realized
by raw configurations; unrealized values of `Signature` are not added. -/
def RealizedSignature := Set.range P.readback

/-- Canonical equivalence between the kernel quotient and the range of the
declared public readback. -/
noncomputable def publicWorldEquivRealizedSignature :
    P.PublicWorld ≃ P.RealizedSignature where
  toFun q := ⟨P.publicSignature q, by
    refine Quotient.inductionOn q ?_
    intro x
    exact ⟨x, rfl⟩⟩
  invFun s := P.toPublicWorld (Classical.choose s.property)
  left_inv q := by
    apply P.publicSignature_injective
    exact Classical.choose_spec (show ∃ x, P.readback x = P.publicSignature q from by
      refine Quotient.inductionOn q ?_
      intro x
      exact ⟨x, rfl⟩)
  right_inv s := by
    apply Subtype.ext
    exact Classical.choose_spec s.property

/-- Every declared finite presentation has an inhabited public quotient. -/
theorem publicWorld_nonempty : Nonempty P.PublicWorld := by
  obtain ⟨x⟩ := P.configNonempty
  exact ⟨P.toPublicWorld x⟩

/-- The public quotient is finite because every class has a raw
representative in the finite configuration carrier. -/
noncomputable instance publicWorldFinite : Finite P.PublicWorld :=
  Finite.of_surjective P.toPublicWorld (by
    intro q
    obtain ⟨x, hx⟩ := Quotient.exists_rep q
    exact ⟨x, hx⟩)

end PublicWorldPresentation

namespace ConsensusTower

variable {ι : Type u} [Preorder ι] (T : ConsensusTower ι)

/-- Complete observer-indexed public-record signature at one regulator. -/
abbrev PublicSignature (r : ι) := T.Observer r → T.Record r

/-- Attach an explicit finite raw configuration/readback presentation to one
regulator of an A3 `ConsensusTower`.  The tower supplies the types of public
records; the caller supplies the configurations and the readback map. -/
@[reducible] def atRegulator (r : ι) (Config : Type v)
    [Fintype Config] [Nonempty Config]
    (readback : Config → T.PublicSignature r) :
    PublicWorldPresentation (T.PublicSignature r) where
  Config := Config
  configFintype := inferInstance
  configNonempty := inferInstance
  readback := readback

/-- A finite synthetic presentation with one hidden private bit.  The seed is
used only to certify that the observer-indexed signature type is inhabited. -/
@[reducible] noncomputable def hiddenBitAtRegulator
    (r : ι) (seed : T.PublicSignature r) :
    PublicWorldPresentation (T.PublicSignature r) where
  Config := T.PublicSignature r × Bool
  configFintype := by classical exact inferInstance
  configNonempty := ⟨(seed, false)⟩
  readback := Prod.fst

/-- The hidden-bit adaptor has genuine raw redundancy: changing only the
private bit changes the raw configuration but not its public-world class. -/
theorem hiddenBit_distinct_but_publicly_equal (r : ι)
    (seed : T.PublicSignature r) :
    (seed, false) ≠ (seed, true) ∧
      (T.hiddenBitAtRegulator r seed).toPublicWorld (seed, false) =
        (T.hiddenBitAtRegulator r seed).toPublicWorld (seed, true) := by
  constructor
  · intro h
    exact Bool.false_ne_true (congrArg Prod.snd h)
  · exact ((T.hiddenBitAtRegulator r seed).toPublicWorld_eq_iff _ _).2 rfl

end ConsensusTower

-- Axiom audit: only Lean/Mathlib quotient principles are used.
#print axioms OPH.Tower.PublicWorldPresentation.toPublicWorld_eq_iff
#print axioms OPH.Tower.PublicWorldPresentation.publicSignature_injective
#print axioms OPH.Tower.PublicWorldPresentation.publicWorldEquivRealizedSignature
#print axioms OPH.Tower.PublicWorldPresentation.publicWorld_nonempty
#print axioms OPH.Tower.ConsensusTower.hiddenBit_distinct_but_publicly_equal

end OPH.Tower
