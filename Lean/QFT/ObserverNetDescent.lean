import QFT.FiniteCausalObserverNet

/-!
# Finite overlap descent for observer nets

This module gives the exact restriction-gluing interface for the finite causal
observer nets in `QFT/FiniteCausalObserverNet.lean`.  The API name
`FiniteCover` denotes only a nonempty finite family of declared subregions of
`W`: the abstract region preorder has no join or pointwise union operation, so
no joint-coverage law is present.  A local family is
compatible when its restrictions from any two cover regions agree on their
declared overlap.  `HasUniqueDescent` says exactly that every compatible
family is the restriction of a unique element of the local algebra of `W`.

Descent is not inferred from isotony, locality, or finiteness.  An
`ObserverNetDescent` explicitly marks the subregion families for which the
receipt is available and requires at least one marked family for every region.  The
negative-control theorem `no_descent_of_indistinguishable_global_sections`
shows that a cover cannot have unique descent when its restrictions fail to
separate two global sections.

The singleton identity family is proved to satisfy descent for every net.  It
is a nonvacuity/type-checking witness only; it supplies no nontrivial gluing
or coverage claim.  Nontrivial declared families must carry their own receipt;
a genuine cover additionally needs an independently defined coverage law.
-/

namespace OPH.QFT

open OPH.Tower

universe u

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι} (N : FiniteCausalObserverNet T)

namespace FiniteCausalObserverNet

/-- A nonempty finite family of declared subregions of one ambient region.
Despite the API name, no join/union coverage law is available in the abstract
region preorder. -/
structure FiniteCover (r : ι) (W : N.Region r) where
  regions : Finset (N.Region r)
  nonempty : regions.Nonempty
  subregion : ∀ U, U ∈ regions → N.regionLE r U W

namespace FiniteCover

variable {N}

/-- A family of local observables agreeing after restriction to every
pairwise declared overlap. -/
structure CompatibleFamily {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W) where
  localSection : ∀ U, U ∈ C.regions → N.localAlgebra r U
  overlap_compatible : ∀ U (hU : U ∈ C.regions)
      V (hV : V ∈ C.regions),
    N.restrict r (N.overlap_le_left r U V) (localSection U hU) =
      N.restrict r (N.overlap_le_right r U V) (localSection V hV)

/-- A global local-algebra element realizes a compatible family when all of
its declared restrictions are the given local sections. -/
def Globalizes {r : ι} {W : N.Region r} (C : N.FiniteCover r W)
    (F : C.CompatibleFamily) (X : N.localAlgebra r W) : Prop :=
  ∀ U (hU : U ∈ C.regions),
    N.restrict r (C.subregion U hU) X = F.localSection U hU

/-- Exact unique restriction-gluing proposition for one declared finite
subregion family. -/
def HasUniqueDescent {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W) : Prop :=
  ∀ F : C.CompatibleFamily,
    ∃! X : N.localAlgebra r W, C.Globalizes F X

/-- Restricting one global section to a cover produces a compatible family. -/
noncomputable def ofGlobal {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W) (X : N.localAlgebra r W) :
    C.CompatibleFamily where
  localSection U hU := N.restrict r (C.subregion U hU) X
  overlap_compatible U hU V hV := by
    rw [N.restrict_trans, N.restrict_trans]

@[simp]
theorem globalizes_ofGlobal {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W) (X : N.localAlgebra r W) :
    C.Globalizes (C.ofGlobal X) X := by
  intro U hU
  rfl

/-- Unique descent forces the family of cover restrictions to be jointly
injective on global sections. -/
theorem jointly_injective_of_unique_descent {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W) (hdescent : C.HasUniqueDescent)
    {X Y : N.localAlgebra r W}
    (hXY : ∀ U (hU : U ∈ C.regions),
      N.restrict r (C.subregion U hU) X =
        N.restrict r (C.subregion U hU) Y) :
    X = Y := by
  obtain ⟨Z, hZ, hunique⟩ := hdescent (C.ofGlobal X)
  have hX : C.Globalizes (C.ofGlobal X) X := C.globalizes_ofGlobal X
  have hY : C.Globalizes (C.ofGlobal X) Y := by
    intro U hU
    exact (hXY U hU).symm
  exact (hunique X hX).trans (hunique Y hY).symm

/-- **Negative control.** If two distinct global sections have identical
restrictions to every member of a cover, that cover cannot satisfy unique
descent.  Isotony and overlap compatibility alone therefore do not supply
the descent receipt. -/
theorem no_descent_of_indistinguishable_global_sections
    {r : ι} {W : N.Region r} (C : N.FiniteCover r W)
    {X Y : N.localAlgebra r W} (hne : X ≠ Y)
    (hXY : ∀ U (hU : U ∈ C.regions),
      N.restrict r (C.subregion U hU) X =
        N.restrict r (C.subregion U hU) Y) :
    ¬ C.HasUniqueDescent := by
  intro hdescent
  exact hne (C.jointly_injective_of_unique_descent hdescent hXY)

end FiniteCover

/-- A proof-carrying declaration of the subregion families on which unique
restriction gluing is available.  `cover_exists` rules out the vacuous packet
with no declared families at all; it does not assert joint coverage. -/
structure ObserverNetDescent where
  declaredCover : ∀ {r : ι} {W : N.Region r},
    N.FiniteCover r W → Prop
  cover_exists : ∀ r (W : N.Region r),
    ∃ C : N.FiniteCover r W, declaredCover C
  unique_descent : ∀ {r : ι} {W : N.Region r}
    (C : N.FiniteCover r W), declaredCover C → C.HasUniqueDescent

namespace ObserverNetDescent

variable {N}

/-- The unique glued section selected by a declared descent receipt. -/
noncomputable def glue (D : N.ObserverNetDescent)
    {r : ι} {W : N.Region r} (C : N.FiniteCover r W)
    (hC : D.declaredCover C) (F : C.CompatibleFamily) :
    N.localAlgebra r W :=
  Classical.choose (D.unique_descent C hC F)

/-- The selected glue restricts to every member of the compatible family. -/
theorem glue_restrict (D : N.ObserverNetDescent)
    {r : ι} {W : N.Region r} (C : N.FiniteCover r W)
    (hC : D.declaredCover C) (F : C.CompatibleFamily)
    (U : N.Region r) (hU : U ∈ C.regions) :
    N.restrict r (C.subregion U hU) (D.glue C hC F) =
      F.localSection U hU :=
  (Classical.choose_spec (D.unique_descent C hC F)).1 U hU

/-- Any other global section with the same restrictions is the selected
glue. -/
theorem glue_unique (D : N.ObserverNetDescent)
    {r : ι} {W : N.Region r} (C : N.FiniteCover r W)
    (hC : D.declaredCover C) (F : C.CompatibleFamily)
    (X : N.localAlgebra r W) (hX : C.Globalizes F X) :
    X = D.glue C hC F :=
  (Classical.choose_spec (D.unique_descent C hC F)).2 X hX

end ObserverNetDescent

/-! ## A universal identity-family receipt -/

/-- The singleton family containing only the ambient region. -/
noncomputable def singletonCover (r : ι) (W : N.Region r) :
    N.FiniteCover r W where
  regions := {W}
  nonempty := ⟨W, by simp⟩
  subregion U hU := by
    simp only [Finset.mem_singleton] at hU
    subst U
    exact N.regionLE_refl r W

/-- Every singleton identity family has exact unique restriction gluing. -/
theorem singletonCover_hasUniqueDescent (r : ι) (W : N.Region r) :
    (N.singletonCover r W).HasUniqueDescent := by
  intro F
  let hW : W ∈ (N.singletonCover r W).regions := by
    simp [singletonCover]
  refine ⟨F.localSection W hW, ?_, ?_⟩
  · intro U hU
    have hUW : U = W := by
      simpa [singletonCover] using hU
    subst U
    simpa using N.restrict_refl r W (F.localSection W hW)
  · intro X hX
    have hx := hX W hW
    simpa using (N.restrict_refl r W X).symm.trans hx

/-- The universal descent packet declaring precisely the singleton identity
families.  This proves the type is inhabited without asserting nontrivial
overlap gluing or a separate coverage law. -/
noncomputable def singletonDescent : N.ObserverNetDescent where
  declaredCover := fun {r} {W} C => C = N.singletonCover r W
  cover_exists := by
    intro r W
    exact ⟨N.singletonCover r W, rfl⟩
  unique_descent := by
    intro r W C hC
    subst C
    exact N.singletonCover_hasUniqueDescent r W

/-! ## Nontrivial-family receipt for the finite partition model -/

section PartitionPublicModel

variable {n k : ℕ}

/-- Every nonempty subregion family in the degenerate partition-public model has unique
descent.  All regional algebras and all restriction maps in that model are
definitionally the same, so pairwise overlap compatibility forces all local
sections to be equal. -/
noncomputable def partitionPublicAllCoversDescent
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    (partitionPublicCausalNet part rho).ObserverNetDescent where
  declaredCover := fun _ => True
  cover_exists := by
    intro r W
    exact ⟨(partitionPublicCausalNet part rho).singletonCover r W, trivial⟩
  unique_descent := by
    intro r W C hC F
    obtain ⟨U, hU⟩ := C.nonempty
    refine ⟨F.localSection U hU, ?_, ?_⟩
    · intro V hV
      have hcompat := F.overlap_compatible U hU V hV
      simpa [partitionPublicCausalNet] using hcompat
    · intro X hX
      have hx := hX U hU
      simpa [partitionPublicCausalNet] using hx

/-- A two-member subregion family of the full two-label region in the finite
partition-public model. -/
noncomputable def partitionPublicTwoRegionCover
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    (partitionPublicCausalNet part rho).FiniteCover ()
      ({0, 1} : Finset (Fin 2)) where
  regions :=
    ({({0} : Finset (Fin 2)), ({1} : Finset (Fin 2))} :
      Finset (Finset (Fin 2)))
  nonempty := by
    refine ⟨({0} : Finset (Fin 2)), ?_⟩
    change ({0} : Finset (Fin 2)) ∈
      ({({0} : Finset (Fin 2)), ({1} : Finset (Fin 2))} :
        Finset (Finset (Fin 2)))
    simp
  subregion (U : Finset (Fin 2)) _ := by
    change U ⊆ ({0, 1} : Finset (Fin 2))
    intro x hx
    fin_cases x <;> simp

/-- The explicit two-member family has exact unique restriction gluing.  This is stronger
than the universal singleton receipt but remains a degenerate common-algebra
model, not a formal coverage theorem or spatial tensor-factor reconstruction. -/
theorem partitionPublicTwoRegionCover_hasUniqueDescent
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    (partitionPublicTwoRegionCover part rho).HasUniqueDescent :=
  (partitionPublicAllCoversDescent part rho).unique_descent
    (partitionPublicTwoRegionCover part rho) trivial

end PartitionPublicModel

end FiniteCausalObserverNet

-- Axiom audit: only Lean/Mathlib extensionality and choice principles occur.
#print axioms OPH.QFT.FiniteCausalObserverNet.FiniteCover.jointly_injective_of_unique_descent
#print axioms OPH.QFT.FiniteCausalObserverNet.FiniteCover.no_descent_of_indistinguishable_global_sections
#print axioms OPH.QFT.FiniteCausalObserverNet.ObserverNetDescent.glue_restrict
#print axioms OPH.QFT.FiniteCausalObserverNet.singletonCover_hasUniqueDescent
#print axioms OPH.QFT.FiniteCausalObserverNet.singletonDescent
#print axioms OPH.QFT.FiniteCausalObserverNet.partitionPublicTwoRegionCover_hasUniqueDescent

end OPH.QFT
