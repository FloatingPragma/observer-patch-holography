import Mathlib.Algebra.FreeAbelianGroup.Finsupp
import Mathlib.GroupTheory.Subgroup.Centralizer

namespace ObserverPatchHolography.A2EndpointCommutator

/-!
# A2 endpoint descent and the six-axis abelian factorization

This file isolates the positive conditional theorem used by the directed
twelve-port source lane.  An accepted primitive step on a presentation `X`
descends through an interpretation `q : X → Q` when it has a meaning-side map
making the naturality square commute.  If two ordered presentation paths have
the same interpreted endpoint for every starting state, surjectivity of `q`
then forces the two meaning-side maps to commute.

For six positive port-axis generators, the fifteen strict-index commutators
are enough to make every generator pair commute.  Adding the six reverse
orientations as formal inverses makes every signed word factor through the
free abelian group on `Fin 6`, equivalently its integer exponent vector.

The hypotheses are ordinary theorem inputs.  This file does not assert that
the simulator source emits the twelve accepted step maps or the fifteen
endpoint-diamond tables.  It also does not identify the meaning object `Q`
with spatial support or prove that the resulting free-abelian action is
faithful.
-/

universe u v

/-- A presentation step `T` descends through `q` to the meaning-side map
`τ`.  In the OPH application this is the exact naturality square supplied by
A2 only after A1 has typed `T` as an accepted data-access map. -/
def Descends {X : Type u} {Q : Type v}
    (q : X → Q) (T : X → X) (τ : Q → Q) : Prop :=
  ∀ x, q (T x) = τ (q x)

/-- **Surjective endpoint-descent lemma.**  All-state equality of the two
interpreted `AB` and `BA` endpoints forces the descended meaning-side maps to
commute.  Equality at only one starting state would not suffice. -/
theorem commute_of_surjective_descents
    {X : Type u} {Q : Type v}
    {q : X → Q} {A B : X → X} {a b : Q → Q}
    (hq : Function.Surjective q)
    (hA : Descends q A a) (hB : Descends q B b)
    (hend : ∀ x, q (A (B x)) = q (B (A x))) :
    Function.Commute a b := by
  intro y
  obtain ⟨x, rfl⟩ := hq y
  calc
    a (b (q x)) = a (q (B x)) := congrArg a (hB x).symm
    _ = q (A (B x)) := (hA (B x)).symm
    _ = q (B (A x)) := hend x
    _ = b (q (A x)) := hB (A x)
    _ = b (a (q x)) := congrArg b (hA x)

/-- The six positive port axes. -/
abbrev Axis := Fin 6

/-- A signed port is one of the six axes together with its orientation.
`reversed = false` is the positive representative and `true` its inverse. -/
structure SignedPort where
  axis : Axis
  reversed : Bool
deriving DecidableEq, Repr

/-- The positive orientation of an axis. -/
def forward (i : Axis) : SignedPort :=
  ⟨i, false⟩

/-- The antipodal orientation of an axis. -/
def reverse (i : Axis) : SignedPort :=
  ⟨i, true⟩

/-- Reverse the orientation while retaining the axis. -/
def antipode (p : SignedPort) : SignedPort :=
  ⟨p.axis, !p.reversed⟩

@[simp]
theorem antipode_forward (i : Axis) : antipode (forward i) = reverse i := rfl

@[simp]
theorem antipode_reverse (i : Axis) : antipode (reverse i) = forward i := rfl

@[simp]
theorem antipode_involutive (p : SignedPort) : antipode (antipode p) = p := by
  cases p with
  | mk i reversed => cases reversed <;> rfl

/-- The signed basis element in the universal six-axis abelian group. -/
def signedBasis (p : SignedPort) : FreeAbelianGroup Axis :=
  if p.reversed then -FreeAbelianGroup.of p.axis else FreeAbelianGroup.of p.axis

/-- The exponent vector of a signed word, represented intrinsically as an
element of the free abelian group on the six axes. -/
def exponentVector : List SignedPort → FreeAbelianGroup Axis
  | [] => 0
  | p :: word => signedBasis p + exponentVector word

/-- The same exponent vector in the concrete finitely supported integer
coordinate model.  Since `Axis = Fin 6`, this is the explicit `ℤ⁶` normal
form. -/
noncomputable def integerExponentVector (word : List SignedPort) : Axis →₀ ℤ :=
  FreeAbelianGroup.toFinsupp (exponentVector word)

@[simp]
theorem signedBasis_forward (i : Axis) :
    signedBasis (forward i) = FreeAbelianGroup.of i := rfl

@[simp]
theorem signedBasis_reverse (i : Axis) :
    signedBasis (reverse i) = -FreeAbelianGroup.of i := rfl

@[simp]
theorem exponentVector_nil : exponentVector [] = 0 := rfl

@[simp]
theorem exponentVector_cons (p : SignedPort) (word : List SignedPort) :
    exponentVector (p :: word) = signedBasis p + exponentVector word := rfl

/-- Each of the six inverse pairs cancels in the exponent vector. -/
@[simp]
theorem exponentVector_forward_reverse (i : Axis) :
    exponentVector [forward i, reverse i] = 0 := by
  simp [exponentVector]

/-- The fifteen strict-index conditions, one for each unordered pair of six
positive generators. -/
def FifteenCommute {G : Type*} [Group G] (g : Axis → G) : Prop :=
  ∀ i j, i < j → Commute (g i) (g j)

/-- The fifteen strict-index conditions imply commutation for every ordered
pair, including the diagonal and reversed index order. -/
theorem commute_all_of_fifteen
    {G : Type*} [Group G] {g : Axis → G}
    (h : FifteenCommute g) (i j : Axis) : Commute (g i) (g j) := by
  rcases lt_trichotomy i j with hij | hij | hij
  · exact h i j hij
  · subst j
    exact Commute.refl (g i)
  · exact (h j i hij).symm

/-- The subgroup generated by the six positive generators. -/
def generatedSubgroup {G : Type*} [Group G] (g : Axis → G) : Subgroup G :=
  Subgroup.closure (Set.range g)

/-- A positive generator as an element of the generated subgroup. -/
def generatorInClosure {G : Type*} [Group G] (g : Axis → G) (i : Axis) :
    generatedSubgroup g :=
  ⟨g i, Subgroup.subset_closure (Set.mem_range_self i)⟩

/-- Interpret a signed port in the subgroup generated by the six positive
generators. -/
def signedGenerator {G : Type*} [Group G] (g : Axis → G)
    (p : SignedPort) : generatedSubgroup g :=
  if p.reversed then (generatorInClosure g p.axis)⁻¹
  else generatorInClosure g p.axis

/-- Evaluate a signed word in its original order. -/
def evalWord {G : Type*} [Group G] (g : Axis → G) :
    List SignedPort → generatedSubgroup g
  | [] => 1
  | p :: word => signedGenerator g p * evalWord g word

@[simp]
theorem signedGenerator_forward {G : Type*} [Group G]
    (g : Axis → G) (i : Axis) :
    signedGenerator g (forward i) = generatorInClosure g i := rfl

@[simp]
theorem signedGenerator_reverse {G : Type*} [Group G]
    (g : Axis → G) (i : Axis) :
    signedGenerator g (reverse i) = (generatorInClosure g i)⁻¹ := rfl

@[simp]
theorem evalWord_forward_reverse {G : Type*} [Group G]
    (g : Axis → G) (i : Axis) :
    evalWord g [forward i, reverse i] = 1 := by
  simp [evalWord]

/-- The range of a fifteen-commuting generator family commutes pairwise. -/
private theorem range_commutes
    {G : Type*} [Group G] {g : Axis → G}
    (h : FifteenCommute g) :
    ∀ x ∈ Set.range g, ∀ y ∈ Set.range g, x * y = y * x := by
  rintro _ ⟨i, rfl⟩ _ ⟨j, rfl⟩
  exact (commute_all_of_fifteen h i j).eq

/-- **Six-axis word factorization.**  Six inverse pairs and the fifteen
positive-generator commutators produce a homomorphism from the free abelian
group.  Every signed word evaluates through its exponent vector. -/
theorem exists_freeAbelian_word_factorization
    {G : Type*} [Group G] (g : Axis → G)
    (h : FifteenCommute g) :
    ∃ φ : FreeAbelianGroup Axis →+
        Additive (generatedSubgroup g),
      (∀ i, φ (FreeAbelianGroup.of i) =
        Additive.ofMul (generatorInClosure g i)) ∧
      ∀ word,
        Additive.toMul (φ (exponentVector word)) = evalWord g word := by
  letI : CommGroup (generatedSubgroup g) :=
    Subgroup.closureCommGroupOfComm (range_commutes h)
  let φ : FreeAbelianGroup Axis →+ Additive (generatedSubgroup g) :=
    FreeAbelianGroup.lift fun i => Additive.ofMul (generatorInClosure g i)
  refine ⟨φ, ?_, ?_⟩
  · intro i
    dsimp [φ]
    exact FreeAbelianGroup.lift_apply_of _ _
  · intro word
    induction word with
    | nil => simp [exponentVector, evalWord]
    | cons p word ih =>
        rcases p with ⟨i, reversed⟩
        have hφi :
            φ (FreeAbelianGroup.of i) =
              Additive.ofMul (generatorInClosure g i) := by
          dsimp [φ]
          exact FreeAbelianGroup.lift_apply_of _ _
        cases reversed with
        | false =>
            change Additive.toMul
                (φ (FreeAbelianGroup.of i + exponentVector word)) =
              generatorInClosure g i * evalWord g word
            rw [map_add, toMul_add, hφi, toMul_ofMul, ih]
        | true =>
            change Additive.toMul
                (φ (-FreeAbelianGroup.of i + exponentVector word)) =
              (generatorInClosure g i)⁻¹ * evalWord g word
            rw [map_add, map_neg, toMul_add, toMul_neg,
              hφi, toMul_ofMul, ih]
            rfl

/-- Under the fifteen commutators, equal exponent vectors imply equal word
actions.  This is the exact sense in which the ordered words factor through
the universal abelian port module. -/
theorem evalWord_eq_of_exponentVector_eq
    {G : Type*} [Group G] (g : Axis → G)
    (h : FifteenCommute g) {left right : List SignedPort}
    (hexponent : exponentVector left = exponentVector right) :
    evalWord g left = evalWord g right := by
  obtain ⟨φ, _hgenerator, hword⟩ :=
    exists_freeAbelian_word_factorization g h
  calc
    evalWord g left = Additive.toMul (φ (exponentVector left)) :=
      (hword left).symm
    _ = Additive.toMul (φ (exponentVector right)) := by rw [hexponent]
    _ = evalWord g right := hword right

/-- Equality of the explicit integer exponent vectors is sufficient for
equality of the signed-word actions. -/
theorem evalWord_eq_of_integerExponentVector_eq
    {G : Type*} [Group G] (g : Axis → G)
    (h : FifteenCommute g) {left right : List SignedPort}
    (hexponent : integerExponentVector left = integerExponentVector right) :
    evalWord g left = evalWord g right := by
  apply evalWord_eq_of_exponentVector_eq g h
  exact (FreeAbelianGroup.equivFinsupp Axis).injective hexponent

/-- The free-abelian factorization is unique once its values on the six
positive generators are fixed. -/
theorem freeAbelian_factorization_unique
    {G : Type*} [Group G] (g : Axis → G)
    (h : FifteenCommute g)
    (φ ψ : FreeAbelianGroup Axis →+ Additive (generatedSubgroup g))
    (hφ : ∀ i, φ (FreeAbelianGroup.of i) =
      Additive.ofMul (generatorInClosure g i))
    (hψ : ∀ i, ψ (FreeAbelianGroup.of i) =
      Additive.ofMul (generatorInClosure g i)) :
    φ = ψ := by
  letI : CommGroup (generatedSubgroup g) :=
    Subgroup.closureCommGroupOfComm (range_commutes h)
  exact FreeAbelianGroup.lift_ext φ ψ fun i => (hφ i).trans (hψ i).symm

/-- Fifteen all-state interpreted endpoint diamonds turn six descended
meaning-side permutations into a `FifteenCommute` family. -/
theorem quotient_fifteen_commute
    {X : Type u} {Q : Type v}
    {q : X → Q} (T : Axis → X → X) (τ : Axis → Equiv.Perm Q)
    (hq : Function.Surjective q)
    (hdesc : ∀ i, Descends q (T i) (τ i))
    (hdiamonds : ∀ i j, i < j → ∀ x,
      q (T i (T j x)) = q (T j (T i x))) :
    FifteenCommute τ := by
  intro i j hij
  have hfun : Function.Commute (τ i : Q → Q) (τ j : Q → Q) :=
    commute_of_surjective_descents hq (hdesc i) (hdesc j)
      (hdiamonds i j hij)
  show τ i * τ j = τ j * τ i
  ext y
  simpa only [Equiv.Perm.mul_apply] using hfun y

/-- Combined endpoint-to-word theorem.  A source packet containing the six
accepted positive step maps and all fifteen endpoint diamonds supplies the
free-abelian word factorization on the descended quotient permutations. -/
theorem exists_quotient_word_factorization
    {X : Type u} {Q : Type v}
    {q : X → Q} (T : Axis → X → X) (τ : Axis → Equiv.Perm Q)
    (hq : Function.Surjective q)
    (hdesc : ∀ i, Descends q (T i) (τ i))
    (hdiamonds : ∀ i j, i < j → ∀ x,
      q (T i (T j x)) = q (T j (T i x))) :
    ∃ φ : FreeAbelianGroup Axis →+
        Additive (generatedSubgroup τ),
      (∀ i, φ (FreeAbelianGroup.of i) =
        Additive.ofMul (generatorInClosure τ i)) ∧
      ∀ word,
        Additive.toMul (φ (exponentVector word)) = evalWord τ word :=
  exists_freeAbelian_word_factorization τ
    (quotient_fifteen_commute T τ hq hdesc hdiamonds)

end ObserverPatchHolography.A2EndpointCommutator

/- Axiom audit: ordinary function, group, quotient, and free-abelian-group
reasoning only. -/

#print axioms ObserverPatchHolography.A2EndpointCommutator.commute_of_surjective_descents
#print axioms ObserverPatchHolography.A2EndpointCommutator.exists_freeAbelian_word_factorization
#print axioms ObserverPatchHolography.A2EndpointCommutator.exists_quotient_word_factorization
