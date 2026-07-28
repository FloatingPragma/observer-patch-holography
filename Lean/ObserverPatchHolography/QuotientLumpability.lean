import Mathlib

namespace OPH.QuotientLumpability

open Finset

/-! # Exact quotient lumpability: the receipt and its certified checker

Issue #592 companion lane.  The fractional quotient-sector sandbox
(`code/particles/fractional/build_fractional_quotient_receipts.py`) records a
`QUOTIENT_LUMPABILITY` readiness gate.  This file supplies the decision
procedure that gate names: a machine-checked definition of strong
lumpability for exact finite kernels, the quotient-kernel construction, the
push-forward commutation theorem it licenses, and a decidable checker that
is sound, complete, and returns an explicit failing representative pair with
the fibre where the weights differ whenever it rejects.

CONTENT.

* `fiberWeight K q x c`: total `K`-weight from state `x` into the `q`-fibre
  of `c`.  Weights live in an arbitrary additive commutative monoid; the
  controls instantiate at `ℤ` (kernel-reducible exact arithmetic) and the
  probability layer at `ℚ`.  No floating point exists anywhere in the file.
* `StronglyLumpable K q`: whenever `q x = q x'`, the fibre weights from `x`
  and `x'` agree on every fibre.
* `quotientKernel K q rep`: the induced kernel on the quotient, through a
  section `rep`; `quotientKernel_section_independent` and
  `quotientKernel_well_defined` show lumpability makes the choice immaterial.
* `push_step_comm`: one presentation-level step commutes with push-forward
  under lumpability; `push_iterate_comm` extends this to `n` steps.
* `quotientKernel_isMarkov`: row-stochasticity is preserved (this needs only
  the section, not lumpability — without lumpability the quotient kernel is
  still stochastic but does not govern the pushed dynamics).
  `quotientKernel_nonneg` is the separate statement for the simulator's
  nonnegative-rate convention: nonnegativity is preserved; no row-sum claim
  is made or true in that convention.
* `lumpabilityWitness K q l lc`: the certificate-producing checker over
  explicit enumerations `l`, `lc` of the state and class types.  `none`
  means lumpable; `some (x, x', c)` hands back a concrete pair with
  `q x = q x'` whose weights into the fibre of `c` differ.
* `lumpabilityWitness_sound`, `lumpabilityWitness_none`,
  `checkLumpable_iff`: the checker is sound and complete whenever the
  enumerations are complete.
* `stronglyLumpable_map`: lumpability transports along any additive monoid
  homomorphism of weights — this carries the `ℤ`-certified control to its
  `ℚ`-probability normalisation.
* Concrete controls on `Fin 4 → Fin 2`: the rate kernel `Wgood` is accepted
  and its quotient kernel computed; `Wbad` is rejected with the named pair
  `(0, 1)` at fibre `0` and weights `12 ≠ 6`.  Both directions run by
  kernel `decide` — no `native_decide`.  `KgoodQ := Wgood / 24` is the
  Markov normalisation, certified lumpable by transport and row-stochastic
  by `norm_num`.
* `fractionalSandboxKernel` is the exact-rational encoding of the concrete
  three-state fractional quotient sandbox in
  `muellerberndt/oph-physics-sim@87767593`, `oph_fractional/compare.py:47-61`.
  The checker accepts it and computes its two-state quotient kernel.

BOUNDARY.  Finite state and quotient types only, exact arithmetic only.
`push_step_comm` is a single-step statement and `push_iterate_comm` its
finite iteration; nothing here concerns continuous time, infinite
presentations, or floating-point kernels.  This file gives the
`QUOTIENT_LUMPABILITY` receipt a decision procedure and certifies the
simulator's concrete fractional sandbox instance below.  That instance is
explicitly a diagnostic toy presentation, not a material Hamiltonian, a
physical H3/KMS cell system, or the physical campaign of issue #592. -/

variable {S Q : Type*} [Fintype S] [DecidableEq Q]

/-- Total `K`-weight from `x` into the `q`-fibre of `c`. -/
def fiberWeight {R : Type*} [AddCommMonoid R]
    (K : S → S → R) (q : S → Q) (x : S) (c : Q) : R :=
  ∑ y ∈ univ.filter (fun y => q y = c), K x y

/-- Strong lumpability: states identified by `q` place equal total weight
into every quotient fibre. -/
def StronglyLumpable {R : Type*} [AddCommMonoid R]
    (K : S → S → R) (q : S → Q) : Prop :=
  ∀ x x', q x = q x' → ∀ c, fiberWeight K q x c = fiberWeight K q x' c

/-- Summing any function fibre-by-fibre recovers the total sum. -/
theorem sum_fibres [Fintype Q] {R : Type*} [AddCommMonoid R]
    (q : S → Q) (f : S → R) :
    ∑ c, ∑ x ∈ univ.filter (fun x => q x = c), f x = ∑ x, f x := by
  simp only [Finset.sum_filter]
  rw [Finset.sum_comm]
  simp

/-- The quotient kernel through a section `rep : Q → S`. -/
def quotientKernel {R : Type*} [AddCommMonoid R]
    (K : S → S → R) (q : S → Q) (rep : Q → S) : Q → Q → R :=
  fun c c' => fiberWeight K q (rep c) c'

/-- Under lumpability the quotient kernel does not depend on the section. -/
theorem quotientKernel_section_independent {R : Type*} [AddCommMonoid R]
    {K : S → S → R} {q : S → Q} (hL : StronglyLumpable K q)
    (rep rep' : Q → S) (hrep : ∀ c, q (rep c) = c) (hrep' : ∀ c, q (rep' c) = c) :
    quotientKernel K q rep = quotientKernel K q rep' := by
  funext c c'
  exact hL _ _ ((hrep c).trans (hrep' c).symm) c'

/-- Under lumpability the quotient kernel row at `c` computes the fibre
weights of every state in the fibre of `c`, not just the representative. -/
theorem quotientKernel_well_defined {R : Type*} [AddCommMonoid R]
    {K : S → S → R} {q : S → Q} (hL : StronglyLumpable K q)
    {rep : Q → S} (hrep : ∀ c, q (rep c) = c)
    {x : S} {c : Q} (hx : q x = c) (c' : Q) :
    quotientKernel K q rep c c' = fiberWeight K q x c' :=
  hL _ _ ((hrep c).trans hx.symm) c'

/-- One evolution step of a (row-vector) measure under a kernel. -/
def step {α R : Type*} [Fintype α] [NonUnitalNonAssocSemiring R]
    (K : α → α → R) (μ : α → R) : α → R :=
  fun y => ∑ x, μ x * K x y

/-- Push a measure forward along the quotient map. -/
def push {R : Type*} [AddCommMonoid R] (q : S → Q) (μ : S → R) : Q → R :=
  fun c => ∑ x ∈ univ.filter (fun x => q x = c), μ x

/-- **Push-forward commutes with one step.** Under strong lumpability,
stepping on the presentation and pushing forward equals pushing forward and
stepping with the quotient kernel. -/
theorem push_step_comm [Fintype Q] {R : Type*} [NonUnitalNonAssocSemiring R]
    {K : S → S → R} {q : S → Q}
    (hL : StronglyLumpable K q) {rep : Q → S} (hrep : ∀ c, q (rep c) = c)
    (μ : S → R) :
    push q (step K μ) = step (quotientKernel K q rep) (push q μ) := by
  funext c'
  have expand : push q (step K μ) c' = ∑ x, μ x * fiberWeight K q x c' := by
    simp only [push, step, fiberWeight]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun x _ => (Finset.mul_sum _ _ _).symm
  have fibres : ∑ x, μ x * fiberWeight K q x c'
      = ∑ c, ∑ x ∈ univ.filter (fun x => q x = c), μ x * fiberWeight K q x c' :=
    (sum_fibres q _).symm
  rw [expand, fibres]
  simp only [step, quotientKernel, push, Finset.sum_mul]
  refine Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun x hx => ?_
  have hxc : q x = c := (Finset.mem_filter.mp hx).2
  rw [hL x (rep c) (hxc.trans (hrep c).symm) c']

/-- Push-forward commutes with `n` steps. -/
theorem push_iterate_comm [Fintype Q] {R : Type*} [NonUnitalNonAssocSemiring R]
    {K : S → S → R} {q : S → Q}
    (hL : StronglyLumpable K q) {rep : Q → S} (hrep : ∀ c, q (rep c) = c)
    (μ : S → R) (n : ℕ) :
    push q ((step K)^[n] μ) = (step (quotientKernel K q rep))^[n] (push q μ) := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [Function.iterate_succ_apply', Function.iterate_succ_apply',
        push_step_comm hL hrep, ih]

/-- Row-stochastic exact kernel: nonnegative entries, rows summing to `1`.
Stated for exact rationals — the probability convention. -/
def IsMarkov {α : Type*} [Fintype α] (K : α → α → ℚ) : Prop :=
  (∀ x y, 0 ≤ K x y) ∧ ∀ x, ∑ y, K x y = 1

/-- Nonnegative rate kernel: the simulator's weaker convention. -/
def IsNonnegKernel {α β R : Type*} [Zero R] [LE R] (K : α → β → R) : Prop :=
  ∀ x y, 0 ≤ K x y

/-- The quotient of a Markov kernel through any section is Markov.  Note
this needs only the section, not lumpability; lumpability is what makes the
quotient kernel govern the pushed dynamics (`push_step_comm`). -/
theorem quotientKernel_isMarkov [Fintype Q] {K : S → S → ℚ} {q : S → Q}
    (hK : IsMarkov K) (rep : Q → S) :
    IsMarkov (quotientKernel K q rep) := by
  constructor
  · intro c c'
    exact Finset.sum_nonneg fun y _ => hK.1 (rep c) y
  · intro c
    calc ∑ c', quotientKernel K q rep c c'
        = ∑ y, K (rep c) y := sum_fibres q (K (rep c))
      _ = 1 := hK.2 (rep c)

/-- The quotient of a nonnegative rate kernel is nonnegative.  This is the
whole of what the rate convention preserves: no row-sum statement holds or
is claimed. -/
theorem quotientKernel_nonneg {R : Type*}
    [AddCommMonoid R] [Preorder R] [AddLeftMono R]
    {K : S → S → R} {q : S → Q} (hK : IsNonnegKernel K) (rep : Q → S) :
    IsNonnegKernel (quotientKernel K q rep) :=
  fun c _c' => Finset.sum_nonneg fun y _ => hK (rep c) y

/-- Fibre weights transport along additive monoid homomorphisms. -/
theorem fiberWeight_map {R R' : Type*} [AddCommMonoid R] [AddCommMonoid R']
    (f : R →+ R') (K : S → S → R) (q : S → Q) (x : S) (c : Q) :
    fiberWeight (fun a b => f (K a b)) q x c = f (fiberWeight K q x c) :=
  (map_sum f _ _).symm

/-- Strong lumpability transports along additive monoid homomorphisms of
the weight monoid — e.g. from an integer rate kernel to its rational
probability normalisation. -/
theorem stronglyLumpable_map {R R' : Type*} [AddCommMonoid R] [AddCommMonoid R']
    (f : R →+ R') {K : S → S → R} {q : S → Q} (hL : StronglyLumpable K q) :
    StronglyLumpable (fun a b => f (K a b)) q := by
  intro x x' hq c
  rw [fiberWeight_map, fiberWeight_map, hL x x' hq c]

/-! ## The certified checker -/

/-- Certificate-producing lumpability checker over explicit enumerations
`l : List S` and `lc : List Q`.  Returns `none` when no violation is
listed, and otherwise `some (x, x', c)`: an explicit pair of states in one
fibre whose total weights into the fibre of `c` differ. -/
def lumpabilityWitness {R : Type*} [AddCommMonoid R] [DecidableEq R]
    (K : S → S → R) (q : S → Q) (l : List S) (lc : List Q) :
    Option (S × S × Q) :=
  (l.flatMap fun x => l.flatMap fun x' => lc.map fun c => (x, x', c)).find?
    fun t => decide (q t.1 = q t.2.1) &&
      decide (fiberWeight K q t.1 t.2.2 ≠ fiberWeight K q t.2.1 t.2.2)

/-- Soundness: a returned witness really is a failing pair with its fibre.
No completeness hypothesis on the enumerations is needed. -/
theorem lumpabilityWitness_sound {R : Type*} [AddCommMonoid R] [DecidableEq R]
    {K : S → S → R} {q : S → Q} {l : List S} {lc : List Q} {x x' : S} {c : Q}
    (h : lumpabilityWitness K q l lc = some (x, x', c)) :
    q x = q x' ∧ fiberWeight K q x c ≠ fiberWeight K q x' c := by
  have hp := List.find?_some h
  simpa only [Bool.and_eq_true, decide_eq_true_eq] using hp

/-- Completeness: if the enumerations cover both types and the checker
finds nothing, the kernel is strongly lumpable. -/
theorem lumpabilityWitness_none {R : Type*} [AddCommMonoid R] [DecidableEq R]
    {K : S → S → R} {q : S → Q} {l : List S} {lc : List Q}
    (hl : ∀ x : S, x ∈ l) (hlc : ∀ c : Q, c ∈ lc)
    (h : lumpabilityWitness K q l lc = none) :
    StronglyLumpable K q := by
  intro x x' hq c
  by_contra hne
  have hmem : (x, x', c) ∈
      (l.flatMap fun x => l.flatMap fun x' => lc.map fun c => (x, x', c)) :=
    List.mem_flatMap.mpr ⟨x, hl x,
      List.mem_flatMap.mpr ⟨x', hl x', List.mem_map.mpr ⟨c, hlc c, rfl⟩⟩⟩
  have hnp := List.find?_eq_none.mp h _ hmem
  simp only [Bool.and_eq_true, decide_eq_true_eq] at hnp
  exact hnp ⟨hq, hne⟩

/-- Boolean face of the checker. -/
def checkLumpable {R : Type*} [AddCommMonoid R] [DecidableEq R]
    (K : S → S → R) (q : S → Q) (l : List S) (lc : List Q) : Bool :=
  (lumpabilityWitness K q l lc).isNone

/-- **The checker decides strong lumpability**: sound and complete over
complete enumerations. -/
theorem checkLumpable_iff {R : Type*} [AddCommMonoid R] [DecidableEq R]
    (K : S → S → R) (q : S → Q) {l : List S} {lc : List Q}
    (hl : ∀ x : S, x ∈ l) (hlc : ∀ c : Q, c ∈ lc) :
    checkLumpable K q l lc = true ↔ StronglyLumpable K q := by
  constructor
  · intro h
    exact lumpabilityWitness_none hl hlc (Option.isNone_iff_eq_none.mp h)
  · intro hL
    cases hw : lumpabilityWitness K q l lc with
    | none => simp [checkLumpable, hw]
    | some t =>
        obtain ⟨x, x', c⟩ := t
        have hs := lumpabilityWitness_sound hw
        exact absurd (hL x x' hs.1 c) hs.2

/-- Strong lumpability is decidable outright on finite data. -/
instance {R : Type*} [AddCommMonoid R] [DecidableEq R] [Fintype Q]
    (K : S → S → R) (q : S → Q) : Decidable (StronglyLumpable K q) := by
  unfold StronglyLumpable
  infer_instance

/-! ## Controls

`Fin 4` states, `Fin 2` quotient classes: `{0, 1} ↦ 0`, `{2, 3} ↦ 1`.
Rate kernels over `ℤ` so every claim runs by kernel `decide`; the Markov
normalisation over `ℚ` is derived by homomorphism transport. -/

/-- Quotient map folding pairs of states. -/
@[reducible] def qHalf : Fin 4 → Fin 2 := ![0, 0, 1, 1]

/-- Section for `qHalf`. -/
@[reducible] def repHalf : Fin 2 → Fin 4 := ![0, 2]

example : ∀ c, qHalf (repHalf c) = c := by decide

/-- Explicit enumeration of the state type. -/
@[reducible] def states4 : List (Fin 4) := List.finRange 4

/-- Explicit enumeration of the class type. -/
@[reducible] def classes2 : List (Fin 2) := List.finRange 2

theorem states4_complete : ∀ x : Fin 4, x ∈ states4 := by decide

theorem classes2_complete : ∀ c : Fin 2, c ∈ classes2 := by decide

/-- A strongly lumpable nonnegative rate kernel: rows within one fibre
differ entrywise but agree on fibre totals `(12, 12)` and `(6, 18)`. -/
@[reducible] def Wgood : Fin 4 → Fin 4 → ℤ :=
  ![![4, 8, 6, 6],
    ![8, 4, 12, 0],
    ![3, 3, 9, 9],
    ![0, 6, 12, 6]]

/-- Non-lumpable perturbation: row `1` moves fibre-`0` weight into
fibre `1`, breaking agreement with row `0`. -/
@[reducible] def Wbad : Fin 4 → Fin 4 → ℤ :=
  ![![4, 8, 6, 6],
    ![3, 3, 18, 0],
    ![3, 3, 9, 9],
    ![0, 6, 12, 6]]

example : IsNonnegKernel Wgood ∧ IsNonnegKernel Wbad := by
  unfold IsNonnegKernel
  constructor <;> decide

/-- Positive control: the checker accepts `Wgood`. -/
theorem wgood_accepted : checkLumpable Wgood qHalf states4 classes2 = true := by
  decide

/-- The positive control is strongly lumpable, through the certified
checker. -/
theorem wgood_lumpable : StronglyLumpable Wgood qHalf :=
  lumpabilityWitness_none states4_complete classes2_complete (by decide)

/-- The computed quotient kernel of the positive control. -/
theorem wgood_quotient :
    ∀ c c', quotientKernel Wgood qHalf repHalf c c'
      = ![![(12 : ℤ), 12], ![6, 18]] c c' := by decide

/-- Negative control: the checker rejects `Wbad` and names the failing
representative pair `(0, 1)` and the fibre `0` where the weights differ. -/
theorem wbad_rejected :
    lumpabilityWitness Wbad qHalf states4 classes2 = some (0, 1, 0) := by
  decide

/-- The named pair really fails: `q`-equal states, unequal fibre weights.
Extracted through the soundness theorem, not recomputed. -/
theorem wbad_witness_valid :
    qHalf 0 = qHalf 1 ∧
      fiberWeight Wbad qHalf 0 0 ≠ fiberWeight Wbad qHalf 1 0 :=
  lumpabilityWitness_sound wbad_rejected

/-- The failing weights themselves: `12` from state `0`, `6` from state `1`,
into fibre `0`. -/
example : fiberWeight Wbad qHalf 0 (0 : Fin 2) = 12
    ∧ fiberWeight Wbad qHalf 1 (0 : Fin 2) = 6 := by
  constructor <;> decide

/-- The negative control is not strongly lumpable. -/
theorem wbad_not_lumpable : ¬ StronglyLumpable Wbad qHalf :=
  fun hL => wbad_witness_valid.2 (hL 0 1 wbad_witness_valid.1 0)

theorem wbad_check_false : checkLumpable Wbad qHalf states4 classes2 = false := by
  decide

/-! ### The Markov normalisation over `ℚ` -/

/-- Division by the total rate `24`, as an additive monoid homomorphism
`ℤ →+ ℚ`. -/
def div24 : ℤ →+ ℚ where
  toFun n := (n : ℚ) / 24
  map_zero' := by norm_num
  map_add' a b := by push_cast; ring

/-- The probability normalisation of the positive control. -/
def KgoodQ : Fin 4 → Fin 4 → ℚ := fun x y => div24 (Wgood x y)

/-- The `ℤ`-certified lumpability transports to the `ℚ` normalisation. -/
theorem kgoodQ_lumpable : StronglyLumpable KgoodQ qHalf :=
  stronglyLumpable_map div24 wgood_lumpable

/-- The normalised control is a genuine Markov kernel, so the controls
exercise the probability convention, not vacuous cases. -/
theorem kgoodQ_isMarkov : IsMarkov KgoodQ := by
  constructor
  · intro x y
    fin_cases x <;> fin_cases y <;> norm_num [KgoodQ, div24, Wgood]
  · intro x
    fin_cases x <;>
      norm_num [KgoodQ, div24, Wgood, Fin.sum_univ_four,
        Matrix.cons_val_two, Matrix.cons_val_three]

/-- The checker accepts the `ℚ` control as well, through the
soundness/completeness theorem rather than kernel computation on `ℚ`. -/
theorem kgoodQ_accepted : checkLumpable KgoodQ qHalf states4 classes2 = true :=
  (checkLumpable_iff KgoodQ qHalf states4_complete classes2_complete).mpr
    kgoodQ_lumpable

/-- The quotient of the normalised control is Markov. -/
theorem kgoodQ_quotient_isMarkov :
    IsMarkov (quotientKernel KgoodQ qHalf repHalf) :=
  quotientKernel_isMarkov kgoodQ_isMarkov repHalf

/-! ## The OPH fractional quotient sandbox

This is a literal exact-rational transcription of the simulator data at
`muellerberndt/oph-physics-sim@877675938812e26417bb006d5bf3752301f8a3f8`,
`oph_fractional/compare.py:47-61`:

* representatives, in source order: `a0`, `a1`, `vac`;
* quotient map: `a0 ↦ anyon_e_over_3`, `a1 ↦ anyon_e_over_3`,
  `vac ↦ vacuum`;
* transition rows: `a0 ↦ {a1: 0.5, vac: 0.5}`,
  `a1 ↦ {a0: 0.5, vac: 0.5}`, `vac ↦ {vac: 1.0}`.

The omitted dictionary entries are zero, exactly as the simulator's
`quotient_lumpability` implementation initializes missing sector weight to
zero (`oph_fractional/quotient.py:64-77`).  Decimal `0.5` is exactly `1 / 2`
as a Python binary float, and `1.0` is exactly `1`; those are encoded in `ℚ`
below.  This is the repository-linked diagnostic
`twisted_tmd_fractional_sandbox` at regulator `12`.  It is not a real material
sample and not the physical H3/KMS presentation requested by issue #592. -/

inductive FractionalSandboxState
  | a0
  | a1
  | vac
  deriving DecidableEq, Fintype

inductive FractionalSandboxSector
  | anyonEOver3
  | vacuum
  deriving DecidableEq, Fintype

/-- The simulator's `canonical`/`quotient_map` dictionary. -/
@[reducible] def fractionalSandboxQuotient :
    FractionalSandboxState → FractionalSandboxSector
  | .a0 => .anyonEOver3
  | .a1 => .anyonEOver3
  | .vac => .vacuum

/-- The simulator's three explicit transition rows, with absent entries zero. -/
@[reducible] def fractionalSandboxKernel :
    FractionalSandboxState → FractionalSandboxState → ℚ
  | .a0, .a1 => 1 / 2
  | .a0, .vac => 1 / 2
  | .a1, .a0 => 1 / 2
  | .a1, .vac => 1 / 2
  | .vac, .vac => 1
  | _, _ => 0

@[reducible] def fractionalSandboxStates : List FractionalSandboxState :=
  [.a0, .a1, .vac]

@[reducible] def fractionalSandboxSectors : List FractionalSandboxSector :=
  [.anyonEOver3, .vacuum]

/-- Auxiliary Lean section for displaying the quotient kernel.  The simulator
does not select a section: `a0` is chosen here for the anyon fibre and `vac`
for the singleton vacuum fibre.  Lumpability makes the result independent of
this auxiliary choice. -/
@[reducible] def fractionalSandboxRep :
    FractionalSandboxSector → FractionalSandboxState
  | .anyonEOver3 => .a0
  | .vacuum => .vac

theorem fractionalSandboxStates_complete :
    ∀ x : FractionalSandboxState, x ∈ fractionalSandboxStates := by
  decide

theorem fractionalSandboxSectors_complete :
    ∀ c : FractionalSandboxSector, c ∈ fractionalSandboxSectors := by
  decide

theorem fractionalSandboxRep_isSection :
    ∀ c, fractionalSandboxQuotient (fractionalSandboxRep c) = c := by
  decide

/-- Exact fibre sums show strong lumpability of the simulator kernel. -/
theorem fractionalSandbox_lumpable :
    StronglyLumpable fractionalSandboxKernel fractionalSandboxQuotient := by
  classical
  have huniv : (univ : Finset FractionalSandboxState) =
      {.a0, .a1, .vac} := by
    decide
  have hsum (f : FractionalSandboxState → ℚ) :
      ∑ x, f x = f .a0 + (f .a1 + f .vac) := by
    rw [huniv]
    rw [Finset.sum_insert
      (by decide : .a0 ∉ ({.a1, .vac} : Finset FractionalSandboxState))]
    rw [Finset.sum_insert
      (by decide : .a1 ∉ ({.vac} : Finset FractionalSandboxState))]
    simp
  intro x x' hq c
  fin_cases x <;> fin_cases x' <;> fin_cases c <;>
    simp_all [fiberWeight, Finset.sum_filter, fractionalSandboxKernel,
      fractionalSandboxQuotient]

/-- `QUOTIENT_LUMPABILITY` for the simulator's concrete fractional sandbox:
the certified checker returns `true` on its exact data. -/
theorem fractionalSandbox_accepted :
    checkLumpable fractionalSandboxKernel fractionalSandboxQuotient
      fractionalSandboxStates fractionalSandboxSectors = true :=
  (checkLumpable_iff fractionalSandboxKernel fractionalSandboxQuotient
    fractionalSandboxStates_complete fractionalSandboxSectors_complete).mpr
      fractionalSandbox_lumpable

/-- The source data is a genuine Markov kernel, not merely nonnegative weights. -/
theorem fractionalSandbox_isMarkov : IsMarkov fractionalSandboxKernel := by
  classical
  have huniv : (univ : Finset FractionalSandboxState) =
      {.a0, .a1, .vac} := by
    decide
  have hsum (f : FractionalSandboxState → ℚ) :
      ∑ x, f x = f .a0 + (f .a1 + f .vac) := by
    rw [huniv]
    rw [Finset.sum_insert
      (by decide : .a0 ∉ ({.a1, .vac} : Finset FractionalSandboxState))]
    rw [Finset.sum_insert
      (by decide : .a1 ∉ ({.vac} : Finset FractionalSandboxState))]
    simp
  constructor
  · intro x y
    fin_cases x <;> fin_cases y <;> norm_num [fractionalSandboxKernel]
  · intro x
    fin_cases x <;>
      norm_num [hsum, fractionalSandboxKernel]

/-- The quotient transition matrix reported by exact fibre summation:
the anyon sector stays anyon with probability `1/2` and moves to vacuum with
probability `1/2`; vacuum is absorbing. -/
theorem fractionalSandbox_quotientKernel :
    ∀ c c', quotientKernel fractionalSandboxKernel fractionalSandboxQuotient
      fractionalSandboxRep c c' =
        match c, c' with
        | .anyonEOver3, .anyonEOver3 => 1 / 2
        | .anyonEOver3, .vacuum => 1 / 2
        | .vacuum, .anyonEOver3 => 0
        | .vacuum, .vacuum => 1 := by
  classical
  have huniv : (univ : Finset FractionalSandboxState) =
      {.a0, .a1, .vac} := by
    decide
  have hsum (f : FractionalSandboxState → ℚ) :
      ∑ x, f x = f .a0 + (f .a1 + f .vac) := by
    rw [huniv]
    rw [Finset.sum_insert
      (by decide : .a0 ∉ ({.a1, .vac} : Finset FractionalSandboxState))]
    rw [Finset.sum_insert
      (by decide : .a1 ∉ ({.vac} : Finset FractionalSandboxState))]
    simp
  intro c c'
  fin_cases c <;> fin_cases c' <;>
    simp [quotientKernel, fiberWeight, Finset.sum_filter, hsum,
      fractionalSandboxKernel, fractionalSandboxQuotient, fractionalSandboxRep]

#print axioms push_step_comm
#print axioms push_iterate_comm
#print axioms quotientKernel_section_independent
#print axioms quotientKernel_isMarkov
#print axioms quotientKernel_nonneg
#print axioms stronglyLumpable_map
#print axioms lumpabilityWitness_sound
#print axioms lumpabilityWitness_none
#print axioms checkLumpable_iff
#print axioms wgood_accepted
#print axioms wgood_lumpable
#print axioms wgood_quotient
#print axioms wbad_rejected
#print axioms wbad_witness_valid
#print axioms wbad_not_lumpable
#print axioms kgoodQ_lumpable
#print axioms kgoodQ_isMarkov
#print axioms kgoodQ_quotient_isMarkov
#print axioms fractionalSandboxStates_complete
#print axioms fractionalSandboxSectors_complete
#print axioms fractionalSandboxRep_isSection
#print axioms fractionalSandbox_accepted
#print axioms fractionalSandbox_lumpable
#print axioms fractionalSandbox_isMarkov
#print axioms fractionalSandbox_quotientKernel

end OPH.QuotientLumpability
