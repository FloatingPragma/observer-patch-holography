import Mathlib

namespace OPH.QuotientLumpability

open Finset

/-! # Exact quotient receipts and their certified checkers

Issue #592 companion lane.  The fractional quotient-sector sandbox
(`code/particles/fractional/build_fractional_quotient_receipts.py`) records a
`QUOTIENT_LUMPABILITY` readiness gate and its sibling quotient-correctness
gates.  This file supplies their source-grounded decision procedures,
beginning with a machine-checked definition of strong
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
* `canonicalizerWitness` and `representativeWitness`: certificate-producing
  checkers for the sibling simulator receipts `CANONICALIZER_IDEMPOTENCE`
  and `REPRESENTATIVE_INVARIANCE`.  The latter follows the simulator source:
  an explicitly supplied observable must be constant on each quotient fibre.
* `orbitSizeBiasWitness`: the exact-rational counterpart of
  `NO_ORBIT_SIZE_BIAS`.  It rejects either a nonpositive declared orbit size
  or sector weights proportional to genuinely varying hidden orbit sizes,
  and distinguishes those two failure witnesses.
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
presentations, or floating-point kernels.  This file gives
`QUOTIENT_LUMPABILITY`, `CANONICALIZER_IDEMPOTENCE`,
`REPRESENTATIVE_INVARIANCE`, and `NO_ORBIT_SIZE_BIAS` decision procedures and
certifies the simulator's concrete fractional sandbox instance below.
`MATERIAL_QUOTIENT_NORMAL_FORM_RECEIPT` is intentionally not reinterpreted:
its four repair-status inputs are absent from this Lean instance.  The
instance here is
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

/-! ## Certified checkers for the sibling quotient receipts -/

/-! ### Canonicalizer idempotence -/

/-- A canonicalizer is idempotent when a second canonicalization changes
nothing.  This is the mathematical property tested by the simulator's
`canonicalizer_idempotence`. -/
def CanonicalizerIdempotent {α : Type*} (canon : α → α) : Prop :=
  ∀ x, canon (canon x) = canon x

/-- A state on which canonicalizing twice differs from canonicalizing once. -/
def canonicalizerWitness {α : Type*} [DecidableEq α]
    (canon : α → α) (l : List α) : Option α :=
  l.find? fun x => decide (canon (canon x) ≠ canon x)

/-- A returned canonicalizer witness is a genuine idempotence failure. -/
theorem canonicalizerWitness_sound {α : Type*} [DecidableEq α]
    {canon : α → α} {l : List α} {x : α}
    (h : canonicalizerWitness canon l = some x) :
    canon (canon x) ≠ canon x := by
  have hp := List.find?_some h
  simpa only [decide_eq_true_eq] using hp

/-- A complete enumeration with no witness proves idempotence. -/
theorem canonicalizerWitness_none {α : Type*} [DecidableEq α]
    {canon : α → α} {l : List α} (hl : ∀ x : α, x ∈ l)
    (h : canonicalizerWitness canon l = none) :
    CanonicalizerIdempotent canon := by
  intro x
  by_contra hne
  have hnp := List.find?_eq_none.mp h x (hl x)
  exact hnp (by simpa only [decide_eq_true_eq] using hne)

/-- Boolean face of the canonicalizer checker. -/
def checkCanonicalizerIdempotent {α : Type*} [DecidableEq α]
    (canon : α → α) (l : List α) : Bool :=
  (canonicalizerWitness canon l).isNone

/-- The canonicalizer checker is sound and complete on a complete list. -/
theorem checkCanonicalizerIdempotent_iff {α : Type*} [DecidableEq α]
    (canon : α → α) {l : List α} (hl : ∀ x : α, x ∈ l) :
    checkCanonicalizerIdempotent canon l = true ↔
      CanonicalizerIdempotent canon := by
  constructor
  · intro h
    exact canonicalizerWitness_none hl (Option.isNone_iff_eq_none.mp h)
  · intro hcanon
    cases hw : canonicalizerWitness canon l with
    | none => simp [checkCanonicalizerIdempotent, hw]
    | some x =>
        exact absurd (hcanon x) (canonicalizerWitness_sound hw)

/-! ### Representative invariance -/

/-- An observable is representative-invariant when it is constant on every
fibre of the quotient map.  This is exactly the generic predicate implemented
by the simulator's `representative_invariance(schema, observable)`. -/
def RepresentativeInvariant {α β Ω : Type*}
    (q : α → β) (observable : α → Ω) : Prop :=
  ∀ x x', q x = q x' → observable x = observable x'

/-- A pair in one quotient fibre on which the observable differs. -/
def representativeWitness {α β Ω : Type*}
    [DecidableEq β] [DecidableEq Ω]
    (q : α → β) (observable : α → Ω) (l : List α) : Option (α × α) :=
  (l.flatMap fun x => l.map fun x' => (x, x')).find? fun p =>
    decide (q p.1 = q p.2) && decide (observable p.1 ≠ observable p.2)

/-- A returned pair is a genuine representative-invariance violation. -/
theorem representativeWitness_sound {α β Ω : Type*}
    [DecidableEq β] [DecidableEq Ω]
    {q : α → β} {observable : α → Ω} {l : List α} {x x' : α}
    (h : representativeWitness q observable l = some (x, x')) :
    q x = q x' ∧ observable x ≠ observable x' := by
  have hp := List.find?_some h
  simpa only [Bool.and_eq_true, decide_eq_true_eq] using hp

/-- A complete enumeration with no violating pair proves invariance. -/
theorem representativeWitness_none {α β Ω : Type*}
    [DecidableEq β] [DecidableEq Ω]
    {q : α → β} {observable : α → Ω} {l : List α}
    (hl : ∀ x : α, x ∈ l) (h : representativeWitness q observable l = none) :
    RepresentativeInvariant q observable := by
  intro x x' hq
  by_contra hne
  have hmem : (x, x') ∈ (l.flatMap fun x => l.map fun x' => (x, x')) :=
    List.mem_flatMap.mpr
      ⟨x, hl x, List.mem_map.mpr ⟨x', hl x', rfl⟩⟩
  have hnp := List.find?_eq_none.mp h _ hmem
  simp only [Bool.and_eq_true, decide_eq_true_eq] at hnp
  exact hnp ⟨hq, hne⟩

/-- Boolean face of the representative-invariance checker. -/
def checkRepresentativeInvariant {α β Ω : Type*}
    [DecidableEq β] [DecidableEq Ω]
    (q : α → β) (observable : α → Ω) (l : List α) : Bool :=
  (representativeWitness q observable l).isNone

/-- The representative-invariance checker is sound and complete on a
complete list. -/
theorem checkRepresentativeInvariant_iff {α β Ω : Type*}
    [DecidableEq β] [DecidableEq Ω]
    (q : α → β) (observable : α → Ω) {l : List α}
    (hl : ∀ x : α, x ∈ l) :
    checkRepresentativeInvariant q observable l = true ↔
      RepresentativeInvariant q observable := by
  constructor
  · intro h
    exact representativeWitness_none hl (Option.isNone_iff_eq_none.mp h)
  · intro hinv
    cases hw : representativeWitness q observable l with
    | none => simp [checkRepresentativeInvariant, hw]
    | some p =>
        obtain ⟨x, x'⟩ := p
        have hs := representativeWitness_sound hw
        exact absurd (hinv x x' hs.1) hs.2

/-! ### No hidden-orbit-size bias -/

/-- Sector weight per hidden representative, using exact rational arithmetic.
Orbit sizes are integers so malformed nonpositive declarations remain
expressible and rejectable. -/
def perRepresentativeWeight {β : Type*}
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) (c : β) : ℚ :=
  sectorWeight c / orbitSize c

/-- Exact form of "sector weights are proportional to hidden orbit sizes":
some declared orbit sizes genuinely vary, at least two sectors have positive
per-representative weight, and every positive per-representative weight is
the same.  The four existential fields retain explicit witnesses for both
the size variation and the nonvacuous positive comparison. -/
def TracksHiddenOrbitCount {β : Type*}
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) : Prop :=
  ∃ sizeLeft sizeRight positiveLeft positiveRight,
    orbitSize sizeLeft ≠ orbitSize sizeRight ∧
    positiveLeft ≠ positiveRight ∧
    0 < perRepresentativeWeight orbitSize sectorWeight positiveLeft ∧
    0 < perRepresentativeWeight orbitSize sectorWeight positiveRight ∧
    ∀ c, 0 < perRepresentativeWeight orbitSize sectorWeight c →
      perRepresentativeWeight orbitSize sectorWeight c =
        perRepresentativeWeight orbitSize sectorWeight positiveLeft

/-- Exact no-orbit-size-bias property implemented by the simulator: every
declared orbit size is positive, and the sector weights do not track a
varying hidden count. -/
def NoOrbitSizeBias {β : Type*}
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) : Prop :=
  (∀ c, 0 < orbitSize c) ∧ ¬ TracksHiddenOrbitCount orbitSize sectorWeight

/-- Explicit failure certificate for `NoOrbitSizeBias`. -/
inductive OrbitSizeBiasWitness (β : Type*)
  | nonpositive (sector : β)
  | tracksHiddenCount
      (sizeLeft sizeRight positiveLeft positiveRight : β)
  deriving DecidableEq

/-- The mathematical proposition certified by each orbit-size-bias witness. -/
def ValidOrbitSizeBiasWitness {β : Type*}
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) :
    OrbitSizeBiasWitness β → Prop
  | .nonpositive c => ¬ 0 < orbitSize c
  | .tracksHiddenCount sizeLeft sizeRight positiveLeft positiveRight =>
      orbitSize sizeLeft ≠ orbitSize sizeRight ∧
      positiveLeft ≠ positiveRight ∧
      0 < perRepresentativeWeight orbitSize sectorWeight positiveLeft ∧
      0 < perRepresentativeWeight orbitSize sectorWeight positiveRight ∧
      ∀ c, 0 < perRepresentativeWeight orbitSize sectorWeight c →
        perRepresentativeWeight orbitSize sectorWeight c =
          perRepresentativeWeight orbitSize sectorWeight positiveLeft

/-- Search the four explicit fields of a hidden-count witness. -/
def hiddenOrbitCountWitness {β : Type*} [Fintype β] [DecidableEq β]
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) (l : List β) :
    Option (β × β × β × β) :=
  (l.flatMap fun sizeLeft =>
    l.flatMap fun sizeRight =>
    l.flatMap fun positiveLeft =>
    l.map fun positiveRight =>
      (sizeLeft, sizeRight, positiveLeft, positiveRight)).find? fun t =>
        decide (
          orbitSize t.1 ≠ orbitSize t.2.1 ∧
          t.2.2.1 ≠ t.2.2.2 ∧
          0 < perRepresentativeWeight orbitSize sectorWeight t.2.2.1 ∧
          0 < perRepresentativeWeight orbitSize sectorWeight t.2.2.2) &&
        l.all fun c => decide (
          0 < perRepresentativeWeight orbitSize sectorWeight c →
          perRepresentativeWeight orbitSize sectorWeight c =
            perRepresentativeWeight orbitSize sectorWeight t.2.2.1)

/-- A returned hidden-count tuple proves the global hidden-count property. -/
theorem hiddenOrbitCountWitness_sound {β : Type*} [Fintype β] [DecidableEq β]
    {orbitSize : β → ℤ} {sectorWeight : β → ℚ} {l : List β}
    {sizeLeft sizeRight positiveLeft positiveRight : β}
    (hl : ∀ c : β, c ∈ l)
    (h : hiddenOrbitCountWitness orbitSize sectorWeight l =
      some (sizeLeft, sizeRight, positiveLeft, positiveRight)) :
    ValidOrbitSizeBiasWitness orbitSize sectorWeight
      (.tracksHiddenCount sizeLeft sizeRight positiveLeft positiveRight) := by
  have hp := List.find?_some h
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at hp
  refine ⟨hp.1.1, hp.1.2.1, hp.1.2.2.1, hp.1.2.2.2, ?_⟩
  intro c
  exact hp.2 c (hl c)

/-- Complete lists make the hidden-count tuple search complete. -/
theorem hiddenOrbitCountWitness_none {β : Type*} [Fintype β] [DecidableEq β]
    {orbitSize : β → ℤ} {sectorWeight : β → ℚ} {l : List β}
    (hl : ∀ c : β, c ∈ l)
    (h : hiddenOrbitCountWitness orbitSize sectorWeight l = none) :
    ¬ TracksHiddenOrbitCount orbitSize sectorWeight := by
  rintro ⟨sizeLeft, sizeRight, positiveLeft, positiveRight, hvalid⟩
  have hmem : (sizeLeft, sizeRight, positiveLeft, positiveRight) ∈
      (l.flatMap fun sizeLeft =>
        l.flatMap fun sizeRight =>
        l.flatMap fun positiveLeft =>
        l.map fun positiveRight =>
          (sizeLeft, sizeRight, positiveLeft, positiveRight)) :=
    List.mem_flatMap.mpr ⟨sizeLeft, hl sizeLeft,
      List.mem_flatMap.mpr ⟨sizeRight, hl sizeRight,
        List.mem_flatMap.mpr ⟨positiveLeft, hl positiveLeft,
          List.mem_map.mpr ⟨positiveRight, hl positiveRight, rfl⟩⟩⟩⟩
  have hnp := List.find?_eq_none.mp h _ hmem
  apply hnp
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true]
  exact ⟨⟨hvalid.1, hvalid.2.1, hvalid.2.2.1, hvalid.2.2.2.1⟩,
    fun c _hc => hvalid.2.2.2.2 c⟩

/-- Certificate-producing checker for orbit-size bias.  Nonpositive sizes
take precedence; otherwise a four-sector tuple certifies hidden-count
tracking. -/
def orbitSizeBiasWitness {β : Type*} [Fintype β] [DecidableEq β]
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) (l : List β) :
    Option (OrbitSizeBiasWitness β) :=
  match l.find? fun c => decide (¬ 0 < orbitSize c) with
  | some c => some (.nonpositive c)
  | none =>
      match hiddenOrbitCountWitness orbitSize sectorWeight l with
      | some (sizeLeft, sizeRight, positiveLeft, positiveRight) =>
          some (.tracksHiddenCount sizeLeft sizeRight positiveLeft positiveRight)
      | none => none

/-- Every returned orbit-size-bias witness is valid. -/
theorem orbitSizeBiasWitness_sound {β : Type*} [Fintype β] [DecidableEq β]
    {orbitSize : β → ℤ} {sectorWeight : β → ℚ} {l : List β}
    {w : OrbitSizeBiasWitness β}
    (hl : ∀ c : β, c ∈ l)
    (h : orbitSizeBiasWitness orbitSize sectorWeight l = some w) :
    ValidOrbitSizeBiasWitness orbitSize sectorWeight w := by
  unfold orbitSizeBiasWitness at h
  split at h
  next c hc =>
    have hp := List.find?_some hc
    cases h
    simpa only [ValidOrbitSizeBiasWitness, decide_eq_true_eq] using hp
  next hnonpos =>
    split at h
    next sizeLeft sizeRight positiveLeft positiveRight ht =>
      cases h
      have htracks := hiddenOrbitCountWitness_sound
        (orbitSize := orbitSize) (sectorWeight := sectorWeight)
        (sizeLeft := sizeLeft) (sizeRight := sizeRight)
        (positiveLeft := positiveLeft) (positiveRight := positiveRight) hl ht
      exact htracks
    next hnone => simp at h

/-- A valid explicit witness refutes `NoOrbitSizeBias`. -/
theorem validOrbitSizeBiasWitness_not_free {β : Type*}
    {orbitSize : β → ℤ} {sectorWeight : β → ℚ}
    {w : OrbitSizeBiasWitness β}
    (h : ValidOrbitSizeBiasWitness orbitSize sectorWeight w) :
    ¬ NoOrbitSizeBias orbitSize sectorWeight := by
  intro hfree
  cases w with
  | nonpositive c => exact h (hfree.1 c)
  | tracksHiddenCount sizeLeft sizeRight positiveLeft positiveRight =>
      apply hfree.2
      exact ⟨sizeLeft, sizeRight, positiveLeft, positiveRight, h⟩

/-- On a complete list, no witness proves no orbit-size bias. -/
theorem orbitSizeBiasWitness_none {β : Type*} [Fintype β] [DecidableEq β]
    {orbitSize : β → ℤ} {sectorWeight : β → ℚ} {l : List β}
    (hl : ∀ c : β, c ∈ l)
    (h : orbitSizeBiasWitness orbitSize sectorWeight l = none) :
    NoOrbitSizeBias orbitSize sectorWeight := by
  unfold orbitSizeBiasWitness at h
  split at h
  next c hc => simp at h
  next hnonpos =>
    constructor
    · intro c
      have hnp := List.find?_eq_none.mp hnonpos c (hl c)
      simpa only [decide_eq_true_eq, not_not] using hnp
    · split at h
      next t ht => simp at h
      next hhidden =>
        exact hiddenOrbitCountWitness_none hl hhidden

/-- Boolean face of the no-orbit-size-bias checker. -/
def checkNoOrbitSizeBias {β : Type*} [Fintype β] [DecidableEq β]
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) (l : List β) : Bool :=
  (orbitSizeBiasWitness orbitSize sectorWeight l).isNone

/-- The orbit-size-bias checker is sound and complete on a complete list. -/
theorem checkNoOrbitSizeBias_iff {β : Type*} [Fintype β] [DecidableEq β]
    (orbitSize : β → ℤ) (sectorWeight : β → ℚ) {l : List β}
    (hl : ∀ c : β, c ∈ l) :
    checkNoOrbitSizeBias orbitSize sectorWeight l = true ↔
      NoOrbitSizeBias orbitSize sectorWeight := by
  constructor
  · intro h
    exact orbitSizeBiasWitness_none hl (Option.isNone_iff_eq_none.mp h)
  · intro hfree
    cases hw : orbitSizeBiasWitness orbitSize sectorWeight l with
    | none => simp [checkNoOrbitSizeBias, hw]
    | some w =>
        exact absurd hfree
          (validOrbitSizeBiasWitness_not_free
            (orbitSizeBiasWitness_sound hl hw))

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

/-! ### `CANONICALIZER_IDEMPOTENCE`

The pinned simulator checks `canonicalize(canonicalize(state)) =
canonicalize(state)` (`oph_fractional/quotient.py:39-47`) and applies that
check to this schema (`oph_fractional/compare.py:132`).  Its strings serve as
both states and canonical sector labels.  In this typed transcription states
and sectors are separate, so the same-type normal-form map is the section
after the quotient map, `rep ∘ q`. -/

/-- Same-type normal-form map corresponding to the simulator canonicalizer. -/
@[reducible] def fractionalSandboxNormalForm :
    FractionalSandboxState → FractionalSandboxState :=
  fun x => fractionalSandboxRep (fractionalSandboxQuotient x)

/-- A deliberately non-idempotent canonicalizer:
`a0 ↦ a1 ↦ vac`. -/
@[reducible] def fractionalSandboxBadCanonicalizer :
    FractionalSandboxState → FractionalSandboxState
  | .a0 => .a1
  | .a1 => .vac
  | .vac => .vac

/-- Negative control: the checker names `a0` as the idempotence failure. -/
theorem fractionalSandboxBadCanonicalizer_rejected :
    canonicalizerWitness fractionalSandboxBadCanonicalizer
      fractionalSandboxStates = some .a0 := by
  decide

/-- The returned canonicalizer witness is valid, through checker soundness. -/
theorem fractionalSandboxBadCanonicalizer_witness_valid :
    fractionalSandboxBadCanonicalizer
        (fractionalSandboxBadCanonicalizer .a0) ≠
      fractionalSandboxBadCanonicalizer .a0 :=
  canonicalizerWitness_sound fractionalSandboxBadCanonicalizer_rejected

/-- The constructed canonicalizer really is not idempotent. -/
theorem fractionalSandboxBadCanonicalizer_not_idempotent :
    ¬ CanonicalizerIdempotent fractionalSandboxBadCanonicalizer :=
  fun h => fractionalSandboxBadCanonicalizer_witness_valid (h .a0)

/-- The Boolean checker rejects the constructed canonicalizer. -/
theorem fractionalSandboxBadCanonicalizer_check_false :
    checkCanonicalizerIdempotent fractionalSandboxBadCanonicalizer
      fractionalSandboxStates = false := by
  decide

/-- `CANONICALIZER_IDEMPOTENCE` verdict for the exact sandbox data: accept. -/
theorem fractionalSandboxCanonicalizer_accepted :
    checkCanonicalizerIdempotent fractionalSandboxNormalForm
      fractionalSandboxStates = true := by
  decide

/-- The accepted sandbox normal-form map is mathematically idempotent. -/
theorem fractionalSandboxCanonicalizer_idempotent :
    CanonicalizerIdempotent fractionalSandboxNormalForm :=
  (checkCanonicalizerIdempotent_iff fractionalSandboxNormalForm
    fractionalSandboxStates_complete).mp
      fractionalSandboxCanonicalizer_accepted

/-! ### `REPRESENTATIVE_INVARIANCE`

The simulator's generic predicate requires an observable to be constant on
each canonical sector (`oph_fractional/quotient.py:50-61`).  The concrete
sandbox call supplies `lambda state: schema.canonicalize(state)`
(`oph_fractional/compare.py:133`), so the observable certified here is exactly
`fractionalSandboxQuotient`.  This receipt is distinct from auxiliary-section
independence of `quotientKernel`, already proved above. -/

/-- A deliberately representative-dependent observable. -/
@[reducible] def fractionalSandboxBadObservable :
    FractionalSandboxState → FractionalSandboxState :=
  id

/-- Negative control: `a0` and `a1` share a sector but have different
observable values. -/
theorem fractionalSandboxBadObservable_rejected :
    representativeWitness fractionalSandboxQuotient
      fractionalSandboxBadObservable fractionalSandboxStates =
        some (.a0, .a1) := by
  decide

/-- The returned representative pair is valid, through checker soundness. -/
theorem fractionalSandboxBadObservable_witness_valid :
    fractionalSandboxQuotient .a0 = fractionalSandboxQuotient .a1 ∧
      fractionalSandboxBadObservable .a0 ≠
        fractionalSandboxBadObservable .a1 :=
  representativeWitness_sound fractionalSandboxBadObservable_rejected

/-- The constructed observable really is not representative-invariant. -/
theorem fractionalSandboxBadObservable_not_invariant :
    ¬ RepresentativeInvariant fractionalSandboxQuotient
      fractionalSandboxBadObservable :=
  fun h => fractionalSandboxBadObservable_witness_valid.2
    (h .a0 .a1 fractionalSandboxBadObservable_witness_valid.1)

/-- The Boolean checker rejects the constructed observable. -/
theorem fractionalSandboxBadObservable_check_false :
    checkRepresentativeInvariant fractionalSandboxQuotient
      fractionalSandboxBadObservable fractionalSandboxStates = false := by
  decide

/-- `REPRESENTATIVE_INVARIANCE` verdict for the exact sandbox call: accept. -/
theorem fractionalSandboxRepresentativeInvariant_accepted :
    checkRepresentativeInvariant fractionalSandboxQuotient
      fractionalSandboxQuotient fractionalSandboxStates = true := by
  decide

/-- The sandbox canonical-sector observable is constant on quotient fibres. -/
theorem fractionalSandboxRepresentativeInvariant :
    RepresentativeInvariant fractionalSandboxQuotient
      fractionalSandboxQuotient :=
  (checkRepresentativeInvariant_iff fractionalSandboxQuotient
    fractionalSandboxQuotient fractionalSandboxStates_complete).mp
      fractionalSandboxRepresentativeInvariant_accepted

/-! ### `NO_ORBIT_SIZE_BIAS`

The simulator rejects nonpositive orbit sizes and, when sizes vary, rejects
sector weights proportional to those hidden sizes
(`oph_fractional/quotient.py:89-118`).  Its sandbox supplies orbit sizes
`(1, 1)` and sector weights `(1.0, 1.0)`
(`oph_fractional/compare.py:54-61,135`).  Here those values are exact integers
and rationals; the numerical tolerance in the floating-point implementation
is replaced by exact equality and strict positivity. -/

@[reducible] def fractionalSandboxOrbitSize :
    FractionalSandboxSector → ℤ
  | .anyonEOver3 => 1
  | .vacuum => 1

@[reducible] def fractionalSandboxSectorWeight :
    FractionalSandboxSector → ℚ
  | .anyonEOver3 => 1
  | .vacuum => 1

/-- Proportional-to-hidden-count negative control from the simulator test:
orbit sizes `(2, 1)` and sector weights `(2, 1)`. -/
@[reducible] def fractionalSandboxBiasedOrbitSize :
    FractionalSandboxSector → ℤ
  | .anyonEOver3 => 2
  | .vacuum => 1

@[reducible] def fractionalSandboxBiasedSectorWeight :
    FractionalSandboxSector → ℚ
  | .anyonEOver3 => 2
  | .vacuum => 1

/-- Negative control: the checker returns the varying-size pair and the two
positive sectors whose per-representative weights are both exactly `1`. -/
theorem fractionalSandboxBiasedOrbit_rejected :
    orbitSizeBiasWitness fractionalSandboxBiasedOrbitSize
      fractionalSandboxBiasedSectorWeight fractionalSandboxSectors =
        some (.tracksHiddenCount .anyonEOver3 .vacuum .anyonEOver3 .vacuum) := by
  norm_num [orbitSizeBiasWitness, hiddenOrbitCountWitness, List.find?,
    perRepresentativeWeight, fractionalSandboxBiasedOrbitSize,
    fractionalSandboxBiasedSectorWeight, fractionalSandboxSectors]
  simp

/-- The hidden-count witness is valid, through checker soundness. -/
theorem fractionalSandboxBiasedOrbit_witness_valid :
    ValidOrbitSizeBiasWitness fractionalSandboxBiasedOrbitSize
      fractionalSandboxBiasedSectorWeight
        (.tracksHiddenCount .anyonEOver3 .vacuum .anyonEOver3 .vacuum) :=
  orbitSizeBiasWitness_sound fractionalSandboxSectors_complete
    fractionalSandboxBiasedOrbit_rejected

/-- The proportional-weight control fails the mathematical property. -/
theorem fractionalSandboxBiasedOrbit_not_free :
    ¬ NoOrbitSizeBias fractionalSandboxBiasedOrbitSize
      fractionalSandboxBiasedSectorWeight :=
  validOrbitSizeBiasWitness_not_free
    fractionalSandboxBiasedOrbit_witness_valid

/-- The Boolean checker rejects the proportional-weight control. -/
theorem fractionalSandboxBiasedOrbit_check_false :
    checkNoOrbitSizeBias fractionalSandboxBiasedOrbitSize
      fractionalSandboxBiasedSectorWeight fractionalSandboxSectors = false := by
  simp [checkNoOrbitSizeBias, fractionalSandboxBiasedOrbit_rejected]

/-- A second negative route exercises malformed nonpositive orbit sizes. -/
@[reducible] def fractionalSandboxNonpositiveOrbitSize :
    FractionalSandboxSector → ℤ
  | .anyonEOver3 => 0
  | .vacuum => 1

theorem fractionalSandboxNonpositiveOrbit_rejected :
    orbitSizeBiasWitness fractionalSandboxNonpositiveOrbitSize
      fractionalSandboxSectorWeight fractionalSandboxSectors =
        some (.nonpositive .anyonEOver3) := by
  decide

theorem fractionalSandboxNonpositiveOrbit_witness_valid :
    ValidOrbitSizeBiasWitness fractionalSandboxNonpositiveOrbitSize
      fractionalSandboxSectorWeight (.nonpositive .anyonEOver3) :=
  orbitSizeBiasWitness_sound fractionalSandboxSectors_complete
    fractionalSandboxNonpositiveOrbit_rejected

/-- `NO_ORBIT_SIZE_BIAS` verdict for the exact sandbox data: accept. -/
theorem fractionalSandboxNoOrbitSizeBias_accepted :
    checkNoOrbitSizeBias fractionalSandboxOrbitSize
      fractionalSandboxSectorWeight fractionalSandboxSectors = true := by
  decide

/-- The accepted verdict implies the exact no-bias mathematical property. -/
theorem fractionalSandboxNoOrbitSizeBias :
    NoOrbitSizeBias fractionalSandboxOrbitSize
      fractionalSandboxSectorWeight :=
  (checkNoOrbitSizeBias_iff fractionalSandboxOrbitSize
    fractionalSandboxSectorWeight fractionalSandboxSectors_complete).mp
      fractionalSandboxNoOrbitSizeBias_accepted

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
#print axioms canonicalizerWitness_sound
#print axioms canonicalizerWitness_none
#print axioms checkCanonicalizerIdempotent_iff
#print axioms representativeWitness_sound
#print axioms representativeWitness_none
#print axioms checkRepresentativeInvariant_iff
#print axioms hiddenOrbitCountWitness_sound
#print axioms hiddenOrbitCountWitness_none
#print axioms orbitSizeBiasWitness_sound
#print axioms validOrbitSizeBiasWitness_not_free
#print axioms orbitSizeBiasWitness_none
#print axioms checkNoOrbitSizeBias_iff
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
#print axioms fractionalSandboxBadCanonicalizer_rejected
#print axioms fractionalSandboxBadCanonicalizer_witness_valid
#print axioms fractionalSandboxBadCanonicalizer_not_idempotent
#print axioms fractionalSandboxBadCanonicalizer_check_false
#print axioms fractionalSandboxCanonicalizer_accepted
#print axioms fractionalSandboxCanonicalizer_idempotent
#print axioms fractionalSandboxBadObservable_rejected
#print axioms fractionalSandboxBadObservable_witness_valid
#print axioms fractionalSandboxBadObservable_not_invariant
#print axioms fractionalSandboxBadObservable_check_false
#print axioms fractionalSandboxRepresentativeInvariant_accepted
#print axioms fractionalSandboxRepresentativeInvariant
#print axioms fractionalSandboxBiasedOrbit_rejected
#print axioms fractionalSandboxBiasedOrbit_witness_valid
#print axioms fractionalSandboxBiasedOrbit_not_free
#print axioms fractionalSandboxBiasedOrbit_check_false
#print axioms fractionalSandboxNonpositiveOrbit_rejected
#print axioms fractionalSandboxNonpositiveOrbit_witness_valid
#print axioms fractionalSandboxNoOrbitSizeBias_accepted
#print axioms fractionalSandboxNoOrbitSizeBias
#print axioms fractionalSandbox_accepted
#print axioms fractionalSandbox_lumpable
#print axioms fractionalSandbox_isMarkov
#print axioms fractionalSandbox_quotientKernel

end OPH.QuotientLumpability
