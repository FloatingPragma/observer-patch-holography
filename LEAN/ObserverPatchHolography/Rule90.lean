import Mathlib
import ObserverPatchHolography.Primitives

/-!
# A cellular-automaton OPH carrier: Rule 90 (issue #304, `Hfib` half)

This module builds a **non-degenerate** OPH carrier out of an additive
cellular automaton — Wolfram's **Rule 90** (`cell' = left XOR right`) — and
discharges the boundary-fiber singleton hypothesis `Hfib` on it for a
*proper information-set* boundary, together with the matching **failure
witness** and a **non-trivial gauge** exhibit.

## Why Rule 90, and why it is not a degenerate demo

`demoCarrier` (in `Primitives.lean`) has a single edge whose only consistent
records are globally constant, so its only `Hfib` witnesses are the trivial
top (`obsMap` itself, `demoCarrier_Hfib_holds_finerB`) and the seed cell
(`demoCarrier_Hfib_holds_seed`) — where the seed already carries the whole
observable. There is no room for a boundary that is *strictly between* "reads
nothing" and "reads everything" yet still identifies the fiber, and none that
can *fail*.

Rule 90 supplies exactly that room. Encode one CA time-step as a two-patch
carrier: patch `false` holds the seed row, patch `true` the next row, and the
single edge's `src` projection applies the Rule-90 map while `tgt` is the
identity, so **edge-consistency ⇔ "the bottom row is the Rule-90 image of the
seed" ⇔ a valid CA diagram**. Because Rule 90 is *linear over `𝔽₂`*, the
consistent set is a genuine linear code and the boundary question becomes the
classical **erasure-correction / information-set** question. On width 3 with
zero boundary the CA map sends a seed `(a,b,c)` to the row `(b, a⊕c, b)`,
whose outer cells are **equal** — a built-in redundancy. Hence:

* `rule90_Hfib_good` — reading bottom-row cells `{0,1}` (`rule90BoundaryGood`)
  **identifies** the fiber: `Hfib` holds. This is a proper information set —
  coarser than the full observable, and it genuinely consumes consistency
  (the unread outer cell is forced to follow cell 0 via the CA redundancy).
* `rule90_Hfib_bad_fails` — reading cells `{0,2}` (`rule90BoundaryBad`) reads
  the *same* bit twice (`b`,`b`), misses the middle `a⊕c`, and **fails**
  `Hfib`: two consistent diagrams agree on the boundary yet differ observably.
* `rule90_gauge_nontrivial` — the CA map is *non-injective* (its `𝔽₂`-kernel
  is 1-dimensional: `(a,b,c)` and `(a⊕1,b,c⊕1)` share an image), so two
  records with **different seeds** are `gaugeEquiv`. The gauge is real: it
  contains the CA map's kernel (a kernel-pair exhibit — that inclusion
  direction is what is machine-checked here).
* `rule90_no_frustrationFree_repair` — **no** frustration-free local repair
  (the `H1 ∧ H2 ∧ H3` binder forms of `Primitives.lean`'s
  `LocalRepairDynamics` section) exists on this carrier: a bottom row with
  unequal outer cells lies outside the Rule-90 image for *every* seed, yet
  `H2` forces the seed-patch move to fire and `H1` pins the bottom row, so
  `H3` is unsatisfiable.

## Honest scope

This is the **`Hfib` half** of #304 on a real CA carrier — strictly richer
than `demoCarrier`'s seed tautology (a proper, failable information set + a
non-trivial gauge). It does **not** yet supply the `HB` (repair-preserved
boundary) premise, so it does not by itself instantiate the *joint* `HB ∧
Hfib` witness that `boundary_fiber_observer_unique`'s doc-comment flags as
future work. Indeed it never can: `rule90_no_frustrationFree_repair` below
proves that **no** frustration-free local repair (`H1 ∧ H2 ∧ H3`) exists on
this carrier at all — the repair route is provably closed here, and the
injectivity reading of `Hfib` (#304) is the one that survives. No
`sorry`, no new axioms (the theorems reduce pointwise; no `native_decide`).
-/

namespace OPH

/-- One step of **Rule 90** on a width-3 tape with zero boundary conditions,
    a triple `(a,b,c) = (cell₀,cell₁,cell₂)`. With the outside cells fixed to
    `0`, `cell'ᵢ = cellᵢ₋₁ ⊕ cellᵢ₊₁` gives
    `(a,b,c) ↦ (0⊕b, a⊕c, b⊕0) = (b, a⊕c, b)`. The outer cells of the image
    coincide (`= b`); that redundancy is what makes `{0,1}` an information set
    and `{0,2}` deficient. -/
def rule90t (t : Bool × Bool × Bool) : Bool × Bool × Bool :=
  (t.2.1, xor t.1 t.2.2, t.2.1)

/-- The **Rule-90 carrier**: two patches (seed row `false`, next row `true`),
    one edge, interface = a width-3 row (`Bool × Bool × Bool`). The `src`
    projection applies the CA map; the `tgt` projection is the identity — so
    edge-consistency says exactly `next row = Rule90(seed)`, i.e. the record
    is a valid one-step CA diagram. Unit weight, discrete row metric. -/
def rule90Carrier : OPHCarrier where
  Patch := Bool
  State := fun _ => Bool × Bool × Bool
  Edge := Unit
  src := fun _ => false
  tgt := fun _ => true
  Iface := fun _ => Bool × Bool × Bool
  projSrc := fun _ s => rule90t s
  projTgt := fun _ s => s
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases h : a = b
    · rw [if_pos h]; exact ⟨fun _ => h, fun _ => rfl⟩
    · rw [if_neg h]; exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 h⟩

/-- **The information-set boundary.** Reads bottom-row cells `{0,1}` — cell 0
    and the middle cell. Strictly coarser than the full observable `obsMap`. -/
def rule90BoundaryGood : Records rule90Carrier → Bool × Bool :=
  fun x => ((x true).1, (x true).2.1)

/-- **The deficient boundary.** Reads bottom-row cells `{0,2}` — the two
    *outer* cells, which the CA redundancy forces equal, so this reads one bit
    twice and never sees the middle `a⊕c`. -/
def rule90BoundaryBad : Records rule90Carrier → Bool × Bool :=
  fun x => ((x true).1, (x true).2.2)

/-- **`Hfib` HOLDS for the information-set boundary `{0,1}` (issue #304).**
    Any two consistent CA diagrams with equal `rule90BoundaryGood` are
    `gaugeEquiv`. The proof genuinely uses consistency: the unread outer cell
    (cell 2) is not in the boundary, but edge-agreement forces it to equal
    cell 0 (the Rule-90 redundancy `image = (b, a⊕c, b)`), so the two read
    cells pin the whole bottom row and hence the observable. -/
theorem rule90_Hfib_good :
    ∀ x y : Records rule90Carrier, rule90BoundaryGood x = rule90BoundaryGood y →
      Consistent rule90Carrier x → Consistent rule90Carrier y → gaugeEquiv rule90Carrier x y := by
  intro x y hB hcx hcy
  rw [consistent_iff_edgeConsistent] at hcx hcy
  -- Consistency: the bottom row IS the Rule-90 image of the seed.
  have hx : rule90t (x false) = x true := hcx ()
  have hy : rule90t (y false) = y true := hcy ()
  -- Boundary gives cell 0 and cell 1 of the bottom row.
  simp only [rule90BoundaryGood, Prod.mk.injEq] at hB
  obtain ⟨hB0, hB1⟩ := hB
  -- CA redundancy: bottom cell 2 = bottom cell 0 (both are seed cell 1).
  have ex : (x true).2.2 = (x true).1 := by rw [← hx]; rfl
  have ey : (y true).2.2 = (y true).1 := by rw [← hy]; rfl
  have hB2 : (x true).2.2 = (y true).2.2 := by rw [ex, ey, hB0]
  -- Hence the whole bottom rows agree.
  have hxtrue : x true = y true :=
    Prod.ext_iff.mpr ⟨hB0, Prod.ext_iff.mpr ⟨hB1, hB2⟩⟩
  -- Observable = (Rule90(seed), bottom row) = (bottom row, bottom row) once consistent.
  show obsMap rule90Carrier x = obsMap rule90Carrier y
  funext e
  cases e
  show (rule90t (x false), x true) = (rule90t (y false), y true)
  rw [hx, hy, hxtrue]

/-- **`Hfib` FAILS for the deficient boundary `{0,2}` — explicit witness.**
    Seed `(0,0,0)` gives bottom `(0,0,0)`; seed `(0,0,1)` gives bottom
    `(0,1,0)`. Both are consistent, both read `(0,0)` on cells `{0,2}`, yet
    their observables differ (`(0,0,0) ≠ (0,1,0)` on the second component), so
    they are not `gaugeEquiv`. The boundary is not an information set: it
    erases the only non-redundant cell. -/
theorem rule90_Hfib_bad_fails :
    ∃ x y : Records rule90Carrier,
      rule90BoundaryBad x = rule90BoundaryBad y ∧
      Consistent rule90Carrier x ∧ Consistent rule90Carrier y ∧
      ¬ gaugeEquiv rule90Carrier x y := by
  refine ⟨fun _ => (false, false, false),
          fun b => bif b then (false, true, false) else (false, false, true), ?_, ?_, ?_, ?_⟩
  · rfl
  · rw [consistent_iff_edgeConsistent]; intro e; cases e; rfl
  · rw [consistent_iff_edgeConsistent]; intro e; cases e; rfl
  · intro hg
    have h2 : (fun _ => (false, false, false) : Records rule90Carrier) true
            = (fun b => bif b then (false, true, false) else (false, false, true) : Records rule90Carrier) true :=
      congrArg Prod.snd (congrFun hg ())
    exact absurd h2 (by decide)

/-- **The gauge is non-trivial — it contains the CA map's kernel.** Rule 90 on
    width 3 is non-injective: seeds `(0,0,0)` and `(1,0,1)` both map to the
    zero row (`a⊕c = 0`, `b = 0`). So two consistent records with **different
    seeds** expose the same observable and are `gaugeEquiv`. The unobservable
    part of the seed — the erasure that `Hfib` correctly quotients away —
    contains `ker(Rule90)` (a kernel-pair exhibit: that inclusion direction
    is what is machine-checked here). -/
theorem rule90_gauge_nontrivial :
    ∃ x y : Records rule90Carrier,
      x false ≠ y false ∧
      Consistent rule90Carrier x ∧ Consistent rule90Carrier y ∧
      gaugeEquiv rule90Carrier x y := by
  refine ⟨fun _ => (false, false, false),
          fun b => bif b then (false, false, false) else (true, false, true), ?_, ?_, ?_, ?_⟩
  · show ((false, false, false) : Bool × Bool × Bool) ≠ (true, false, true); decide
  · rw [consistent_iff_edgeConsistent]; intro e; cases e; rfl
  · rw [consistent_iff_edgeConsistent]; intro e; cases e; rfl
  · show obsMap rule90Carrier _ = obsMap rule90Carrier _
    funext e; cases e; rfl

/-- Both outer cells of any Rule-90 image row coincide — each equals the middle
    seed cell — so a bottom row with **unequal** outer cells lies outside the
    image of `rule90t` for *every* seed. -/
theorem rule90t_outer_eq (s : Bool × Bool × Bool) :
    (rule90t s).1 = (rule90t s).2.2 := rfl

/-- **No frustration-free local repair exists on this carrier.** There is no
    local move `lr` satisfying the `LocalRepairDynamics` hypotheses of
    `Primitives.lean` (binder forms verbatim): `H1` — firing at `i` changes
    patch `i` only; `H2` — the move at `i` fires **iff** some edge incident to
    `i` is inconsistent; `H3` — after firing at `i`, all edges incident to `i`
    are consistent. Reason: the record whose bottom row is `(0,0,1)` has
    unequal outer cells, hence its only edge is broken for **every** seed row
    (`rule90t_outer_eq`). `H2` then forces the *seed*-patch move to fire, `H1`
    pins the bottom row, and `H3` demands a Rule-90 preimage of `(0,0,1)` —
    which cannot exist. The repair route is closed on this carrier; the
    injectivity reading of `Hfib` (#304) is the one that survives. -/
theorem rule90_no_frustrationFree_repair :
    ¬ ∃ lr : rule90Carrier.Patch → Records rule90Carrier → Records rule90Carrier,
      (∀ (i : rule90Carrier.Patch) (x : Records rule90Carrier) (j : rule90Carrier.Patch),
          j ≠ i → (lr i x) j = x j) ∧
      (∀ (i : rule90Carrier.Patch) (x : Records rule90Carrier),
          lr i x ≠ x ↔
            ∃ e : rule90Carrier.Edge,
              (rule90Carrier.src e = i ∨ rule90Carrier.tgt e = i) ∧
                ¬ edgeConsistentAt e x) ∧
      (∀ (i : rule90Carrier.Patch) (x : Records rule90Carrier),
          lr i x ≠ x →
            ∀ e : rule90Carrier.Edge,
              (rule90Carrier.src e = i ∨ rule90Carrier.tgt e = i) →
                edgeConsistentAt e (lr i x)) := by
  rintro ⟨lr, H1, H2, H3⟩
  -- It suffices to refute any record whose bottom row is the out-of-image (0,0,1).
  suffices main : ∀ x : Records rule90Carrier, x true = (false, false, true) → False by
    exact main (fun b => bif b then (false, false, true) else (false, false, false)) rfl
  intro x hxt
  -- The single edge is broken at `x`, whatever the seed row is.
  have hbroken : ¬ edgeConsistentAt (C := rule90Carrier) () x := by
    intro h
    have h' : rule90t (x false) = x true := h
    have h02 := rule90t_outer_eq (x false)
    rw [h', hxt] at h02
    exact absurd h02 (by decide)
  -- H2 (⇐): a broken edge incident to the seed patch forces its move to fire.
  have hfire : lr false x ≠ x := (H2 false x).mpr ⟨(), Or.inl rfl, hbroken⟩
  -- H3: after firing at the seed patch, the edge must be consistent …
  have hcons : rule90t ((lr false x) false) = (lr false x) true :=
    H3 false x hfire () (Or.inl rfl)
  -- … while H1 keeps the bottom row untouched (`true ≠ false`).
  have htgt : (lr false x) true = (false, false, true) :=
    (H1 false x true (fun h => Bool.noConfusion h)).trans hxt
  -- So (0,0,1) would be a Rule-90 image — its outer cells would coincide.
  have h02 := rule90t_outer_eq ((lr false x) false)
  rw [hcons, htgt] at h02
  exact absurd h02 (by decide)

-- Axiom audit: these must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms rule90_Hfib_good
#print axioms rule90_Hfib_bad_fails
#print axioms rule90_gauge_nontrivial
#print axioms rule90_no_frustrationFree_repair

end OPH
