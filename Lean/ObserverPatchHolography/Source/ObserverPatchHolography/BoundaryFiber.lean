import Mathlib
import ObserverPatchHolography.Primitives
import ObservableNormalForms.ObserverConfluence

/-!
# The declared boundary map `B_OPH` and its injectivity modulo gauge (issue #304)

This module discharges the concrete application gate of issue #304 on the
declared consistent-state domain. It supplies the three objects the gate
names, all on the verified rooted-tree packet-net domain of *Reality as a
Consensus Protocol* (Definition `def:tree-packet-domain`):

1. **The domain.** `TreePacketNet` encodes the declared fixed-cutoff patch
   class verbatim: a finite rooted tree of patches, per-patch state
   `S_i = A × K_i` (public packet times hidden gauge label), every interface
   projection reading the packet component only. `TreePacketNet.carrier`
   realises it as an `OPHCarrier`, so `Records`, `Consistent`, `obsMap`, and
   `gaugeEquiv` are the repository's own primitives, not fresh look-alikes.
   The hidden labels `K_i` are the declared redundancy data: they are read by
   no interface projection, so relabeling them moves inside one `gaugeEquiv`
   class (`fourVertexNet_gauge_nontrivial` exhibits this concretely).

2. **The map.** `TreePacketNet.BOPH` is the declared physical boundary/sector
   map: the root-packet readback `x ↦ (x root).1`. This is the "root-packet
   map" entry of the paper's protected-data taxonomy
   (Definition `def:finite-quotient-repair-presentation`), defined
   concretely, not left as a binder.

3. **The theorem.** `TreePacketNet.BOPH_injective_modulo_gauge` proves the
   gate statement in the issue's own shape:
   `u, v ∈ C_OPH` and `B_OPH(u) = B_OPH(v)` imply `u ~gauge v`,
   for every net in the class (every finite rooted tree, every packet
   alphabet, every hidden-label family). The proof genuinely consumes
   consistency: edge agreement propagates the root packet down the tree by
   strong induction on depth, so the unread bulk packets are forced while the
   hidden labels stay free. The conclusion is `gaugeEquiv`, not record
   equality; `fourVertexNet_gauge_nontrivial` shows the quotient is proper.

## Connection to the generic identifiability layer

`TreePacketNet.BOPH_boundaryIdentifiesModulo` restates the theorem as the
`BoundaryIdentifiesModulo` premise of the substrate-neutral package under
`Proofs/ObservableNormalForms/`. The two generic hypotheses of the neutral
endpoint theorem are then discharged for the domain's own repair law
(`treeRepairAt`, the paper's `T_i`):

* `treeStep_observationPreserving` — accepted tree repairs never write the
  root, so `B_OPH` is preserved (hypothesis `ObservationPreserving`);
* `treeStep_completeFor` — a record is a `treeStep` normal form iff it is
  `Consistent` (hypothesis `CompleteFor`).

`BOPH_observerEndpointUnique` composes the three through
`boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo`: any two records
with equal root-packet readback settle, under any maximal tree-repair
schedule, to gauge-equal endpoints. This is the physical payoff form of #304
on the declared domain, with no confluence input anywhere.

## Scope

The machine-checked theorem covers the declared verified domain (packet-copy
constraints along a rooted tree). The paper's Theorem
`thm:functional-selected-fiber` extends the same propagation argument to
general single-parent deterministic extension maps; that generalization is
paper-level, with the packet-copy rule as its realized Lean instance here.
The theorem does not assert
identifiability for an arbitrary finite patch net with an arbitrary boundary
map: that statement is false, with machine-checked countermodels in the tree
(`demoCarrier_Hfib_fails`, `rule90_Hfib_bad_fails`). For any claim
surface whose patch net leaves the declared class, boundary identifiability
stays a named per-net premise. Edge weights are set to `1` here because the
identifiability statement never reads them; the paper's weight ordering
matters for repair descent, not for the fiber question. No `sorry`, no new
axioms; the audit block at the end must report only standard axioms.
-/

namespace OPH

open Relation

/-- The verified rooted-tree packet-net domain (*Reality as a Consensus
    Protocol*, Definition `def:tree-packet-domain`): a finite rooted tree of
    patches with per-patch state `A × K i`. `parent` and `depth` encode the
    tree; `parent_progress` forces every non-root ancestor chain to descend
    to the root, so the type carries exactly the rooted-tree nets. -/
structure TreePacketNet where
  /-- Patch vertices of the tree. -/
  V : Type
  /-- The vertex set is finite. -/
  [vFintype : Fintype V]
  /-- Vertices have decidable equality. -/
  [vDecEq : DecidableEq V]
  /-- The public packet alphabet `A`, `|A| ≥ 1` suffices for the statement. -/
  A : Type
  /-- Packets have decidable equality (used by the discrete edge metric). -/
  [aDecEq : DecidableEq A]
  /-- Per-patch hidden gauge labels `K_i` (declared redundancy data). -/
  K : V → Type
  /-- The root patch (the screen port). -/
  root : V
  /-- Parent map of the rooted tree (value at `root` is irrelevant). -/
  parent : V → V
  /-- Tree depth used as the descent measure of the parent map. -/
  depth : V → ℕ
  /-- Every non-root vertex has a strictly shallower parent. -/
  parent_progress : ∀ i : V, i ≠ root → depth (parent i) < depth i

attribute [instance] TreePacketNet.vFintype TreePacketNet.vDecEq TreePacketNet.aDecEq

namespace TreePacketNet

variable (T : TreePacketNet)

/-- The rooted-tree packet net as an `OPHCarrier`: one edge `{parent i, i}`
    per non-root vertex `i`, interface alphabet `A`, both projections reading
    the packet component only (the hidden labels `K i` are invisible to every
    interface), unit weights, discrete edge metric. -/
def carrier : OPHCarrier where
  Patch := T.V
  State := fun i => T.A × T.K i
  Edge := {i : T.V // i ≠ T.root}
  src := fun e => T.parent e.val
  tgt := fun e => e.val
  Iface := fun _ => T.A
  projSrc := fun _ s => s.1
  projTgt := fun _ s => s.1
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases h : a = b
    · rw [if_pos h]; exact ⟨fun _ => h, fun _ => rfl⟩
    · rw [if_neg h]; exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 h⟩

/-- **The declared physical boundary/sector map `B_OPH`**: the root-packet
    readback `x ↦ (x root).1`. This is the "root-packet map" entry of the
    protected-data taxonomy in Definition
    `def:finite-quotient-repair-presentation`, made concrete on the verified
    domain. It reads the public packet of the screen port and nothing else:
    no bulk packet and no hidden label. -/
def BOPH (x : Records T.carrier) : T.A := (x T.root).1

/-- Unique continuation on the declared domain: on an edge-consistent record
    every packet equals the root packet. Strong induction on `depth`; the
    parent edge of each non-root vertex forces its packet to follow its
    parent, and `parent_progress` drives the chain to the root. -/
theorem packet_eq_root_of_edgeConsistent
    {x : Records T.carrier} (hx : EdgeConsistent T.carrier x) :
    ∀ i : T.V, (x i).1 = (x T.root).1 := by
  have main : ∀ (n : ℕ) (i : T.V), T.depth i ≤ n → (x i).1 = (x T.root).1 := by
    intro n
    induction n with
    | zero =>
        intro i hle
        by_cases hroot : i = T.root
        · rw [hroot]
        · exact absurd (lt_of_lt_of_le (T.parent_progress i hroot) hle) (Nat.not_lt_zero _)
    | succ n ih =>
        intro i hle
        by_cases hroot : i = T.root
        · rw [hroot]
        · have hedge : (x (T.parent i)).1 = (x i).1 := hx ⟨i, hroot⟩
          have hp : T.depth (T.parent i) ≤ n :=
            Nat.lt_succ_iff.mp (lt_of_lt_of_le (T.parent_progress i hroot) hle)
          exact hedge.symm.trans (ih (T.parent i) hp)
  exact fun i => main (T.depth i) i le_rfl

/-- **THEOREM — the #304 application gate on the declared domain.** In the
    issue's own shape: `u, v ∈ C_OPH` and `B_OPH(u) = B_OPH(v)` imply
    `u ~gauge v`, for every rooted-tree packet net. Consistency propagates
    the root packet through the bulk (`packet_eq_root_of_edgeConsistent`), so
    equal root readback forces equal packets on every interface; the residual
    freedom is exactly the hidden labels, which lie inside `gaugeEquiv`.
    The conclusion is observable equality, not record equality:
    `fourVertexNet_gauge_nontrivial` below shows the gauge quotient is
    proper, so `B_OPH` is injective on `C / ~gauge` and not on `C`. -/
theorem BOPH_injective_modulo_gauge
    {u v : Records T.carrier}
    (hu : Consistent T.carrier u) (hv : Consistent T.carrier v)
    (hB : T.BOPH u = T.BOPH v) :
    gaugeEquiv T.carrier u v := by
  rw [consistent_iff_edgeConsistent] at hu hv
  have hu' := T.packet_eq_root_of_edgeConsistent hu
  have hv' := T.packet_eq_root_of_edgeConsistent hv
  show obsMap T.carrier u = obsMap T.carrier v
  funext e
  obtain ⟨i, hne⟩ := e
  show ((u (T.parent i)).1, (u i).1) = ((v (T.parent i)).1, (v i).1)
  rw [hu' (T.parent i), hu' i, hv' (T.parent i), hv' i]
  show (T.BOPH u, T.BOPH u) = (T.BOPH v, T.BOPH v)
  rw [hB]

/-- The gate restated in the generic vocabulary of the substrate-neutral
    package: `B_OPH` satisfies `BoundaryIdentifiesModulo` on the consistent
    set, modulo `gaugeEquiv`. This is the exact premise the neutral endpoint
    theorem leaves open, discharged for the declared domain. -/
theorem BOPH_boundaryIdentifiesModulo :
    ObservableNormalForms.BoundaryIdentifiesModulo
      {x : Records T.carrier | Consistent T.carrier x} T.BOPH (gaugeEquiv T.carrier) := by
  intro u v hu hv hB
  exact T.BOPH_injective_modulo_gauge hu hv hB

/-! ## The domain's own repair law and the endpoint payoff

The paper's repair `T_i` sets the packet of a non-root vertex to its parent's
packet and leaves the hidden label alone. The two generic hypotheses of the
neutral endpoint theorem hold for it by direct computation, so the physical
payoff form of #304 (same boundary readback, gauge-equal endpoints under any
maximal schedule) follows on the declared domain with no confluence input. -/

/-- The paper's tree repair `T_i`: copy the parent's packet into vertex `i`,
    keep the hidden label. The root move is the identity (the guard in
    `treeStep` fires only at non-root vertices). -/
def treeRepairAt (i : T.V) (x : Records T.carrier) : Records T.carrier :=
  Function.update x i ((x (T.parent i)).1, (x i).2)

/-- One accepted asynchronous tree-repair step: some non-root vertex's copy
    move changes the record. -/
def treeStep (x y : Records T.carrier) : Prop :=
  ∃ i : T.V, i ≠ T.root ∧ y = T.treeRepairAt i x ∧ T.treeRepairAt i x ≠ x

/-- Accepted tree repairs never write the root, so the declared boundary map
    is preserved: the `ObservationPreserving` hypothesis of the neutral
    endpoint theorem, for `B_OPH`. -/
theorem treeStep_observationPreserving :
    ObservableNormalForms.ObservationPreserving T.treeStep T.BOPH := by
  intro x y hxy
  obtain ⟨i, hne, rfl, _⟩ := hxy
  show (x T.root).1 = ((T.treeRepairAt i x) T.root).1
  have hupd : (T.treeRepairAt i x) T.root = x T.root :=
    Function.update_of_ne (Ne.symm hne) _ _
  rw [hupd]

/-- A record is a `treeStep` normal form iff it is `Consistent`: the
    `CompleteFor` hypothesis of the neutral endpoint theorem. Every edge of
    the net is the parent edge of exactly one non-root vertex, and the copy
    move at `i` fires iff that parent edge is broken. -/
theorem treeStep_completeFor :
    ObservableNormalForms.CompleteFor T.treeStep
      {x : Records T.carrier | Consistent T.carrier x} := by
  intro x
  rw [Set.mem_setOf_eq, consistent_iff_edgeConsistent]
  constructor
  · intro hnf e
    obtain ⟨i, hne⟩ := e
    show (x (T.parent i)).1 = (x i).1
    by_contra hbad
    apply hnf (T.treeRepairAt i x)
    refine ⟨i, hne, rfl, fun heq => ?_⟩
    have hi := congrFun heq i
    rw [treeRepairAt, Function.update_self] at hi
    exact hbad (congrArg Prod.fst hi)
  · intro hcons y hstep
    obtain ⟨i, hne, rfl, hfire⟩ := hstep
    apply hfire
    have hedge : (x (T.parent i)).1 = (x i).1 := hcons ⟨i, hne⟩
    funext j
    rcases eq_or_ne j i with rfl | hji
    · rw [treeRepairAt, Function.update_self]
      exact Prod.ext_iff.mpr ⟨hedge, rfl⟩
    · rw [treeRepairAt, Function.update_of_ne hji]

/-- **Endpoint payoff on the declared domain.** Any two records with the same
    `B_OPH` readback settle, along any maximal accepted tree-repair
    schedules, to gauge-equal normal forms. Composition of the discharged
    `BoundaryIdentifiesModulo` premise with the neutral equivalence; the
    boundary does the work, confluence never enters. -/
theorem BOPH_observerEndpointUnique :
    ObservableNormalForms.ObserverEndpointUniqueModulo
      T.treeStep T.BOPH (gaugeEquiv T.carrier) :=
  (ObservableNormalForms.boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    T.treeStep_observationPreserving T.treeStep_completeFor).mp
    T.BOPH_boundaryIdentifiesModulo

end TreePacketNet

/-! ## The paper's exported four-vertex instance

The verified domain record cited by Theorem `thm:tree-packet-domain`: vertices
`r, a, b, c` with edges `r–a`, `a–b`, `a–c`, packet alphabet `ℤ₃`, hidden
labels `ℤ₂`. It witnesses that the class is non-degenerate (multi-patch,
multi-edge, branching) and that the gauge quotient is proper: `B_OPH` is
injective modulo gauge but genuinely not injective on records, and a readback
that misses the protected packet fails to identify. -/

/-- The four vertices of the exported instance. -/
inductive FourV : Type
  | r | a | b | c
  deriving DecidableEq, Fintype

/-- The exported four-vertex net: root `r`, edges `r–a`, `a–b`, `a–c`,
    packets `ℤ₃` (as `Fin 3`), hidden labels `ℤ₂` (as `Fin 2`). -/
def fourVertexNet : TreePacketNet where
  V := FourV
  A := Fin 3
  K := fun _ => Fin 2
  root := FourV.r
  parent := fun i => match i with
    | .r => .r
    | .a => .r
    | .b => .a
    | .c => .a
  depth := fun i => match i with
    | .r => 0
    | .a => 1
    | .b => 2
    | .c => 2
  parent_progress := by
    intro i h
    cases i
    · exact absurd rfl h
    · decide
    · decide
    · decide

/-- The record with packet `p` everywhere and hidden label `h` everywhere. -/
def fourVRecord (p : Fin 3) (h : Fin 2) : Records fourVertexNet.carrier :=
  fun _ => (p, h)

theorem fourVRecord_consistent (p : Fin 3) (h : Fin 2) :
    Consistent fourVertexNet.carrier (fourVRecord p h) := by
  rw [consistent_iff_edgeConsistent]
  intro _
  rfl

/-- **The gauge quotient is proper on the declared domain.** Two consistent
    records with the same `B_OPH` readback that differ as records (distinct
    hidden labels) yet are `gaugeEquiv`. So the conclusion of
    `BOPH_injective_modulo_gauge` cannot be strengthened to record equality:
    the theorem is genuinely a statement modulo the declared redundancy. -/
theorem fourVertexNet_gauge_nontrivial :
    ∃ u v : Records fourVertexNet.carrier,
      u ≠ v ∧
      fourVertexNet.BOPH u = fourVertexNet.BOPH v ∧
      Consistent fourVertexNet.carrier u ∧ Consistent fourVertexNet.carrier v ∧
      gaugeEquiv fourVertexNet.carrier u v := by
  refine ⟨fourVRecord 0 0, fourVRecord 0 1, ?_, rfl,
    fourVRecord_consistent 0 0, fourVRecord_consistent 0 1, ?_⟩
  · intro h
    have hr := congrFun h FourV.r
    have hsnd : (0 : Fin 2) = 1 := congrArg Prod.snd hr
    exact absurd hsnd (by decide)
  · show obsMap fourVertexNet.carrier _ = obsMap fourVertexNet.carrier _
    funext _
    rfl

/-- A deficient readback on the same instance: reading the root's hidden
    label instead of its packet. Protected data misdeclared as gauge and
    gauge data misdeclared as protected. -/
def fourVHiddenReadback (x : Records fourVertexNet.carrier) : Fin 2 :=
  (x FourV.r).2

/-- **Failure witness: a readback missing the protected packet does not
    identify.** Two consistent records with equal hidden-label readback whose
    observables differ. The identifiability of `BOPH_injective_modulo_gauge`
    is a property of the declared boundary map reading the protected packet,
    not of readback maps in general; the coarse-boundary countermodels
    elsewhere in the tree (`demoCarrier_Hfib_fails`, `rule90_Hfib_bad_fails`)
    make the same point on other carriers. -/
theorem fourVHiddenReadback_not_identifying :
    ∃ u v : Records fourVertexNet.carrier,
      fourVHiddenReadback u = fourVHiddenReadback v ∧
      Consistent fourVertexNet.carrier u ∧ Consistent fourVertexNet.carrier v ∧
      ¬ gaugeEquiv fourVertexNet.carrier u v := by
  refine ⟨fourVRecord 0 0, fourVRecord 1 0, rfl,
    fourVRecord_consistent 0 0, fourVRecord_consistent 1 0, ?_⟩
  intro hg
  have he := congrFun hg ⟨FourV.a, fun h => FourV.noConfusion h⟩
  have hfst : (0 : Fin 3) = 1 := congrArg Prod.fst he
  exact absurd hfst (by decide)

/-- The gate theorem instantiated on the exported four-vertex record, so the
    class-level statement has a concrete, fully evaluated member. -/
theorem fourVertexNet_BOPH_injective_modulo_gauge
    {u v : Records fourVertexNet.carrier}
    (hu : Consistent fourVertexNet.carrier u)
    (hv : Consistent fourVertexNet.carrier v)
    (hB : fourVertexNet.BOPH u = fourVertexNet.BOPH v) :
    gaugeEquiv fourVertexNet.carrier u v :=
  fourVertexNet.BOPH_injective_modulo_gauge hu hv hB

/-! ### Axiom audit

The gate theorem, its generic-vocabulary form, the endpoint payoff, and the
witnesses must depend only on the standard axioms (`propext`,
`Classical.choice`, `Quot.sound`), never on `sorryAx`: none of this module
touches the three declared admissions in `Primitives.lean`. -/
#print axioms TreePacketNet.BOPH_injective_modulo_gauge
#print axioms TreePacketNet.BOPH_boundaryIdentifiesModulo
#print axioms TreePacketNet.treeStep_observationPreserving
#print axioms TreePacketNet.treeStep_completeFor
#print axioms TreePacketNet.BOPH_observerEndpointUnique
#print axioms fourVertexNet_gauge_nontrivial
#print axioms fourVHiddenReadback_not_identifying
#print axioms fourVertexNet_BOPH_injective_modulo_gauge

end OPH
