import OphGap.Checker

/-!
# Transport: from the BFS-relabelled certificate back to the ORIGINAL graph

`OphGap/Checker.lean` certifies a graph whose vertices are `0, …, n-1` with the spanning
tree placed first. OPH's graph is not stored that way: its vertices are sparse carrier ids
and its edges are in OPH's own order and orientation. This file proves that positive
definiteness of `Dᵀ D` is carried back from the relabelled graph to the original one,
under hypotheses that are themselves decided by the Lean kernel on the literal data:

* a vertex relabelling `new : ℕ → ℕ`, injective on the original carrier list (witnessed by
  a left inverse `oldT`), landing in `0, …, n-1`;
* an edge reordering `perm : ℕ → ℕ` (with left inverse `invPerm`) carrying each original edge,
  relabelled by `new` and *canonicalised* (endpoints sorted, so orientation is forgotten), onto
  the equal canonicalised edge of the relabelled list — decided by one linear pass over the
  lists with `O(log n)` table lookups.

The mathematical content is the quadratic form `x ⬝ Dᵀ D x = ∑_e (σ_e x(tgt e) - x(src e))²`:
each summand is invariant under reversing an edge whenever `σ_e = ±1`, and a sum over
`Fin m` is invariant under an index bijection. So the form of the original graph at `x` equals the
form of the relabelled graph at `x ∘ new⁻¹`; injectivity of `D` transports.

No unproved placeholders, no axioms beyond `propext`, `Classical.choice`, `Quot.sound`.
-/

open Matrix

namespace GapScout

/-! ## Lookup tries: `O(log k)` table lookup that the kernel can evaluate cheaply

A binary trie keyed by the bits of `k` from the least significant end. The kernel's
accelerated `Nat.mod`/`Nat.div` make each step a handful of reductions, whereas indexing
into a `List` of 16,000 entries would cost thousands of reductions per lookup. -/

inductive Trie (α : Type) where
  | leaf : Trie α
  | node (v : Option α) (l r : Trie α) : Trie α

def Trie.get {α : Type} : Trie α → Nat → Option α
  | .leaf, _ => none
  | .node v l r, k => if k = 0 then v else if k % 2 = 1 then l.get (k / 2) else r.get (k / 2)

/-- Overwrite one key (fuel-bounded; used only by negative controls). -/
def Trie.set {α : Type} : Nat → Trie α → Nat → Option α → Trie α
  | 0, t, _, _ => t
  | f + 1, .leaf, k, v =>
    if k = 0 then .node v .leaf .leaf
    else if k % 2 = 1 then .node none (Trie.set f .leaf (k / 2) v) .leaf
    else .node none .leaf (Trie.set f .leaf (k / 2) v)
  | f + 1, .node w l r, k, v =>
    if k = 0 then .node v l r
    else if k % 2 = 1 then .node w (Trie.set f l (k / 2) v) r
    else .node w l (Trie.set f r (k / 2) v)

abbrev NTree := Trie Nat

/-- The relabelling function read off a tree (`0` outside the table; the checks below make
sure every carrier is in the table). -/
def newF (t : NTree) (a : Nat) : Nat := (t.get a).getD 0

/-! ## Canonical edges, relabelling, and the kernel-checkable transport certificate -/

/-- Forget orientation: put the smaller endpoint first. -/
def canon (e : Edge) : Edge := if e.1 ≤ e.2.1 then e else (e.2.1, e.1, e.2.2)

/-- Relabel both endpoints and canonicalise. -/
def relab (new : Nat → Nat) (e : Edge) : Edge := canon (new e.1, new e.2.1, e.2.2)

/-- `a` is labelled: `newT` sends it below `n` and `oldT` sends the label back to `a`. -/
def labelled (n : Nat) (newT oldT : NTree) (a : Nat) : Bool :=
  match newT.get a with
  | some v => decide (v < n) && decide (oldT.get v = some a)
  | none => false

def nodesOk (n : Nat) (newT oldT : NTree) (l : List Nat) : Bool := l.all (labelled n newT oldT)

def endpointsOk (n : Nat) (newT oldT : NTree) (es : List Edge) : Bool :=
  es.all fun e => labelled n newT oldT e.1 && labelled n newT oldT e.2.1 && decide (e.1 ≠ e.2.1)

/-- Strictly increasing (hence duplicate-free). -/
def sortedOk : List Nat → Bool
  | [] => true
  | [_] => true
  | a :: b :: l => decide (a < b) && sortedOk (b :: l)

/-- `t` tabulates `canon` of the list, starting at key `k`: `t.get (k+i) = some (canon l[i])`. -/
def edgesTOk (t : Trie Edge) : Nat → List Edge → Bool
  | _, [] => true
  | k, e :: es => decide (t.get k = some (canon e)) && edgesTOk t (k + 1) es

/-- `perm` sends original index `i` to an index `j < m` of the relabelled list, `invPerm` sends
it back, and the tabulated relabelled edge `j` is the relabelled original edge `i`. -/
def permOk (new : Nat → Nat) (perm invPerm : NTree) (t : Trie Edge) (m : Nat) :
    Nat → List Edge → Bool
  | _, [] => true
  | i, e :: es =>
    (match perm.get i with
      | some j => decide (j < m) && decide (invPerm.get j = some i) &&
          decide (t.get j = some (relab new e))
      | none => false) && permOk new perm invPerm t m (i + 1) es

/-- THE TRANSPORT CERTIFICATE. `nodes`/`origEs` are the original graph, `es` the relabelled
one, `newT`/`oldT` the vertex relabelling and its inverse, `perm`/`invPerm` the edge
reordering and its inverse, `edgesT` a lookup table for `es.map canon`. Every table is a
hint that is re-checked here; none is trusted. -/
def transportCheck (n : Nat) (newT oldT : NTree) (nodes : List Nat) (origEs es : List Edge)
    (perm invPerm : NTree) (edgesT : Trie Edge) : Bool :=
  decide (nodes.length = n) && sortedOk nodes && nodesOk n newT oldT nodes &&
    endpointsOk n newT oldT origEs && decide (origEs.length = es.length) &&
    edgesTOk edgesT 0 es && permOk (newF newT) perm invPerm edgesT es.length 0 origEs

/-! ### Unpacking the checks -/

theorem labelled_spec {n : Nat} {newT oldT : NTree} {a : Nat} (h : labelled n newT oldT a = true) :
    newF newT a < n ∧ oldT.get (newF newT a) = some a := by
  unfold labelled at h
  split at h
  · rename_i v hv
    simp only [Bool.and_eq_true, decide_eq_true_eq] at h
    simp [newF, hv, h]
  · simp at h

theorem nodesOk_spec {n : Nat} {newT oldT : NTree} {l : List Nat} (h : nodesOk n newT oldT l = true) :
    ∀ a ∈ l, newF newT a < n ∧ oldT.get (newF newT a) = some a :=
  fun a ha => labelled_spec (List.all_eq_true.mp h a ha)

theorem endpointsOk_spec {n : Nat} {newT oldT : NTree} {es : List Edge}
    (h : endpointsOk n newT oldT es = true) :
    ∀ e ∈ es, labelled n newT oldT e.1 = true ∧ labelled n newT oldT e.2.1 = true ∧ e.1 ≠ e.2.1 := by
  intro e he
  have := List.all_eq_true.mp h e he
  simpa [Bool.and_eq_true, decide_eq_true_eq, and_assoc] using this

theorem sortedOk_tail : ∀ {a : Nat} {l : List Nat}, sortedOk (a :: l) = true →
    sortedOk l = true ∧ ∀ x ∈ l, a < x
  | _, [], _ => ⟨rfl, by simp⟩
  | a, b :: l, h => by
    simp only [sortedOk, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨h1, h2⟩ := h
    refine ⟨h2, ?_⟩
    intro x hx
    simp only [List.mem_cons] at hx
    rcases hx with rfl | hx
    · exact h1
    · exact lt_trans h1 ((sortedOk_tail h2).2 x hx)

theorem sortedOk_pairwise : ∀ {l : List Nat}, sortedOk l = true → l.Pairwise (· < ·)
  | [], _ => List.Pairwise.nil
  | _ :: _, h => List.Pairwise.cons (sortedOk_tail h).2 (sortedOk_pairwise (sortedOk_tail h).1)

theorem sortedOk_nodup {l : List Nat} (h : sortedOk l = true) : l.Nodup :=
  (sortedOk_pairwise h).imp Nat.ne_of_lt

theorem edgesTOk_spec {t : Trie Edge} : ∀ {k : Nat} {l : List Edge}, edgesTOk t k l = true →
    ∀ i (hi : i < l.length), t.get (k + i) = some (canon l[i])
  | _, [], _, i, hi => absurd hi (Nat.not_lt_zero _)
  | k, e :: es, h, i, hi => by
    simp only [edgesTOk, Bool.and_eq_true, decide_eq_true_eq] at h
    cases i with
    | zero => simpa using h.1
    | succ i =>
      have := edgesTOk_spec h.2 i (by simpa using hi)
      simpa [Nat.add_assoc, Nat.add_comm 1 i] using this

theorem permOk_spec {new : Nat → Nat} {perm invPerm : NTree} {t : Trie Edge} {m : Nat} :
    ∀ {k : Nat} {l : List Edge}, permOk new perm invPerm t m k l = true →
    ∀ i (hi : i < l.length), ∃ j, perm.get (k + i) = some j ∧ j < m ∧
      invPerm.get j = some (k + i) ∧ t.get j = some (relab new l[i])
  | _, [], _, i, hi => absurd hi (Nat.not_lt_zero _)
  | k, e :: es, h, i, hi => by
    simp only [permOk, Bool.and_eq_true] at h
    cases i with
    | zero =>
      obtain ⟨h1, _⟩ := h
      split at h1
      · rename_i j hj
        simp only [Bool.and_eq_true, decide_eq_true_eq] at h1
        exact ⟨j, by simpa using hj, h1.1.1, by simpa using h1.1.2, by simpa using h1.2⟩
      · simp at h1
    | succ i =>
      obtain ⟨j, hj⟩ := permOk_spec h.2 i (by simpa using hi)
      refine ⟨j, ?_⟩
      simpa [Nat.add_assoc, Nat.add_comm 1 i] using hj

/-! ### Slicing the checks (each slice is one kernel evaluation in its own module) -/

theorem edgesTOk_append (t : Trie Edge) (l₁ l₂ : List Edge) : ∀ k,
    edgesTOk t k (l₁ ++ l₂) = (edgesTOk t k l₁ && edgesTOk t (k + l₁.length) l₂) := by
  induction l₁ with
  | nil => intro k; simp [edgesTOk]
  | cons e es ih =>
    intro k
    simp only [List.cons_append, edgesTOk, ih, List.length_cons, Bool.and_assoc]
    rw [show k + 1 + es.length = k + (es.length + 1) by omega]

theorem permOk_append (new : Nat → Nat) (perm invPerm : NTree) (t : Trie Edge) (m : Nat)
    (l₁ l₂ : List Edge) : ∀ k,
    permOk new perm invPerm t m k (l₁ ++ l₂) =
      (permOk new perm invPerm t m k l₁ && permOk new perm invPerm t m (k + l₁.length) l₂) := by
  induction l₁ with
  | nil => intro k; simp [permOk]
  | cons e es ih =>
    intro k
    simp only [List.cons_append, permOk, ih, List.length_cons, Bool.and_assoc]
    rw [show k + 1 + es.length = k + (es.length + 1) by omega]

theorem edgesTOk_split (t : Trie Edge) (k a k' : Nat) (l : List Edge) (la : (l.take a).length = a)
    (hk : k + a = k') (h1 : edgesTOk t k (l.take a) = true) (h2 : edgesTOk t k' (l.drop a) = true) :
    edgesTOk t k l = true := by
  rw [← List.take_append_drop a l, edgesTOk_append, la, hk, h1, h2]; rfl

theorem permOk_split (new : Nat → Nat) (perm invPerm : NTree) (t : Trie Edge) (m k a k' : Nat)
    (l : List Edge) (la : (l.take a).length = a) (hk : k + a = k')
    (h1 : permOk new perm invPerm t m k (l.take a) = true)
    (h2 : permOk new perm invPerm t m k' (l.drop a) = true) :
    permOk new perm invPerm t m k l = true := by
  rw [← List.take_append_drop a l, permOk_append, la, hk, h1, h2]; rfl

theorem nodesOk_split (n : Nat) (newT oldT : NTree) (a : Nat) (l : List Nat)
    (h1 : nodesOk n newT oldT (l.take a) = true) (h2 : nodesOk n newT oldT (l.drop a) = true) :
    nodesOk n newT oldT l = true := by
  rw [← List.take_append_drop a l, nodesOk, List.all_append, ← nodesOk, ← nodesOk, h1, h2]; rfl

theorem endpointsOk_split (n : Nat) (newT oldT : NTree) (a : Nat) (l : List Edge)
    (h1 : endpointsOk n newT oldT (l.take a) = true)
    (h2 : endpointsOk n newT oldT (l.drop a) = true) :
    endpointsOk n newT oldT l = true := by
  rw [← List.take_append_drop a l, endpointsOk, List.all_append, ← endpointsOk, ← endpointsOk, h1, h2]
  rfl

/-! ### Refusal through a single failing conjunct (for the controls) -/

theorem transportCheck_false_of_len {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {perm invPerm : NTree} {edgesT : Trie Edge}
    (h : decide (nodes.length = n) = false) :
    transportCheck n newT oldT nodes origEs es perm invPerm edgesT = false := by
  unfold transportCheck; rw [h]; rfl

theorem transportCheck_false_of_nodes {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {perm invPerm : NTree} {edgesT : Trie Edge}
    (h : nodesOk n newT oldT nodes = false) :
    transportCheck n newT oldT nodes origEs es perm invPerm edgesT = false := by
  unfold transportCheck; rw [h]; simp

theorem transportCheck_false_of_perm {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {perm invPerm : NTree} {edgesT : Trie Edge}
    (h : permOk (newF newT) perm invPerm edgesT es.length 0 origEs = false) :
    transportCheck n newT oldT nodes origEs es perm invPerm edgesT = false := by
  unfold transportCheck; rw [h]; simp

/-- A labelled id is one of the original carriers: `i ↦ new nodes[i]` is an injection
`Fin nodes.length → Fin n`, hence (same cardinality) a surjection, and `oldT` then
identifies `a` with a carrier. No scan of the 8,662-element list is evaluated for this. -/
theorem mem_of_labelled {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    (hlen : nodes.length = n) (hnd : nodes.Nodup) (hN : nodesOk n newT oldT nodes = true)
    {a : Nat} (ha : labelled n newT oldT a = true) : a ∈ nodes := by
  classical
  obtain ⟨hlt, hold⟩ := labelled_spec ha
  have spec := nodesOk_spec hN
  let g : Fin nodes.length → Fin n :=
    fun i => ⟨newF newT nodes[i.1], (spec _ (List.getElem_mem i.2)).1⟩
  have hg : Function.Injective g := by
    intro i j hij
    have h1 := (spec _ (List.getElem_mem i.2)).2
    have h2 := (spec _ (List.getElem_mem j.2)).2
    have hv : newF newT nodes[i.1] = newF newT nodes[j.1] := congrArg Fin.val hij
    rw [hv, h2] at h1
    exact Fin.ext ((hnd.getElem_inj_iff).mp (Option.some.inj h1).symm)
  have hs : Function.Surjective g := (Finite.injective_iff_surjective_of_equiv (finCongr hlen)).mp hg
  obtain ⟨i, hi⟩ := hs ⟨newF newT a, hlt⟩
  have hv : newF newT nodes[i.1] = newF newT a := congrArg Fin.val hi
  have h1 := (spec _ (List.getElem_mem i.2)).2
  rw [hv, hold] at h1
  rw [Option.some.inj h1]
  exact List.getElem_mem i.2

/-! ## The original graph as a matrix indexed by its own carrier ids -/

/-- Vertices of the original graph: the carrier ids that appear in OPH's node list. -/
abbrev Carrier (nodes : List Nat) := {v : Nat // v ∈ nodes}

/-- Every edge joins two distinct listed carriers. -/
def EndpointsIn (nodes : List Nat) (es : List Edge) : Prop :=
  ∀ e ∈ es, e.1 ∈ nodes ∧ e.2.1 ∈ nodes ∧ e.1 ≠ e.2.1

def origSrc (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) (e : Fin es.length) :
    Carrier nodes :=
  ⟨(es.get e).1, (h _ (List.get_mem es e)).1⟩

def origTgt (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) (e : Fin es.length) :
    Carrier nodes :=
  ⟨(es.get e).2.1, (h _ (List.get_mem es e)).2.1⟩

/-- The signed derivative of the ORIGINAL graph: rows are OPH's edges in OPH's order and
orientation, columns are OPH's carrier ids. -/
def origD (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) :
    Matrix (Fin es.length) (Carrier nodes) ℝ :=
  signedDerivative (origSrc nodes es h) (origTgt nodes es h) (σOf es)

theorem orig_loop_free (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) :
    ∀ e, origSrc nodes es h e ≠ origTgt nodes es h e :=
  fun e heq => (h _ (List.get_mem es e)).2.2 (congrArg Subtype.val heq)

/-! ## Quadratic-form bookkeeping -/

/-- One summand of the quadratic form, on a potential given by raw id. -/
def term (f : Nat → ℝ) (e : Edge) : ℝ := (sgn e.2.2 * f e.2.1 - f e.1) ^ 2

/-- Reversing an edge does not change its summand (this is where `σ = ±1` is used). -/
theorem term_canon (f : Nat → ℝ) (e : Edge) : term f (canon e) = term f e := by
  unfold canon
  split
  · rfl
  · obtain ⟨a, b, s⟩ := e
    cases s <;> simp [term, sgn] <;> ring

theorem sum_map_eq_sum_fin {α : Type} (l : List α) (g : α → ℝ) :
    (l.map g).sum = ∑ i : Fin l.length, g (l.get i) := by
  have h1 : l.map g = List.ofFn (g ∘ l.get) := by rw [← List.map_ofFn, List.ofFn_get]
  rw [h1, List.sum_ofFn]
  rfl

/-! ## The transport theorem -/

/-- TRANSPORT (kernel-triviality form). If the relabelled graph `es` passes `check`, and the
original graph `(nodes, origEs)` is carried onto it by a vertex map `new` that is injective
on `nodes`, with the same multiset of canonical signed edges, then the original `D` has
trivial kernel. -/
theorem origD_injective {n : Nat} {nodes : List Nat} {origEs es : List Edge}
    {cyc : List (Nat × Bool)} (new : Nat → Nat)
    (hcheck : check n es cyc = true) (hends : EndpointsIn nodes origEs)
    (hlt : ∀ a ∈ nodes, new a < n)
    (hinj : ∀ a ∈ nodes, ∀ b ∈ nodes, new a = new b → a = b)
    (hσ : ∃ σ : Fin origEs.length ≃ Fin es.length,
      ∀ i, relab new (origEs.get i) = canon (es.get (σ i)))
    (x : Carrier nodes → ℝ) (hx : origD nodes origEs hends *ᵥ x = 0) : x = 0 := by
  obtain ⟨σ, hσ⟩ := hσ
  classical
  have hok := check_edgesOk hcheck
  let φ : Carrier nodes → Fin n := fun a => ⟨new a.1, hlt a.1 a.2⟩
  have hφ : Function.Injective φ := by
    intro a b hab
    exact Subtype.ext (hinj a.1 a.2 b.1 b.2 (congrArg Fin.val hab))
  let y : Fin n → ℝ := Function.extend φ x 0
  have hy : ∀ a, y (φ a) = x a := fun a => hφ.extend_apply x 0 a
  let fy : Nat → ℝ := fun v => if h : v < n then y ⟨v, h⟩ else 0
  let fx : Nat → ℝ := fun a => if h : a ∈ nodes then x ⟨a, h⟩ else 0
  have hfy : ∀ a (ha : a ∈ nodes), fy (new a) = x ⟨a, ha⟩ := by
    intro a ha
    simp only [fy, dif_pos (hlt a ha)]
    exact hy ⟨a, ha⟩
  have hfx : ∀ a (ha : a ∈ nodes), fx a = x ⟨a, ha⟩ := fun a ha => by simp only [fx, dif_pos ha]
  -- the original quadratic form vanishes at x
  have hrow : ∀ e : Fin origEs.length,
      term fx (origEs.get e) = (origD nodes origEs hends *ᵥ x) e ^ 2 := by
    intro e
    rw [origD, signedDerivative_mulVec_apply _ _ _ (orig_loop_free nodes origEs hends)]
    have hm := hends _ (List.get_mem origEs e)
    simp only [term, σOf, origSrc, origTgt, hfx _ hm.1, hfx _ hm.2.1]
  have hSorig : (origEs.map (term fx)).sum = 0 := by
    rw [sum_map_eq_sum_fin]
    apply Finset.sum_eq_zero
    intro e _
    rw [hrow e, hx]
    simp
  -- it equals the relabelled quadratic form at y
  have hrel : ∀ e ∈ origEs, term fy (relab new e) = term fx e := by
    intro e he
    obtain ⟨h1, h2, _⟩ := hends e he
    rw [relab, term_canon]
    simp only [term]
    rw [hfy _ h1, hfy _ h2, hfx _ h1, hfx _ h2]
  have hSrel : (es.map (term fy)).sum = 0 := by
    rw [sum_map_eq_sum_fin]
    calc ∑ j : Fin es.length, term fy (es.get j)
        = ∑ j : Fin es.length, term fy (canon (es.get j)) := by simp only [term_canon]
      _ = ∑ i : Fin origEs.length, term fy (canon (es.get (σ i))) :=
          (Equiv.sum_comp σ (fun j => term fy (canon (es.get j)))).symm
      _ = ∑ i : Fin origEs.length, term fx (origEs.get i) := by
          apply Finset.sum_congr rfl
          intro i _
          rw [← hσ i, hrel _ (List.get_mem origEs i)]
      _ = 0 := by rw [← sum_map_eq_sum_fin]; exact hSorig
  -- so every relabelled row vanishes at y, and the certified graph forces y = 0
  have hD' : signedDerivative (srcOf n es hok) (tgtOf n es hok) (σOf es) *ᵥ y = 0 := by
    funext e
    have hnn : ∀ i ∈ (Finset.univ : Finset (Fin es.length)), 0 ≤ term fy (es.get i) :=
      fun i _ => sq_nonneg _
    have hsum : ∑ i : Fin es.length, term fy (es.get i) = 0 := by
      rw [← sum_map_eq_sum_fin]; exact hSrel
    have hall := (Finset.sum_eq_zero_iff_of_nonneg hnn).mp hsum e (Finset.mem_univ _)
    have hm := edgesOk_mem hok _ (List.get_mem es e)
    rw [signedDerivative_mulVec_apply _ _ _ (loop_free n es hok)]
    have h1 : fy (es.get e).2.1 = y (tgtOf n es hok e) := by
      simp only [fy, tgtOf, dif_pos hm.2.1]
    have h2 : fy (es.get e).1 = y (srcOf n es hok e) := by
      simp only [fy, srcOf, dif_pos hm.1]
    simp only [term] at hall
    rw [σOf, ← h1, ← h2, Pi.zero_apply]
    exact (pow_eq_zero_iff two_ne_zero).mp hall
  have hy0 : y = 0 := kernel_trivial_of_check hcheck y hD'
  funext a
  rw [← hy a, hy0]
  rfl

/-- TRANSPORT (positive-definiteness form): the ORIGINAL graph's `Dᵀ D` is positive
definite. -/
theorem origGram_posDef_of_transport {n : Nat} {nodes : List Nat} {origEs es : List Edge}
    {cyc : List (Nat × Bool)} (new : Nat → Nat)
    (hcheck : check n es cyc = true) (hends : EndpointsIn nodes origEs)
    (hlt : ∀ a ∈ nodes, new a < n)
    (hinj : ∀ a ∈ nodes, ∀ b ∈ nodes, new a = new b → a = b)
    (hσ : ∃ σ : Fin origEs.length ≃ Fin es.length,
      ∀ i, relab new (origEs.get i) = canon (es.get (σ i))) :
    ((origD nodes origEs hends)ᵀ * origD nodes origEs hends).PosDef := by
  apply gram_posDef_of_injective
  intro x1 x2 h12
  dsimp only at h12
  have h0 : origD nodes origEs hends *ᵥ (x1 - x2) = 0 := by rw [mulVec_sub, h12, sub_self]
  exact sub_eq_zero.mp (origD_injective new hcheck hends hlt hinj hσ _ h0)

/-! ## From the Boolean certificate to the hypotheses -/

theorem transportCheck_spec {n : Nat} {newT oldT : NTree} {nodes : List Nat} {origEs es : List Edge}
    {perm invPerm : NTree} {edgesT : Trie Edge}
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    nodes.length = n ∧ sortedOk nodes = true ∧ nodesOk n newT oldT nodes = true ∧
      endpointsOk n newT oldT origEs = true ∧ origEs.length = es.length ∧
      edgesTOk edgesT 0 es = true ∧
      permOk (newF newT) perm invPerm edgesT es.length 0 origEs = true := by
  unfold transportCheck at hT
  simp only [Bool.and_eq_true, decide_eq_true_eq] at hT
  exact ⟨hT.1.1.1.1.1.1, hT.1.1.1.1.1.2, hT.1.1.1.1.2, hT.1.1.1.2, hT.1.1.2, hT.1.2, hT.2⟩

/-- The edge tables give an index equivalence carrying each original edge to its relabelled
copy. (`perm` is injective because `invPerm` is a left inverse; equal lengths then make it a
bijection — no surjectivity pass is evaluated.) -/
theorem index_equiv_of_check {n : Nat} {newT oldT : NTree} {nodes : List Nat} {origEs es : List Edge}
    {perm invPerm : NTree} {edgesT : Trie Edge}
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    ∃ σ : Fin origEs.length ≃ Fin es.length,
      ∀ i, relab (newF newT) (origEs.get i) = canon (es.get (σ i)) := by
  classical
  obtain ⟨_, _, _, _, hlen, hE, hP⟩ := transportCheck_spec hT
  have specP := permOk_spec hP
  have specE := edgesTOk_spec hE
  let f : Fin origEs.length → Fin es.length := fun i =>
    ⟨((perm.get i.1).getD 0), by
      obtain ⟨j, hj, hjm, _, _⟩ := specP i.1 i.2
      rw [Nat.zero_add] at hj
      rw [hj]
      exact hjm⟩
  have hf : ∀ i : Fin origEs.length, ∃ j, perm.get i.1 = some j ∧ (f i).1 = j ∧
      invPerm.get j = some i.1 ∧ edgesT.get j = some (relab (newF newT) (origEs.get i)) := by
    intro i
    obtain ⟨j, hj, _, hinv, ht⟩ := specP i.1 i.2
    rw [Nat.zero_add] at hj hinv
    exact ⟨j, hj, by simp [f, hj], hinv, by simpa using ht⟩
  have hinj : Function.Injective f := by
    intro i i' hii'
    obtain ⟨j, hj, hfj, hinv, _⟩ := hf i
    obtain ⟨j', hj', hfj', hinv', _⟩ := hf i'
    have : j = j' := by rw [← hfj, ← hfj', hii']
    subst this
    rw [hinv] at hinv'
    exact Fin.ext (Option.some.inj hinv')
  have hbij : Function.Bijective f :=
    ⟨hinj, (Finite.injective_iff_surjective_of_equiv (finCongr hlen)).mp hinj⟩
  refine ⟨Equiv.ofBijective f hbij, ?_⟩
  intro i
  obtain ⟨j, hj, hfj, _, ht⟩ := hf i
  have hE' := specE (f i).1 (f i).2
  rw [Nat.zero_add] at hE'
  have hfi : edgesT.get (f i).1 = some (relab (newF newT) (origEs.get i)) := by rw [hfj]; exact ht
  rw [hfi] at hE'
  rw [Equiv.ofBijective_apply, List.get_eq_getElem]
  exact Option.some.inj hE'

theorem endpointsIn_of_check {n : Nat} {newT oldT : NTree} {nodes : List Nat} {origEs es : List Edge}
    {perm invPerm : NTree} {edgesT : Trie Edge}
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    EndpointsIn nodes origEs := by
  obtain ⟨hlen, hsorted, hN, hE, _⟩ := transportCheck_spec hT
  intro e he
  obtain ⟨h1, h2, h3⟩ := endpointsOk_spec hE e he
  exact ⟨mem_of_labelled hlen (sortedOk_nodup hsorted) hN h1,
    mem_of_labelled hlen (sortedOk_nodup hsorted) hN h2, h3⟩

theorem hinj_of_check {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    (hN : nodesOk n newT oldT nodes = true) :
    ∀ a ∈ nodes, ∀ b ∈ nodes, newF newT a = newF newT b → a = b := by
  intro a ha b hb hab
  have h1 := (nodesOk_spec hN a ha).2
  rw [hab, (nodesOk_spec hN b hb).2] at h1
  exact (Option.some.inj h1).symm

/-- THE BRIDGE, fully checked: a passing `check` on the relabelled data plus a passing
`transportCheck` give positive definiteness for the ORIGINAL graph. -/
theorem origGram_posDef_of_transportCheck {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {cyc : List (Nat × Bool)} {perm invPerm : NTree} {edgesT : Trie Edge}
    (hcheck : check n es cyc = true)
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    ((origD nodes origEs (endpointsIn_of_check hT))ᵀ *
      origD nodes origEs (endpointsIn_of_check hT)).PosDef := by
  obtain ⟨_, _, hN, _, _⟩ := transportCheck_spec hT
  exact origGram_posDef_of_transport (newF newT) hcheck _ (fun a ha => (nodesOk_spec hN a ha).1)
    (hinj_of_check hN) (index_equiv_of_check hT)

/-! ## Opaque statement wrappers (same device as `GramPosDef`: keeps the elaborator from
unfolding an 11,816-edge literal while elaborating the statement). -/

/-- `Dᵀ D` of the ORIGINAL graph `(nodes, es)` is positive definite. -/
def OrigGramPosDef (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) : Prop :=
  ((origD nodes es h)ᵀ * origD nodes es h).PosDef

/-- Every eigenvalue of the ORIGINAL graph's `Dᵀ D` is strictly positive (the eigenvalues
are Mathlib's `IsHermitian.eigenvalues`, taken from the positive-definiteness proof `hpd`). -/
def OrigEigenvaluesPos (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es)
    (hpd : OrigGramPosDef nodes es h) : Prop :=
  ∀ i : Carrier nodes, 0 < hpd.1.eigenvalues i

/-- The ORIGINAL graph's signed derivative has trivial kernel. -/
def OrigKernelTrivial (nodes : List Nat) (es : List Edge) (h : EndpointsIn nodes es) : Prop :=
  ∀ x : Carrier nodes → ℝ, origD nodes es h *ᵥ x = 0 → x = 0

theorem origEigenvaluesPos_of_transportCheck {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {cyc : List (Nat × Bool)} {perm invPerm : NTree} {edgesT : Trie Edge}
    (hcheck : check n es cyc = true)
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    OrigEigenvaluesPos nodes origEs (endpointsIn_of_check hT)
      (origGram_posDef_of_transportCheck hcheck hT) :=
  fun i => (origGram_posDef_of_transportCheck hcheck hT).eigenvalues_pos i

theorem origKernelTrivial_of_transportCheck {n : Nat} {newT oldT : NTree} {nodes : List Nat}
    {origEs es : List Edge} {cyc : List (Nat × Bool)} {perm invPerm : NTree} {edgesT : Trie Edge}
    (hcheck : check n es cyc = true)
    (hT : transportCheck n newT oldT nodes origEs es perm invPerm edgesT = true) :
    OrigKernelTrivial nodes origEs (endpointsIn_of_check hT) := by
  obtain ⟨_, _, hN, _, _⟩ := transportCheck_spec hT
  exact origD_injective (newF newT) hcheck _ (fun a ha => (nodesOk_spec hN a ha).1)
    (hinj_of_check hN) (index_equiv_of_check hT)

/-! ## VACUITY CONTROLS on a graph that is genuinely gapless

The all-negative 4-cycle (`square_not_posDef` in `SignedGap.lean`) has a nontrivial kernel, so
by the transport theorem NO `(newT, oldT, es, cyc)` with `check` and `transportCheck` both
true can exist for it. We exhibit the two natural attempts and let the kernel refuse them:

* `sqRelabBij`: the identity relabelling (a genuine bijection of the right vertex set) onto a
  certified 4-vertex graph with a *different* edge — refused by the multiset check;
* `sqRelabCollapse`: a non-injective "relabelling" that folds vertex 3 onto vertex 2 —
  refused by `nodesOk` (no left inverse).

And a 3-vertex certified graph cannot receive the 4-vertex square at all (length check). -/

def sqNodes : List Nat := [0, 1, 2, 3]
def sqEs : List Edge := [(0, 1, true), (1, 2, true), (2, 3, true), (3, 0, true)]

/-- Triangle with a pendant vertex: certified (spanning tree + negative triangle). -/
def triPendEs : List Edge := [(0, 1, true), (1, 2, true), (2, 3, true), (0, 2, true)]
def triPendCyc : List (Nat × Bool) := [(0, true), (1, true), (3, false)]

theorem triPend_check : check 4 triPendEs triPendCyc = true := by decide

/-- identity tables on `{0,1,2,3}` -/
def idT4 : NTree :=
  (((Trie.leaf.set 8 0 (some 0)).set 8 1 (some 1)).set 8 2 (some 2)).set 8 3 (some 3)
/-- folds `3 ↦ 2` -/
def collapseT4 : NTree :=
  (((Trie.leaf.set 8 0 (some 0)).set 8 1 (some 1)).set 8 2 (some 2)).set 8 3 (some 2)
/-- table of `triPendEs.map canon` -/
def triPendT : Trie Edge :=
  (((Trie.leaf.set 8 0 (some (0, 1, true))).set 8 1 (some (1, 2, true))).set 8 2
    (some (2, 3, true))).set 8 3 (some (0, 2, true))

/-- Right vertex bijection, identity edge order, wrong edge multiset: REFUSED
(`sqEs[3] = (3,0)` is not `triPendEs[3] = (0,2)`). -/
theorem sqRelabBij_refused :
    transportCheck 4 idT4 idT4 sqNodes sqEs triPendEs idT4 idT4 triPendT = false := by decide

/-- Non-injective vertex map: REFUSED. -/
theorem sqRelabCollapse_refused :
    transportCheck 4 collapseT4 idT4 sqNodes sqEs triPendEs idT4 idT4 triPendT = false := by decide

/-- The certified triangle on 3 vertices cannot receive a 4-carrier original: REFUSED. -/
def triEs : List Edge := [(0, 1, true), (1, 2, true), (0, 2, true)]
def triT : Trie Edge :=
  ((Trie.leaf.set 8 0 (some (0, 1, true))).set 8 1 (some (1, 2, true))).set 8 2 (some (0, 2, true))
theorem tri_check : check 3 triEs [(0, true), (1, true), (2, false)] = true := by decide
theorem sqToTri_refused :
    transportCheck 3 idT4 idT4 sqNodes sqEs triEs idT4 idT4 triT = false := by decide

/-- Positive control for the machinery itself: the same triangle, relabelled by a genuine
permutation `0↦2, 1↦0, 2↦1`, with edges reordered (`0↦0, 1↦2, 2↦1`) and two of them reversed,
is ACCEPTED. -/
def permT3 : NTree := ((Trie.leaf.set 8 0 (some 2)).set 8 1 (some 0)).set 8 2 (some 1)
def permInvT3 : NTree := ((Trie.leaf.set 8 0 (some 1)).set 8 1 (some 2)).set 8 2 (some 0)
def swapT3 : NTree := ((Trie.leaf.set 8 0 (some 0)).set 8 1 (some 2)).set 8 2 (some 1)
def triOrigEs : List Edge := [(2, 1, true), (1, 0, true), (0, 2, true)]
theorem triPerm_accepted :
    transportCheck 3 permT3 permInvT3 [0, 1, 2] triOrigEs triEs swapT3 swapT3 triT = true := by
  decide

/-- And the same with a wrong edge permutation (not injective: `1↦0, 2↦0`): REFUSED. -/
def badPermT3 : NTree := ((Trie.leaf.set 8 0 (some 0)).set 8 1 (some 0)).set 8 2 (some 0)
theorem triBadPerm_refused :
    transportCheck 3 permT3 permInvT3 [0, 1, 2] triOrigEs triEs badPermT3 swapT3 triT = false := by
  decide

end GapScout
