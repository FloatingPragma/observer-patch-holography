import ObserverPatchHolography.Primitives

namespace OPH.Locality

open OPH

/-!
# A certified finite dependency-cone bound for local repair

The constructed single-site repair move reads only the closed edge
neighborhood of its firing site and writes only the site itself. This
module turns that locality into an exact finite upper-bound statement: after
any schedule of `n` local repair moves, the records on a region depend
only on the initial records inside the `n`-fold grown neighborhood of
that region. Changing initial data outside that finite cone cannot
change the readout, for every schedule of that length.

The theorem quantifies over every fixed exogenous repair word shared by the
two compared inputs.  It therefore does not cover an adaptive scheduler whose
next site is chosen from the current global state.  It also does not prove
that this upper bound is minimal or that influence reaches its boundary.  No
probability, continuum limit, propagation speed, or physical distance enters;
the cone is the combinatorial shadow of single-site transactional repair, and
its physical reading requires separate scheduler-locality, length, and clock
attachments.
-/

variable {C : OPHCarrier}

/-- The closed edge neighborhood of a site: the site itself together
with every endpoint of an edge incident to it. This is exactly the
data one repair move at the site reads. -/
def closedNbhd (i : C.Patch) : Set C.Patch :=
  {j | j = i ∨ ∃ e : C.Edge,
    (C.src e = i ∨ C.tgt e = i) ∧ (C.src e = j ∨ C.tgt e = j)}

theorem self_mem_closedNbhd (i : C.Patch) : i ∈ closedNbhd (C := C) i :=
  Or.inl rfl

theorem endpoint_mem_closedNbhd {i : C.Patch} {e : C.Edge}
    (hinc : C.src e = i ∨ C.tgt e = i) {j : C.Patch}
    (hend : C.src e = j ∨ C.tgt e = j) :
    j ∈ closedNbhd (C := C) i :=
  Or.inr ⟨e, hinc, hend⟩

/-- One growth step of a region: the region together with the closed
neighborhoods of all its sites. -/
def grow (S : Set C.Patch) : Set C.Patch :=
  S ∪ ⋃ i ∈ S, closedNbhd i

theorem subset_grow (S : Set C.Patch) : S ⊆ grow S :=
  Set.subset_union_left

theorem closedNbhd_subset_grow {S : Set C.Patch} {k : C.Patch}
    (hk : k ∈ S) : closedNbhd (C := C) k ⊆ grow S :=
  fun _p hp => Or.inr (Set.mem_iUnion₂.mpr ⟨k, hk, hp⟩)

/-- The `n`-fold grown neighborhood of a region: a certified dependency-cone
upper bound for `n` sequential repair moves. -/
def ball (S : Set C.Patch) : ℕ → Set C.Patch
  | 0 => S
  | n + 1 => grow (ball S n)

/-- Dependent-update evaluation agrees on records that agree at the
evaluation site or evaluate at the updated site. -/
theorem update_eval_agree {x y : Records C} (i : C.Patch)
    (s : C.State i) (j : C.Patch) (hj : j = i ∨ x j = y j) :
    Function.update x i s j = Function.update y i s j := by
  by_cases hji : j = i
  · subst hji
    simp only [Function.update_self]
  · simp only [Function.update_of_ne hji]
    rcases hj with hj | hj
    · exact absurd hj hji
    · exact hj

/-- Edge consistency at an incident edge is a function of the records
on the closed neighborhood of the site. -/
theorem edgeConsistentAt_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p)
    {e : C.Edge} (hinc : C.src e = i ∨ C.tgt e = i) :
    edgeConsistentAt e x ↔ edgeConsistentAt e y := by
  unfold edgeConsistentAt
  rw [hagr (C.src e) (endpoint_mem_closedNbhd hinc (Or.inl rfl)),
    hagr (C.tgt e) (endpoint_mem_closedNbhd hinc (Or.inr rfl))]

/-- Post-update edge consistency at an incident edge is a function of
the records on the closed neighborhood of the site. -/
theorem edgeConsistentAt_update_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p) (s : C.State i)
    {e : C.Edge} (hinc : C.src e = i ∨ C.tgt e = i) :
    edgeConsistentAt e (Function.update x i s)
      ↔ edgeConsistentAt e (Function.update y i s) := by
  unfold edgeConsistentAt
  rw [update_eval_agree i s (C.src e) ?hsrc,
    update_eval_agree i s (C.tgt e) ?htgt]
  case hsrc =>
    by_cases h : C.src e = i
    · exact Or.inl h
    · exact Or.inr
        (hagr (C.src e) (endpoint_mem_closedNbhd hinc (Or.inl rfl)))
  case htgt =>
    by_cases h : C.tgt e = i
    · exact Or.inl h
    · exact Or.inr
        (hagr (C.tgt e) (endpoint_mem_closedNbhd hinc (Or.inr rfl)))

/-- The firing trigger reads only the closed neighborhood. -/
theorem localTrigger_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p) :
    LocalTrigger i x ↔ LocalTrigger i y :=
  exists_congr fun _e => and_congr_right fun hinc =>
    not_congr (edgeConsistentAt_agree hagr hinc)

/-- Transactional solvability reads only the closed neighborhood,
pointwise in the candidate state. -/
theorem solvesAt_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p) (s : C.State i) :
    SolvesAt i x s ↔ SolvesAt i y s :=
  forall_congr' fun _e => imp_congr_right fun hinc =>
    edgeConsistentAt_update_agree hagr s hinc

theorem locallySolvable_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p) :
    LocallySolvable i x ↔ LocallySolvable i y :=
  exists_congr (solvesAt_agree hagr)

/-- `Classical.choose` picks the same witness from pointwise-equal
predicates. -/
private theorem choose_eq_of_pred_iff {α : Sort*} {p q : α → Prop}
    (hpq : ∀ a, p a ↔ q a) (hp : ∃ a, p a) (hq : ∃ a, q a) :
    Classical.choose hp = Classical.choose hq := by
  have hpq' : p = q := funext fun a => propext (hpq a)
  subst hpq'
  rfl

/-- **One repair move reads one closed neighborhood.** If two records
agree on the closed neighborhood of the firing site, the repaired
records agree at the site and at every other point where the inputs
agree. -/
theorem localRepair_agree {x y : Records C} {i : C.Patch}
    (hagr : ∀ p ∈ closedNbhd (C := C) i, x p = y p)
    (j : C.Patch) (hj : j = i ∨ x j = y j) :
    localRepair C i x j = localRepair C i y j := by
  by_cases hx : LocalTrigger i x ∧ LocallySolvable i x
  · have hy : LocalTrigger i y ∧ LocallySolvable i y :=
      ⟨(localTrigger_agree hagr).mp hx.1,
        (locallySolvable_agree hagr).mp hx.2⟩
    have hs : Classical.choose hx.2 = Classical.choose hy.2 :=
      choose_eq_of_pred_iff (solvesAt_agree hagr) hx.2 hy.2
    rw [localRepair_of_fire C i x hx, localRepair_of_fire C i y hy,
      ← hs]
    exact update_eval_agree i (Classical.choose hx.2) j hj
  · have hy : ¬ (LocalTrigger i y ∧ LocallySolvable i y) := fun hy' =>
      hx ⟨(localTrigger_agree hagr).mpr hy'.1,
        (locallySolvable_agree hagr).mpr hy'.2⟩
    rw [localRepair_of_quiescent C i x hx,
      localRepair_of_quiescent C i y hy]
    rcases hj with hj | hj
    · subst hj
      exact hagr j (self_mem_closedNbhd j)
    · exact hj

/-- A repair schedule: a finite word of firing sites, applied left to
right. -/
noncomputable def applyWord (w : List (Site C)) (x : Records C) :
    Records C :=
  w.foldl (fun z k => localRepair C k z) x

@[simp] theorem applyWord_nil (x : Records C) :
    applyWord ([] : List (Site C)) x = x := rfl

theorem applyWord_cons (k : Site C) (w : List (Site C))
    (x : Records C) :
    applyWord (k :: w) x = applyWord w (localRepair C k x) := rfl

/-- One growth step transports agreement through one repair move: if
the inputs agree on `grow T`, the outputs of any single move agree on
`T`. -/
theorem localRepair_agree_on_grow {x y : Records C} (k : Site C)
    (T : Set C.Patch)
    (hagr : ∀ p ∈ grow (C := C) T, x p = y p) :
    ∀ j ∈ T, localRepair C k x j = localRepair C k y j := by
  intro j hj
  by_cases hjk : j = k
  · subst hjk
    have hk : j ∈ T := hj
    exact localRepair_agree
      (fun p hp => hagr p (closedNbhd_subset_grow hk hp))
      j (Or.inl rfl)
  · rw [localRepair_apply_of_ne C k x j hjk,
      localRepair_apply_of_ne C k y j hjk]
    exact hagr j (subset_grow T hj)

/-- **Finite dependency-cone bound.** For every repair schedule `w` and
region `S`, agreement of the initial records on the cone
`ball S w.length` forces agreement of the final records on `S`. The
readout after `n` sequential repair moves on a region is a function of the initial
data inside the `n`-fold grown neighborhood, for every fixed exogenous word of
length `n`. State-dependent schedule selection, a converse, and minimality are
not claimed. -/
theorem applyWord_agree_on (w : List (Site C)) (S : Set C.Patch)
    {x y : Records C}
    (hagr : ∀ p ∈ ball (C := C) S w.length, x p = y p) :
    ∀ j ∈ S, applyWord w x j = applyWord w y j := by
  induction w generalizing x y with
  | nil =>
      intro j hj
      simpa using hagr j hj
  | cons k w' ih =>
      intro j hj
      rw [applyWord_cons, applyWord_cons]
      exact ih
        (localRepair_agree_on_grow k (ball S w'.length)
          (by simpa [ball, List.length_cons] using hagr))
        j hj

/-- **No influence outside the cone.** Changing the initial record at
a site outside `ball S n` cannot change the readout on `S` after any
repair schedule of length `n`. -/
theorem no_influence_outside_ball (w : List (Site C))
    (S : Set C.Patch) (x : Records C) {p : C.Patch}
    (hp : p ∉ ball (C := C) S w.length) (s : C.State p) :
    ∀ j ∈ S, applyWord w (Function.update x p s) j
      = applyWord w x j :=
  applyWord_agree_on w S fun q hq =>
    Function.update_of_ne (fun h => hp (by rw [← h]; exact hq)) s x

end OPH.Locality

#print axioms OPH.Locality.localRepair_agree
#print axioms OPH.Locality.applyWord_agree_on
#print axioms OPH.Locality.no_influence_outside_ball
