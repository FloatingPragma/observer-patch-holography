import ObserverPatchHolography.ObservablePriority.ObservedChannel

/-! Convert from the convenient chronological indexing in the local multiplication
to an actual field of physical collar labels. Distinct orderings can then be
compared entry by entry without either copying or observing the original input. -/
namespace StrictQuantum.Addressed
open StrictQuantum
open scoped BigOperators
noncomputable section

variable {E V A : Type} [DecidableEq V]

abbrev Collar := Option (Sign × Bool)

def readField [DecidableEq E] (es : List E) (rs : List A) (e : E) : Option A :=
  (es.zip rs).lookup e

theorem field_readback [DecidableEq E] (es : List E) (rs : List A)
    (hn : es.Nodup) (hl : es.length = rs.length) :
    es.map (readField es rs) = rs.map some := by
  induction es generalizing rs with
  | nil => cases rs <;> simp_all
  | cons a es ih =>
    cases rs with
    | nil => simp at hl
    | cons r rs =>
      obtain ⟨hna,hnd⟩ := List.nodup_cons.mp hn
      have hlen : es.length = rs.length := by simpa using hl
      simp only [List.map_cons]
      congr 1
      · simp [readField]
      · calc
          _ = es.map (readField es rs) := by
            apply List.map_congr_left
            intro b hb
            have hab : a ≠ b := fun h => hna (h ▸ hb)
            simp [readField, List.lookup_cons, beq_eq_false_iff_ne.mpr (Ne.symm hab)]
          _ = _ := ih rs hnd hlen

theorem scan_length (q : Screen E V) (es : List E) (x : Loads V) :
    (scan q es x).2.length = es.length := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih => simp [scan, ih]

theorem labels_length (rs : List (Option Sign)) (c : Coins rs) : (labels rs c).length = rs.length := by
  simpa using congrArg List.length (labels_signs rs c)

/-- The only physical output registers used here are port loads and addressed
collars. Blank registers outside this execution window keep their supplied value. -/
def embed [DecidableEq E] (q : Screen E V) (es : List E)
    (z : (x : Loads V) × Coins (scan q es x).2) : Loads V × (E → Collar) :=
  ((scan q es z.1).1, fun e => (readField es (labels (scan q es z.1).2 z.2) e).join)

theorem embed_readback [DecidableEq E] (q : Screen E V) {es : List E} (hn : es.Nodup)
    (x : Loads V) (c : Coins (scan q es x).2) :
    es.map (embed q es ⟨x,c⟩).2 = labels (scan q es x).2 c := by
  have h := field_readback es (labels (scan q es x).2 c) hn
    (by rw [labels_length, scan_length])
  have hw := congrArg (List.map Option.join) h
  simpa [List.map_map, Function.comp_def, embed] using hw

theorem embed_injective [DecidableEq E] (q : Screen E V) {es : List E} (hn : es.Nodup) :
    Function.Injective (embed q es) := by
  rintro ⟨x,c⟩ ⟨y,d⟩ h
  apply physical_injective q es
  apply Prod.ext
  · have hh := congrArg Prod.fst h
    exact hh
  · change labels (scan q es x).2 c = labels (scan q es y).2 d
    rw [← embed_readback q hn x c, ← embed_readback q hn y d, h]

theorem labels_complete (rs : List (Option Sign)) (ys : List Collar)
    (h : ys.map (Option.map Prod.fst) = rs) :
    ∃ c : Coins rs, labels rs c = ys := by
  induction ys generalizing rs with
  | nil => subst rs; exact ⟨(), rfl⟩
  | cons y ys ih =>
    subst rs
    obtain ⟨c, hc⟩ := ih _ rfl
    cases y with
    | none => exact ⟨c, congrArg (List.cons none) hc⟩
    | some l => exact ⟨(l.2,c), by simp [labels,hc]⟩

/-- Matrix entry indexed by a physical bra (y, coll) and a blank-input ket x.
No copy of the input is present in that physical bra. -/
def element (q : Screen E V) (es : List E) (M : E → Bool) (x y : Loads V)
    (coll : E → Collar) : ℝ := by
  classical
  exact if (scan q es x).1 = y ∧
      (scan q es x).2 = es.map (fun e => (coll e).map Prod.fst) ∧
      (∀ e, e ∉ es → coll e = none)
    then Native.matrixCoefficient q M es x else 0

theorem embed_outside [DecidableEq E] (q : Screen E V) {es : List E} (x : Loads V)
    (c : Coins (scan q es x).2) {e : E} (he : e ∉ es) : (embed q es ⟨x,c⟩).2 e = none := by
  have hnone : readField es (labels (scan q es x).2 c) e = none := by
    apply List.lookup_eq_none_iff.mpr
    rintro ⟨a,b⟩ hb
    have ha := (List.of_mem_zip hb).1
    have hne : e ≠ a := fun h => he (h ▸ ha)
    simp [hne]
  simp [embed, hnone]

theorem supported_is_embed [DecidableEq E] (q : Screen E V) {es : List E} (hn : es.Nodup)
    (x y : Loads V) (coll : E → Collar)
    (hl : (scan q es x).1 = y)
    (hr : (scan q es x).2 = es.map (fun e => (coll e).map Prod.fst))
    (hout : ∀ e, e ∉ es → coll e = none) :
    ∃ c : Coins (scan q es x).2, embed q es ⟨x,c⟩ = (y,coll) := by
  obtain ⟨c,hc⟩ := labels_complete (scan q es x).2 (es.map coll)
    (by simpa [List.map_map, Function.comp_def] using hr.symm)
  refine ⟨c, Prod.ext hl ?_⟩
  funext a
  by_cases ha : a ∈ es
  · exact List.map_inj_left.mp ((embed_readback q hn x c).trans hc) a ha
  · exact (embed_outside q x c ha).trans (hout a ha).symm

/-- Every sparse coordinate drawn from the recursive calculation has exactly
its prescribed amplitude at the actual material/coin registers. -/
theorem element_at_embed [DecidableEq E] (q : Screen E V) {es : List E} (hn : es.Nodup)
    (M : E → Bool) (x : Loads V) (c : Coins (scan q es x).2) :
    element q es M x (embed q es ⟨x,c⟩).1 (embed q es ⟨x,c⟩).2 =
      Native.matrixCoefficient q M es x := by
  classical
  have hl := congrArg (List.map (Option.map Prod.fst)) (embed_readback q hn x c)
  simp only [List.map_map, labels_signs] at hl
  unfold element
  exact if_pos ⟨rfl, hl.symm, fun _ he => embed_outside q x c he⟩

theorem weight_perm {as bs : List (Option Sign)} (p : as.Perm bs) :
    weight as = weight bs := by
  induction p with
  | nil => rfl
  | cons a p ih => cases a <;> simp [weight, ih]
  | swap a b l => cases a <;> cases b <;> rfl
  | trans _ _ h₁ h₂ => exact h₁.trans h₂

theorem map_eq_transfer {D : Type} {as bs : List E} (p : as.Perm bs) (f g : E → D) :
    as.map f = as.map g ↔ bs.map f = bs.map g := by
  simp only [List.map_inj_left]
  exact ⟨fun h e he => h e (p.mem_iff.mpr he), fun h e he => h e (p.mem_iff.mp he)⟩

/-- The only scheduling interface is sharp port-load agreement and the same
signed record at each address for this input. These data are classical and
contain no quantum parameters, matrix elements or superpositions. -/
theorem element_schedule_agreement (q : Screen E V) {as bs : List E} (p : as.Perm bs)
    (M : E → Bool) (x y : Loads V) (coll : E → Collar) (record : E → Option Sign)
    (hload : (scan q as x).1 = (scan q bs x).1)
    (ha : (scan q as x).2 = as.map record) (hb : (scan q bs x).2 = bs.map record) :
    element q as M x y coll = element q bs M x y coll := by
  classical
  have hw := weight_perm (p.map record)
  have hm := map_eq_transfer p (fun e => (record e).isSome) M
  have hr := map_eq_transfer p record (fun e => (coll e).map Prod.fst)
  have hout : (∀ e, e ∉ as → coll e = none) ↔ (∀ e, e ∉ bs → coll e = none) := by
    simp only [p.mem_iff]
  simp only [element, Native.matrix_outcome_projection, ha, hb, flags, List.map_map,
    Function.comp_def, hload, propext hm, propext hr, propext hout, hw]

#print axioms embed_injective
#print axioms supported_is_embed
#print axioms element_at_embed
#print axioms element_schedule_agreement

end
end StrictQuantum.Addressed
