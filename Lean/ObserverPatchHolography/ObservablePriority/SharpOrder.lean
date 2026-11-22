import ObserverPatchHolography.ObservablePriority.Sharp

/-! Equivalence under every reordering that preserves the priority of incident
    seams. The local commutation proof acts on physical loads and signed collars. -/
namespace ObserverPatchSharp

private theorem disjoint_facts {a b : Edge} (hd : meets a b = false) :
    a ≠ b ∧ incident b a.lo = false ∧ incident b a.hi = false := by
  exact ⟨fun he => by subst b; simp at hd, Bool.or_eq_false_iff.mp hd⟩

theorem disjoint_write_ready (s : Material) {a b : Edge} (hd : meets a b = false) :
    ready (write s b) a = ready s a := by
  obtain ⟨hne, hl, hh⟩ := disjoint_facts hd
  simp [ready, initialReady, write, hne, hl, hh]

private theorem disjoint_write_commute (s : Material) {a b : Edge} (hd : meets a b = false) :
    write (write s a) b = write (write s b) a := by
  obtain ⟨hab, bl, bh⟩ := disjoint_facts hd
  obtain ⟨hba, al, ah⟩ := disjoint_facts ((meets_symm ..).trans hd)
  unfold write
  congr 1
  · funext p
    cases ha : incident a p <;> cases hb : incident b p <;> simp [ha, hb]
  · funext f
    by_cases hb : f = b
    · subst f
      simp [hba, sign, al, ah]
    · by_cases ha : f = a
      · subst f
        simp [hab, sign, bl, bh]
      · simp [hb, ha]

/-- A locally evaluated strict write at its service slot. -/
def service (e : Edge) (s : Material) : Material := if ready s e then write s e else s

theorem disjoint_service_commute (s : Material) {a b : Edge} (hd : meets a b = false) :
    service a (service b s) = service b (service a s) := by
  have hr := disjoint_write_ready s hd
  have hs := disjoint_write_ready s ((meets_symm ..).trans hd)
  cases ha : ready s a <;> cases hb : ready s b <;>
    simp [service, ha, hb, hr, hs, (disjoint_write_commute s hd).symm]

theorem serviceFold_eq_material (L : Source) (xs : List Edge) :
    xs.foldr service (initial L) = material L (accepted (initialReady L) meets xs) := by
  induction xs with
  | nil => exact initial_material L
  | cons a prior ih =>
    simp only [List.foldr_cons, ih, service, accepted]
    rw [← ready_material]
    split
    · exact write_material _ _ _ (by assumption)
    · rfl

private theorem pull_forward (s : Material) (a : Edge) (front tail : List Edge)
    (h : ∀ b ∈ front, meets a b = false) :
    (front ++ a :: tail).foldr service s = (a :: (front ++ tail)).foldr service s := by
  induction front with
  | nil => rfl
  | cons b rest ih =>
    change service b ((rest ++ a :: tail).foldr service s) =
      service a (service b ((rest ++ tail).foldr service s))
    rw [ih (fun d hd => h d (List.mem_cons_of_mem b hd))]
    exact (disjoint_service_commute _ (h b (by simp))).symm

private theorem legal_folds_agree (rank : Edge → Nat) (s : Material)
    {first second : List Edge} (p : first.Perm second)
    (hf : Ordered rank first) (hs : Ordered rank second) :
    first.foldr service s = second.foldr service s := by
  induction first generalizing second with
  | nil =>
    have h : second = [] := p.symm.eq_nil
    subst second
    rfl
  | cons a rest ih =>
    have ha : a ∈ second := p.mem_iff.mp List.mem_cons_self
    obtain ⟨front, tail, rfl⟩ := List.mem_iff_append.mp ha
    have prest : rest.Perm (front ++ tail) := (p.trans List.perm_middle).cons_inv
    obtain ⟨hfront, htail, hcross⟩ := List.pairwise_append.mp hs
    have hs' : Ordered rank (front ++ tail) :=
      List.pairwise_append.mpr ⟨hfront, (List.pairwise_cons.mp htail).2,
        fun b hb d hd => hcross b hb d (List.mem_cons_of_mem a hd)⟩
    have hd : ∀ b ∈ front, meets a b = false := by
      intro b hb
      apply Bool.eq_false_iff.mpr
      intro hm
      have hab := (List.pairwise_cons.mp hf).1 b (prest.mem_iff.mpr (by simp [hb]))
      have hba := hcross b hb a (by simp)
      have hr1 := hab.2 hm
      have hr2 := hba.2 ((meets_symm ..).trans hm)
      omega
    rw [pull_forward s a front tail hd]
    exact congrArg (service a) (ih prest (List.pairwise_cons.mp hf).2 hs')

/-- Two schedulers that finish the same finite screen and respect incident
    priorities give identical final loads and every signed collar. -/
theorem sharpRun_eq_of_legal_perm (rank : Edge → Nat) (L : Source)
    {first second : List Edge} (p : first.Perm second)
    (hf : Ordered rank first) (hs : Ordered rank second) :
    sharpRun rank first L = sharpRun rank second L := by
  rw [sharpRun_eq_material rank first L hf, sharpRun_eq_material rank second L hs,
    ← serviceFold_eq_material, ← serviceFold_eq_material]
  exact legal_folds_agree rank (initial L) p hf hs

def AllSchedulesAgree (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool) :=
  ∀ (full localOrder : List Edge), full.Perm screen → localOrder.Perm (screen.filter keep) →
    Ordered rank full → Ordered rank localOrder →
    ∀ (L : Source) (e : Edge), e ∈ screen → keep e = true →
      (sharpRun rank full L).collar e = (sharpRun rank localOrder L).collar e

theorem predecessorClosed_iff_allSchedulesAgree (rank : Edge → Nat) (screen : List Edge)
    (keep : Edge → Bool) (ho : Ordered rank screen) :
    PredecessorClosed rank screen keep ↔ AllSchedulesAgree rank screen keep := by
  have hfilter := List.Pairwise.filter keep ho
  constructor
  · intro hc full localOrder pf pp hfull hlocal L e he hk
    rw [sharpRun_eq_of_legal_perm rank L pf hfull ho,
        sharpRun_eq_of_legal_perm rank L pp hlocal hfilter]
    exact closed_agrees rank screen keep ho hc L e he hk
  · intro hh
    exact (predecessorClosed_iff_agrees rank screen keep ho).mpr
      (hh screen (screen.filter keep) .rfl .rfl ho hfilter)

#print axioms disjoint_service_commute
#print axioms serviceFold_eq_material
#print axioms sharpRun_eq_of_legal_perm
#print axioms predecessorClosed_iff_allSchedulesAgree

end ObserverPatchSharp
