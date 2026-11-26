import Mathlib.Tactic

namespace IntervalEcho.AdaptiveRepair

/-- At any history the next offer is selected before its fresh fault.
The three continuations may make wholly different later choices.
The index counts the exposures still owed to one given position. -/
inductive SiteTree : ℕ → Type
  | done : SiteTree 0
  | miss {r : ℕ} (correct hold wrong : SiteTree r) : SiteTree r
  | hit {r : ℕ} (correct hold wrong : SiteTree r) : SiteTree (r+1)

def finalBit {r : ℕ} (P : SiteTree r) (x : ℚ) : ℚ :=
  match P with
  | .done => x
  | .miss a b c => (12*finalBit a x + 3*finalBit b x + finalBit c x)/16
  | .hit a b c => (12*finalBit a 0 + 3*finalBit b x + finalBit c 1)/16

def potential (r : ℕ) (x : ℚ) := 1/13 + (3/16)^r * (x-1/13)

theorem every_adapted_site {r : ℕ} (P : SiteTree r) (x : ℚ) :
    finalBit P x = potential r x := by
  induction P generalizing x with
  | done => simp [finalBit, potential]
  | miss a b c ha hb hc =>
    simp only [finalBit, ha, hb, hc]
    ring
  | hit a b c ha hb hc =>
    simp only [finalBit, ha, hb, hc, potential, pow_succ]
    ring

def covered (offer : Fin 8) (site : Fin 6) : Bool :=
  decide (site.val / 3 = offer.val % 2) &&
    match offer.val / 2 with
    | 0 => decide (site.val % 3 ≠ 2)
    | 1 => decide (site.val % 3 ≠ 1)
    | 2 => decide (site.val % 3 ≠ 0)
    | _ => true

/-- Every chamber of a no-replacement traversal has exactly three hits per site. -/
theorem three_visits : ∀ i : Fin 6,
    (Finset.univ.filter (fun o : Fin 8 => covered o i)).card = 3 := by decide

/-- A controller chooses a remaining label before drawing the current fault.
The continuation can separately depend on each possible outcome. -/
inductive Schedule : Finset (Fin 8) → Type
  | done : Schedule ∅
  | offer {rest : Finset (Fin 8)} (o : Fin 8) (unused : o ∉ rest)
      (correct hold wrong : Schedule rest) : Schedule (insert o rest)

def visits (remaining : Finset (Fin 8)) (i : Fin 6) :=
  (remaining.filter (fun o => covered o i)).card

def projectSite {remaining : Finset (Fin 8)} (plan : Schedule remaining)
    (i : Fin 6) : SiteTree (visits remaining i) :=
  match plan with
  | .done => by simpa [visits] using SiteTree.done
  | @Schedule.offer rest o unused a b c =>
    if h : covered o i = true then by
      simpa [visits, Finset.filter_insert, h, unused] using
        (SiteTree.hit (projectSite a i) (projectSite b i) (projectSite c i))
    else by
      simpa [visits, Finset.filter_insert, h] using
        (SiteTree.miss (projectSite a i) (projectSite b i) (projectSite c i))

theorem whole_schedule (plan : Schedule Finset.univ) (initial : Fin 6 → ℚ) :
    ∑ i, finalBit (projectSite plan i) (initial i) =
      27/4096 * (∑ i, initial i) + 1878/4096 := by
  have term (i : Fin 6) : finalBit (projectSite plan i) (initial i) =
      27/4096 * initial i + 313/4096 := by
    rw [every_adapted_site]
    have three : visits Finset.univ i = 3 := three_visits i
    rw [three]
    unfold potential
    ring
  simp only [term, Finset.sum_add_distrib, ← Finset.mul_sum]
  norm_num

/-- A prefix may reserve different next offers on different noise branches.
The remaining upper exposure budget is common, and may be slack at a leaf. -/
inductive Prefix : ℚ → ℕ → Type
  | stop (d : ℚ) (budget : ℕ) : Prefix d budget
  | next {d a b c : ℚ} {k n : ℕ}
      (drift : (12*a+3*b+c)/16 ≤ d + (k:ℚ)/16)
      (good : Prefix a n) (hold : Prefix b n) (wrong : Prefix c n) : Prefix d (k+n)

def prefixMean {d : ℚ} {n : ℕ} (P : Prefix d n) : ℚ :=
  match P with
  | .stop x _ => x
  | .next _ a b c => (12*prefixMean a + 3*prefixMean b + prefixMean c)/16

theorem prefix_bound {d : ℚ} {n : ℕ} (P : Prefix d n) :
    prefixMean P ≤ d + (n:ℚ)/16 := by
  induction P with
  | stop d n =>
    simp only [prefixMean, le_add_iff_nonneg_right]
    positivity
  | next drift a b c ha hb hc =>
    simp only [prefixMean, Nat.cast_add]
    linarith

theorem native_drift (d k bad : ℚ) (nonnegative : 0 ≤ bad) :
    (12*(d-bad)+3*d+(d+k-bad))/16 ≤ d + k/16 := by linarith

def cost (o : Fin 8) := (Finset.univ.filter (fun i : Fin 6 => covered o i)).card

theorem exposure_budget (s : Finset (Fin 8)) (proper : s ≠ Finset.univ) :
    ∑ i ∈ s, cost i ≤ 16 := by
  have exists_missing : ∃ a, a ∉ s := by
    by_contra! h
    apply proper
    ext a
    simp [h]
  obtain ⟨a, ha⟩ := exists_missing
  have sub : s ⊆ Finset.univ.erase a := by
    intro x hx
    simp only [Finset.mem_erase, Finset.mem_univ, and_true]
    exact fun eq => ha (eq ▸ hx)
  have bound : ∀ o : Fin 8, ∑ i ∈ Finset.univ.erase o, cost i ≤ 16 := by decide
  exact (Finset.sum_le_sum_of_subset_of_nonneg sub (by intros; omega)).trans (bound a)

theorem strict_prefix {d : ℚ} (P : Prefix d 16) : prefixMean P ≤ d + 1 := by
  simpa using prefix_bound P

end IntervalEcho.AdaptiveRepair
