import Mathlib
import A5OPH

namespace OPH.TrichotomyCases

/-! # The assembled trichotomy case enumeration

Issue #568.  The compact-Lie trichotomy proof enumerates centre dimension
against semisimple dimension.  This file assembles the finite receipts into
one statement whose signature carries every mathematical hypothesis
explicitly, so nothing is hidden inside executable definitions:

* `hdz` — the centre dimension lies in `{0, 1, 5, 6, 7, 11, 12}`: the
  rational-centre subset-sum receipt (`Screen/A5CharacterField.lean`
  proves this list from Galois stability; the torus/cocharacter step that
  makes the centre rational is a declared classical input).
* `hm` — every simple-ideal dimension lies in `{3, 8, 10}`: the declared
  compact-simple dimension list below twelve (classification of compact
  simple Lie algebras, declared classical input).
* `hsum` — centre and semisimple dimensions total twelve: the port-module
  dimension.

`centre_case_enumeration` returns exactly four arithmetic survivors.  The
`dz = 6` survivor is the dimension-six branch, eliminated at the module
level by `Screen/A5SixAxes.lean` (`no_three_plus_three_split`); feeding
that elimination in as the explicit premise `hsix` gives the trichotomy
(`trichotomy_of_branch_elimination`).  Under the inner-action premises —
`dz ≤ 1` from the S5 centre receipt (`Screen/A5PortModule.lean`) and
`dz ≠ 0` from the `su(2)⁴` fixed-space exclusion (declared, Lie-theoretic)
— the Standard-Model case is the unique survivor
(`inner_action_selection`).

BOUNDARY.  All statements are dimension arithmetic over declared lists;
the bridge from Lie algebras to these dimension data (reductive
decomposition, characteristic centre, exponential surjectivity) stays a
declared classical input, and no physical receipt is formalized. -/

open OPHSelection (SimpleDims)

theorem three_le_of_mem {d : ℕ} (hd : d ∈ SimpleDims) : 3 ≤ d := by
  fin_cases hd <;> omega

/-- Semisimple dimension zero forces the empty ideal multiset. -/
theorem sum_eq_zero (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hs : m.sum = 0) : m = 0 := by
  by_contra h
  obtain ⟨a, ha⟩ := Multiset.exists_mem_of_ne_zero h
  have h3 : 3 ≤ a := three_le_of_mem (hm a ha)
  have hle : a ≤ m.sum :=
    Multiset.single_le_sum (fun x _ => Nat.zero_le x) a ha
  omega

/-- A multiset of allowed compact-simple dimensions summing to six is
exactly `{3, 3}`: the dimension-six branch content `su(2) ⊕ su(2)`. -/
theorem sum_eq_six (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hs : m.sum = 6) : m = {3, 3} := by
  have h3 : ∀ d ∈ m, 3 ≤ d := fun d hd => three_le_of_mem (hm d hd)
  have hcard : m.card ≤ 2 := by
    have := Multiset.card_nsmul_le_sum h3
    simp only [smul_eq_mul] at this
    omega
  interval_cases hc : m.card
  · rw [Multiset.card_eq_zero] at hc
    rw [hc] at hs
    simp at hs
  · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp hc
    have haM := hm a (by rw [ha]; simp)
    rw [ha] at hs
    simp at hs
    fin_cases haM <;> omega
  · obtain ⟨a, b, hab⟩ := Multiset.card_eq_two.mp hc
    have haM := hm a (by rw [hab]; simp)
    have hbM := hm b (by rw [hab]; simp)
    rw [hab] at hs ⊢
    simp at hs
    fin_cases haM <;> fin_cases hbM <;> simp_all

/-- The four arithmetic survivors of the centre-versus-semisimple
enumeration: abelian, dimension-six branch, Standard-Model type, and
`su(2)⁴`. -/
theorem centre_case_enumeration
    (dz : ℕ) (m : Multiset ℕ)
    (hdz : dz ∈ ({0, 1, 5, 6, 7, 11, 12} : Finset ℕ))
    (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hsum : dz + m.sum = 12) :
    (dz = 12 ∧ m = 0) ∨ (dz = 6 ∧ m = {3, 3}) ∨
      (dz = 1 ∧ m = {3, 8}) ∨ (dz = 0 ∧ m = {3, 3, 3, 3}) := by
  have hm' : ∀ d ∈ m, d ∈ OPHTrichotomy.SimpleDims := by
    intro d hd
    have h := hm d hd
    simp only [SimpleDims, Finset.mem_insert, Finset.mem_singleton] at h
    simp only [OPHTrichotomy.SimpleDims, Finset.mem_insert,
      Finset.mem_singleton]
    exact h
  fin_cases hdz
  · exact Or.inr (Or.inr (Or.inr
      ⟨rfl, OPHSelection.sum_eq_twelve m hm (by omega)⟩))
  · exact Or.inr (Or.inr (Or.inl
      ⟨rfl, OPHSelection.sum_eq_eleven m hm (by omega)⟩))
  · exact absurd (by omega : m.sum = 7) (OPHTrichotomy.sum_ne_seven m hm')
  · exact Or.inr (Or.inl ⟨rfl, sum_eq_six m hm (by omega)⟩)
  · exact absurd (by omega : m.sum = 5) (OPHTrichotomy.sum_ne_five m hm')
  · exact absurd (by omega : m.sum = 1) (OPHTrichotomy.sum_ne_one m hm')
  · exact Or.inl ⟨rfl, sum_eq_zero m hm (by omega)⟩

/-- The trichotomy: with the dimension-six branch eliminated (the module
argument of `Screen/A5SixAxes.lean`, fed in as the explicit premise
`hsix`), exactly the three paper outcomes survive. -/
theorem trichotomy_of_branch_elimination
    (dz : ℕ) (m : Multiset ℕ)
    (hdz : dz ∈ ({0, 1, 5, 6, 7, 11, 12} : Finset ℕ))
    (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hsum : dz + m.sum = 12)
    (hsix : dz ≠ 6) :
    (dz = 12 ∧ m = 0) ∨ (dz = 1 ∧ m = {3, 8}) ∨
      (dz = 0 ∧ m = {3, 3, 3, 3}) := by
  rcases centre_case_enumeration dz m hdz hm hsum with h | h | h | h
  · exact Or.inl h
  · exact absurd h.1 hsix
  · exact Or.inr (Or.inl h)
  · exact Or.inr (Or.inr h)

/-- The inner-action selection: with the S5 centre bound `dz ≤ 1`
(`Screen/A5PortModule.lean`) and the declared `su(2)⁴` fixed-space
exclusion `dz ≠ 0`, the Standard-Model dimension data is the unique
survivor. -/
theorem inner_action_selection
    (dz : ℕ) (m : Multiset ℕ)
    (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hsum : dz + m.sum = 12)
    (hle : dz ≤ 1)
    (hne : dz ≠ 0) :
    dz = 1 ∧ m = {3, 8} := by
  have h1 : dz = 1 := by omega
  subst h1
  exact ⟨rfl, OPHSelection.sum_eq_eleven m hm (by omega)⟩

end OPH.TrichotomyCases

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.TrichotomyCases.sum_eq_six
#print axioms OPH.TrichotomyCases.centre_case_enumeration
#print axioms OPH.TrichotomyCases.trichotomy_of_branch_elimination
#print axioms OPH.TrichotomyCases.inner_action_selection
