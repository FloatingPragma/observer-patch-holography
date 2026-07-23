import Mathlib

/-!
# Synchronous settling on ranked dependency systems

This module supplies a synchronous counterpart to the manuscript's
asynchronous ranked-functional theorem.  A site's generator may inspect the
whole configuration syntactically, but `causal` proves that its value depends
only on strictly lower-ranked sites.  Consequently one additional rank
settles at every synchronous round.
-/

namespace ObservableNormalForms

universe u v

/-- A ranked dependency system with heterogeneous site-value types. -/
structure RankedSynchronousSystem where
  Site : Type u
  Value : Site → Type v
  rank : Site → ℕ
  generate : (s : Site) → ((t : Site) → Value t) → Value s
  causal : ∀ (s : Site) (x y : (t : Site) → Value t),
    (∀ t : Site, rank t < rank s → x t = y t) →
      generate s x = generate s y

namespace RankedSynchronousSystem

variable (S : RankedSynchronousSystem.{u, v})

abbrev Configuration : Type (max u v) :=
  (s : S.Site) → S.Value s

/-- One synchronous round: rank-zero boundary sites are protected; every
other site recomputes from the pre-round configuration. -/
def synchronousStep (x : S.Configuration) : S.Configuration :=
  fun s =>
    if _h : S.rank s = 0 then x s else S.generate s x

/-- Iteration with the convention that round zero is the input. -/
def synchronousEvolve : ℕ → S.Configuration → S.Configuration
  | 0, x => x
  | n + 1, x => S.synchronousStep (synchronousEvolve n x)

/-- The candidate extension satisfies every generated (positive-rank)
equation. -/
def IsGeneratedExtension (e : S.Configuration) : Prop :=
  ∀ s : S.Site, S.rank s ≠ 0 → S.generate s e = e s

/-- Two configurations agree on the protected rank-zero boundary. -/
def SameBoundary (x e : S.Configuration) : Prop :=
  ∀ s : S.Site, S.rank s = 0 → x s = e s

/-- After `n` rounds, every site of rank at most `n` agrees with the generated
extension.  This is the rank-induction core of synchronous DAG settling. -/
theorem synchronousEvolve_agrees_through_rank
    {x e : S.Configuration}
    (hboundary : S.SameBoundary x e)
    (hext : S.IsGeneratedExtension e) :
    ∀ (n : ℕ) (s : S.Site), S.rank s ≤ n →
      S.synchronousEvolve n x s = e s := by
  intro n
  induction n with
  | zero =>
      intro s hs
      apply hboundary s
      omega
  | succ n ih =>
      intro s hs
      simp only [synchronousEvolve, synchronousStep]
      split
      next hs0 =>
        exact ih s (by simp [hs0])
      next hs0 =>
        calc
          S.generate s (S.synchronousEvolve n x)
              = S.generate s e := by
                apply S.causal
                intro t ht
                apply ih t
                omega
          _ = e s := hext s hs0

/-- A generated extension is a fixed point of the synchronous update. -/
theorem synchronousStep_fixed
    {e : S.Configuration} (hext : S.IsGeneratedExtension e) :
    S.synchronousStep e = e := by
  funext s
  simp only [synchronousStep]
  split
  next _ => rfl
  next hs0 => exact hext s hs0

theorem synchronousEvolve_fixed
    {e : S.Configuration} (hext : S.IsGeneratedExtension e) (n : ℕ) :
    S.synchronousEvolve n e = e := by
  induction n with
  | zero => rfl
  | succ n ih =>
      simp only [synchronousEvolve, ih]
      exact S.synchronousStep_fixed hext

/-- Synchronous depth-settling corollary: if every site has rank at most
`depth`, exactly `depth` rounds reach the generated extension. -/
theorem synchronous_depth_settling
    {x e : S.Configuration} {depth : ℕ}
    (hboundary : S.SameBoundary x e)
    (hext : S.IsGeneratedExtension e)
    (hdepth : ∀ s : S.Site, S.rank s ≤ depth) :
    S.synchronousEvolve depth x = e := by
  funext s
  exact S.synchronousEvolve_agrees_through_rank hboundary hext depth s (hdepth s)

/-- Generated extensions are unique once the boundary and a common finite
rank bound are fixed. -/
theorem generatedExtension_unique
    {e₁ e₂ : S.Configuration} {depth : ℕ}
    (hext₁ : S.IsGeneratedExtension e₁)
    (hext₂ : S.IsGeneratedExtension e₂)
    (hboundary : S.SameBoundary e₁ e₂)
    (hdepth : ∀ s : S.Site, S.rank s ≤ depth) :
    e₁ = e₂ := by
  calc
    e₁ = S.synchronousEvolve depth e₁ :=
      (S.synchronousEvolve_fixed hext₁ depth).symm
    _ = e₂ := S.synchronous_depth_settling hboundary hext₂ hdepth

end RankedSynchronousSystem

end ObservableNormalForms
