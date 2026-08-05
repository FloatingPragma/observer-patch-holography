import FiniteConditionalRepair

namespace OPH.Thermodynamics

/-!
# Stationary source dynamics: the second law without reversibility

The conditional-resampling kernel is reversible, but the transition-side
second law needs less.  Any finite stochastic kernel that preserves a faithful
reference contracts relative entropy to that reference.  This is the exact
interface needed to interpret a source-counted, possibly nonreversible repair
chain without silently importing microscopic reversibility.

The lazy directed three-cycle is an exact rational separation witness: it is
stochastic, preserves the faithful uniform law, and violates detailed balance.
The stationary second law follows.  Detailed fluctuation relations and
Onsager reciprocity therefore remain separate receipts.

No source construction, physical energy, clock, or common-reference
identification is asserted here.
-/

variable {Omega : Type*} [Fintype Omega] [DecidableEq Omega]

/-- **Stationary-kernel second law.** A faithful stationary reference is
enough for relative-entropy descent.  Detailed balance is not a premise. -/
theorem stationary_secondLaw (pi p : Omega → Real)
    (K : Omega → Omega → Real)
    (hp : ∀ x, 0 ≤ p x) (hpi : ∀ x, 0 < pi x)
    (hK0 : ∀ x y, 0 ≤ K x y)
    (hK1 : ∀ x, (Finset.univ.sum fun y => K x y) = 1)
    (hstationary : push pi K = pi) :
    kl (push p K) pi ≤ kl p pi := by
  have h := kl_push_le p pi K hp (fun x => le_of_lt (hpi x))
    (fun x hx => absurd hx (ne_of_gt (hpi x))) hK0 hK1
  rwa [hstationary] at h

/-- Expected value of an observable after one transition from `x`. -/
noncomputable def rowExpectation (K : Omega → Omega → Real)
    (f : Omega → Real) (x : Omega) : Real :=
  Finset.univ.sum fun y => K x y * f y

/-- Every stochastic kernel fixes constant observables.  This is useful as a
control, but it is not evidence for a nontrivial protected charge. -/
theorem constantObservable_fixed {Omega' : Type*} [Fintype Omega']
    (K : Omega' → Omega' → Real)
    (hK1 : ∀ x, (Finset.univ.sum fun y => K x y) = 1)
    (c : Real) (x : Omega') :
    rowExpectation K (fun _ => c) x = c := by
  rw [rowExpectation, ← Finset.sum_mul, hK1 x, one_mul]

/-! ## Exact nonreversible separation witness -/

/-- A lazy clockwise walk on three states. -/
noncomputable def directedLazy3 : Fin 3 → Fin 3 → Real := fun x =>
  ![![1 / 2, 1 / 2, 0],
    ![0, 1 / 2, 1 / 2],
    ![1 / 2, 0, 1 / 2]] x

/-- The faithful uniform reference on the three-cycle. -/
noncomputable def uniform3 : Fin 3 → Real := fun _ => 1 / 3

theorem directedLazy3_nonneg (x y : Fin 3) :
    0 ≤ directedLazy3 x y := by
  fin_cases x <;> fin_cases y <;> norm_num [directedLazy3]

theorem directedLazy3_row_sum (x : Fin 3) :
    (Finset.univ.sum fun y => directedLazy3 x y) = 1 := by
  fin_cases x <;> norm_num [directedLazy3, Fin.sum_univ_succ]

theorem uniform3_pos (x : Fin 3) : 0 < uniform3 x := by
  norm_num [uniform3]

theorem directedLazy3_stationary :
    push uniform3 directedLazy3 = uniform3 := by
  funext y
  fin_cases y <;>
    norm_num [push, uniform3, directedLazy3, Fin.sum_univ_succ]

/-- The stationary witness is not in detailed balance: the clockwise edge
`0 → 1` has positive equilibrium flux and its reverse has zero flux. -/
theorem directedLazy3_not_detailedBalance :
    uniform3 0 * directedLazy3 0 1 ≠
      uniform3 1 * directedLazy3 1 0 := by
  norm_num [uniform3, directedLazy3]

/-- Relative entropy to the uniform law descends for the nonreversible exact
witness. -/
theorem directedLazy3_secondLaw (p : Fin 3 → Real)
    (hp : ∀ x, 0 ≤ p x) :
    kl (push p directedLazy3) uniform3 ≤ kl p uniform3 := by
  exact stationary_secondLaw uniform3 p directedLazy3 hp uniform3_pos
    directedLazy3_nonneg directedLazy3_row_sum directedLazy3_stationary

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.stationary_secondLaw
#print axioms OPH.Thermodynamics.directedLazy3_stationary
#print axioms OPH.Thermodynamics.directedLazy3_not_detailedBalance
#print axioms OPH.Thermodynamics.directedLazy3_secondLaw
