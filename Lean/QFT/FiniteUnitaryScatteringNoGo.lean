import Mathlib

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A finite-unitary asymptotic-limit no-go

This module isolates one obstruction relevant to the scattering part of the
QFT lane.  In a Hausdorff topological group, ordinary convergence of the
iterates `g ^ n` forces the one-step evolution `g` to be the identity.  The
finite-dimensional unitary corollary therefore says that a nontrivial exact
unitary update cannot itself settle to an operator-norm (equivalently, in
finite dimension, entrywise) large-time limit at fixed cutoff.

The result is deliberately narrower than scattering theory.  A Møller
operator is a limit of a *relative* evolution such as `U(-t) U₀(t)`, not a
limit of either evolution separately.  The final theorem supplies an exact
relative-evolution control: two identical nontrivial dynamics cancel to the
constant identity sequence.  Thus this module neither rules out scattering
nor constructs it.  It proves that a direct fixed-cutoff limit of one
nontrivial exact unitary is not an adequate scattering carrier.  Viable
routes may instead use, for example, comparison dynamics, selected projected
scalar or observable readouts, an infinite-dimensional weak limit reached
through a continuum or infinite-volume passage, nonunitary/open evolution,
or a finite-time recurrence-aware operational readout.  Merely replacing
norm convergence by full weak-operator convergence at the same fixed finite
dimension does not evade the result: the usual finite-dimensional operator
topologies coincide.
-/

namespace OPH.QFT.FiniteUnitaryScatteringNoGo

open Filter Topology

/-! ## The abstract topological-group obstruction -/

/-- If the natural powers of an element of a Hausdorff topological group
converge, then that element is the identity.

The proof compares the same shifted sequence in two ways.  Continuity of
right multiplication sends `g ^ n → a` to `g ^ n * g → a * g`, while
cofinality of `n ↦ n + 1` sends `g ^ n → a` to `g ^ (n + 1) → a`.
Since `g ^ (n + 1) = g ^ n * g`, Hausdorff uniqueness gives `a = a * g`,
and cancellation gives `g = 1`. -/
theorem tendsto_powers_forces_identity
    {G : Type*} [TopologicalSpace G] [Group G] [IsTopologicalGroup G]
    [T2Space G] (g a : G)
    (h : Tendsto (fun n : ℕ => g ^ n) atTop (nhds a)) :
    g = 1 := by
  have hshift : Tendsto (fun n : ℕ => g ^ (n + 1)) atTop (nhds a) := by
    simpa only [Function.comp_apply] using
      h.comp (tendsto_add_atTop_nat 1)
  have hright : Tendsto (fun n : ℕ => g ^ n * g) atTop (nhds (a * g)) :=
    h.mul tendsto_const_nhds
  have hag : a = a * g := by
    apply tendsto_nhds_unique hshift
    simpa only [pow_succ] using hright
  have hone : (1 : G) = g := by
    apply mul_left_cancel (a := a)
    simpa only [mul_one] using hag
  exact hone.symm

/-- Contrapositive form: a nonidentity group element has no ordinary limit
for its sequence of natural powers. -/
theorem nontrivial_powers_have_no_limit
    {G : Type*} [TopologicalSpace G] [Group G] [IsTopologicalGroup G]
    [T2Space G] (g : G) (hg : g ≠ 1) :
    ¬ ∃ a : G, Tendsto (fun n : ℕ => g ^ n) atTop (nhds a) := by
  rintro ⟨a, ha⟩
  exact hg (tendsto_powers_forces_identity g a ha)

/-! ## Fixed-cutoff unitary specialization -/

/-- A nontrivial finite-dimensional unitary matrix has no ordinary
large-time limit as a sequence of unitary operators.  The topology on the
unitary subgroup is inherited from the finite matrix space. -/
theorem finite_unitary_powers_have_no_limit
    {d : Type*} [Fintype d] [DecidableEq d]
    (U : unitary (Matrix d d ℂ)) (hU : U ≠ 1) :
    ¬ ∃ L : unitary (Matrix d d ℂ),
      Tendsto (fun n : ℕ => U ^ n) atTop (nhds L) :=
  nontrivial_powers_have_no_limit U hU

/-- The ambient-matrix form of the fixed-cutoff result.  It does not evade
the obstruction by allowing a putative limit outside the unitary subtype:
the unitary matrices form a closed set, so any ambient limit of unitary
powers is itself unitary. -/
theorem finite_unitary_ambient_powers_have_no_limit
    {d : Type*} [Fintype d] [DecidableEq d]
    (U : unitary (Matrix d d ℂ)) (hU : U ≠ 1) :
    ¬ ∃ L : Matrix d d ℂ,
      Tendsto (fun n : ℕ => ((U ^ n : unitary (Matrix d d ℂ)) : Matrix d d ℂ))
        atTop (nhds L) := by
  rintro ⟨L, hL⟩
  have hLunitary : L ∈ unitary (Matrix d d ℂ) := by
    apply isClosed_unitary.mem_of_tendsto hL
    exact Filter.Eventually.of_forall fun n => (U ^ n).property
  let Lu : unitary (Matrix d d ℂ) := ⟨L, hLunitary⟩
  have hsubtype : Tendsto (fun n : ℕ => U ^ n) atTop (nhds Lu) := by
    apply tendsto_subtype_rng.mpr
    exact hL
  exact hU (tendsto_powers_forces_identity U Lu hsubtype)

/-! ## An explicit boundary: relative evolution is not excluded -/

/-- Relative comparison of a dynamics with itself is identically one, even
when the individual power sequence has no limit.  This elementary control is
included to prevent the no-go from being misquoted as a no-go for Møller
operators or scattering theory. -/
theorem identical_relative_evolution_is_constant
    {G : Type*} [Group G] (g : G) :
    (fun n : ℕ => (g ^ n)⁻¹ * g ^ n) = fun _ => 1 := by
  funext n
  simp

/-- Consequently the identical relative evolution converges, regardless of
whether the underlying one-step evolution is nontrivial. -/
theorem identical_relative_evolution_tendsto
    {G : Type*} [TopologicalSpace G] [Group G] (g : G) :
    Tendsto (fun n : ℕ => (g ^ n)⁻¹ * g ^ n) atTop (nhds 1) := by
  rw [identical_relative_evolution_is_constant]
  exact tendsto_const_nhds

#print axioms tendsto_powers_forces_identity
#print axioms nontrivial_powers_have_no_limit
#print axioms finite_unitary_powers_have_no_limit
#print axioms finite_unitary_ambient_powers_have_no_limit
#print axioms identical_relative_evolution_is_constant
#print axioms identical_relative_evolution_tendsto

end OPH.QFT.FiniteUnitaryScatteringNoGo
