# Bridge Boundary Theorem Index

This table maps the finite Lean results used by the cross-paper bridge audit.
The default Lake build includes the source library, the finite event algebra,
the observable-normal-form package, and the Screen modules.

| Paper statement | Lean declaration | Default build | Support status | Physical receipt outside Lean |
|---|---|---:|---|---|
| Strict descent and atomic transactions do not imply local confluence | `OPH.BridgeBoundaries.composed_writes_change_protected`, `peak_endpoints_differ` | yes | Exact finite counterexample | Protected-support-complete reads, protected-conflict-complete components, and a checked quotient local diamond |
| Compatible one- and two-party marginals do not determine a union-collar payload | `OPH.BridgeBoundaries.pair_marginals_do_not_determine_global_law` | yes | Exact classical counterexample | Aligned exact-Markov coherence, a canonical recovery with commuting-square and pentagon checks, or a primitive aggregate payload |
| A transitive 12-port set has no invariant literal `8+3+1` partition | `OPH.oph_defect_generator_not_A5_equivariant` | yes | Exact finite anti-bridge | A full-rank physical current map; compact commutator closure; common inner `A5` action; refinement naturality |
| The twelve-dimensional compact bracket exists on an abstract coefficient algebra | declarations in `Screen/Compact12.lean` | yes | Exact algebraic witness | Identification of that coefficient algebra with physical currents |
| Accepted repair, equilibrium relaxation, and modular flow are different maps | not formalized in this module | n/a | Typed paper boundary | Source-derived reversible proposal rates, an equilibrium generator, a physical noncommutative cap algebra, and the matching cap state |
| The overlap-generated geometric observables form a von Neumann subalgebra | not formalized in this module | n/a | Corrected operator-algebra definition | A normal conditional expectation, if that optional projection is used |
| The hyperbolic three-space is an observer-frame fiber | not formalized in this module | n/a | Corrected geometric type | Affine translation generators, semantic event ancestry, coincidence, and an affine overlap cocycle |
| A conformal celestial sphere does not determine a Lorentz quadratic cone | not formalized in this module | n/a | Explicit non-implication in the papers | Independently inferred quadratic form, eigenvalue gap, inertia `(1,3)`, fit residual, and wrong-signature controls |

The anti-bridges prevent finite bookkeeping from being promoted into a
physical implication. They do not weaken the positive conditional theorems;
they identify the hypotheses those theorems consume.
