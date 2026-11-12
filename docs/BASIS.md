# The OPH basis

Observer Patch Holography proposes a closure hypothesis, an operational
agreement condition and an executable specification. It also adopts
state-selection and response-endogeneity principles, motivated by the
requirement of internal realization but not deduced from closure alone.
The conditional theorems of the [theorem library](../paper/README.md) keep
their displayed premises; the [premise map](../model/PREMISES.md) records the
connections between those premises and the executable specification and
which premises remain declared inputs.

## The closure hypothesis

The Universe is a timeless, self-referential structure. Its observers recover
the specification of their own world from completed records and build
hardware whose execution realizes that same history, including its builders
and their hardware. The structure reads itself into existence.

Write `U` for the complete inhabited history, `D(U)` for an account of its
architecture recoverable through records inside it, and `R` for an admissible
realization of that account constructed inside `U`. The closure condition is
`R(D(U)) ≅ U`, where `≅` is agreement in the observable and interventional
structure. As a fixed point of a universe-level closure operator it reads
`T(U) = U`. No external host, external time or global clock is supplied;
construction and execution happen in emergent, mutually calibrated local
times. Fixedness of the complete history permits changing states within it.

The equation is schematic until a candidate specifies the history domain and
the action of `T`. Existence and identification of our Universe with such a
fixed point require evidence for the specified candidate. Different
architectures may realize the same observable universe; architectural
uniqueness is not a completion requirement. The
[working terms](INTERPRETATION_AND_TERMS.md) distinguish complete history,
instantaneous fixed state, stationary distribution and continuing regime.

## The no-outside rule

A self-contained realization accounts for its operations and information
access within the structure. This motivates the following justification
policy and two additional constitutive principles. Their mathematical
content is stated separately: absence of an external implementer does not
by itself select a state or the connected group that implements rechartings.

- **Justification policy.** A premise is justified operationally when observers recover it and the
  specification that implements it runs and reproduces their measurements.
  The carrier, the observer algebra, the repair law and the closure
  coordinates P and N must meet that test; a construction choice needs no
  uniqueness theorem unless uniqueness is part of its claim.
- **State-selection principle (A3).** On the declared feasible state space,
  exact reference, observer cover and positive weights, adopt the least
  informative state relative to that reference under every observer-visible
  constraint. This is a constitutive selection law. A fixed-point equation
  alone supplies neither its relative-entropy objective nor its reference.
- **Response-endogeneity principle (part of A2).** On the complete-response
  branch, require each proper carrier recharting to be implemented by the
  connected reversible response group, up to a unitary centralizer acting
  trivially on the response algebra. This is stronger than implementation by
  some internal operation: an internally available disconnected operation
  can act nontrivially on an abelian current algebra, whose inner action is
  trivial.

The conditional entropy-selection and gauge-type theorems retain these
principles as premises. Closure and operational agreement do not prove their
necessity; a physical application must also realize them on the same source.

## The agreement condition

Observers assign the same operational meaning to accepted data on shared
boundaries, naturally under every restriction, recharting and refinement. A
record is public when it survives every authorized comparison, and "our
Universe" is the world that public record describes. This is a definition and
an acceptance criterion. The [interpretation contract](INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity)
states what inter-observer agreement must preserve.

## The executable specification

The specification consists of the twelve-port carrier, the native repair
law, the observer algebra of states, preparations, measurements, retained
records and feedback, the declared services, and the closure inputs P and N.
One carrier supports the local observer, clock and geometry results; the
spherical federation, whose cross-carrier routing is a declared convention,
is the construction used for larger observer domains and refinement. The
specification is implemented and run in the
[model directory](../model/README.md). Its status is a candidate, judged by
the [acceptance criteria](../model/CONSTRUCTION.md#completion-target-and-acceptance-criteria)
of the completion target.

## Entry point and order

The exposition enters at the executable specification, because it is
concrete and because observers stand there, holding records and inferring.
The hypothesis specifies the internal-realization contract; the adopted
principles and executable choices then have to satisfy their stated
mathematical and empirical tests. The structure itself has no first element
in either sense, and there is no starting point in a timeless loop.

## Theorem-library premises

The library papers state their theorems with three labelled premises. They
name parts of the basis above.

| Library label | Registry name | Part of the basis |
| --- | --- | --- |
| A1 | `observer_screen_architecture` | The carrier, federation and complete-response layer of the executable specification |
| A2 | `observer_meaning_agreement` | The agreement condition, together with the adopted response-endogeneity principle |
| A3 | `conditional_maximum_randomness` | The adopted state-selection principle on its declared reference and constraint space |

The [premise reference](AXIOM_REFERENCE.md) gives the formal statements of
A1, A2 and A3; it is retained byte-identical for pinned receipts, and its
"three core axioms" wording is superseded by this basis. The machine-readable
source is
[the premise registry](../claims/axiom_registry.yaml). A theorem that uses them
states an implication; the [premise map](../model/PREMISES.md) records the
connections to the executable specification and which antecedents remain
declared inputs.

## The paper set

- The [flagship synthesis](../flagship/from_observer_consensus_to_standard_physics.tex)
  gives the account, states the closure and executable-model hypotheses and
  walks through the tower of results layer by layer.
- The [executable observer paper](../paper/executable_observer_construction.tex)
  states the specification, the executed components and the ordered joins
  that remain.
- The theorem library holds the long proofs: consensus, observers and records, gauge structure, matter, screen microphysics, and spacetime with Einstein dynamics.
