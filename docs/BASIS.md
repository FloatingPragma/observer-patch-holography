# The OPH basis

Observer Patch Holography proposes a closure hypothesis, an operational
agreement condition and an executable specification. It also adopts
state-selection and response-endogeneity principles, motivated by the
requirement of internal realization but not deduced from closure alone.
The conditional theorems of the [theorem library](../paper/README.md) keep
their displayed premises; the [premise register](registers/PREMISE_REGISTER_V3.md)
records the connections between those premises and the executable
specification.

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

The executable specification is the OPH law, stated in the
[simulator specification](../simulator/SPEC.md); the
[model page](../model/README.md) lists what is proven and measured on it. Its carrier is a torus with
three spatial axes and a private cyclic register at every address. Each
constituent moves by reflected routing; two constituents that share an
address in matched slots receive a contact coin and step together as a bound
pair whose two rails form a clock; a reversible capture vertex turns a
head-on collision into such a pair plus a marker that leaves on its own wire;
and a newborn pair carries an age counter of 8192 states and opens its hop
by a reversible partial swap on the profile `(1 - cos(pi u))/2` over 128
updates, then hops fully. The creator may supply this architecture;
observers need not invent its law. On the OPH law the field of one
constituent has a proved 3+1D Klein–Gordon limit, moving clocks keep proper
time at long wavelengths, and formed matter obeys the long-wave relation
`omega^2 = m^2 + c^2 |k|^2`.

The carrier is hardware of the computation and no physical space. Physical
space, time, fields and particles are effective descriptions that the
observers inside the computation use for the records they agree on, and every
test compares distributions of admissible observer records with what
observers in our world record. The theory proves that observers who inhabit
the simulation can describe their experiences as relativistic 3+1D continuum
spacetimes. The model supplies a field with a proved 3+1D continuum limit,
complete instruments and the count clocks of its observers; the readouts and
the refinement family that the theorem consumes in addition are work in
progress.
Other finite constructions retain the hypotheses of their own theorems.

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
| A1 | `observer_screen_architecture` | The declared screen/carrier and complete-response hypotheses of the theorem; their source-specific realization requires a checked map |
| A2 | `observer_meaning_agreement` | The agreement condition, together with the adopted response-endogeneity principle |
| A3 | `conditional_maximum_randomness` | The adopted state-selection principle on its declared reference and constraint space |

The [premise reference](AXIOM_REFERENCE.md) gives the formal statements of
A1, A2 and A3, which this basis places among the principles above. The
machine-readable source is
[the premise registry](../claims/axiom_registry.yaml). A theorem that uses them
states an implication; the [premise register](registers/PREMISE_REGISTER_V3.md)
records the connections to the executable specification.

## The paper set

- The [flagship synthesis](../flagship/from_observer_consensus_to_standard_physics.tex)
  gives the account, states the closure and executable-model hypotheses and
  walks through the tower of results layer by layer.
- The theorem library holds the long proofs: consensus, observers and records, gauge structure, matter, screen microphysics, and spacetime with Einstein dynamics.
