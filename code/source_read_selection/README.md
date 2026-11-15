# Constrained information projection and the M1 read-law boundary

This package establishes a criterion for when finite classical A3 selection
can guarantee a read. It covers arbitrary convex constraints, correlated
schedules and faithful nonuniform references, and distinguishes entropy on
complete histories from entropy on a state-determining observer cover.
The theorem is connected to actual scalar transitions and source
interventions. It does **not** derive M1 or give a complete A1--A3 countermodel.
The exact objective and exit are in [CONTRACT.md](CONTRACT.md).

## General constrained selection criterion

Let K be any convex family of nonnegative finite vectors, let r have strictly
positive coordinates, and let every entropy weight be strictly positive.
If p* attains the minimum of the weighted classical KL objective on K, then
every coordinate positive in any feasible q is positive in p*. The reference
need not be feasible. The theorem does not require independent events,
uniform weights, a full simplex, affine constraints or an interior optimizer
supplied as a premise.

This is the classical I-projection support property rather than an entropy
principle introduced here; see Csiszar and Shields,
[Information Theory and Statistics: A Tutorial, Theorem 3.1](https://www.renyi.hu/~csiszar/Publications/Information_Theory_and_Statistics%3A_A_Tutorial.pdf).
The Lean proof uses convexity of x log x and the exact negative t log t
mixture term at a zero coordinate. `weightedKL` uses Mathlib's `klFun`, and
its entropy-plus-linear expansion is proved, not postulated.

Two consequences are kernel checked:

* For complete-history atoms, a failure event has zero probability at p*
  **if and only if every feasible law gives that event zero mass**.
* The same statement holds for any failure readout that is a nonnegative
  linear combination of the scored cover atoms. With local classical
  entropies, flattening the cover gives those atoms; injectivity alone does
  not make every history event such a positive readout.

`SourceConstrainedRead.lean` composes the support theorem with the native
cut theorem: if any feasible schedule assigns positive weight to a word that
never crosses the source cut, the selected complete-history KL optimizer
cannot guarantee the remote read under both source interventions. Conversely,
a guaranteed read forces **every feasible schedule** to exclude every such
word. This is a necessary source-grammar condition, not sufficient routing.
Unlike the uniform-word bound, it permits correlated schedules and
nonuniform references; it asserts positive obstruction mass, not the same
numerical 29% bound for all these families.

The exact affine certificate theorem additionally proves
`D(q||r) = D(q||p*) + D(p*||r)` when the displayed log-ratio annihilates all
feasible moment differences. Positivity and the equality case give a unique
minimum over the entire family. No numerical optimization or finite sampling
stands in for that universal statement.

## Why the actual A3 cover matters

Use the convex family `p=(t,1/2,1/2-t)`, `0 <= t <= 1/2`, and one faithful
reference `r=(1/4,1/4,1/2)`. Both the identity cover and the coarse map
`C(p)=(p0+p1,p2)` determine p on this family. The coarse reference is the
actual pushforward `C(r)=(1/2,1/2)`; neither the reference nor the constraints
is changed between the two problems.

| Entropy scored | Proved unique minimizer | Probability of history atom 0 |
| --- | --- | --- |
| Full-history KL to r | `(1/6,1/2,1/3)` | `1/6` |
| Coarse-cover KL to C(r) | `(0,1/2,1/2)` | `0` |

The full-support feasible witness `(1/4,1/2,1/4)` survives in both families.
There is no contradiction with the support theorem: the coarse selected
atoms both have positive probability. Recovering history atom 0 requires
`p0=C(p)0-1/2`, with a negative offset. State determination does not turn this
into a nonnegative observable of the scored atoms.

Lean proves convexity, compatible faithful reference, injectivity on the
whole feasible family, both unique minima, and the support difference. This
blocks a false extension of the support theorem to every A3 cover. It is a
finite classical control of that inference, not a realization of all A1--A3
clauses or a declaration that these covers are source selected.

A separate two-edge word example fixes the first-letter probability at 3/4
and uses reference `(1/8,3/8,1/4,1/4)`. The unique optimizer over all such laws
is `(3/16,9/16,1/8,1/8)`, in AA, AB, BA, BB order. It is correlated. Only AB
communicates to the remote receiver by the two-step deadline, so the specified
sign decoder fails with probability 7/16. The exact group-ratio certificate
proves the optimizer for every feasible law, including boundary laws.

The nonconvex control has just two feasible distributions: a vertex and the
uniform law. A faithful reference selects the vertex; the exact ratio of
exponentiated four-times-KL scores is 343/256 > 1. This tests the convexity
boundary. A separate successful singleton-face control forces AB by adding
constraints excluding all other words; it is explicitly classified as an
assumed-success constraint, not a source derivation.

## Uniform finite-word quantitative specialization

Let the candidate alphabet have M distinct scalar seams. Let c of them cross
a source region's boundary. Prepare two different payloads inside that
region, with identical initial data everywhere outside it. Use the same
schedule law, read interface and metadata for both interventions.

Every word avoiding the cut leaves **the complete exterior readout history**
identical under the two interventions. This holds for every word length and
every initial real-valued state, not just the controls below. A decoder may
know the whole word and retain every exterior sample; it cannot return
both different payloads correctly on that word. Consequently any accepted
history claimed correct for both inputs must contain a crossing event.
That condition is necessary, not sufficient for transport.

Under `RepairWordSchedule.FullWordKLProjection`, every length-H word has
positive mass M^(-H), and the mass of cut-avoiding words is exactly

```
p_avoid(H) = ((M-c)/M)^H.
```

An independently sampled balanced choice of the two payloads therefore has
decoding error at least `p_avoid(H)/2`. Equivalently, at least one input has
error at least that large; it need not be both. The Lean two-input inequality
adds the two error probabilities before dividing by two. This prevents an
incorrect per-input bound when a decoder always guesses one payload.

When a noncrossing move exists, every prescribed finite horizon has a
positive-probability obstructing word. The theorem is about a finite deadline.
It neither proves failure of eventual delivery nor bounds the runtime of a
source-selected stopping policy. If c>0, the avoidance probability tends to
zero as H grows; crossing alone does not prove correct decoding.

## Captured wiring and exact controls

The existing pinned W12 level-three capture contains 1,280 carriers and
15,360 scalar ports. Expanding its intra-carrier and glued seams gives
**46,050 distinct undirected edges**. Carrier 0's actual cut has **11 edges**;
the verifier does not assume a twelve-edge boundary or fabricate a regular
graph. The candidate alphabet here is these undirected real pair means, not
the directed integer-completion alphabet of the repair-law RFC.

At H=2,279, exact integer arithmetic certifies

```
2491767692 / 4294967296 <= p_avoid <= 2491767693 / 4294967296
balanced-input error >= 622941923 / 2147483648 > 29/100.
```

These are operation-count bounds under the specified random schedule, not
physical probabilities measured in hardware. H=2,279 is also the event count
of one native-update experiment in `code/source_native_updates/`, which
deliberately schedules transport;
its successful readouts are not contradicted by this calculation.

The captured positive control transports a balanced +/-1 payload from ports
0,9 through pairs (1,5), (14,40), (23,38), to receiver ports 45,39 in eight
scalar means, with receiver amplitude +/-1/16. The negative control executes
sixteen real supported means on both sides of the cut without crossing it.
Both source interventions give identical receiver samples. Every consumed
value, writer identity and resulting value is retained for these four
histories. Their ten active ports start at baseline 2 except the two source
rails; inactive captured ports have the same baseline for both inputs.

On the separate chain 0--1--2--3, the producer and verifier independently
reconstruct every word of lengths zero through six under both +/-1 source
interventions: **2,186 histories and 12,030 means**. At length three, only
word (0,1),(1,2),(2,3) succeeds. Conditioning on successful words removes
26/27 of the schedule mass. It is a different optimization problem; omitted
attempts and the acceptance mechanism need M1 resource accounting.

The toy histories require 24,060 scalar reads and writes each, 8,744 initial
preparation writes and 2,186 terminal receiver samples. The four captured
controls use 48 means, 96 scalar reads and writes each, 40 active preparation
writes and eight terminal receiver samples. These describe the small controls,
not the total cost/capacity of a full M1 service or of preparing the entire
captured federation. Evidence construction, exhaustive verifier work and
immutable controller storage are separate from those executed operation
counts. Compact toy commitments are checked by full replay every time.

## Scope against the full axioms

The [canonical reference](../../docs/AXIOM_REFERENCE.md) states the A1--A3
clauses against which this package is scoped. The
[source-selection guide](../../model/SOURCE_SELECTION.md) distinguishes those
clauses from supplied candidate and interface inputs.

| Required object | What is present | Required for full-axiom promotion |
| --- | --- | --- |
| A1 architecture and response | A pinned finite support and real scalar transition action | Complete compatible local algebras, source-generated response/completeness, operational interfaces and the refinement tower for this same model |
| A2 accepted data and endogenous transport | Same-word source interventions and exact exterior locality | A1-typed accepted-record domain, all meaning/naturality diagrams, complete same-response endogenous overlap implementers |
| A3 feasible grammar | General convex-KL support theorem and necessary cut-word exclusion condition; uniform specialization | A1-generated history observables and a factorization proof covering every A2-visible temporal constraint |
| A3 reference, cover and output | Positive-readout criterion and an exact same-reference cover sensitivity control | Source derivation of the actual cover/reference and its failure-readout map, followed by a justified map to accepted metric-neighbour reads |
| M1 implementation | A necessary crossing condition and finite positive/negative controls | Native record service, full q=13/q=21 executions and intervention-preserving refinement with complete cost accounting |

`CoreAxioms.lean` supplies typed shadows; inhabiting them would not discharge
the full response, accepted-domain and constraint-completeness clauses. A
two-word scalar schedule example is likewise not a countermodel satisfying
all of A1--A3.
The support theorem does not derive complete primitive or temporal grammars
from A1--A3 or insert the desired metric read menu as an A3 constraint.

The next selection obligation is specific: derive the complete temporal
constraint grammar and the actual entropy cover/readout from the source
architecture. In the complete-history or positively failure-visible branch,
guaranteed finite-deadline reads require the grammar itself to exclude
cut-avoiding histories. In a coarse-cover branch, optimizer selection requires
its own source-bound cover and output theorem; state determination is
insufficient. An accepted-record stopping/publication policy is another open
route whose correctness and cost must be proved. A policy merely declared
to accept successful routes retains the M1 premise.
The current theorem does not choose between these alternatives, and the
axioms' entailment or non-entailment of the full M1 law is unresolved.

[Temporal acceptance](../source_temporal_acceptance/README.md) derives the
complete linear observation grammar from native prefixes. It characterizes
the individual-record menu, strengthens the crossing test to an exact
observability condition, and proves an arbitrary-continuation erasure
obstruction. Selection of the metric meaning and full A1--A3 interface is
outside those conditional theorems.

[Native read acceptance](../source_read_acceptance/README.md) derives a maximal
sound local partial decoder for a specified uncertainty model, with native
error and stopping-work theorems. Its preparation, paired alphabet and A3
schedule premises are supplied.

## Reproduction

```
python code/source_read_selection/build.py
python code/source_read_selection/verify.py --write-receipt
python code/source_read_selection/verify.py
python -m pytest -q code/source_read_selection
cd Lean
lake build Geometry.SourceReadSelectionAxiomAudit Geometry
```

Forty-eight Lean theorem declarations pass a transitive gate allowing only
`propext`, `Classical.choice` and `Quot.sound`. The verifier does not import
the producer: it reconstructs adjacency independently, propagates exact
source coefficients rather than payload states, and encloses probabilities
by integer Euclidean division. Constrained minima are checked by exact
group-ratio/Pythagorean certificates; cover injectivity by an independently
computed nonzero determinant; and the nonconvex comparison by an exact
rational log-score ratio. All source files in the local proof dependency
closure, including the imported audit command, are pinned. Adversarial tests
reseal modified commitments, invoke the real CLI, corrupt source files and
test missing/commented CI trust gates.
