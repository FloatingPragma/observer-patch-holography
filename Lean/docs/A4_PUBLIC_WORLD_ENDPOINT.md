# A4 finite public-world endpoint and E2 handoff

Issue: `#699`

Files:

- `Tower/PublicWorldQuotient.lean`
- `Tower/FixedPointEndpoint.lean`
- `Tower/AdaptiveFixedPointEndpoint.lean`
- `ObserverPatchHolography/Execution/AdaptiveRunStratification.lean`
- `ObserverPatchHolography/Execution/AdaptiveRunCounterexamples.lean`

## Audited dependency boundary

A3's `ConsensusTower` supplies finite observer/record labels, public
algebras, and refinement maps. It does **not** supply an `OPHCarrier`, raw
configurations, a public readback from those configurations, a repair family,
or a source-selected physical state. A4 therefore exposes the missing
realization data instead of pretending that A3 determines it.

The repair proofs reuse `ObserverPatchHolography/Primitives.lean`:
`OPHCarrier`, `Records`, `obsMap`, `gaugeEquiv`, `localRepair`,
`acceptedStep`, `acceptedStepLR`, `Repair`, termination,
`repair_respects_gauge`, and the conditional H1--H3 completeness theorem.
Confluence is never imported silently. The existing `demoLR` is a concrete
termination-plus-completeness witness but is provably nonconfluent.

## Literal quotient

`PublicWorldPresentation Signature` contains a finite inhabited raw type and
an explicit `readback`. Its public equivalence is `Setoid.ker readback`, and
`PublicWorld` is the literal Lean quotient. Machine-checked results include:

- `toPublicWorld_eq_iff`: quotient equality iff readback equality;
- `publicSignature_injective`: descended readback is a complete invariant;
- `publicWorldEquivRealizedSignature`: the quotient is exactly the realized
  readback range;
- finite inhabitedness of the public quotient;
- `ConsensusTower.atRegulator`, universe-polymorphic in the raw
  configuration type;
- a hidden-bit nonvacuity example with distinct raw representatives in one
  public class.

No quotient is called physical merely because it exists as a finite type.

## Generic endpoint theorem

`FiniteRepairSystem P` is a finite family of labelled repair maps. Its
accepted `step` relation records only nonstuttering moves. The relevant
interfaces are deliberately separate:

- `CompletedSchedule schedule start` means only that the finite endpoint is
  fixed by every labelled move;
- `CompleteFor consistent` means fixed iff the declared semantic consistency
  predicate holds;
- `QuotientCompatible` contains both repair-output congruence and enabledness
  congruence under public equivalence.

`QuotientCompatible` is stronger than output congruence alone because
completion must not depend on a hidden representative. It yields descended
labelled `publicMove`s, preservation of public equivalence by finite runs, and
representative-independent completion.

The main generic theorem is
`public_endpoint_exists_unique_on_public_class`. Under explicit termination,
confluence, quotient compatibility, and `CompleteFor`, it constructs a finite
completed schedule whose endpoint is consistent and proves that every
completed schedule from every raw representative of the same public class
has the same public endpoint. Same-start equality is stronger: completed raw
endpoints are equal before quotienting.

This is **not** an infinite weak- or strong-fairness theorem. A4 implements an
explicitly named terminal-completion postcondition. Termination proves that at
least one such finite schedule is reachable. Repair labels and histories are
proof witnesses, not a clock, worldline, preferred foliation, or component of
world identity.

## Canonical adaptive attempt packet

The repository's `adaptiveRun` is a state-dependent sequence of attempted
single-site repairs and may stutter when the selected site does not fire. The
adaptive packet separates three claims that a completed finite schedule does
not express:

- every adaptive attempt stream is eventually constant, because every genuine
  change is an accepted repair step and strictly lowers `mismatchCount`;
- `PathwiseWeakFair` rules out a reducible constant tail and therefore makes
  the stable record a normal form;
- the stronger `WorkConserving` condition gives a normal and stable index no
  larger than the initial mismatch count.

`alwaysProbe` supplies a positive work-conserving witness. The reducible
`remoteReader flagUp` run is constant under an unfair scheduler and shows why
eventual constancy does not imply normality.

`AdaptiveFixedPointEndpoint.lean` is the downstream consumer. Under the
separately named completeness and confluence hypotheses, every pathwise
weak-fair adaptive run reaches a stable consistent raw endpoint equal to the
canonical `Repair` endpoint and hence the same `PublicWorld` fixed-point
object. The packet constructs no stochastic policy, active-source law,
hitting probability, rate, physical scheduler, clock, or refinement theorem.

## Cumulative attempt-capacity packet

`CumulativeAttemptCapacity.lean` charges every chosen-site invocation, including
an equality stutter, and keeps genuine-change and stutter counts separate. The
initial mismatch rank bounds genuine changes but cannot bound attempts: the
committed `delayThenProbe` family holds the initial mismatch at one, delays the
first genuine repair for any prescribed finite prefix, and then normalizes.

`BoundedWaste C q sigma` is the exact quantitative premise that closes this
gap. From every scheduler index and every reducible record, one genuine repair
must occur among the next `q + 1` attempts. It gives a stable normal index no
larger than `(q + 1) * mismatchCount initial`; work conservation is the `q = 0`
special case. The arbitrary-width independent-defect carrier attains that
work-conserving threshold exactly, and the different-patch-cardinality TwoCell
source is a second sharp instance.

`CumulativeCapacityEndpoint.lean` is the load-bearing Tower consumer. A
sufficient attempt budget and bounded waste reach the existing
`FixedPointObject` before exhaustion. Completeness and confluence then identify
the raw record with canonical `Repair` and its image with the same public
fixed-point object. No physical time, rate, energy, bandwidth, hardware, fee,
stochastic-policy, infinite-state, or cross-regulator result is supplied.

## Direct OPH primitive packet

`OPHPrimitiveEndpoint.presentation C seed` uses exactly
`Config := OPH.Records C` and `readback := OPH.obsMap C`. Patch-state
finiteness is explicit in `FinitePatchStates C`.

For the repository's constructed repair:

- `localSystem` is exactly `OPH.localRepair`, and its `step` is definitionally
  `OPH.acceptedStep`;
- `local_terminating` reuses `OPH.termination_holds`;
- `localQuotientCompatible` reuses the existing single-move observation and
  firing congruence theorems;
- the recursive `OPH.Repair` is a fixed normal form and is idempotent;
- `publicRepair` descends `OPH.Repair` to the gauge quotient using
  `OPH.repair_respects_gauge`;
- `publicRepair_idempotent` and
  `canonicalEndpoint_rep_independent` prove an idempotent canonical public
  endpoint map independent of raw gauge representative;
- every canonical endpoint belongs to the explicit `FixedPointObject`.

The constructed `localRepair` is not claimed complete on frustrated carriers
and is not claimed confluent. Accordingly, its canonical `Repair` is a
choice-selected normal-form map, not a proof that all schedules agree or that
every normal form is overlap-consistent.

For an arbitrary frustration-free `lr`, H1--H3 supply termination and
`fixed iff OPH.Consistent`. `lr_public_endpoint_exists_unique_on_gauge_class`
adds explicit confluence, gauge congruence, and enabledness congruence and then
proves the full consistent endpoint theorem across schedules and gauge
representatives. No premise is hidden in the theorem name.

## A3 consensus-tower adaptor

`ConsensusTower.OPHPrimitiveReadbackAdapter T r C` is the typed handoff at a
regulator. It carries:

- a raw OPH seed; and
- an injective encoding
  `OPH.Obs C -> T.PublicSignature r`.

The encoding is a required realization receipt, not data manufactured by A3.
Its injectivity makes equality in the attached tower quotient exactly
`OPH.gaugeEquiv`. The adaptor then provides:

- `primitivePresentation`, literally through
  `ConsensusTower.atRegulator`;
- the actual attached local-repair system and its quotient compatibility;
- `primitivePublicRepair`, a well-defined idempotent endomorphism of the
  tower-regulator public quotient;
- an inhabited public fixed-point object, using existing OPH termination;
- `primitiveLR_endpoint_exists_unique_on_gauge_class`, the full conditional
  termination/confluence/completeness/congruence/completion result stated in
  the A3 regulator's public-world type.

## Negative controls

Three finite countermodels isolate the indispensable assumptions.

1. `forkSystem` terminates and has two completed schedules from one start,
   but is nonconfluent and reaches distinct public endpoints.
2. `settleSystem` terminates and is confluent, but its empty truncated
   schedule is not completed and reports a different public value from the
   genuine fixed endpoint.
3. `representativeSystem` terminates, is confluent, and has completed runs,
   yet two publicly equivalent starts settle to different public bits because
   the repair reads a hidden bit. `representative_no_descended_repair` proves
   that no quotient endomorphism can represent this repair.

Thus same-start confluence cannot replace quotient congruence, and a finite
prefix cannot be relabelled fair merely because the system terminates.

## Closure status and E2 obligations

The bounded mathematical deliverables of `#699` are closure-ready: literal
gauge/public quotients, descended repair maps, idempotence, fixed-point
objects, inhabitedness, conditional schedule-and-representative independence,
the OPH primitive adaptor, the A3 regulator adaptor, and all three negative
controls are machine checked without `sorry`.

This closure is conditional finite structure, **not** a source-realized
physical world or a prediction. In particular, no exhibited
example is asserted to combine H1--H3 completeness, confluence, and all gauge
compatibility premises; the known `demoLR` fails confluence. E2 must:

1. construct the readback adaptor on one inhabited enriched tower selected by
   the intended source/model;
2. supply one repair family satisfying the full premise packet if it wants a
   consistent schedule-independent endpoint rather than only the canonical
   choice-selected `OPH.Repair` endpoint;
3. define raw-configuration refinement maps and prove readback and repair
   naturality across regulators.

A4 proves no refinement-independent limit, locally covariant net, physical
clock, repair rate, continuum existence, QFT reconstruction, Einstein
equation, or empirical prediction.
