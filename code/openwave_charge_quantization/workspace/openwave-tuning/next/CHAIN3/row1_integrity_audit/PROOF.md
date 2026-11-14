# Independent record-only row1 integrity audit

This audit concerns the finite native source with total one or two, known
actual preparation, conserved total, and faithful closing E records. It
validates conditional simulation-generated charge records at L1. It does
not identify an observed particle or award an official row. The all31 goal
requires 31 validated obligations; a prepared packet alone is not validation.

The producer is frozen by `MANIFEST.json` SHA256
`1f8c91d4820fa2fe5f439f0a8e344645b4c29b11e3541ac32b7c52a71f1d73c4`.
`PINNED_INPUTS.json` binds the producer, its receipts, all 54 runtime and
external frozen entries, and the primary/AU5 provenance. The audit writes
only its own directory and checks all input hashes before and after.

## Independent expected values

The audit does not invoke the producer's `run.py` or its experiment helpers.
Its independent native transition splits the pair total: for seam `(a,b)`,
`s=n[a]+n[b]`, set `n'[a]=floor(s/2)+coin*(s mod 2)` and `n'[b]=s-n'[a]`.
The current is `n[a]-n'[a]`. Since `s` and `n[a]-n[b]` have the same parity,
this equals the declared current formula. It is compared with the native
kernel and both frozen AU5 modules for all 5,400 branches in the admitted
90-state source space.

Routes are recomputed by queuing complete paths in seam order. A separate
controller transition consumes one native event and explicitly tracks the
occupied port, stage, remaining targets, retries already spent, and signed
measurement steps. Discovery advances the seam after a complete block of
holds; a move ends discovery and identifies the reached endpoint. Routing
and measurement advance only on a move. A block of holds in either stage
terminates with the corresponding timeout. Exhausting the six-edge cycle
targets returns the signed winding; an empty measurement target list still
requires discovery and routing.

The absorbing recurrence advances this controller event by event, assigning
each coin half the incoming exact rational mass. It merges only identical
future controller states and preserves elapsed attempts. Terminals retain
reason, winding or null, partial signed steps, and attempts. Every event
either advances a finite controller location, spends a bounded retry, or
terminates, so the recurrence is finite and has mass one. Its complete laws
are compared against the producer for every initial column and command
`-2..2` at retries 1, 2, and 3. Native executor fixtures are generated from
this separate controller, including all failure stages and a timeout after
three completed signed measurement steps. Their entire records and terminal
results must equal the actual frozen executor's results.

For an initial port `x`, only its five incident discovery blocks can move
the unit. If `h=2^-R`, the first moving block has weight `h^k(1-h)`, and the
following `d` routing moves succeed with weight `(1-h)^d`. This derives

`P_reset(x)=sum_k h^k(1-h)^(1+d(y_k,0))`.

Completing the command multiplies by `(1-h)^(6|w|)`. Discovery timeout has
mass `h^5`; routing and measurement timeouts retain all remaining mass.
The event recurrence verifies these formulas without using them to produce
its expected laws. An unknown initial classical distribution gives the
convex mixture of these twelve unconditional laws. No initial prior is
implied by admission, and no successful-run normalization is used.

`branch_tree_counts.py` separately propagates integer path multiplicities
through the independent controller. Coin outcomes with identical load and
current are identified exactly as in the producer's quotient tree. It
reconstructs all sixty reported leaf counts and summed leaf-record lengths,
including the totals 2,970 and 65,910. This is a combinatorial check; it does
not claim to re-execute all sixty complete executor trees. The native
execution evidence is the independent actual-record comparison above.

For all 768 frozen joint instruments, an independent integer recurrence
propagates the two prefix events and thirty closing events. Its state keeps
both the actual L0 load and the current load, so each source setting accounts
for all `2^32` coin words. Only after propagation are regional L0 and L1
windings projected. Multiplication by the independent event-level pointer
law gives the full joint outcome map. Every frozen joint-map digest,
support size, outcome mass, expectation, and maximum attempt count must
match exactly. This avoids rerunning the producer's expensive experiment
pipeline while independently reconstructing the entire claimed laws.

## Source integrity and latch meaning

Let `delta = actual - reconstructed` at L0. Under faithful closing reception,
both states are advanced by the same incidence current at each closing
event, so `delta` is constant throughout the sweep. The contrast difference
on seam `(a,b)` is `delta[a]-delta[b]`. Vanishing differences on a connected
graph imply constant `delta`. Equal conserved totals then imply `delta=0`.
Conversely, zero discrepancy remains zero. This argument covers arbitrary
prefix histories, all sweep orders, and every coin word in the stated
domain. The executable additionally checks every equal-total pair of source
and reconstruction states and both directions of a lost moving current on
every seam. Each single-current discrepancy affects exactly nine contrasts.

Counters advance only from received currents. Induction gives
`c[v]=(G[v]-n[v]) mod 6` for a faithfully reconstructed history. For total
`L` in `{1,2}`, every regional winding is in `[-L,L]`, on which reduction
modulo six is injective. The pointer command comes from the L1 counters.
The thirty closing events can disturb the source; the executable retains
a witness with different L0 and L1 winding and rejects a stale L0 readout.
Source audit snapshots are deliberately poisoned without affecting admission.

The full frozen AU5 validator requires total two and a `final_n` field.
The producer supplies a synthetic reconstruction from G and received
currents. Its terminal comparison verifies replay consistency, not source
authenticity. Total one uses the frozen pointer validator and the explicit
source checks. An empty shortened history with synthetic `final_n` passes
AU5 alone, while a faithful closing sweep rejects an omitted moving current.

Opposite currents and a complete moving six-edge cycle can be omitted while
leaving `B deltaF=0`. The audit executes both cases: the resulting accepted
claim has correct current reconstruction, but its retained attempt count
and fine-record weight differ from the complete emitted history. Therefore
neither historic authenticity, original per-seam winding, nor elapsed clock
of omitted events follows from the closing theorem. An authenticated attempt
ledger or schedule would be additional input.

Replacing every E record by the reconstruction's own native prediction can
produce acceptance with wrong reconstruction. Six explicit counterexamples
are retained for both totals and three orders. Faithful closing reception
is consequently a premise, not an output of admission. A constant nonzero
discrepancy is contrast-invisible; equal total is essential. That last
counterexample is algebraic and does not execute an out-of-domain load.

## Confirmed pointer snapshot defect and additive correction

The raw producer checks every proposed initial pointer column through the
frozen validator, which also compares `meter.final_loads`. On a discovery
timeout this field equals the hidden initial state. A concrete R=1 all-hold
pointer tape is compatible with initial ports 0 and 9. With exactly that
observed tape, setting only `final_loads` to unit 0 produces compatible set
`[0]`; setting it to unit 9 produces `[9]`. Removing the snapshot rejects the
claim. Thus raw admission's compatible set describes tape plus diagnostic
final state, rather than the tape alone. Both accepted variants still return
null charge, so this does not refute the per-column probability laws or
conditional source reconstruction theorem.

`observer_only.admit_record_only` is the additive public admission facade.
It discards the supplied final pointer snapshot, replays the native events
for every hypothesized initial column, and validates each tape with that
synthetic final state. It invokes unchanged producer admission for every
distinct replay-compatible final placeholder and unions the compatible
initial columns. This pruning is compared with the slower enumeration of
all twelve possible final states on 21 representative cases: both load
totals, success, all three timeout stages, partial signed measurement, and
the ambiguous `[0,9]` tape. It requires
all admitted candidates to agree on the observable verdict, including null
timeouts and retained attempt count. It invents no distribution over those
columns and does not sum candidate probabilities.

The returned schema exposes the twelve-entry compatible-initial indicator,
the fine-record weight conditional on a compatible column, and the diagonal
classical effect with that weight in compatible columns and zero elsewhere.
An explicit prior evaluates this effect by a weighted sum. The audit groups
the complete event recurrence into observable outcome effects and checks
that their sum is the twelve-dimensional identity for each of fifteen
command/retry instruments. All failure effects participate. The exact joint
source/pointer recurrences likewise have unit mass separately in every
initial column. Thus refinement to fine tapes and coarse observation both
retain a normalized classical instrument without a chosen prior.

For a fixed initial column and a finite observed tape, native replay fixes
the final pointer state uniquely. That state is among the replay-derived
placeholders. The frozen validator admits that pair exactly when
the tape obeys its policy and claimed outcome. Hence existential elimination
of the placeholder produces exactly the tape-compatible initial columns.
The construction cannot introduce a non-native pointer event or bypass any
source, clock, weight, or probability-law check. The facade's union is also
compared with an independent event-controller replay that never reads
`final_loads`.

The executable checks final-snapshot removal, arbitrary malformed values,
and contradictory valid snapshots for complete traversals and all timeout
stages across all twelve initial columns and all five commands. It retains
the `[0,9]` witness and rejects command echoes, malformed source/closing
records, wrong totals, stale latches, false weights, false conditional laws,
and guessed zero after a partial measurement timeout. The raw frozen package
is left intact. Public claims of record-only admission must use the facade.

## Resource and scope boundary

For a complete emitted fine history with prefix length `N` and pointer
attempt count `A`, the weight is `2^-(N+30+A)`. Even when the two coins give
the same state, they remain distinct fine records and their multiplicities
are included. Discovery costs at most `30R`, routing at most `3R`, and
measurement at most `6LR`, so total attempts are bounded by
`N+30+(33+6L)R`. Holds, sweep events, and failed attempts are retained.

The fine-record weight is conditional on a fixed compatible pointer initial
column, or describes a fine history that includes that column. For a prior
`pi` and an observed tape compatible with columns `S`, its marginal mass is
`2^-(N+30+A) * sum(pi[x] for x in S)`. No prior is given by this construction.
Consequently the twelve coin-unconditional column laws are exact, while a
single prior-independent marginal probability for an observed tape is not
claimed. In the `[0,9]` ambiguity witness the likelihood is the same in both
columns, but compatibility itself does not select their prior weights.

No reserved acceptance configuration or total-four source is executed.
Counters, scheduling, the charge unit, native law, and apparatus remain
declared architecture. The audit proves no autonomous counter sensor,
electron identification, spin, force, or continuum interpretation.
