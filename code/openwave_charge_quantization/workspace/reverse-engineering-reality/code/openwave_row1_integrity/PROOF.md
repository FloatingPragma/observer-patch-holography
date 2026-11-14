# Charge readout from a closing contrast sweep

The construction is a finite classical native instrument on twelve ports and
thirty oriented seams. The source contains one or two nonnegative integer
loads. Its preparation $G$ is known, its total load is conserved, and the
closing contrast records reach the observer faithfully. The result is a
conditional statement about simulation-generated records. The charge unit,
current law, scheduler, finite counters and pointer carrier are declared
architecture. The construction establishes no electron, spin, force or
continuum identification.

## Frozen executable inputs

`INPUTS.json` pins the two unmodified AU5 modules, the portable primary native
kernel and pointer, their geometry packet, and provenance. `run.py` checks
these bytes before and after execution. When their workspace paths exist it
also checks all 36 primary files and the frozen AU5 and route inputs. The
primary manifest digest is
`d5e692b0597dbbeafabb88e6f7d417da8f7d7d1fd39a334c08fa73ffe40e8d90`.

The copied `strict_tape_validator.py` digest is
`6c2f1501a4596ba1a09e1e203ea9c0f2411f33286e8b77cc8467aaeacd0c576a`.
The copied `closing_sweep.py` digest is
`0f1ff01daf8aa2afa1c351eb93cdc6243bb1cbb6cced82b9be950baff47155ce`.
Their callable signatures and scientific premises appear in
`provenance/AU5_FREEZE.md` and `provenance/AU5_PROOF.md`.

The public native pointer imports its portable kernel and geometry packet.
Its functions `execute`, `route`, `reset_success`, `favorable_coins` and
`independent_reset_distribution` have identical abstract syntax trees to the
original frozen policy. No pulse, ideal terminal projector or replacement
pointer participates in this instrument. `run.py` compares all 5,400
one-event branches from the 90 total-one and total-two preparations across
the kernel and both AU5 modules before instrument tests.

## Source and observer interfaces

For seam $e=(a,b)$ and coin $c\in\{0,1\}$, the native event reads
$d=n_a-n_b$, transfers

\[
 j(d,c)={d-(2c-1)(|d|\bmod2)\over2},\qquad n'=n-j e_a+j e_b.
\]

A C record retains the seam, coin and current. An E record also retains the
contrast. `source_execute` applies the native source events, then calls the
frozen `closing_sweep` for thirty additional E events, each seam once. It
returns received records separately from its audit dictionary. The audit
holds source snapshots and emitted prefix records for the scientific tests.
Those fields are absent from the arguments of `observe_source` and `admit`.

`observe_source` receives known $G$, received C/E records and the declared
total. It forms \(\widehat n=G-B\widehat F\), compares received contrasts,
and updates its reconstruction and counters using received currents. It
checks source-sector closure, the native current alphabet, the complete
closing seam order and all thirty E records. No source terminal load enters
this calculation. A source-snapshot field added to an otherwise valid claim
has no effect on admission, including a deliberately false field.

The frozen full AU5 validator admits only total two and requires a
`final_n` field. `_au5_adapter` supplies this field as `reconstruct(G,
received_records)`. It never copies an executor audit field. This use checks
the frozen source replay, counter, pointer, fine-weight and conditional-law
rules; its synthetic final-load comparison supplies no independent
information about the source. Integrity comes from the received E sweep
under the premises below. Total one uses the same frozen
`validate_pointer_tape` and the explicit source-record checks.

## Closing-sweep theorem

Let $B$ have column $e_a-e_b$ for seam $(a,b)$. Let $n_0$ be the
actual source at latch L0, \(\widehat n_0=G-B\widehat F\) the observer's
reconstruction, and \(\delta=n_0-\widehat n_0\). Assume:

1. The swept seams connect every port.
2. Every closing seam, coin, contrast and current is received faithfully.
3. The retained preparation $G$ is the actual preparation.
4. The source evolves only through conserved seam currents.

Then every closing contrast agrees with its reconstruction if and only if
\(\widehat n_0=n_0\). A pass also implies equality after every closing
event and at latch L1.

Proof. Both source and reconstruction change by the same received incidence
current during the sweep, so their discrepancy is constant. The contrast
difference on $(a,b)$ is \(\delta_a-\delta_b\). If every comparison
agrees, this difference vanishes on all swept edges. Connectedness makes
\(\delta\) constant, while equal conserved totals make its sum zero;
hence \(\delta=0\). Conversely zero discrepancy gives equal contrasts
and persists under the common updates. This is the frozen AU5 theorem,
with its comparison isolated from executor diagnostics.

The theorem covers arbitrary source histories in the stated domains and
all sweep coin words. It is algebraic. The package's bounded experiments
use exact weighted finite-state recursion, and claim no literal enumeration
of all $2^{30}$ sweep words.

A single lost moving current $j=\pm1$ changes the reconstruction by a
signed seam incidence vector. Exactly nine carrier seams have nonzero
contrast discrepancy. The test checks both signs on every seam. The
mismatching set is independent of sweep order and coin word.

Equal-load reconstructions can have different tapes: omitted opposite
currents or a closed current cycle have $B\Delta F=0$. The closing theorem
certifies load reconstruction and regional charge, not complete historical
reception, winding of every individual register, or the elapsed clock of
omitted cancelling events. The full instrument's resource theorem below
uses its complete emitted tape. A received clock for an adversarially
shortened prefix requires an authenticated schedule or attempt ledger.

A channel that replaces each E record by the reconstruction's own native
prediction produces an accepted tape with wrong source reconstruction.
`run.py` executes that counterexample for both permitted totals and three
orders. Faithful reception is a premise, not a fact established by matching
record checks. A constant discrepancy is invisible to all contrasts;
equal totals exclude it. The constant-offset counterexample is stated
algebraically, without executing an out-of-domain source.

## Latches and finite counters

The observer starts with twelve zero counters and updates

\[
 c_a\leftarrow(c_a+j)\bmod6,\qquad
 c_b\leftarrow(c_b-j)\bmod6.
\]

By induction, for a complete native history
$c_v=(G_v-n_v)\bmod6$. For a region $A$, its counter class is
\(\sum_{v\in A}c_v\bmod6\). At total $L\in\{1,2\}$, the integer
\(w_A=G(A)-n(A)\) lies in $[-L,L]$. Reduction modulo six is injective
there, so `lift` obtains $w_A$ from the finite counter alone. The admitted
charge is $Q_A=n(A)-G(A)=-w_A$, in the declared charge units.

L0 precedes the sweep and L1 follows its thirty events. Both counter
snapshots and reconstructions are retained. The pointer command uses L1,
after every sweep current has advanced the counters. The reported charge
therefore refers to L1. The receipt includes a concrete case in which L0
and L1 classes differ. Complement regions give opposite windings because
the total relative charge is zero. These counters and updates are supplied
observer memory; no autonomous physical counter sensor is derived.

## Actual pointer and unknown initial state

The unmodified native pointer starts with one load on an arbitrary unknown
port. Its controller scans the thirty seams for a nonzero current,
identifies the reached endpoint from that current, routes to port zero,
and attempts $6|w_A|$ steps along the oriented six-edge cycle. A cursor
advances only on a moving native record. Every addressed edge gets at most
$R$ attempts. Discovery, routing and measurement timeouts retain their
full tapes and have no charge value, including no guessed zero.

Admission calls the frozen `validate_pointer_tape` for every initial port.
It retains the compatible columns rather than accepting a supplied hidden
initial-position label. Every emitted pointer event is replayed according
to the actual native law and discovery/routing/measurement policy. Its
winding is the signed moving cycle-step sum divided by six. Empty tapes,
command echoes, deleted moving records, corrupt currents and guessed
timeout readouts are rejected. A completed zero command must still perform
native discovery and routing.

For a fixed initial column $x$, define $h=2^{-R}$, and let $y_k$ be the
other endpoint of the $k$-th incident seam in scan order, for $k=0,\ldots,4$.
If $d(y_k,0)$ is the frozen route length,

\[
 P_{\rm reset}(x)=\sum_{k=0}^4h^k(1-h)^{1+d(y_k,0)},\qquad
 P_{\rm success}(x,w)=P_{\rm reset}(x)(1-h)^{6|w|}.
\]

The remainder is retained failure probability. The observer's unknown
initial state can have any classical distribution over the twelve columns;
its law is the corresponding mixture. No initial-state prior or
normalization on successful runs is imposed.

`pointer_distribution` computes the complete exact law by stage, signed
completed step count and number of attempts, using the actual native kernel
and frozen routing function. Discovery traverses each actual seam in
order; both native outcomes are propagated. Routing and measurement use
their native move/hold branches, with timeout weight $2^{-R}$. Induction
on these finite controller stages gives the same branch law as `execute`.
The script also explores all sixty actual-executor branch trees at retry one
(twelve initial columns and five commands from minus two to plus two) and
compares their entire rational outcome maps with that computation. Equal
native successor states combine their multiplicities; this quotient keeps
all fair-coin mass. Every resulting actual pointer tape is passed through
the frozen validator.

## Complete bounded instruments and clock

The exact joint-law domain comprises twelve one-load preparations and four
two-load preparations (units on \(\{0,0\},\{0,1\},\{0,2\},\{0,3\}\));
two prefix events on seams 0 and 5, recorded as C then E; both fair coins;
a thirty-event closing sweep in identity or reverse order; the region
\(\{0,1,2\}\) and its complement; all twelve unknown pointer columns;
and $R=1$. There are $16\cdot2\cdot2\cdot12=768$ complete instruments.
Each source setting covers $2^{32}$ coin words through integer
multiplicity recursion. The result retains the joint L0 and L1 winding,
actual pointer outcome, partial signed step count and total attempts.

For each source branch, multiplying by the conditional pointer branch law
gives its joint rational weight. Each source distribution and each pointer
distribution has mass one. Their product law therefore has mass one,
including all three timeout outcomes. `RESULTS.json` retains a digest,
normalization, every outcome's unconditional mass, attempt expectation and
support bound for each complete instrument. `run.py` regenerates the full
maps used by those digests.

A fine complete emitted source-plus-pointer record has weight
$2^{-(N+30+A)}$, where $N$ is prefix length and $A$ pointer attempts.
The fair coin is drawn even when both coins cause the same native state
update; recursion merges their equal successors with the correct
multiplicity. In the bounded joint law $N=2$. Discovery costs at most
$30R$, routing $3R$, and measurement $6LR$, so

\[
 T\le N+30+(33+6L)R.
\]

All holds and the thirty disturbing source E events count. The complete
executor examples use $R=2$; exact joint laws use $R=1$. These are
development domains. No total-four configuration, withdrawn acceptance
case or other reserved configuration is executed, and no official row
status changes.
