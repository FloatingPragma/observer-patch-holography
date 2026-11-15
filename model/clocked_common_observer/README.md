# The common work observer at its first native current

This construction lets the native current record determine when the common
work observer receives its current-controlled kick. The work register evolves
through every waiting slot, then continues through its four reads and feedback
without another preparation. The retained output includes the native state,
the disturbed work state, the first-current sign and slot, every work record,
and the unresolved timeout branch.

The construction joins the existing
[native first-current instrument](../native_relational_readout/README.md) and
[common source work program](../common_source/COMPOSITION.md). Native stopping
is executed from the actual finite native Kraus maps. Work operations are
compiled from the existing source pulse words and executed by the shared
[density-matrix runtime](../primitive_observer/runtime.py). These two execution
contracts are explicit: the native clock uses finite channel semantics,
not compiled native pulse words.

## One preparation and a retained stopping record

The native register is the invariant binary total-two carrier of dimension
66. Its initial state is the existing mixture of ten coherent phase pairs
`rho_+`. The work register has two scalar sites of occupation dimension two
and the four-dimensional gauge/matter sector, for total dimension 16. Its
initial state is the declared source-prepared finite work vacuum. The two
registers are prepared independently.

Choose native law C or E and a positive integer cap `N`. In each controller
slot the native instrument chooses one of the thirty source seams uniformly
and uses the fair placement coin. Its signed visible current is
`j in {-1,0,+1}`. In that same slot, work executes the source-compiled waiting
map `A_num`, the finite product evolution for Hamiltonian model coordinate
`delta=1/4`. There are 128 product steps in this evolution and in each later
Hamiltonian interval.

The first nonzero current stops the native waiting process. If it occurs at
slot `n`, work has undergone `A_num^n`. The current kick is then applied,
followed by the existing complete continuation:

```text
prepare native and work once
repeat native attempt + work A_num until first current j, or cap N
on first current at n:
    kick(j)
    evolve(1), gauge read r0, feedback(r0), evolve(1/2)
    smooth read r1, evolve(1/2), smooth read r2, feedback(r2)
    evolve(1/2), gauge read r3
on timeout:
    retain the native no-current state and the work state after A_num^N
```

The kick occurs **after** the waiting evolution. The suffix consumes the
actual waiting density matrix, including its coherences. Starting that suffix
from its default preparation would define a different experiment.

Successful records are `(n,j,(r0,r1,r2,r3))`; timeout is `(N,0,())`.
Consequently a cap `N` has `32*N+1` complete output records. The pair `(n,j)`
encodes the entire stopped current word `0^(n-1),j`; timeout encodes `0^N`.
The seam, coin and candidate repair label are deliberately probability-summed
inside each native current outcome. No claim is made that these forgotten
labels remain available as a larger public archive. Every declared current,
work and timeout outcome remains, including zero-probability outcomes.

All output blocks are unnormalized. There is no division by the probability
of obtaining a current before the cap. The native register remains quantum;
it is not replaced by a sampled classical configuration or by its current
probabilities.

## Complete output and composition

Let `T_j` denote the exact native current subchannel, let
`U(rho)=exp(-i*delta*H_2) rho exp(i*delta*H_2)` be the exact waiting channel,
and let `V_j,w` denote the exact work suffix branch with work record `w`.
The target successful block is

```text
T_j T_0^(n-1)(rho_+) tensor V_j,w(U^n(rho_work)),
                 1 <= n <= N, j in {-1,+1}.
```

The timeout block is

```text
T_0^N(rho_+) tensor U^N(rho_work).
```

These tensor factors represent the same experiment, not separately prepared
experiments at each possible stopping slot. Simulation can share prefixes
because each outcome asks for a different branch of that one declared
instrument. Such arithmetic reuse does not provide a physical copying
operation or remove repeated physical operation costs.

The construction validates the entire ordered history family against its
source programs, reconstructs every resource ledger, and checks the same
initial work/reference state throughout. Output iterators snapshot their
factors before returning. In-process integrity seals do not authenticate
an entirely fabricated and resealed history; the shared runtime's importer
replays saved untrusted checkpoints before accepting them.
Certificate inputs must use base NumPy arrays; array subclasses are refused
before their methods can alter support, entries or numeric diagnostics.

The native Kraus families are the existing partial injections, split by the
actual signed current. Their effects sum to the identity on all 66 native
basis states. The work suffix keeps all sixteen records and its final quantum
state. Its complete instrument is trace preserving in the ideal model.
Independence of the initial registers, completeness of both instruments and
the causal selection of the suffix therefore give total output trace one.
No theorem asserting equality of C and E quantum outputs is used: they have
the same signed waiting probabilities but different retained native
coherences.

The existing [finite Kraus criterion and CPTP composition result](../../Lean/Dynamics/ChoiCPTP.lean)
supplies the corresponding complete-positivity argument, including an
untouched reference. The exact hypotheses here are the finite native Kraus
families, unitary waiting operation and complete work instruments. The growing
record space is a direct sum; a literal application of a fixed square-matrix
formalization additionally requires its declared finite tape embedding.
Floating arrays are compared with these exact target maps by a separate
numerical certificate; numerical near-completeness alone is not treated as
an exact CPTP proof.

## Finite numerical certificate

The work certificate acts on all sixteen orthonormal input columns. For
each column, it reuses the existing finite polynomial Hamiltonian and
smooth-pointer enclosures, carrying their remainders through the complete
record instrument. At dimension two these finite-column estimates require
no vacuum or initial occupation-moment premise. The original binary source
branch matrices are compared with those exact polynomial columns by
rational arithmetic, including the rounding used to compose source matrices.

If `eta` bounds the difference of the stacked numerical and ideal branch
maps, the ideal complete stack is an isometry. Consequently

```text
half-diamond error <= eta*(2+eta)/2,
numerical CP trace-norm gain <= (1+eta)^2.
```

These bounds hold on arbitrary work inputs and finite reference extensions.
They can therefore be composed through waiting and the suffix even though
the incoming suffix state differs from the prepared vacuum. The numerical
gain must remain in the error propagation: the binary arrays are not assumed
to define an exactly trace-preserving map.

The executed density matrices require another check. Exact integer arithmetic
compares each returned binary matrix with the action of the retained binary
operators; a directed Frobenius bound supplies a nuclear-norm bound. This
includes small non-Hermitian arithmetic defects without symmetrizing or
renormalizing the reported matrix. The preparation is compared once with the
exact finite dimension-two preparation. The retained numerical execution has
no work reference; the all-input map theorem and separate reference tests
must not be mistaken for a certificate of every later floating application.

More explicitly, let `a` and `b_j` be the waiting and signed-suffix stacked
operator errors. Let `P_n` be the stored numerical waiting matrix product,
and let `r_n` bound the actual multiplication defect in
`P_n - A_num*P_(n-1)`. The producer propagates

```text
h_0 = 0,
h_n <= (1+a)*h_(n-1) + a + r_n,
d_n,j = (1+b_j)*h_n + b_j.
```

The waiting timeout uses `d_n,0=h_n`. If `e_prep` is the preparation error
and `q_n,j` is the sum of actual density-arithmetic errors over the complete
work menu, the work contribution is

```text
w_n,j <= d_n,j*(2+d_n,j)/2 + (1+d_n,j)^2*e_prep + q_n,j.
```

All intermediate widenings are upward rational operations. This propagates
the one preparation error through the numerical norm gain and accounts for
every reused matrix product and every returned density, without charging a
new preparation at the stopping event.

For each native stopping outcome, the native arrays are independently
compared with the exact rational `T_j T_0^(n-1)` target, including timeout.
If the native half nuclear-norm error is `e_r`, its reported nuclear norm is
at most `p_r+2*e_r`, where `p_r` is the exact target mass. If `w_r` bounds the
complete work output error at that outcome, one tensor triangle gives

```text
complete finite output error <= sum_r e_r + sum_r (p_r+2*e_r)*w_r.
```

For a success, `w_r` sums all sixteen work records. There is no extra
factor equal to the number of records. For timeout it concerns the retained
waiting state. This composition does not assume positivity or exact unit
trace of the numerical factors. A comparison with the unlimited exact
finite-occupation observer adds the exact unresolved native mass once.

An effect-probability interval allows twice this half nuclear-norm bound:
the original binary output is not assumed to have exactly unit trace.
Differences between two observers include both errors before asserting a
separation. The certificate covers the exact tensor products of the retained
binary factors; any later floating multiplication forming a large joint
matrix would need its own arithmetic bound.

## Refining the cap

Every successful block with `n<=N` is identical in the cap `N` and cap `M>N`
experiments. The only further execution starts from the unresolved cap `N`
block. It continues both retained registers and stops at the next actual
nonzero current. Thus cap refinement changes the finite execution allowance,
not the preparation or the earlier stopped outputs.

For the supplied phase preparation, the exact unresolved mass is

```text
p_timeout(N) = (2/5)*(13/15)^N + (3/5)*(5/6)^N.
```

The general native input bound is `p_timeout(N)<=(13/15)^N`. It follows from
the diagonal effects of the zero-current native subchannel, including an
untouched reference. Each exact waiting step and complete suffix preserves
trace, so the same tail bounds the complete native/work refinement error.
After embedding outputs into a common record space with separate timeout
labels, the half nuclear-norm discrepancy between a cap and its unlimited
continuation is at most this unresolved mass. A finite cap never claims to
execute the infinite history.

Merely forgetting later stopping labels does **not** recover the cap `N`
timeout quantum block: later native and work operations have acted.
Agreement of coarsened classical survival probabilities is a weaker claim.
The quantum compatibility statement is forward continuation from the
retained timeout state.

## Clock meaning and controls

At full activity the slot `n` counts native attempts. Work's model coordinate at the
end of a successful path is `n/4+5/2`; timeout reaches `N/4`. These quantities
are neither source pulse counts nor laboratory durations. The Hamiltonian,
per-slot model increment, preparation, coherent source services and uniform
native scheduler are declared inputs.

The existing [record-count identities](../../Lean/QFT/SourceRecurrenceClock.lean)
apply to `omega(k)=abs(j_k)` with counted symbol one. They supply additivity,
monotonicity and the bound by the number of slots. Their separately specified
two-state return law is not imported as this native waiting law. The
[count-clock enclosure](../../Lean/Time/SourceCountClockEnclosure.lean)
requires additional count-volume and proper-time premises that this packet
does not supply.

An independent identity idle, with probability `1-q` in each native slot,
replaces `T_0` by `(1-q)Id+q*T_0` and `T_j` by `q*T_j`. For the native register
alone, its unlimited first-current instrument is invariant under this change
after forgetting the stopping-slot record and retaining the sign and native
state. The full slot-labelled native instrument changes with `q`.
Here work advances by `A_num` in every slot, with exact target `U`.
Its stopping-time state
therefore generally changes. The joint observer cannot inherit the
native-only idling invariance. The retained half-activity control uses this
distinction explicitly.

The retained reset control replaces the continued work by a fresh vacuum at
the event. Its first-slot gauge observation differs from the continued-state
prediction. Both controls keep the timeout population; neither compares
renormalized successful samples. The kick's declared position after `A_num^n`
is also bound by the source program and its complete output certificate.

The work marginal is the same under C and E for this independently prepared
experiment: it depends only on their common `(n,j)` probabilities and the
same work continuation. The joint quantum outputs can differ through
the native factor. A work-only readout is not asserted to select between
the native candidates.

## Resource and remaining-obligation boundaries

Work preparation is charged once. At full native activity, a successful path
with record `(n,j,w)` uses
`n` native attempts, `n` waiting work evolutions and one current kick with
one complete four-read suffix. A timeout path uses `N` native attempts and
`N` waiting work evolutions and performs no suffix. Every repeated source
word invocation is kept in its path count even if its contracted matrix is
cached. Costs of a realized path are distinct from the host work needed to
enumerate all possible paths and from the classical storage for its records.
For the idling control, the observed slot count does not determine how many
hidden native activity draws led to actual native attempts. Per-record
attempt costs are therefore upper bounds; exact expected active attempts
are `q` times the expected controller-slot count. Work incurs
one waiting evolution in each such slot.

At cap 128 the declared record alphabet contains 4,097 outcomes. A fixed
encoding can use nine bits for the native stopping outcome and four work
bits, with the timeout determined by the native outcome. The packet bounds
these retained record widths; it does not assign a physical implementation
cost to memory, clock counters or native hidden labels.

This is a finite occupation construction. The retained target and numerical
certificate do not turn its arbitrary incoming work states into the
vacuum occupation theorem. Infinite-occupation comparison, physical source
and storage errors, a physical halt device, and calibrated clock durations
are separate obligations. There is no empirical agreement claim or
construction of the complete consensus tower. A later repair selected by a
work record is not part of this stopping-and-work output; connecting that
additional feedback is a distinct extension.

## Demonstrated output and budgets

The retained [receipt](results.json) executes both C and E at caps 4, 16, 64
and 128, with source split count 128. Each error below bounds the half
nuclear norm of the complete classical/quantum output, including all work
records and the native state. Eventual comparison refers to unlimited
stopping in the same finite occupation model.

| Cap | Complete records | Finite output error upper | Eventual output error upper |
| ---: | ---: | ---: | ---: |
| 4 | 129 | 0.000009 | 0.515028 |
| 16 | 513 | 0.000018 | 0.072993 |
| 64 | 2,049 | 0.000019 | 0.000066 |
| 128 | 4,097 | 0.000019 | 0.000019 |

The displayed bounds hold for both laws. At cap 128 the exact timeout mass
is approximately `4.48134e-9`, and the eventual final-gauge-one probability
lies in `[0.649195, 0.649269]`. The outward component bounds are `0.000001`
for one waiting step and `0.000016` for either complete signed suffix in
half diamond norm. Native and work density arithmetic are each checked
against `1e-8`, as is the complete mass defect. These arithmetic checks are
included in the complete output budget, rather than added informally to
the reported intervals.

The retained-minus-reset difference for first-slot success with final gauge
bit one lies in `[-0.002136, -0.001999]`. Keeping the actual waiting state
therefore produces a certified absolute difference of at least `0.001999`.
The half-activity native-idling control has finite output error at most
`0.000022` and eventual error at most `0.000089`; its eventual final-gauge-one
probability exceeds the baseline by `[0.005698, 0.006125]`. The latter interval
includes both finite arithmetic budgets and both unresolved tails, and
therefore witnesses failure of native-only idling invariance for this common
work observer.

The cap-128 maximum path costs 2,906,209 work source pulses, with one work
preparation and four fresh pointer blanks. Its expected work source pulse
upper bound is below 233,466.599; the exact rational bound is retained in
the receipt. The expected number of full-activity native attempts is
approximately 6.5999999665, with at most 128. These are mathematical source
operation budgets, with native channel attempts counted separately from
the work source pulses; they do not supply physical service error bounds.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.clocked_common_observer --check model/clocked_common_observer/results.json
python -m pytest -q -W error model/clocked_common_observer
```

Replay regenerates the native source branches, source work operations and
complete retained output summaries. Published probabilities include the
timeout; exact native calculations and directed work error bounds are kept
separate from floating diagnostics.
