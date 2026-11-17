# Native selection and count-clock rest energy

Let L be the complete240-position native graph Laplacian and
T=I-L/1260. The native endpoint census supplies T, including every seam,
coin, current and no-operation record. The preparation program is
U_prep(t)=exp(-itH_prep), H_prep=-1260 log T. Its complete M6 native implementation is the
source-word construction bound in PREREG.json, with its full processor and
source-padding failure instruments. This is a specified coherent program;
unread classical endpoint repair is a different instrument.

Let P project onto the first twelve-port cell. Two integer vectors are

    f=e4-e5-e6+e7,
    g=-e5-e6+e8+e11.

Their Gram matrix is [[4,2],[2,4]]. The exact full-coordinate equalities
Lf=6f and Lg=6g hold, with zero exterior residual. Write D for their span
and R for its orthogonal projector. R is only used in the proof and terminal
diagnostic. The selection controller contains no R, f, g or energy filter.

For one preparation epoch put U=U_prep(4), Q=PUP on the twelve-dimensional cell, and
sigma=I_cell/12. The controller applies U, records inside or escape without
resolving the inside coordinate, and replaces an escaped source by sigma.
The complete reduced map is

    Phi(rho)=Q rho Q* + Tr((I-Q*Q)rho) sigma.             (1)

The second term is a measure-and-prepare channel whose complete escape and
reinjection records are retained. It is not normalization of a success
branch. Equation(1) is trace preserving and positive on arbitrary source
and reference inputs. Its field location may be selected by an earlier
classical source record. An actual unit at any cell address is an admitted
initial condition. The known-total-one source bridge supplies such states;
the coherent address encoding is an explicit apparatus isometry.

## The dark sector is selected by the process

On D, U=e^(-4iE_log)I, where E_log=-1260 log(209/210). Thus D is reducing for Q,
Q*Q is identity on D and the escape effect annihilates D. Both coordinates
and off-diagonal density operators in D are preserved by(1). The stationary
subspace is a single eigenspace both for H_prep and for the count-clock
generator L, whose eigenvalue is6, even when its density is mixed.

Put B=I_cell-R. The adjoint channel is

    Phi*(A)=Q* A Q + Tr(A sigma)(I-Q*Q).                (2)

It maps positive operators to positive operators. Since D is reducing and
the escape effect annihilates D, Phi*^n(B) is supported in B for every n.
The independent checker bounds the exact ideal-source operator
K=Phi*^64(B) by comparing its stored raw-Q reconstruction with the full
source-word and spectral certificates. It retains binary64 multiplication
error and a separate conservative bound for the64 matrix iterations. A
strict bound ||K||<1/4 gives the Loewner inequality

    0 <= K <= B/4.

Positivity then gives Phi*^(64k)(B)<=4^(-k)B. Therefore every initial cell
density, including every basis unit and any state with an untouched quantum
reference, obeys

    Tr(B Phi^(64k)(rho)) <= 4^(-k).                    (3)

The numerical tail from a stored Q need not itself annihilate D exactly.
That numerical matrix is never iterated as if it had the exact support
property. The comparison first bounds the exact ideal-source tail; (3)
then uses exact invariance. Cross terms between D and its complement do not
enter either the population or the escape probability. A single-step
singular norm need not be strictly below1, which is why the certificate
uses64 complete epochs.

This is source-driven spectral purification from a coordinate unit or a
mixed seed. The monitor specifies a region, but cannot choose the two dark
modes or their nonzero common energy. All80 lower native modes are retained
in the source construction and are included in the contraction test. A
HOLD source with the identical monitor has Q=I: the seed remains I/12 with
dark weight1/6. Coordinate dephasing after each epoch gives a separate
complete control of the coherent cancellation mechanism.

The number80 has an exact independent proof in
../../reviews/2026-11-15_codex_rotation_energy_review/LOWER_MODES_RESULT.json.
The20 identical body Laplacians have eigenvalues0(1),5-sqrt(5)(3),6(5),
5+sqrt(5)(3). Their integer polynomial and traces are checked directly. The
30 exterior edges form a matching, so their positive Laplacian is at most2I.
Minmax puts the first80 eigenvalues at most7-sqrt(5)<6 and the remaining160
at least6. The pump certificate disposes of every component orthogonal to D;
it neither discards these lower modes nor pretends6 is the ground energy.

The exact selected sector is held by source U_prep(t) and by exp(-itL) after the monitor and
reinjection are disabled. Its exterior amplitude is zero at every real t.
Perturbations with support in the cell can be subjected to the same complete
controller, and(3) applies to the perturbed density. Arbitrary unrecorded
external disturbances and autonomous hardware error correction are not
assumptions of this assertion.

## Complete native apparatus and finite implementation

The coarse region predicate is the Boolean union of address00000xxx and
address000010xx. A native CCX/flip network computes this flag using eight
source bits, one pointer bit and four scratch bits, then uncomputes the
scratch. Its truth table covers all512 source/flag inputs. Full13-qubit
raw source-word witnesses and the arbitrary-input operator bound retain
dirty-scratch and processor failures. The native pointer reads the one
region flag; measuring every source-coordinate bit would destroy the
required inside coherence and is not used as the monitor.

On escape, eight native bit reset operations prepare the source blank.
Four native coins propose an address in0 through15. Values12 through15
reject the draw; sixteen draws retain timeout probability4^-16. Every
accepted address has the same probability. Native flip words route the
blank to that address, producing I/12 when the private selection record is
discarded. The record itself remains in the archive. The certificate
charges all21 source/flag/scratch reset bits every epoch, even if the branch
does not need all of them. Each bit keeps the inherited2678-attempt reset
contract and timeout35/2^31. Reset execution ancestry is stated separately
from the64 fresh selector calls and the fresh coherent source words.

Let delta_flow be the complete half-diamond error per source flow and
delta_app the complete monitor/reset/reinjection error per epoch. Contractivity
of the comparison channels and complete first-failure composition give

    delta_512 <=512(delta_flow+delta_app).             (4)

Equation(4) is the comparison between exact complete physical instruments.
Stored floating maps can have gain above1. The final receipt therefore also
uses x/(1-2x), x=512(delta_flow+delta_app), before composing the two pumps.
If step differences are bounded by a_j in half-diamond norm, the complete
operator norm gain is at most1+2a_j. Expanding products and using
prod(1+2a_j)<=1/(1-2 sum a_j) gives total half-diamond difference at most
S/(1-2S), where S=sum a_j<1/2. The final routed receipt applies this guard
again to the preparation, reader and route sum. This conservative double
guard includes tensor and downstream cross terms; no exact unit gain is
assumed for a stored numerical center.

Every term includes invalid labels, source-padding aborts, processor
complements, reset timeouts and selector timeouts. The actual compiled
map need not have an exact dark eigenspace. Combining(3) and(4), its full
output has selected-sector probability at least

    1-4^-8-delta_512                                  (5)

from every admitted cell density. The failure symbols are unresolved
outcomes, not assigned energy zero. Since0<=L<=12 on the valid source chart,
the unconditioned valid-output second moment about6 is bounded by
36(4^-8+delta_512). Conditional on a valid output, divide this bound by
its certified probability lower bound1-delta_512. Thus exact sharpness
belongs to the ideal source program, and the finite native implementation
has an explicit source-energy window. Observed small numerical leakage is
reported beside this complete bound, never in its place.

## The count clock fixes the generator and energy origin

The sealed count-clock packet in ../2026-11-15_codex_poisson_source_clock
executes242175 native calls, including the full source columns of both compact
one-charge sectors and all four products in their separated signed-pair sector.
For one source its transfer eigenvalue is1254/1260. For the signed pair it is
1248/1260. Every exterior and vacuum residual in these equations is zero.
The second cell uses its own native compact basis:

    r=-e12-e15+e17+e18,
    s=-e12-e15+e16+e19.

Blindly shifting the first-cell basis would touch an external seam and would
not satisfy the source equation. The second-cell preparation uses its actual
basis and a separately certified source subblock. From coordinate14, its
64-epoch tail bound is less than1/16, so the ideal512-epoch bad weight is at
most2^-32. Its native region predicate and all reinjection routes execute.

One count-clock tick has duration1/1260 and a native-coin Poisson(1) number
of source events. The exact ideal event mixture is

    S=exp(T-I),       H_clock=-1260 log S=1260(I-T).     (6)

This follows from the complete native event series; it does not insert a
Hamiltonian eigenvalue chosen from a measured mass. On the actual signed-pair
source use its full transient restriction and retained vacuum branch. The
same mixture gives its Dirichlet generator. The source census therefore
supplies6 on each single-charge compact sector and12 on the separated signed
pair. The nonlinear finite-attempt logarithm is not used as the holding-mode
energy. Its independently compiled comparison is retained below.

The finite clock has a128-bit count decision and32-event cap. It keeps every
unresolved count and address draw. Its complete time4 error is bounded by
5040 count-timeout probabilities plus161280 address-timeout probabilities;
no resolved sample is renormalized into an exact Poisson law. Equation(6)
is the exact ideal-clock theorem, with this explicit finite realization.

The physical zero-load reference is checked on all1260 native seam/coin
cases: every load and current remains zero. Thus T has eigenvalue1 on that
reference, S has eigenvalue1 and H_clock has eigenvalue0. This fixes a common
energy origin. It introduces no separately chosen constant per particle.
The energy-mode chart includes the physical vacuum as label240 and retains
labels241 through255 as invalid. The unit-only preparation mode retains its
original padding failures; a vacuum is never silently admitted to that mode.
A complete native address read distinguishes the vacuum reference from an
unresolved charged state. The idle interferometer branch is an apparatus
reference and is not declared to be the vacuum merely by naming it.

## A coupled source and energy instrument

The declared coherent holding mode implements exp(-it H_clock) on the measured
compact sectors. This continuation is a specified quantum mode of the same
source, realized by native M6 words. It is not the unread Poisson jump channel.
The operator to be realized is fixed by (6) and the exact native source census,
rather than by a target energy or a511keV input. The source's available native
control algebra supplies the implementation. Architectural uniqueness is not
assumed. All proper-code load/export and controller precision premises are
listed with the program.

A native coarse-sector instrument first tests the complete two-dimensional
sector. It uses two controlled rank-one reflections, an auxiliary flag and
the native pointer. It keeps coherent superpositions within the sector. Its
complement is an explicit unresolved output. Both cells' full12-column maps,
processor complements and source flag effects execute and are retained. This
measurement is a terminal energy assay, not the operation that selects the
state: the source pump reaches that sector before it is used.

The observed source flag controls a native two-bit phase coupling. For phase
angle theta, the actual I/J source word realizes

    exp(i theta I) exp(-i theta diag(1,1,1,0))
       =diag(1,1,1,exp(i theta))                       (7)

on the flag/reference proper code. The pair computes the conjunction of its
two source flags with native CCX, applies(7), and uncomputes the scratch.
Thus an off-sector input does not receive an energy-success record. The
complete cq instrument is the composition of the actual source flag maps,
classical record encoding, native flag-controlled word and native reference
pointer. Every final source flag, reference bit and scratch outcome remains
in the record. A scalar phase calibration without this source connection
is stored only as a calibration component.

The complete output is evaluated on the512-epoch source states, selected-sector
calibrations, HOLD's unchanged coordinate0, off-sector inputs, physical vacuum
and invalid padding. The latter inputs have energy-success probability below
the numerical1e-12 gate and retain unresolved/vacuum/invalid records. Their
finite complete-instrument error is also stated; no exact zero for a stored
floating matrix is substituted for that bound. The vacuum's native generator
zero is a separate source result.

The native X/Y pointer supplies a complex accepted fringe z. At duration
1/10, equal to126 count-clock ticks, the source degree bound implies
0<=H_single<=12 and0<=H_pair<=24, so every admissible phase magnitude is below
pi. The decoder is E=-atan2(Im z,Re z)/t. It uses that source-derived bound
and the measured records; it does not use an expected6 or12 to choose a
winding. The decoded values are6 and12 to the recorded numerical precision,
with pair/single ratio2. Time4 is an independent longer-duration consistency
check of the same instrument and clock.

The numerical central agreement is not the reported physical precision.
For complete output error epsilon and an acceptance lower bound p_min,
each conditional X/Y mean changes by at most2epsilon/(p_min-epsilon).
Consequently the complex fringe radius is at most
2sqrt(2)epsilon/(p_min-epsilon). For radius r<1 the angular error is at most
r/(1-r). The short-time interval determines a unique winding of the long
measurement; no expected6 or12 selects that winding. Dividing the long-phase
error by4 gives the final energy interval. Taking the worst numerator and
denominator endpoints propagates the ratio and the single511keV anchor.
The final routed ratio is enclosed in
[1.999883624587696,2.0001163823007317], giving approximately
[1021.94053,1022.05948]keV. The old pair excess is outside this finite interval.

The single permitted511keV anchor converts the measured single-source energy
unit. The separate native signed-pair generator and coupled reader then give
1022keV in the exact ideal construction. This is the second job of the anchor;
no second scale or pair fit enters. The pair value is not generated by setting
E_pair=2 E_single in the reader. It is compiled from the separately enumerated
signed-pair source generator. Under the same source sector, source flags and
reader, the finite-attempt logarithmic clock yields a time4 nonadditivity
fringe of about0.1153193. Its original2.45085keV excess remains in the record.
The changed count clock, rather than a hidden recalibration, removes that
finite-attempt composition effect.

## The same prepared source and complete error budget

For routed coordinate0 in the first cell, R rho_initial=0. Since Q reduces D
and reinjection I/12 is block diagonal, every ideal preparation output has

    rho_N=(1-w_N) R/2 + sigma_bad,    Tr sigma_bad=w_N.  (8)

The second cell has the same property from coordinate14. Thus its trace
distance to the target sector density is w_N, without a square-root loss.
For arbitrary admitted inputs a gentle-projection comparison remains valid;
(8) is the stronger bound for the actual routed inputs. Native charge
conjugation maps the second preparation to the negative-charge sector with
coin/current relabelling. The separate negative-unit handoff and signed-pair
loader certify that source connection; a positive preparation is never
silently treated as the signed pair.

The loader and current-controlled routing packets are respectively
../2026-11-15_codex_signed_pair_loading and
../2026-11-16_codex_signed_pair_routing. All12 actual native continuations
end at delta0-delta14, with real0-to1-to0 current records defining the frame,
without a source reset. Their failure alphabet includes route timeout; its
per-trial bound1-(1-2^-32)^8 is included. The separate address memories are
prepared by the two unit-source pump programs and then exported to the
disjoint signed-source chart p*256+n. This is the declared loader boundary:
the two global unit pumps are not claimed to be the interacting signed-pair
evolution during preparation. Accepted loading intertwines the actual native
signed-source action, and all contact, outside, collision and vacuum types
remain in its separate source certificate.

The nonrouted all12 actual-input energy records also remain in the coupled
receipt. Equation(8) does not apply to their arbitrary source coordinates.
For those inputs the energy-output bound follows instead because the ideal
coarse D measurement removes D-to-complement coherences and every D state
has the same measured phase. Its entire classical output law differs from
the accepted ideal energy law by at most the bad-sector probability. This
argument does not assert trace distance w to R/2 for an arbitrary input.

The exact largest measured full6x6 primitive error is enclosed by
1.411037534545332e-15. The3,750,720-pulse preparation word has complete
half-diamond bound9.849889300102983e-8, or5.043143321652727e-5 for512 uses.
This is an operator bound from the executed primitive catalogue, with all
padding and processor outcomes, not a small state residual. Each cell's
512-step region/reset/reinjection apparatus contributes about0.0001825.
Together with(8), the clock caps and the complete sector/flag/phase reader,
the routed joint source-preparation and energy-output comparison, including
the final global gain guard, is at most0.0004841318708469099.
The component rational terms and hashes are in `coupled_energy/RECEIPT.json`;
the route and final global-gain terms are in `routed_join/RECEIPT.json`.

The output records retain upstream sector-processor failures before the
flag-controlled reader, downstream failures, and the separate arithmetic
normalization residual. Pair upstream first failure is f0+s0*f1, with s0
the first source's successful sector-instrument mass. Signed numerical
centers as small as-2.56e-13 are preserved without clipping. They are centers
of a certified comparison, not physical probabilities. The exact complete
instrument has a positive normalized output law within the stated radius.
On the unresolved branch the actual native pointer measures all eight
source-address bits; the resulting address distinguishes charged complement,
physical vacuum240 and invalid241..255. These are measured outcomes, not
classes assigned from the input's filename.

The raw floating success maps are kept without state renormalization. A tight
1e-7 diagnostic on their accumulated unconditional trace failed after512
uses; their one-word trace drift predicts a deficit of about1.2e-7. Its
original inputs and failure are retained. A separately frozen replay records
the noncompact population and numerical trace deficit separately. The exact
source theorem and the complete operator error budget carry the physical
process comparison; they do not infer a negative source probability from
rounding, or hide the failed numerical diagnostic.

## What the energy interpretation asserts

Within this declared finite construction, rest energy is the positive native
count-clock generator above the verified zero-load reference, read by a
source-coupled coherent instrument on a localized, source-selected state. The
same operator holds that state and supplies the separated signed-pair test.
This gives a concrete operational meaning to the calibrated quantity. A
transfer eigenvalue, a standalone programmed phase or a modular record cost
by itself would not provide that meaning.

Preparation, sector diagnosis, reset, reference preparation and coherent
control are apparatus operations. Their native event/pulse counts, memory,
timeouts and complete output errors are included. They are not claimed to
perform zero work, nor are their record counts converted into joules. Both
sides of the separated-source comparison use the same declared apparatus
and source clock. The quantity reported is the held source's excitation
energy relative to its common reference; it is not a thermodynamic balance
sheet for the complete measuring device. No inertial, relativistic or
laboratory proper-time identity is inferred from this static experiment.
The allowed anchor supplies the laboratory energy unit.

This theorem is a candidate-level construction and a quantified native
implementation. Checklist closure, source-pin adoption and official row
acceptance are recorded by their actual reviewers and owners. Same-group
execution and review confer no other-group G6 or maintainer award.
