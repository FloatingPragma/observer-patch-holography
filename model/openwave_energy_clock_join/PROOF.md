# One transfer-energy scale does not predict two reflection frequencies

## Source and claim

The finite diagnostic `static-two-mode-arc-diagnostic-v1` consumes the
240-port, 630-seam graph at
`code/openwave_static_signed_source/evidence/openwave/shared/localized_rest_join/evidence/graph.json`,
SHA-256 `2e8be4c030e8e04ab903adb5a8a786f11654ecca0380ee54d17509a707f49085`,
in OPH commit `b37ba4dae0e069dc030ceb0f130e52905fbb5465`.
The calculator reads its seams and port count, and reconstructs the endpoint
law, transfer and modes. It does not consume the file's stored rays or
endpoint table, another producer, a spectrum or a success receipt.

Its parent is `A-signed-compact-arc-clock-v1`, described by
[partial_records](../../code/openwave_partial_records/PROOF.md#row-3-de-broglie-clock).
The diagnostic adds a second mode and specifies the same ideal preparation
recipe for both. It is an apparatus extension, with no claim that its added
preparation is compiled by that parent. The official M10.2 candidate
`A-static-signed-source-v1` at `c330c6e7f23506b5f87c339cd73906928c1a7cd9`
is a separate contract. No future candidate selection is inferred.

The repository's official executable model, [CCF-LT-2](../SPEC.md), has a
distinct source and Hamiltonian. The spectral and instrument statements here
apply only to the declared static diagnostic; no source or observation map
to CCF-LT-2 is supplied.

The tested hypothesis is the positive-branch relation

    E_i = k theta_i,   k = E_f/theta_f,

with one zero reference, one scale and two distinct nonzero modes. A different
scale on the comparison mode is a recalibration. The constant graph mode has
E=theta=0 and fixes an additive offset to zero. Fitting two independent affine
parameters to the two nonzero modes is outside this hypothesis.

**Result:** this scalar relation fails on the two specified modes. The exact
nonlinear relation `E=-1260 log(cos(theta))` holds on their positive-transfer
branch. A physical energy-clock relation on the original apparatus is not
established by either statement. Review status is provisional; the internal
opposition calculation has no independent reviewer.

## Native transfer, memory and padding

One positive unit at p undergoes the native integer endpoint rule

    d=q_a-q_b,
    j=(d-(2c-1)(abs(d) mod 2))/2,
    q'=q-j e_a+j e_b.

A fair seam and coin have 1260 alternatives, including held events. Exactly
one coin moves the unit along each incident seam. With graph Laplacian L,
the mean one-attempt transfer is `T=I-L/1260`. The edge-square formula gives
`0<=L<=12I`, since maximum degree is six. Thus T is positive on the entire
240-position domain. Energy means the declared transfer decoder
`H_transfer=-1260 log T`, using the real logarithm.

The coherent register has 256 position addresses and 2048 event addresses,
dimension 524288. On physical positions define

    (Af)(p,e)=f_p/sqrt(1260),   e<1260,

with zero on all other coordinates. The disjoint columns have unit norm,
so `A*A=I_240`. The permutation S sends a moving arc to the reverse endpoint
with complementary coin, and fixes held and padding arcs. Direct endpoint
enumeration gives `S²=I` on all 524288 columns and `A*SA=T` on all 240 input
columns. The checker derives those counts independently of L.

R is `2AA*-I` on valid positions, including their event padding, direct sum
identity on all invalid positions. This is the extension used in the parent
recompute. R is an involution; `W=SR` is unitary. Physical event arcs form a
reducing subspace. On event padding W is -I; on invalid positions W is I.
Those conventions do not extend the admitted preparation domain.

W is a declared coherent memory word. It is distinct from the committed
repair instrument and from `exp(-it H_transfer)`. No passive-dynamics or
architecture-uniqueness claim is required or made.

## Two actual modes

The integer vector f has entries `(1,-1,-1,1)` at ports `(4,5,6,7)` and zero
elsewhere. Its full-graph equation is `Lf=6f`, with `||f||²=4`.

For the second mode put `a=5+sqrt(5)` and let B be any basis of the
three-dimensional a-eigenspace of the first cell Laplacian. Each of the
twenty cells has three ports on external seams. For its ordered port set b,
define its profile by

    g_cell = B (B[b,:])^-1 (1,1,1)^T.

The four boundary sets are `(1,3,10)`, `(2,8,10)`, `(0,8,10)` and `(0,2,8)`.
Every boundary matrix is invertible in Q(sqrt(5)). For example, the first
profile is

    (-1,1,-1,1,r-2,r-2,2-r,2-r,2-r,-1,1,r-2),  r=sqrt(5).

At both ends of every external seam g equals one. Hence the external
Laplacian annihilates g, and the internal equations prove `Lg=a g` on all
240 positions. In particular `g_0=-1`, `g_5=sqrt(5)-2`, `sum g=0`, and
`||g||²=880-288sqrt(5)>0`. The script verifies the full equation exactly.
This mode is extended over twenty cells; no second localized particle is
asserted. Normalization changes neither eigenvalue nor spectral decoder.

The corresponding transfer multipliers are

    lambda_f=209/210,
    lambda_g=1-(5+sqrt(5))/1260.

They belong to this full graph, rather than a restricted-cell scheduler.

## Same reference and complete instrument

For each normalized position mode f_i set `u_i=Af_i`,
`s_i=sqrt(1-lambda_i²)` and `w_i=(Su_i-lambda_i u_i)/s_i`.
Compression gives an orthonormal invariant pair with

    W u_i=lambda_i u_i+s_i w_i,
    W w_i=-s_i u_i+lambda_i w_i.

The oriented triangle `(0,1,5)` defines v: +1/sqrt(6) on its three forward
moving arcs and -1/sqrt(6) on their reversals. Port sums vanish and `Sv=-v`,
so `A*v=0`, `Wv=v` and v is orthogonal to both invariant pairs. The same v
is used for both modes, without a phase drive between words.

The supplied initial density is `tau=I_body/302400`, where body means all
physical position/event arcs, including native holds. For the requested mode
i, one ideal filter has effects `P_chi_i` and `I-P_chi_i` on that body,
where `chi_i=(u_i+v)/sqrt(2)`. Its conditional successful density is
`P_chi_i`, with probability `p=1/302400`. Failure is terminal and retains the
subnormalized branch density `(I-P_chi_i)/302400`, whose trace is `1-p`.
Conditioned on failure, the normalized density is `(I-P_chi_i)/302399`.
The complete instrument retains the subnormalized branch and its original
weight. The filter is tried once.
The input, reference, preparation recipe and detector are fixed, with i a
mode selection. Implementing the second filter is a declared coherent service,
with no native-pulse compilation or physical fidelity certificate.

On success apply n complete W words, with requested integer `0<=n<=16`,
then the full coordinate PVM on all 524288 addresses. The complete alphabet
contains one preparation-failure flag and every coordinate z. Its law is

    Pr(fail)=1-p,
    Pr(z)=p/2 |cos(n theta_i) u_i(z)+sin(n theta_i) w_i(z)+v(z)|²,
    theta_i=acos(lambda_i) in (0,pi/2).

These are original-input probabilities; their sum is one. The failure source
is retained and receives no assigned energy. The bounded ideal sequence has
no unresolved waiting branch. No stochastic retry or successful-only timing
is substituted. A request beyond the specified word range is outside this
protocol; no hardware deadline law is supplied.

At the same arc `z=(5,3)`, `v(z)=1/sqrt(6)`. Both u_i(z) are nonzero, so the
cos(theta_i) interference coefficient `v(z)u_i(z)` is nonzero. Without v,
the u-return probability is `cos²(n theta_i)` and has harmonic 2theta_i.
The measured menu here includes the stationary reference and exposes theta_i.
An integer-word 2pi alias is unobservable; the principal branch and positive
real log domain specify which frequency-energy comparison is tested.

The transfer energy labels the underlying position mode f_i. It is not an
energy measurement on chi_i or a transfer-Hamiltonian assignment to v.
The comparison tests a proposed conversion of this label into the declared
arc observation; it does not identify the two physical processes by definition.
The complete-output claim is for these fixed preparations. It is not a
channel-distance estimate for arbitrary inputs or external references.

## Certified discrepancy and observable witness

Exact rational series bounds give the following outward decimal enclosures:

| Quantity | Enclosure |
| --- | --- |
| E_f | [6.014331228348, 6.014331228349] |
| theta_f | [0.097628775039, 0.097628775040] |
| E_g | [7.256925920238, 7.256925920239] |
| theta_g | [0.107223284457, 0.107223284458] |
| k | [61.604083692711, 61.604083692712] |
| k theta_g | [6.605392189505, 6.605392189506] |
| E_g-k theta_g | [0.651533730732, 0.651533730733] |
| E_g/k, the predicted second phase | [0.117799429603, 0.117799429604] |

The code emits finer rational endpoints. For an exact explanation of the
sign, `(-log cos theta)/theta` is strictly increasing on `(0,pi/2)`:
its derivative has numerator `theta tan theta+log cos theta`, zero at zero
and with derivative `theta sec² theta>0`. Since theta_g>theta_f, the anchored
linear conversion necessarily underpredicts E_g.

The discrepancy affects a complete observable, beyond assigning two spectral
numbers. Replace theta_g by the single-scale prediction E_g/k in the same
invariant pair, preparation and coordinate PVM. At n=1, z=(5,3), actual minus
predicted original-input probability lies strictly in

    [2.994e-10, 2.995e-10].

The factor 1/302400 is included. The calculator encloses both amplitudes
using the exact g profile and norm, rational square roots and alternating
sine/cosine bounds. This is a probability-law distinction; no finite-shot
measurement campaign or sample sufficiency is inferred from its small size.

For logarithms the positive atanh series has geometric tail
`2 z^(2N+1)/((2N+1)(1-z²))`, `z=(1-x)/(1+x)`. For acos the positive arcsine
series uses `acos x=2 asin sqrt((1-x)/2)` and tail at most
`2 y^(2N+1)/(1-y²)`. Integer square roots and rational outward rounding bound
all inputs. Increasing floating precision is not the certification method.

## Controls, execution and boundary

The complete-law predicate accepts literal W execution against the analytic
law on all coordinates for each n from zero through sixteen. It rejects a
doubled frequency, reversed phase, dephased or removed reference, quarter-turn
reference phase and removal of the preparation failure. Replacing the second
scale by its separately fitted value fails the shared-parameter predicate.
Dropping R from the word changes the same full-law predicate after two words.
Complex logarithm branches fail the real-energy domain predicate; padding
inputs fail the preparation-domain predicate. Changing an internal seam in
one cell is rejected by the structural requirement of identical ordered
internal cells. Moving an external seam from `(1,82)` to `(4,82)` preserves
those cells and passes source admission, then fails the actual compact
anchor eigen-equation. Both tests require their specific rejection message;
neither uses a stored digest rejection.

Positive controls retain two orthogonal L=6 compact rays under the same
calibration, the constant zero mode and the exact compressed Chebyshev identity
`A*W²A=2T²-I`. The known nonlinear relation is compatible with both modes.

A 2pi phase alias is explicitly non-discriminant. Relabeling one W word as
one native attempt cannot be detected from these probability laws alone:
a native implementation/resource-duration map is absent. The executed service
counts are one filter attempt, n reflections, n reversals and one coordinate
measurement on success, with zero words on preparation failure. No duration,
native attempt count or primitive pulse cost is assigned to those services.
Substituting a continuous Dirac source for an ordered macro belongs to the
[separate exact counterexample](../openwave_macro_generator_identity/PROOF.md),
with a different source; that control is not imported into this diagnostic.

`check.py` uses a 524288-entry permutation and vectors, plus the 240-position
integer count matrix and exact cell equations. It constructs no dense arc
operator. Exact compression and the reducing subspaces prove all integer
words; floating array comparisons for n=0..16 are software diagnostics.
`crosscheck.py` reconstructs the second mode through a cell spectral
projector and an edgewise full-graph action. It encloses log by its power
series and acos by rational cosine bisection, importing neither calculator
nor bound kernel. These are different implementations by the same author,
sharing conceptual source dependencies. No independent review is claimed.

The [WQ-03 lepton timeline](../lepton_clock_energy/PROOF.md),
[elementary density/wave comparison](../observer_time_fields/REFINEMENT.md),
[clock/reference tests](../clock_defect_energy/CLOCK_WAVE.md) and
[typed chain checks](../../code/openwave_chain_checks/provenance/PROOF.md)
have their own source and observation domains. Their results supply no silent
transfer to this graph apparatus. The finite scalar counterexample implies
no OPH-wide obstruction, physical mass prediction, OpenWave cell award,
absolute-unit defect or closure of a physical completion gate.
