# Recovery evidence must describe the supplied state and counts

This repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
was developed on main `3d848b7f` and rebased onto `555852f8` after #1050/#1051
merged. It extends the finite-state validation
work and the squared-fidelity correction in #1038 to the remaining Stage 1
recovery and tomography pipeline. It supplies no physical model, recovery
instrument, or evidence distinguishing OPH from ordinary quantum mechanics.

## Reproduced defects

The first commit, `4d129953`, records **27 failing cases and four passing
existing controls**, before any implementation change.

| Input | Previous result | Correct treatment |
| --- | --- | --- |
| Zero matrix compared with `I/8` | Fidelity 1, obtained by manufacturing `I/8` | Reject: a zero matrix is not a normalized state |
| `diag(1-epsilon,0,0,0,0,0,0,epsilon)`, `epsilon=2^-40` | The recovered rare sector was zero | Preserve both resolved sectors; this classical Markov state recovers exactly |
| `diag(1,1e-100)` within the declared trace tolerance | Zero entropy | Retain the positive eigenvalue contribution, without normalizing the input |
| Empty tomography data | `I/8`, with apparently perfect recovery | Reject missing measurements |
| Empty, negative, fractional, Boolean or malformed counts | Zero or otherwise plausible expectations | Require matching binary outcomes and nonnegative integer counts, with positive shots per setting |
| Three estimates of XI with signed shot totals `(1,-100,1)` and shot counts `(1,100,1)` | `1/3` | The pooled estimate is `-98/102` |

The old inverse-root cutoff `1e-10` and subsequent output normalization were
separate errors. For the Markov example with `0<epsilon<=1e-10`, deleting the
rare B direction leaves a trace-deficient recovered operator with trace
`1-epsilon`. Renormalizing it yields `|000><000|`. Its trace distance from the
true state is epsilon, its squared fidelity is `1-epsilon`, and the Z read
changes by `2 epsilon`. A state normalization check alone cannot detect this
loss. The conditional sector is erased completely.

## One finite channel, with its support domain made explicit

Let tau be a positive operator on BC with positive trace, `t=Tr_C(tau)`, and Pi the
support projector of t. Let T be the Moore--Penrose inverse square root of t,
with every positive eigenvalue inverted. Define

```text
K_c = sqrt(tau) (I_B tensor |c>) T,
R_0(X) = sum_c K_c X K_c^*.
```

Direct multiplication gives `sum_c K_c^* K_c = T t T = Pi`. Thus R_0 is
completely positive and trace preserving on the supported input algebra;
the pseudoinverse alone does not make it trace preserving on all B inputs.

The new `petz_recovery_kraus` adds the explicit completion

```text
R(X) = R_0(X) + Tr[(I-Pi) X] tau/Tr(tau).
```

For an orthonormal kernel basis `{v_r}`, its additional Kraus operators are
`L_jr = sqrt(tau/Tr(tau))|j><v_r|`. Their adjoint products sum to `I-Pi`.
This agrees with replacement by tau for a normalized reference and preserves
unit trace on kernel inputs for accepted source trace drift. Consequently the
exact fixed map is completely positive and
trace preserving on the full B algebra. There is no division by the trace
of its output and no selection depending on the state being processed.
The replacement on the unused kernel is a disclosed completion choice,
not a derived physical law.

Positivity and `Tr[t(I-Pi)]=0` imply that tau is supported in
`Pi tensor I_C`: a positive operator with zero trace on a subspace
annihilates that subspace. Hence `R(t)=tau`. If rho_AB has marginal t,
the same support argument shows that its kernel-completion term vanishes.
Thus, for normalized rho_ABC, using tau=rho_BC recovers a normalized ABC state with the original
A marginal and BC marginal, and agrees with the usual supported Petz formula.
This proves the exact mathematical statements for finite states, including
singular references. The executable implementation evaluates them numerically.

The tests check every Choi entry through matrix-unit action, trace preservation,
reference recovery, arbitrary positive inputs on the completed kernel, and
65-digit independent matrix-function evaluations in dimensions 2x2, 2x3 and
3x2. A symbolic Bell-reference calculation proves that its channel replaces
any B input X with `Tr(X) P_Bell`; a valid recovery channel need not leave its
input marginal unchanged.

The general channel and fidelity definitions are the standard ones in
[Watrous, chapters 2 and 3](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
They are used here as software controls, not new OPH assumptions.

## Complete classical conclusion, and the quantum distinction

For any classical distribution p_abc, the recovered distribution is

```text
q_abc = p_ab p_bc / p_b     when p_b > 0,
q_abc = 0                 when p_b = 0.
```

Summing the expression proves normalization and preservation of both AB
and BC. It is Markov, and a second reconstruction gives q again. For any
candidate Markov distribution `r_abc=s_b u_a|b v_c|b`, expansion of the
logarithm gives

```text
D(p || r) = I_p(A:C|B) + D(p_B || s_B)
            + sum_b p_b D(p_A|b || u_A|b)
            + sum_b p_b D(p_C|b || v_C|b).
```

The identity uses the usual support convention: a positive numerator over a
zero denominator gives infinity, and zero-p_b terms contribute zero.
Nonnegativity of the three added divergences proves that q minimizes
`D(p||r)` over the entire classical Markov family, with minimum
`D(p||q)=I_p(A:C|B)`. Equality fixes the candidate distribution q, including
its zero-weight sectors. This does not turn the generally nonconvex Markov
family into an affine constraint set. Nor is optimization over the second
argument of D the A3 optimization over its first argument.

An exact rational distribution is constructed on **each of the 255 nonempty
supports** of three classical bits. Independent rational marginal arithmetic
checks the recovered q; high-precision scalar logarithms check the entropy
identity. Separate Markov controls retain rare masses down to `1e-310`.

The corresponding quantum map need not preserve AB, be idempotent, or return
a Markov state. Noncommuting controls retain all three failures while checking
the A and BC preservation identities. The Fawzi--Renner theorem guarantees
the existence of a suitable recovery channel, not a fidelity promise for this
particular unrotated map. Its squared-fidelity bound in bits is `2^-I`.
[Fawzi and Renner](https://arxiv.org/abs/1410.0664)

The regression also replays the full-rank example obtained from Table I and
equation (20) of [Zhang's published counterexample](https://yuxuanzhang1995.github.io/assets/pdf/agentic/petz-cmi.pdf),
with maximally mixed admixture `10^-4`. A separate 65-digit implementation
reconstructs the exact decimal data, contracts indices and evaluates the
matrix functions. It gives `I + log2(F_Petz)` approximately `-0.00180240517`
bits. The ordinary Petz fidelity is below the optimal-recovery bound even
for this strictly positive state. The benchmark preserves that result rather
than adjusting its fidelity to satisfy an inapplicable inequality. This is
an attributed external control, not a new OPH result; the regression is a
high-precision numerical replay, not a replacement for that work's certificate.

## Missing measurements can hide one full bit

The tomography contract is a complete set of local Pauli settings, with
strict binary outcomes in Qiskit's count order. Repeated estimates of one
observable pool signed counts and shots before division. With equal shots,
this is the old mean; unequal shots retain their actual weight.

Completeness is essential for this unconstrained linear inversion. The states

```text
rho_0 = I/8,
rho_1 = (I + Z tensor Z tensor Z)/8
```

give identical uniform outcome probabilities in every local-Pauli setting
except ZZZ. Nevertheless `I(A:C|B)` is zero for rho_0 and one bit for rho_1.
For rho_1 every proper marginal is maximally mixed, while its full entropy
is two bits. An independent product-projector calculation checks all 26
remaining settings. Filling the missing ZZZ expectation with zero therefore
chooses a Markov answer that the data cannot distinguish from a non-Markov
state. Rejecting incomplete data closes that ambiguity; it does not imply
that tomography with a separately declared restricted model is impossible.

Shot noise can make complete linear inversion nonpositive. The existing
normalized-positive-part estimator remains an explicit tomography operation.
It is neither a channel nor a claim to be the nearest trace-one PSD matrix.
New reports expose its raw minimum eigenvalue, removed negative spectral
mass, Frobenius correction, per-setting shots and estimator name. They also
retain the complete count dictionaries for later replay. They do not report
those corrections as confidence intervals or certify statistical error.
The analysis and state-distance APIs never call this estimator.

### Preserve positivity through the returned matrix

The cross-PR review exposed an integration defect with the exact information
validator in #1049. Four independent three-qubit Born-count experiments,
each with 8192 shots per Pauli setting, produced rounded estimates with
negative eigenvalues near `1e-17` after the nominal positive-part projection.
Exact characteristic coefficients confirmed that the returned entries were
indefinite. The information API correctly refused them. Its exact-state
contract must not be weakened to conceal this estimator defect.

The estimator tests its returned Hermitian entries by rational Schur
elimination. For a positive pivot, positivity is equivalent to positivity
of its Schur complement; a zero diagonal requires the corresponding row
to vanish. This handles singular states without a tolerance. A candidate
that passes is retained, including exact rank-deficient candidates.

When ordinary reconstruction loses positivity, form the numerical factor
`B = V sqrt(p)` from the normalized positive spectral part. Treat the supplied
entries of B as exact binary rationals, and compute `G = B B*` exactly.
Round each off-diagonal to binary64 once, with conjugate symmetry. Let
`e_ij` be that exact rounding error. Choose each output diagonal by upward
rounding of

```text
G_ii + sum_(j != i) (|Re e_ij| + |Im e_ij|).
```

For the resulting M, `E=M-G` is Hermitian diagonally dominant with nonnegative
diagonal, hence positive semidefinite. Thus M is positive semidefinite as
the actual returned binary64 matrix. Upward rounding of `Tr(E)` bounds both
its trace norm and operator norm. There is no fixed positive floor. This
certificate concerns rounding the numerical factor; it does not bound the
factor's eigensolver error or experimental statistical uncertainty.

The final state retains the shared trace-one tolerance (`1e-12`); no extra
normalization follows the certified rounding. Normalizing those entries
again would require a new positivity check. `gram_rounding_trace_bound` is
null when the original candidate passes exact positivity, and is the bound
for the fallback Gram reconstruction otherwise. The report includes the
normalization tolerance and the complete correction from raw inversion.
Exact characteristic-polynomial controls independently verify both M and
`M-B B*` across real/complex factors and extreme scales. Complete noisy-count
replay exercises reconstruction, information and recovery together.

The audit implication is specific: ideal Born-probability reconstruction
does not exercise finite-shot spectral projection. Retain complete sampled
counts through the full analysis path, combine stricter downstream validators
with their producers, and certify properties of the returned entries rather
than the eigenspectrum intended before reconstruction.
The same count fixture can require outward rounding on one LAPACK build and
already be exactly positive on another. Controls assert exact positivity in
both cases and check any emitted bound; they do not require the fallback to
run when its precondition is absent. The first Ubuntu CI replay exposed and
corrected an overly strict assertion about that optional bound.

### Circuit states must satisfy the same contract

A further audit with the actual Qiskit 2.5.2 SDK found that the ordinary
`DensityMatrix(Statevector(...))` constructor also loses exact positivity.
With #1049's strict information validator, the structured theta=.6 and 1
circuits and depth-three random circuits with seeds 0 through 3 all failed
before reporting CMI. The theta=0 and GHZ controls alone did not expose it.
Six retained SDK-array-contract regressions fail on `72ce9e1d`; they check
the returned entries independently using characteristic coefficients.

The circuit constructor now reverses the statevector tensor axes into
q0-first order, then applies the same outward Gram rounding directly to
that numerical vector. It does not project or renormalize a supplied state.
The bound is relative to the numerical amplitudes, and does not certify
simulation error relative to ideal gates. Exactly representable products
receive no artificial floor. Invalid, non-normalized or incorrectly sized
amplitudes fail explicitly.

Every circuit-reference report and random-control candidate now retains a
`circuit_state` record with the construction, rounding bound, scope and trace
tolerance. Legacy fields `exact_analysis` and `exact_cmi_bits` remain for
compatibility: they denote the numerical circuit reference, not an exact
ideal-circuit certificate. The tests replay the complete JSON report path
and check that these qualifications survive serialization. The real-SDK
replay also covers all ten default random candidates and the five selected
catalog circuits; their complete CMI/recovery analysis succeeds with #1049.

This adds an audit requirement to the preceding lesson: follow every state
producer into the consumer, including optional SDK adapters and the actual
command-line report. A successful tomography reconstruction test does not
validate the separate circuit-reference constructor.

### Unresolved recovery must not erase acquired evidence

The maintainer's review of `d5e30c7c` found a separate reporting defect:
strict support refusal aborted the run before any counts were saved. The
complete 27-setting fixture assigns outcome `k` the count
`{'X': (8,2), 'Y': (9,1), 'Z': (5,5)}[basis[1]][(k>>1)&1]`.
Each setting has 40 shots. Reconstruction returns an exactly PSD state near
`I_A/2 tensor rho_B tensor I_C/2`, with B Bloch vector `(0.6,0.8,0)` and
CMI zero. Its roundoff-scale dense support remains unresolved for inversion.
Qiskit 2.5.2's depth-three seed-1 reference independently reaches that same
refusal. Neither result justifies removing positive spectral mass.

The low-level inverse still refuses that support, now with the specific
`UnresolvedPetzSupport` subclass of `ValueError`. `analyze_state` first
computes and validates CMI, then handles recognized recovery precision errors.
The report retains CMI and the optimal-recovery bound, sets `petz_status`
to `unresolved_support`, records the reason, and serializes fidelity,
trace distance and observable mismatch as JSON `null`. Available recovery
has status `available` and a null reason. The recovery-ordering check is
null if any of its three required fidelities is unavailable; the two
CMI comparisons remain evaluated. Unavailable is neither pass nor fail.
Invalid states, invalid CMI and unrelated analysis errors still raise.

Immediately after the sampler returns, before count regrouping or reference
and tomography analysis, `acquired_counts.json` saves all returned counts,
the state/basis-to-circuit mapping, catalog and run metadata. It can be
replayed even if subsequent analysis raises. Successful runs also retain
counts in `summary.json` as before. Circuit candidate CMI is necessarily
computed before acquisition to select the random control; it does not
require Petz inversion. This change neither submits a cloud job nor alters
the numerical support, information or state-construction contracts.

Seven failing report/evidence regressions were retained in `c1459805`
before repair. Mandatory tests include the exact count fixture, the
unmodified seed-1 numerical amplitudes, saved-count replay and forced
post-acquisition errors. Two optional cases use actual Qiskit circuits and
CLI parsing with `--random-seeds 1`, replacing only the sampler with a
local count fixture. The fixtures test reporting, not hardware behavior.

The custody audit of `431556fe` also reproduced a sequential-rerun defect:
after a successful 135-record run, incomplete acquisition in the same output
directory left 134 replacement count records beside the old successful
summary. Four failing controls are retained in `49ccbf69`. The runner now
refuses existing Stage 1 output files before sampling; it neither overwrites
nor removes them. Separate runs require separate output locations. Fresh
directories retain the post-acquisition error behavior described above.

Additional independent controls retain nonzero information through unresolved
recovery: a shared classical AC bit and a Bell AC pair, each independent of
a dense B state with eigenvalue `2^-50`, have CMI one and two bits. Their
reports retain those values and bounds `1/2` and `1/4` while recovery is
unavailable. Complete integer-count fixtures with AC correlations
`c=0,0.6,1` independently give `I=1-h2((1+c)/2)` and Petz fidelity
`(1+sqrt(1-c*c))/2`; permuting these fixtures checks both true and false
fingerprint outcomes. Zero-only unavailable cases cannot detect an
implementation that erases valid CMI, and equal-fidelity cases cannot
validate a nontrivial ordering.
Three injected faults passed the former SDK-free Stage 1 test file at
`431556fe` (79 passes, two optional skips): zeroing unavailable CMI and its
bound, forcing the recovery ordering to succeed, and forcing both CMI
comparisons to fail. The strengthened controls reject each fault. Removing
the output preflight also fails all four custody controls. The twelve
mathematical and seven earlier reporting mutations were rerun and rejected.

### Source support, square roots and reported distances

The follow-up reviewed main `84025a97491f6e1173cc29d72822c989e4ef6ed5`.
Passing trace preservation did not establish correct channel action:

| Original input and independent result | Defective result |
| --- | --- |
| `P=vv*/4`, `v=(1,i,-1,-i)`: `P^2=P`, so `sqrt(P)=P` and `F(P,I/4)=1/4`. | Fidelity `0.2500000061422713`; spurious positive eigenvalues were amplified by square roots. |
| `P=vv*/4`, `v=(1,i,1,-i)`: `Tr_C(P)=I/2`, so the Petz map is `X -> Tr(X)P`. | Action on the first-coordinate projector wrong by Frobenius norm `9.26e-9`, despite passing trace preservation. |
| `tau=[[.5,.3,0],[.3,.18,0],[0,0,.32]]`: the exact leading minor of the supplied binary64 entries is positive, hence rank three. With C dimension one its Petz channel is the identity. | The eigensolver and rounded kernel residual both indicated a false zero. A fabricated kernel completion changed the first-coordinate projector by norm `0.706341`, while remaining numerically trace preserving. |
| Two independent off-diagonal source entries equal to the least subnormal `s`, with their transposed partners zero and positive surrounding diagonal blocks. | Hermitian averaging before subtraction erased both `s/2` coherences; trace distance was zero instead of `s`. |

The shared exact pivoted LDL factorization determines positivity and rank of
the original Hermitian average before binary64 rounding. It also replaces the
separate Schur positivity checker used by tomography; that producer's outward
Gram-rounding contract is unchanged. An exactly indefinite source is invalid,
even below the density validator's historical negative-eigenvalue tolerance.

For `A=XX*`, `B=YY*`, squared fidelity is `||X*Y||_1^2`. The LDL factors give
the exact rational overlap and its rank before numerical singular values.
Rank zero proves zero fidelity; rank one gives exactly `Tr(AB)`. At rank two,
`F=t+2 sqrt(e2)`, where `t=||X*Y||_F^2` and `e2` is the sum of squared 2-by-2
minors. Cauchy--Binet expresses both in rational LDL columns and positive
pivots. Rational square roots give exact rounding, including both directions
of a subnormal tie. Higher ranks retain every certified positive singular
value and refine precision. No positive eigenvalue cutoff is introduced.

Roots retain exact support and refine the function before output conversion.
A zero component needs an exact certificate: vanishing components of
`A,...,A^n`, a verified positive root on an invariant block, or rational spectral
projectors with exact square-root cancellations. For the last certificate,
rational eigenvalues are grouped by rational square class; distinct classes
have linearly independent square roots over the rationals. Controls include
ordinary full-rank matrices whose irrational root has an exact zero, and nearby
nonzero components. Inverse roots retain the conservative dense-support guard,
additionally requiring the numerical positive count to match the source rank.
The false-kernel example raises `UnresolvedPetzSupport`; true singular
completion and tiny coordinate support remain usable.

Petz marginals and trace-distance differences retain original entries until
output conversion. The public recovered state evaluates the supported map
directly from the original marginals, before rounding its root and inverse
factors. Exact diagonal marginals use the rational conditional law. General
quantum marginals use refined matrix functions and contractions in a private
multiprecision context. Exact component zeros can be certified by the real
polynomial span of the two spectral functions: each unordered Hermitian pair
`V_i (AB tensor I_C) V_j* + V_j (AB tensor I_C) V_i*` must vanish in that
component, with diagonal pairs included once. Here the `V_i` span the lifted
products `BC^r (B^s tensor I_C)`. Both spectral functions vanish at zero,
so positive powers suffice. This sufficient certificate has a bounded span;
an unresolved cancellation raises rather than becoming a guessed zero.

Exact marginal products simplify the map before numerical evaluation:
if `BC Tr(rho)=B tensor C`, its output is `AB tensor C/Tr(rho)`;
if `AB Tr(rho)=A tensor B`, its output is `A tensor BC/Tr(rho)`.
These identities apply even when the original ABC state has correlations
absent from those two marginals. They cannot justify returning that input.
For a remaining numerical cancellation on a state of dimension at most 16,
the existing exact finite Markov criterion supplies a second certificate.
A proposed output is accepted only if it is exactly PSD, its exact AB and BC
marginals equal the original marginals, and it is exactly Markov. Such a state
is its own Petz recovery, and its recovery map and input marginal are the
specified ones; therefore it equals the required output. This proves the
entire candidate, including any provisionally zero components. Merely small
entries or matching marginals without the Markov certificate are insufficient.
The source itself can also be certified as an exact Markov fixed point when
the numerical candidate is unresolved. Both checks retain the
channel-domain checks. The dimension
bound limits the cost of the relative-modular moment calculation; it is not
a physical restriction or a tolerance for approximate Markov behavior.
Neighboring non-Markov controls must recover the marginal-determined state,
not inherit the fixed-point result. The independent audit retained full-rank
complex products and a nontrivial `AB_L tensor B_R C` Markov state that the
polynomial zero certificate alone could not resolve.

The separately returned binary64 Kraus operators are approximations to that
source map. A necessary additional channel
check compares its action on the reference marginal with each nonzero source
component on that component's own scale. For B dimension one and reference
`[[.5,i*s],[-i*s,.5]]`, correctly rounded root components still produce an
inaccurate replacement channel: its binary64 application doubles the tiny
coherence. This case now raises; the neighboring diagonal `(1,s)` replacement
remains available. Reference preservation and completeness are necessary
checks, not a proof of all channel action. Independent matrix-unit and
closed-form channel controls remain part of the tests.

The author audit reproduced a separate downstream error: with
`rho=diag(.25,.25,.25,.25,0,5s,0,0)`, exact conditional populations lie just
below and above `2.5s`, respectively. The old rounded-factor contraction
returned `(3s,3s)` instead of `(2s,3s)`. Independent rational block-conditional
controls include coherent real and complex BC references, so a diagonal-only
repair fails. Kernel replacement controls independently check both signs of
accepted trace drift while preserving the original supported reference.

If rounding the recovered matrix makes it exactly indefinite, recovery raises
`UnresolvedStatePrecision` instead of adding a positivity floor. Stage 1 catches
recognized precision failures from recovery and its three metrics, reports
`petz_status=unresolved_precision` (or `unresolved_support` for its subclass),
and retains CMI, its optimal-recovery bound and acquired counts. Invalid CMI
and unrelated exceptions still propagate. Stage 1 also passes original entries
to CMI: its former pre-averaging wrapper could erase a positive coherence and
report zero information. Controls retain representable tiny CMI, exact zero,
and unrepresentable positive CMI that must raise rather than become zero.

Stage 1 Pauli contractions use the original exact Hermitian average. Observable
mismatch forms differences before rounding either expectation, then averages
them before its single output conversion. For diagonal states with common
unit mass at 000 and an added mass m at 001 versus 010, four of 36 low-weight
Paulis differ by `2m`; the mean is `2m/9`. Separate float expectations erased
this positive signal for ordinary tiny m and for `m=9s`. A mean below the
representable range is unavailable, not zero. A complex dyadic control checks
all 64 Pauli labels against independent literal Kronecker matrices.

Additional exact paths resolve rational answers with irrational intermediate
spectra: arbitrary-rank diagonal fidelity uses the cross-root sum, and trace
distance uses one-sign rational characteristic factors and Vieta's formula.
Mixed-sign irreducible factors remain numerical. Normal and subnormal ties
exercise both rounding directions. Source Hermiticity and trace checks compare
exact squared defects against the shared tolerance, avoiding a preliminary
floating norm or eigensolve that can fail under the caller's underflow policy.

A source with a rank-one A or C marginal is an exact product across that cut,
so its Petz recovery equals itself. The public recovery uses this fixed-point
identity after the channel-domain checks, preventing factor rounding from
making a valid complex product state unavailable. Exact source rank selects
the identity; numerical near-purity does not. The neighboring pure state with
amplitudes `(1,1,1,1,1,1,1,i)/sqrt(8)` has CMI
`h2(1/2+sqrt(10)/8)=0.48376694420142585` bits and remains a nontrivial available
recovery control. These are finite matrix identities, not a new physical
Markov premise or a source-selection claim for OPH.

The initial implementation audit rejected twelve isolated mutations through
retained scientific assertions:
old roots, old fidelity, missing inverse-rank validation, pre-averaged source
entries, constant-zero fidelity, blanket root refusal, missing exact spectral
zero certificates, missing reference-component checks, missing semidefinite
distance ties, pre-averaged CMI, an overly broad exception handler, and a
missing recovered-positivity guard. Positive controls pass in the same copied
packages; collection failures are excluded. This bounded set tests specific
failure mechanisms, not completeness of the numerical domain.

The author review adds twelve isolated mutations: floating prevalidation,
missing arbitrary-rank diagonal fidelity, checking only one state's
diagonality, missing source tolerances, rounded Pauli contraction, rounded
expectations before subtraction, rational-eigenvalue-only distance ties,
diagonal-only recovery repair, missing whole-output certificate, omitted
Markov identity, omitted exact marginal equality, and unnormalized kernel
replacement. Each fails a scientific assertion or a required valid result;
the copied unmodified controls pass. The coherent recovery mutation reports
`3s` instead of the independently derived `2s`. Independent dense quantum
oracles, invalid candidates, non-Markov neighbors and trace-drift controls
guard against a repair that merely recognizes the retained example.

## Numerical and downstream boundaries

Strict numeric conversion rejects masked data, mixed Booleans and lossy
binary64 conversions before the shared density validator. Normalization and
Hermiticity retain the shared `1e-12` acceptance convention; roots, distances
and recovery additionally require exact PSD of the supplied Hermitian average.
Inputs are not rescaled and positive eigenvalues are not floored.
Inverse roots reject unresolved dense support instead of guessing rank;
coordinate-diagonal tiny positive entries remain explicitly supported.
Non-rational matrix-function values remain numerical calculations, with
refinement through 1280/1320 decimal digits, support checks and output-scale
and midpoint guards. Exhausted precision and nonzero output underflow raise
`UnresolvedStatePrecision`; agreement is not an interval enclosure. These
diagnostics are not clipped to make them satisfy a theorem. A finite returned
channel still has binary64 factor error; independent action controls remain
essential alongside trace preservation.

The new code is shared with the finite quantum-information package. The
Stage 1 wrappers retain bit units and the legacy squared-fidelity result
fields. The existing optional-hardware imports remain lazy: numerical tests
need no Qiskit or provider access. No quantum jobs are submitted by this repair.

The frozen Stage 1 summaries contain neither complete counts nor reconstructed
states, so their measured recovery metrics cannot be recomputed from this
archive. Their bytes remain unchanged, and their reading notice records the
additional limitations. The structured circuit reference family can be rebuilt
analytically: its CMI and trace-distance values are retained to roundoff;
fidelity changes by about `2.1e-8` and `2.4e-8` at theta=.6 and 1, respectively,
consistent with the old inner-square-root evaluation's rank sensitivity.
This does not revalidate the unreplayable measured metrics.

The source-support follow-up changes live root, recovery and distance values
and their Stage 1 availability reporting. No additional live receipt producer
uses this API. Frozen summaries still lack the inputs needed for replay, so
the correction cannot refresh or promote their measured recovery metrics.

The separately frozen issue-509 benchmark, other receipts, registry payloads,
papers, book, simulator, physical model gates #1025/#1026 and byte-pinned
mandatory runner are unchanged. Existing papers already distinguish supported
Petz action from a full channel and recovery existence from a chosen map.
The full-domain completion here makes that existing distinction executable.

## Reproduction and hostile controls

```text
PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
```

The existing quantum-information workflow runs the first command on both
platforms; mandatory shard zero already runs the Stage 1 file on both.
At the initial reviewed head, the affected core suite passed 927 tests on
Linux and 921 on Windows, with six extended-precision skips on Windows;
the Stage 1 file passed 49 tests, for 95 added cases across both suites.
The positivity follow-up adds 38 independent Gram/Schur controls and five
tomography controls. Its core suite passes 965 tests on Linux and 959 on
Windows (six platform skips); all 54 Stage 1 tests pass on both platforms.
The same Stage 1 tests pass with the strict #1049 information validator.
Warnings are treated as errors.
The circuit follow-up adds thirteen cases (six pre-fix reproductions, five
invalid-amplitude controls, an exact product and a complete report replay).
All 67 Stage 1 cases pass on Windows and Linux in the combined four-PR tree;
that tree's core suite passed 1367 tests on Linux and 1361 on Windows with
six expected skips before this consumer-only correction. Core sources are
unchanged by the circuit follow-up. The PR adds 151 cases in total.
The tests include independent rational and symbolic controls, high-precision
matrix functions, complete Born-count reconstruction, invalid input rejection
under optimized Python, and a guard forbidding state repair during analysis.

Twelve isolated deliberate faults were detected: restored positive support
cutoff, omitted kernel completion, wrong inverse power, lost complex adjoint,
root fidelity reported as squared, constant perfect fidelity, zero-state
replacement, equal setting weights, old incomplete tomography, negative
counts, fabricated zero-shot data, and fidelity forced to meet the
optimal-recovery bound. Each ran in a separate copied package tree; collection
or import errors were not counted as detection. These controls exercise
incorrect scientific answers, not merely an implementation's success flag.

All twelve controls were rerun after the positivity repair. Seven additional
mutations are rejected: skip that repair, omit row compensation, round
diagonals to nearest, use the wrong complex Gram product, report a false zero
rounding bound, ignore a nonzero row at a zero Schur pivot, and return a
constant zero Gram matrix. The latter controls use exact returned-entry
inequalities and the full noisy-count path.

Four circuit-specific mutations are also rejected: restore the ordinary
rounded outer product, omit qubit-axis reversal, falsify the rounding bound
as zero, and drop the circuit certificate from the saved report. These
controls complement the actual SDK replay without making the mandatory
numerical tests depend on optional circuit or cloud packages.
