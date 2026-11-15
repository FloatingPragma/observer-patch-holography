# Finite-state entropy and the role of sector labels

The [Gibbs response contract](../maxent/GIBBS_RESPONSE_AUDIT.md) distinguishes
resolved binary64 states from resolved weighted responses. State producers
conservatively reject subnormal thermal populations, including joint sector
weights. The public MaxEnt covariance retains populations and original
observable entries in multiprecision until the final response is rounded.

This shared implementation supplies density-matrix, entropy,
partial-trace and modular-splitting operations to the MaxEnt, collar-alignment,
null-net and Einstein-closure evidence. It uses natural logarithms and the
ordinary matrix trace. It implements finite numerical diagnostics, not a
physical state selection or a proof of exact rank from approximate matrices.

The [entropy-response identities](../geometry/ENTROPY_FIRST_LAW_AUDIT.md) separate
tangent first laws from finite entropy changes with relative-entropy
remainders. `entropy_response` tests independently supplied generators on
all normalized tangents, supplies central probability-transfer witnesses,
and computes Gibbs tangent slopes separately from finite secants.

The [algebra and repair results](ALGEBRAS_AND_REPAIR.md) give regional
separation and basis normalization through reduced spectra, and classify
exactly when finite relative-entropy completion is an affine quantum repair.
They include independent recognition
of supplied conditional expectations and a noncommuting repair convergence
proof. The [operator-span contract](ALGEBRAS_AND_REPAIR.md#preserve-the-supplied-span-before-computing-its-action)
distinguishes exact annihilators from ill-conditioned basis artifacts and retains exact
source intersections before numerical conversion. Independent original-input
controls also cover valid complex algebras and extreme basis units.

The [recovery evidence](RECOVERY_EVIDENCE.md#source-support-square-roots-and-reported-distances)
preserves original-entry support in square roots and state fidelity, and checks
Petz action beyond trace preservation. These APIs require exact positivity of
the supplied Hermitian average; unresolved numerical functions raise instead
of clipping a spectrum or manufacturing a recovery score. Stage 1 preserves
independent CMI and acquired counts when recovery is numerically unavailable.

The [fixed-point and relaxation results](REPAIR_FIXED_POINTS.md) prevent almost coincident
repair algebras from manufacturing protected observables and an artificially
large relaxation gap. Intersections use exact retained input spans; channel
and gap calculations are numerical diagnostics. The note derives the
complete decay law of the counterexample and separates permanent records
from records whose drift stays small over a finite observation time.

`expectations.repair_generator` accepts rates in any common clock unit.
Its identity and spectral-resolution checks use the generator divided by
`max(rates)`; the returned generator, gap and defect magnitudes retain the
input rate units. The report's `rate_scale` is this divisor, and `tolerance`
applies to defect magnitudes divided by it. Unresolved relative gaps, lost
positive rates, and unrepresentable rescaled generators still raise; a small
absolute rate alone does not imply poor numerical resolution.

<a id="audit-findings-and-corrections"></a>

## Input contracts and independent controls

| Input or mathematical distinction | Required implementation behavior |
| --- | --- |
| A multi-sector Markov state requires its central label in the conditioning system. | Construct a labelled direct sum. Its CMI is zero to numerical precision; erasing the label need not preserve that conclusion. |
| Empty or zero-weight block families, NaNs, negative states and inconsistent dimensions do not establish alignment. | Require a nonempty normalized block family and finite normalized positive states before evaluating the criterion. |
| An all-observable modular-invariance criterion cannot be certified by random probes. | Test every matrix unit of the left algebra; linearity supplies completeness. |
| Orthogonal pure states have infinite relative entropy. A negative density is invalid. | Validate the states and preserve reference zeros. Detected support escape returns infinity; full-algebra modular logarithms reject singular states. |
| Gibbs constraints, multipliers and central-sector energies define one complete optimization problem. | Validate complete pairings, independence modulo identity, faithful targets and convergence. Use one common energy shift across sectors and reject spectral underflow instead of flooring it. |

Positive eigenvalues are not discarded below `1e-14` or replaced by a
`1e-300` floor. The shared `1e-12` tolerance only admits floating-point error
in Hermiticity, normalization and negative eigenvalues. Accepted Hermiticity
roundoff is symmetrized and negative spectral roundoff is set to zero for
spectral calculations. Inputs are never implicitly normalized. No positive
support leakage is rounded away to make relative entropy finite. Relative
entropy reevaluates ambiguous reference support from the original Hermitian
representative using exact positivity, rank and kernel checks. A rounded
zero eigenvalue or matrix/kernel product alone cannot justify infinity.
An unresolved finite evaluation raises. A diagonal reference supplies its
spectrum directly; explicit positive entries as small as `1e-310` remain
supported. The exact checks certify these supplied entries, not the ranks
of arbitrary approximate or physically uncertain matrices.

Additional controls separate spectral support from numerical resolution:

- A real Hadamard change of basis applied to spectra `(.5,.25,.25,0)` and
  `(.25,.5,.25,0)` preserves common support and gives `D=.25*log(2)`.
  Kernel roundoff cannot justify infinity. The routine returns the
  correct value when support is resolved or raises an explicit support error.
  A companion state with mass in the missing direction must give infinity
  or the same explicit uncertainty error, never a finite divergence.
- A smallest positive binary64 sector weight multiplying `I_2/2` underflows.
  Direct sums reject component underflow; explicit zero sectors and
  representable tiny sectors remain valid.
- Complex central energies or inverse temperatures cannot supply real sector
  probabilities. The scalar validators reject complex, Boolean, nonfinite and
  nonscalar inputs, including convergence and alignment tolerances.
- For `H=30 X`, both Gibbs eigenweights are positive, but reconstruction of
  the dense matrix erases the smaller direction. Gibbs states and sector
  blocks validate the reconstructed faithful state. Diagonal controls
  preserve the same small probability when representable. The public
  Duhamel response instead evaluates the original thermal family before
  output rounding; its resolved covariance can be returned without a
  representable binary64 state. See the [response contract](../maxent/GIBBS_RESPONSE_AUDIT.md).

## Preserve the source before evaluating entropy responses

Surviving positive populations and a correct derivative ratio do not suffice
to validate a finite entropy response. The source entries, their small
differences and the final contractions each require their own checks.

For `T = diag(0, 2^-60, 1)` and multiplier 700, subtracting the binary64
mean from T can erase the low-energy doublet splitting. The resulting
variance is `4.929838271879886e-305`, whereas the original scalar partition
function gives `1.88079096131566e-37`. Both derivative components change by
the same factor, so `dS/dt = 700` does not detect the error. All three thermal
populations are normal binary64 values; no population underflow is required.

The entropy-response API shares the original-input spectral machinery with
the MaxEnt covariance calculation. Supplied entries, energy differences,
populations and response contractions remain in multiprecision
arithmetic until output conversion. The returned state also retains small
real or imaginary coherences from the source eigenprojectors. For this API,
subnormal populations and unresolved faithful output states are refused;
the covariance-only API has its separately documented, broader domain.
Repeated-precision agreement is a numerical control, not an interval proof.

Relative entropy has a separate range-sensitive control. Set `p=2^-40`, `t=2^-538`,
`rho=[[p,t],[t,1-p]]` and `sigma=diag(p,1-p)`. These supplied matrices have
exact trace one and positive determinants. Their divergence rounds to
`3.5e-323`, but squaring a spectral overlap before multiplying its logarithmic
factor can produce zero. A correct support decision does not prevent this
arithmetic error.

Range-sensitive contractions trigger reevaluation from the original
Hermitian entries, with exact positivity, rank and support checks, increased
numerical precision and a second-precision comparison. An unrepresentable
nonzero result raises instead of becoming equality. Exact equality returns
zero; detected support escape retains infinity. Accepted trace offsets are
not normalized away. The clipped-negative spectral convention is retained
on the ordinary faithful-reference path; exact positivity is required when
source-level refinement or support checks are needed. These checks do not certify
arbitrary floating-point eigensolves. Cancellation against an accepted trace
offset also triggers reevaluation: `rho=(.5-t^2) I+t X`, `t=2^-27`, against
`I/2` has positive divergence `2.2597578014143572e-32`; a negative numerical
value cannot establish equality. Exact unequal-state cancellations
that the bounded numerical refinement cannot resolve raise explicitly.

Support is checked before classifying a divergence as infinite. For
`sigma=[[p,x+iy],[x-iy,1-p]]`, with `p=2^-20`, `x=2^-14` and binary64
`y=0.0009746528185015864`, the exact determinant is positive, although both
the numerical smallest eigenvalue and matrix/kernel product can round to
zero. The original determinant gives `D(I/2||sigma)=25.904130790941124`.
Exact kernel containment also retains genuine complex singular examples.

The finite-balance consumer retains original Hermitian differences for its
trace and modular pairings as well as its divergence. Otherwise accepted
one-sided subnormal entries can disappear during symmetrization before the
divergence calculation is called. A reference whose averaging loses component
precision is refused before taking its logarithm; ordinary Hermitian
roundoff remains admissible. When cancellation makes `D-Tr(rho-sigma)`
unreliable, that remainder is evaluated before either term is rounded. For
`rho=(.5+3*2^-53) I` against `I/2`, subtracting correctly rounded terms had
understates the remainder by one ninth. The logarithm is a numerical
calculation.

Near a subnormal rounding midpoint, two precision levels can agree while
both miss the term that decides the correct neighbor. Gibbs state components
whose rounding side is unresolved are refused. Neighboring representable
subnormals and algebraically certified zeros remain accepted. This resolution
gate is a numerical safeguard, not a general correct-rounding certificate.
The relative-entropy fallback uses the same midpoint gate and can increase
precision to resolve the side. Rational conversion of the final MP value
avoids a second rounding during conversion to a binary64 subnormal.
Two high-precision eigensolves increase runtime materially. This API serves
small finite diagnostics, not a cost-equivalent replacement for large dense solves.

The retained controls independently evaluate scalar partition derivatives
and qubit characteristic roots. They check variance and both derivatives,
complex state coherences, opposite signs, changed energy units, exact zeros,
support boundaries and lost-output refusals. Public MaxEnt and serialized
finite-entropy-balance consumers exercise the same defects. Existing fixed
receipt producers are subject to their numerical tolerances. These numerical
checks do not promote an experimental or physical claim.

## One normal form for the finite A3 objective

The existing [A3 specification](../../docs/AXIOM_REFERENCE.md) supplies fixed
positive weights `w_P`, local states `rho_P`, references `tau_P`, a compatible
feasible family K, and declared local traces. Set `W=sum_P w_P`, `p_P=w_P/W`,
and form the labelled direct sums

```text
R(rho) = direct_sum_P p_P rho_P,
T      = direct_sum_P p_P tau_P.
```

Use each declared local trace on its summand. On a positive block,
`log(p_P rho_P)=log(p_P) I+log(rho_P)` on its support. Thus

```text
sum_P w_P D(rho_P || tau_P) = W D(R(rho) || T).
```

Both sides are infinite in exactly the same support-escape cases. This is
an identity for the existing finite A3 data; it adds no new physical premise.
For a nonempty cover W is positive. If an empty cover is permitted, its
injectivity already forces K to be a singleton and the objective is zero;
there is then no remaining optimization to encode.
R is affine, so optimizing on K is equivalent to optimizing on its image
R(K). The A3 cover's injectivity makes the correspondence one-to-one.
The weights, references and compatible feasible image must still be retained;
optimization over the entire direct-sum state space is a different problem.

The label is the cover index, not a constructed joint state on the physical
union of patches. Three pair distributions that make AB, BC and CA each
perfectly anticorrelated have compatible uniform singleton marginals but no
joint classical ABC extension. Their labelled direct sum is nevertheless a
valid encoding of the A3 objective. The test suite checks this obstruction
by complete enumeration of the eight assignments.

This normal form also explains the bulk/edge accounting. For normalized
sector probabilities p,

```text
S(direct_sum_j p_j rho_j) = H(p) + sum_j p_j S(rho_j).
D(direct_sum_j p_j rho_j || direct_sum_j q_j tau_j)
  = D(p || q) + sum_j p_j D(rho_j || tau_j).
```

The second formula has the usual support convention: a positive p block
with zero q has infinite divergence. When p=q, its classical term vanishes.
When p varies, omitting that term changes the selection problem. For
`rho_j=rho_bulk,j tensor I_edge,j/d_j`, the first formula adds
`sum_j p_j log d_j`, with H(p) counted once in the bulk entropy. This is the
identity already used by the Einstein-closure code. An auxiliary cover label
does not acquire physical edge entropy merely by being introduced in software.

These identities follow by evaluating block logarithms and traces; they use
the finite entropy definitions in [Watrous, section 5.2](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
They reorganize existing OPH mathematics without selecting its physical inputs.

## Why a central label cannot be replaced by a mixture

The [quantum Markov structure theorem](https://arxiv.org/abs/quant-ph/0304007)
uses a decomposition of the conditioning system B:

```text
H_B = direct_sum_j (H_bL,j tensor H_bR,j),
rho_ABC = direct_sum_j p_j (rho_A,bL,j tensor rho_bR,C,j).
```

The sector j is retained in B. Every entropy in
`I(A:C|B)=S(AB)+S(BC)-S(B)-S(ABC)` contains the same H(p); these terms cancel.
The within-sector product terms cancel too, giving zero. Replacing the direct
sum by an unlabelled mixture does not preserve this argument.

The smallest retained control is an equal classical mixture of `A=J=C=0`
and `A=J=C=1`. With J in the collar, `I(A:C|J)=0`. After J is erased,
`I(A:C)=log(2)`. No approximation, simulation fit or optimizer decides this
difference. In multiple sectors, modular splitting must likewise be tested
inside each central block; zero CMI does not imply an unconditioned product
across one globally chosen left/right tensor cut. The numerical cross-lane
test verifies both assertions on the same labelled state.

The API consequence is explicit: `markov_chain_state(da,dl,dr,dc,weights)`
returns tensor dimensions `(da,len(weights)*dl,dr,dc)`, ordered
`A,(sector,bL),bR,C`. The single-sector interface is unchanged. An application
without access to that sector label must use the erased state and test its
actual CMI; it cannot assume the label exists as a physical repair.

## Downstream scope and reproduction

| Consumer | Current contract |
| --- | --- |
| Null-net standardness | Retains the central label in multi-sector states; single-sector modular-locality and nonlocal-Gibbs controls have their declared domains. |
| Collar alignment | Uses a complete basis and rejects vacuous or invalid families. The Bell counterexample is Markov but not aligned. |
| MaxEnt closure | Enforces input and convergence premises for the multiplier counts, nonclosure example and exact product-subfamily control. |
| Einstein closure | Uses the same central-label entropy identity, with separate first-law and normalization controls. |

The shared helpers are one implementation, not four independent calculations.
The [MaxEnt coordinate contract](../maxent/PROJECTION_COORDINATES.md)
accounts for observable units, basis and scalar energy origins, and
connects its invariant stopping error to a bound on the true optimum. It
retains the nonzero coarse-graining closure defect.

The [sector Gibbs construction](../collar_alignment/GIBBS_SECTOR_AUDIT.md)
applies the energy-origin contract to collars. It shares observable decomposition
and thermal diagonalization with MaxEnt, retains relative sector partition
functions, and prevents a large scalar origin from manufacturing alignment.
Its finite interaction/log-state identity separates
center-probability error from conditional-state error in relative entropy.
Independent controls use explicit index contraction, scalar closed forms,
SciPy matrix logarithms, complete classical enumeration and known product/Bell
states. Tests include malformed matrices, missing sectors, invalid tensor
partitions, singular references, tiny positive eigenvalues, actual label
erasure and execution under `python -O`.

From the repository root with the pinned requirements and `PYTHONPATH=code`:

```text
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent
```

These finite constructions do not amend A1–A3, prove a physical Markov collar,
identify gravity or select the physical observer model. This shared module
governs the live consumers named above.

The [collar information-budget identities](../collar_alignment/INFORMATION_BUDGET.md)
exclude entropy-subtraction false zeros in mutual and conditional mutual
information. Exact reductions, stable classical remainders and guarded
quantum evaluation preserve small correlations. `is_markov_exact` decides
zero CMI for the supplied finite matrix using exact relative-modular moments
when simpler certificates do not suffice. The same construction resolves fixed-cut
alignment into four nonnegative terms and proves that their weighted total
is the relative-entropy projection cost onto the declared aligned family.
For trace-roundoff inputs, MI uses the homogeneous convention `T I(rho/T)`;
unresolved precision, indefinite states and nonzero output underflow raise.
