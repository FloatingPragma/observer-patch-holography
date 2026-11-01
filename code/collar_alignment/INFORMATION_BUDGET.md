# Resolve the full collar information budget

Maintainer follow-up (#1049): support containment in the relative-modular
moment test is now compared in the canonical Gaussian rational field, using
a zero residual. Symbolic expression-tree equality falsely rejected a
full-rank complex conditional-product state and a nearby non-Markov state.
The new regression exercises exact Markov, public CMI, collar CMI and the
alignment budget, with both full-rank and singular conditional products;
real and complex support escapes must still raise. This changes acceptance
of valid inputs, not the moment theorem, frozen evidence or paper claims.

The [review summary](https://github.com/FloatingPragma/observer-patch-holography/pull/1049#pullrequestreview-5423084144)
and [inline comment](https://github.com/FloatingPragma/observer-patch-holography/pull/1049#discussion_r4190904723)
identify this same defect. The eleven direct controls in
`test_complex_markov_support.py` cover both entries' requirements; five fail
against the reviewed `10fffe19` implementation. The broader follow-up retains
24 independently constructed complex conditional-product families, with
unequal sector weights and complex states at both endpoints. In each case,
an explicit within-sector correlation preserves both marginals and destroys
conditional independence. Exact Markov, public CMI and collar CMI must accept
the product and distinguish the perturbed state. Twenty-three of these 24
controls fail on the reviewed implementation; all 35 controls pass with
canonical support equality. The branch includes main `7f99c8cc`.

This is a focused evidence repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033),
built on main `3d848b7f`. It repairs false zero information scores and makes
the existing alignment criterion more informative. It builds on the
sector-label and Gibbs audits without depending on the first-law repair
#1048 (merged). No physical model or new physical hypothesis is selected.

## Reproduced failures

The first commit, `568cc679`, contains fourteen tests that fail before the
implementation changes. The previous shared routines computed information
by subtracting binary64 entropies. With ordinary trace and natural logarithms:

| Supplied state | Previous result | Independent result |
| --- | --- | --- |
| `diag(1/4+e,1/4-e,1/4-e,1/4+e)`, `e=2^-30` | MI = 0; collar alignment accepted at tolerance `e^2` | MI = `6.938893903907228e-18`, greater than the tolerance |
| Classical three-bit parity perturbation `p=1/8+e*(-1)^(a+b+c)`, `e=2^-30` | CMI = 0 | CMI = `2.7755575615628914e-17` |
| `I_8/8+e X tensor Y tensor X`, `e=1e-150` | CMI = 0 | CMI approximately `3.2e-299` |
| `I_4/4+1e-200 X tensor X` | MI = 0 | Positive, below binary64 range; the numerical API must raise |
| A masked density entry | Hidden payload accepted | Missing evidence must be rejected |

The scalar controls use the separate formula
`[(1+x)log(1+x)+(1-x)log(1-x)]/2`, where `x=d*e`.
Further controls use independent 2-by-2 characteristic roots and matrix
logarithms. They also retain the fourth-order example
`rho=I_4/4+e*(X tensor I+I tensor X)`, whose MI is
`128 e^4+O(e^6)`. Its second-order entropy terms cancel; resolving only the
leading state perturbation does not resolve its information.

## Four terms, with no additional state-selection premise

Fix the already-declared collar tensor factors `A,L,R,D`. For every
normalized finite state, including singular states, the chain rule gives

```text
I(AL:RD) = I(L:R) + I(A:R|L) + I(L:D|R) + I(A:D|LR).             (1)
```

To prove this, expand `I(AL:RD)` first across `R,D`, then expand each term
across `L,A`:

```text
I(AL:R)   = I(L:R)   + I(A:R|L),
I(AL:D|R) = I(L:D|R) + I(A:D|LR).
```

Every term is nonnegative by subadditivity or strong subadditivity. These
are established finite quantum-information results; see
[Watrous, chapter 5](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
Equation (1) therefore gives an exact equivalence: the specified cut is a
product precisely when all four terms vanish. It also gives a quantitative
budget: the whole defect is at most the sum of any four valid upper bounds.
No faithfulness, Gibbs model, small-coupling expansion or continuum limit
is required for this identity.

The usual collar CMI is only the last term. Each term can be the sole
obstruction: correlate just the corresponding pair of classical bits and
leave the other two bits independent. The retained Bell counterexample
has zero collar and middle terms but two endpoint terms of `2 log(2)` each,
giving total alignment defect `4 log(2)`. Thus its obstruction is located,
rather than merely labelled nonaligned.

`alignment_information_budget` reports all four terms for every positive
sector, their weighted values, the existing maximum-sector defect, and
the numerical chain-rule residual. A small weight can make the average
small while the maximum remains large. Zero-weight sectors have no effect.

## The whole budget is a variational distance

Let the retained center and cuts be fixed, with
`rho=direct_sum_j p_j rho_j` and normalized conditional states. Compare with
an arbitrary aligned state on this same declared algebra,
`tau=direct_sum_j q_j (sigma_L,j tensor sigma_R,j)`. Block logarithms and
partial traces give the exact decomposition

```text
D(rho || tau)
 = D(p || q)
   + sum_j p_j [ I(L_j:R_j)_rho_j
                 + D(rho_L,j || sigma_L,j)
                 + D(rho_R,j || sigma_R,j) ].                 (2)
```

Support escape has its usual infinite-divergence meaning. All added terms
are nonnegative. Consequently the minimum is
`J=sum_j p_j I(L_j:R_j)`, attained by retaining `p_j` and replacing each
conditional state by its own marginal product. Equality fixes these
weights and marginals on positive sectors; choices on absent sectors are
irrelevant. This proves the minimizer directly, without an optimizer.
For the collar, `L_j` means `A bL,j` and `R_j` means `bR,j D`, so (1)
resolves this same minimum into four contributions.

Quantum Pinsker gives
`||rho-rho_aligned||_1/2 <= sqrt(J/2)`. The implementation evaluates this
inequality numerically, not as an interval certificate. For total mass `T`
the homogeneous version is `sqrt(T*J/2)`, capped at `T`; the report retains
that mass and evaluates the square roots before multiplication to avoid
intermediate underflow.

This projection is onto the fixed center and cut. It does not optimize
over all possible quantum Markov decompositions. In particular, a quantum
CMI alone is not generally the relative-entropy distance to that larger
family; see [Ibinson, Linden and Winter](https://arxiv.org/abs/quant-ph/0611057).
The repair does not use that false identification.

## Numerical implementation and exact zero decisions

The information evaluator reuses the shared input conversion and dimension
checks. It rejects masked, Boolean, nonfinite and lossy-conversion inputs.
Accepted Hermitian roundoff is averaged with exact rational arithmetic.
Partial traces then sum the supplied entries exactly, including rare-event
mass that would disappear when added to a unit diagonal entry.

Positivity is checked algebraically on that Hermitian representative.
Its characteristic polynomial supplies the elementary symmetric functions
of its real eigenvalues. They are all nonnegative precisely when the
matrix is PSD: otherwise `det(t I+A)` has a positive root, whereas a
polynomial with nonnegative coefficients and positive leading coefficient
cannot vanish for positive `t`. Trailing zero coefficients determine exact
nullity. A negative eigenvalue admitted by the older roundoff tolerance
is therefore rejected instead of clipped into apparently valid evidence.

For a classical joint distribution, define
`q_abc=p_ab*p_bc/p_b`, omitting zero-mass conditioning blocks. This is a
nonnegative distribution with the same total mass as `p`. Its divergence
is exactly the CMI. The evaluator sums
`p log(p/q)-p+q`, using exact probability differences and a convergent
`log1p` remainder near equality. These terms cannot cancel a positive
correlation against an order-one entropy.

General quantum combinations use private 80-digit arithmetic, escalating
to 400 digits when needed, and must agree with a separate evaluation at
40 additional digits. This is a precision diagnostic, not an interval
bound on a numerical eigensolve. Exact nullity distinguishes zero
eigenvalues from unresolved positive spectral mass. Nonzero output that
underflows or suffers excessive subnormal quantization raises explicitly.
The global mpmath context is unchanged.

No uncertain entropy cancellation is floored to zero. Exact product checks
handle common zero cases without requiring resolved numerical eigenvalues:
positivity has already been proved exactly. Dimension-one conditioning
factors do not add distinct product splits. Only the nontrivial factors
are enumerated, so the number of splits is at most the conditioning-space
dimension, regardless of how many trivial factors are declared.
An exact isospectral check also recognizes hidden
single-product decompositions: compare the power traces of
`rho_AB tensor rho_BC` and `rho_B tensor rho_ABC` through their common
dimension. Newton identities determine their spectra, and the tensor
entropy identity then gives zero CMI. This shortcut is sufficient but is
not used as a complete Markov criterion.

The final exact test closes that recognition gap. For source `r` and
reference `s` containing its support, set

```text
M_k(r,s) = Tr(r^k s^(1-k)),             k=1,2,...,
```

where negative powers act on the exact support of `s`. Compare the pairs

```text
(rho_ABC, rho_A tensor rho_BC),
(rho_AB,  rho_A tensor rho_B).
```

**Finite zero test.** CMI is zero if and only if these moments agree for
`k=1,...,N`, where `N=n_ABC^2+n_AB^2`.

For sufficiency, the positive relative modular operator
`Delta=L_r R_(s inverse)` has at most `n^2` positive spectral points.
With vector `sqrt(s)`, its positive spectral measure has moments `M_k`
and its `t log(t)` integral is `D(r||s)`. Subtract the two measures and
combine coincident positive spectral points. There are at most `N` points.
The first `N` positive moments determine every remaining weight by an
invertible Vandermonde system; the points are nonzero, so starting at
power one is sufficient. Zero spectral points do not contribute to
`t log(t)`. Thus the divergences agree. Their difference is precisely CMI.

For necessity, use the established equality structure
`rho_ABC=direct_sum_j p_j rho_A,Lj tensor rho_Rj,C`, with the label retained
in B ([Hayden, Jozsa, Petz and Winter](https://arxiv.org/abs/quant-ph/0304007)).
In each sector both moments reduce to
`p_j Tr[rho_A,Lj^k (rho_A tensor rho_Lj)^(1-k)]`;
the right factor contributes trace one. This also holds on singular
supports. Summing proves equality for every `k`. A common nonunit source
mass multiplies both moment sequences by the same factor at each power.

All these matrix products, inverses on support, traces and comparisons are
rational operations for the supplied binary64 entries. `is_markov_exact`
exposes this decision independently of any floating information value;
it can reject a non-Markov perturbation even below binary64 output range.
The exact moment path is intended for small finite witnesses and can be
more expensive than numerical evaluation. It proves a property of the
specified matrix, not of an unknown noisy preparation.

### Additional audit findings

Commit `4fe2afd1` reproduces four failures in the first version of this PR.
An exactly independent state and a conditional Markov state were rejected
because their positive eigenvalues could not be resolved at 400 digits.
Both use the six-dimensional integer Gram matrix
`G=(I+2^26 S)^T (I+2^26 S)`, where `S` is the unit superdiagonal, then
multiply it by the smallest positive binary64 value.
Every input entry is representable, although a positive eigenvalue is below
`1e-400`. Requiring a numerical entropy before honoring an algebraic zero
was unnecessary. Exact product zeros now return after exact PSD validation;
other unresolved spectral cases reach the exact zero test before refusal.
A nearby state with an additional classical conditional correlation still
fails explicitly, so this fallback cannot simply declare unresolved input
Markov.

The other two failures concern the numerical and exact APIs with 24
dimension-one conditioning factors. They previously attempted `2^24`
copies of the same product check on a four-dimensional state. The shared
split helper now tests it once. Regression controls also retain nontrivial
conditioning partitions and exercise complex singular quantum supports
through the full modular-moment test, bypassing its shortcuts.

## Compatibility and downstream impact

- Existing MI/CMI function signatures and nats remain unchanged. CMI keeps
  its four-entropy definition. MI now uses the same formula with an empty
  conditioning system, whose entropy is `-T log T`. Consequently MI is
  `T I(rho/T)` for accepted trace-roundoff inputs; this removes an artificial
  trace-offset contribution. No matrix is silently renormalized.
- The collar wrapper preserves entries until exact Hermitian averaging.
  Its weighted CMI is accumulated as an exact sum of products before final
  conversion. It cannot turn an unrepresentable nonzero total into zero.
- The live null-net standardness witness inherits the corrected CMI. Its
  Markov comparison remains below its stated tolerance, and its interacting
  Gibbs comparison remains non-Markov. The separate Gaussian null-net
  receipt producer does not call this information evaluator.
- MaxEnt and Einstein receipts continue to use their existing entropy and
  relative-entropy interfaces. No frozen registration, pinned receipt,
  claim payload, paper, book or release artifact changes. Exact normalized
  formulas and physical qualification in those surfaces are unchanged.
  The public distinctions between a declared finite collar and a physical
  source remain in place; this work does not select #1025's model or satisfy
  #1026's start gate.
- The existing finite-quantum-information workflow collects the new tests
  on Linux. The byte-frozen mandatory runner is unchanged.

## Reproduction and independent checks

```text
PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
```

Controls include independent scalar spectra through fourth order, separate
matrix logarithms, noncommuting nonuniform Markov perturbations, exact
tensor permutations, rare-event probabilities, zero and singular supports,
missing data, negative states, underflow and weighted-sector distinctions.
The exact modular-moment verifier is checked against **all 255 nonempty
classical three-bit supports**, using integer conditional-independence
identities as the independent answer. Quantum controls include a basis
change hiding a product collar and a multi-sector state with different
left/right decompositions inside the conditioning space.

The variational test evaluates an independently chosen trial family with
wrong marginal states and wrong sector weights, retaining every extra term
in (2). Four separate witnesses isolate each term of (1). The proof concerns
all finite states in its stated algebra; these controls check the executable
diagnostic and its input boundary rather than replacing the proof.

Initial local validation before maintainer review: **859 tests passed on
Linux; 853 passed on Windows with six existing extended-precision skips**,
with warnings treated as errors. The initial implementation added **77
cases**, including fourteen original pre-fix
failures and four audit regressions committed before their repair.
The null-net Gibbs witness has CMI `0.006842238011420933`; its floating
Markov construction has residual CMI `1.9828816605061793e-31`. The latter
is below the existing numerical threshold, not an exact-zero certificate
for the rounded constructor output.

After the maintainer correction, the additional 35 support controls and
the merge of main `7f99c8cc`, the same complete local suite passes **989
tests on Linux and 983 on Windows**, with six expected Windows
extended-precision skips and warnings treated as errors. The PR now adds
112 cases in total. Replaying the reviewed implementation fails five of
the eleven direct controls and 23 of the 24 broader complex-state controls;
the canonical-arithmetic correction passes all 35. Claim registry, axiom
inventory, reader style and whitespace checks also pass.

Sixteen isolated mutations were tested against 75 of the new controls
(the full 255-support enumeration and slower multi-sector replay were
separately included in the complete platform suites). Every copy contained
its own source and tests outside the repository import path:

| Deliberate defect | Failing controls |
| --- | ---: |
| Erase quantum coherences | 30 |
| Round conditional reference probabilities first | 1 |
| Omit total-mass correction | 2 |
| Accept nonzero output underflow | 1 |
| Skip exact positivity | 1 |
| Check only the first modular moment | 7 |
| Reduce quantum precision | 9 |
| Suppress one endpoint obstruction | 6 |
| Discard later weighted sectors | 7 |
| Trust a product without comparing it | 29 |
| Return bits through a nats interface | 21 |
| Transpose the relative-modular reference | 2 |
| Require resolved eigenvalues for an exact product | 1 |
| Repeat splits over dimension-one factors | 2 |
| Reject unresolved spectra before exact Markov verification | 1 |
| Declare unresolved spectra Markov without verification | 1 |

The transposed-reference mutation initially survived the fast subset.
The retained complex BC-state control closes that coverage gap; its B
marginal is real, so the wrong reference cannot cancel between both sides.
The initial fast-test selector also unintentionally excluded four
`small_classical` cases; the final replay uses complete function names to
exclude only the two expensive tests. All cases run in the platform suites.
The numerical output is still a checked approximation. The PSD, product,
support and finite-moment decisions use exact rational arithmetic.
