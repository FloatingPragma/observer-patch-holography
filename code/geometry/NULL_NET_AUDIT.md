# Audit of the finite Gaussian null-net report

This is a focused repair of **finite diagnostics promoted to stronger
null-net receipts**, under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
It concerns the supplied anti-periodic, half-filled free-fermion ring in
[null_net_receipts.py](null_net_receipts.py). It does not construct a physical
null net or change the axioms. The report now gives full-subspace leakage,
an actual bound for the particular expectation it measures, and complete
finite momentum-envelope sums. These replace selected-packet, finite-trend
and truncated-shape evidence that had been assigned stronger meanings.

## Reproductions before the repair

On main `7f5805a2`, the four regression tests committed before the repair
fail with these results:

* The full N=16 covariance has eigenvalues between approximately
  **-0.140729 and 1.140729**, so it cannot be a fermionic covariance. Its
  projector defect `||C^2-C||_F` is approximately 0.527622.
* `instrument_null_net((16,))` witnesses half-sided compression although
  the maximum leakage probability over the B subspace is **0.653606 in
  both time directions**. The selected packet's leakage ratio is 9.5669.
* That same single-stage report witnesses mixed-GNS Cauchy convergence and
  percent-level modular Lie closure, despite having no refinement difference
  and **no Lie samples**.
* `lie_closure_receipt(32, rmax=0)` returns a NaN residual after a warning.
  It does not reject the empty signal. NaN itself did not pass the old
  percent-level comparison; the empty-family `all(...)` did.

The old docstring also asserted that a nonzero imaginary antisymmetric
Hermitian commutator lies approximately in a span of real symmetric
generators. Section 4 gives an exact obstruction. Finite residual decrease
was incorrectly labeled a certified convergence rate.

## 1. A valid full covariance, and why finite arcs are faithful

Put `k_t=(2t+1) pi/N` and fill `cos(k_t)>0`, with N divisible by four.
The normalized discrete Fourier columns form a unitary matrix. If V
contains its occupied columns, then `C=V V*`, so `C^2=C`, `C*=C`, and
`rank C=N/2`. In the chosen gauge its entries are real:

\[
C_{ij}=\frac1N\sum_{\cos k_t>0}\cos(k_t(i-j)).
\]

The anti-periodic seam has `C(N-d)=-C(d)`. Replacing distance d by
`min(d,N-d)`, as the old helper did, therefore changes signs. The repair
uses the actual distance `|i-j|`. Independent complex Fourier projectors,
including the complementary unoccupied projector, check the whole matrix.
Local arc entries of size at most N/2 were unaffected by this seam bug;
the old interior bond expectation and arc modular generators do not change
for this reason.

There is also a structural proof of finite arc faithfulness, independent of
the small computed eigenvalue gaps. For an m-site consecutive arc with
`m<=N/2`, restrict V to those rows. Any m distinct occupied columns give a
Vandermonde matrix with distinct nodes `exp(i k_t)`, up to nonzero column
factors. It has rank m. Hence `C_arc=V_arc V_arc* > 0`. The same argument
using unoccupied columns gives `I-C_arc>0`. Every occupation eigenvalue
therefore lies strictly between zero and one. The reduced Gaussian density
has positive eigenvalues, products of occupation eigenvalues and their
complements, and is faithful.

The numerical gap calculation uses its own 120-digit context and rejects
an unresolved gap. The Vandermonde proof supplies strict positivity for
every mathematical finite stage, **not** a uniform lower bound or a
cyclicity/standardness result on a common limiting GNS space. The finite
relative commutant dimension `4^(|A|-|B|)` remains valid: these are full
matrix factors. In CAR notation, complement generators need the B-parity
twist to commute with B; bare odd operators anticommute. Site coverage
likewise remains a finite algebra-generation fact.

## 2. One selected observable has a proved limit

Keep exactly the original mathematical observable, with the seam bond
excluded:

\[
E_N=\frac1N\sum_{j=0}^{N-2}f(j/N)
\langle c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j\rangle,
\qquad f(x)=e^{-200(x-1/4)^2}.
\]

The finite geometric sum over occupied momenta gives
`C_N(1)=1/(N sin(pi/N))`. Set `S_N=(1/N) sum_{j=0}^{N-2} f(j/N)`,
`I=integral_0^1 f`, and `c_N=2/(N sin(pi/N))`. Thus `E_N=c_N S_N`
and the limit is `L=(2/pi) I`, approximately **0.07978843321**.

For completeness, `0<=f<=1` and its total variation on [0,1] is less
than 2. On each mesh cell, its left-rectangle quadrature error is bounded
by the cell width times its variation on that cell. Summing gives
`|(1/N) sum_{j=0}^{N-1} f(j/N)-I|<=2/N`. Removing the last term gives
`|S_N-I|<=3/N`.

For `x=pi/N<=pi/8`, `x-x^3/6<=sin x<=x`. Therefore

\[
|E_N-L|\leq\frac{2}{\pi}
 \frac{3/N+x^2/6}{1-x^2/6}
\leq B_N:=\frac23\frac{3/N+a_N}{1-a_N},
\qquad a_N=\frac{(22/(7N))^2}{6}.
\]

The second inequality uses `3<pi<22/7` and monotonicity in the nonnegative
argument a. `B_N` is an exact rational, tends to zero, and is less than
`3/N` for N>=8. The latter inequality is equivalent to
`a_N(2N/3+3)<1`; using `a_N<2/N^2` bounds its left side by
`4/(3N)+6/N^2`, already below one at N=8 and decreasing thereafter.
Consequently `|E_N-E_M|<=B_N+B_M<3/N+3/M`. This controls every
allowed pair, not only a sampled decreasing sequence. In particular,
stages above `6/epsilon` are Cauchy to epsilon for this observable.

The report stores exact rational bounds. Its expectation and limit floats
are numerical evaluations, **not** outward-rounded enclosures. Separate
Fourier sums and adaptive quadrature test the values and bounds. This
proves convergence of the single supplied observable. It does not prove
the mixed inner products, compatible embeddings and observable-family
convergence that the mixed-GNS receipt consumes.

## 3. Why packet asymmetry cannot certify half-sided inclusion

Let P project onto a finite nonempty proper subspace B, `Q=I-P`, and
`V(t)=exp(i h_A t)`. The internal factor `W(t)=exp(-i h_B t)` preserves
B, so

\[
\sup_{\psi\in B,\ \|\psi\|=1}\|QV(t)W(t)\psi\|^2
=\|QV(t)P\|^2.
\]

Write `A=P V(t) P` as a square matrix on B. Unitarity gives

\[
(QVP)^*(QVP)=I_B-A^*A,
\qquad
(PVQ)(PVQ)^*=I_B-AA^*.
\]

The matrices `A*A` and `AA*` have the same eigenvalues, including zero
multiplicities, since A is finite and square. The two leakage blocks have
the same singular values apart from harmless rectangular zeros. Since
`V(-t)=V(t)*`, the maximum leakage probabilities are identical for t and
-t. This holds also for complex Hermitian generators and unequal B and
complement dimensions; it is stronger than real-matrix time reversal.

Moreover, if a unitary sends finite B into B, its image has the same
dimension and equals B. Proper one-sided subspace inclusion is impossible.
The same elementary dimension argument excludes proper inclusion of a
finite-dimensional algebra into itself by conjugation. This does **not**
exclude half-sided modular inclusions of infinite-dimensional algebras or
of real standard subspaces in their proper setting. It explains why this
particular finite test cannot serve as their receipt.

An exact three-dimensional counterexample makes the packet problem sharp.
Let U cycle `e0 -> e1 -> e2 -> e0` and `B=span(e0,e1)`. The packet e0 has
zero forward leakage and unit backward leakage, while the worst-case
leakage is one in both directions. U has a Hermitian logarithm, so it is a
member of a continuous finite unitary flow. Exact symbolic matrix checks
and a numerical realization test this example independently.

On the actual supplied rings at model time 0.12:

| N | Selected packet leakage ratio | Maximum leakage probability, either sign |
|---|---:|---:|
| 16 | 9.5669 | 0.653606 |
| 32 | 38.2244 | 0.983561 |
| 64 | 2941.09 | 0.999994 |

The large packet ratio survives as a directional transport diagnostic.
It cannot discharge half-sided inclusion. Zero leakage ratios are recorded
as undefined where appropriate, rather than divided by an invented floor.
The flow's numerical unitarity defect is reported and must be <=1e-10;
that tolerance is a numerical resolution guard, never a proof of inclusion.

## 4. What the commutator envelope can and cannot establish

For real symmetric H and K, `D=i[H,K]` is imaginary antisymmetric and
Hermitian. Every real symmetric S is orthogonal to D in the real
Hilbert-Schmidt inner product. Thus

\[
\|D-S\|_F^2=\|D\|_F^2+\|S\|_F^2.
\]

For nonzero D, the best relative distance to the entire real symmetric
space is **one**. The old stated span of real arc modular generators,
number operators and identity is contained in this space and contains
zero, so its best distance is also one. An exact Pauli X/Z calculation
checks `||i[X,Z]-[[a,b],[b,c]]||_F^2=8+a^2+2b^2+c^2`.
This obstruction is to that particular claimed span; it does not rule
out a larger conformal algebra with the required momentum generators.

The code actually measures the real antisymmetric commutator C through
the scalar envelope

\[
p(j)=\sum_{r\geq2,\ r\ \mathrm{even}}
 \frac r2(-1)^{(r-2)/2}C_{j+1-r/2,\,j+1+r/2},
\]

retaining terms whose indices fit. This map is not injective even among
even-range commutators. For a 6x6 antisymmetric C with `C[1,3]=2` and
`C[0,4]=1` (and their negative transposes), p is identically zero although
`||C||_F^2=10`. It is itself a commutator: take `H=diag(0,...,5)` and
the symmetric `K[i,j]=C[i,j]/(i-j)` off diagonal. Multiplying this invisible
component by an arbitrary constant leaves every envelope sample fixed.
An accurate envelope fit therefore cannot bound full operator error.

The repair sums all fitting even ranges by default. For a supplied cutoff
R it also reports the triangle bound

\[
|p_{\rm full}(j)-p_R(j)|\leq
 \sum_{r>R,\ r\ \mathrm{even, fitting}}\frac r2|C_{j+1-r/2,\,j+1+r/2}|.
\]

This is an exact finite inequality for exact entries. Reported evaluation
of that bound is floating-point, not an interval certificate. Independent
exact entry-by-entry enumeration checks it, including R beyond the matrix.
The historical R=8 fit remains visible:

| N | Old R=8 shape residual | Complete finite sum residual |
|---|---:|---:|
| 32 | 0.00165175 | 0.000903682 |
| 64 | 0.00663635 | 0.000369533 |

Removing a fixed-range truncation improves this supplied-family diagnostic.
Neither two decreasing residuals nor fitted normalization prove a continuum
rate or an operator relation. Empty sample families cannot pass, fewer than
two fit coordinates or zero/nonfinite signals are rejected, and no finite
trend is promoted to a convergence-rate certificate.

## Numerical follow-up audit

An adversarial review of PR head `e732046f` found two additional failure
classes, reproduced by six failing regression cases before their repair:

* Dividing a complex packet by a subnormal scale, such as `1e-320`, can
  overflow NumPy's internal reciprocal and return NaN leakage. The same
  operation in the Hermitian check produced NaN comparisons, admitting
  `[[0,1e-320],[0,0]]` as a complex Hermitian matrix.
* Input promotion converted the unequal integers `2^53` and `2^53+1` into
  equal doubles, yielding a falsely perfect scalar fit. Squaring the
  residual of `measured=[1,1e-200]`, `target=[1,0]` also returned zero
  although its relative norm is representable near `1e-200`.

The repair scales the real and imaginary components separately. It supports
packets from the smallest binary64 subnormal through components of magnitude
`1.7e308`, including complex magnitudes that overflow if formed directly.
Unresolvable component dynamic ranges fail explicitly. The Hermitian guard
uses a `64*eps` tolerance after component scaling and the flow retains its
separate numerical unitarity check; neither tolerance proves an exact identity.

Input validation and rational-norm conversion reuse the audited primitives
in `null_tomography.py`. They inspect original sequence elements before
promotion, reject booleans, nonfinite values and lossy conversion, and retain
representable small results. For the one-scalar shape fit, compute the exact
rational least-squares coefficient of the supplied binary samples, round
that coefficient once, and replay the **returned** coefficient against the
original samples. Both squared residual and squared measured norm are
computed exactly. The diagnostic threshold compares their ratio directly
to `(1/50)^2`; it never compares an underflowed or rounded displayed norm.
The displayed norm is still a numerical square root. The exact ratio concerns
the given binary samples, not eigensolver, envelope-construction or continuum
error. An unrepresentable nonzero coefficient or norm is rejected.

Independent 400-digit projection controls cover twenty scaled data sets;
analytic Pauli evolution checks the subnormal complex flow. The rounded
subnormal-coefficient regression separately checks that representation error
is included in the returned residual. The default ring conclusions and
reported precision-level agreement remain unchanged after these repairs.

## Contract, custody and downstream impact

The public function names and report path remain for existing callers, but
the report is **schema 2**. Legacy `lie_closure_*` verdicts refer explicitly
to scalar envelopes; the corresponding stronger `receipts_witnessed` flag
is false. Mixed-GNS and HSM witness flags are also false. The finite
algebra facts and the new single-observable limit are positive results.
`compressing_sign` becomes `lower_packet_leakage_sign`; packet ratio and
direction verdicts are named as packet diagnostics. This is an intentional
meaning/schema correction. No in-repository consumer depends on the removed
keys beyond the updated tests.

Ring families must be nonempty, valid and strictly increasing. Full report
rings are divisible by eight; the standalone covariance and bond helpers
also accept multiples of four. Lie fits require N>=32 and divisibility by
eight. Invalid or unresolved numerical inputs fail; JSON rejects NaN and
infinity. Matrix, packet, time and fit inputs must be exactly representable
as finite real or complex binary64 components. Modular generators retain the supplied model and 120-digit
eigendecomposition, rounded to binary64. No physical clock is attached.

`claims/claim_registry.yaml` cites this report under
`OPH-GR-E2E-BRANCH-ENTRY`, already classified as a conditional composition
with open nonemptiness. Its null-net claim explicitly says finite
one-particle diagnostics do not instantiate the scaling-limit receipts.
The spacetime paper (including Theorem 4.3d in the shared technical fragment),
flagship and book use common-GNS
standardness, cofinal relative-commutant cyclicity, actual half-sided
inclusions and modular-intersection premises. This repair supplies none of
those stronger premises and changes none of their theorem statements.
Current public surfaces already retain the physical/nonemptiness boundary;
no paper/PDF, claim classification, observation rung, premise, Lean theorem
or axiom is promoted or rewritten.

The generated report at its existing registered path is deliberately
replaced because its witness flags were wrong. The prior report remains
recoverable from main `7f5805a2`; the corrected report retains the packet
numbers and historical R=8 values for comparison. No frozen registration,
source routing tape or hash-pinned evidence is changed. The registry points
to a mutable generated report, not a pinned hash. There is no M1 or
physical-model completion claim, and standing issue #1033 remains open.

## Reproduction and independent controls

```sh
OPENBLAS_NUM_THREADS=1 python code/geometry/null_net_receipts.py
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/geometry/test_null_net_receipts.py -W error
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
```

The existing Linux finite-geometry and quantum-information workflows
run these tests. They include complex Fourier projectors, independent
spectral evolution, exact permutation and Pauli controls, an invisible
commutator, exact tail sums, separate quadrature, adversarial empty/zero/
nonfinite inputs, and full schema/numerical replay of the generated report.
The proofs above supply the general statements; finite tests support their
implementation and do not substitute for those proofs.

Local validation after the follow-up audit: **710 Linux tests passed;
705 Windows tests passed with five expected extended-exponent platform
skips**, with warnings treated as errors. The null-net module has 93 tests.
Claim-registry, axiom-consistency and reader-style checks pass. The unchanged
Lean trust-inventory checker passes on Linux; on Windows it reports a path
separator mismatch (`Screen\\...` versus `Screen/...`), not a changed proof
inventory. No Lean source or inventory checker is edited here.
