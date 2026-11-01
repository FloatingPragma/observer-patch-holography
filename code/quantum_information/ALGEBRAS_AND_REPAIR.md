# When finite entropy completion is a repair channel

For a finite observable algebra and a faithful reference state, relative-entropy
completion of all its marginal states is affine **exactly when** the algebra
is invariant under the reference's modular flow. In that case the completion
is the dual of the unique reference-preserving conditional expectation.
Outside that case, no single linear channel implements all those completions.
This classifies the whole finite algebraic problem, including boundary states;
it is not a fit to selected successful inputs.

The result connects the finite A3 objective, the collar alignment criterion
and the conditional expectations in the
[current source-selection contract](../../model/SOURCE_SELECTION.md).
It also identifies which premises are doing work. An observable span need
not be an algebra, a weighted projection need not be positive, and a selected
reference need not be compatible with a proposed retained algebra.
These are three different tests.

The implementation uses the shared [entropy and sector-label conventions](README.md),
ordinary matrix trace and natural logarithms. The finite mathematical results are proved below;
the Python checks are floating-point diagnostics supplemented by exact
rational counterexamples. None of these checks identifies a simulator's
transition law merely by constructing the desired channel.

## Correct subsystem and separation calculations

Let \(H=H_R\otimes H_S\), with dimensions \(D=d_Rd_S\), where the tensor
slots in R can be noncontiguous. On the full Hilbert-Schmidt space \(HS(H)\),
use left multiplication \(L(a)v=av\) and the vector \(\Omega=\sqrt{\rho}\).
Then

\[
\|L(a_R\otimes I_S)\Omega\|^2
 =\operatorname{Tr}(\rho_R a_R^*a_R),\qquad
\|L(a_R\otimes I_S)\|_{\mathrm{HS}}^2
 =D d_S\operatorname{Tr}(a_R^*a_R).
\]

Consequently the separation modulus with the *ambient Hilbert-Schmidt
operator norm* used by the standardness helper is

\[
\inf_{\|L(a)\|_{\mathrm{HS}}=1}\|L(a)\Omega\|
 =\sqrt{\frac{\lambda_{\min}(\rho_R)}{D d_S}}.
\]

The infimum follows by choosing a rank-one \(a_R^*a_R\) in a minimum
eigenspace. In particular, the represented regional algebra separates
\(\Omega\) exactly when its actual reduced state is faithful. For the full
algebra, both cyclicity on full \(HS(H)\) and separation hold exactly when
\(\rho\) is faithful: right multiplication by \(\sqrt{\rho}\) is then
invertible. If its rank is r, its cyclic span instead has dimension Dr.

For singular \(\rho\), that proper cyclic subspace is the actual GNS
representation of the state. The GNS vector is cyclic there by construction.
Calling all of \(HS(H)\) the state's GNS space and then calling its GNS vector
noncyclic confuses two representations. The faithful case agrees with both
descriptions.

For a general independent operator basis \(b_i\), set
\(G_{ij}=\operatorname{Tr}(b_i^*b_j)\) and
\(H_{ij}=\langle b_i\Omega,b_j\Omega\rangle\).
The squared modulus is the least generalized eigenvalue of \(Hc=\lambda Gc\).
It is not generally the least eigenvalue of H divided by one basis norm.
Orthonormalizing the operator span gives the same calculation without assuming
orthogonality or equal norms.

| Independent control | Required conclusion |
| --- | --- |
| For \(\rho=\mathrm{diag}(1/2,0,1/2,0)\), dimensions (2,2), region [1], the nonzero operator \(I\otimes|1\rangle\langle1|\) annihilates \(\Omega\). | Use the actual region's partial trace. The A region separates; the B region does not. The reversed example is also tested. |
| Bases [I,Z] and [I,I+Z] span the same diagonal qubit algebra at \(\Omega=(1,1)/\sqrt2\). | Both give \(1/\sqrt2\), also under unequal rescaling. An independent generalized eigensolver checks the result. |
| The explicit spectrum (1,1.5e-20) is positive. | Cyclicity and separation use the same resolved spectral support; no arbitrary 1e-10 threshold discards a positive eigenvalue. |
| A dense full-standardness representation requires \(D^2\) left-action matrices of size \(D^2\) by \(D^2\). At D=64 their complex entries alone need one tebibyte. | Use the D by D state and its spectrum. A D=64 test disables both dense representation constructors and checks full and noncontiguous regional standardness. |

The state validator rejects invalid densities. Dense spectra whose
support cannot be resolved at floating-point precision raise an explicit
error; a numerical answer does not certify exact rank. Explicit diagonal
positive entries, including the smallest binary64 subnormal, remain positive.

## Preserve the supplied span before computing its action

Ordinary floating-point normalization and QR can change the
operator span while still returning an orthonormal basis. For the exactly
representable inputs

```text
A = [[1, 2], [3, 4]],
B = [[1, 2 + 2^-40], [3, 4 - 2^-40]],
Omega = (1, 0),
```

`B-A` is nonzero and annihilates Omega. The separation modulus is therefore
exactly zero. A rounded chart can instead give a positive modulus above the
`1e-10` comparison used by the null-net commutant diagnostic.
For the same supplied span and `Omega=(0,1)`, the independent generalized
Gram pencil gives modulus `sqrt(9/14)`, approximately `0.801783725737273`.
Neither an orthonormality check nor a blanket zero result catches both errors.

The same numerical issue can obscure a valid algebra. With
`N=[[1,2],[2,-1]]`, the basis `[3I+N, 3I+N+2^-30 I]` spans exactly
`span{I,N}` and `N^2=5I`, but its numerical chart could fail the identity
membership check. Complex Hermitian variants also expose loss when an exact
intersection is rounded before its ill-conditioned presentation is removed.
Squaring basis entries to find their norms separately can cause underflow or
overflow for finite rescalings of ordinary matrix units.

The shared implementation treats the **supplied binary64 entries** as
its source. It does not infer a different ideal algebra from a nearby
rounded presentation:

1. Divide real and imaginary components by their common maximum magnitude
   before computing each column norm. This avoids both squaring extreme
   physical units and forming a complex reciprocal of a subnormal scale.
2. Use numerical modified Gram-Schmidt with two orthogonalization passes
   only when the normalized source columns have condition number at most
   32. This is an accuracy guard, not a rank proof.
   Otherwise perform row reduction of the transposed source columns over
   the Gaussian rationals and apply the same guard to this canonical chart.
3. If the canonical chart still cannot be resolved, perform exact Hermitian
   Gram-Schmidt before conversion and scalar normalization. Exact intersections use the same
   ordering of operations and retain their rational span for subsequent
   conservation decisions. Zero/dependent bases remain errors. Normalization
   and exact-to-float conversion refuse lost nonzero components; the numerical
   chart also refuses an entirely lost source coordinate row.

The exact preconditioning preserves the mathematical span: row operations on the
transpose are invertible changes of the independent column basis, and
Gram-Schmidt subtracts exact combinations of preceding directions. Its
nonzero squared norms are positive rationals. The numerical path has a
bounded input condition number; the exact orthogonal fallback needs only
scalar normalization. Neither inherits the original presentation's bad
conditioning.
The returned chart, GNS expectation construction, closure-test products
and ordinary singular values still use floating-point arithmetic; this
is not an interval certificate for their final error. HS projection and
distance instead use the original exact span as described below.

Trying the canonical chart first also avoids an artificial refusal. Put
`U=[[1,0],[epsilon,0]]`, `V=[[0,1],[epsilon,0]]`, `epsilon=2^-600`.
The shear `[U,U+2^-10 V]` spans the same space as `[U,V]` and the exactly
orthogonal pair `[U-V,U+V]`. Exact Gram-Schmidt of the canonical chart
introduces a coordinate of order `epsilon^2`, which underflows on conversion;
the canonical chart itself is well-conditioned and preserves its small
coordinates numerically. Its action on `(0,1)` has the exact annihilator U.
The retained real and complex controls in `test_operator_span_recharting.py`
check this zero and small projector entries on their own scale. Refusal of
one attempted chart does not prove that the supplied span is unrepresentable.
Coordinate permutations also expose partial loss within an occupied row:
for `U=[[epsilon,1],[0,0]]`, `V=[[epsilon,0],[1,0]]`, Householder QR can erase
one of the two small projector entries even though the input condition
number is near one. Reorthogonalized modified Gram-Schmidt preserves the
tested entries; checking only whether an entire coordinate row vanished
would miss the defect.

Small operator actions require a separate check. Their resolution alarm is
measured against the unit operator/reference scale, even if the entire
action nearly cancels. Below that alarm, exact rank of the action on the
original entries decides whether a nonzero annihilator exists. A kernel
returns zero. For full-rank actions whose exact orthogonal operator
directions also have orthogonal images, the squared singular values are
exact ratios of image and operator squared norms; taking their square
roots with an exact squared-midpoint comparison preserves directly evaluable
tiny positives. The final binary64 value is rounded once, with ties to even;
a rounded intermediate root can otherwise land on the wrong subnormal
neighbor. Orthogonality compares against zero in the same Gaussian-rational
domain, including for multiple nonzero action columns.
For a single operator, the same exact norm ratio is used at every scale.
Other unresolved positive actions raise an explicit error, as do positive
outputs that underflow. A small computed singular value alone is never
silently changed to zero. Masked and Boolean reference entries are rejected
before coercion.

The retained controls in `test_operator_span_charts.py` use exact
original-entry Gram pencils and projectors, together with independently
known annihilators and positive neighbors. They cover real and complex
shears, reordered bases, finite extreme scales, conditional expectations,
repeated intersections and repair gaps. In particular, nearly aligned
but distinct algebras still have their actual scalar intersection; an
unresolved slow gap cannot become a fast gap by inventing a conserved
observable. These controls validate calculations on the supplied inputs;
they do not select a physical repair law.

## The complete finite completion theorem

Let \(B\subset M_D(\mathbb C)\) be a unital complex *-subalgebra and
\(\rho>0\) a density matrix. A state b on B means **all** its observable
expectations, not just a chosen nonclosed list of observables. Define

\[
C_\rho(b)=\arg\min_{\tau\ge0,\ \operatorname{Tr}\tau=1,\ \tau|_B=b}
                   D(\tau\|\rho).
\]

The feasible set here is the entire fiber of states with B data b.
Every such fiber is nonempty and compact. Since \(\rho\) is faithful, the
objective is continuous and strictly convex, so its minimizer exists and
is unique, including on boundary fibers.

**Theorem.** The following are equivalent:

1. \(C_\rho\) is affine on the entire state space of B.
2. B is invariant under \(x\mapsto\rho^{it}x\rho^{-it}\) for every real t.
3. There is a \(\rho\)-preserving unital completely positive projection E
   with range B, fixing every element of B.
4. In Wedderburn coordinates
   \[
   H=\bigoplus_\alpha(H_{L,\alpha}\otimes H_{R,\alpha}),\qquad
   B=\bigoplus_\alpha(M_{n_\alpha}\otimes I_{m_\alpha}),
   \]
   the reference has the form
   \[
   \rho=\bigoplus_\alpha
            p_\alpha\eta_\alpha\otimes\xi_\alpha,
   \quad p_\alpha>0,\quad \eta_\alpha>0,\quad\xi_\alpha>0,
   \]
   with normalized sector weights and normalized factor states.

When they hold, E is unique and

\[
E(x)=\bigoplus_\alpha
  \operatorname{Tr}_{R,\alpha}[(I\otimes\xi_\alpha)x_{\alpha\alpha}]
            \otimes I_{R,\alpha},
\qquad
E_*(\sigma)=\bigoplus_\alpha b_\alpha\otimes\xi_\alpha,
\]

where \(b_\alpha=\operatorname{Tr}_{R,\alpha}\sigma_{\alpha\alpha}\) are
subnormalized left marginals. E acts on observables; its Hilbert-Schmidt
adjoint \(E_*\) is the trace-preserving state channel. E itself need not
preserve ordinary trace. The completion is \(C_\rho(\sigma|_B)=E_*(\sigma)\).

### Proof of necessity, including affinity

Any positive affine assignment of full states to **every** state of B,
consistent with those B data, must have the product form
\(b\mapsto\bigoplus_\alpha b_\alpha\otimes\xi_\alpha\) with fixed right
states. To see this, first take data supported in one central sector.
Positivity forces the output to have zero support in the other sectors,
including cross-sector coherences. A pure left marginal then forces the
output to be \(|u\rangle\langle u|\otimes\xi_u\).

For orthogonal u,v, compare the two decompositions of
\((|u\rangle\langle u|+|v\rangle\langle v|)/2\), using u,v and
\((u+v)/\sqrt2,(u-v)/\sqrt2\). Affinity and the diagonal blocks give
\(\xi_u=\xi_v=(\xi_++\xi_-)/2\); the off-diagonal block gives
\(\xi_+=\xi_-\). Within any orthonormal basis all right states therefore
agree. The assigned state of \(I/n_\alpha\) is independent of its
decomposition, so this common right state is independent of the basis.
Every pure state lies in such a basis; spectral decomposition gives the
formula for every mixed left state. For \(n_\alpha=1\) there is only one
normalized left state. Finally, affinity over central weights gives the
formula on all of B's state space.

Apply this to \(C_\rho\). At the reference's own data it returns \(\rho\),
since \(D(\tau\|\rho)\ge0\), with equality only at \(\tau=\rho\).
Thus affinity forces exactly the block product reference in item 4.
This argument assumes positivity and full-domain affinity, not complete
positivity. It is the finite assignment argument underlying
[Pechukas's assignment-map result](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.73.1060);
the extension over the central sectors is explicit here.
The arbitrary-system assignment result is also proved by
[Rodriguez-Rosario, Modi and Aspuru-Guzik](https://arxiv.org/abs/0910.5568).

For modular invariance, a continuous automorphism group of B fixes each
minimal central projection: it cannot continuously permute a finite set.
Therefore \(\log\rho\) has no cross-sector blocks. In a sector, its
commutator defines a derivation of \(M_n\otimes I_m\). Such a derivation is
inner on \(M_n\): equivalently, examining matrix-unit commutators shows that
off-diagonal m by m blocks of \(\log\rho\) are scalar, and differences
between its diagonal blocks are scalar. Hence
\((\log\rho)_\alpha=K_{L,\alpha}\otimes I+I\otimes K_{R,\alpha}\).
Exponentiation and normalization give item 4. The converse follows by
direct conjugation. This is the finite form of
[Takesaki's conditional-expectation theorem](https://www.sciencedirect.com/science/article/pii/0022123672900043).

For completeness, item 3 also forces item 4 without assuming the theorem.
Since E fixes B pointwise and is unital completely positive, equality in
the Schwarz inequality places B in its multiplicative domain:
\(E(bxc)=bE(x)c\). Its dual consequently depends only on the restriction
of the input state to B, and is a positive affine consistent assignment
on the entire state space. The preceding assignment argument applies.
Preservation of \(\rho\) then gives item 4.

### Sufficiency, uniqueness and exact entropy accounting

For the product reference, the displayed E is a composition of central
pinching, positive-state partial evaluation and inclusion. It is unital and
completely positive, fixes B, and preserves \(\rho\). For any expectation
with these properties, bimodularity gives

\[
\operatorname{Tr}(\rho b^*E(x))
 =\operatorname{Tr}(\rho E(b^*x))
 =\operatorname{Tr}(\rho b^*x),\qquad b\in B.
\]

It is therefore the orthogonal projection onto B for the faithful GNS
inner product \(\langle x,y\rangle_\rho=\operatorname{Tr}(\rho x^*y)\),
which proves uniqueness.

Write \(\Gamma=E_*(\sigma)\). For faithful marginals,
\(\log\Gamma-\log\rho\) belongs to B: the right-factor logarithms cancel.
Since \(\sigma\) and \(\Gamma\) have the same B expectations,

\[
D(\sigma\|\rho)
 =D(\sigma\|\Gamma)+D(\Gamma\|\rho).
\tag{1}
\]

The same equation holds on boundary fibers. Let Q be the direct sum of
the supports of \(b_\alpha\), tensored with each full right factor.
Positivity and the zero marginal weights force \(\sigma=Q\sigma Q\),
including when \(\sigma\) has central coherences. Since every \(\xi_\alpha\)
is faithful, Q is exactly the support of \(\Gamma\). On Q, the right
logarithms still cancel, and the compressed left logarithms define an
observable in B after extending it by zero outside Q. Equality of its
expectations proves (1). No logarithm of zero or assumed invertibility of
the marginal is needed. The reference term is \(Q(\log\rho)Q\); it must
not be replaced by the logarithm of \(Q\rho Q\). The two differ when Q
does not commute with \(\rho\), a case included in the regression suite.

Every \(\tau\) in the same B fiber has the same \(\Gamma\). Equation (1)
then proves that \(\Gamma\) is its unique minimum, since
\(D(\tau\|\Gamma)\ge0\) with equality only at \(\tau=\Gamma\). This proves
item 4 implies item 1 and the stated completion formula. It also shows
exactly how much relative entropy one completed repair removes. Conditional
expectations and their entropy identities are treated in
[Bardet, Capel and Rouze, sections 2 and 3](https://link.springer.com/article/10.1007/s00023-021-01088-3).

The [labelled direct sum](README.md#one-normal-form-for-the-finite-a3-objective)
turns the fixed-weight finite A3
objective into one relative entropy. The present theorem applies to such
an objective when its feasible set is precisely a complete B-data fiber.
An arbitrary OPH feasible family K may impose further constraints or lack
a joint matrix-state realization. This theorem neither removes those
constraints nor constructs that realization. In particular, the finite
classification does not derive a complete repair grammar, retained algebra or
reference state from the [canonical premises](../../docs/AXIOM_REFERENCE.md).

The [native fixed-seam attachment](../../model/primitive_observer/DERIVATION_COMPLETION.md)
identifies E with this expectation on the declared 66-dimensional binary
shell and tracial reference. C has the same completion endpoint but different
finite-event coherence. This does not extend the theorem to all carrier
sectors or select a complete recorded instrument.

## Three exact controls that distinguish the premises

**Nonlinear completion.** Take B to be the diagonal qubit algebra and
\(\rho=\begin{pmatrix}1/2&1/4\\1/4&1/2\end{pmatrix}\).
The diagonal data (1,0) and (0,1) each have only one positive completion.
Their mean completion is \(I/2\). At their mean data, however, the unique
relative-entropy minimum is \(\rho\), because it already has those data.
The affine candidate loses by
\(D(I/2\|\rho)=\tfrac12\log(4/3)>0\). Thus a reference alone does not make
A3 completion a linear repair channel. There is no optimizer in this test.

**A weighted projection need not be positive.** Let
\(\rho=(I+\tfrac12 X\otimes X)/4\) and \(B=M_2\otimes I_2\).
The \(\rho\)-GNS orthogonal projection is unital, idempotent, fixes B,
has range B, preserves \(\rho\), and is GNS self-adjoint. Yet it sends the
positive input \(I+Z\otimes X\), whose eigenvalues are 0 and 2, to
\(I+\tfrac12 ZX\otimes I\), which is not even Hermitian. Here X and Z are
Pauli matrices. For the projected left factor y, the orthogonality equations
are \(y\rho_L=\operatorname{Tr}_R(x\rho)\); here \(\rho_L=I/2\), giving
the stated image. The tests check those equations with exact rational
matrices against a full Pauli basis. This prevents treating weighted
least squares as an automatic quantum channel.

**An operator system need not support a quantum repair.** Retain only the
span of I,X,Z with tracial reference \(I/2\). Entropy maximization sets the
unconstrained Y coefficient to zero. The resulting map
\(T(a)=(a+a^{\mathsf T})/2\) is affine, positive, unital and idempotent,
but its unnormalized Choi eigenvalues are
\((-1/2,1/2,1/2,3/2)\). It is not completely positive. Its range is closed
under adjoint but not multiplication, so it is not B as required above.
The algebra validator rejects that range before any repair construction.

## Noncommuting primitive repairs still have a common limit

Let \(E_m\) be finitely many reference-preserving expectations on one
common finite workspace, with strictly positive supplied rates \(\gamma_m\).
Define

\[
\mathcal L=\sum_m\gamma_m(E_m-I),\qquad B_\cap=\bigcap_m B_m.
\]

Each \(E_m\) is a GNS orthogonal projection. Therefore

\[
\langle x,-\mathcal Lx\rangle_\rho
 =\sum_m\gamma_m\|(I-E_m)x\|_\rho^2,\qquad
\ker\mathcal L=B_\cap.
\tag{2}
\]

Positivity of every rate makes the kernel equality exact. No commutation
between different \(E_m\) is needed. Put
\(\gamma=\sum_m\gamma_m\) and
\(P=\sum_m(\gamma_m/\gamma)E_m\). Then
\(e^{t\mathcal L}=e^{-\gamma t}\sum_{k\ge0}(\gamma t)^kP^k/k!\) is unital,
completely positive and reference-preserving. The finite spectral theorem
and (2) show convergence to the GNS projection \(E_\cap\) onto \(B_\cap\).
That intersection is modular invariant, so the projection is its unique
reference-preserving conditional expectation.

When the orthogonal complement is nonzero, let \(\lambda>0\) be the least
eigenvalue of \(-\mathcal L\) on it. Then
\(\|(e^{t\mathcal L}-E_\cap)x\|_\rho
\le e^{-\lambda t}\|(I-E_\cap)x\|_\rho\).
When the intersection is the full algebra there is no positive gap to
report. The dual semigroup is a state channel, so relative entropy to
\(\rho\) is nonincreasing by data processing. This does not identify the
GNS gap with an entropy decay constant.

An exact control uses qubit dephasing in Z and
\(N=(3X+4Z)/5\), each at rate 1/2 with \(\rho=I/2\).
These expectations do not commute. Their common algebra is the scalars
and the generator spectrum is \(0,-1/10,-9/10,-1\).
SymPy verifies the rational characteristic spectrum independently of the
numerical gap routine. Explicit matrix exponentials check complete
positivity, entropy decrease and the corresponding norm bound.
This supplies a concrete noncommuting example, not just a commuting
special case.

## Implementation, independent checks and downstream effects

[algebras.py](algebras.py) validates a complete independent basis, its
identity, adjoints and every pairwise product. It embeds operators in the
declared tensor slots, computes intersections, and tests modular invariance
by all commutators with \(\log\rho\). This finite generator criterion is
equivalent to invariance under all modular times.

[expectations.py](expectations.py) separates the weighted projection
candidate from a validated expectation. Its recognizer accepts the **full
linear superoperator**, not a callable tested on a few random states.
It independently recomputes the full Choi matrix, unitality, reference
preservation, range, fixation of B, idempotence and GNS self-adjointness.
It never decides acceptance by comparing a supplied map with the producer's
projection formula. The observable map and its dual use explicit row-major
vectorization; complex unitary examples test the convention.

The numerical tolerance is finite, positive and bounded; it cannot be set
to infinity or a Boolean to disable checks. Invalid states, singular
references, dependent bases, empty families, omitted/zero/negative rates,
nonfinite matrices, diagnostic overflow and unresolved Gram matrices fail
explicitly. Validation remains active under optimized Python.
Vector norms are evaluated in complex128 arithmetic before checking
normalization. The int64 vector \((2^{32},1)\)
has a wrapped squared norm of one in int64 arithmetic; float32
could also round the squared norm of \((1,10^{-4})\) to one. Regression
controls reject both, along with unsigned-integer and complex64 variants.
These guards are not interval arithmetic: deviations below the declared
tolerance can pass. Channel conditions and spectral gaps remain numerical
diagnostics. The [fixed-point calculation](REPAIR_FIXED_POINTS.md) uses
exact elimination on the retained supplied spans. Two distinct, almost
parallel qubit repairs have only a scalar conserved algebra. Their true
slow mode has rate of order the squared angle and must not be deleted from
the complement. Unresolved positive gaps still raise. The calculation gives
the complete finite spectrum and record lifetime, and keeps exact span
data separate from its numerical basis. A numerical construction is not
a source-law receipt.

### Preserve observable actions and residuals

The supplied span also determines projection, not just intersection. If
\(C\) has the independent original vectorized operators as columns, the
ordinary-trace Hilbert--Schmidt projector is
\(P=C(C^*C)^{-1}C^*\). Indeed, \(Px\) belongs to the span and
\(C^*(x-Px)=0\), which uniquely characterizes the orthogonal projection.
[linear_actions.py](linear_actions.py) evaluates this factored expression
over the exact rational real and imaginary parts of the supplied binary64
entries. It caches the coordinate factor, rather than a dense ambient
projector, and refreshes that cache when an intersection supplies new exact
columns. Clearing a common denominator from \(C\) leaves \(P\) unchanged;
the Gram inverse is stored as an integer numerator and denominator. This
keeps repeated contractions over Gaussian integers and postpones rational
division until output, avoiding repeated fraction reduction in closure
checks. Dense exact calculations remain more expensive than floating chart
calculations; their cost depends on both span dimension and integer bit
length. The full matrix algebra needs no factorization.

With \(B=\mathbb C I_4\), \(\rho=I_4/4\) and
\(X=\operatorname{diag}(2^{56},-2^{56},1,0)\), the HS projection and
the expectation must give \(I_4/4\). The associated rate-one repair action
must retain its final diagonal response \(1/4\). Numerical contraction can
lose these responses even when the map passes the complete channel recognizer.
Real and imaginary controls retain them; finite results also survive overflowing
or underflowing intermediate products.

Distance is computed as \(\sqrt{(x-Px)^*(x-Px)}\), with subtraction and the
sum of squares exact until the final scalar rounding. Subtracting a rounded
projection is insufficient: for \(B=\mathbb C I_3\) and
\(X=\operatorname{diag}(1,1+2^{-52},1)\), the projection rounds to \(I_3\),
but the distance is \(2^{-52}\sqrt{2/3}\), not \(2^{-52}\). A large scalar
part plus a tiny off-diagonal perturbation must likewise retain the tiny
distance; QR roundoff must not supply a much larger fictitious residual.
Algebra validation compares the exact squared distance with its squared
tolerance, so a residual below the floating-point range can still be
correctly accepted. This preserves the existing tolerance-based closure
contract; it does not certify exact algebra closure.

For an independently supplied superoperator, `apply_map` preserves that
map's entries and the row-major/adjoint convention. It uses a floating dot
product only when separate real and imaginary cancellation guards and
range checks permit it; otherwise it evaluates the original entries
exactly. Upper-range magnitude estimates also trigger exact evaluation:
even an all-positive sum can round down to the largest finite float when
the original sum lies at the overflow midpoint. Real/imaginary and dual
controls check both sides of this boundary, retaining representable results.
The ordinary floating path is not a correctly-rounded or interval
certificate. Projection components and exact fallback map components are
rounded once; distance uses an integer midpoint comparison to round its
rational square root without an intermediate floating square. A final
nonzero component or distance that rounds to zero, or a result that
overflows on final rounding, raises explicitly. A scalar distance can
therefore succeed even when individual projected components cannot be
represented.

The controls in [test_linear_actions.py](test_linear_actions.py) and
[test_projection_residuals.py](test_projection_residuals.py) use independent
rational contractions, analytic projections and high-precision scalar
norms. They cover mixed scales, complex dual maps, retained intersections,
and final-range refusals. Observable evaluation, channel recognition and
spectral support are distinct checks; none selects a physical repair law.

| Component | Current contract |
| --- | --- |
| Null-net standardness | Uses the actual region, a basis-independent modulus and resolved spectral support. Unresolved support raises instead of silently deciding rank. |
| Collar modular alignment | Its complete-basis criterion is checked against existence of the finite expectation. Aligned product blocks pass and generic incompatible blocks fail. |
| Finite MaxEnt/A3 | Identifies exactly which full-marginal problems give affine repair and proves the Pythagorean accounting; arbitrary feasible families retain their own constraints. |
| Source selection | Supplies finite existence, recognition and noncommuting convergence for specified algebras and reference. Protected records, conserved functionals, checkpoint instruments, grammar completeness, refinement and clock attachment require their own evidence. |
| Classical consensus projection | A separately specified row-stochastic three-state transition table reproduces weighted fiber resampling. Replacing stationary weights by uniform weights fails. The existing demand for independent transition extraction is retained. |
| Source operator-join evidence | Its exact source-specific verifier remains independent. A generic constructed channel does not replace its source identification or Lean proof. |

The one-particle zero-extension helper in the null-net receipts is also a
different operation from tensor-identity embedding and remains separate.
These finite results do not identify a physical model.

The regression suite uses explicit tensor indices, independent generalized
Gram eigenvalues, weighted partial traces, independent Choi blocks, SciPy
matrix logarithms, exact rational controls and adversarial supplied maps.
A central-sector test includes a pure state with cross-sector coherence
and a singular completed marginal, checking (1) against scalar logarithms.
Separate controls treat a singular marginal whose support does not commute
with the reference, and noncommuting repairs preserving a complex nontracial
reference after a global unitary change of coordinates.
Run the live consumers and their tests with the pinned requirements:

~~~text
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  code/quantum_information code/geometry code/collar_alignment code/maxent
~~~

The dedicated finite-quantum-information workflow runs these tests on
Linux.
