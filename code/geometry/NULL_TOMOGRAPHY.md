# Finite null tomography: reconstruction, consistency and resolution

This audit repairs the numerical realization of the existing finite
null-tomography theorem. The old solver returned a zero-residual tensor from
one direction or twelve copies of it, even when a nonzero trace-free tensor
was completely invisible. It also accepted timelike designs. Normalizing a
spatial vector by `sqrt(dot(d,d))` made `(1e200,0,0)` produce the **timelike**
four-vector `(1,0,0,0)` and made `(1e-200,0,0)` produce nonfinite data.

The repair in [null_tomography.py](null_tomography.py) separates three facts:
the design resolves nine components; the sampled charges lie within a declared
distance of the tensor image; and the inverse has a measured sensitivity.
Rank and sensitivity cannot be inferred just from a small residual. This is a focused contribution
to [the existing-theory audit, #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).

## The finite theorem and its angular meaning

Work in a fixed observer frame with `eta = diag(-1,1,1,1)` and `k=(1,n)`,
`|n|=1`. A real symmetric tensor has null charge

\[
q_T(n)=T_{00}+2T_{0i}n_i+T_{ij}n_i n_j.
\]

Write the spatial block as `S = a I + Q`, with `tr Q = 0`. Then
`q_T = (T00+a) + 2 T0i n_i + n^T Q n`. These are the degree-zero, degree-one,
and degree-two spherical harmonics, of dimensions `1+3+5=9`. They are
independent: the odd part fixes the linear term; averaging the even part over
the sphere fixes its constant; a trace-free quadratic form constant on the
sphere vanishes. Thus the only invisible symmetric tensors are `phi*eta`.
The gauge `-T00+T11+T22+T33=0` removes exactly that ambiguity.

For a finite normalized family, let `A` evaluate a Frobenius-orthonormal basis
of this nine-dimensional space. At rank nine:

* The unique least-squares coordinates are `x_hat=A^+q`.
* A sampled family has an exact tensor representative iff `w^T q=0` for
  every `w` satisfying `A^T w=0` (the dependent-family relations).
* With `r=q-A x_hat`, orthogonality gives
  `||q-Ax||^2 = ||r||^2 + ||A(x-x_hat)||^2`. Hence `||r||` is the exact
  distance to the sampled tensor image. For nonzero residual,
  `w=r/||r||` is a unit incompatibility witness, with `w^T q=||r||`.
* If `q=A x_true+e` and `||e|| <= delta`, then
  `||T_hat-T_true||_F <= delta/s_min(A)` in the trace-free gauge. This is
  sharp: noise along the least singular left vector attains equality.

These statements follow from orthogonal projection and the singular-value
decomposition. The tests verify the equality case of the bound, as well as
an independent exact incompatibility witness; they do not merely check that
a producer accepts its own data. The error bound concerns charge error with
the directions fixed, and excludes error in the geometry or in the numerical
solve. The norm and conditioning refer to the declared observer frame;
Frobenius norms are not Lorentz invariant.

## Reusing the proved frame and exposing a sampling blind spot

The deterministic frame is the one already formalized in
[Tensor.lean](../../Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean):
six spatial axes and the three body diagonals
`(1,1,1)/sqrt(3)`, `(1,1,-1)/sqrt(3)`, `(1,-1,1)/sqrt(3)`.
The Lean `tomographyDecoder_charge` and
`nine_charge_metric_ambiguity` theorems already prove its inverse and metric
ambiguity. An independent SymPy calculation builds the charge map from the
bilinear form, obtains determinant `8192/27` in the existing nine-coordinate
basis, inverts it exactly, and tests the numerical solver on every coordinate.
The numerical implementation uses an orthonormal basis instead, so that its
singular values measure Frobenius error directly. No new formal axiom is used.

Nine independent rows interpolate **any** nine scalar readings: they impose
no compatibility test on the readings themselves. An explicit obstruction is
the cubic harmonic

\[
f(n)=n_x(n_y^2-n_z^2).
\]

It vanishes on all nine frame directions but takes value `-2/9` at
`n=(2/3,1/3,2/3)`. Its first nine readings reconstruct the zero tensor; its
tenth excludes every symmetric tensor. The tests compute the tenth row,
its unique dependent-family relation, its nonzero violation, and the exact
least-squares distance using symbolic arithmetic before comparing the
numerical fit and witness. This counterexample identifies the logical
boundary of finite sampling. Even many sampled directions do not by themselves
establish tensoriality on all unsampled directions. The paper's theorem about
a whole directional family requires its dependent-family relations across
that family, not just interpolation on one nine-element subset.

## Numerical contract

`fit_null_charges(q, rays)` accepts finite real binary64 inputs, one scalar
charge per nonzero Minkowski-null four-vector, and a resolved rank-nine design.
It accepts future or past representatives. It checks the spatial norm of
`k/k0` against one within `64*eps` and rejects a smallest singular value at or
below `64*eps*max(n,9)*s_max`. This is a declared numerical resolution cutoff,
not a new geometric theorem. Small departures from exact nullness within that
validation tolerance remain ordinary input roundoff.

The raw inputs obey `(k,q) -> (c*k,c^2*q)`. The solver normalizes both to
`(k/k0,q/k0^2)`, so arbitrary ray representatives cannot alter rank, residual,
conditioning or consistency. Normalization uses exact rational division of
the supplied binary numbers followed by one binary64 rounding. Charge
evaluation likewise sums exact quadratic products before rounding; this
preserves `2e-200` when `T00=-T22=1e200`, `T02=T20=1e-200`, and `k=(1,0,1,0)`.
It also avoids forming overflowing or underflowing ray squares. This modest
exact work is appropriate for the finite audit, not intended as a large-data
simulation kernel.

`fit.require_consistent(error_budget=...)` replays the returned tensor against
the **original** binary inputs. It compares the exact rational sum
`sum(((q_i-T(k_i,k_i))/k_i0^2)^2)` with the exact squared **explicit** finite
nonnegative error budget. Thus charge/ray normalization and solve roundoff
cannot cause an out-of-budget tensor to be accepted. The displayed
`residual_norm` is a numerical square root of this sum; acceptance does not
rely on that rounded number. There is no hidden absolute tolerance. NaN or
infinite budgets are errors. A rejection can reflect insufficient numerical
resolution and alone does not prove that no better tensor exists. The result
also supplies singular values, the sharp data-noise amplification estimate,
a unit residual witness, its charge contraction, and its design defect
`||A^T w||`. These SVD diagnostics are floating-point estimates, not interval
certificates; a witness whose design defect is appreciable is not an exact
dependent-family proof. Exact proofs in the tests and Lean have separate
authority. In particular, roundoff on a consistent family can produce a tiny
residual with an unreliable normalized witness. The exact replay establishes
the returned tensor's finite error budget; it does not certify SVD optimality,
an inverse error bound or the physical interpretation. Result arrays are
backed by immutable bytes and exposed through separate array views so neither
their values nor shape/dtype metadata can be edited under an already checked
residual. Acceptance and the legacy wrapper return independent writable copies.

Malformed data, complex or nonfinite entries, lossy conversion to binary64,
zero directions, nonsymmetric
tensors, insufficient or unresolved designs, unrepresentable nonzero charges,
and rescaling that erases a nonzero input component fail closed. In particular,
two unequal integer charges above `2^53` cannot silently become equal through
float conversion, including mixed Python integer/float sequences. Residual norms
use a scaled norm so a violation of size `1e-200` does not become zero by
squaring. Reconstruction remains an ordinary binary64 solve; large error
amplification is reported, not interpreted as physical evidence.

The maintainer-style follow-up audit reproduced two further defects before
repair: the intermediate solver residual understated the returned tensor's
residual (`4.0283e-15` versus `4.7972e-15` in the retained control), and mutable
NumPy arrays allowed an accepted zero tensor to be edited while keeping its
zero residual. Exact replay and immutable result storage close these paths.
Controls compare both sides of the acceptance boundary, down to adjacent
binary64 budgets, against an independent exact contraction.
The same audit checks mixed-type conversion before NumPy promotion and
computes the metric projection exactly before rounding: dividing each trace
summand first previously erased a minimum-subnormal metric coefficient.
A separate 120-digit QR solver checks thirty families across three sample
counts, three charge scales and individually rescaled past/future rays. It
builds the coordinate matrix directly from the bilinear form, without the
producer's basis, SVD or charge evaluator.

The existing `reconstruct_from_charges` remains a **diagnostic** returning
`(tensor, residual)` and can report an inconsistent full-rank family. Its
residual now uses normalized rays; earlier callers all used `k0=1`.
`design_matrix` and `charges_of` retain raw quadratic weights.
The compatibility wrapper now rejects incomplete designs instead of silently
returning an underdetermined pseudoinverse.

## Impact and reproduction

The Einstein closure receipt now uses the existing exact frame plus three
additional directions, enforces a stated numerical budget on its consistent
family, and reports sensitivity and a witness for its inconsistent family.
The entropy, central first-law, MaxEnt and baseline countermodel code in the
same file is unchanged. The claim registry, paper statements, Lean results,
frozen predictions and pinned historical receipts are unchanged: the repair
makes the live numerical implementation respect the existing finite theorem.
It supplies no local-stress identification, Ward identity, universal coupling,
or physical Einstein branch.

From the repository root, with the pinned requirements installed:

```sh
python -m pytest -q code/geometry/test_null_tomography.py code/geometry/test_einstein_closure_receipts.py
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent
```

The existing finite-quantum-information workflow runs the full second command
on Linux. The focused controls include the exact frame and
cubic obstruction, worst-case noise, Lorentz covariance modulo the metric,
ray rescaling and permutation, deficient and nearly collapsed designs,
extreme units, cancellation, and malformed inputs. The extended-exponent
conversion control runs where `longdouble` has a wider exponent range;
Windows skips that control because its two floating types share the range.
