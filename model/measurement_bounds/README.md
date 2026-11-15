# Measurement bounds that do not assume a Gaussian preparation

This packet gives predictions for a complete native readout interface that
hold for **every positive state**, including non-Gaussian states and states
after arbitrary dynamics. It separates a declared detector-gain problem from
the preparation-family problem examined by the
[six-probe diagnostic](../empirical_diagnostics/README.md).

This packet uses no additional experimental outcomes. The anharmonic comparison and its
score are those of the anharmonic packet. The content here is a state-independent calibration
control, a necessary relation between two characteristic-function reads, and
an outward bound on the relation's defect for the actual finite Kraus matrices.

## An additional prediction without choosing a preparation family

Let `U` be the unitary measured by one complex characteristic read. Measuring
the same quadrature at twice its coordinate gives `U^2`. For any normalized
positive state, write `z=<U>` and `w=<U^2>`. Then

\[
 |w-z^2|\le 1-|z|^2.
\]

Thus measuring both quadratures of `z` restricts the next complex read to a
disk with center `z^2` and radius `1-|z|^2`. No displacement, temperature,
Gaussian family, Hamiltonian, oscillator occupation cutoff or physical clock
appears in this inequality. The common quadrature and relative doubling are
interface premises. The relation holds for the unrestricted oscillator and
for each finite Hermitian quadrature separately; it needs no canonical
commutation identity at a cutoff.

For illustration, the unchanged native baseline gives
`z=0.6253565912-0.6523263793i`. Conditional on that first complex read, **every
state** must satisfy approximately

```text
-0.999272631 <= Im chi(0.8i) <= -0.632473773.
```

The baseline prediction `-0.666242523` lies inside this interval. The existing
non-Gaussian control has a different first complex read and a different
allowed interval; it also passes. These numbers illustrate the bound with
source-generated expectations, not measured frequencies. Real finite-shot
readings require the uncertainty allowance below.

This bound is weaker than the exact Gaussian identity
`chi(2ih)=|chi(ih)|^2 chi(ih)^2`, but it requires much less knowledge of the
preparation. Failure of the Gaussian identity can indicate an inadequate
preparation family. Failure of the universal bound, after justified errors,
cannot be repaired merely by choosing another state or Hamiltonian for the
same repeated experiment. It instead challenges some joint premise: stable
preparation across probe settings, their common unitary/power relation,
detector calibration, error assumptions, or the positive-state description.
The test alone cannot select which premise failed.

## The finite theorem attachment and compiled-word correction

The existing
[`EventAlgebra.expectation_nonneg` and `expectation_one`](../../Lean/EventAlgebra/StateExpectation.lean)
supply positivity and normalization of the finite trace functional. Positivity
of `<X*X>` implies Cauchy-Schwarz for its Gram form. With
`X=U*-conj(z)I` and `Y=U-zI`, it yields the displayed moment inequality.
This short specialization is a direct mathematical argument; no Lean
theorem is claimed.

The executable does **not** assume its floating source matrices satisfy the
ideal group law exactly. From the actual source-contracted pointer Kraus
matrices, including the spin record, it forms the exact-binary operators

```text
M = E0-E1,    A=M_real(h)+i M_imag(h),
B=M_real(2h)+i M_imag(2h).
```

Here each effect sums both retained spin branches. Outward interval matrix
arithmetic bounds

```text
delta = max(||A* A-I||, ||A A*-I||),
eta   = ||B-A^2||.
```

For any normalized positive input, Cauchy-Schwarz applied to
`X=A*-conj(<A>)I`, `Y=A-<A>I` gives

\[
 |\langle B\rangle-\langle A\rangle^2|
 \le 1-|\langle A\rangle|^2+\delta+\eta.
\]

Both noncommuting norm defects are included. The implementation also bounds
each instrument's effect-sum defect, the null controls and phase-reversal
identities. Frobenius enclosures bound the required operator norms. Every
elementary operation and square root is widened outward; BLAS is not used
inside these certificates. Exact-rational tests check the enclosure kernel,
and a deliberately wrong doubled-coordinate instrument produces a large
detected defect.

For the retained 32-level execution, the complete spin-plus-pointer correction
`delta+eta` is below `8.31e-13`; all instrument effect-sum defects are below
`7.76e-14`. The separately contracted noisy probability laws agree with the
full CP execution within `3.34e-16`. These numerical scales describe the
finite implementation, not experimental uncertainties.

These certificates describe the exact binary matrices returned in the fresh
execution. They do not certify the floating propagation of density matrices,
the implementation of an ideal source word in hardware, or physical service
errors. The small effect-sum defects are retained rather than silently
renormalizing outcomes. Independent matrix contractions check the complete
executed probabilities separately. A replay recomputes its own enclosures;
stored radii are not attached to a different binary execution.

Finite normalization must also be kept explicit. For exact binary effects
with `M=E0+E1`, the noisy mean is
`b<M>+g<E0-E1>`, not literally `b+g<E0-E1>` when `M` differs from identity.
The reported effect-sum bounds control that difference. The `delta+eta`
certificate applies to the raw actual `A,B`; using null-calibrated floating
means additionally requires charging effect-sum, null-control and calibration
errors, as well as propagation error. The receipt's calibrated moment-witness
values are therefore **diagnostics**, not an end-to-end certified rejection
test. Its ideal calibration and sampling formulas below assume a complete
ideal interface; physical use must add the justified interface errors.

## Separate a stationary bit-readout channel from preparation

The declared detector error model has probabilities `e0` for reporting one
when the ideal bit is zero, and `e1` for reporting zero when it is one. It is
a supplied, state-independent classical channel, constant across settings and
independent of the recorded spin. Its mean obeys

```text
m_observed = b + g m_ideal,
b=e1-e0,    g=1-e0-e1.
```

At zero probe coordinate the two opposite phase settings have ideal means
`+1` and `-1`, regardless of the oscillator state. Their observed means
therefore identify `b=(m_plus+m_minus)/2` and `g=(m_plus-m_minus)/2`.
Positive nonzero gain is required for the chosen orientation. The imaginary
null setting must give `b`. Opposite phases at any nonzero coordinate must
sum to `2b`. These redundant controls can reject the stationary-channel
assumption; they do not establish a universal hardware calibration model.

The retained example supplies `e0=.05`, `e1=.10`. Actual complete source
execution gives null means `.9` and `-.8`, recovering `b=.05` and `g=.85`
for both the original Gaussian and non-Gaussian preparations. A second
executed control changes the detector offset by `.04` only for the two
imaginary phase settings. Their reversal sum acquires residual `.08`, exposing
that setting dependence without changing preparation or dynamics.

This identifies an affine binary-readout gain and offset under the stated
conditions. It does **not** identify the physical displacement scale. Replacing
`h` by `s*h` preserves both the unitary moment inequality and its doubling
relation. The probe-scale gauge and the archived `.636` versus
`sqrt(2)*.4` provenance issue are therefore unresolved. Raw counts, phase
settings and independently authenticated coordinate metadata are needed.

## What executes

The [protocol](protocol.json) fixes nine settings: three zero-coordinate
controls, the real and imaginary `h=.4` probes with their phase reversals,
and both quadratures at `2h=.8`. Every setting uses a separate fresh
preparation; probes are not sequential destructive reads of one state.

For each example, the computation executes the unchanged native Fock and
displacement preparation, the native three-pulse spin pointer, and the native
oscillator pointer. The supplied bit-error channel is implemented by summing
its Kraus alternatives as CP branches. The common primitive runtime retains
all four `(spin, observed detector bit)` outcomes and their terminal quantum
states. There is no postselection. The non-Gaussian example is the same fixed
positive mixture used by the [empirical diagnostic](../empirical_diagnostics/README.md); no target values or later
trajectory observations enter its construction.

The host computes complete instrument laws rather than fabricating shots.
The receipt binds preparation provenance, emitted word hashes, source
dependencies, complete probabilities and the actual operator certificates.
The examples have zero evolution duration; the mathematical measurement
constraints hold for any state supplied to that same interface.

## Sampling errors and prospective use

For `n` independent fresh binary trials per setting and `k` settings, Hoeffding
and a union bound give a simultaneous raw-mean radius

\[
 r=\sqrt{2\log(2k/\alpha)/n}
\]

with failure probability at most `alpha`. If the two null calibrations are
also within this radius, then for each corrected mean
`mu_hat=(m_hat-b_hat)/g_hat`,

```text
|mu_hat-mu| <= 3r/|g_hat|.
```

The calibration uncertainty is charged, rather than treating the inferred
gain as exact. For a positive true gain certificate require `g_hat>r`.
Two quadratures give complex-moment radius `sqrt(2)*3r/|g_hat|`.
These deliberately conservative bounds assume the stationary bit channel,
stable repeated preparation and independent trials; real apparatus error
needs its own additional allowance.

With `k=9`, `alpha=.01` and observed gain `.85`, the corrected complex radii
are approximately `1.367`, `.4321` and `.1367` at 200, 2,000 and 20,000 trials
per setting. The paper's stated 200 shots would be insufficient for a useful
confidence-level version of this conservative calibration-and-bound test,
even if authenticated counts were available. This is a different, more
model-independent task than choosing between the two locked
preparation candidates.

For complex estimates `z_hat,w_hat` with justified radii `epsilon1,epsilon2`,
the executable uses the robust necessary condition

\[
 |\hat w-\hat z^2|+|\hat z|^2-1
 \le \delta+\eta+\epsilon_2
       +4|\hat z|\epsilon_1+2\epsilon_1^2.
\]

The scalar formula and Hoeffding guarantee are mathematical bounds; their
displayed logarithm evaluations are floating approximations, separately
labeled from the outward matrix certificates. No statistical inference is
made from source-generated expectations or the deposited derived positions.

The next useful measurement is the complete null/phase-reversal calibration
and both complex probes on contemporaneous, identically prepared trials.
Passing these checks supports a consistent measurement interface. It does
not recover the physical action, derive the prepared state, explain the anharmonic
trajectory discrepancy or distinguish OPH from ordinary quantum mechanics.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.measurement_bounds --check model/measurement_bounds/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/measurement_bounds
```

Use `--output <path>` to write a new receipt. The output identifies the
remaining preparation, interface and detector-channel premises explicitly.
