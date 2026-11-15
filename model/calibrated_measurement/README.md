# From complete counts to a calibrated measurement prediction

This packet closes the calibration and finite-sampling gap identified in
[measurement bounds](../measurement_bounds/README.md). It gives a conditional
confidence test for the same preparation-independent moment relation using
**six measured settings**, and can predict the last two settings from the
first four. The calculation charges finite instrument completeness,
imperfect null controls, calibration uncertainty and declared apparatus
errors together. It uses exact rational scalar arithmetic and outward square
root/logarithm bounds.

It uses no additional laboratory outcomes. The anharmonic comparison,
its unresolved coordinate provenance and its RMSE are those of that packet. This
result is a usable prospective calculation, not a resolution of that discrepancy
or an empirical distinction between OPH and ordinary quantum mechanics.

## What the next measurement can test

Measure complete `(spin, detector bit)` counts at these fixed settings:

| Setting | Role |
| --- | --- |
| `null_plus`, `null_minus` | Zero-coordinate, opposite-phase gain and offset calibration |
| `h_real`, `h_imag` | First complex moment at the common coordinate `h` |
| `2h_real`, `2h_imag` | Second complex moment at the doubled coordinate |

Each setting uses independently prepared trials, not successive destructive
reads of one state. The four non-null settings must concern the same prepared
state. The null controls are bounded uniformly over states and need not use
that particular state. The detector's unknown classical confusion channel
must stay the same across settings and be independent of the spin record.
All four joint outcomes are counted; summing the spin rows is a marginal,
not postselection. Fixed trial budgets are required; optional stopping or
selecting the most favorable probe after seeing counts needs another bound.

The three imaginary-null/phase-reversal controls are useful checks
on stationarity. They are not required for the six-setting implication, and
omitting them does not establish stationarity. The producer compiles
their operators to reuse the existing nine-setting operator certificate;
only six settings enter this measurement family or its trial budget.

As a concrete illustration, rounding the unchanged native baseline's means
to a **software count fixture** at 20,000 trials per setting gives:

| Setting | `(spin 0, bit 0)` | `(spin 0, bit 1)` | Both spin-1 counts |
| --- | ---: | ---: | --- |
| `null_plus` | 19,000 | 1,000 | 0, 0 |
| `null_minus` | 2,000 | 18,000 | 0, 0 |
| `h_real` | 15,816 | 4,184 | 0, 0 |
| `h_imag` | 4,955 | 15,045 | 0, 0 |

If actual independent counts had these values, the stated stationary-channel
premises held, and service errors were zero for the named complete mathematical
interface below, the prospective 99% family guarantee would imply

```text
-0.853255 <= expected observed mean at 2h_imag <= -0.034184
```

The displayed endpoints are rounded outward. This is a bound on the future
setting's expected detector mean, not on one future bit or its sample mean.
Future counts receive their own simultaneous sampling allowance. The fixture
is not an observation or a random sample, so no experimental confidence
claim is made from it. Its assigned spin-zero counts are bookkeeping for a
software control, not an independently certified simulated trajectory.

At 200 or 2,000 trials per setting the conservative interval is too broad to
give this sign prediction. At 20,000 trials, the first complex-moment error
radius is about **0.132906**, including the inferred calibration. The
nine-setting design has corresponding radius about 0.1367; this packet
contributes the complete error accounting and executable count interface.

## The actual computation and its finite comparison

Fresh execution contracts the existing native source words for the oscillator
pointer and the three-pulse spin pointer on the 32-level oscillator code.
Each detector outcome has both spin Kraus branches, acting on a 64-dimensional
carrier. The unchanged outward matrix kernel encloses the exact-binary effects

```text
E_(s,b) = sum_spin K_(s,b,spin)* K_(s,b,spin),
T_s = E_(s,0)+E_(s,1),       M_s = E_(s,0)-E_(s,1),
A = M_h_real + i M_h_imag,   B = M_2h_real + i M_2h_imag.
```

The replay checks fixed rational guarantees on those fresh arrays:
`||T_s-I|| <= 10^-12`, both null-observable errors at most `10^-12`, and

```text
max(||A* A-I||, ||A A*-I||) + ||B-A^2|| <= 10^-11.
```

These are conservative exported guarantees; the existing sharper raw
enclosures are recomputed rather than read from a receipt. They certify the
**returned binary Kraus matrices**. The existing detector and spin-pointer
producers use their guarded polar representation of contracted native words.
This packet adds no projection or normalization of those returned matrices;
the raw-word-to-representation telescope is a floating diagnostic, not
a rigorous source-implementation bound. That discrepancy requires a justified
additional service allowance if the raw pulse execution is the target.
Native word metadata and source hashes bind the fresh construction. Base-array validation prevents
an array subclass from concealing entries through overridden operations.

Binary source effects need not sum to identity exactly. To define a genuinely
complete linear comparison instrument, let

```text
L_(s,b,spin) = K_(s,b,spin) T_s^(-1/2).
```

For `c_s=||T_s-I||<1`, positive finite-dimensional spectral calculus gives
`T_s>0` and `sum L*L=I`. All outcomes are kept. This comparison is a linear
quantum instrument, not state-dependent probability renormalization or
conditional success selection. Its effects are used as a mathematical target;
the code neither executes an extra native inverse-square-root operation nor
claims a source implementation of that correction.

Let its binary observable be `Mtilde_s=T_s^(-1/2) M_s T_s^(-1/2)`.
Positivity of `E0,E1` gives `||M_s||<=1+c_s`, while spectral calculus gives
`||T_s^(-1/2)||<=1/sqrt(1-c_s)` and
`||T_s^(-1/2)-I||<=1/sqrt(1-c_s)-1`. Therefore

```text
||Mtilde_s-M_s|| <= kappa_s = c_s(1+c_s)/(1-c_s).
```

This comparison error is included explicitly. Any actual apparatus discrepancy
from this complete comparison, including physical implementation errors, must
enter through a justified additional observed-mean bound `beta_s`.
A complete-output total-variation bound `epsilon_s` can supply
`beta_s=2 epsilon_s`; missing bounds are `null`, never zero.

The producer also executes the baseline preparation and the six
complete readouts in the shared primitive runtime. Its floating expectations
are retained as **diagnostics only**, used to make clearly labeled rounded
count fixtures. No floating propagation result enters the confidence theorem.
The theorem covers every normalized positive state on the finite carrier,
including a state correlated with an untouched reference after reduction.
An unrestricted oscillator experiment needs a justified interface/cutoff
comparison in its service bounds; this packet does not supply one.

## Joining the existing theorems

The finite positive state expectation in
[`StateExpectation.lean`](../../Lean/EventAlgebra/StateExpectation.lean)
and Gram Cauchy-Schwarz give the moment bound. For the actual binary
operators, with `z=<A>`, `w=<B>` and the defect bound `d` above,

```text
|w-z^2| <= 1-|z|^2+d.
```

The complete Kraus interface in
[`ChoiCPTP.lean`](../../Lean/Dynamics/ChoiCPTP.lean)
provides the compatible finite instrument contract. The attachment is
the explicit completion/calibration/sampling error chain for this concrete
source interface. The spectral estimate and scalar argument here are direct
mathematics checked by executable controls, not new Lean theorems.

Write the stationary bit-channel mean as `b+g <Mtilde_s>`, where
`b=e1-e0`, `g=1-e0-e1>0`, and both unknown flip probabilities are in `[0,1]`.
For fixed, independent binary trials within each setting, Hoeffding's bounded
sum inequality and a union bound give

```text
r_s = sqrt(2 log(12/alpha)/n_s),
Pr(all six sample means are within r_s of their expectations) >= 1-alpha.
```

Independence between settings is not required by the union bound. This uses
[Hoeffding (1963), Theorem 2](https://www.cs.rpi.edu/academics/courses/spring06/random/hoefding.pdf).
The code bounds the logarithm by range reduction and a positive atanh series
with an explicit tail, then bounds its square root using integer arithmetic.
No floating logarithm decides acceptance or rejection.

Let `d_plus,d_minus` be the null-observable errors and let `beta_s` bound the
observed-mean service discrepancy. Set

```text
g_hat = (m_hat_plus-m_hat_minus)/2,
b_hat = (m_hat_plus+m_hat_minus)/2,
t = (r_plus+r_minus+d_plus+d_minus
     +kappa_plus+kappa_minus+beta_plus+beta_minus)/2.
```

On that simultaneous event, both `|g_hat-g|` and `|b_hat-b|` are at most `t`.
Require `g_hat>t`; otherwise the result is inconclusive. Since
`|<M_s>|<=1+c_s`, the corrected mean
`x_hat_s=(m_hat_s-b_hat)/g_hat` satisfies

```text
|x_hat_s-<M_s>| <= [r_s+kappa_s+beta_s+(2+c_s)t]/g_hat.
```

This charges the finite normalization error even though calibration was
inferred from the same count table. Combining each pair of quadrature radii
in Euclidean norm gives `epsilon1,epsilon2`. The exact rejection condition is

```text
|w_hat-z_hat^2|+|z_hat|^2-1
  > d+epsilon2+4|z_hat| epsilon1+2 epsilon1^2.
```

The left side is bounded downward and the right side upward. A rejection
falsifies the **joint declared premises** with family error at most `alpha`;
it does not identify which premise failed. Non-rejection establishes neither
a Gaussian preparation nor an OPH-specific law. The four-setting prefix
uses the same inequalities to bound the future `B` moment and then propagates
the gain/offset and service intervals into the two future observed means.
Those predictions do not depend on the later settings' count values.

## Run a fresh certificate or supply counts

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.calibrated_measurement --check model/calibrated_measurement/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/calibrated_measurement
python -m model.calibrated_measurement --counts /path/to/counts.json
```

The count file contains `counts`, `service_mean_errors` and
`failure_probability`. A count row is `[n00,n01,n10,n11]`, in spin/bit order.
Supply exactly the first four settings for a prospective prediction, or all
six for a comparison. Service entries are exact rational strings or `null`;
all six service obligations must be present even for a prediction prefix.
For example, `failure_probability` is `"1/100"`. A full JSON input example
with unresolved physical service bounds is in [counts_example.json](counts_example.json).

Every invocation recompiles and certifies its own source interface. Zero
service errors in the retained software controls name the complete
mathematical target; they are not a hardware calibration result. The open
physical requirement is acquisition of contemporaneous raw counts and justified
phase, coordinate, stationarity and apparatus bounds. Absolute displacement
scale, physical clocks, the interacting action and the anharmonic mismatch
are separate obligations.
