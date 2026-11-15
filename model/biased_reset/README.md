# Continuing the fixed biased moving-code process

This CP-03 / [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
component extends the existing [moving code](../moving_defect/README.md)
through its antipodal boundary. Every native flag is retained. A supplied
coherent service decodes a known boundary image and re-encodes the unknown
logical qubit at core `(0,1)`; ordinary legal core motion is preserved.
Both before and after the reset the physical carrier is the same binary
two-load source. Classical core metadata is not coherent spatial position.

The source weights are `(25,25,8,...,8)`, with uniform seam attempts and
the complete C/E placement coins. The code family is exactly the 60
encoders declared by `moving_defect`. This is a continuing
controller with explicit services, separate from the stopped receipt.
A core reset is not the number of quantum recovery operations:
active recovery follows every attempted repair.

## Exact optimal asymptotic reset rate

Let `Q` be the complete substochastic kernel on the 60 legal cores, keeping
the six forbidden destinations as distinct boundary outputs. Solve exactly

```text
h(c) = 1 + sum_legal Q(c,d) h(d),  h(forbidden)=0.
```

The positive two-attempt absorption bound establishes finite excursions.
The unique maximum of `h` is at `(0,1)`, with

```text
H = 97909756773858989776401485931436511086245
    /482253783605192403282205691964646745488
  = approximately 203.0253781358 native attempts.
```

This is a count of attempted repairs, including any declared holds; it is
not a physical lifetime. The canonical experiment has repair probability one.
Each attempt lowers conditional expected `h` by one. Every counted reset,
including a preemptive or randomized reset, raises `h` by at most `H`.
For `N_n` resets in `n` completed attempts, telescoping gives

```text
E[h(X_n)] - h(X_0) = -n + E[sum_reset delta h],
E[N_n] >= max(0, (n-h(X_0))/H).
```

Resetting only on reaching a forbidden core, always to `(0,1)`, attains
mean renewal length `H`. Its full 60-core stationary law is independently
checked against the continued kernel, and its reset rate is exactly `1/H`,
approximately `0.00492549261172` per native attempt. Consequently this rate
is globally minimal **within the declared fixed code family and source law**.
The argument permits adaptive C/E choices because their complete core kernels
are identical. It permits arbitrary adaptive/randomized legal-core reset
destinations. It does not permit skipping attempts, halting at the boundary,
discarding outcomes, or substituting another code family without a new map.

For finite horizons, the actual Bellman recursion includes a voluntary reset
before each attempt as well as the compulsory boundary reset afterward:

```text
B_n(c) = sum_live Q(c,d) V_(n-1)(d) + r(c)[1+min_d V_(n-1)(d)],
V_n(c) = min(B_n(c), 1+min_d B_n(d)),  V_0=0.
```

The retained exact calculation covers horizons 0 through 64. Core `(0,1)`
minimizes each of those values, and voluntary resets give no improvement
there. The value from `(0,1)` at horizon64 is approximately `0.1911928840`
resets. These checks are **not** a closed-form proof that the same policy
minimizes every possible finite horizon. The asymptotic proof above does
cover the stated unrestricted horizon limit.

## A complete quantum service, not a reset of the unknown state

For each positive flagged image, its two logical columns have equal size and
disjoint native addresses. The full-register inverse permutation and inverse
uniform mixer decode it; the existing target encoder then produces

```text
R_reset K_flag V_core = sqrt(p_flag) V_(0,1).
```

This identity holds on every logical matrix unit and arbitrary untouched
reference. No branch probability is divided out. The all43,200-slot census
includes 120 positive boundary resets, 11,460 ordinary continuing flags,
17,220 valid zero slots and 14,400 raw abort slots. Abort and zero outputs
remain present; the additional service is identity there. Sixty inverse
Gray schedules are checked on all512 original addresses, including padding.

The fresh floating source replay executes a moving C branch, an E boundary
reset, then another moving C branch. At every declared prefix it retains all
8/16/8 raw outcomes and the unnormalized source/reference amplitudes.
An operator-column test covers all logical inputs, including a complex
reference extension. It does not claim a numerical replay of all43,200
branches or a physical processor-failure certificate. The fixed numerical
admission cap is `1e-8`; it is not an outward hardware-error bound.

The added boundary service costs at most102,326 source pulses. The receipt
also recomputes the original ordinary-recovery/native costs, the complete
stationary mean pulse cost and the largest complete attempt, charging every
invocation even when compiler matrices are cached. The growing archive uses
25 or26 fixed-format bits per C or E attempt with the format in the receipt.
Physical storage, initialization, export and classical-control apparatus
retain their declared service obligations.

## Open requirements

The source law differs from the fair law of
[PR1153](https://github.com/FloatingPragma/observer-patch-holography/pull/1153).
Its fair-control stopping means are150 and120, and a reset to any edge
attains rate1/150. This package does not import that PR's rotation-equivariant
encoder family or its every-code finite-horizon theorem into the biased law.
The [compatibility note](../completion/SAME_PROCESS_INTERFACES.md) records
the exact boundaries.

Continuation does not fix the reference-like exterior obstruction. Nor does
it make the fixed modular flow preserve all moving codes: only19 of60
declared code subspaces are invariant. A reset count, source pulse budget or
active logical qubit is therefore not assigned a topological particle,
physical spin, rest energy or mass. CP-03 and the integrated completion
gates are open.

The optimal reset code itself is invariant under that fixed clock, but one
native attempt can leave the invariant subset. The new stationary receipt
retains the exact positive mass on incompatible codes and the exact
mean coefficient of `sin^2(s log(25/8)/2)` in clock leakage. Optimal resetting
therefore does not secretly repair the moving-clock compatibility premise.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.biased_reset --check model/biased_reset/results.json
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/biased_reset
```
