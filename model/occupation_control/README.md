# Controlled occupation removal for the oscillator apparatus

The [occupation-error derivation](DERIVATION.md) bounds a continuous, unrestricted oscillator's records using an
interval-enclosed finite computation. The bound includes thermal preparation,
preparation displacement, the complete pulse history, the pointer detector,
time discretization and evaluator roundoff.

The result is useful early in the evolution and increasingly loose later.
For case `a` of the supplied [anharmonic apparatus](../anharmonic_transfer/README.md),
the retained 64-level run has complete half-trace error below `1.32e-8` at
0.4 ms and `0.002041` at 2 ms. At 10 ms it is about `0.19946`, too loose to
justify a precision claim. The 32-level bound exceeds one at 10 ms and does
not certify that the selected reference branch even has positive mass. The
32-level conditional estimator certificate is unavailable at both 5 ms and
10 ms. Differences between two cutoffs alone do not establish these guarantees.

| Published coordinate | 32-level half-trace radius | 64-level half-trace radius |
|---|---:|---:|
| 0 ms | 6.75e-13 | 9.51e-13 |
| 0.4 ms | 0.009815 | 1.319e-8 |
| 2 ms | 0.247765 | 0.002041 |
| 5 ms | 0.688188 | 0.026703 |
| 10 ms | 1.674805 | 0.199459 |

These are errors relative to the declared apparatus law, not experimental
error bars. They are conservative integrated residual bounds. A large bound
does not demonstrate a large actual truncation error. The long-window bound
is dominated by occupation removal; making the time integrator more accurate
alone will not fix it. The sparse-cutoff receipt covers case `a`; the complete changed-potential
comparison below certifies the same 64-level evaluator across all five cases
and every deposited endpoint.

For the 64-level run at 2 ms, the selected two-point estimator is
`-0.1434118369` with certified numerical/occupation radius `0.016264164`.
At 0.4 ms the corresponding radius is below `9.905e-8`. These radii include
terminal conditioning and final scalar evaluation; they are wider than the
complete-state half-trace radii because selection and estimator scaling
amplify error.

## What runs and what is measured

The evaluator receives the fixed preparation parameters
`alpha=-1.0081292326480422`, `nbar=0.13314125172459518`, and the supplied pulse,
noise and detector settings. It never reads empirical outcome files. Its
computation is:

1. Prepare a capped thermal mixture, charge its exact tail and rounding, and
   evolve the real-displacement preparation.
2. Retain the spin while evolving the negative/positive/negative sideband
   pulses with number dephasing and the two intervening spin rotations.
   Finite carrier gaps are included for nonzero phase cases.
3. At each requested endpoint, attach a pointer in `|+>`, execute opposite
   half-displacements and read `-Y_pointer`, keeping both final spin branches.
4. Report the spin-selected detector mean and its two-point position
   estimator `Im chi(i*0.4)/(sqrt(2)*0.4)`. This is not exactly `<x>`.
5. Bound every polynomial segment's full differential-equation defect and
   endpoint discontinuity with outward arithmetic, including the single
   occupation edge crossing out of the retained space.

The retained paths use order-12 polynomials, five cells per 50 microsecond
sideband, eight cells for each displacement, and two for each spin rotation.
The integration parameter is supplied drive elapsed duration; displacement
and instantaneous rotation segments use their own normalized control
parameters. Host runtime has no physical-clock interpretation.

**This is an independent semantic evaluator and error certificate.** It does
not execute native observer source words. The apparatus, preparation family,
quantum structure, noise and clock are supplied. Source-derived action
selection lies outside this packet. This packet does not replace native execution
receipts. For orientation, its case-`a` final estimator differs from the
corresponding source-executed finite prediction by about `2.4e-7` at both cutoffs; that numerical
agreement is separate from the error certificate.

## Guarantee and assumptions

The [derivation](DERIVATION.md) maps the existing
[Etienney–Robin–Rouchon theorem](https://arxiv.org/html/2501.09607v6) to this
apparatus and constructs its conservative unrestricted propagator. The
finite polynomial need not be positive. Its full continuous defect, rather
than a collection of sampled boundary populations, controls the comparison.

The interval kernel assumes ordinary IEEE-754 binary64 round-to-nearest
arithmetic, gradual underflow and correctly rounded elementary operations
and square root. It avoids BLAS in the enclosure calculation, widens each
operation, and fails on nonfinite bounds. Target numeric parameters are the
exact binary values in the receipt, not exact symbolic multiples of pi.

The complete error is `half_trace_error_upper`. Final selected records have
separate probability, mean and estimator radii, including scalar evaluation
rounding. Conditional certificates are unavailable when the bound cannot
guarantee positive reference acceptance. The packet additionally reports
`unrestricted_vs_exact_finite_half_trace_error_upper` for composition with
independent finite-source execution bounds.

The certificate excludes experimental parameter uncertainty, detector
normalization provenance, preparation-transfer validity, the laboratory's
unmodeled selection dwell, native compilation error and physical service
error. It establishes no different empirical score or OPH-specific law.

## Scope and reproduction

Use the early-window bound when testing proposed source laws; there is
an explicit numerical uncertainty to keep below the discriminator size.
For longer windows, enlarge the retained occupation space adaptively or use
a tighter boundary-row nuclear-norm enclosure. The current boundary is
computed only at the last occupied level, so adaptive refinement can reuse
the same theorem and sparse solver. Native compilation can be reserved for
the candidate protocols that survive these inexpensive semantic tests.
Larger semantic cutoffs do not by themselves extend the current native
compiler's 64-occupation limit.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/occupation_control
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.occupation_control --check model/occupation_control/results.json
```

To generate a fresh result, replace `--check` with `--output <path>`.
The two retained windows take a few minutes on the development host. Replay
recomputes every bound; it does not infer a certificate from agreement with
the stored rounded numbers. Each receipt's certificate belongs to its own
computed binary path. Cross-platform replay tolerance is a reproducibility
check, not permission to attach a stored radius to different computed means.
Tests include exact rational arithmetic checks,
independent high-precision generator enclosures, actual occupation escape,
a deliberately inaccurate integration step, detector phase conventions and
failure to certify terminal selection.

## Complete changed-potential comparison

The [transfer evaluator](transfer.py) applies the same 64-level apparatus
certificate to all five cases and all 25 positive-time coordinates in the
[published deposit](https://zenodo.org/records/18752287). Its compressed
[receipt](transfer_results.json.gz) retains every full error budget, both
spin branches, all 125 comparisons and five separate zero-time diagnostics.
The calibrated preparation, potential phases, noise, probe `h=.4`, pulse
durations and finite carrier intervals are unchanged. There is no new fit.
Every case is computed and sealed before the scoring code opens the measured
values. This is a retrospective test of supplied apparatus, not prospective
evidence or an OPH-specific prediction. Each coordinate is a terminal read on
an independently prepared trial; storing several computational endpoints
does not impose a sequence of physical measurements on one evolving state.

The conditioning bound uses the approximate output's recorded trace defect.
Write `Delta=rho_exact-rho_approx`, with `||Delta||_1 <= 2 epsilon` and
`|Tr Delta| <= d`. For the accepted detector effect `0 <= A <= Q <= I`,
let `p=Tr(A rho_exact)/Tr(Q rho_exact)`. The spectrum of `A-pQ` lies in
`[-p,1-p]`. Subtracting that interval's midpoint proves

```text
|Tr[(A-pQ) Delta]| <= epsilon + d/2 = e.
```

The same bound for `Q` certifies positive physical acceptance when the
approximate accepted mass `q` exceeds `e`. Then the conditional probability
error is at most `e/q`, the spin-mean error at most `2e/q`, and the position
estimator error at most `2e/(q sqrt(2) h)`, plus the separately enclosed
scalar evaluation errors. This needs Hermiticity, not positivity of the
approximate path. The exact trace defect is enclosed by widening its single
binary rounding. The original conservative radii remain alongside these
tighter radii; neither is inferred from agreement between cutoffs.

For a deposited decimal `y`, certified centre `c` and radius `r`, the
comparison reports the exact rational floor `max(0, |y-c|-r)` on all
remaining sampling, calibration or apparatus/model error. A positive floor
means numerical and occupation error alone cannot cover the discrepancy.
A zero floor does not establish agreement within experimental uncertainty.

All 125 endpoints have certified conditional intervals; **60 measured central
values lie outside them**. The largest required remaining allowance in each
case is shown below, rounded down. No endpoint is excluded.

| Case | Positive discrepancy floors / 25 | Largest floor | Published time |
| --- | ---: | ---: | ---: |
| a | 14 | 0.60059383 | 8.4 ms |
| b | 11 | 0.42598315 | 1.6 ms |
| c | 14 | 0.46868933 | 1.2 ms |
| d | 12 | 0.36578125 | 2.8 ms |
| e | 9 | 0.33039699 | 3.6 ms |

Floors and radii use dimensionless two-point position-estimator units.
Some late-time intervals remain wide. The numerical enclosures do not
establish agreement at the points they contain. Independent checks compare
all 130 endpoint centres, including zero time, with the corresponding
64-level native execution, and test all asymmetric schedules against a
separate dense master-equation solver. This is finite implementation
consistency, not a bound on native compiler or laboratory service error.

The [shared-calibration diagnostic](shared_calibration.py) permits one common
nonnegative output gain `g` and unrestricted offset `b` across all 125 rows.
It minimizes the remaining uniform allowance `e` subject to
`g*lower+b-e <= observed <= g*upper+b+e`, scaling each numerical interval
with the gain. The exact optimum is **0.4916177969622167** estimator units
(rounded here). A feasible rational point and nonnegative dual witness in the
[receipt](shared_calibration_results.json) prove this value against all 252
endpoint and boundary constraints. The dual uses case-a 1.2, 3.2 and 8.4 ms;
every other observation remains checked. The same exact optimum holds when
**each of the five cases has its own nonnegative gain and unrestricted
offset**. The verifier binds all 25 times in each case and checks the same
dual against case a alone. That gives the lower bound even with independent
calibrations; the existing common feasible point gives the matching upper
bound for all five cases. No second optimization or changed prediction is
needed. Requiring shared calibration across different potentials is therefore
not the source of this residual floor.

The witness also constrains heterogeneous nonnegative remaining allowances
`e_1.2`, `e_3.2`, `e_8.4` at those three case-a endpoints. With
`D=10195098079083289586`, its exact necessary inequality is

```text
(4428631686509795952*e_1.2 + 668917353031848841*e_3.2
 + 5097549039541644793*e_8.4)/D
 >= 5012091657452653873877253447/10195098079083289586000000000.
```

The weights sum to one and cancel the gain and offset exactly. They are
approximately `0.4343883357`, `0.0656116643` and `0.5`; the receipt retains
their exact values and endpoint binding. This is a necessary error-budget
condition, not an assignment of experimental uncertainties or a sufficient
test that particular error bars fit the data.

Independent interval choices make
this a conservative relaxation of the common apparatus dynamics. The
optimizing gain and offset are diagnostic and are never adopted. This rules
out a purely affine calibration explanation with zero remaining allowance;
changes to the physical probe coordinate, preparation, timing or dynamics are
outside this family. It supplies no statistical rejection or OPH verdict.

The [paper](https://arxiv.org/html/2603.04744v1#S3.SS2) reports 200 nominal
measurements per point and spin-measurement-noise error bars in Figure 5.
The [deposit](https://zenodo.org/records/18752287) supplies no numerical
per-point uncertainties, accepted counts, raw trials, calibration covariance
or resolution of the documented estimator-normalization ambiguity. No
experimental confidence interval, p-value or detector correction is assigned.

Replay recomputes every path and bound; the compressed format only reduces
artifact size. The sparse-cutoff receipt is a separate certificate.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.occupation_control.transfer --check model/occupation_control/transfer_results.json.gz
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/occupation_control/test_transfer.py
python -m model.occupation_control.shared_calibration --check model/occupation_control/shared_calibration_results.json
```

The last command checks exact arithmetic on the authenticated retained
intervals; the transfer command recomputes their numerical trajectories.
