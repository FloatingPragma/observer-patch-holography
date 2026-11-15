# Check the preparation before testing the dynamics

The native observer predicts six additional, fixed preparation readouts.
They expose a concrete ambiguity in the anharmonic comparison: a non-Gaussian
mixture gives the **same original imaginary characteristic read** as the
transferred thermal preparation, yet gives different additional readouts.
The additional outcomes are not measured here. The 125-point
comparison, its RMSE of 0.287643 and its unresolved normalization are those
of the anharmonic packet.

This extends [action identification](../identification_design/README.md) in
a different direction. That packet identifies a specified family of action
corrections. This one tests the preparation and readout assumptions that must
be settled before a trajectory discrepancy can identify an action. The
parent preparation calibration uses a complex scan. What is missing
is an authenticated, contemporaneous scan for each changed-potential run;
transferring a scan from another run does not establish the state in a later run.

The [protocol](protocol.json) fixes the six probes before their contrasts are
computed. Its design was motivated by the retrospective audit; it is
not a blinded analysis of the deposited dataset. The preparation parameters are
the authenticated parent values, with no refit. The producer never opens
the 125 target trajectory observations; it extracts only the calibration
parameters from the authenticated frozen parent receipt.

## What runs and what is measured

Each probe is a separate fresh preparation, with no evolution or clock:

1. Draw the supplied capped thermal occupation and execute its native Fock
   preparation and real displacement on a 32-level code.
2. Initialize the evolution spin in its supplied zero state. Execute its
   three-pulse native pointer copy and retain both spin outcomes.
3. Execute a second native pointer instrument for the requested
   characteristic coordinate and real or imaginary quadrature.
4. Retain all four `(spin, detector)` probabilities and their quantum output
   branches through the common primitive runtime. No branch is discarded.

The baseline is the supplied displaced-thermal family with
`alpha=-1.0081292326480422` and `nbar=0.13314125172459518`.
The control uses a supplied classical mixture of two source-executed coherent
preparations, with displacements `0` and `-1.5`. The latter receives weight
`0.7581835796978423`. This weight is constructed from the baseline's analytic
prediction to demonstrate blindness of the original read; it is not inferred from
an observed value. The predictor contracts native words; analytic formulas
are independent checks and define this counterexample.

Explicitly, the exact ideal weight is
`w=Im chi_gaussian(0.4i)/Im chi_coherent,-1.5(0.4i)`, since the vacuum
component has zero imaginary characteristic value. This gives exact
blindness in the unrestricted analytic model. The implementation uses its
binary floating representation; its two native original-read laws differ by
total variation about `2.6e-15`, a numerical residual rather than an exact
zero certificate for floating source execution.

For `chi(beta)=Tr[rho D(beta)]`, the two detector outcomes have difference
`Re chi(beta)` or `Im chi(beta)`. The finite-difference position conversion
is not used. Reporting raw four-outcome frequencies in a future experiment
would avoid introducing an unauthenticated position-estimator denominator.
It would require calibration of the physical probe coordinate.

| Fixed probe | Baseline `P(detector=0)` | Non-Gaussian control |
| --- | ---: | ---: |
| Real part, `beta=0.4i` | 0.812678296 | 0.738417930 |
| Imaginary part, `beta=0.4i` (original read) | 0.173836810 | 0.173836810 |
| Real part, `beta=0.8i` | 0.485930446 | 0.384809749 |
| Imaginary part, `beta=0.8i` | 0.166878739 | 0.314060429 |
| Real part, `beta=0.4` | 0.951829772 | 0.961558173 |
| Imaginary part, `beta=0.4` | 0.500000000 | 0.500000000 |

These rounded values are unrestricted analytic reference values; the
[receipt](results.json) separately contains native finite predictions,
complete laws, contraction metadata and independently certified semantic
paths. The spin is initially zero and no dynamics is applied, so the other
spin row is zero in the ideal experiment. Its outcomes remain in every law.

## What the additional probes identify

Within the declared real-displacement thermal family,

\[
 \chi(u+iv)=\exp[-(\bar n+1/2)(u^2+v^2)+2i\alpha v].
\]

If the coordinate `h` is trusted, the two quadratures at `ih` determine

\[
 \alpha=\frac{\arg\chi(ih)}{2h},\qquad
 \bar n=-\frac{\log|\chi(ih)|}{h^2}-\frac12.
\]

This inverse requires nonzero characteristic amplitude and a phase branch
`-pi/(2h) < alpha <= pi/(2h)`. Displacements separated by `pi/h` remain
aliased, even after doubling the coordinate. Near a vanishing characteristic
amplitude, the inverse is poorly conditioned. Its exact inverse is not a
small-shot parameter-estimation guarantee.

The other probes test two independent restrictions of that family:

\[
 \chi(2ih)=|\chi(ih)|^2\chi(ih)^2,
 \qquad \chi(h)=|\chi(ih)|.
\]

The constructed mixture passes the original imaginary probe and fails these
additional restrictions. The real probe at `0.4i` separates the two candidate
binary laws by total variation `0.074260366`; the imaginary probe at `0.8i`
separates them by `0.147181690`. Thus a single initial position estimator
cannot establish which state enters the subsequent evolution. This result
does not assert that the experimental preparation was this particular
mixture.

Under independent ideal trials and the two locked candidates, deciding by
the midpoint of their real-`0.4i` detector probabilities has error at most
`exp(-m*Delta^2/2)` for either candidate, where `Delta=0.074260366`.
The conservative bound falls below 1% at **1,671 fresh trials**. If each
candidate probability has an independently justified uncertainty `epsilon`,
replace `Delta/2` by `max(0, Delta/2-epsilon)` in the usual
`exp(-2*m*margin^2)` bound. No physical apparatus uncertainty or statistical
confidence is inferred from the deposited derived values.
For the stronger, also predeclared imaginary-`0.8i` probe, the certified
contrast exceeds `0.14718168968`, and **426 fresh ideal trials** suffice for
the same conservative bound. Choosing a probe from its calculated separation
is prospective experimental design; it is not validation by fresh data.

## Cutoff and numerical control

An independent evaluator uses the existing
[outward occupation-control construction](../occupation_control/README.md).
It evolves the capped thermal state through the displacement, appends the
pointer, and evolves its opposite half-displacements. Polynomial differential
defects, the occupation edge, the exact thermal tail, endpoint arithmetic
and final readout rounding enter its complete probability-error bound.
Both detector quadratures use the same independently evolved pointer state.
The baseline's complete-law total-variation radii are below `7.06e-13`;
the mixture's are below `3.66e-10`. After these errors, the real-`0.4i`
contrast remains at least `0.07426036559`. These finite errors are far smaller
than the proposed discriminator.

Each certificate applies to its own semantic path and unrestricted supplied
oscillator model. Native source-word floating error is separately checked
against that path; its floating residual is not promoted to a rigorous
native-execution error certificate. The packet does not change the parent
occupation-control receipt or attach its radii to a different trajectory.
All paths here concern preparation and readout at zero evolution duration.
No long-time anharmonic bound is claimed.

## The open probe-scale ambiguity

Suppose the actual coordinate is `s*beta` for unknown positive scale `s`.
The entire Gaussian characteristic scan then depends only on
`A=s*alpha` and `W=s^2*(nbar+1/2)`. The transformation

\[
 s\mapsto c s,\quad \alpha\mapsto\alpha/c,\quad
 \bar n+1/2\mapsto(\bar n+1/2)/c^2
\]

preserves every scan value whenever the transformed thermal occupation stays
nonnegative. More Gaussian probes cannot independently determine all three
quantities. An external displacement calibration or independently established
preparation scale is necessary. The deposited `.636` normalization therefore
cannot be fixed merely by fitting a wider Gaussian scan; a position-conversion
factor and a displacement coordinate are also different apparatus quantities.

The [primary experiment](https://arxiv.org/html/2603.04744v1) describes
characteristic-function readout and two-point position estimation. The
[measurement provenance analysis](../identifiability_checks/README.md#measurement-provenance-and-interpretation) records why its printed
formula and the deposited derived values are not authenticated as
one complete measurement protocol.

The required experiment acquires these contemporaneous probes,
retains raw spin/detector counts and phase/coordinate metadata, and separates
preparation calibration from later evolution tests. Passing this finite menu
would support the stated preparation family; it would not establish an
arbitrary state, select the action or distinguish OPH from ordinary quantum
mechanics. The architecture, preparation family, draws, probe coordinates,
pointer blanks and committed read are supplied services.

## Reproduction

The separate [fixed spectral cross-scale calculation](cross_scale.py) tests
Wissner-Gross's exact entangled-pair scalar through the declared common
$P=\varphi+\sqrt{\pi}\,\alpha$ and $n_s=1-P/48$ maps. Exact rational
square-root and Machin-series bounds are computed before opening the
CODATA reference and the retained ACT DR6 joint chains. The inverse coupling
is $0.5983883$ reported standard uncertainties above CODATA 2022; both tilt
enclosure endpoints retain exactly $3114/508985$ ACT lower-tail mass.
The tilt changes by only $2.4709142\times10^{-14}$ from the inherited
CODATA-calibrated branch. This scalar substitution does not remove the
conditional cosmological tension.

Its [receipt](cross_scale_results.json) also checks the actual scope of the
chiral and oscillator consumers: none of the four retained mass prescriptions
passes both conservative lepton-ratio windows, and those consumers do not
establish one common physical action, clock, readout or complete uncertainty
contract. Arithmetic enclosures are not physical model-error bounds; posterior
positions are not p-values. The [empirical status](../../docs/EMPIRICAL_STATUS.md#one-fixed-spectral-candidate-across-scales)
records the dataset and physical-map boundaries. This calculation replays
the scalar and the ACT weighted projection, while reading the retained
upstream particle/oscillator scope rather than reexecuting their histories.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.empirical_diagnostics --check model/empirical_diagnostics/results.json
python -m model.empirical_diagnostics.cross_scale --check model/empirical_diagnostics/cross_scale_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/empirical_diagnostics
```

Replace `--check` with `--output <path>` to generate a fresh receipt. The
receipt hashes its protocol, preparation input, numerical helpers and native
source dependencies. Source custody, replay and finite-model checks do not
authenticate missing laboratory metadata.
