# Interaction recovery with the complete field participating

The observer identifies the local quartic action while **every onsite and
bond term participates in every calibration trial**. This removes the
isolated-generator access used by the [preceding reconstruction](../interacting_reconstruction/README.md).
The inverse receives only complete binary detector laws, opaque labels,
protocol declarations and separately executed clock probabilities.

This is a finite source-executed identification benchmark. The action family,
number preparation, pair detectors, schedule and clock calibration are
supplied. The inference is Python arithmetic on accessible records, not
an arithmetic program compiled into the native observer. No empirical
comparison or interacting continuum construction is claimed.

## Executed computation

The backend uses eight four-level sites, so the calibration field has 65,536
amplitudes. It keeps the same local canonical scalar action and hidden channel
permutation as the preceding packet. For each of two separately read clock
intervals, near `1/16384` and `1/8192`, the source executes:

1. A product number preparation: either all sites in their oscillator blanks,
   or one addressed site in `|1>` and all other sites blank. These are supplied
   oscillator preparations, not the interacting vacuum.
2. One complete symmetric field step: eight onsite half steps, all twelve
   bond steps, then the eight onsite half steps in reverse order. Each local
   matrix is contracted from actual emitted twelve-port source words.
3. A binary source pointer on the addressed site, or on an addressed pair.
   The single-site effects read `(|0> ± i|2>)/sqrt(2)` and
   `(|1> ± i|3>)/sqrt(2)`; pair effects read `(|00> ± i|11>)/sqrt(2)`.
   Each phase uses a fresh trial. Every spectator is traced over, and the
   pointer's entire complement is retained without postselection.

All 28 opaque pairs are probed, including those eventually inferred to have
no direct edge. Pair instruments operate on the full 16-level pair register
plus a binary pointer. Their accessibility is supplied; the construction
does not show how a local apparatus in the recovered metric reaches arbitrary
pairs. The model-coordinate interval scales the same ordered gate schedule
at both durations. It is not a host-runtime measurement.

There are 88 binary trial laws per interval, or 176 altogether. The host
contracts nine independent preparation columns together for efficiency;
these columns are not coherently superposed and do not reduce the declared
fresh-trial resource count. Calibration requests 4,231,808 native pulses in
total. These are computed full probability laws, not physical shot samples.

## Why the short calculation suffices

For a local imaginary-coherence readout, the difference of the two phase
probabilities is an odd function of the interval. Its first derivative
identifies the relevant generator matrix element. The blank spectators have
zero position mean, so external bonds disappear from the first onsite
derivative. For the pair probe, only the direct pair coupling contributes at
first order. Higher-order interactions and indirect paths remain in the
actual probabilities.

The two-duration Richardson combination cancels the cubic term. It can be
applied **directly to the complete product schedule**, even though one
product step is not exactly `exp(-itH)`. A declared upper bound of 96 on
the sum of centered factor-generator norms controls the fifth derivative.
The final backend audit verifies that this sum is below 87. Consequently
calibration needs no separately amplified splitting-error allowance.

The factor sequence must stay fixed across durations. The runtime rejects
mismatched step counts or missing whole-field declarations. A counterexample
in the tests shows that mixing step counts invalidates the cancellation.
A bound on the norm of the summed Hamiltonian alone would also be insufficient;
the [derivation](DERIVATION.md) supplies the product-rule proof and counterexample.

Once identified, the coefficients define the same local quartic family at
the two withheld occupation cutoffs. The independent predictor computes
complete two-read histories before their corresponding source executions:
coupled evolution, first measurement, first-bit-controlled feedback,
continued evolution and second measurement. Both native pointer phases and
all 1,056 joint native/quantum outcomes are retained at each cutoff.
These are the parent's predeclared measurement protocols, withheld from this
calibration; they are not prospective experimental data.

## Retained result

The [receipt](results.json) recovers all twelve edges and the rank-three cube
Gram, conditional on the same unit-edge convention and topology-separation
premise as before. The recovered coefficients are approximately:

```text
bare mass squared = 1.00000000027
bond weight       = 0.99999999997
quartic g         = 0.04166666664   (target 1/24)
```

The bare mass and bond/speed coefficients describe the supplied action;
they are not measurements of the interacting excitation gap or relativistic
propagation speed. A positive finite-duration interval resolves `g>0` within
the declared canonical polynomial family.

| Complete two-read comparison | 8 × 2⁷ levels | 8 × 4⁷ levels |
| --- | ---: | ---: |
| Numerical prediction/execution TV | <1.2e-13 | <6.8e-13 |
| Ideal coefficient + evolution-splitting TV bound | <2.97e-6 | <2.14e-5 |
| Matched quartic-deletion discrepancy | >3.90e-5 | >3.90e-5 |

Numerical agreement is not a roundoff certificate. The ideal budgets exclude
clock uncertainty, floating errors and physical service errors. Existing
occupation bounds transfer under the declared cube relabeling; they are
larger than this particular interaction-deletion contrast, so the contrast
is not certified after occupation removal.

The calibration is observably different from isolated probing. As a separate
test at model interval `1/8`, the complete-field onsite probabilities differ
from the isolated ones by more than `1e-4` (about `0.0053` in the retained
run); some nonedge pair phase asymmetries exceed `1e-5` (about `0.001`). The
receipt does not store these two values; the tests assert the thresholds. These reflect higher-order coupled effects. Independent
full-Hamiltonian comparisons, spectator-trace tests and source-corruption
tests check that the executable does not substitute isolated laws.

## Remaining work

The improvement is removal of generator isolation. The pair instruments,
number preparations, fixed action family, separated-zero topology assumption
and independent calibrated clock are premises. Short-time inference
amplifies probability noise; one evolution step does not solve physical
sampling cost. No finite-shot experiment was run.

The next useful extension is an economical observation protocol with fewer
trusted controls and a clock included in the continuing process, followed by
interacting observables that retain their meaning over a common spatial,
occupation and energy refinement domain. Action selection by native agreement,
native/field geometry identification, interacting relativistic covariance
and a controlled continuum theory are open.

## Reproduce

From the repository root with `model/requirements.txt` installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.coupled_reconstruction --check model/coupled_reconstruction/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/coupled_reconstruction
python -m model.verify --hashes-only
```

The [manifest](manifest.json) binds all package artifacts. The result binds
the imported producers, including the frozen parent's prediction and
occupation-certificate code. Existing evidence packets are unchanged.
Replay allows `2e-8` floating variation in inferred coefficient diagnostics,
because the short-time inverse amplifies rounding; raw laws, clocks and
predictions retain `2e-10`. Structures, counts, hashes and fixed proof strings
are compared exactly. These numerical replay tolerances are not physical error bounds.
