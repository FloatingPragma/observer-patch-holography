# Two-order characteristic constraints from executable probes

This packet makes an axis-to-interior measurement constraint without a
Gaussian preparation or an oscillator evolution law. It executes both
orders of two noncommuting source-compiled unitary probes, then attaches the
standard quantum Gram-positivity constraint to the published oscillator
characteristic measurements. The laboratory comparison assumes the reported
Weyl displacement instruments and a common normalized preparation. The
randomized future protocol instead permits arbitrary state changes between
trials, provided its setting is drawn independently after preparation.

The result is a **state-independent constraint, a descriptive
consistency diagnostic, and a future randomized test allowing preparation drift**, not a fitted trajectory or an OPH-specific
empirical success. Both unresolved coordinate extents and both operator
orders are retained. The Gaussian and anharmonic comparisons retain their separate scores. The input means have sampling and apparatus errors; they are not
asserted to be exact characteristic values.

## Positivity theorem and instrument attachment

For any normalized positive state functional and unitaries `U,V`, write
`a=<U>`, `b=<V>`. Positivity of the operator Gram matrix gives

```text
|<UV>-ab| <= sqrt[(1-|a|²)(1-|b|²)],
|<VU>-ab| <= sqrt[(1-|a|²)(1-|b|²)].
```

For the first inequality take the three operators `I,U*,V`. Their Gram
matrix is

```text
G = [[1, conjugate(a), b],
     [a, 1, <UV>],
     [conjugate(b), conjugate(<UV>), 1]].
```

Its determinant is `(1-|a|²)(1-|b|²)-|<UV>-ab|²`. Exchange `U,V`
for the second inequality. Keeping the same `a,b` matrix convention and
replacing `<UV>` by `<VU>` also gives a positive matrix: it is the
entrywise conjugate of the Gram matrix for `I,U,V*`. Thus either order
supports the same witness-coordinate convention. This is a small quantum Bochner constraint;
it reuses positivity, not a new axiom about observers. A primary reference
for Weyl relations and quantum Bochner positivity is Lami, Sabapathy and
Winter, [*All phase-space linear bosonic channels are approximately Gaussian
dilatable*, Lemma 2 and Appendix A](https://arxiv.org/pdf/1806.11042).
The finite algebra prerequisite is represented by
[StateExpectation](../../Lean/EventAlgebra/StateExpectation.lean). This
packet supplies a mathematical attachment and executable checks, not a new
Lean formalization of the continuous Weyl theorem.

For the laboratory convention `D(alpha)=exp(alpha a*-conjugate(alpha)a)`,
put `U=D(u)`, `V=D(iv)` for real `u,v`. Then

```text
UV = exp(-iuv) D(u+iv),   VU = exp(+iuv) D(u+iv).
```

Consequently `chi(u+iv)` lies in the **intersection** of disks centered at
`exp(+iuv) chi(u)chi(iv)` and `exp(-iuv) chi(u)chi(iv)`, each with the radius
above. Changing the unresolved coordinate direction swaps the two disks;
it cannot justify selecting the order giving a smaller residual. This
extends the [collinear moment disk](../measurement_bounds/README.md)
to two noncommuting displacement directions. The disks are generally broad;
they do not determine a unique interior mean.

## What the executable observer actually runs

The finite source demonstration chooses `U=exp(i pi Z/7)` and
`V=exp(i pi X/5)` on a declared proper two-dimensional code. The existing
source compiler constructs these words and their controlled versions from
native pulses. A pointer starts at zero, receives a Hadamard, a controlled
word, a phase shift, and another Hadamard. Its two complete measurement
maps are `(I+exp(i phase)W)/2` and `(I-exp(i phase)W)/2` in exact arithmetic.
Zero phase reads the real mean; phase `-pi/2` reads the imaginary mean.

The runtime executes all eight word/quadrature settings, retains both
outcomes, and checks entangled-reference execution and serialized checkpoint
continuation. The matrix products are read from executed source words;
target Pauli matrices are not substituted for those words. Chronological
`V` then `U` implements `UV`, and both orders are run. Raw numerical maps
are neither polarized nor renormalized. Their mass/unitarity defects and
Ramsey residuals are reported separately from the exact theorem.

The chosen code and coherent loading/control/export services are interface
inputs. This two-dimensional demonstration verifies the general unitary
instrument construction. It does **not** implement the exact continuous
Weyl relations: finite matrices cannot represent the oscillator canonical
commutation relation exactly. The laboratory Weyl identification therefore
is an explicit apparatus premise; it is not derived from this qubit
experiment. No physical service-fidelity certificate is available here.

## Public data and results analysis

The data are the authenticated, unchanged zero-duration means from
[McGarry et al.'s deposition](https://zenodo.org/records/18752287), described
in their [measurement protocol](https://arxiv.org/html/2603.04744v1).
The existing loader verifies the deposited hashes. Exact decimal strings
are then read as rational numbers. The procedure uses the 21 axis records
and 100 off-axis complex records of one quadrant; exact mirrored copies
are not additional observations. The measured origin is preserved, not
used to rescale the data. Unit Gram diagonals express normalized-state
mathematics; they do not replace the noisy origin measurement.

For each extent, predictions consume only the axes and coordinate indices;
their output is sealed before the scorer consumes interior responses.
This is a retrospective split by a declared rule. The public values were
available and inspected during development, so this is not a
preregistered or prospective validation. No outcome is fitted and no
coordinate extent is selected by its score.

At extent **4**, one of the 200 order-specific plug-in Gram matrices is
certifiably nonpositive. At extent **5**, three are. In both cases the
strongest retained witness is the `(1,1)` interior point in the positive
product-order convention. A rational witness proves that making its three
complex means jointly physical requires changing at least one by at least
**0.02439** for extent 4, or **0.03129** for extent 5, in complex absolute
value. These are conservative lower bounds on deterministic correction,
not measured detector biases or statistical significance levels.

The certificates evaluate `v*Gv` with rational arithmetic and outward
Taylor intervals for the Weyl phase. If all three complex entries changed
by at most `epsilon`, then `|v*DeltaG v| <= 2 epsilon ||v||²`. A negative
certified quadratic form therefore gives the reported correction floor.
Floating eigenvalues only provide diagnostics. No shot counts are inferred,
no p-value is assigned, and ordinary measurement noise could account for
the incompatibilities. The deposited aggregate means do not provide the
authenticated accepted counts/raw trials needed to decide that question.

## Randomized prediction with memory and preparation drift

A fixed Gram witness can be tested without identical preparations or
independent trials. Write its quadratic form as

```text
q = v*Gv = n + sum_j c_j mu_j >= 0,
n = ||v||², L = sum_j |c_j|,
p_j = |c_j|/L, Y = n + L sign(c_j) X_j.
```

Here `mu_j` denotes one of the six real/imaginary probe means and `X_j` is
its complete binary outcome. Zero-coefficient settings are omitted. Draw
setting `j` with the **frozen** probability `p_j` after the trial's state
has been prepared, independently of that state and its reference memory.
Then `E[Y | past] = q >= 0` for every conditional quantum state. The state
can depend on all previous records and interventions; no reset or stationary
state is required. A preparation success condition may be applied before
setting choice, but readout outcomes cannot be discarded.

For conditional detector mean errors bounded by `epsilon_j`, set
`B = algebra_error + sum_j |c_j| epsilon_j`. Errors must hold uniformly
for the possible pasts and states; marginal calibration alone does not
establish them. The rational Weyl-phase coefficient enclosure contributes
`algebra_error`. The fixed-horizon prediction is

```text
Pr[mean(Y) + B + L sqrt(2 log(1/alpha)/N) < 0] <= alpha.
```

To see why, `Y` lies in an interval of width `2L` and has conditional mean
at least `-B`. Conditional Hoeffding's lemma makes
`exp(-lambda sum(Y+B) - N lambda² L²/2)` a nonnegative
supermartingale. Markov's inequality and optimization over `lambda > 0`
give the displayed bound. This is an attachment of the usual bounded
martingale argument, not a new statistical axiom; see Howard et al.,
[*Time-uniform Chernoff bounds via nonnegative supermartingales*](https://doi.org/10.1214/18-PS321).
The implementation uses a fixed horizon and does not assert an optional
stopping guarantee. Logarithms and square roots use the existing
[exact scalar interval bounds](../calibrated_measurement/README.md).
Missing service bounds produce `unbounded_services`, not a rejection.

The executable source demonstration fixes `v=(-1-i,1,i)`, which gives
`n=4`, coefficients `(-2,2,-2,2,0,-2)`, and five equally likely settings.
It executes the random choice, source Ramsey readout, outcome-dependent
source feedback, and a second random choice/readout on the continuing
carrier. All **100** two-trial records and their final quantum states are
retained, including an entangled reference; checkpoint replay is exact at
the reported floating precision. Conditional second-trial witness means
are positive despite varying with the first record. A control choosing
preparation *after* the setting gives a mean near **-2.833**. This demonstrates
why fresh setting independence is a real premise.

The exact binary source detector arrays are additionally compared with a positive
square `A* A`, where `A = v0 I + v1 U* + v2 V`. Rational matrix
arithmetic and outward square-root bounds prove a uniform score defect at most
**1e-13**, including the correction from each raw detector to its
mathematical trace-preserving completion `K T^(-1/2)`. Each execution recomputes exact rational bounds from its actual binary
arrays and fails if they exceed the fixed rational certificate envelopes.
This permits different compiler/BLAS last bits without comparing exact
platform-dependent fractions as if they were universal values. Per-run
floating diagnostics are not inputs to certification. Execution retains
the raw arrays and reports its roughly **2.22e-15** mass defect; it does not
silently normalize them. The randomizer uses equal binary amplitudes, so
its mathematical trace-preserving completion has exactly the declared
five equal probabilities; its raw mass defect is reported separately.
Feedback can change the state arbitrarily without changing this conditional
score theorem. The score certificate is not a bound on the full executed
output law or its accumulated floating error. It does not certify physical
gate fidelity, randomizer independence, or an exact continuous Weyl representation.
The decision API verifies that its coefficients and phase-error allowance
come from a Gram vector; arbitrary unsigned positivity claims are rejected.

The public means select and freeze one candidate witness per coordinate
extent for a **future** experiment. Under zero physical service error,
**26,103** trials for extent 4 or **14,914** for extent 5 suffice for the
fixed-horizon rejection threshold *if the future mean score is at least as
negative as the deposited plug-in score*. Each test receives `alpha=1/200`.
These are threshold-crossing budgets, not a power calculation or a claim
that those future means will recur. Nonzero service bounds increase the
budgets and may make a crossing impossible. To exclude both unresolved
coordinate scenarios, both relevant hypotheses must be rejected;
rejecting one does not select or reject the other.

No such test is applied to the deposited means: they lack authenticated
ordered random-setting records, complete accepted counts, and the required
conditional service bounds. They are descriptive evidence only.

## What it means and next steps

The preparation-specific Gaussian relation can be checked alongside a
weaker constraint applying to every state measured by the declared Weyl
instruments. The continuing randomized protocol turns positivity into a
finite-sample prediction without assuming a common state across trials,
while exposing setting independence and apparatus accuracy as the actual
remaining experimental obligations. It does not establish that the dataset
is unreliable, select an OPH repair law, or settle the coordinate extent.

The next empirical step is fresh random settings drawn from frozen
probabilities, with all outcomes retained, coordinate calibration and conditional detector-error
bounds. Its witness can be tested despite preparation drift and reference
memory. Public aggregate means cannot substitute for those records. The
independent [closed native-law protocol](../closed_prediction_protocol/README.md)
provides a more discriminating future experiment once a physical mapping
of its services is available. The Gram test is shared by ordinary quantum
models, so agreement would validate the joint instrument/state premises,
not independently confirm OPH or its effective spacetime description.

```sh
python -m model.characteristic_constraints --check model/characteristic_constraints/results.json
python -m pytest -q -W error model/characteristic_constraints
```
