# IBM Quantum Cloud: Frozen Hardware-Fidelity Benchmark Archive

> **Frozen on 11 July 2026.** This program is a hardware-fidelity benchmark
> suite. Prepared-state ratio tests cannot discriminate OPH from standard
> quantum mechanics: QM predicts every prepared state, dynamic circuit,
> measurement record, and feedback operation used here, so each run measures
> how faithfully IBM devices execute the programmed circuits and nothing more.
> The runs are retained as engineering receipts. No further quantum-cloud
> campaign should be presented as OPH evidence until OPH supplies a
> quantitative observable on which its prediction differs from QM.

## Informal story

If OPH is right, then the local quantum structure of reality has a specific
form. Small toy sectors provide controlled places to test implementations of
that structure on IBM quantum hardware. The first campaign checked whether
prepared OPH-motivated states survived hardware execution and readout. A later
blinded campaign tested a stronger, generative question: whether a bounded
patch could read a local boundary record and apply the preregistered repair
move on one development and one held-out backend/layout role.

- Stage 1 asks a structural question: do OPH-like local states look more Markovian and more recoverable than generic or strongly non-Markovian controls?
- `Z_3` asks whether the simplest reduced-sector encoding is internally self-consistent.
- `Z_5` asks whether the measured ratio reproduces a programmed golden-ratio target.
- `S_3` asks the same kind of question in the first nonabelian case, where representation dimensions matter and the target ratio is `2`.
- The blinded generative benchmark asks whether a measured record, rather than
  an initial amplitude alone, selects the correct finite repair against frozen
  open-loop, delayed, shuffled, inverted, state-preparation, noise, label, and
  decoy controls.

The first campaign is state-preparation and readout consistency evidence. The
blinded campaign is implementation evidence for an observer-like self-reading
patch with local state, a measured boundary, a record, conditional feedback,
final readback, append-only receipts, and a public evidence bundle. Neither
campaign supplies evidence for OPH over unrestricted quantum mechanics.

## What these tests actually are

The names `Z_3`, `Z_5`, and `S_3` come from simple symmetry groups.

- `Z_3` means a 3-step cyclic symmetry: imagine a clock with three positions.
- `Z_5` means a 5-step cyclic symmetry: the same idea, but with five positions.
- `S_3` means the six permutations of three objects: it is the smallest genuinely nonabelian case, so order matters.

These are small toy symmetry systems. They are useful for checking that the
OPH-motivated encodings and analysis are internally coherent. Because the
target patterns are directly prepared with standard quantum gates, success is
not evidence that OPH is favored over quantum mechanics.

## What the hardware is doing, informally

The original exact-ratio experiments have the same broad shape:

1. pick a tiny reduced sector where the OPH construction specifies a probability structure;
2. encode that sector into 2 or 3 qubits;
3. prepare the corresponding state on a real IBM device;
4. measure many shots;
5. reconstruct the probabilities from the counts;
6. compare the measured pattern with the programmed OPH-motivated target.

The blinded generative experiment instead prepares balanced finite inputs,
measures a decision record, conditionally applies repair, and retains every
seven-bit joint outcome. Standard quantum mechanics predicts both the prepared
states and this programmed dynamic instrument. The scientific question is
therefore whether the engineered self-reading implementation passes its frozen
reduced-controller nulls, not whether IBM hardware violates or outperforms
standard quantum theory.

## The basic math

For the exact-ratio experiments, the reduced-sector probabilities are modeled in heat-kernel form:

```math
p_R(t) \propto d_R e^{-t \lambda_R}.
```

Here:

- `R` labels the reduced sector or representation;
- `d_R` is its dimension;
- `\lambda_R` is the Laplacian eigenvalue for that sector;
- `t` is the effective heat-kernel time.

The intuition is straightforward:

- sectors with larger `\lambda_R` are more suppressed;
- sectors with larger dimension `d_R` get more multiplicity weight;
- `t` controls how strongly the suppression acts.

So the experiment checks whether the relative measured weights reproduce the
exponential pattern encoded by the circuit. Standard quantum mechanics makes
the same prediction for that preparation.

The measured probabilities are then turned into log-gaps

```math
\Delta_R = -\log(p_R/p_0),
```

or, in the `S_3` nonabelian case,

```math
\Delta_{\mathrm{std}} = -\log\!\left(\frac{p_{\mathrm{std}}/2}{p_0}\right),
```

because the standard representation has dimension `2`.

Why take logs at all? Because an exponential law becomes linear after the log. If

```math
p_R \propto e^{-t\lambda_R},
```

then the corresponding gap behaves like

```math
\Delta_R \sim t \lambda_R.
```

That means ratios of gaps should reproduce ratios of eigenvalues. `Z_5` and `S_3` are useful because they turn the test into a clean number-comparison problem.

The concrete tests are:

### Stage 1: Markov / recoverability

For a 3-qubit local state `\rho_{ABC}`, we reconstruct the density matrix from Pauli tomography and evaluate the conditional mutual information

```math
I(A:C|B)=S_{AB}+S_{BC}-S_B-S_{ABC},
```

along with Petz-recovery quality. OPH favors the direction "lower CMI, better recoverability" for more structured local states.

Informally, conditional mutual information asks: conditioned on the middle subsystem `B`, how much extra correlation remains between the two ends `A` and `C`? Smaller values mean the state is closer to the kind of locally consistent, recoverable structure OPH expects.

### `Z_3`

`Z_3` is the sanity check. The two nontrivial sectors have the same eigenvalue, so extracting `t` from either one should agree:

```math
t_1 = -\frac{1}{3}\log(p_1/p_0), \qquad
t_2 = -\frac{1}{3}\log(p_2/p_0).
```

If `t_1` and `t_2` agree, the reduced-sector preparation and readout path is behaving coherently.

`Z_3` is the pipeline sanity check. It tells us whether the reduced-sector preparation and readout path is internally consistent before the sharper ratio claims are interpreted.

### `Z_5`

`Z_5` is the sharp abelian ratio test. The two nontrivial gap scales satisfy

```math
\frac{\Delta_2}{\Delta_1} = \frac{\lambda_2}{\lambda_1} = \varphi^2 \approx 2.618,
```

where `\varphi` is the golden ratio.

That makes `Z_5` a sharper engineering check than a vague "looks about right"
test. Obvious wrong alternatives such as `2` and `4` are well separated from
the programmed target, although QM predicts that target too.

`Z_5` is useful as an engineering check because it asks for a specific
nontrivial number with little room for interpretive adjustment.

### `S_3`

`S_3` is the first nonabelian exact-ratio test. The gauge-invariant one-plaquette sector is 3-dimensional, with representations

- trivial, dimension `1`
- sign, dimension `1`
- standard, dimension `2`

and OPH-normalized Laplacian eigenvalues

```math
\lambda_{\mathrm{triv}} = 0,\qquad
\lambda_{\mathrm{sign}} = 6,\qquad
\lambda_{\mathrm{std}} = 3.
```

So the parameter-free prediction is

```math
\frac{\Delta_{\mathrm{sign}}}{\Delta_{\mathrm{std}}} = 2.
```

`S_3` is the first nonabelian case in this ladder. It is the first test where the representation dimension plays an essential role in the observable.

## Frozen public code and data

The runnable code, representative legacy derived summaries, and complete
blinded-run evidence are stored in:

- `code/ibm_quantum_cloud/programs/`
- `code/ibm_quantum_cloud/qc_data/`
- `code/ibm_quantum_cloud/runs/issue_509_20260711_v2/`

The public data bundle includes the main Stage 1, `Z_3`, `Z_5`, and `S_3`
hardware outputs, the `S_3` layout diagnostic, and the complete blinded
generative run: manifest, lock, reveal, provider job IDs, raw joined counts,
calibration receipts, packet, blind report, hashes, and replay commands.

## Frozen archived results

In the tables below, `Prepared t` is the operating point used to prepare the state. The main tested observable depends on the experiment.

- In `Z_3`, the tested observable is the extracted heat-kernel time, so the comparison is prepared `t` versus measured `t`.
- In `Z_5`, the tested observable is the ratio `Δ₂/Δ₁`, so the comparison is expected `φ²` versus measured `Δ₂/Δ₁`.
- In `S_3`, the tested observable is the ratio `Δ_sign/Δ_std`, so the comparison is expected `2` versus measured ratio.

### Stage 1: programmed ordering reproduced on two backends

Hardware runs:

- `ibm_marrakesh`, job `d6t4da6sh9gc73di7720`
- `ibm_fez`, job `d6t4ejngtkcc73cm8l6g`

Programmed target and reason:

OPH favors locally structured states with lower conditional mutual information and better recoverability. For this benchmark, that gives three direct expectations:

- `CMI(structured_theta_0.00) < CMI(random_seed_2)`
- `CMI(structured_theta_0.00) < CMI(ghz_control)`
- `Petz fidelity(structured_theta_0.00) > Petz fidelity(random_seed_2) > Petz fidelity(ghz_control)`

Measured comparison:

| Programmed target | Expected relation | `ibm_marrakesh` | `ibm_fez` | Result |
| --- | --- | --- | --- | --- |
| structured state has lower CMI than random control | true | `0.2309 < 0.3890` | `0.1498 < 0.4992` | pass |
| structured state has lower CMI than GHZ control | true | `0.2309 < 0.9474` | `0.1498 < 0.9166` | pass |
| recoverability ranks structured > random > GHZ | true | `0.9297 > 0.8649 > 0.6066` | `0.9479 > 0.8117 > 0.6449` | pass |

Representative reconstructed values:

| Backend | State | CMI (bits) | Petz fidelity |
| --- | --- | ---: | ---: |
| `ibm_marrakesh` | `structured_theta_0.00` | `0.2309` | `0.9297` |
| `ibm_marrakesh` | `random_seed_2` | `0.3890` | `0.8649` |
| `ibm_marrakesh` | `ghz_control` | `0.9474` | `0.6066` |
| `ibm_fez` | `structured_theta_0.00` | `0.1498` | `0.9479` |
| `ibm_fez` | `random_seed_2` | `0.4992` | `0.8117` |
| `ibm_fez` | `ghz_control` | `0.9166` | `0.6449` |

On both backends, the most structured state reconstructs below both controls in CMI, and recoverability worsens as structure is reduced. That is the OPH-favored direction.

The Stage 1 export also preserves the fixed random-control choice recorded in the run catalog: `random_seed_2` was selected from a small candidate seed set by exact CMI before the hardware run.

### `Z_3`: internal-consistency pass

Hardware run:

- `ibm_marrakesh`, job `d6t4ic790okc73et0n4g`

Programmed target and reason:

In `Z_3`, the two nontrivial sectors have the same eigenvalue `lambda = 3`, so both channels should return the same heat-kernel time. The expected measured value is therefore the prepared `t` itself.

| Prepared `t` | Expected measured `t` | Mean measured `t` | Abs. error | Internal abs. `t_q1 - t_q2` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0.30` | `0.3000` | `0.3162` | `0.0162` | `0.0142` | `0.10%` |
| `0.60` | `0.6000` | `0.6147` | `0.0147` | `0.0034` | `0.07%` |
| `0.90` | `0.9000` | `0.9153` | `0.0153` | `0.0043` | `0.05%` |

The extracted `t` value stays within about `0.02` of the prepared `t` at all three operating points, and the two independent extractions agree closely. This is a clean internal-consistency success.

### `Z_5`: prepared-state readings around the programmed ratio

Hardware runs:

- `ibm_marrakesh`, job `d6t4iiv90okc73et0nbg`
- `ibm_marrakesh`, job `d6t4rkjbjfas73fp65m0`
- `ibm_fez`, job `d6t4ipvgtkcc73cm8qg0`
- focused high-shot `ibm_marrakesh` rerun, job `d6t50e790okc73et186g`

Programmed target and reason:

In `Z_5`, the gap ratio is fixed by the programmed Laplacian spectrum:

`Delta2 / Delta1 = lambda2 / lambda1 = phi^2 ~= 2.6180`

So the expected value at every prepared `t` is the same number, `2.6180`.

Main `ibm_marrakesh` sweep:

| Prepared `t` | Expected `Delta2 / Delta1` | Measured `Delta2 / Delta1` | Abs. ratio error | Relative error | Internal abs. `t1 - t2` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0.30` | `2.6180` | `2.5974` | `0.0207` | `0.8%` | `0.0025` | `0.29%` |
| `0.60` | `2.6180` | `2.7383` | `0.1203` | `4.6%` | `0.0268` | `0.39%` |
| `0.90` | `2.6180` | `2.5786` | `0.0395` | `1.5%` | `0.0137` | `0.24%` |

Second `ibm_marrakesh` sweep:

| Prepared `t` | Expected `Delta2 / Delta1` | Measured `Delta2 / Delta1` | Abs. ratio error | Relative error | Internal abs. `t1 - t2` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0.30` | `2.6180` | `2.4532` | `0.1649` | `6.3%` | `0.0199` | `0.22%` |
| `0.60` | `2.6180` | `2.8778` | `0.2597` | `9.9%` | `0.0556` | `0.22%` |
| `0.90` | `2.6180` | `2.6488` | `0.0308` | `1.2%` | `0.0104` | `0.39%` |

`ibm_fez` replication:

| Prepared `t` | Expected `Delta2 / Delta1` | Measured `Delta2 / Delta1` | Abs. ratio error | Relative error | Internal abs. `t1 - t2` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0.30` | `2.6180` | `2.7262` | `0.1082` | `4.1%` | `0.0127` | `1.03%` |
| `0.60` | `2.6180` | `2.4068` | `0.2113` | `8.1%` | `0.0521` | `0.66%` |
| `0.90` | `2.6180` | `2.4610` | `0.1570` | `6.0%` | `0.0537` | `0.54%` |

Focused high-shot follow-up at `t = 0.90`:

| Backend | Expected `Delta2 / Delta1` | Measured `Delta2 / Delta1` | Abs. ratio error | Relative error | Internal abs. `t1 - t2` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ibm_marrakesh` | `2.6180` | `2.5498` | `0.0682` | `2.6%` | `0.0238` | `0.31%` |

One verdict governs this lane: the high-precision run excludes the programmed
target. The focused high-shot rerun measured `2.5498` with interval
`[2.4919, 2.6023]`, and `phi^2 = 2.6180` lies outside that interval. The
lower-shot sweeps scatter on both sides of the target and carry no verdict of
their own. In every case, QM predicts the prepared circuit rather than a
competing value, so the exclusion is a device-fidelity statement, not an
OPH-versus-QM statement.

### `S_3`: preregistered test failed; layout-sensitive post-hoc follow-up

Hardware diagnostic bundle:

- base layout `ibm_marrakesh`, job `d6t5c47gtkcc73cm9pd0`
- reversed layout `ibm_marrakesh`, job `d6t5cdvgtkcc73cm9pr0`
- direct reversed-layout repeat, job `d6t5csmsh9gc73di8e8g`

Programmed target and reason:

In `S_3`, the programmed one-plaquette reduced sector has `lambda_sign = 6`
and `lambda_std = 3`, so the target is:

`Delta_sign / Delta_std = 2.0000`

At `t = 0.60`:

| Configuration | Expected ratio | Measured ratio | Abs. ratio error | Relative error | `Δt/t` | Leakage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| base layout `[7, 8]`, raw | `2.0000` | `1.8536` | `0.1464` | `7.3%` | `-0.0760` | `0.08%` |
| base layout `[7, 8]`, mitigated | `2.0000` | `1.8724` | `0.1276` | `6.4%` | `-0.0659` | `0.06%` |
| reversed layout `[8, 7]`, raw | `2.0000` | `2.0348` | `0.0348` | `1.7%` | `0.0172` | `0.04%` |
| reversed layout `[8, 7]`, mitigated | `2.0000` | `2.0299` | `0.0299` | `1.5%` | `0.0148` | `0.01%` |

Direct repeat of the selected mapping:

- `ibm_marrakesh`, seed `17`, job `d6t5csmsh9gc73di8e8g`
- expected ratio `2.0000`
- measured ratio `2.0657`
- abs. ratio error `0.0657`
- relative error `3.3%`
- `Δt/t = 0.0323`
- leakage `0.01%`

The preregistered test is recorded as failed: the base-layout runs measured
`1.8536` (raw) and `1.8724` (mitigated) against the frozen target `2.0000`.
That verdict is permanent under the failure-permanence rule.

The physical mapping clearly changes the extracted value. Because the reversed
mapping was selected after the first result was seen, its near-target repeat is
a useful device diagnostic, not a blinded confirmation of the ratio. The
reversed-layout observation is a new hypothesis, not a rescue of the failed
test; it requires a fresh preregistered run across multiple layouts and
backends, with confidence intervals, before it carries any weight.

### Frozen blinded generative repair: preregistered engineering pass

The production-random `issue509-20260711-v2` run used `ibm_fez` for the
development role and submitted the held-out `ibm_kingston` role only after the
development jobs had completed. The run contained 3,104 circuits, 606,208
shots, and 14 completed IBM Runtime jobs. No shots were dropped or
postselected.

The frozen blind report passed every calibration, leakage, completeness, and
shot-conservation gate. Record-gated repair exceeded the preregistered
conditional-likelihood threshold against every required null on each backend
and pooled. The weakest comparison was the 256-component label/layout mixture,
with likelihood ratio approximately `256`; all other primary comparisons had
per-backend log-likelihood ratios above `77,000`. The committed reference
label/layout was the unique top component among 256 possibilities.

After the blind report was file-hashed, the reveal commitment and all 3,104
normalized-OpenQASM logical-circuit hashes verified. The raw receipts and exact
replay are public in
`code/ibm_quantum_cloud/runs/issue_509_20260711_v2/`.

## Frozen assessment and interpretation boundary

The engineering outcomes are real:

- the legacy runs show how well IBM devices prepared and read the selected
  finite-sector states, including visible layout sensitivity;
- Stage 1's offline reconstruction preserved its intended ordering on two
  backends; and
- the blinded run strongly identifies the programmed record-gated controller
  against its restricted null family and provides unusually complete receipts.

Their evidentiary weight for **OPH versus standard quantum mechanics is zero**.
The state experiments directly encode the target amplitudes. The dynamic
experiment uses standard measurement, a classical record, and standard
conditional feedback. Its likelihood ratios are controller-versus-controller
comparisons, not OPH-versus-QM comparisons, and they are not Bayes factors for
OPH.

More shots, qubits, IBM backends, or blinding would improve engineering
replication without repairing that missing theory contrast. The archive is
therefore frozen until OPH produces a source-closed observable with a numerical
QM disagreement, predicted effect size, nuisance model, audit against existing
bounds, and preregistered decision rule. The detailed reactivation gate is in
[`../code/ibm_quantum_cloud/README.md`](../code/ibm_quantum_cloud/README.md).
