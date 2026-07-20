# Frozen IBM Quantum Cloud Engineering Bundle

> **Status: frozen on 11 July 2026.** No experiment in this bundle
> distinguishes OPH from standard quantum mechanics. The code, derived legacy
> summaries, and complete issue `#509` receipts are retained as an engineering
> and reproducibility archive. New OPH-evidence runs on quantum-cloud hardware
> are paused until OPH supplies a quantitative prediction that differs from
> standard quantum mechanics.

This directory is the public code-and-data archive for the completed IBM
Quantum Cloud OPH-motivated engineering experiments.

It contains:

- `programs/`: runnable preflight, sealing, execution, harvesting, and analysis scripts
- `qc_data/`: representative derived JSON summaries from the legacy hardware runs
- `protocols/`: frozen or pre-hardware experimental protocol specifications
- `tests/`: mandatory finite-kernel tests and optional Qiskit/Aer circuit tests

## Scope

The bundle contains two scientifically different kinds of experiment.

1. The original Stage 1, `Z_3`, `Z_5`, and `S_3` runs prepare states whose
   mathematical pattern is already specified, then ask whether a real QPU can
   preserve and read that pattern. These are state-preparation, readout, and
   layout-consistency tests.
2. The issue `#509` benchmark starts from ordinary basis states and generates
   its output through a mid-circuit record and conditional feedback. It was
   preregistered, blinded, run on a development QPU and a held-out QPU, and
   analyzed before the semantic reveal was opened. Its complete receipts are
   under [`runs/issue_509_20260711_v2/`](runs/issue_509_20260711_v2/).

The second kind is stronger evidence for the engineered controller claim
because the bounded patch reads its own local state, writes a record, uses that
record to choose a repair move, and publishes the evidence needed to audit the
process. It is not evidence for OPH over unrestricted quantum mechanics.
Standard quantum mechanics predicts every programmed circuit in this bundle.

The archived legacy summaries describe 10 provider jobs, 294 circuit
instances, and 400,384 shots. The blinded benchmark adds 14 jobs, 3,104
circuits, and 606,208 shots. Those totals do not imply equal evidential weight:
the legacy `qc_data/` files mostly preserve derived metrics rather than every
provider count, whereas the blinded run preserves the complete raw joined
counts, job journals, commitments, and replayable analysis.

## How evidence strength is graded here

“Strong evidence” is meaningful only after stating what the evidence is for.
This README uses four practical levels.

- **Sanity check:** verifies that an encoding, analysis, or readout behaves as
  expected. It is useful for finding mistakes but weak as a theory test.
- **Supportive hardware consistency:** repeats a predicted pattern on real
  hardware, sometimes on more than one chip, but the target was prepared in
  the state or an important choice was made after looking at data.
- **Confirmatory implementation test:** freezes the target, controls, layouts,
  analysis, and failure rule before a blinded held-out run. This can strongly
  test the stated device or controller claim.
- **Theory-discriminating test:** competing physical theories make different
  predictions for the same intervention. None of the experiments in this
  directory reaches this level for OPH versus standard quantum mechanics.

The immediate engineering claim and the broad OPH claim therefore receive
different assessments:

| Experiment | Immediate outcome | Strength for the immediate claim | Strength for OPH as a physical theory |
| --- | --- | --- | --- |
| Stage 1 recoverability | Structured states retained lower conditional dependence and higher recovery fidelity than random and GHZ controls on two QPUs | Supportive hardware consistency | None as a discriminator; the ordering is designed into QM-prepared states |
| `Z_3` heat-kernel check | Two independently decoded nontrivial states gave nearly the same heat time at all three prepared points | Strong sanity check | None as a discriminator; QM predicts the directly prepared distribution |
| `Z_5` golden-ratio check | Several runs landed near `phi^2`, but the values spread from `2.407` to `2.878`, and the high-shot point was below the exact target | Mixed supportive consistency | None as a discriminator; QM predicts the encoded amplitude ratio |
| `S_3` nonabelian check | The base layout missed `2`; reversing the physical qubits moved the result near `2`, and a repeat stayed near `2` | Useful layout diagnostic | None as a discriminator; QM predicts the prepared state and measurement |
| Blinded record-gated repair | Passed every frozen primary null on both the development and held-out QPUs | Strong confirmatory evidence for the programmed self-reading repair implementation | None as a discriminator; QM predicts measurement, records, and feedback |
| Secondary dimension-exponent test | Identifiable nonabelian spectra favored `kappa = 1`; abelian negative controls correctly had zero separation | Supportive diagnostic evidence | None as a discriminator; it checks a programmed stationary-law convention |

## Experiments and outcomes in detail

### 1. Stage 1: Markov structure and recoverability

#### What it asks

Imagine three neighboring pieces of a system, called A, B, and C. If B carries
almost all of the information needed to connect A with C, then the full state
is close to a local Markov chain. The benchmark measures this with:

- **conditional mutual information (CMI):** lower is more locally Markovian;
- **Petz recovery fidelity:** closer to `1` means a standard recovery map,
  computed from the reconstructed state, rebuilds the full three-part state
  more accurately from overlapping local pieces.

The circuit family contains three deliberately structured states, one
fixed-depth random control, and one GHZ control. Full three-qubit Pauli
tomography reconstructs each density matrix. Each tomography circuit uses 512
shots. The Petz recovery is an offline calculation on that reconstruction; the
QPU does not execute a Petz-recovery circuit. The random control, seed 2, was
chosen from ten candidate random circuits because it had the largest ideal
CMI. The test was run on `ibm_marrakesh` and independently on `ibm_fez`.

#### What happened

The simplest structured state, `structured_theta_0.00`, had the lowest measured
CMI and the highest measured Petz fidelity on both chips:

| Backend | CMI: structured | CMI: random | CMI: GHZ | Recovery fidelity: structured | Random | GHZ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ibm_marrakesh` | `0.231` | `0.389` | `0.947` | `0.930` | `0.865` | `0.607` |
| `ibm_fez` | `0.150` | `0.499` | `0.917` | `0.948` | `0.812` | `0.645` |

All three preregistered qualitative fingerprint checks returned `true` on both
backends. The provider jobs are `d6t4da6sh9gc73di7720` and
`d6t4ejngtkcc73cm8l6g`; the public derived reconstruction summaries are in
[`qc_data/stage1/`](qc_data/stage1/).

#### How strong is it?

This is useful replicated evidence that IBM hardware and the tomography
pipeline preserve the intended recoverability ordering. It is not evidence for
OPH over quantum mechanics, which predicts the prepared states. The structured
circuits were constructed to have lower ideal CMI, generic quantum information
theory already connects low CMI with recovery, and the random control was one
selected member of a candidate family rather than a held-out control. The
public legacy files do not provide the complete tomography count tables, and
the reported hardware metrics have no formal uncertainty test. Projecting a
noisy reconstruction onto a physical density matrix can also bias nonlinear
quantities such as CMI and fidelity. The result is best called **supportive
hardware consistency**, not a blinded confirmation.

### 2. `Z_3`: the three-state internal-consistency check

#### What it asks

`Z_3` is a three-position cycle. The experiment prepares a heat-kernel-shaped
probability distribution at three declared heat times, `t = 0.3`, `0.6`, and
`0.9`. Its two nontrivial states should encode the same heat time. Reading each
one independently therefore gives a simple cross-check: if the two inferred
times disagree badly, the encoding or readout is not trustworthy.

The run used 4,096 shots per prepared time on `ibm_marrakesh`, with one unused
two-qubit code retained as leakage rather than discarded.

#### What happened

| Prepared `t` | Inferred from state 1 | Inferred from state 2 | Absolute difference | Leakage |
| ---: | ---: | ---: | ---: | ---: |
| `0.3` | `0.3233` | `0.3091` | `0.0142` | `0.098%` |
| `0.6` | `0.6164` | `0.6130` | `0.0034` | `0.073%` |
| `0.9` | `0.9174` | `0.9131` | `0.0043` | `0.049%` |

The two decodings agree closely at all three points. The provider job is
`d6t4ic790okc73et0n4g`, with derived probabilities and bootstrap intervals in
[`qc_data/z3/`](qc_data/z3/).

#### How strong is it?

This is a **strong sanity check** for the small encoding and analysis path. It
is not evidence for OPH over QM: the heat-kernel probabilities were placed in
the initial state, both inferred times come from that same prepared state,
and the check used one backend. The analysis renormalizes the probabilities
inside the three valid codes after separately reporting the very small invalid-
code leakage. That is reasonable for checking the encoded sector, but it is
another reason not to read this as an independently generated law of nature.

### 3. `Z_5`: the golden-ratio state-preparation check

#### What it asks

`Z_5` is a five-position cycle encoded in three qubits. For the chosen spectrum,
two logarithmic probability gaps are supposed to have the ratio
`phi^2 = 2.6180339887...`. The experiment prepares the corresponding amplitudes,
measures them, and reconstructs that ratio. A sixth, seventh, or eighth code is
counted as leakage.

#### What happened

The target neighborhood was reproduced several times, but not with the
stability expected of a new exact-number measurement:

| Run | Shots per point | Ratio at `t=0.3` | `t=0.6` | `t=0.9` |
| --- | ---: | ---: | ---: | ---: |
| `ibm_marrakesh`, seed 11 | 4,096 | `2.597` | `2.738` | `2.579` |
| `ibm_marrakesh`, seed 17 | 4,096 | `2.453` | `2.878` | `2.649` |
| `ibm_fez`, seed 17 | 4,096 | `2.726` | `2.407` | `2.461` |
| `ibm_marrakesh`, high-shot `t=0.9` | 32,768 | n/a | n/a | `2.550` |

The observed ratios range from `2.407` to `2.878` around the target `2.618`.
The high-shot `t=0.9` bootstrap interval was approximately `[2.492, 2.602]`,
which sits slightly below the exact target. Leakage remained between about
`0.22%` and `1.03%`. Seven of the nine ordinary 4,096-shot intervals included
the target; the separate high-shot interval did not. Provider job IDs and the
full per-point intervals are in [`qc_data/z5/`](qc_data/z5/).

#### How strong is it?

The repeated near-target values are **supportive hardware consistency** for
preparing and reading the intended five-state family. The spread, including a
high-shot interval that excludes the exact value, makes the evidence mixed for
an exact hardware ratio. More importantly, the target is encoded directly in
the prepared amplitudes. The bootstrap intervals quantify finite-shot
variation, not calibration drift, model choice, or other systematic error.
This is not an independent derivation of `phi^2` by the QPU and provides no
OPH-versus-QM evidence.

### 4. `S_3`: the nonabelian ratio and layout diagnostic

#### What it asks

`S_3`, the six permutations of three objects, is the smallest nonabelian group.
Its three reduced representation sectors are encoded in two qubits. After
accounting for the two-dimensional standard representation, the two measured
logarithmic gaps should have ratio `2`. This is again a prepared-state test,
but it also probes whether assigning the two logical sector bits to different
physical qubits biases the answer.

#### What happened

At prepared `t = 0.6`, 16,384-shot jobs on `ibm_marrakesh` gave:

| Physical layout | Raw ratio | Readout-mitigated ratio | 95% bootstrap interval after mitigation |
| --- | ---: | ---: | ---: |
| base `[7, 8]` | `1.854` | `1.872` | `[1.807, 1.928]` |
| reversed `[8, 7]` | `2.035` | `2.030` | `[1.966, 2.100]` |

The base layout missed the target. Reversing the same two physical qubits moved
the ratio across the target neighborhood. A separate seed-17 confirmation on
the reversed layout measured `2.066`, with interval `[1.995, 2.127]` and
`0.012%` leakage. The jobs and assignment matrices are in
[`qc_data/s3/`](qc_data/s3/).

#### How strong is it?

This is a valuable **hardware-bias diagnostic**: it shows that the sector-to-
qubit map materially changes the extracted ratio, and the archived summaries
make that dependence visible. It is not an OPH-versus-QM test. The
reversed layout was investigated after the base-layout bias was observed, so
landing near `2` is partly a post-observation layout result. The later repeat
supports the diagnosis but does not retroactively blind the layout choice. All three
jobs used the same IBM backend, the same physical-qubit pair, and essentially
the same operating point. The bootstrap also treats the fitted assignment
matrix as fixed, so its uncertainty is not propagated into the interval.

### 5. Blinded generative record-gated repair

#### What it asks

The older tests ask whether a target already placed in an initial state survives
the hardware. This benchmark asks whether the device can generate a repair
decision from a record created during the circuit. “Generative” is used in
that control-system sense: the output is generated by measurement and
feedback. It does not refer to generative AI.

Each dynamic circuit uses three qubits for a finite mismatch state and one
decision qubit. It:

1. starts from an ordinary computational-basis state;
2. applies a balanced lazy disturbance;
3. reads the disturbed state into a mid-circuit record;
4. decides whether a proposed move reduces the mismatch;
5. records that decision;
6. applies the repair only when the decision says it should; and
7. records the final state.

That is the observer-like technology claim in concrete form: a bounded patch
has local state and a boundary, reads itself, writes a record, feeds the record
back into a repair move, and emits a public audit trail.

The primary repair arm is compared with matched lazy heat that ignores the
record, a delayed record, a shuffled record, an inverted decision, a
state-preparation-only explanation, a calibrated-noise model, and a mixture of
256 possible label/layout maps. Invalid codes and failed feedback stay in the
outcome table; they are never postselected away.

The state-preparation-only null is especially important. It preserves the
final-state histogram and the joint disturbed-state/decision histogram, but
breaks the temporal link connecting the record, decision, and later repair.
Beating it therefore shows that matching static before-and-after histograms is
not enough to explain the observed joint process.

#### Why the blinding matters

Before hardware submission, the catalog froze every circuit, shot count,
backend policy, physical layout, null model, leakage rule, analysis threshold,
and exclusion rule. Semantic group labels, controller identities, backend
roles, and layout identities were hidden behind opaque IDs. `ibm_fez` was run
as the development slot. Only after all its jobs completed was the held-out
`ibm_kingston` slot submitted. The blind report was written and file-hashed
before the private reveal was opened.

#### What happened

- 14 IBM Runtime jobs completed: 3,104 circuits and 606,208 shots.
- Every submitted job and expected row was included, with no dropped or
  postselected shots.
- After unblinding, the reveal commitment and all 3,104 logical-circuit hashes
  verified.

The primary statistic compares the probability of all observed joint outcomes
under record-gated repair with their probability under each frozen null. A
positive log-likelihood ratio favors repair; larger values mean greater
separation. The required thresholds were `log(10)` per backend and `log(100)`
pooled.

| Frozen null | `ibm_fez` log ratio | held-out `ibm_kingston` | Pooled |
| --- | ---: | ---: | ---: |
| lazy heat, which ignores the record | `136,738.537` | `141,443.849` | `278,182.386` |
| delayed record | `77,145.222` | `78,954.396` | `156,099.617` |
| shuffled record | `93,252.706` | `96,396.875` | `189,649.582` |
| inverted decision | `275,461.330` | `284,142.447` | `559,603.778` |
| state preparation only | `99,891.886` | `103,598.342` | `203,490.228` |
| calibrated noise | `171,961.767` | `180,409.442` | `352,371.209` |
| mixture of 256 label/layout maps | `5.545` | `5.545` | `5.545` |

The smallest entry, `5.545`, corresponds to an ordinary likelihood ratio of
about `256`. The preregistered reference mapping was also the unique best one
among the 256 candidates. The much larger process ratios accumulate evidence
over many shots and many deliberately distinct circuit variants; they should
not be compared directly with a single small laboratory measurement.

The run also passed all frozen quality gates:

- the largest four-bit basis error was `3.516%` on `ibm_fez` and `2.930%` on
  `ibm_kingston`, below the `15%` limit;
- pooled Cayley leakage was `3.362%` and `1.189%`, below the `10%` limits;
- the worst individual circuit had `19.792%` leakage, retained in the data and
  below the separate `25%` catastrophic limit;
- every comparison still passed when the frozen contamination fraction was
  varied from `4%` through `15%`, around the central `8%` model; and
- the complete primary S3 transition law had zero violations in a simultaneous
  99% prediction envelope covering 55,296 locked outcome cells.

That last number counts the many cells in the locked transition table, not
55,296 independent experiments. “Zero violations” means that none fell outside
the tolerance calibrated for the whole family at once.

The full result, jobs, raw counts, checksums, and exact replay commands are in
[`runs/issue_509_20260711_v2/`](runs/issue_509_20260711_v2/).

The very large log-likelihood ratios need context. They accumulate over many
shots and many sharply separated, nearly deterministic circuit variants. They
show that the hardware data match the programmed conditional controller far
better than the listed frozen controllers. They are not odds that OPH is true,
and they are not Bayes factors.

#### How strong is it?

This is the strongest evidence in the directory for its **immediate
implementation claim**. It is preregistered, blinded, multiplicity-corrected,
held out on a second QPU, fail-closed on calibration and leakage, and fully
receipt-backed. It strongly supports the statement that IBM hardware realized
the programmed finite self-reading repair instrument and that this instrument
outperformed the restricted nulls. The held-out chip is a real replication
barrier, but both chips came from one provider and one campaign, so this does
not yet establish robustness across providers, hardware generations, or time.

It does not support the stronger statement that OPH beat standard quantum
mechanics. Unrestricted quantum mechanics includes measurement, classical
records, and feedback, so it predicts the same programmed instrument. The run
also says nothing direct about OPH relativity, Standard Model, cosmology, or
ontology claims.

### 6. Secondary dimension-exponent diagnostic

#### What it asks

The companion Metropolis--Hastings circuits test three programmed stationary-
law conventions. `kappa = 1` gives one factor of representation dimension,
`kappa = 0` removes that factor, and `kappa = 2` uses a squared or Plancherel-
like factor. S3, A4, and a seeded random spectrum contain unequal dimensions,
so the choices can be distinguished. Z3 and Z5 have only one-dimensional
sectors, so all three choices must be identical; those are negative controls.

#### What happened and how strong is it?

On the two backend slots, the identifiable S3, A4, and seeded-random families
favored the programmed `kappa = 1` law over `kappa = 0` and `kappa = 2`. The Z3
and Z5 negative controls returned exactly zero log-likelihood separation, as
required. Fourteen Holm-corrected secondary goodness-of-fit tests produced no
rejection at family alpha `0.01`.

| Pooled family | log ratio: `kappa=1` over `kappa=0` | log ratio: `kappa=1` over `kappa=2` |
| --- | ---: | ---: |
| S3 | `46.216` | `83.635` |
| A4 nonabelian decoy | `28.276` | `14.092` |
| seeded random nonabelian decoy | `189.026` | `926.073` |
| Z3 abelian negative control | `0.000` | `0.000` |
| Z5 abelian negative control | `0.000` | `0.000` |

This is **supportive diagnostic evidence** that the circuit family and analysis
correctly distinguish the dimension conventions when they are identifiable
and refuse to do so when they are not. It did not select the primary verdict,
and it does not derive nature's stationary law or independently confirm OPH.

## Why this program is frozen

For every experiment above, standard quantum mechanics predicts the same
probabilities as the advertised OPH-motivated target:

- Stage 1, Z3, Z5, and S3 place the target structure into a state prepared with
  standard quantum gates, then measure it with standard quantum readout.
- The generative benchmark is a standard dynamic quantum circuit: measurement
  creates a classical record and standard conditional feedback acts on that
  record.
- Its impressive likelihood ratios compare the implemented controller with
  deliberately restricted controller and noise nulls. They do not compare an
  OPH probability law with a different quantum-mechanical probability law.

This is a missing-model-contrast problem, not a statistical-power problem.
More shots, qubits, backends, or elaborate blinding would make the engineering
replication stronger but would not create evidence that separates OPH from
quantum mechanics. Hardware noise and provider calibration would matter once a
specific deviation existed, but they are not the reason for the present
freeze: there is currently no different OPH number to measure.

The hardware program may be reopened only after a discriminator proposal
provides all of the following before platform selection or data collection:

1. a source-closed observable for which OPH and standard quantum mechanics
   give different numerical predictions;
2. the sign and magnitude of the predicted difference as a function of the
   controlled experimental parameters;
3. a nuisance and calibration model showing that the difference is identifiable
   above the proposed platform's noise floor;
4. an audit against existing experimental bounds, especially matter-wave,
   optomechanical, collapse-model, gravitational-response, and spontaneous-
   emission constraints when relevant; and
5. a preregistered decision rule, power calculation, and public evidence plan.

Until that gate is met, maintenance, security fixes, and exact replay of the
archive are welcome, but no new IBM Quantum run should be presented or funded
as evidence for OPH over standard quantum mechanics.

## Where the code and data live

- `programs/ibm_runtime_common.py`
  Shared IBM Runtime helper utilities.
- `programs/check_ibm_setup.py`
  Auth and backend-discovery verifier.
- `programs/discrete_heatkernel_test.py`
  `Z_3`, `Z_5`, and `S_3` reduced-sector runner.
- `programs/s3_diagnostic_bundle.py`
  `S_3` layout and readout diagnostic bundle.
- `programs/stage1_markov_fingerprint.py`
  3-qubit Markov / recoverability benchmark.
- `programs/generative_repair_kernel.py`
  Pure finite-matrix implementation of the record-gated Cayley repair kernel,
  the matched open-loop heat null, and the dimension-exponent diagnostic.
- `programs/record_gated_cayley_circuits.py`
  Local dynamic-circuit preflight for the self-reading repair instrument. It
  contains no IBM submission path.
- `programs/generative_repair_circuits.py`
  Secondary Metropolis--Hastings `kappa = 0, 1, 2` diagnostic circuits,
  including mandatory abelian non-identifiability controls.
- `programs/blind_preregister.py`
  Opaque circuit catalog, private commit/reveal payload, and digest verifier.
- `programs/blind_preregister_orchestrate.py`
  Two-phase production seal: catalog precommitment followed by an analysis-lock
  binding and atomically committed public/private artifacts.
- `programs/cayley_blind_likelihood_analysis.py`
  Frozen joint-process conditional likelihood ratios, contamination
  sensitivity, calibration/leakage gates, and blinded verdict logic.
- `programs/record_gated_cayley_runtime.py`
  Explicit role-scoped IBM plan, submit, and harvest path with append-only
  hash-chained journals.
- `programs/runtime_analysis_packet.py`
  Fail-closed conversion of the complete two-role runtime journals into the
  blinded analysis packet.
- `protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md`
  Claim boundary, exact predictions, blinding contract, nulls, and failure rule
  for the issue `#509` benchmark.
- `runs/issue_509_20260711_v2/`
  Public manifest, committed reveal, raw provider count journals, sealed blind
  packet/report, job receipts, checksums, and the completed hardware result.
- `qc_data/README.md`
  Index of the representative legacy derived summaries.

## Technical details and local verification

The exact frozen contract is in
[`protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md`](protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md),
and the public run README gives the hashes and end-to-end replay commands. Each
dynamic variant used 192 shots per backend slot; separate four-bit calibration
circuits used 512 shots per prepared code. The calibration circuits are
validity diagnostics, not a factorized correction applied to the feedback
data. Development had to finish before held-out submission, and both roles
remained bound to the same sealed catalog and analysis lock.

### Local preflight

Create the pinned optional environment:

```bash
python3 -m venv /tmp/oph-ibm-quantum
/tmp/oph-ibm-quantum/bin/python -m pip install -r code/ibm_quantum_cloud/requirements-ibm.txt
```

Run the pure finite-kernel and optional circuit tests:

```bash
python3 -m pytest -q code/ibm_quantum_cloud/tests/test_generative_repair_kernel.py

/tmp/oph-ibm-quantum/bin/python -m pytest -q \
  code/ibm_quantum_cloud/tests/test_generative_repair_kernel.py \
  code/ibm_quantum_cloud/tests/test_record_gated_cayley_circuits.py
```

Generate local validation receipts without credentials or hardware submission:

```bash
python3 code/ibm_quantum_cloud/programs/generative_repair_kernel.py \
  --json-out /tmp/oph-generative-kernel-validation.json

/tmp/oph-ibm-quantum/bin/python \
  code/ibm_quantum_cloud/programs/record_gated_cayley_circuits.py \
  --model s3 \
  --protocol record_gated \
  --shots 192 \
  --json-out /tmp/oph-s3-record-gated-preflight.json
```

These local receipts validate an engineered record-conditioned process against
frozen open-loop and stationary-law nulls. The hardware workflow additionally
requires the two-phase seal, role-scoped sequential execution, complete
hash-chained journals, blind conditional-likelihood report, and verified
reveal required by the issue `#509` acceptance protocol.

A passing hardware result supports only the realization of a blinded finite
self-reading repair instrument against its frozen controller, label, noise,
and decoy nulls. Standard quantum mechanics predicts the programmed instrument,
so the result cannot distinguish OPH from unrestricted quantum mechanics and
does not establish an OPH vacuum, relativity, Standard Model structure, or
cosmology.

## Archived rerun instructions

These commands are retained for reproducibility and maintenance. Their
presence does not reopen the frozen OPH-evidence campaign.

The credential file is not part of this repo. The scripts accept either a raw
IBM Cloud API key or a local note file with a line of the form:

```text
IBM cloud API key: <YOUR_API_KEY>
```

Pass it explicitly with `--credentials-file` if it is not in the current working directory.

Example commands:

```bash
python3 code/ibm_quantum_cloud/programs/check_ibm_setup.py --credentials-file /path/to/IBM_cloud.txt

python3 code/ibm_quantum_cloud/programs/discrete_heatkernel_test.py \
  --group z3 \
  --mode hardware \
  --backend ibm_marrakesh \
  --shots 4096 \
  --outdir code/ibm_quantum_cloud/qc_data/z3/rerun_example \
  --credentials-file /path/to/IBM_cloud.txt

python3 code/ibm_quantum_cloud/programs/s3_diagnostic_bundle.py \
  --mode hardware \
  --backend ibm_marrakesh \
  --t-value 0.60 \
  --shots 8192 \
  --use-best-pair-layout \
  --outdir code/ibm_quantum_cloud/qc_data/s3/rerun_example \
  --credentials-file /path/to/IBM_cloud.txt
```

## Notes

- The legacy JSON exports keep original job IDs and the reported derived
  metrics, but account-specific instance CRNs are redacted in the public copy.
- Complete provider count dictionaries were not archived for most Stage 1,
  Z3, Z5, and direct S3 runs. The S3 layout summaries retain target counts and
  fitted assignment matrices, but not every calibration count dictionary.
- The issue `#509` v2 directory is different: it retains the complete raw
  joined counts for all 3,104 circuits, provider metadata, hash-chained
  journals, the blind packet and report, and everything needed for exact
  public replay.

## License And Patent Policy

This code and data bundle is part of the OPH public repository. See the main
[LICENSE](../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../PATENTS.md).
