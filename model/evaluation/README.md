# Evaluating a candidate

This folder grades candidate laws. The program-wide status of every observable and theorem premise is the [completion ledger](../../docs/registers/README.md).

The [scorecard](SCORECARD.md) is the official evaluation: it compares every registered law against every requirement. Each criterion ID opens its acceptance test, the reason for each grade and the original evidence in the [evidence ledger](EVIDENCE.md). Overlapping observers own their local physical registers. A material-clock test specifies which patch reads each register and when.

Run these commands from the repository root with the model dependencies installed.

```sh
python -m model evaluate
python -m model evaluate --criterion OBS-02
python -m model evaluate --criterion PR-54
python -m model evaluate --criterion JOINT-patch-measurement
python -m model evaluate --law interval_echo --criterion G2
python -m model evaluate --write
python -m model evaluate --check
python -m pytest -q model/evaluation
```

## What is executed

On every evaluation the common harness constructs each candidate’s cell generator and checks hermiticity, native coupling supports, the action of the 60 proper rotations, and the multiplicity of the first nonconstant response level. It builds the candidate on the 20-, 80- and 320-cell federations where implemented and validates basis coverage and coupling locality. It checks parameter declarations and reads the headline numerical quantities from the measured files. The JSON receipts retain the numerical values, standard errors when supplied, and exact source pointers.

A candidate can deposit measured premise assessments in `premise_results.json` beside its law module. Each entry gives the premise ID, status, evidence paths relative to the candidate, and scope in the same format as `PREMISES`. To revise an existing law-module assessment, include `supersedes_sha256` with the SHA-256 of the prior entry serialized as compact, key-sorted JSON (ASCII escaping, UTF-8 bytes). The evaluator retains the earlier evidence and scope and rejects ambiguous or stale revisions. This allows a run to pin the law module while subsequent measurements grade its premises.

Each candidate retains its recorded criterion and premise assessments. The programme register has its own grades. When a candidate supplies a structured full-credit certificate, the evaluator checks the evidence, stated scope and listed premises. The rendered first label always gives the full criterion grade. Additional measurements appear after it with their apparatus and evidence.

The scalar and contrast-memory results have separate apparatuses. The finite band test uses the candidate's declared register representation and the first distinct spectral level above its ground space. A sharp port probability has the same value under H and −H*. On a real uniform hopping this prevents a sign or oriented slow-energy assignment from positions alone. A calibrated phase-sensitive preparation can separate them. [Port reconstruction and coherence control](patch_spatial_readout.json)

The remaining judgments are evidence-backed assessments. The evaluator checks their vocabulary, identifiers, premise crosswalk and evidence paths. The linked proofs and experimental receipts establish their scientific scope. Reproducing an experimental result requires running its instrument. For example:

```sh
python -m model run --law interval_echo --full
python -m model run --law contrast_memory_spherical_v1
python -m model run --law scalar_channel --extended
python -m model run --law private_record_walk --check
python -m model.experiments.cage
python -m model.evaluation.repair_band
python -m model.evaluation.patch_spatial_readout
python -m model.evaluation.observer_encounter
```

The band check measures the operator that carries the three slow directions. Spectral multiplicities use energy differences and a tolerance relative to the centered operator. The rotation check compares products on rays, so a common phase cancels and relative phases between registers enter the comparison.

The patch readout reconstructs the slow directions from relabelled, exact port-detection probabilities and tests two held exposures. The records admit both $H$ and $-\overline H$ because their propagators are complex conjugates in the port basis. With a supplied phase reference, the preparation $(|0\rangle+i|1\rangle)/\sqrt{2}$ distinguishes them. At twice the calibration exposure the two predicted probabilities at port zero are 0.353601 and 0.585897. The reference belongs to this control. The original instrument uses sharp port preparations. [Patch readout](patch_spatial_readout.py)

The encounter instrument evolves the complete two-load sector of the 20-cell sphere from a uniform ensemble of separated starts. It reports both the first detection rate and persistence conditional on that detection for the native contrast law and an exchange-only control.

The archived eight-slot routing experiments supply the historical comparison. On the macOS arm64 replay, two tests still report seeded count differences in the electron-wall and loaded-clock archives. The historical rows identify the published results and apparatus. The registered candidate replays are independent of those count archives. [Electron-wall replay](../benchmarks/eight_slot_torus/electron_wall/test_electron_wall.py) · [Loaded-clock replay](../benchmarks/eight_slot_torus/formed_clocks/test_loaded_clock_archive.py)

## What is scored

Each grade applies to the named registered law. A row can link multiple comparisons, each with its own construction, finding and evidence. The photon row gives separate entries for the written field and the controlled quantized Maxwell box. The evidence ledger gives its law, its finding and its source alongside the registered result. For Interval Echo the written-field comparison traces a signed load, a field received in a distant cell and the orbital states that can accompany the load. The pinned-law dark-contact check uses the original 99-state operator. [Electron and related row studies](../candidate_laws/interval_echo/openwave/electron_field/README.md) · [Pinned dark-contact calculation](../candidate_laws/interval_echo/openwave/electron_field/dark_contact.py)

Where a measurement demonstrates one part of a criterion, its scorecard cell names that measurement and links to its experiment. The separate full-criterion grade determines the category totals.

The dark-sector test starts at the Einstein tower. One finite experiment now varies a signed modular charge at fixed full stored-energy distribution. Use those held preparations and a fully recovered control for the curvature comparison. The meter must extract area and proper volume from a shared vacuum and calibrated clocks while retaining access to interior modes invisible to a single-load exterior history. On the same vacuum and matter perturbation, measure the regional bulk entropy, area and proper-volume changes in ten independent timelike directions. Test whether the residual divided by the fourth power of the proper radius tends to zero as the regions shrink. Calibrate the gravitational coupling on ordinary matter, freeze it and predict the curvature change for the independently measured collar charge. Rotation curves and shear can then test a profile with frozen normalization. To identify the origin of an excess, this comparison also needs an indicator of overlap recoverability determined independently of the gravity data. OW-19 asks for a neutral massive excitation. Its binding measurement has a separate physical interpretation.

| Source | Coverage |
| --- | --- |
| [Theorem criteria](theorem_criteria.json) | A1–A3, finite repair, consensus, quantum probabilities, current and matter premises, G1–G6 for the ten-timelike-observer Einstein theorem |
| [Observer criteria](observer_criteria.json) | Seven historical field, clock and pair experiments, plus finite patch measurements with their preparations and transfer limits |
| [OpenWave criteria](openwave_criteria.json) | One row for each of 31 phenomena in the pinned OpenWave table, with its upstream number, dynamics, accepted outcome and local candidate evidence |
| [Gravitational particle](gravitational_particle.json) | Separate assessment of the hypothetical graviton. The native contact and ideal tidal readers test different parts of the assignment |
| [Particle criteria](particles.json) | A discriminator for each of 17 Standard Model particle types. A candidate can declare a species assignment in `PARTICLE_IDENTIFICATIONS` with a grade, evidence and explanation |
| Premise rows of the [observables register](../../docs/registers/OBSERVABLES.json) and [roles](roles.json) | Complete formal assumption inventory and its law, designer or import owner. Each plugin declares evidence in `PREMISES` and same-process claims in `JOINT` |

The six gravity grades use one source, common proper balls, and ten independent timelike observers. The three small-ball errors concern bulk entropy and independently measured matter stress, area and curvature, and entropy stationarity at fixed proper volume. Each error must vanish faster than the fourth power of the proper radius. Comparing opposite light rays can separately test pressure, flux and focusing on a calibrated optical instrument. Those comparisons have no score of their own. [Ten-observer theorem and normalizations](../../code/matched_einstein_limit/PROOF.md) · [Optional light-ray measurements and controls](null_direction_comparisons.json)

`--write` generates the scorecard, evidence ledger, candidate index table and per-candidate `EVALUATION.json`. `--check` regenerates those views in memory and compares their contents. A new plugin defaults to open on reviewed criteria until it has evidence in the correct physical domain. Add narrowly scoped positive results where the data support them, and distinguish calibrations, held tests, earlier benchmark constructions and physical identifications.
