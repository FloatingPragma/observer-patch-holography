# Evaluating a candidate

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

The remaining judgments are evidence-backed assessments. The evaluator checks their vocabulary, identifiers, premise crosswalk and evidence paths. The linked proofs and experimental receipts establish their scientific scope. Reproducing an experimental result requires running its instrument. For example:

```sh
python -m model run --law interval_echo --full
python -m model run --law contrast_memory_spherical_v1
python -m model run --law scalar_channel --extended
python -m model.experiments.cage
python -m model.evaluation.repair_band
python -m model.evaluation.patch_spatial_readout
python -m model.evaluation.observer_encounter
```

The band check measures the operator that carries the three slow directions. The patch readout reconstructs them from relabelled, exact port-detection probabilities and tests two held exposures. The encounter instrument evolves the complete two-load sector of the 20-cell sphere from a uniform ensemble of separated starts. It reports both the first detection rate and persistence conditional on that detection for the native contrast law and an exchange-only control.

The archived eight-slot routing experiments supply the historical comparison. On the macOS arm64 replay, two tests still report seeded count differences in the electron-wall and loaded-clock archives. The historical rows identify the published results and apparatus. The three current candidate replays are independent of those count archives. [Electron-wall replay](../benchmarks/eight_slot_torus/electron_wall/test_electron_wall.py) · [Loaded-clock replay](../benchmarks/eight_slot_torus/formed_clocks/test_loaded_clock_archive.py)

## What is scored

Where a measurement demonstrates one part of a criterion, its scorecard cell names that measurement and links to its experiment. The separate full-criterion grade determines the category totals.

The dark-sector test starts at the Einstein tower. On one finite overlap, vary the signed modular charge, compare with a fully recovered control and independently measure area and proper volume using the same vacuum and clocks. Repeat the curvature measurement with ordinary matter and a second radius before transferring the coupling to galaxies. Rotation curves and shear can then test a profile with frozen normalization. To identify the origin of an excess, this comparison also needs an indicator of overlap recoverability determined independently of the gravity data. OW-19 asks for a neutral massive excitation. Its binding measurement has a separate physical interpretation.

| Source | Coverage |
| --- | --- |
| [Theorem criteria](theorem_criteria.json) | A1–A3, finite repair, consensus, quantum probabilities, current and matter premises, G1–G6 and their joint null balance |
| [Observer criteria](observer_criteria.json) | Seven historical field, clock and pair experiments, plus finite patch measurements with their preparations and transfer limits |
| [OpenWave criteria](openwave_criteria.json) | One row for each of 31 phenomena in the pinned OpenWave table, with its upstream number, dynamics, accepted outcome and local candidate evidence |
| [Particle criteria](particles.json) | A discriminator for each of 17 Standard Model particle types. A candidate can declare a species assignment in `PARTICLE_IDENTIFICATIONS` with a grade, evidence and explanation |
| [Premise register](../../tracking/premise_register.json) and [roles](roles.json) | Complete formal assumption inventory and its law, designer or import owner. Each plugin declares evidence in `PREMISES` and same-process claims in `JOINT` |

`--write` generates the scorecard, evidence ledger, candidate index table and per-candidate `EVALUATION.json`. `--check` regenerates those views in memory and compares their contents. A new plugin defaults to open on reviewed criteria until it has evidence in the correct physical domain. Add narrowly scoped positive results where the data support them, and distinguish calibrations, held tests, earlier benchmark constructions and physical identifications.
