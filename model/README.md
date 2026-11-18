# The executable OPH model

The executable model is **polynomial-sv**. Its law, the contrast-memory law, is a candidate for the native
repair law: twelve-port cells whose seams keep a small private memory of how much their two ends disagree.
The [README](polynomial_sv/README.md) walks through it, the [specification](polynomial_sv/SPEC.md) fixes the
law, the carriers and the records, and the [proof](polynomial_sv/PROOF.md) shows why a bound pair stays a pair.

Run it with `pip install -r model/requirements.txt`, then `python -m model` (add `--full` to recompute
everything) and `python -m pytest -q model/polynomial_sv`.

The [construction manual](construction_manual/README.md) shows the hardware and the two laws under test on one
sheet and describes them part by part, with the reason for each design choice and every equation in plain
words.

Two experimental packages extend the model, each with its own README, proof and tests. The
[scalar channel](scalar_channel/README.md) executes the second law of the manual on one cell: local preparation
and observer instruments. The [finite clock](contrast_clock/README.md) runs a fixed program of contrast-memory
events and two observer readings under one time-independent law, with a clock register as the sequencer of the
computation. Their results stand apart from polynomial-sv and support no claim about it without an explicit
implementation proof.

A few small theory packages sit next to it (`support_refinement`, `matched_einstein_limit`,
`defect_particle_classes`, `carrier_selection_reduction`, `finite_instrument_calibration`,
`maxwell_local_selection`, `maxwell_weak_field`); each has its own README and tests.
