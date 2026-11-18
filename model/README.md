# The executable OPH model

The executable model is **polynomial-sv**. Its law, the contrast-memory law, is a candidate for the native
repair law: twelve-port cells whose seams keep a small private memory of how much their two ends disagree.
The [README](polynomial_sv/README.md) walks through it, the [specification](polynomial_sv/SPEC.md) fixes the
law, the carriers and the records, and the [proof](polynomial_sv/PROOF.md) shows why a bound pair stays a pair.

The experimental [finite clock](contrast_clock/README.md) executes native
contrast-memory gates and retained observer readings with a fixed controller
Hamiltonian. Its source mapping and finite observation interval are explicit;
the controller is an additional declared architecture.

Run it with `pip install -r model/requirements.txt`, then `python -m model` (add `--full` to recompute
everything) and `python -m pytest -q model/polynomial_sv`.

A few small theory packages sit next to it (`support_refinement`, `matched_einstein_limit`,
`defect_particle_classes`, `carrier_selection_reduction`, `finite_instrument_calibration`,
`maxwell_local_selection`, `maxwell_weak_field`); each has its own README and tests.
