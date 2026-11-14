# Local construction integration gates

Run from any directory:

```sh
python /path/to/bundle/integration/run.py
```

Requires NumPy and SciPy. Portable geometry and native events come from the candidate package; no RER archive is required.
Use `--output PATH` to save machine-readable results. The fixed random seed and imported source hashes make reruns comparable.

Three issues closed here:

1. Executable algebra gates for the massive port/seam generator and independent-load composition.
2. Exact signed-current Gauss identity and charge conservation through 10,000 native events on the canonical ring federation.
3. Reproducible result emission with provenance and a guard against promoting running/interim work into any of the 31 official cells.

These are local regression checks on already explored configurations. They do not consume frozen acceptance configurations and cannot confer an OpenWave pass. Fable's separate acceptance freeze is maintained independently. `local_result.json` records a local run only.

## All-row dispatch and source-bound receipts

Run from the bundle root (`code/openwave_candidate` in RER or `openwave-tuning` in the development workspace).

```sh
python integration/assess.py --output integration/assessment.json
python integration/assess.py --verify integration/assessment.json
```

`assess.py` dispatches all 31 frozen criteria. Three finite-construction probes give partial support for integer charge, internal oscillations and the massive wave relation. The joined controls supply diagnostics for source persistence (row 4), supplied static flux and rotation (rows 6 and 7), selected record energy (row 20), and integer string energy (row 22). The other 23 rows report missing physical observables. None of these diagnostics identifies a physical particle or closes its complete criterion. Every row retains its official status and `full_criterion_met: false`. In particular the internal oscillation is at 2E; it is not credited as the requested electron rest clock.

Each probe also rejects a deliberately broken control: a missing recorded current, an incorrect oscillation frequency, and an incorrect rest-block sign. The receipt binds the executing script, implementation, frozen platform packet, portable geometry and joint record implementation by SHA-256. Verification recomputes all probes and rejects any changed source, row, measurement or incomplete source list. It requires the recorded Python/NumPy environment and exact deterministic results; a different environment should generate its own receipt for comparison. Persisted integrity tests run with `python -m unittest discover -s integration -p test_assess.py`. This does not replace independent acceptance on fresh configurations.

The regional apparatus, canonical general-cut identity and local-coarse modular
control are retained under `record_extensions` and linked from row1. They
close operational prerequisites while preserving the physical-charge and
admitted-sector gate. The modular candidate has an explicit C/E instrument and
per-seam rate; its results do not inherit integer-charge or global-monitoring
claims. Both source domains are recorded for comparison.
