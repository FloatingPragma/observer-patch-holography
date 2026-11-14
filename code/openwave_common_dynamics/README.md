# Shared matter, field and clock records

This package regenerates finite simulation-generated records from one frozen
massive matter-field source with a native clock and complete position/field
readout. It also executes separately typed quantum-cycle, reference-preparation
and two-charge phase checks. These experiments have different preparations
and apparatus domains; their inclusion does not establish one interacting
continuum theory or validate the 31 physical criteria.

Run with **Python 3.13** and the pinned dependencies in
[requirements.txt](requirements.txt). The numerical libraries are NumPy,
SciPy, SymPy and mpmath. The frozen carrier also imports NetworkX and schema
validation libraries through its original package initializers.

```sh
python3 code/openwave_common_dynamics/replay.py --output /tmp/common-replay.json
```

Each experiment runs in a fresh temporary workspace under `-I -S`, with one
numerical thread and an independent file-access guard. The runner preserves
the selected interpreter path, adds only explicit runtime site directories,
and executes no `.pth` startup files. The repository's Python 3.12 virtual
environment is outside this replay contract: `-S` changes its site discovery.

The [contract](CONTRACT.json) distinguishes the source, state, clock and
observable domains. The [proof](PROOF.md) explains the execution and error
certificates. [SNAPSHOT_MANIFEST.json](SNAPSHOT_MANIFEST.json) binds a flat
workspace with the original relative paths and exact source bytes. No live
project source or symlink fallback is used. Exact repeated-block certificates
replace the expansion of millions of primitive pulses.

Use `--mode integration`, `quantum`, `reference` or `two_charge` to replay one
experiment. The compact two-charge export excludes the large mobile-trajectory
archives listed in the snapshot manifest. Its independent matrix replay
regenerates all 64 phase-instrument outcomes.
