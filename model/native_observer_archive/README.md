# Autonomous observers with finite memory

A moving clock packet drives two native source events, writes observer A's
record, carries it to observer B, and lets B's received bit select the second
event's address. One time-independent Hamiltonian performs the whole program.
Uniform local couplings and a supplied compact packet give more than 99.9%
readiness throughout a finite interval, including every unfinished outcome.
The [proof](PROOF.md), [experiment](protocol.json) and
[source mapping](provenance.json) specify this experimental construction.
Physical particle identification and an indefinitely continuing observer
require additional source-specific construction.

`OPH-NOA-1` uses the exact native event gates of
[`OPH-NMR-1`](../native_memory_clock/README.md), with a different clock,
preparation, apparatus and adaptive address rule. The official calibrated
CCF-LT-2 model and its evidence are unchanged. The observers copy orthogonal
bit records; their copies do not clone arbitrary quantum states. Source and
reference coherences within a public pointer outcome are retained.

## Offline replay

With `model/requirements.txt` installed, run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.native_observer_archive.experiment
python -B -O -m unittest discover -s model/native_observer_archive -t .
```

The compressed `results.json.gz` contains every clock probability and all
twelve address branches at each program prefix. Their product reconstructs
the complete diagnostic law without repeating identical register tables.
The public pointer-only measurement is a separate instrument, checked with
retained quantum references. No random samples or successful-arrival
normalization replace the complete laws. Controls disable communication,
publication and feedback, reverse the packet, or change the native event law.
Late readouts retain the loss of the records during reverse passage.

The optional extended command propagates the full 49,152-dimensional prepared
span of the same 4,096-site construction through literal Hamiltonian columns:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.native_observer_archive.experiment --extended
```

The headline includes a separate small literal-column composition check.
The uniform interval bound is analytic; numerical time points do not prove it.
Packet representation has an exact error bound. SciPy propagation residuals
are implementation checks, with no certified floating-point exponentiation
claim. Verification is read-only. `--record` refuses to overwrite evidence;
source amendments require a separate freeze. All 31 complete OpenWave
physical rows are unvalidated; candidate registration and adversarial row
validation belong in Bernhard's OpenWave fork.

The verifier checks each saved scientific gate before numerical comparison.
Versioned producer receipts and their freeze are retained under `receipts/v1`;
their source revision and hashes are specified in the source mapping.
