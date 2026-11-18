# Finite memory and an autonomous event clock

This experimental construction carries the directed-bit source along a finite
memory rail. One time-independent Hamiltonian performs its reversible updates
and stores each address and overwritten bit. Conditioning on a measured clock
position gives the exact repeated native record law, including no-operations.
The same dynamics reverses and uncomputes its records. This is a finite
computation with an internal ordinal clock; physical particle identification
and a permanent observer archive require further construction.

The [proof](PROOF.md), [specified experiment](protocol.json) and
[complete records](results.json) define this result. `OPH-NMR-1` is a separate
experimental construction. The official calibrated CCF-LT-2 source and evidence
are unchanged. Source ancestry and theorem boundaries are in
[provenance.json](provenance.json).

## Offline replay

From the repository root, with `model/requirements.txt` installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.native_memory_clock.experiment
python -B -O -m unittest model.native_memory_clock.test_clock model.native_memory_clock.test_experiment
```

The headline calculates all four two-event source/background cases, every
terminal outcome and 3,072 sampled outcome indices. These samples are separate
preparations, with shared random uniforms for matched controls. The fixed
simultaneous concentration bound is reported beside each table. A table's
index identifies the complete record, including an unfinished clock.
`native.decode_history` reconstructs the primary seam/coin/current records
and the initial word from the endpoint and archive, without an initial-word
input. The classical archive readout and the inverse-decoder interference
readout are alternative instruments; decoding coherently uncomputes the
archive.

Extended checks are separate:

```sh
python -B -O -m unittest model.native_memory_clock.test_native
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.native_memory_clock.experiment --extended
```

They enumerate the full dirty-memory gate domain and calculate the three- and
four-event rails. The compressed extended receipt contains every computed
terminal record law; decompression needs only Python's standard library.
The full-carrier spectral theorem retains all modes and their exact
multiplicities. Numerical propagation uses a verified reducing span, with no
selected low-energy window or empirical mass input.

Verification is read-only. `--record` creates evidence only when its target
does not exist and the frozen sources match. Source amendments require a
separate freeze; failed evidence is preserved. Passing these finite checks
does not validate any of OpenWave's 31 complete physical rows. Candidate
registration and platform validation belong in Bernhard's OpenWave fork.
