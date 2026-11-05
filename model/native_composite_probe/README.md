# A detector for charged composite excitations

A coherent native detector measures the cubic particle response that is
orthogonal to all linear fermion probes on the specified interacting
eight-site source. Its 1,008 native pulses preserve the complete field,
all selector outputs and every abort. The measured success probability
agrees with the independent certified 138-line spectral weight.

The [proof](PROOF.md) states the source controls and complete error bound;
the [receipt](results.json) binds the executed detector and lossless output
archive. Resolving individual energies and identifying observed particles
require further physical attachment.

Run `.venv/bin/python -m model.native_composite_probe --check` for a cold
replay, or `--write` to retain its receipt and archive.
