# Faster complete particle compilation

Exact native momentum factors give a complete-output compiler error that
scales linearly with volume. The same localized particle energy and detector
records therefore require far fewer certified cycles, while retaining all
ancillary outputs and quantum references. The calculation uses the specified
coherent controls and does not provide an autonomous microscopic clock.

See the [derivation](PROOF.md) and [replay](results.json).
Run `python3 -m model.openwave_factorized_compiler --check`.
