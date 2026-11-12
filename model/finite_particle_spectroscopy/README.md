# Finite phase-pointer mass records

A complete 32-outcome detector measures the existing calibrated source through
its actual energy evolution, retaining every probe return, phase bin,
inconclusive result and many-excitation complement. It preserves color
and unresolved family coherence. Its mass confidence sets have ideal
coverage at least 65/72; source calibration errors and hardware errors
remain separate.

An exactly equivalent adaptive implementation uses one reusable quantum
qubit, five retained classical bits and five sequential source evolutions.
Every original outcome and source/reference state is preserved.

The [receipt](results.json) supplies the complete law, a 17-type conditional
GeV menu and sufficient ideal precision/exposure budgets. The large
planned registers are not executed, and the executed 32-outcome pointer does not
resolve the neutrino mass species. See the [proof](PROOF.md) and
[code](spectroscopy.py).

Run `python3 -m model.finite_particle_spectroscopy --check` from the repository root.
