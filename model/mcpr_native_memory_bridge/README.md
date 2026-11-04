# Pixel-dependent charged masses in native memory

The conditional pixel mass law can be realized by the same twelve-port
memory and native gauge-current probes. Independent probe durations and a
declared feedback amplitude produce a matrix with exactly the three
charged-lepton coefficients of that law. Every repair record is retained.
A bound on the full correlated memory proves the corresponding Hamiltonian
in the weak-pulse limit.

The architecture is a chosen simulator design with target-informed mass-law
provenance. Its consequences are source-computable; identifying its
coefficients with physical interacting particle masses requires additional
physics. See [the proof](PROOF.md) and [the receipt](results.json).

Run `python -m model.mcpr_native_memory_bridge --check` and
`pytest model/mcpr_native_memory_bridge` from the repository root.
