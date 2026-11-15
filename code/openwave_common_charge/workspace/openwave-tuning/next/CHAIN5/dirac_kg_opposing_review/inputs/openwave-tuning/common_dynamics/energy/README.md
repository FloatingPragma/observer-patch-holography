# Energy and spectral records

The native walk has complete phase-estimation records and a conserved positive
spectral observable with local density and flux. Its quasienergy is measured
relative to the shared clock, with a declared logarithm branch. This construction
supplies neither a unique vacuum nor a finite positive continuum rest mass.

Run `python energy.py` for the source-bound calculation and `python freeze.py`
for transitive custody verification. The [proof](PROOF.md) gives the complete
instrument and its scope. [ENERGY_CONTRACT.json](ENERGY_CONTRACT.json) contains
all phase probabilities, native fine-pointer factors, source error and reset
budgets, local currents and rejecting controls.

The main interfaces are `phase_estimation(state,bits)`,
`quasienergy_labels(bits,h,tau0)`, `qpe_contract(bits)` and
`positive_witness(state)`. The full source array is always128 by16. A phase bin
is a finite-resolution random record, not an exact eigenvalue assertion.
Field and many-body extensions must control their entire augmented evolution
and preserve its own reference, wire allocation and complete outcomes.
