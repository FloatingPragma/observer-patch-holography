# Native transfer and finite energy gaps

The actual repair process has a positive stationary population transfer.
Its path probabilities satisfy finite reflection positivity, and its
logarithm defines a finite Hamiltonian with explicit rational gap bounds.
The computation uses the same 240-, 960- and 3,840-port graphs and the
66-state twelve-port exclusion sector. It requires no fitted spectrum or
dense eigensolver. A physical clock and continuum field limit are separate
requirements.

Run `python3 -m model.openwave_native_euclidean_transfer --check`.
The [proof](PROOF.md) states the exact population marginal and its relation
to the complete observer instrument. The [receipt](results.json) contains
all finite bounds and the failing nonlazy control.
