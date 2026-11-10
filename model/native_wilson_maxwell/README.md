# A finite charged-current and field interaction

Actual native charged states supply a complete fermion occupation code.
An explicit spinor current couples it to a finite Maxwell field on an
open cube. The calculation bounds a nonzero pair-plus-field record from
a prepared free state under the coupled Hamiltonian, retaining every
charge, field and preparation complement. It is a finite nonstationary
interaction, without a claim of relativistic scattering or observed
particle identity.

The [proof](PROOF.md) specifies the common source and complete records.
[The current control](CURRENT.md) explains why a scalar-band Peierls
coupling cannot replace this spinor current. From the repository root run
`python -B -m model.native_wilson_maxwell --write`,
`python -B -m model.native_wilson_maxwell --check` and
`python -m pytest -q model/native_wilson_maxwell` with the scientific
environment. The receipt replays both independent current calculations
and the shared finite field provider.
