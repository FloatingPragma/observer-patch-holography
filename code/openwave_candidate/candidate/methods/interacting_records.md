# Selected record energy in the native monitored instrument

`InteractingRecords` subclasses `JointRecords`; its native event target map,
C/E records, Gauss sector, integer cutoff and `step` channel are inherited.
The supplied configuration specifies energy kind, coefficient g and spectral
energy mass m_energy. The coherent Hamiltonian is

```text
H_selected = H_joint + g diag(U_phi(E(z,F)));
E_e(z,F) = (F_ea - F_eb)/2.
```

The record-energy functional and transit extension are defined in
[record_operators.md](record_operators.md). This selected diagonal operator
preserves the Gauss sector and changes the generator. The free massive-square
and bare-clock identities are not transferred to H_selected.

`coherent_diagonal_addition` exposes the selected extra diagonal to the existing
`ObserverProtocol` sparse builder. Dense `JointRecords.step` and sparse observer
evolution therefore use the same selected H. The regression compares both
matrices and their evolved pure states and rejects adding the energy twice.
At zero coupling the selected finite channel agrees exactly with the free
`JointRecords.step` channel.

For each step, `step_with_energy` retains the initial energy, energy after the
coherent interval, and final energy. With U=exp(-it H_selected), the coherent
interval preserves Tr(rho H_selected). The monitored instrument generally
changes it. For each resolved committed outcome i, with probability p_i and
basis target t_i, the energy ledger assigns

```text
weighted_output_energy_i = p_i (H_selected)_(t_i,t_i).
```

The no-commit outcome contributes `(1-p) Tr(U rho U^dagger H_selected)`.
Summing all ledger entries reconstructs the final expectation exactly to
numerical tolerance. Monitoring energy change is the final expectation minus
the expectation after coherent evolution. This includes measurement
backaction and native placement. No heat/work split, apparatus energy model,
or monitored energy-conservation claim is made. Boundary rejections and
transit noops retain their outcomes and energy contributions.

The reproduced variants use the oriented triangle, integer cutoff 1, initial
port 0 with zero records, propagation mass 0.23, interval 0.37, commit probability
0.29, energy coupling 0.29 and spectral energy mass 0.23. Energy kinds are
`first_order` and `massive`. The equal numerical propagation and energy masses
are explicit choices with different operator roles. The first-order kind
ignores m_energy and has zero cycle-energy weight; the massive kind has
positive cycle-energy weight m_energy/2.

The returned measurements include complete observer position probabilities,
port/transit and boundary-attempt readouts, change from the zero-interaction
channel, all energy expectations, trace and ledger probability. Rejected
controls test the false conservation of energy through monitoring and the
false reuse of the free channel at nonzero coupling.

This supplies an executable declared interaction on the same finite instrument.
The default candidate has no selected record energy. The finite variants do
not establish a physical charge dictionary, calibrated electric coupling,
force law between moving charges, confinement tension, or string breaking.
Rows 20/22 receive interaction-domain tests; row 30 requires validation of the
selected interacting generator. No full physical criterion follows.

Implementation: [interacting_records.py](../oph_openwave/interacting_records.py).
Checks: [test_interacting_records.py](../tests/test_interacting_records.py).
