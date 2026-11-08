# Source energy, action and repair effort

Energy is measured from the specified source clock, using the same selected
eight-spin cube state as `protected_defect_catalog`. Complete native charge
creation and energy measurements give positive finite excitation energies.
The [energy table](ENERGY_READOUT.md) retains both charge weights and their
unconditioned energy moments. These are conditional mean energies in source
units J=1; they identify no particle masses or laboratory species.

The construction supplies two independent readouts. A 64-setting Pauli
instrument measures the first moment of the actual H_B Hamiltonian. A
controlled H_B clock measures its phase derivative, with the common ground
phase included. Both start from blank, prepare the same selected ground,
create the coherent local charge channels and retain every native failure.

There is a precise positive relation between energy and signed action:
under this common clock, relative action accumulated for duration T equals
T times the excitation energy. The native port action is computed from all
twelve source generators. It is state dependent and signed. Under H_B
evolution, energy is constant and ongoing system power is zero. Bare C/E
repairs act identically on every protected spin state, so inserting them
requires no system-energy input. The controlled clock's implementation cost
is separate; unread C/E repair is not identified with H_B evolution.

Repair records alone cannot determine this energy. The complete protected
repair histories are identical for source states of energies 0 and -6, and
for the two created charge channels with different mean energies. The
proof gives the exact measurement sufficiency criterion and constructs a
measurement that satisfies it. A reversible native control cycle also has
positive pulse cost and zero ideal energy work. Neither pulse count nor
record entropy is a universal energy functional.

The actual regional source state supplies a modular first-law check without
a Gibbs reference: local Z creation changes physical energy and modular
relative entropy while preserving regional entropy exactly. A reservoir
temperature, apparatus energy and memory-reset heat are separate premises.
The [proof](PROOF.md) states these boundaries and the complete error budgets.

Run `python -m model.source_energy_effort --check
model/source_energy_effort/results.json` and
`python -W error -m pytest -q model/source_energy_effort`.
