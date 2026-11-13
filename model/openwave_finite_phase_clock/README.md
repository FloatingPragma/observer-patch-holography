# Finite resources for a native particle clock

A finite controller moves through preparation and readout regions under
one time-independent Hamiltonian. All three native charged masses retain
their free source frequency under two preparation strengths. The 4096-state
calculation includes every controller position, reference record and energy
cost. An exact local commutator bounds contamination from the interaction
regions. The supplied moving controller and phaseful references distinguish
this finite experiment from a self-starting particle clock.
[Proof](PROOF.md), [source](source.py) and [complete records](results.json)
state the numerical and calibration domains.

Run `python3 -m model.openwave_finite_phase_clock --write`, followed by
`python3 -m model.openwave_finite_phase_clock --check` and
`python3 -m pytest -q model/openwave_finite_phase_clock`.
