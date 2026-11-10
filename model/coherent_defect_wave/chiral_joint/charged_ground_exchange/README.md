# Charged responses and spin exchange

The atomic ground spins interact through virtual hopping in the complete
quantum-link source. This calculation evaluates their second-order exchange
coefficient with its actual charged energy denominator. Separate lepton-spin
conservation makes the quark spin response vanish at every intermediate time.
The two lepton responses retain the quantum scalar and the full maximal quark
module; their exact algebraic resolvents determine the exchange.

The [proof](PROOF.md) and [machine receipt](results.json) retain the complete
source and the scalar energy shift. The exchange enters the controlled
small-hopping expansion; the calculation supplies no stationary-ground
transfer to unit hopping or particle-mass identification.

Run `python -m model.coherent_defect_wave.chiral_joint.charged_ground_exchange --check`
for an independent exact replay. The mutually exclusive `--write` mode
generates the receipt.
