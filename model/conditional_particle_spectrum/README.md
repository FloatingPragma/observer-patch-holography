# Conditional particle energies

One calibrated quadratic action produces finite energy records for all nine
charged fermions, three neutrino mass eigenfields, W, Z, Higgs, photon and eight
gluon components. Its native matter tensors, Higgs currents and canonical
kinetic terms give 180 fermionic modes and 56 physical bosonic modes at two
opposite momenta. Complete probes retain antiparticles, antineutrino flavour
analysers, zero outcomes and the bosonic no-operation record.

The GeV coefficients inherit empirical Yukawa, gauge and Higgs inputs; the
optional neutrino extension also assumes normal ordering and zero lightest
mass. The quark and charged-lepton values are running action parameters, and
the boson values are tree comparison coefficients. Physical interacting
particle masses require source, clock and pole matching. Separate composite
and native-source candidates are classified without adding them to this action.

The [proof](PROOF.md) specifies the finite domain and readout;
the [receipt](results.json) contains every mode and record.

```sh
python3 -m unittest model.conditional_particle_spectrum.test_spectrum
python3 -m model.conditional_particle_spectrum --write
python3 -m model.conditional_particle_spectrum --check
```
