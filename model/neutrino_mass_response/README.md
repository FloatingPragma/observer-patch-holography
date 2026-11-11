# Neutral propagation and coupled sterile neutrinos

The native neutral Higgs-fermion composite propagates in an executed
three-site source. Its complete local ground response also gives an exact
selfenergy through eight retained complementary modes. Lepton number
forbids a Majorana term in that original law.

A separate unit-coupling extension adds a sterile right-handed neutrino
through its left-handed conjugate. Its complete lepton and sterile control
has a unique ground and neutral energy gaps 1 and 5, with all particle and
antiparticle records retained. A unit sterile Majorana interaction gives a
separately calculated ground and spectrum. An active tau antineutrino is a
weak-flavor label, distinct from the sterile coordinate.

The [proof](PROOF.md) states the source domains. The [receipt](results.json)
includes a conditional common-scale GeV example with its empirical input.
These finite gaps have no proved relativistic rest-mass attachment.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/neutrino_mass_response
OPENBLAS_NUM_THREADS=1 python -m model.neutrino_mass_response --write
OPENBLAS_NUM_THREADS=1 python -m model.neutrino_mass_response --check
```
