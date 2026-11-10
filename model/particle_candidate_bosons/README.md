# Which boson detector can see the source excitation?

One complete loop apparatus on the specified primary source has nonzero
original-state records in six interacting gauge multiplets. An even scalar
or a vector detector cannot see the isolated A1u transition from the gauge
ground; an imaginary scalar loop does. The apparatus also sees distinct
scalar- and vector-type source channels. All probabilities retain the rare
empty-matter sector weight and the full reference.

The [proof](PROOF.md) gives the full-complement residue bound and exact
detector zeros. These cube labels do not identify photon, gluon, W, Z or
Higgs particles. Their physical field, vacuum, dispersion and mass
attachments are separate requirements.

```sh
python -m model.particle_candidate_bosons --write --check
python -m pytest -q model/particle_candidate_bosons
```
