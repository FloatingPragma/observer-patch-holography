# Original energy records and their conserved labels

The actual selected polynomial state suppresses negative-energy reverse
records strongly enough that every complete Walsh spatial row has positive
creation-energy probability greater than1/24. Its G1 and H parts exceed
5/384 and11/384. Complete total-cube and native-number detectors refine
these records without changing their original energy marginal.

[The proof](PROOF.md) also identifies explicit original G1g/Hg final-energy
branches with a small authenticated density floor. The actual original transfer measures have nonzero width from an exact
positive residual norm. The bound is small; it neither locates their lines
nor selects an isolated mass from the quantum labels. [The producer](model.py) retains this distinction.

```sh
python -m model.particle_primary_original_branch --check
python -m pytest model/particle_primary_original_branch/test_model.py -q -W error
```
