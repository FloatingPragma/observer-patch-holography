# A complete finite ground and its detector records

The same finite scalar interaction has a unique ground on its entire
786,432-state carrier. An exact small-block certificate locates that state,
proves its gauge and spin singlet character, and gives complete scalar
records on one clock. Its energy is about −7.9463053 in the declared source
units. The 307-state sector ground at −7 is a prepared sector state rather
than the carrier vacuum. Relativistic particles require spatial propagation
and charged probes about this full ground.

The [proof](PROOF.md) explains the complete-carrier reduction, preparation
and detector law. The [code](ground.py) admits rational inequalities and
exact algebraic roots; floating calculations only propose checked witnesses.

```sh
python -m model.particle_scalar_ground
python -m pytest -q model/particle_scalar_ground/test_ground.py
```
