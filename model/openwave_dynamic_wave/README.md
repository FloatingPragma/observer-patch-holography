# Autonomous geometry and wave records

One12800-state source retains native8×8 matter, two quantum geometry modes,
all49 planar wave modes and the field vacuum. Exact harmonic extension on
the same7×7 graph determines geometry profiles, stiffness and inertia.
The geometry evolves during propagation, including field backreaction.

At source time4, the native-source transverse shift is certified in
[−0.001846753,−0.001846538] sites. A complete scheduled arrival instrument
gives delay[0.010034899,0.010075386], with hit and timeout probabilities
retained. Coherent native-clock/geometry/wave observations share the same
Hamiltonian. A cap2 comparison preserves both response signs.

The architecture is a declared scalar plane with two geometry coordinates,
κ=100 and a finite product oscillator cap. No physical G, observed lens
angle, Einstein closure or universal clock calibration is inferred.
`PROOF.md` gives the exact block reduction, numerical arithmetic contract,
complete instruments and scope.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_dynamic_wave --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_dynamic_wave --check
python3 -m unittest model.openwave_dynamic_wave.test_wave
```
