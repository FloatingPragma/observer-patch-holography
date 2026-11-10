# Charged energies with the transverse field retained

The fixed [mobile affine-Gauss source](../affine_gauss_mobile/PROOF.md) has a
70,928-dimensional one-charge sector. This package encloses its thirteen
lowest energy levels without discarding the transverse field. It certifies
a unique charged-sector ground, a positive gap above it, and a strict
energy reduction caused by field dressing. A complete compensated-charge
pump and energy instrument prepares that actual ground with bounded
success probability; both failure outcomes remain.

The actual neutral ground, Hamiltonian, reference and clock are those of
the parent. These are finite source energies in a problem with a fixed
boundary reference. Neither an isolated relativistic rest mass nor a
Standard Model species or GeV conversion is inferred.

```sh
python -m model.affine_gauss_spectrum --write
python -m model.affine_gauss_spectrum --check
python -m pytest -q -W error model/affine_gauss_spectrum
```

The [proof](PROOF.md) controls the entire photon complement analytically.
Outward interval LDL inertia certifies the small compression; no numerical
eigensolver or fitted oscillator tail establishes an enclosure.
