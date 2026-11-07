# Matching a field cone to observer directions

A field's propagation geometry and an observer's causal geometry must agree
in the same coordinates and clock. The exact native twelve-direction decoder
reduces this matching problem to twelve scalar tests. The criterion includes
an explicit positive example and is covariant under a common change of frame.
It is conditional on the field tensor, spatial readback and clock conversion.

The [proof](PROOF.md) distinguishes this principal-cone test from equality of
complete field observations. A two-radius control explains why the existing
unwarped radial product cannot be identified with flat polar cones through
the fixed map `x=r*n` and one global clock scaling. It does not exclude that
product's curved Lorentzian description or other coordinate dictionaries.
No source action is changed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/field_cone_compatibility
python -m model.field_cone_compatibility --check
```

`cone.certificate()` replays the native Gram and tensor decoder and returns
exact examples. The proof's general real-matrix statements are proved in
prose; the executable works over exact `Q(sqrt(5))` data.
