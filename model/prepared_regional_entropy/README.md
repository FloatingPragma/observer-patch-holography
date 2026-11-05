# Finite entropy changes of the prepared interacting field

The literal prepared Weyl/contact plaquette supports a two-sided family of
finite regional entropy variations, with all sixteen regional occupation
states included. Exact rational bounds give a positive relative-entropy
quadratic form and control the finite entropy error of the regional modular
polynomial over all 255 mathematical Hermitian trace-zero directions,
including the complete 127-dimensional parity-preserving subspace. The reference is
the actual subnormalized success block of the
[certified source preparation](../certified_reference/README.md).

For a Hermitian trace-zero regional change `X` with `||X|| <= 1/40`, the
changed density is positive and has the same success probability. With
`S(sigma)=-Tr(sigma log sigma)`, the exact identity and uniform bounds are

```
Delta S = Tr(X (-log sigma)) - D(sigma+X || sigma),
(1/2) ||X||_HS^2 <= D(sigma+X || sigma) <= 20 ||X||_HS^2,
|Delta S-Tr(X P8)| <= 2e-9 ||X||_1 + 20 ||X||_HS^2.
```

The abort complement is unchanged; success is never renormalized. The
[proof](PROOF.md) also gives a positive full-field marginal lift of every
such variation and a symmetric entropy derivative estimator with an
explicit error quadratic in the finite amplitude. No commutation with
the reference or particle-number restriction is needed.

This supplies a finite matter-entropy input for the interacting reference.
The variations and their full-field lifts are mathematical state families;
their record-resolved native preparation, a geometric modular charge,
area normalization and joint gravitational stationarity require separate
constructions.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.prepared_regional_entropy --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/prepared_regional_entropy
```

`--write` regenerates the exact receipt. Replay verifies the inherited
source hashes and recomputes the spectral and modular error bounds. It
does not repeat the native control certificate. Independent matrix tests
check noncommuting variations, the complete signed full-field lift,
finite entropy differences, and failure under unsupported normalization
or altered certificate inputs.
