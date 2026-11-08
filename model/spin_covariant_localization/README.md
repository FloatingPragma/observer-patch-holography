# Covariant finite source localization

The continuum field of the selected spin ground has an exact covariant
position observable on its field-generated cyclic subspace. Its momentum is
the existing spatial translation generator, so the canonical commutator and
uncertainty inequality hold on an explicit common core. The
[proof](CONTINUUM_POSITION.md) derives the required spectral density from the
source field bounds. Identifying this generally nonlocal Gram-whitened
measurement with the native carrier position record requires an operational
attachment; the full vacuum-containing space requires a no-position outcome.

The same selected finite spin ground and local probe determine eight orthogonal
spatial-character modes. Their Fourier transform gives a complete covariant
coordinate detector with an explicit dark complement on the full source.
The detector probabilities and source-energy leakage are certified from the
complete finite witness. This gives a source-defined localization measurement;
its whitening is not proved local, and its probe subspace is not preserved
by the Hamiltonian. The [finite proof](PROOF.md) distinguishes these finite graph
coordinates from particle position, continuous momentum and relativistic mass.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.spin_covariant_localization --check
python -W error -m pytest -q model/spin_covariant_localization
```
