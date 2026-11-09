# A common physical observable limit

The selected Wilson source supports one subsequential bounded observable
C* algebra with a stationary state, a point-norm continuous time group,
and a nonnegative canonical GNS energy generator. The construction uses
the same finite source, physical Gauss state and clock at every level.
It includes gauge-invariant endpoint-charge projections and normalized
Wilson plaquettes at centered rational positions, with fixed triangular
time filters. One subsequence carries every finite star-word moment.
Complete two-outcome effect instruments built from these observables have
convergent fixed finite history probabilities and finite pointer density
matrices on that same subsequence, including every complementary outcome.

The [proof](PROOF.md) builds a separable algebra of bounded operator
sequences, quotients norm-null sequences, extracts its state, and applies
the existing finite downward-weight bound to every cyclic vector.
The existing Wilson cone gives exact limiting commutation outside its
conservative cone, enlarged by the temporal widths of both probes.
No exact state embeddings between finite levels are required.

The common limit can have only scalar or zero-energy responses. The
construction supplies no positive lower bound on a particle signal,
surviving interaction, Lorentz symmetry, particle masses or gravity.
The countable observable menu is not claimed to generate all physical
observables. An infinite fixed-charge menu requires charge tightness to
inherit completeness from the complete finite charge PVMs. The retained
idle-reference purification is not inserted into the stationary source
algebra.

Run the finite premise and adversarial checks with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_field_algebra
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_primary_field_algebra
```

`algebra.py` authenticates actual selected members and concrete probes,
computes rational time-word, negative-energy and filtered-cone bounds,
and executes a separate exact two-state matrix control. The source plans
do not allocate the full spatial Hamiltonian. The compactness theorem
does not give an effective subsequence or correlation convergence rate.
The finite pointer consequence is a marginal statement; no limit of the
growing canonical reference or efficient compiled apparatus is claimed.
