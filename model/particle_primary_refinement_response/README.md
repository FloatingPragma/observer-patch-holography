# Original finite responses across the primary capacity family

Every selected finite Wilson member has a populated neutral collective
response on its unchanged Gauss polynomial state and source clock.
The complete instrument retains the idle reference and its failure
outcome. Its original transfer weight is at least3/100 in
`[1/(100 beta),1446/a]`; the corresponding final-excitation weight is
at least1/50 in `[1/(100 beta),1446/a+100 epsilon]`. Complete finite
energy clocks retain at least1/50 in the displayed enlarged record band.

The proof combines a product of physical two-cell Ritz trials, an exact
256-state local CAR Schur certificate, and the polynomial passivity
lemma. The full-source beta-defect density is less than97/10 uniformly
over the selected family. A separate native-block retention gives a
positive q=6 matched-pair response at level one. These statements concern
averages over the source box, without a per-site uniformity assertion.

The energy band changes with capacity. No fixed finite-energy continuum
response, isolated branch, physical momentum or mass identification is
inferred. The trial state is used only for a full-ground upper bound;
the selected source density is unchanged. The large spatial density
matrix and clock histogram are not allocated by this certificate.

See [the proof](PROOF.md) for the source domains, finite-theorem consumers
and the separate noncollapse premises. Run:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_refinement_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/oph-completion-venv/bin/python -m model.particle_primary_refinement_response --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/oph-completion-venv/bin/python -m model.particle_primary_refinement_response --check
```

`bind_member(member)` and `verify_bound(member,receipt)` regenerate the
actual source, state, reference and clock. `certificate()['members']`
contains the receipts for levels zero, one and two; the proof of the
uniform bound covers every finite nonnegative level.
