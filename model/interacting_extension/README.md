# Interacting fields and complete detector records

The charged-lepton Higgs vertex acts with fermionic hopping and two native
quantum gauge links on complete finite physical sectors of dimensions 438 and
1596. Six compensated creation probes give a complete Ramsey experiment,
including blocked creation, both detector outcomes and their back-action.
The graph, scalar cutoff, coupling placement, prepared ground state and clock
are declared inputs. This construction supplies no physical mass prediction.

The [Higgs calculation](higgs_join/PROOF.md) verifies the full generated
complement and nonzero propagator correction. The
[record proof](PROOF.md) gives the measurement law and an exact finite-shot
indistinguishability bound. The [native reference calculation](source_selection/PROOF.md)
derives event-conditioned references and the correlated movement-count law.

The [spatial control](chiral_lattice/PROOF.md) supplies a free propagating
overlap branch and a separate interacting abelian flux ring. The
[electric calculation](chiral_electric/PROOF.md) verifies the connection and
positive geometric term for a changing chiral subspace. These are separate
finite controls; they do not establish one interacting chiral continuum or
an Einstein coarse-graining limit.

Source provenance is commit `3cf3c14810671a16e61b7915f69079a1773a12a2`
(r2062). The maintained implementation checks source bytes, including its
project import dependencies, against `source_manifest.json`. Descendant
checkouts are supported when the consumed source bytes agree.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/interacting_extension
OPENBLAS_NUM_THREADS=1 python -m model.interacting_extension.certificate --check
OPENBLAS_NUM_THREADS=1 python -m model.interacting_extension.higgs_join.certificate --check
OPENBLAS_NUM_THREADS=1 python -m model.interacting_extension.source_selection.certificate --check
OPENBLAS_NUM_THREADS=1 python -m model.interacting_extension.chiral_lattice.certificate --check
OPENBLAS_NUM_THREADS=1 python -m model.interacting_extension.chiral_electric.certificate --check
```

The full signed spectral measure, probability arrays, postmeasurement
vectors and 12,681,216 simulated bits can be exported with
`python -m model.interacting_extension.certificate --check --export /tmp/oph-records`.
The export is generated from the complete Hamiltonians, without a stored
eigenbasis or a measured-mass table. Precision-ladder and optional spatial
exports are documented in their owning packages.
