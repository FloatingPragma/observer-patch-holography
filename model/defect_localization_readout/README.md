# A created local spin excitation protected through repair

A two-pulse local spin operation on the prepared source cube creates an
opposite spin-parity sector while leaving its initial exterior state unchanged.
The same source's complete detector distinguishes it from the reference.
Every repair branch preserves this symmetry label. A bound on the full
recorded dynamics keeps a local detector signal above `2/5` through ten
repairs up to time `1/16` on the specified strong-coupling trajectory.

The [proof](PROOF.md) distinguishes this finite symmetry label from a
topological defect, a permanently localized state and a particle mass.
The parity is conserved by the chosen Hamiltonians and every C/E repair;
coherent creation controls can change it. It is not an invariant under all
admitted observer interventions or an asserted infinite-volume superselection
charge. The bare occupation records cannot read the internal spin.

The two complete programs prepare the actual eight-site selected ground,
apply either no creation or the local pulse, then read X. They execute
1,090 and 1,092 source pulses on all 131,072 work rows. All first-failure
slots are retained. The ten-repair lifetime is an analytic all-record
bound, not an executed long tape. Preparation, coherent controls and the
pointer consume the existing declared service interfaces.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_localization_readout --check
python -m pytest -q -W error model/defect_localization_readout
```

The receipt replays both finite programs and checks the actual Hamiltonian,
complete raw branches and certified ground proposal. Its physical mass table
is empty.
