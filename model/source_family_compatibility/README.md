# A common clock for spin and fermionic source sectors

The full-ground strong-coupling theorem imposes an exact compatibility
constraint on particle reconstruction. Throughout its stated domain,
keeping the ordered spin field's ballistic clock nonzero sends every
occupation-charged source excitation to infinite energy. A chemical
reference shift cannot retain both conjugate charges on that clock.
This applies to the specified projection family, not to every OPH source.

The [proof](PROOF.md) gives the complete energy-weight bound, source grading,
an exact finite countercontrol with both spin and odd CAR excitations,
and the precise finite gauge, chirality and statistics attachments.
It also replays an existing gauge-invariant experiment whose two admitted
operation parameters give different complete records. No new Hamiltonian
or fitted mass is introduced.

```sh
python -W error -m pytest -q model/source_family_compatibility
python -W error -m model.source_family_compatibility --check model/source_family_compatibility/results.json
```

The viable next construction must retain the required field sectors on
one source, state and clock. It can leave the certified spin-projection
domain while keeping the full CAR source, or construct different emergent
sectors with their own justified grading. Neither option is supplied by
renaming the spin fields or an internal representation.
