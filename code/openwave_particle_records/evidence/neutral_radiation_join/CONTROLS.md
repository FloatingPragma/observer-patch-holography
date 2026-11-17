# Same-reader controls

All controls use the preregistered characteristic-field reader.

| Control | Required discrimination | What it tests |
| --- | --- | --- |
| No initial pair | No emitted field or remote inverse response | Field originates in contact, not the clock or reader |
| Local swaps disabled | Unit-speed characteristic update fails; no antipodal response | Ballistic relation depends on actual source transport |
| Detector disabled | Field norm remains2 and no pair is recreated | Reciprocal measurement is an actual operation |
| One current tag reversed | Same field norm, mismatched metadata, no inverse response | Norm alone does not supply the detector dictionary |
| Nine planted local amplitude pairs | Correct E/B for all pairs | Reader distinguishes signs and lane orientation |
| All single-slot impulse columns | Correct transfer and all finite phases | Frequencies come from the executed map |
| Native checkpoint write disabled | Conclusive symbol0 and changed nonzero field | Native readout depends on the actual write |
|32 holding coins on a write edge | Explicit write timeout | A failed store is not decoded as a success |
| Four all-holding read scans | Explicit seed timeout with no label | An inconclusive census gives no guessed symbol |

The transport-off control is expected to fail the unit-speed transfer law;
its failure is the required discrimination, not a suppressed bad source trial.
The detector-off and tag-scramble controls preserve the dimensionless field
norm but fail the original pair response. All their events are retained.

Attempt01's recorder alias and serialization error are implementation failures,
not controls or evidence. They remain in the packet's execution history.
