# H3/KMS campaign status

The repaired 4,096-patch cell is a valid instrument run with exact disk replay.
It is not a physical H3/KMS emergence receipt.

The frozen status is:

```text
instrument_status = VALID_PASS
campaign_complete = false
cell_scientific_status = INCOMPLETE
physical_promotion_allowed = false
postrun_scientific_failures = []
replay_scientific_status = NOT_EVALUATED_BY_DISK_REPLAY
```

The campaign has four declared rungs, 4,096, 16,384, 65,536, and 262,144,
and three declared seeds per rung. Only the repaired 4,096-patch cell has this
frozen replay bundle. One cell cannot complete the campaign.

Scientific evaluation remains unavailable for seven source-typing reasons:

1. The nested icosahedral support tower is an independent regulator. Source
   repair normal forms do not produce it.
2. The finite M4 density and modular generator are constructed from one source
   snapshot without an independent physical state producer.
3. The cross-ratio coordinate is target-blind after capture, but it lacks an
   independently normalized geometric-flow observable.
4. The BW payload is a structural postcapture reconstruction without a
   cofinal physical state tower and independent modular/geometric pair.
5. Source capture lacks independent modular-transport time and geometric-flow
   parameter pairs.
6. The H3, S2, E3, and E4 features are model-specific sensitivity diagnostics,
   not equal-footing physical controls.
7. Event positions, boxes, cones, and frame coordinates are postcapture
   diagnostics without an independent event-manifold producer.

The compact source ledger is
[`runs/h3_kms_repaired_4k_status.json`](runs/h3_kms_repaired_4k_status.json).
It records the upstream receipt hashes and the fail-closed scientific status.
The full replay bundle remains in the workspace campaign archive; it is not
duplicated here because its event report is about 13 MB and this ledger does
not claim an independent replay.
