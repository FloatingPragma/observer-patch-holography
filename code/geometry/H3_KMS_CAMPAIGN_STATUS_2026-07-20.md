# H3/KMS campaign status

The frozen 4,096-carrier acceptance bundle,
`run_4k_acceptance_20260720_v5`, has a
valid replay instrument and exact disk replay. Its physical campaign preflight
returns `NO_GO`. It is not a physical H3/KMS emergence receipt.

The frozen status is:

```text
instrument_status = VALID_PASS
campaign_complete = false
cell_scientific_status = INCOMPLETE
physical_promotion_allowed = false
postrun_scientific_failures = []
replay_scientific_status = NOT_EVALUATED_BY_DISK_REPLAY
physical_preflight_verdict = NO_GO
physical_preflight_receipt = false
physical_stage_status = NOT_EVALUATED
```

The campaign has four declared rungs, 4,096, 16,384, 65,536, and 262,144,
and three declared seeds per rung. Only the primary 4,096-carrier cell has this
frozen replay bundle. One cell cannot complete the campaign.

The final v5 accounting is exact: all 21 top-level replay receipts pass, the
physical preflight lists 11 typed blockers, all nine physical stages are
`NOT_EVALUATED`, and the twelve-cell matrix contains one `INCOMPLETE` cell and
eleven `NOT_EVALUATED` cells. Promotion and branch retirement are both
forbidden.

Version v5 has a standalone numerical-runtime receipt that binds the Python implementation, ABI,
compiler, NumPy/SciPy stack, floating ABI, BLAS/SIMD metadata, byte order, and
five thread caps before source construction. A mismatched runtime is rejected
and a near-equal output has no replay status. The replay-manifest SHA-256 is
`3e05b72b326cf6266ff3c30880ea2c5dbe60289ee6585e0378a8a6b6f40782f4`;
the numerical-runtime receipt SHA-256 is
`ca84a0b44801ed1c9aa3ba4d6a73d9f499a1199ff62c3abcb540f25a3599b0e1`;
the source-capture semantic SHA-256 is
`b87f880350bf488b0b756899d54211402f8e9a2f6202a2ede81d5e4640b98f31`;
and the frozen campaign-plan SHA-256 is
`b8f3e870ad7b40266b8b7745a79bc5618d3787d2f728addfecb76a362c3342b7`.
These are replay-provenance receipts, not physical-emergence receipts.

The postcapture sensitivity diagnostic reported lower held-out loss for E4
than for S2, H3, or E3:

```text
E4 = 1.102824e-16
S2 = 0.0076832337
H3 = 0.0076971990
E3 = 0.0094118124
```

H3 therefore did not win this diagnostic comparison. This is not a physical
H3 failure: the feature maps were constructed after capture and are not
equal-footing outputs of an independent event producer. No clock candidate was
scored either, because the source lacks an independently produced modular-time
and geometric-flow pair.

Scientific evaluation is unavailable for seven source-typing reasons:

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

The diagnostic objects do not close those gates. The finite M4 state selects
the first sixteen presentation-ordered amplitudes and fails its source-selection
A5-invariance receipt. Its algebraic mixed-GNS construction repeats the same
constructed density matrix across support fibers; the paper's multiresolution
certificate and conditional-expectation receipt are false. The native BW
payload parses, but only BW02 passes its diagnostic predicate and the overall
BW receipt is false. The semantic event packet passes E1 population only; E2,
E3, E4, and the Lorentzian-cone receipt are false.

The compact source ledger is
[`runs/h3_kms_repaired_4k_status.json`](runs/h3_kms_repaired_4k_status.json).
It records the upstream receipt hashes and the fail-closed scientific status.
The full replay bundle is stored in the workspace campaign archive. It is not
duplicated here, and this compact ledger does not claim an independent replay.
