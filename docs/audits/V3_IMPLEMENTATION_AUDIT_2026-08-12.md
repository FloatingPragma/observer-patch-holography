# Completion plan V3 implementation audit

Audit window: 2026-08-12 through 2026-08-13. Baseline commit: `283e0b286a35d9dcb03dac52bb7e8b9178fc523e`. Audited Fable head: `1381ff871ee5a06edc962a649e62f1efe6e2f23b`. Repair commit: `3ab5bc2064235a740bb5574ea165564e43046bca`.

This was an adversarial implementation and epistemic-status audit, not a new scientific comparison. Three separated reviewer contexts inspected the plan, live issues, all post-baseline commits, tracking and claim registries, generated documentation, formal composition surfaces, validators, and the sibling OL-A1 simulator campaign. The machine-readable custody record is `AUD-V3-2026-08-12` in `tracking/audit_register.json`.

## Result

The V3 adequacy-first premise is sound, but the first implementation was not safe to accept unchanged. It mixed exact conditional mathematics, formal precursors, simulation evidence, and physical adequacy in several rows; omitted major physics domains; lacked a common-world and architecture-version contract; and allowed multiple registries to pass while incomplete or stale. Those defects were corrected without promoting any empirical prediction.

The observation ledger now separates four rungs: formal precursor, structural, emergent, and predictive. Only four rows remain `attained`, each within a narrow conditional structural boundary: OL-C1 (Born/Lueders under PR-02 and PR-03), OL-E2 (exact fluctuation identities under PR-07), OL-E3 (finite Landauer bound under PR-07 and PR-15), and OL-G3 (one-generation representation/anomaly package under PR-12). None establishes that the architecture selects its carrier or realizes the result in one inhabited physical world.

## Critical corrections

- The OL-A1 replication is recorded as `FAILED`. At 65,536 carriers the (1,3) signature reproduced in five of five replicates, but the 16,384 rung reproduced in only two of five; the registered ratio band missed in all five; and the ancestry-permutation null reproduced the verdict in fourteen of fifteen cells. The campaign cannot support an ancestry-emergence claim.
- The premise register was expanded to name physical time/dynamics, subsystem factorization, variational data, global-form attachment, chirality, Higgs/Yukawa structure, common-load and spacetime/photon/gauge-action bridges, and baryon-operator/proton attachments. Ledger rows now distinguish consumed premises from still-open premises.
- Spacetime, causality, no-signalling, Schrödinger dynamics, Maxwell propagation, family/chirality results, and several constants claims were narrowed or demoted where their cited theorems provide only formal or conditional precursors.
- The Standard Model table is a fixed nineteen-term correspondence, not an assembled Lagrangian. It loads the canonical premise register and exposes missing scalar, Yukawa, global-form, and physical-action attachments.
- Duplicate claim-registry keys, stale live-owner prose, obsolete issue URLs, incomplete surface inventory, and placeholder open-problem rows were repaired and given fail-closed gates.

## Completion-plan assessment

V3 is a materially better program than V2 because it permits a theorem of adequacy under an explicit premise vector before attempting premise discharge. That is the appropriate ordering for a reverse-engineering program. Three additions were necessary:

1. One common-world integration lane (#740), so separate conditional models cannot be counted as one realized universe.
2. Append-only architecture/protocol versions (#741), so changing the simulated substrate invalidates dependent promotions rather than silently changing their meaning.
3. Four missing physics lanes: cosmology/astrophysics (#742), interacting QFT/RG/scattering (#743), QCD through atomic physics (#744), and electroweak/flavor/neutrino physics (#745).

The phrase “P and N are the only numerical inputs” was also too strong. The corrected contract says P and N are the only *proposed fundamental free numerical parameters*; calibrations, measured comparisons, fitted anchors, and imported transport data remain numerical inputs with explicit ancestry.

## Simulator custody

The sibling `oph-physics-sim` campaign is pinned at commit `42aa96607f40dae94c5a8b65e9dd8e71e5b6434e`, including the exact manifest and fifteen receipts. A separate standard-library replay checked the exact inventory, all hashes and configurations, eigenvalue-derived fields, controls, and the recomputed FAILED verdict. The commit is reachable on its configured remote.

Raw feature matrices were not retained. Consequently this audit could not independently reconstruct observables from raw captures; it could only independently replay the retained summaries and arithmetic. Any future promotion-grade empirical campaign should preserve raw or fully reproducible source material and use a verifier independent of the producer.

## Residual scientific obligations

AV-0 is deliberately `exploratory_uninhabited`. It is a versioned conditional architecture contract, not a witness that any model satisfies every registered interface. Issue #740 owns that common-world witness. Issues #742 through #745 own substantial physics that the original fourteen-lane plan did not cover. Issue #739 remains the premise-discharge queue. Standing custody #738 must continue to prevent frozen or retrospective artifacts from being promoted without an eligible comparison contract.

No V3 issue required reopening: #726 through #739 were all still open at audit time. Audit amendments were added to those contracts, and the missing issues #740 through #745 were created instead of hiding their work inside unrelated lanes.
