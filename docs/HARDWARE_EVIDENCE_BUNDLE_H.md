# Evidence-Bundle Sufficiency for Hardware Claim Class H

Proof packet for issue #325. This document defines when a bundle of receipts
is sufficient evidence for a hardware claim, and when it is not. Hashes and
signatures prove identity and integrity relative to a trust root; they never
prove the physical truth of a claim. This packet defines the additional
structure that a physical claim requires.

## 1. Claim class H

A claim of class H asserts that a physical device produced a measured effect
under stated conditions. Formally, a class-H claim is a tuple

```text
H = (device D, protocol Pi, conditions C, effect statement E, magnitude M, uncertainty U)
```

where `E` names an observable, `M` a measured value or bound, and `U` a stated
error model. Examples in this repository: an energy-balance claim for a
resonator cell, a candidate-enrichment claim for an optical sampler, a lift or
effective-weight claim for a driven frame. Non-examples: simulation outputs,
design targets, and theory-side theorems; those live in other claim classes.

An OPH technology claim has an additional subject boundary. It concerns an
observer-like self-reading system: a bounded physical or software patch with
local state, ports or boundaries, readback, durable records, feedback or repair
moves, and a public evidence bundle. A generic optics, vibration, artificial
intelligence, mining, or engineering result does not become an OPH result by
using the same hardware.

## 2. Evidentiary predicates

A bundle `B` for a class-H claim consists of typed records. Each predicate
below is either satisfied by explicit artifacts in `B` or the bundle is
insufficient.

- **Raw capture.** Unprocessed instrument output for every channel named in
  `Pi`, captured before any analysis, with device identity and firmware state.
- **Calibration chain.** For every instrument, a calibration record tracing to
  a declared reference, dated inside the calibration validity window, with the
  transformation from raw units to reported units stated explicitly.
- **Custody.** A continuous record of who or what held the device and data
  between capture and publication, sufficient to exclude substitution.
- **Controls.** Sham runs, detuned twins, or blanks executed under the same
  protocol, interleaved with live runs, with their raw captures included.
- **Analysis binding.** The analysis code, its inputs, and the claim text are
  content-addressed, and the published `(E, M, U)` is reproducible from the
  bundle by an independent party running only bundle contents.
- **Completeness.** The bundle declares the full population of runs executed
  under `Pi`, including failures and aborts. A bundle that reports a favorable
  subset is insufficient regardless of the integrity of each member.

## 3. Threat model

Sufficiency is always relative to a threat model. The default model for this
repository assumes any of the following may occur and must be excluded or
detected by the bundle:

- **Signer compromise.** A valid signature by a compromised key. Mitigation:
  independent attestation (Section 4), never signature count.
- **Replay.** Re-presentation of an old valid bundle for a new claim.
  Mitigation: per-run nonces or timestamps bound into raw captures.
- **Selective reporting.** Publishing the runs that succeeded. Mitigation: the
  completeness predicate plus preregistered run schedules where feasible.
- **Device substitution.** Measuring a different device than the one named.
  Mitigation: custody records plus physical identity marks in raw captures.
- **Analysis degrees of freedom.** Post-hoc analysis choices that manufacture
  the effect. Mitigation: analysis code frozen before unblinding, recorded in
  the bundle.

## 4. Attestation rule

A class-H bundle is **sufficient** only when, in addition to every predicate
in Section 2, at least one of the following holds:

1. an independent party reproduced the effect from a fresh device or a fresh
   run series, with its own bundle; or
2. an independent party with no stake in the outcome witnessed the protocol
   end to end and attests to the custody and control records; or
3. the threat model for the specific claim is explicitly argued down (for
   example, a null result or an upper bound), and that argument is part of
   the bundle.

Claims of extraordinary physical effects (energy gain, lift, computational
advantage beyond classical baselines) always require rule 1 or rule 2.

## 5. Rejected counterexamples

The following bundles are integrity-valid and evidentially insufficient. Each
is rejected by a named predicate.

- A signed, hashed video of a device operating: fails raw capture and
  controls.
- A complete raw dataset with valid hashes whose calibration records are
  expired: fails the calibration chain.
- A perfectly bound analysis over five runs selected from forty: fails
  completeness.
- Two independent signatures over the same captured dataset: signatures do
  not compose into attestation; fails Section 4.
- A reproducible simulation matching the claimed effect: simulations are not
  class H; the bundle proves the model, not the device.

## 6. Executable v1 contract

The machine-readable scaffold consists of:

- `schemas/hardware_evidence_bundle_h_v1.schema.json`, the closed v1 type;
- `tools/verify_hardware_evidence_bundle_h.py`, an independent fail-closed
  verifier;
- `tools/hardware_evidence_external.py`, the operator-pinned trust,
  provenance, preregistration, analysis-replay, and attestation verifier;
- `code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/`, a
  tiny schema-valid fixture that explicitly says no physical device or
  measurement exists; and
- `code/audit/test_hardware_evidence_bundle_h.py`, the fast adversarial suite
  executed by the normal mandatory audit lane.

The verifier recomputes every artifact SHA-256, requires the canonical root to
cover every artifact, rejects unsafe paths, and cross-checks the bound run
schedule, successes, failures and aborts, raw run/device/nonce records,
calibration validity windows, control coverage, device identity, custody
intervals, analysis inputs, protocol and claim text. When replay is in scope,
an identified external nonce-registry snapshot is mandatory. Values in
`producer_assertions` are reported as ignored keys and never authorize a
predicate.

Run the reference fixture from the repository root:

```bash
python3 tools/verify_hardware_evidence_bundle_h.py \
  code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/bundle.json \
  --replay-registry \
  code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/replay_registry.json
```

The command returns `INSUFFICIENT` and exits nonzero. That is the required
result: this fixture tests the contract and makes no class-H claim. The suite
also rejects a seen replay nonce, omitted scheduled run, stale calibration,
device substitution, custody break, changed analysis or claim text,
compromised signer, and multiple signatures over one capture masquerading as
independent attestation. Another control gives an attacker every self-authored
declaration, lets the attacker rename itself as an independent party,
fabricate a Moon-levitation claim, and recompute all hashes. The packet cannot
promote without roots pinned by the verifier operator.

### Authenticated decision path

The producer bundle cannot nominate its own trust roots. The verifier operator
supplies two separate inputs:

1. a closed trust policy containing Ed25519 public keys, party and
   organization identities, roles, validity intervals, and revocation state;
2. an external-evidence packet containing signatures over the bundle root and
   artifact hashes, signed pre-run commitments, and a signed replay-registry
   snapshot.

The verifier checks actual Ed25519 signature bytes. It requires role-separated
provenance for raw captures, controls, calibration, device identity, custody,
analysis, protocol, claim text, and attestation. The run schedule and analysis
recipe must be signed by the preregistration authority before the first
capture. Those commitments bind the policy, bundle, claimant, device, and
protocol, so they cannot be transplanted to another claim. Every nonce must
have one signed pre-run assignment and post-capture consume receipt for the
same bundle and canonical root. The replay authority signs the complete
registry state, including prior seen nonces. An end-to-end witness must belong
to a different party and organization from the claimant and must sign both the
attestation artifact and the canonical bundle root.

Analysis replay uses a closed declarative recipe rather than executing
producer-supplied code. The implemented
`mean_max_deviation.v1` operation reads the bound raw samples and recomputes
the effect statement, unit, magnitude, and maximum-deviation uncertainty
exactly with rational arithmetic. A mismatch fails the analysis predicate.
Additional analysis operations require separately reviewed verifier code.

`test_authenticated_contract_has_a_nonvacuous_sufficient_path` constructs
independent ephemeral keys and a synthetic decision-procedure packet. The
verifier returns `SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL` and the CLI
exits zero. The fixture makes no claim about an actual device; it proves that
the policy is satisfiable rather than rejecting every input. Mutation tests
reject a bad signature, witness reuse after bundle rebinding, a witness
administered by the claimant's organization, a cross-bundle or post-capture
schedule commitment, unsigned replay-state drift, a nonce consumed for another
bundle, stale authority snapshots, malformed nested inputs, and a fully
re-signed raw mutation that disagrees with deterministic analysis replay.

The six named external gates remain machine-visible whenever their evidence
is absent or fails. They close per bundle when the operator-pinned trust,
independent witness, deterministic analysis, preregistration, provenance, and
replay-authority checks pass. `INVALID` has exit code `2`, `INSUFFICIENT` has
exit code `1`, and `SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL` has exit code
`0`.

## 7. Exact claim boundary

The test suite constructs a synthetic satisfiability fixture. The repository
ships no promoted physical claim. The issue #509 IBM bundle is an independently
replayable engineering specimen against its frozen controller nulls, but this
verifier does not bind or promote it, and its programmed circuit is
non-discriminating between OPH and standard quantum mechanics. The application
concepts in [APPLICATIONS.md](APPLICATIONS.md) are design documents.

A sufficient verifier result establishes that the named evidence predicates
hold relative to the operator's pinned trust policy and the declared threat
model. It does not prove the physical truth of the effect, protect against
collusion by every independent authority, or convert a generic hardware
result into an OPH result. OPH attribution additionally requires the
observer-like self-reading structure stated in Section 1.
