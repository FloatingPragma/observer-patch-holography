# Evidence-Bundle Sufficiency for Hardware Claim Class H

Proof packet for issue #325. This document defines when a receipt bundle is
sufficient evidence for a hardware claim. Hashes and signatures establish
identity and integrity relative to named trust roots. Physical truth requires
the additional predicates and trust assumptions stated here.

## 1. Claim class H

A claim of class H asserts that a physical device produced a measured effect
under stated conditions. Formally, a class-H claim is a tuple

```text
H = (device D, protocol Pi, conditions C, effect statement E, magnitude M, uncertainty U)
```

`E` names an observable, `M` gives a measured value or bound, and `U` states
the error model. Examples include an energy-balance claim for a resonator cell,
a candidate-enrichment claim for an optical sampler, and a lift or
effective-weight claim for a driven frame. Simulation outputs, design targets,
and theory-side theorems belong to other claim classes.

An OPH technology claim has an additional subject boundary. It concerns an
observer-like self-reading system: a bounded physical or software patch with
local state, ports or boundaries, readback, durable records, feedback or repair
moves, and a public evidence bundle. Shared hardware alone does not turn a
generic optics, vibration, artificial-intelligence, mining, or engineering
result into an OPH result.

## 2. Evidentiary predicates

A bundle `B` for a class-H claim consists of typed records. Every predicate
below must pass.

- **Raw capture.** The executable v1 protocol names one measurement channel,
  its raw unit, the reported unit, the device, and a firmware hash. Each
  scheduled run carries exactly one raw artifact with a nonempty sample
  population. The raw artifact binds the run, device, protocol, capture nonce,
  capture time, sequence index, channel, raw unit, firmware hash, and hashed
  physical serial mark. A protocol with several channels requires a separately
  reviewed operation and schema version.
- **Calibration chain.** Each calibration names the channel and units, its
  validity interval, a reference identifier, a reference URI, and the SHA-256
  digest of the included reference-certificate artifact. The
  raw-to-reported transformation is an explicit affine map with rational scale
  and offset. The verifier applies this map to every raw and control sample,
  and the validity window must cover both captures. The declared calibration
  identifiers must equal the complete calibration population used by the run
  schedule; unused calibration rows are rejected.
- **Custody.** Continuous time segments cover capture through publication for
  the named device. The signed custody record also lists the complete bound
  data population, which must equal every bundle artifact other than the
  custody record itself.
- **Controls.** Every scheduled run has exactly one blank, sham, or detuned
  control. Its sequence index immediately precedes the live capture. Its
  timestamp falls after the preceding live capture, if one exists, and before
  its paired live capture. Device mark, firmware, protocol, channel, and raw
  unit must agree across the pair. Raw and control captures carry distinct
  registry-managed nonces.
- **Analysis binding.** The protocol, analysis recipe, inputs, and exact
  structured claim are content-addressed. The declared input population must
  equal the bound evidence population. The run schedule, protocol, and analysis
  recipe carry signed pre-run commitments from the pinned preregistration
  authority. The verifier runs a closed declarative operation;
  producer-supplied code and producer-selected replay commands cannot promote a
  bundle.
- **Completeness.** The ordered preregistered run schedule equals the ordered
  reported population, including successful, failed, and aborted runs. The
  declared raw and control artifacts equal the artifacts consumed by those
  runs and by the analysis.

## 3. Threat model

Sufficiency is relative to a threat model. The closed v1 schema fixes the
following modes as in scope; a producer cannot disable one in its bundle:

- **Signer compromise.** A compromised key can produce a valid signature. A
  physical promotion therefore requires one independently administered
  witness key to sign the closed attestation, the canonical bundle root, and
  the replay-registry digest.
  The witness party and organization must differ from the claimant and every
  measurement, calibration, device, custody, preregistration, and replay
  authority.
- **Replay.** An old valid packet can be presented under a different claim.
  A signed pre-run reservation artifact lists every ordered raw and control
  capture nonce. Each nonce has one post-capture consume receipt bound to the
  bundle identifier and canonical root. The replay authority signs the complete
  registry snapshot, and the independent witness signs its exact file digest.
- **Selective reporting.** A favorable subset can hide failed or aborted
  runs. The preregistered schedule, the typed run population, the raw capture
  population, and the control population must agree exactly.
- **Device substitution.** A different device can be measured under the
  declared identifier. Device-authority records, raw captures, controls,
  calibration, custody, firmware, and hashed physical marks must agree, and
  the independent witness covers device identity and capture.
- **Analysis degrees of freedom.** Post-capture choices can manufacture an
  effect. The protocol and analysis recipe are committed before capture. A
  closed verifier operation recomputes the result from every paired sample.

The model assumes the operator selected the trust policy independently of the
producer. It also assumes that at least one required independent authority
remains outside a coalition. Cryptographic receipts cannot detect a coalition
that controls every pinned authority and signs a mutually consistent fiction.

## 4. Attestation rule

A class-H bundle is sufficient only when every predicate in Section 2 passes
and one of these evidence routes is resolved:

1. an independent party reproduces the effect from a fresh device or run
   series and supplies a separately verified bundle;
2. an independent party with no stake in the outcome witnesses the protocol
   end to end and attests to the schedule, identity, capture, calibration,
   controls, custody, nonce reservation, replay registry, analysis freeze, and
   claim; or
3. the claim-specific threat model excludes reproduction and witnessing, as
   can be appropriate for a null result or upper bound, and the exclusion
   argument is itself bound and reviewed.

Extraordinary-effect claims, including energy gain, lift, and computational
advantage beyond classical baselines, require route 1 or route 2. The
executable v1 sufficient path implements route 2. A declared reproduction
returns `INSUFFICIENT` until a verifier resolves the separate fresh-run bundle.
The executable contract does not promote route 3.

## 5. Rejected counterexamples

These packets can have valid internal hashes while remaining evidentially
insufficient:

- a signed video of a device operating, which lacks raw capture and controls;
- a complete raw dataset with an expired or transform-free calibration;
- five favorable runs selected from a preregistered population of forty;
- two signatures over the same captured dataset, which do not constitute an
  independent witness or reproduction;
- a raw mutation with recomputed hashes and trusted signatures whose
  deterministic replay disagrees with the structured claim;
- a witness from the measurement organization, or separate witnesses that
  divide the attestation-artifact and root signatures between two keys;
- a protocol edited after capture while retaining its pre-run commitment; and
- a reproducible simulation, which is outside claim class H.

The adversarial suite also constructs a Moon-levitation claim whose producer
controls every bundle field and recomputes every internal hash. It fails
because producer-authored keys and independence declarations are not
operator-pinned evidence.

## 6. Executable v1 contract

The machine-readable packet consists of:

- `schemas/hardware_evidence_bundle_h_v1.schema.json`, the closed bundle type;
- `tools/verify_hardware_evidence_bundle_h.py`, the fail-closed bundle
  verifier;
- `tools/hardware_evidence_external.py`, the operator-pinned trust,
  provenance, preregistration, deterministic-analysis, replay, and attestation
  verifier;
- `code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/`, a
  schema-valid synthetic fixture with closed JSON protocol, analysis, and
  structured-claim artifacts; and
- `code/audit/test_hardware_evidence_bundle_h.py`, the adversarial test suite
  executed by the mandatory audit lane.

The verifier reads and hashes every artifact, requires the canonical root to
cover every artifact, rejects unsafe paths, and treats `producer_assertions` as
untrusted audit notes. It checks ordered schedule equality, paired captures,
calibration transforms, custody of the complete data population, exact claim
equality, protocol and analysis preregistration, role-separated Ed25519
signatures, witness independence, and authoritative nonce consumption.

Run the reference fixture from the repository root:

```bash
python3 tools/verify_hardware_evidence_bundle_h.py \
  code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/bundle.json \
  --replay-registry \
  code/audit/fixtures/hardware_evidence_bundle_h/reference_nonphysical/replay_registry.json
```

The command returns `INSUFFICIENT` with exit code `1`. The fixture declares
that no device or physical measurement exists, and it supplies no
operator-pinned trust packet.

### Authenticated decision path

The producer bundle cannot nominate its trust roots. The verifier operator
supplies a closed Ed25519 trust policy and a separate external-evidence packet.
The policy binds public keys to party, organization, role, validity interval,
and revocation state. Every anchor must contain distinct Ed25519 public-key
material, so separate identity labels cannot disguise one signing key as both
an evidence authority and an independent witness. Party and organization
identifiers are nonempty, and a declared witness signer must name the same
party as its pinned anchor. The evidence packet carries signatures, pre-run
commitments, and the signed replay-registry snapshot.

The claimant and independent witness sign the canonical root. Authorities sign
the artifacts assigned to their roles. The preregistration authority signs the
ordered run schedule, protocol, nonce reservation, and analysis recipe before
the first control capture. Each post-capture nonce assignment binds the
policy-independent registry identifier, bundle identifier, and canonical root.
One witness key must sign the closed end-to-end attestation, the root, and the
replay-registry digest. That witness must be separate from every
evidence-producing authority. Its explicit witnessed scope covers the protocol
as well as schedule, device, captures, calibration, controls, custody,
analysis freeze, claim, nonce reservation, and replay registry. The supplied
registry snapshot is a closed object; unsigned extension fields are rejected.

The closed analysis operation takes paired raw and control samples
`(x_k, c_k)`. For the declared affine calibration `T(z) = a z + b`, it computes

```text
d_k = T(x_k) - T(c_k)
M   = mean_k(d_k)
U   = max_k |d_k - M|
```

All arithmetic uses exact rational numbers. The resulting effect statement,
unit, `M`, and `U` must equal the structured claim. The claim artifact must
equal the entire structured claim, including device, protocol, conditions,
effect, magnitude, uncertainty, and extraordinary-effect flag.

The focused suite constructs independent ephemeral keys and a synthetic
policy packet, then obtains
`SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL`. Negative cases cover malformed
schema input, unsafe evidence substitutions, missing trust roots, bad
signatures, compromised signers, same-organization witnessing, split witness
signatures, post-capture or cross-bundle commitments, protocol mutation,
selective reporting, replay-state drift, nonce rebinding, stale snapshots,
calibration windows and reference-certificate drift, device and firmware
substitution, custody gaps at the control or live capture boundary, unpaired
controls, raw or control nonce reuse, extra analysis inputs, altered claim
text, producer-selected replay code, attempts to disable a v1 threat, and
nonce-reservation omissions. It also covers analysis disagreement after fully
re-signed raw, control, or calibration mutations.

The six external gate codes remain visible when their evidence is absent or
fails. `INVALID` has exit code `2`, `INSUFFICIENT` has exit code `1`, and
`SUFFICIENT_RELATIVE_TO_DECLARED_THREAT_MODEL` has exit code `0`.

## 7. Acceptance map

| Required property | Fail-closed enforcement |
|---|---|
| Claim class H and evidentiary predicates | Closed schema, typed artifact formats, predicate report, and exact structured claim |
| Trust roots, compromise, replay, selective reporting, substitution | Operator-supplied policy, role signatures, preregistered nonce population, independently co-signed nonce registry, ordered population checks, and device/firmware/mark agreement |
| Raw, calibration, analysis, and claim binding | Artifact hashes, canonical root, signed provenance, explicit affine transform, closed replay, and exact claim equality |
| Reproduction or attestation where required | Physical promotion requires the resolved witness path; reproduction remains insufficient until its separate bundle is verified |
| Integrity-valid physical falsehood controls | Self-authored false claims, claim/data disagreements, protocol reuse, witness conflicts, and omitted predicates fail |

## 8. Exact claim boundary

The repository ships no promoted physical class-H claim through this verifier.
The positive test establishes that the evidence-policy decision procedure is
satisfiable. The reference fixture explicitly denies any physical experiment.

A sufficient verdict states that the named evidence predicates hold relative
to the operator's pinned trust policy and the declared threat model. It does
not establish the physical truth of the effect, protect against a coalition
controlling every independent authority, or convert a generic hardware result
into an OPH result. OPH attribution also requires the observer-like
self-reading structure stated in Section 1.

The issue-509 IBM bundle is an independently replayable engineering specimen
against its frozen controller nulls. This verifier does not bind or promote
that bundle, and its programmed circuit is non-discriminating between OPH and
standard quantum mechanics. The application concepts in
[APPLICATIONS.md](APPLICATIONS.md) are design documents.
