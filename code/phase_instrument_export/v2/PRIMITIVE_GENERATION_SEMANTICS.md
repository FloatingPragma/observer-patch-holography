# INS-03 v2 primitive generation semantics (proposed, unfrozen)

Semantics id: `oph.ins03.primitive_generation_semantics.v2_proposed`.
Transcript schema id: `oph.sim.ins03_primitive_transcript.v2_proposed`.
Status string: `PROPOSED_UNFROZEN`.

This document is a proposal. It specifies the grammar of primitive
operations for producing a phase-instrument export in the v1 schema
`oph.sim.ins03_phase_instrument_export.v1` from committed source objects, the
canonical machine state and its step digests, and the transcript format that
the executable replay verifier
`code/phase_instrument_export/v2/transcript_replay_verifier.py` checks. The
owner of instrument lane
[#737](https://github.com/FloatingPragma/observer-patch-holography/issues/737)
freezes this document before any producer transcript carries evidential
weight; until that freeze event exists, every transcript under this grammar
is a demonstration of the machinery and nothing else. The freeze discipline
is stated in section H, and the missing custody half is specified in
`AUTHENTICATED_BINDING_SPEC.md` in this directory.

Governing records: the deep audit
`plan/audits/RER_POST_R2020_ACTUAL_COMMITS_DEEP_AUDIT_2026-08-24.md`
(finding F3 and completion item 1), the plan section V3.26 of
`plan/COMPLETION_PLAN_V3.md` (the INS-03 passage), section C.2 (c) of
`plan/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md` (source-operation
origin), and the v1 contract
`code/phase_instrument_export/VALIDATOR_CONTRACT.md`. The audit's exact
judgments are load-bearing here and are repeated where they bind:

* `in_class = false` is not provenance. It states only that an operation is
  outside the enumerated constructor closure of
  `EventAlgebra.SourceReachability.Reachable`; it does not prove
  observer-patch origin.
* Replay verifies derivation, not origin. A verified replay proves that the
  claimed export is computable from the committed starting objects and the
  transcribed arguments under this grammar. It proves nothing about who or
  what produced the transcript, when, or from what physical or simulated
  process.
* Origin requires the authenticated binding of the companion specification:
  owner key custody, a signature over the binding tuple, and independent
  replay. A self-digest is not authenticated custody.

## A. Relation to the v1 checker

The v1 validator `code/phase_instrument_export/ins03_export_validator.py` is
a static committed-fixture checker: it can emit only
`STATIC_COMMITTED_FIXTURE_CONFORMANT`, it hard-pins the committed effect
table, run state, and count literals, and it rejects every producer claim
with `PRODUCER_AUTHENTICATION_UNIMPLEMENTED`. The v2 replay verifier is the
executable derivation half of the separate v2 gate that the register and
plan specify: frozen primitive-generation semantics, a primitive transcript,
machine proof or independent replay, authenticated artifact and
configuration binding, and target-blind deletion and mutation controls.

The division is deliberate:

* v1 binds one static document to the committed fixture literals. It cannot
  accept a lawful departure from the committed run state.
* v2 replay binds a transcript to this grammar and to the committed starting
  objects. It does not pin the derived export to the static fixture, so a
  future frozen semantics can derive lawful departures; the verifier still
  reads the claimed export through the v1 checker as an informational,
  non-gating cross-check.
* Neither half authenticates production. The authenticated binding of the
  companion specification is the only route to a producer verdict, and it is
  unimplemented.

The verifier reuses the exact arithmetic of the v1 module by import: the
rational field `Q`, the real quadratic field `Q(sqrt(3))`, the complexified
field `Q(sqrt(3), i)`, the matrix operations, the encodings, the canonical
serialization, and the committed reference constants are single-sourced in
`ins03_export_validator.py`, and the v2 tests check that the constants seen
through v2 are byte-identical to the v1 literals.

## B. Committed starting objects

Every derivation under this grammar starts from these committed objects and
from nothing else:

1. The committed eight-context effect table: contexts `web_diagonal`,
   `web_conjugated_0` through `web_conjugated_5`, `phase`, with the
   outcome-0 effect matrices of `ins03_export_validator.COMMITTED_EFFECT0`,
   transcribing `Lean/EventAlgebra/LuedersPhaseInstrument.lean`
   (`luedersPhaseInstrument_run_table`, `committedContextEffect_diagonal_eq`,
   `sourcePhaseLift_entries`). The outcome-1 effect is the exact complement
   `1 - E_0`.
2. The committed count literals `(111, 68, 179)`, `(315, 401, 716)`,
   `(179, 179, 358)` per context
   (`ins03_export_validator.COMMITTED_COUNTS`).
3. The committed run-state diagonal `diag(111/179, 68/179)` and the
   off-diagonal modulus bound `7548/32041`
   (`COMMITTED_RUN_STATE_DIAGONAL`, `COMMITTED_OFFDIAG_NORMSQ_BOUND`,
   specializing `EventAlgebra.prep_offdiag_normSq_le`).

The grammar contains no primitive that reads any other data source. A
producer needing operations outside this closure must extend the grammar,
and the extension binds only through a fresh freeze of this document.

## C. Exact arithmetic and encodings

All semantics are exact. No floating-point value enters any computation.

* An exact rational is encoded as the two-integer array
  `[numerator, denominator]` in lowest terms with positive denominator.
* An element of `Q(sqrt(3), i)` is encoded as the four-string list
  `[re_rational, re_sqrt3, im_rational, im_sqrt3]` of canonical rational
  literals (the `C3.encode` convention of
  `oph_fpe/quantum/phase_operation.py`, as transcribed by v1).
* A two-by-two matrix is encoded by the nested-list `encode_matrix` form.
* Canonical serialization of any JSON object is
  `json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False)`.
* Every digest is the lowercase hexadecimal SHA-256 of the UTF-8 bytes of a
  canonical serialization with no trailing line break.

## D. The replay machine state

The machine state is the septuple

```text
(effects, preparation, outcome_maps, summed_channel, readback,
 declarations, export)
```

with initial value `(None, None, {}, {}, {}, {}, None)`. Its canonical
encoding is the JSON object

```json
{
  "schema": "oph.ins03.replay_state.v2_proposed",
  "effects": null | {context: {"0": matrix, "1": matrix}},
  "preparation": null | {"rho_00": [n,d], "rho_11": [n,d],
                          "rho_01": [4 strings]},
  "outcome_maps": {...},
  "summed_channel": {...},
  "readback": {...},
  "declarations": {...},
  "export": null | {...}
}
```

where matrices use `encode_matrix` and the export-shaped blocks use exactly
the v1 export field shapes. The state digest after a step is the SHA-256 of
the canonical serialization of this encoding. Because serialization sorts
keys, the digest is independent of construction order and depends only on
the mathematical content of the state.

## E. The primitives

Six primitives exist. Each is deterministic and total on its stated domain;
an argument outside the domain is a named replay failure, never a silent
default. `where` in an error names the failing step as `steps[k]`.

### E.1 `load_committed_effect_table`

* Arguments: none (`{}`).
* Preconditions: `effects` is `None`; `export` is `None`.
* Semantics: sets `effects` to the committed eight-context table of section
  B item 1, outcome 0 from the committed literal and outcome 1 as the exact
  complement `1 - E_0`, in the committed context order.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on any argument;
  `STEP_PRECONDITION_UNMET` on a second load or after emission.

### E.2 `prepare_state`

* Arguments: `rho_00` and `rho_11` as exact rationals, `rho_01` as one
  exact `Q(sqrt(3), i)` scalar, in the encodings of section C.
* Preconditions: `preparation` is `None`; `export` is `None`.
* Semantics: checks exactly `rho_00 + rho_11 = 1`, `rho_00 >= 0`,
  `rho_11 >= 0`, and the positivity bound
  `|rho_01|^2 <= rho_00 * rho_11` in `Q(sqrt(3))`; on success stores the
  state coordinates. This is state preparation from transcribed literals:
  the grammar records what was prepared, and nothing in this primitive
  makes the literals source-produced.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on a malformed encoding, a trace
  different from one, a negative diagonal entry, or a positivity failure;
  `STEP_PRECONDITION_UNMET` on a second preparation or after emission.

### E.3 `assemble_lueders_kraus`

* Arguments: `context`, one committed context name.
* Preconditions: `effects` is loaded; the context is not yet assembled;
  `export` is `None`.
* Semantics: checks exactly that the outcome-0 effect is Hermitian and
  idempotent (`E = E^dagger`, `E * E = E`); the Lueders singleton assembly
  is defined for projector effects only. Builds the v1 field-1 through
  field-3 block for both outcomes with Kraus families `{E}` and `{1 - E}`,
  the induced effects `sum_k K_k^dagger K_k`, the exactly zero effect
  residuals, and the matrix-unit trace diagnostics, and the v1 field-4
  summed-channel block with the Kraus normalization identity and its
  spanning-set trace checks, all recomputed exactly.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on an unknown context name or a
  non-projector effect; `STEP_PRECONDITION_UNMET` on a missing table, a
  repeated context, or after emission.

### E.4 `readout_table`

* Arguments: `context`, integer `count_0`, integer `count_1`, integer
  `mass`.
* Preconditions: the context is assembled; `preparation` exists; the
  context has no readback yet; `export` is `None`.
* Semantics: checks exactly `mass > 0`, `0 <= count_i <= mass`, and
  `count_0 + count_1 = mass`. Computes the exact Born weight
  `Tr(rho E_{c,i})` from the prepared state and the assembled effects; the
  primitive is defined only where that weight is a plain rational, so that
  the exact compatibility residual
  `Tr(rho E_{c,i}) - count_i / mass` has the v1 rational encoding. Builds
  the v1 field-5 readback block. A nonzero residual is inside the grammar:
  the replay layer transcribes and recomputes it and imposes no tolerance;
  the `TOL_READBACK` band of design section C.4 is a freeze-time gate of
  the decision rule, not of replay.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on malformed integers, count bounds,
  a mass mismatch, or an irrational Born weight;
  `STEP_PRECONDITION_UNMET` on a missing assembly or preparation, a
  repeated readout, or after emission.

### E.5 `bind_declaration`

* Arguments: `field`, one of `carrier_id`, `source_record_ids`,
  `operations`, `provenance`, `labels`, `derived_for_display`; and `value`.
* Preconditions: the field is not yet bound; `export` is `None`.
* Semantics: shape-checks the value (`carrier_id` a nonempty string;
  `source_record_ids` a list of nonempty strings; `operations` a list of
  objects with string `operation`, `carrier_before`, `carrier_after` and
  boolean `in_class`; `provenance` an object whose `simulator_commit` is 40
  lowercase hex digits; `labels` an object with boolean `exploratory`,
  boolean `evidential`, and a nonempty `claim_boundary`;
  `derived_for_display` an object) and stores it verbatim.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on an unknown field or a shape
  failure; `STEP_PRECONDITION_UNMET` on a rebinding or after emission;
  `PRODUCER_AUTHENTICATION_UNIMPLEMENTED` when `labels.evidential` is
  `true`, because an evidential label is a custody claim and no
  authenticator exists.

Declared fields carry no derivational weight. They are freely stipulable
tokens: replay copies them and checks their shape, and their verification
is transcription, not computation. Every claim of source origin, custody,
or evidence routed through a declared field fails closed. The derived
fields (the outcome maps, summed channels, readbacks, state coordinates,
positivity certificate, record-diagonal certificate, and every digest) are
data whose verification is computation, and only they are checked by
recomputation.

### E.6 `emit_export`

* Arguments: none (`{}`).
* Preconditions: `effects` is loaded; `preparation` exists; every loaded
  context is assembled and has readback; `carrier_id`,
  `source_record_ids`, `operations`, `provenance`, and `labels` are bound;
  `export` is `None`.
* Semantics: assembles the complete v1 export object: schema id,
  `provenance_class` fixed to `synthetic` (a transcript cannot emit
  `producer`; that claim fails closed at the transcript level), the
  context list in committed order, carrier dimension two, the accumulated
  field-1 through field-5 blocks, the field-6 preparation block with the
  recomputed positivity certificate, record-diagonal certificate, and
  `preparation_content_sha256` over the canonical preparation body (the
  v1 convention: carrier id, producer commit, source record ids, and state
  coordinates, with no context field), the declared provenance and labels,
  the optional `derived_for_display` block, and the recomputed
  `custody_digest_sha256` under the v1 digest convention. Sets `export`.
  Any step after emission is outside the grammar, so a verified transcript
  ends with exactly one emission.
* Errors: `ARGUMENT_OUTSIDE_GRAMMAR` on any argument;
  `EMIT_PRECONDITION_UNMET` on a missing table, preparation, assembly,
  readout, or declaration; `STEP_PRECONDITION_UNMET` on a second emission.

## F. The transcript format

A transcript is one JSON file in canonical serialization with one trailing
newline, top-level fields exactly:

```json
{
  "schema": "oph.sim.ins03_primitive_transcript.v2_proposed",
  "semantics_id": "oph.ins03.primitive_generation_semantics.v2_proposed",
  "semantics_status": "PROPOSED_UNFROZEN",
  "transcript_class": "synthetic",
  "steps": [
    {"index": 0, "primitive": "...", "args": {...},
     "state_digest_sha256": "..."}
  ],
  "claimed_export": {...},
  "claimed_export_digest_sha256": "...",
  "transcript_digest_sha256": "..."
}
```

* `steps` is a nonempty ordered list; `index` equals the list position;
  `args` carries the exact arguments of section E;
  `state_digest_sha256` is the section D digest of the machine state after
  executing the step. Floats are outside the grammar in every argument
  except inside the `derived_for_display` value of `bind_declaration`.
* `claimed_export` is the complete export object the transcript claims to
  derive, including its own `custody_digest_sha256`.
* `claimed_export_digest_sha256` is the digest of the canonical
  serialization of `claimed_export`.
* `transcript_digest_sha256` is the digest of the canonical serialization
  of the transcript object with that one field removed. This self-digest
  is an integrity convenience for referencing the transcript by content.
  It is not custody: whoever edits the transcript recomputes it, so it
  binds the content to no agent, key, time, or execution, which is the
  audit's ruling on self-digests.
* `transcript_class` is `synthetic` or the reserved `producer`. The
  `producer` class always fails with
  `PRODUCER_AUTHENTICATION_UNIMPLEMENTED` until the authenticated binding
  of the companion specification is implemented and this document is
  frozen.
* One optional field `authenticated_binding` is reserved for the companion
  specification. Its presence invokes the named verification hook, which
  is unimplemented and fails closed.

## G. Replay verification and the verdict grammar

`transcript_replay_verifier.py` re-executes every primitive from the
committed starting objects in exact arithmetic, checking:

1. canonical serialization of the file text
   (`NONCANONICAL_SERIALIZATION`);
2. the schema, semantics id, status string, and transcript class
   (`TRANSCRIPT_SCHEMA_MISMATCH`, `SEMANTICS_ID_MISMATCH`,
   `SEMANTICS_STATUS_INVALID`, `TRANSCRIPT_CLASS_INVALID`); a
   `semantics_status` other than `PROPOSED_UNFROZEN`, including a claimed
   freeze, fails closed because no freeze registration exists to verify
   against;
3. the transcript and claimed-export digests
   (`TRANSCRIPT_DIGEST_MISMATCH`, `CLAIMED_EXPORT_DIGEST_MISMATCH`, with
   the format codes `TRANSCRIPT_DIGEST_FORMAT`,
   `CLAIMED_EXPORT_DIGEST_FORMAT`, `STEP_DIGEST_FORMAT`);
4. step structure and order (`STEP_FORMAT`, `STEP_INDEX_NONCONSECUTIVE`,
   `UNKNOWN_PRIMITIVE`);
5. each step by exact re-execution (`ARGUMENT_OUTSIDE_GRAMMAR`,
   `STEP_PRECONDITION_UNMET`, `EMIT_PRECONDITION_UNMET`), then the
   declared intermediate digest against the recomputed state digest
   (`STEP_DIGEST_MISMATCH` at the named step `steps[k]`);
6. the final replayed export against the claimed export byte-for-byte
   under canonical serialization (`FINAL_EXPORT_MISMATCH`; a transcript
   whose replay produces no export fails with `EMIT_ABSENT`);
7. every producer or custody claim, failing closed with
   `PRODUCER_AUTHENTICATION_UNIMPLEMENTED`: the `producer` transcript
   class, a `producer` provenance class in the claimed export, an
   `evidential` label, or an `authenticated_binding` block.

Verdict grammar, extending v1 exactly:

* `TRANSCRIPT_REPLAY_VERIFIED_UNAUTHENTICATED`: every check passes and the
  transcript class is `synthetic`. The verdict states a derivation fact
  and nothing else. It is explicitly not production, not provenance, and
  not custody; the report carries that boundary verbatim and the class
  marker `SYNTHETIC`.
* `NONCONFORMANT`: any check fails; the report lists every named error.
* `PRODUCER_AUTHENTICATION_UNIMPLEMENTED` is preserved from v1 as the
  fail-closed error code for every producer or custody claim. No producer
  success verdict exists in v2 replay.
* The v1 verdict `STATIC_COMMITTED_FIXTURE_CONFORMANT` is unchanged and
  appears in the v2 report only as the informational, non-gating
  `v1_static_verdict` reading of the claimed export.

## H. Freeze discipline

This document binds nothing until frozen. The freeze event, owned by lane
#737 under the register's artifact rule, consists of: the owner records
this document's SHA-256 and the verifier module's SHA-256 as freeze
artifacts in the instrument register; the status string of registered
transcripts moves from `PROPOSED_UNFROZEN` to a frozen status that names
the semantics digest; and the authenticated binding of the companion
specification is implemented so that a `producer` transcript can carry a
verifiable signature. Until all three exist, no transcript under this
grammar carries evidential weight, and the verifier accepts only the
`PROPOSED_UNFROZEN` status. A change to this document after freeze is a
new semantics id, never an in-place edit.

## What is not proved here

No producer exists, no owner key exists, no run exists, and no freeze
event exists. This document and its verifier certify, at most, that one
transcript's claimed export is derivable from the committed starting
objects through this proposed grammar. That is derivation, not origin:
nothing here establishes source production, provenance, custody, public
outcomes, or run reality, and `in_class = false` in any operations list
establishes only nonmembership in the enumerated Lean constructor closure.
The shipped `sample_synthetic_transcript.json` is marked synthetic and can
never present as production. Register rows PR-03, PR-64, and PR-65 are
open; PR-04 stays a declared row under its recorded disposition; OL-C5
stays `partial`. Nothing here arms an instrument, draws a seed, discharges
or demotes a premise, freezes a rule, scores a comparison, or supports any
physical claim of the OPH program.
