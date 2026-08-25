# INS-03 v1 static committed-fixture checker contract

Schema under validation: `oph.sim.ins03_phase_instrument_export.v1`, the
simulator-side export interface specified in section A of
`plan/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md` (OPH meta planning
surface, outside this repository).  Validator:
`code/phase_instrument_export/ins03_export_validator.py`.  Owning lanes:
issues 730 (OL-C5, quantum) and 737 (instrument register).  The typed Lean
binding interface `Lean/EventAlgebra/SourceBoundInstrumentInterface.lean` is
the typed boundary inside the corpus. The v1 checker is not the receiving half
of a producer handshake: it authenticates no external source, run, provenance,
or custody data. It checks only a synthetic static transcription of the
committed fixture. A finite produced run requires a separately frozen v2
schema, decision rule, and authenticator.

## What the validator certifies

A `STATIC_COMMITTED_FIXTURE_CONFORMANT` verdict certifies, for one synthetic
JSON document, all of the following, checked fail-closed with named error
codes and exact arithmetic (`fractions.Fraction` plus a minimal exact
`Q(sqrt(3))` implementation; no floating-point value enters any check):

1. **Schema id and version.** The `schema` field equals
   `oph.sim.ins03_phase_instrument_export.v1`; the version is the `.v1`
   suffix of that string.
2. **Canonical serialization.** The file text equals
   `json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False)` with at
   most one trailing newline; no absolute path, no timestamp, and no float
   appears outside a `derived_for_display` block.
3. **Exact encodings.** Every exact rational is a two-integer array
   `[numerator, denominator]` in lowest terms with positive denominator
   (`oph.sim.qm_observer_viz.v1` convention, `oph_fpe/qm_observer/DESIGN.md`
   section 7).  Every scalar of `Q(sqrt(3), i)` is the four-string list
   `[re_rational, re_sqrt3, im_rational, im_sqrt3]` written by `C3.encode`
   of `oph_fpe/quantum/phase_operation.py`, each component a canonical
   rational literal.  A matrix of dimension two uses the nested-list
   `encode_matrix` form; any other dimension uses the `encode_matrix_n`
   extension, the same per-entry encoding in row-major order with the
   integer `n` carried beside the matrix.
4. **Required fields per context and outcome**, transcribing design
   section A: `outcome_maps[c][i]` with the ordered Kraus family (field 1),
   `effect_from_kraus`, `declared_effect`, and the exact `effect_residual`
   (field 2), and the legacy-named `trace_nonincreasing` diagnostic values on
   the frozen spanning set of matrix units (field 3); `summed_channel[c]`
   with the Kraus normalization, its residual, and the spanning-set trace checks
   (field 4); `readback[c][i]` with the public outcome symbol, integer
   count, context mass, and exact compatibility residual (field 5);
   `preparation` in the two-by-two coordinatization with the positivity
   certificate, the record-diagonal certificate, the carrier id, the source
   record ids, the `operations` list flagged against the constructor list
   of `EventAlgebra.SourceReachability.Reachable`, and
   `preparation_content_sha256` (field 6); `provenance` with producer
   module hashes, commits, repository URL, run id, input inventory, runtime
   read log, and the import-graph independence report (field 7); `labels`
   with the exploratory and evidential booleans and the verbatim claim
   boundary (field 8).
5. **Dimensional consistency.** Every matrix matches the declared
   `carrier_dimension`; the preparation coordinatization requires dimension
   two, as design field 6 states.
6. **Exact algebraic conformance**, computable from the export alone:
   * each Kraus family's induced effect `sum_k K_k^dagger K_k` equals the
     declared `effect_from_kraus`, which equals `declared_effect`, with the
     residual matrix exactly zero, and the induced effect passes the exact
     Hermitian two-by-two positive-semidefinite test;
   * for each outcome, `1 - E_i` equals the other outcome's induced effect
     and passes the same exact positive-semidefinite test. This is the
     outcome trace-nonincreasing certificate: Kraus form gives `E_i >= 0`,
     while `1 - E_i >= 0` gives `E_i <= 1` and hence
     `Tr(Phi_i(X)) <= Tr(X)` for every positive-semidefinite input `X`;
   * the declared spanning-set trace values equal the exact recomputation
     through the Kraus families, per outcome and for the summed channel.
     These matrix-unit values are diagnostics of the linear trace action,
     not positivity tests, and their differences are not required to vanish;
   * the per-context effects sum exactly to the identity, and the
     summed-channel Kraus normalization equals the identity with residual
     exactly zero;
   * the declared run-state outcome traces `Tr(rho E_{c,i})`, recomputed
     from the exported preparation and effects, equal the declared
     frequencies `count_i / mass` exactly, and each declared compatibility
     residual equals the exact recomputation;
   * the preparation trace is one, the positivity bound
     `|rho_01|^2 <= rho_00 * rho_11` holds exactly in `Q(sqrt(3))`
     arithmetic, and the record-diagonal certificate matches the exported
     off-diagonal coordinate.
7. **The committed-table cross-check.** Where the export declares the
   committed eight contexts (`web_diagonal`, `web_conjugated_0` through
   `web_conjugated_5`, `phase`), its effect table, run state, and count
   table must equal the committed literals, hard-coded in the validator as
   reference constants: the effect table and run-state diagonal
   `diag(111/179, 68/179)` from `Lean/EventAlgebra/LuedersPhaseInstrument.lean`
   (`luedersPhaseInstrument_run_table`, `committedRunState_eq_literal`,
   `committedContextEffect_diagonal_eq`, `sourcePhaseLift_entries`) and
   `Lean/EventAlgebra/SourceBoundInstrumentInterface.lean`
   (`SourceBoundDeterminedData.publicTable_literal`,
   `modelFrequency_phase_zero`); the count literals `(111, 68)` at mass
   179, `(315, 401)` at mass 716, and `(179, 179)` at mass 358; the eight
   outcome-0 entries `111/179`, `111/179`, `111/179`, `315/716`, `315/716`,
   `315/716`, `315/716`, `1/2`; and the off-diagonal modulus bound
   `7548/32041` of design section A field 6.
8. **Digest integrity.** The declared `custody_digest_sha256` equals the
   lowercase hexadecimal SHA-256 of the UTF-8 bytes of the canonical
   serialization of the export object with that field removed, with no
   trailing line break.

Digest convention of record: the design document names the committed
`canonical_sha256` of `oph_fpe/core/charged_response.py` together with the
indent-2 canonical text and the convention of
`plan/SIM_ALIGNMENT_2026-08-20.md` (the JSON text without a trailing
newline).  The committed function re-serializes its argument with compact
separators before hashing, so the two conventions name different byte
streams; this contract binds the digest to the indent-2 canonical text,
the same byte stream the canonical-serialization check pins.

## Verdict grammar and the provenance class

* `STATIC_COMMITTED_FIXTURE_CONFORMANT`: every check passes and
  `provenance_class` is `synthetic`. The shipped
  `sample_conforming_export.json` carries this verdict. It demonstrates
  static fixture/schema conformance and nothing else.
* `NONCONFORMANT`: any check fails.  The report lists every named error.

The `provenance_class` field is a validator-lane addition to the design's
section A field list and is mandatory. `producer` is a reserved input value,
but v1 always rejects it with `PRODUCER_AUTHENTICATION_UNIMPLEMENTED`, even
when the document contains no known synthetic marker. Self-declared module
hashes, commits, run identifiers, runtime logs, and a canonical document
digest are not authentication. V1 intentionally exposes no producer success
or evidential verdict. A later frozen v2 schema may add one only together
with an external verifier of the referenced artifacts and execution.

Named error codes: `FILE_UNREADABLE`, `JSON_PARSE`, `MISSING_FIELD`,
`UNEXPECTED_FIELD`, `SCHEMA_ID_MISMATCH`, `PROVENANCE_CLASS_INVALID`,
`NONCANONICAL_SERIALIZATION`, `NONCANONICAL_PATH`, `NONCANONICAL_TIMESTAMP`,
`FLOAT_OUTSIDE_DISPLAY`, `RATIONAL_ENCODING`, `SCALAR_ENCODING`,
`MATRIX_ENCODING`, `DIMENSION_MISMATCH`, `CARRIER_DIMENSION_UNSUPPORTED`,
`CONTEXT_COVERAGE`, `OUTCOME_COVERAGE`, `SPANNING_SET_INCOMPLETE`,
`EFFECT_FROM_KRAUS_MISMATCH`, `EFFECT_DECLARED_MISMATCH`,
`EFFECT_RESIDUAL_MISMATCH`, `EFFECT_NOT_PSD`, `TRACE_CHECK_MISMATCH`,
`CONTEXT_SUM_NOT_IDENTITY`, `EFFECT_COMPLEMENT_MISMATCH`,
`EFFECT_COMPLEMENT_NOT_PSD`, `KRAUS_NORMALIZATION_MISMATCH`,
`KRAUS_NORMALIZATION_NOT_IDENTITY`, `KRAUS_NORMALIZATION_RESIDUAL`,
`SUMMED_TRACE_MISMATCH`, `READBACK_COUNT_INVALID`, `MASS_INCONSISTENT`,
`COUNT_MASS_MISMATCH`, `BORN_WEIGHT_NOT_RATIONAL`,
`READBACK_TRACE_MISMATCH`, `READBACK_RESIDUAL_MISMATCH`,
`PREP_TRACE_NOT_ONE`, `PREP_DIAGONAL_NEGATIVE`, `PREP_CERTIFICATE_MISMATCH`,
`PREP_POSITIVITY_VIOLATION`, `PREP_RECORD_DIAGONAL_INCONSISTENT`,
`PREP_OPERATION_CLASS_MISMATCH`, `PREP_OPERATION_DESCRIPTION_MISSING`,
`PREP_HASH_MISMATCH`, `PROVENANCE_FORMAT`, `LABELS_INVALID`,
`SYNTHETIC_EVIDENTIAL_CONFLICT`, `SYNTHETIC_MARKER_IN_PRODUCER`,
`PRODUCER_AUTHENTICATION_UNIMPLEMENTED`,
`COMMITTED_EFFECT_MISMATCH`, `COMMITTED_RUN_STATE_MISMATCH`,
`COMMITTED_FREQUENCY_MISMATCH`, `DIGEST_FORMAT`, `DIGEST_MISMATCH`.

## What the validator cannot certify

The validator reads one JSON document.  It cannot certify source
production, provenance, custody, or run reality: whether any run happened,
whether the declared source records exist, whether the declared producer
modules produced the counts, whether the runtime read log is complete, or
whether the import-graph independence report describes any actual
execution.  The Lean binding interface proves only that its current
placeholder fields are freely stipulable and do not authenticate custody:
`committed_corpus_does_not_determine_binding` and
`binding_digest_free_parameter` in
`Lean/EventAlgebra/SourceBoundInstrumentInterface.lean` show that the
placeholder custody fields are freely stipulable over the committed
determined part, so a document passing every check here is constructible
without any source production, as the shipped synthetic sample
demonstrates by existing. A richer in-corpus or external data-bearing
construction remains viable. V1 therefore rejects every `producer` claim rather
than attach an authentic-sounding verdict to unauthenticated declarations.

The frozen decision rule of design section C is outside this checker.
No placeholder of section C.6 is bound here: no `TOL_READBACK`, no
`DELTA_PHASE`, no `EPS_PHASE`, no `EPS_ADD`, no interval construction, and
no verdict of the PASS/FAIL/SOURCE_PRODUCER_MISSING/INCONCLUSIVE grammar.
The readback checks here are the exact-equality reading that a
deterministic committed-table export satisfies; a produced run's sampling
residuals are gated by the frozen `TOL_READBACK` band under the frozen
rule, and the v1 committed-table cross-check binds the preparation and counts
to the old static literals. A produced run, including a route-R4 export whose
preparation lawfully departs from the committed run-state diagonal, cannot be
validated by v1. It requires a separately frozen v2 validator that carries
the finite-run tolerance, per-context preparation/custody records, and an
external provenance authenticator without hard-pinning the produced state or
sampled counts to the static fixture.

## What is not proved here

No run exists, no seed is drawn, no freeze event exists, and no producer
export exists.  The shipped sample is synthetic, marked as such in its
`provenance_class`, its labels, and its placeholder values, and it
discharges nothing.  Register rows PR-03, PR-64, and PR-65 are open;
PR-04 stays a declared row under its recorded disposition.  Nothing in
this directory discharges the source-production row, demotes any row,
freezes any instrument, or scores any comparison.  The `provenance_class`
field and the digest-convention binding above are declared choices of this
lane, recorded here pending any later registration under the owning lanes.

## No producer handshake in v1

The `custody_digest_sha256` check proves only that one JSON document is
self-consistent under the declared canonical serialization. It does not show
that the referenced source records, producer modules, runtime log, custody
manifest, or run exist, and it cannot supply the external custody datum of
`SourceBoundInstrumentBinding` evidentially. A future v2 authenticator must
resolve and hash the referenced artifacts independently, bind them to a
frozen rule and run receipt, and verify the source-to-outcome custody chain.
Until then there is no producer handshake, no producer verdict, and no route
from a v1 pass to a premise, instrument, or observation promotion.

## Files

* `ins03_export_validator.py`: the validator, the committed reference
  constants with cited sources, and the synthetic sample builder.  CLI:
  `python3 ins03_export_validator.py <export.json>`; exit 0 on a
  conformant static synthetic fixture, 1 on a nonconformant or reserved
  producer-class document, 2 on an unreadable or unparseable file.
* `sample_conforming_export.json`: the synthetic conforming sample, built
  from the committed objects only (the Lueders Kraus family of each
  committed context is the singleton projector family, since each
  committed effect is its own square root), with clearly marked
  placeholders in every producer-side field the corpus does not fix.
* `test_ins03_export_validator.py`: pytest coverage for the sample
  verdict, the mutation guards, the committed reference constants, and the
  fail-closed producer boundary.
