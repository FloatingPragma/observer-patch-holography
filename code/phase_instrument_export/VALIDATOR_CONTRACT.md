# INS-03 export validator contract

Schema under validation: `oph.sim.ins03_phase_instrument_export.v1`, the
simulator-side export interface specified in section A of
`plan/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md` (OPH meta planning
surface, outside this repository).  Validator:
`code/phase_instrument_export/ins03_export_validator.py`.  Owning lanes:
issues 730 (OL-C5, quantum) and 737 (instrument register).  The typed Lean
binding interface `Lean/EventAlgebra/SourceBoundInstrumentInterface.lean` is
the receiving half of this interface inside the corpus; this validator is
the receiving half on the export side.

## What the validator certifies

A `SCHEMA_CONFORMANT_*` verdict certifies, for one JSON document, all of the
following, checked fail-closed with named error codes and exact arithmetic
(`fractions.Fraction` plus a minimal exact `Q(sqrt(3))` implementation; no
floating-point value enters any check):

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
   (field 2), and the `trace_nonincreasing` values on the frozen spanning
   set of matrix units (field 3); `summed_channel[c]` with the Kraus
   normalization, its residual, and the spanning-set trace checks
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
     residual matrix exactly zero;
   * the declared spanning-set trace values equal the exact recomputation
     through the Kraus families, per outcome and for the summed channel;
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

* `SCHEMA_CONFORMANT_SYNTHETIC`: every check passes and `provenance_class`
  is `synthetic`.  The shipped `sample_conforming_export.json` carries this
  verdict.  It is a distinct verdict word from the producer one: a
  synthetic export can demonstrate schema conformance and nothing else.
* `SCHEMA_CONFORMANT_PRODUCER`: every check passes and `provenance_class`
  is `producer`.  A producer export must carry no synthetic placeholder
  values (no `SYNTHETIC_PLACEHOLDER` string, no all-zero commit or digest,
  no `example.invalid` URL); any such value fails closed.
* `NONCONFORMANT`: any check fails.  The report lists every named error.

The `provenance_class` field is a validator-lane addition to the design's
section A field list.  It exists to expose the synthetic/producer
distinction the lane requires and is mandatory; a later freeze of the
export schema may adopt or replace it.

Named error codes: `FILE_UNREADABLE`, `JSON_PARSE`, `MISSING_FIELD`,
`UNEXPECTED_FIELD`, `SCHEMA_ID_MISMATCH`, `PROVENANCE_CLASS_INVALID`,
`NONCANONICAL_SERIALIZATION`, `NONCANONICAL_PATH`, `NONCANONICAL_TIMESTAMP`,
`FLOAT_OUTSIDE_DISPLAY`, `RATIONAL_ENCODING`, `SCALAR_ENCODING`,
`MATRIX_ENCODING`, `DIMENSION_MISMATCH`, `CARRIER_DIMENSION_UNSUPPORTED`,
`CONTEXT_COVERAGE`, `OUTCOME_COVERAGE`, `SPANNING_SET_INCOMPLETE`,
`EFFECT_FROM_KRAUS_MISMATCH`, `EFFECT_DECLARED_MISMATCH`,
`EFFECT_RESIDUAL_MISMATCH`, `TRACE_CHECK_MISMATCH`,
`CONTEXT_SUM_NOT_IDENTITY`, `KRAUS_NORMALIZATION_MISMATCH`,
`KRAUS_NORMALIZATION_NOT_IDENTITY`, `KRAUS_NORMALIZATION_RESIDUAL`,
`SUMMED_TRACE_MISMATCH`, `READBACK_COUNT_INVALID`, `MASS_INCONSISTENT`,
`COUNT_MASS_MISMATCH`, `BORN_WEIGHT_NOT_RATIONAL`,
`READBACK_TRACE_MISMATCH`, `READBACK_RESIDUAL_MISMATCH`,
`PREP_TRACE_NOT_ONE`, `PREP_DIAGONAL_NEGATIVE`, `PREP_CERTIFICATE_MISMATCH`,
`PREP_POSITIVITY_VIOLATION`, `PREP_RECORD_DIAGONAL_INCONSISTENT`,
`PREP_OPERATION_CLASS_MISMATCH`, `PREP_OPERATION_DESCRIPTION_MISSING`,
`PREP_HASH_MISMATCH`, `PROVENANCE_FORMAT`, `LABELS_INVALID`,
`SYNTHETIC_EVIDENTIAL_CONFLICT`, `SYNTHETIC_MARKER_IN_PRODUCER`,
`COMMITTED_EFFECT_MISMATCH`, `COMMITTED_RUN_STATE_MISMATCH`,
`COMMITTED_FREQUENCY_MISMATCH`, `DIGEST_FORMAT`, `DIGEST_MISMATCH`.

## What the validator cannot certify

The validator reads one JSON document.  It cannot certify source
production, provenance, custody, or run reality: whether any run happened,
whether the declared source records exist, whether the declared producer
modules produced the counts, whether the runtime read log is complete, or
whether the import-graph independence report describes any actual
execution.  This is exactly the externality the Lean binding interface
proves: `committed_corpus_does_not_determine_binding` and
`binding_digest_free_parameter` in
`Lean/EventAlgebra/SourceBoundInstrumentInterface.lean` show that the
custody data of a source binding is freely stipulable over the committed
determined part, so a document passing every check here is constructible
without any source production, as the shipped synthetic sample
demonstrates by existing.  A `SCHEMA_CONFORMANT_PRODUCER` verdict states
that the document is well formed, algebraically exact, and free of
synthetic markers; it does not state that a producer ran.

The frozen decision rule of design section C is outside this validator.
No placeholder of section C.6 is bound here: no `TOL_READBACK`, no
`DELTA_PHASE`, no `EPS_PHASE`, no `EPS_ADD`, no interval construction, and
no verdict of the PASS/FAIL/SOURCE_PRODUCER_MISSING/INCONCLUSIVE grammar.
The readback checks here are the exact-equality reading that a
deterministic committed-table export satisfies; a produced run's sampling
residuals are gated by the frozen `TOL_READBACK` band under the frozen
rule, and the committed-table cross-check binds an export that declares
the committed contexts to the committed literals, so a produced run under
the frozen rule is evaluated by that rule and not by this validator.  A
route-R4 export whose preparation lawfully departs from the committed
run-state diagonal is likewise a matter for the frozen rule; this
validator's cross-check reports the departure and the frozen rule decides
its meaning.

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

## The handshake with the Lean binding interface

A genuine producer export that passes this validator with
`provenance_class` equal to `producer` supplies one value the corpus
cannot supply for itself: its `custody_digest_sha256`, the custody
artifact digest that the `custodyDigest` field of
`SourceBoundInstrumentBinding` in
`Lean/EventAlgebra/SourceBoundInstrumentInterface.lean` consumes under a
declared mapping of hex digests into the natural numbers.  The Lean
theorems prove that the committed corpus determines every other field of
the binding and leaves exactly that custody data external; this validator
is the conformance gate the external data passes through on its way in.
Registration of such an export, the binding of the section C.6
placeholders, any freeze, and any evidential use are owner and
preregistration actions under the freeze discipline of issue 737 and are
not performed, prepared, or implied by a validation verdict.

## Files

* `ins03_export_validator.py`: the validator, the committed reference
  constants with cited sources, and the synthetic sample builder.  CLI:
  `python3 ins03_export_validator.py <export.json>`; exit 0 on a
  conformant export, 1 on a nonconformant one, 2 on an unreadable or
  unparseable file.
* `sample_conforming_export.json`: the synthetic conforming sample, built
  from the committed objects only (the Lueders Kraus family of each
  committed context is the singleton projector family, since each
  committed effect is its own square root), with clearly marked
  placeholders in every producer-side field the corpus does not fix.
* `test_ins03_export_validator.py`: pytest coverage for the sample
  verdict, the mutation guards, the committed reference constants, and the
  synthetic/producer verdict distinction.
