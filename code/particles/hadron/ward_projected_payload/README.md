# Ward-Projected Hadronic Transport Payload Harness (Generator G1)

The emitter is source-only and unscored. Target validation occurs in a
separate process after the source artifact has been sealed. No real payload is
currently scoreable: the corrective v3 target is permanently inactive, the present
emission is a singleton sampled-grid diagnostic, and `Delta_EW` remains open.

This directory implements the payload harness for the Ward-projected
`U(1)_Q` Thomson transport lane: a deterministic chain from a
spectral-measure module to `Delta_source` in inverse-alpha units, evaluated
at the internal fixed-point root of the Stage-5 chain in the same
subtraction scheme as `a0(P)`. The creation of this directory is the
declared payload-work start under the frozen v2 target's
`payload_work_start_definition`. That chronology is historical only; it does
not make V1 eligible for promotion.

Everything computable from first principles on this machine is computed:
the free-parton one-loop transport with internal Stage-5 masses, the
massless MS-bar R-ratio corrections driven by the chain's own
`alpha_3(m_Z;P*)` and the source-only Lambda_QCD lane, a declared
constituent-dressing branch, and the strict bracket over the declared
IR-cutoff and truncation grid. The remainder is a nonperturbative
Ward-projected spectral measure that no current method can supply at the
historical v2 precision scale; `PAYLOAD_STATUS.md` Section 4 preserves that
comparison as a diagnostic only. V3 registers no active precision or pass/fail
threshold.

## Usage

```bash
# one payload
python3 payload_harness.py --module parton_free

# source-only declared grid at the source-derived internal P
python3 run_bracket.py --output /sealed/source-bracket.json

# separate fail-closed comparison attempt, after sealing the emission
# (currently exits 2 with NOT_EVALUABLE)
python3 score_bracket.py \
  /sealed/source-bracket.json \
  ../../../../falsification/frozen_targets/hadronic_closure_target_2026-07-16_v3.json

# tests
python3 -m pytest test_payload_harness.py test_score_bracket.py -q
```

## Boundaries

- Inputs: Stage-5 internal masses, the chain's `alpha_3(m_Z;P*)`, and the
  dimensionless `Lambda3/mZ` ratio from the source-only transmutation lane.
- Excluded from all computations: CODATA/NIST alpha, measured hadronic cross
  sections, PDG hadron data, the empirical endpoint interval, and everything
  under `falsification/`.
- `run_bracket.py` contains no target value, tolerance, comparison rule, or P
  override. It uses the source-derived internal P and emits distinct total,
  residual, `S_QEW`, and `S_hadronic` coordinate kinds.
- The emitted bracket is explicitly an uncertified sampled-module envelope at
  one P. It is not a function over the P basin and is not an uncertainty
  interval.
- `score_bracket.py` is the only process that reads a target. It rejects an
  inactive target, unknown schemas, coordinate-kind errors, hash failures,
  P-domain mismatches, missing provenance, an open `Delta_EW` gate, and an
  uncertified artifact. It implements no scalar containment shortcut and no
  completed-map solver. Its receipt checks are schema and binding guards, not
  independent signature, timestamp, or interval-proof verification.
- Corrective target
  `falsification/frozen_targets/hadronic_closure_target_2026-07-16_v3.json`
  separates CL-1 and CL-2 but is a permanently inactive, post-target-access
  erratum scaffold. It cannot be made blind by timestamping it later. An
  eligible experiment needs a detached successor whose complete method was
  frozen before genuinely withheld data, or an audited clean-room producer.
- `promotion_allowed = false` on every source-bracket artifact. The corpus no-go
  (Theorem 6, `THOMSON_TRANSPORT_THEOREMS.md`) applies: this bracket does
  not determine the spectral moment and cannot be promoted to an endpoint
  derivation.
- Historical `runtime/ward_projected_payload_bracket_current.json` is retained
  byte-for-byte. It came from the contaminated legacy emitter and cannot be
  promoted or rescored into a promoted result.

See `V1_PROTOCOL_ERRATUM_2026-07-16.md` for the corrected process boundary and
the permanent historical V1 disposition. `PAYLOAD_STATUS.md` preserves the
technical status and labels its old comparison as historical.
