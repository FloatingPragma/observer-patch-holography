# Ward-Projected Hadronic Transport Payload Harness (Generator G1)

LABEL: development bracket, non-blind environment; the protocol pass requires
an isolated re-run.

This directory implements the payload harness for the Ward-projected
`U(1)_Q` Thomson transport lane: a deterministic chain from a
spectral-measure module to `Delta_source` in inverse-alpha units, evaluated
at the internal fixed-point root of the Stage-5 chain in the same
subtraction scheme as `a0(P)`. The creation of this directory is the
declared payload-work start under the frozen v2 target's
`payload_work_start_definition`.

Everything computable from first principles on this machine is computed:
the free-parton one-loop transport with internal Stage-5 masses, the
massless MS-bar R-ratio corrections driven by the chain's own
`alpha_3(m_Z;P*)` and the source-only Lambda_QCD lane, a declared
constituent-dressing branch, and the strict bracket over the declared
IR-cutoff and truncation grid. The remainder is a nonperturbative
Ward-projected spectral measure that no current method can supply at the
required precision; `PAYLOAD_STATUS.md` Section 4 states that wall with
numbers.

## Usage

```bash
# one payload
python3 payload_harness.py --module parton_free

# full declared bracket grid -> runtime/ward_projected_payload_bracket_current.json
python3 run_bracket.py

# tests
python3 -m pytest test_payload_harness.py -q
```

## Boundaries

- Inputs: Stage-5 internal masses, the chain's `alpha_3(m_Z;P*)`, and the
  dimensionless `Lambda3/mZ` ratio from the source-only transmutation lane.
- Excluded from all computations: CODATA/NIST alpha, measured hadronic cross
  sections, PDG hadron data, the empirical endpoint interval, and everything
  under `falsification/`.
- Canon target diagnostics (`S_required`, `c_Q`, `Delta_missing`) appear
  only in compare-only blocks labeled as non-blind development comparisons.
- `promotion_allowed = false` on every artifact. The corpus no-go
  (Theorem 6, `THOMSON_TRANSPORT_THEOREMS.md`) applies: this bracket does
  not determine the spectral moment and cannot be promoted to an endpoint
  derivation.

See `PAYLOAD_STATUS.md` for the implemented contract, the development
bracket, the precision-wall statement, and the blind-run protocol.
