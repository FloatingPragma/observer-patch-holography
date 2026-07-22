# Ward-Projected Hadronic Transport Harness

This directory emits source-side transport diagnostics for the source-chain and gauge-width pixel maps of the
P-closure issue (#545). It contains no comparison target, tolerance, decision rule, scorer,
or promotion mechanism.

The harness computes the declared source quantities at the internal Stage-5
root. The nonperturbative Ward-projected spectral measure, the same-scheme
remainder, and the electroweak transport contribution are work in progress.
The emitted sampled-grid envelope is not a certified interval or a physical
Thomson prediction.

## Usage

```bash
python3 payload_harness.py --module parton_free
python3 run_bracket.py --output /tmp/source-bracket.json
python3 -m pytest test_payload_harness.py -q
```

No CODATA/NIST alpha value, measured hadronic cross section, PDG hadron datum,
or empirical endpoint interval enters the source calculation. Generated
artifacts set `promotion_allowed = false`.

See [`PAYLOAD_STATUS.md`](PAYLOAD_STATUS.md) for the source quantities and
missing physical construction.
