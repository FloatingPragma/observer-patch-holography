# Conditional Ward spectral diagnostics

`run_bracket.py` evaluates 58 declared spectral variants at the internal
Stage-5 evaluation point and writes
`runtime/ward_projected_payload_source_bracket_current.json`. The coupling
uses four-loop running and high-to-low flavor matching. The Lambda grid
matches the conditional transmutation producer's central value and four
sampled threshold corners. The grid envelope is a numerical sensitivity
diagnostic, without a full uncertainty enclosure or a physical prediction.

The harness computes the declared source-chain and gauge-width pixel-map
quantities. The nonperturbative Ward-projected source measure, same-scheme
remainder and electroweak transport contribution are separate physical
premises. Generated diagnostic artifacts set `promotion_allowed=false`.
See [the payload contract](PAYLOAD_STATUS.md) for these requirements.

The retained `runtime/ward_projected_payload_bracket_current.json` is a
fixed comparison-bearing fixture with its own declared inputs and target
access. It is not an output of the source diagnostic producer or a receipt
for its corrected threshold convention. Its stored parameter values and
content digest identify the comparison protocol; replaying that fixture
does not certify the live source calculation.

```sh
python code/particles/hadron/ward_projected_payload/run_bracket.py
python code/particles/hadron/ward_projected_payload/payload_harness.py --module parton_free
python -m pytest -q -W error code/particles/hadron/ward_projected_payload/test_threshold_matching.py
python -m pytest -q code/particles/hadron/ward_projected_payload/test_payload_harness.py
```
