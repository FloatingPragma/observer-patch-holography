# Charge, rest energy and frame on one signed source

One signed source on the twelve-port carrier, run under the native repair law and
read by a complete current pointer, gives three results:

- **Charge comes in whole units.** An integer face certificate on the cell makes
  one signed pair a primitive class, so two sharp half charges lie outside the
  measured source image. The stopped flow E[F] = B^T L^+ q follows from the law:
  all 264 native trials absorb, all 48 changed-law controls do not.
- **A compact charged state has a positive rest energy.** On the full 630-seam
  graph the four-port charged ray has energy 6.014331228348673 in transfer units.
  Fixing the unit with 511 keV (a calibration), a disjoint pair comes out at
  1024.4508487 keV.
- **The current pattern responds to a half-turn.** The measured covariance returns
  under the executed half-turn, the nonzero mean reverses, and a second half-turn
  returns both.

```sh
python3 code/openwave_signed_source/verify.py
python3 tools/build_openwave_signed_source_claim.py --check
python3 -B code/openwave_static_signed_source/replay.py   # full experiment replay
```

The [proof](PROOF.md) states the construction and its scope. `recompute.py`
rederives all three results from the sealed carrier graph and the native rule.
The experiment records, manifests and standalone replay are in
[`code/openwave_static_signed_source/`](../openwave_static_signed_source/).
The claim ledger is [`claims/openwave_signed_source.json`](../../claims/openwave_signed_source.json).
