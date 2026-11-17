# Charged and neutral particle records

One signed source supports integer charge, a compact positive-energy state and
a measured half-turn relation between its mean and covariance. Separate source
extensions supply a localized internal clock, field spin, coherent field return
and neutral record transport. Their physical identifications require explicit
state, rotation, energy and clock maps. The [proof](PROOF.md) states the finite
results; the [claim ledger](../../claims/openwave_particle_records.json) retains
failed comparisons and unexecuted proposals; the [verifier](verify.py) checks
archived evidence without rerunning the scientific experiments.

```sh
python3 code/openwave_particle_records/verify.py
python3 tools/build_openwave_particle_records.py --check
```

Add `--workspace ..` to verify every original sealed leaf in the sibling
OpenWave workspace. `EVIDENCE.json` binds each manifest and selected exact
proof, receipt and code bytes. The `evidence/` tree is an archival selection;
its source scripts require their original dependencies and are not a standalone
simulation. Full replay instructions and dependency boundaries are in
[REPRODUCE.md](REPRODUCE.md).
