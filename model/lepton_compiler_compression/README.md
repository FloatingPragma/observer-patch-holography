# A smaller complete lepton compiler

The fixed neutral compiler in
[lepton_mass_compiler](../lepton_mass_compiler/README.md) has an exact
implementation with **90 physical and 12 ancillary CAR modes per site**.
The original 323,910 ancillary modes are reconstructed by an explicit
number-preserving archive isometry; no failure, quantum reference or
coherence within a measured flag sector is discarded.

The actual source projector selects a two-dimensional address span for
each character. A sharper leakage calculation also replaces the
`108 V² t²/N` complete-output bound by
`2 V (6−2 sum e_k²)t²/N < 12 V t²/N`. The round count and quantum archive
both remain finite resources; the quoted precision campaign is a bound,
not an executed run.

```
VECLIB_MAXIMUM_THREADS=1 python3 -m pytest -q model/lepton_compiler_compression
VECLIB_MAXIMUM_THREADS=1 python3 -m model.lepton_compiler_compression --write
VECLIB_MAXIMUM_THREADS=1 python3 -m model.lepton_compiler_compression --check
```

[PROOF.md](PROOF.md) states the arbitrary-occupation and arbitrary-reference
domain, the exact archive correspondence, the limited minimality claim,
and the block-Duhamel bound. This reduces the declared implementation;
it supplies no passive clock or physical Standard Model identification.
