# Vendored external K7 sources

These files are third-party material, copied byte for byte so that the
[K7 / OPH gauge kinetic cross-control](../k7_gauge_kinetic_cross_control.md)
can be replayed offline. They are negative-control inputs only. They are not
OPH results, OPH evidence or OPH predictions, and the cross-control verdict
`NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES` establishes no
OPH physics.

The upstream PDF contains the illustrations referenced by the Markdown manuscript. [K7 framework 3.5, illustrated edition](https://github.com/Arithmon/K7/blob/210480b2ba8a13c98e6a6f146fbfd965e4904c00/publications/papers/pdf/k7_framework_3_5_main.pdf)

| Vendored file | Upstream repository | Upstream commit | Upstream path |
| --- | --- | --- | --- |
| `k7_gauge/ERRATUM_v3.5.md` | `Arithmon/K7` | `210480b2ba8a13c98e6a6f146fbfd965e4904c00` | `publications/ERRATUM_v3.5.md` |
| `k7_gauge/k7_framework_3_5_main.md` | `Arithmon/K7` | `210480b2ba8a13c98e6a6f146fbfd965e4904c00` | `publications/papers/markdown/k7_framework_3_5_main.md` |
| `k7_gauge/GaugeSector.lean` | `Arithmon/K7-Lean` | `2553ed170d04db8777f84fb2f9ef2b66b386f7fb` | `GIFT/Relations/GaugeSector.lean` |
| `k7_gauge_scale/ERRATUM_v3.5.md` | `Arithmon/K7` | `26d18fdd2fce61cc0231ec3485f8f5fe7efd5fc6` (K7 baseline of the scale contract) | `publications/ERRATUM_v3.5.md` |

The cross-control receipt and the scale contract in [`../runtime/`](../runtime/)
record the Git blob and SHA-256 digest of each file, and their verifiers
reject any drift.

Author of the K7 framework paper: Brieuc de La Fourniere. The vendored texts
carry no licence statement of their own, so their redistribution status is
unknown; the upstream repositories are the authoritative copies. The main
paper contains its own predictions and relative links to upstream files that
are not vendored here; those links do not resolve inside this repository.
