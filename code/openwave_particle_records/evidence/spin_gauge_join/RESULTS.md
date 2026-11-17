# Finite same-field join: passed, with physical response still open

The corrected registered run passed in44.181seconds. It used no measured
nature constants and made no new complete M6 execution claim for the
integer/tag services.

| Check | Result |
| --- | --- |
| CAR/tagged-spin Gauss carrier intertwiner |332,640 coherent basis columns, exact |
| Native C/E coin instrument |665,280 coin columns, exact effect sums and partial injections |
| CAR crossing signs retained |13,888 negative-sign columns |
| Whole-event overflow retained |11,200 columns |
| Missing-sign control |13,888 failing columns |
| Missing-electric-update control |22,400 failing columns |
| Spin-reversal control |22,400 failing columns |
| Proposed-current-on-overflow control |11,200 failing columns |
| Same-field loop |All4 spin inputs; identity on their full spin coherence |
| Complete two-block pointer law, phase0 |5/16,3/16,3/16,5/16 |
| Complete two-block pointer law, phase1 |3/8,1/8,1/8,3/8 |
| Nonselective electric energy, supplied kappa=1/8 |1/16,5/32,11/32 |
| Finite current/path source dimensions |486 columns for each q=0,1; spin-current dimension972 |
| Independent two-step path records |4,374 for each q |
| Complete current-query unitary completions |3 per q, all columns exact |
| Maximum Pauli-square/path reconstruction error |1.11e-16, binary64 diagnostic |

The two-path cross response is nonzero already at zero integer electric
coupling: Frobenius norm42.21493363832878 for q=0, versus28.143289092219053
for q=1. Both have maximum entry0.7817580303394194. In the latter finite
window each seam has216 moving columns and108 overflow columns, versus324
moving and0 overflow at q=0. These norms are not a magnetic-moment estimate:
they include geometry, active-memory characters and finite-boundary effects.

The new result is an exact join of the conserved CAR field spin to the
existing signed-integer Gauss and coherent loop carrier, together with a
complete specified spin-current/path reader on that object. The current
contraction has not been selected as a physical kinetic generator. No
stationary sourced dipole, g or Larmor result is claimed. The local tagged
spin has the existing half-unit field representation; electron, rotation
and unit identification remain open. The entire neutral pair is not
misidentified as a single spin-half particle.

An earlier implementation preflight failed because its selector-to-native-
coin conversion was reversed. It emitted no results. Its exact reconstructed
code, patch, traceback and original DESIGN/TASK-only preregistration are
preserved under `attempts/01_invalid/`. That attempt lacked a preregistered
code/source manifest; this custody limit is explicit. Root verified and
registered the corrected code and192 input hashes at
2026-11-14T17:05:03.858818+00:00, before the successful run.

Reproduce the complete certificate:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave/shared/spin_gauge_join/verify.py
```

This compares the receipt, exact source pins and every retained path array.
`python3 openwave/shared/spin_gauge_join/custody.py` checks sealed file/input
hashes without repeating the calculation. These are internal producer
checks, not an independent G6 review or platform acceptance.
