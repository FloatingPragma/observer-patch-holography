# Preparing the native feedback reference

The actual native source/current-sign controller converges within a certified
finite number of attempts. Its classical reference has support 129 inside the
132-state source/bit carrier. The certificate uses every input, including the
three transient states, and checks the complete C and E record instruments.
It concerns the continuing one-record observer with older records exported,
repair probability one and uniform exploration probability one half.

Let `K` be its exact transition matrix. Each entry has denominator 120.
The executable computes `K^64` with integer arithmetic and the common measure
`a_j = min_i (K^64)_ij`. Its total mass exceeds `2/5`. Normalizing `a` gives
`nu`, and the nonnegative residual has the stochastic decomposition

```
K^64 = (2/5) 1 nu + (3/5) R.
```

For probability laws `p,q`, the constant part cancels in their difference.
A stochastic matrix cannot increase total variation, so
`TV(p K^64, q K^64) <= (3/5) TV(p,q)`. The stored exact reference `pi` is
checked against the actual native matrix, including normalization and every
stationarity equation. Thus for every classical initial law,

```
TV(p K^n, pi) <= (3/5)^floor(n/64).
```

The full record instrument uses only the source and current-sign bit as
inputs to its next step. Let `Q` project the last record to that bit and
`L` write one complete record. Exhaustive endpoint branch checks establish
`Q L = K`; the next full-record marginal is `L Q`. Hence one extra attempt
transfers the same bound to the stationary source/last-record law `L pi`.
Fixed tagged mixtures of C and E follow by linearity. All emitted branches
remain part of the instrument, and every overwritten record is exported.
The total archive grows and is not assigned a stationary distribution.

This supplies preparation accuracy for the finite reference used by the
[theorem-to-observer map](../completion/THEOREM_MAP.md). It does not identify
native attempt counts with physical duration, implement a modular clock, or
prove decay of arbitrary source coherences. The separate
[balanced adaptive controller](../adaptive_reference/README.md) has its own
quantum mixing theorem on a different state space. Its constants do not
transfer to this controller.

Run from the repository root:

```sh
python3 -m model.adaptive_preparation --check
python3 -m pytest -q model/adaptive_preparation
```

The receipt contains the exact common-column measure, conservative rational
budgets and source digests. The argument is an elementary finite stochastic
proof checked by exact arithmetic; no additional Lean formalization is
asserted.
