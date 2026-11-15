# Fewer full-field source steps at the same certified error ceiling

The [retained full-band clock](../band_clock/README.md) uses 32 symmetric
splitting steps. The [contact spectral enclosure](../pair_cyclic_control/README.md)
reduces the clock-rate error term without changing its action or reference.
This separate execution variant uses that improvement to reduce the number
of source steps. It retains the entire 8,128-state field sector, all momenta,
scattering states, Weyl tastes, arbitrary untouched references and all four
occupation-pointer outcomes. The original execution and receipt are retained.

For dimensionless modular parameter `|s|<=1`, the same native count decoder
gives the declared rate window and statistical failure bound. The uniform
complete average half-diamond estimate for `n` splitting steps is

```
epsilon(n,s) = min(1,
    (7/4)^3 (5|s|/4)^3 / (3 n^2)
    + 77|s|/2560 + 381169/10^12).
```

The first term is the retained ideal symmetric-product estimate. The second
uses the proved full-field centered generator norm `77/96` and the same
rate-error bound `3/80`. The last term is the unchanged statistical failure
allowance. This estimate is uniform over the full field with an untouched
reference; contraction covers the complete occupation instrument. Source
floating residuals and physical apparatus errors are separate from this
ideal-word bound.

The original ceiling at `|s|=1` is approximately `0.03856403003`. The improved
estimate is approximately `0.03799047087` with 21 steps, while 20 steps give
`0.03880144725`. Thus 21 is the smallest positive integer certified by this
particular bound under the original uniform ceiling. This is not optimality
over all source programs, other splitting formulas, or tighter estimates.
`minimum_steps` uses exact rational arithmetic and integer square roots;
it can also exploit a smaller declared `|s|`. It does not fit output data.

At `|s|=1`, the complete evolution and the two pointer services require
17,219,334 supplied pulses instead of 26,238,982, saving 9,019,648 pulses
(about 34.4 percent). Routing and every proper-code load/export operation
are counted. Preparation of the band packet and global Gibbs reference,
and the retained native record acquisition, are separate supplied costs.
No pulse count is interpreted as energy or physical time.

The implementation calls the same full-field source evolution with the
smaller integer schedule. The reference is the same supplied global Gibbs
state; the rate-transfer actuator is the same supplied service. Stop records
keep the complete input and reference. Random complex-reference histories
over the full sector are compared with an independently exponentiated
Hamiltonian, and the accounting test checks every pulse category including
both pointer services.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/band_clock_budget
```

This reduces a finite simulation cost. It supplies no native field encoding,
physical clock calibration, particle mass or continuum limit.
