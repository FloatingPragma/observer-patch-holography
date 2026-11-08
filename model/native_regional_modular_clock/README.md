# Regional modular clocks on the fair native source

The same transfer reference `rho=P/Tr(P)` has an exact regional density,
including every supported local occupation sector. After rescaling modular
time by `2M`, its regional generator approaches a local exchange generator
with an explicit density factor and boundary occupation potential.
[PROOF.md](PROOF.md) derives the restriction, signed logarithm remainder,
conditional preparation target and exact fixed-cut obstruction.

Fixed nested regions do not have compatible modular flows. A growing collar
provides controlled interior compatibility instead: [COLLAR.md](COLLAR.md)
derives a finite exchange propagation bound and joins it to the regional
modular estimate. On the authenticated level-14, radius-four source window,
the flow error at absolute exchange time `1/100` and density `1/2` is below
`1.068e-6` per norm-one, number-preserving root observable. The source has 348 window ports,
905 internal seams, 17 crossing seams and port collar distance nine.

The matrices and finite windows are reconstructed from the original native
source. No growing global occupation matrix is allocated. Local ratios use falling
factorials; the arithmetic tests include capacity `2^400`. The owning replay
checks all 4,356 entries of the H66 transfer and all thirty seam restrictions,
plus the finite routed collars. It binds the inherited rooted-source files
and the native transfer source files; it does not rerun every inherited
producer or claim a full Python-import closure.

A scalable global preparation of `rho` remains a premise. The positive
sampler requires a canonical configuration sampler and an addressed
coherent seam mixer; repair-only classical feedback cannot create its
coherences. An infinite port net requires an all-level routing or a
justified rooted subsequence. Growing regions require both the stated
capacity bound and a density rate. This is a graph exchange clock, with
explicit modular time orientation. It does not identify physical Lorentz
boosts, a relativistic vacuum, responsive geometry or Einstein dynamics.

Run the exact owner and its independent tests from the repository root:

```sh
python -m model.native_regional_modular_clock.packet --write
python -m model.native_regional_modular_clock.packet --check
python -W error -m pytest -q model/native_regional_modular_clock
```
