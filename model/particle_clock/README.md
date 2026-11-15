# A native record clock for one interacting energy pair

This finite CP-02/CP-13 attachment identifies the **same two-dimensional
particle subspace** in a declared interacting Weyl/contact action, a faithful
reference recovered from native records, and a complete source-executed
energy interferometer. It adds the operator map missing from a comparison of
unrelated spectral gaps. The field encoding and action parameters are supplied;
physical particle identification and clock calibration remain open.
The shared operator is executed through
[`bound_pair.fiber.clock_subspace`](../bound_pair/fiber.py), rather than
reconstructed from a separate spectrum fixture.

On an even periodic cube of side `L>=4`, take two fermions with
`h(k)=-(1/8) sigma.sin(k)` and contact interaction `U=1`. At total momentum
`K=(pi,pi,pi)`, `h(K-k)=h(k)`. The spin singlet is annihilated by
`sigma_a tensor I + I tensor sigma_a`. Consequently the onsite eta pair `d`
and the zero-mean offsite singlet `r`, whose relative momentum amplitude is
`sqrt(2) sin(k_x)`, span an invariant space `S` with

```
H r=0, H d=U d, H|S=diag(0,U).
```

The supplied isometry maps native occupations `{1,2}` to `r` and `{7,2}` to
`d`. The canonical native product reference assigns these source states
`200/7505` and `64/7505`. Its corner is therefore

```
P pi P = (264/7505) rho_S,
rho_S = diag(25/33,8/33) = exp(-beta H|S)/Z,
beta U = log(25/8).
```

This supplies the faithful reference and exact Gibbs/modular identity on
`B(S)=M2`: Heisenberg model coordinate `t=-beta s`, where `s` is the modular
parameter. It does not make this rank-two state faithful on the full field,
derive a laboratory temperature or second, or equate native repair with the
interacting Hamiltonian. The many-site isometry is a mathematical observer
dictionary; its physical source encoding has not been compiled here.

The estimate comes only from the complete native coin counts. The two selected
ports have distance two in the fixed inference tree. On the inherited
simultaneous confidence event, the error in `beta U` is at most `3/80`; after
subtracting a common phase, the clock half-diamond error is at most `3/160`
per unit modular parameter. The unchanged statistical failure allowance is
added to the complete average error. No rate or mass is inferred by fitting
the field response.

The executed protocol keeps both outputs of a nondemolition corner herald,
applies the existing balanced source probe on success, performs one uniformly
selected native C or E repair, and heralds the corner again. Outside outcomes
halt with their actual source/reference state. On return it applies the
record-derived source clock, the analyzer, and the native `(1,7)` detector.
Coin and current together decode the pre-detector occupation, hence the
declared `0/U` energy label. Source-address export aborts remain separate
outputs. Structurally zero raw slots remain available through `raw_output`.

The coarse herald is an actual source word: reversible CNOT fan-in makes one
work bit free, and a fully executed nine-qubit controlled-X macro copies the
other eight work bits' predicate to a blank pointer. Its 2,068 pulses preserve
coherence inside the accepted corner. Certification checks all 512 original
inputs, including padding, and both pointer outcomes. Preparation, pointer
initialization/readout, coherent services and record storage retain their
explicit apparatus obligations. Proper-code failure mass is included in the
bound; those processor-environment outputs are bounded, not materialized.

The exact surviving corner map has population kernel

```
[[859/990, 5/198],
 [ 4/495,404/495]]
```

and coherence multipliers `alpha_C=1519/1980`, `alpha_E=37/60`. After the
balanced probe, return probability is `85/99` and outside probability is
`14/99`. The energy-`r` record probability, divided here **only to display
the analytic coefficient**, is

```
p_r / (264/7505) = 85/198 + (17/66) alpha_law cos(s log(25/8)).
```

The reported complete instrument retains the initial outside mass, subsequent
leakage, every native record and abort. Its probabilities and error bounds are
never conditioned on a successful herald. Preparation error is charged before
the first herald, so it is not divided by the small corner probability.
Incompatible training or source certification produces an explicit stop on
the same terminal carrier. Fixed admission caps, rather than a realized
transcript's error, give the unconditional average comparison.

Both reference levels have fermion number two: `r` is a zero-contact singlet,
not a vacuum state. The measured `U` is their energy difference in the supplied
action convention. This is a lossy energy interferometer for that pair. It is
not autonomous protection of every superposition, a full interacting Gibbs
preparation, an interacting continuum limit, a rest mass in physical units,
or an accepted integrated completion gate. It strengthens the finite
Gibbs/modular attachment in [the theorem map](../completion/THEOREM_MAP.md)
and the earlier [fixed thermal clock](../common_source/THERMAL_CLOCK.md) by
giving this concrete native preparation/record map and a complete continuing
experiment on one specified interacting excitation subspace.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.particle_clock --write
OPENBLAS_NUM_THREADS=1 python -m model.particle_clock --check
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/particle_clock
```
