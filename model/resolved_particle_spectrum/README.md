# Full resolution on the protected eight-line support

This successor performs the previously prospective sixteen-setting spectral
experiment by **certified accelerated emulation**. Counts at times0,39,...,273
and pointer phases0,pi/2 recover all eight weights on the source-specified
support. Every inferred interval is strictly positive and contains its
independent exact source residue. The calculation does not infer arbitrary
unknown frequencies or identify physical particle species.

The reference is the same excited fixed-N6 state, the local probe is the same
one-site number phase, and the action is the original2x2x2 Weyl/contact cube.
The observable sector is its eight-dimensional inelastic cyclic response,
plus the known elastic reference weight1/16. It is not all119 eta-three
states or all1820 protected states. Its frequencies are E_H-3; changing the
common energy zero leaves the readout invariant.

## Why millions of steps can be evaluated exactly enough

The original chronological source word lists the two spin hoppings on
each edge consecutively. Their four mode addresses are disjoint, so the
two generators commute. Their sum is exactly the original complete Weyl
edge and commutes with physical eta-plus. Grouping them changes no ordering.
Each hopping half-word therefore commutes with eta-plus squared.
All onsite contacts commute with each other, and their complete sum D obeys

`D_N6 T = T D_N2 + 2T`, with `T=(eta-plus)^2`.

Consequently each original ideal symmetric step has the exact intertwiner

`S_N6(h) T = exp(-2ih) T S_N2(h)`,

while the eta-four reference R gets phase exp(-3ih). The finite four-dimer
complete repairs are identities on these ideal states, with all native
zero and outside outcomes retained. Isolated spin entries and partial
contact layers are not asserted to have this intertwiner.

`factorization.py` checks the source ordering, every Gaussian-integer edge
intertwiner, the contact identity and the actual singlet chart. A normalized
36-dimensional singlet chart retains all35 eta-three singlets and the
eta-four parent. In the actual Hubbard onsite frame, let U be the spatial
8x8 forward half-hop product and W=Sym²U on normalized unordered pairs.
Then the complete step is

`S = W^T diag(exp(-ih D)) W`.

Transpose is appropriate because every local Hubbard hopping matrix is
complex symmetric; it reverses their order without reversing time. The
construction includes certified sqrt2 factors for distinct-site singlets.
The probe q is supported on diagonal pair states and has exact integer
norm squared56, so its expectation needs no numerical normalization.

## Exact arithmetic and complete native error

The common step is h=39/2^19. At time39j the word has j*2^19 steps, divided
equally around complete C and E repair. At the longest time273 this is
3,670,016 steps, or38,755,368,960 field pulses before preparation/readout.

All step entries have outward160-bit interval enclosures. The power routine
rounds exact Gaussian-integer products onto128-bit dyadics. For approximants
to unitary targets with errors a,b, each multiplication has error at most

`a+b+ab+d*2^-128`, with d=36.

The error is itself rounded outward. Thus binary powering has a proved
operator error, not a floating matrix residual. The longest power takes24
matrix products and has error249561087/2^128, below7.34e-31.

The actual native word is separate from this ideal compression. Its full
10,560-pulse step is consumed against independent literal-source templates,
including all routed FSWAPs, phases, loads and export complements. The exact
complete step error is below1.323294e-12. Telescoping over every repeated
step includes every possible first-failure branch. Full-Fock splitting adds
`(629/256)*t^3/n^2`. Including initialized source preparation, pointer gates,
all complete repairs, certified power arithmetic and exact binary64
probability-conversion error, every setting has systematic probability
error below8.568121e-6, hence below1/65536.

The long native pulse stream and its individual failure weights are **not
executed or materialized**. `budgets.py` gives the compact repeated-word
description, all pulse/export counts and record addresses. A record address
is `(half,step,pulse,complement)`; the first failure terminates, while outside
repair continues without normalization. The complete-history norm theorem
bounds every native event probability without pretending that a projected
rounded native word has the exact ideal symmetry.

## Sixteen simulated settings and eight inferred weights

The executed computation powers the full36-state ideal word and forms the
two Ramsey probabilities, rather than summing supplied spectral weights.
It samples a binary64 proxy to these certified probabilities using PCG64
seed20261102. Exact rational conversion error is included. These are counts
from the stated surrogate probability law, not exact samples of the
unmaterialized native failure instrument. Their distance from the native
probabilities is bounded as a systematic probability error; no claim of
identical individual native and surrogate records is made.

There are9*2^29=4,831,838,208 independent reset trials per setting, for
77,309,411,328 simulated trials overall. Zero is contrasted with **all other
complete outcomes**, without survival postselection. Under the specified
independent binomial/reset service, Hoeffding and a union bound give
simultaneous sampling error1/32768 except with probability below
`32*(3/8)^9 < .01`. The same conditional recovery theorem applies to actual
native counts satisfying the stated sampling event, using their proved
systematic probability bounds.

`inversion.py` uses only the eight exact source frequency intervals and the
counts. Rational Lagrange polynomials propose a96-bit inverse B of the
sampling matrix V. Outward arithmetic certifies `||B||_infinity<2` and
`||BV-I||_infinity<1e-18`. For `y=(16C-1)/15`, the error is bounded by

`||B||_infinity*(16/5)*(systematic+sampling)+||BV-I||_infinity`.

The residual term uses the independently proved fact that the true weights
form a probability vector, hence max|w_j|<=1. No true numerical weight is
used to construct B or infer the weights. Comparison with the source
residues happens only after reconstruction.

| Line in ascending frequency order | Recovered weight |
| --- | ---: |
| 1 | .002457894673 |
| 2 | .005641397395 |
| 3 | .005890881572 |
| 4 | .007131698528 |
| 5 | .004528321839 |
| 6 | .421578515341 |
| 7 | .415796456292 |
| 8 | .136974834359 |

Each has conditional error below.000258393828. Every interval is positive
and contains the independent exact weight. This closes full weight
resolution on the supplied support in the certified emulator. Native
hardware execution, blind support discovery, relativistic mass shells,
physical charge/spin and Standard Model species identification remain
separate questions.

```sh
.venv/bin/python -W error -m pytest model/resolved_particle_spectrum -q
.venv/bin/python -W error -m model.resolved_particle_spectrum --check model/resolved_particle_spectrum/results.json
```
