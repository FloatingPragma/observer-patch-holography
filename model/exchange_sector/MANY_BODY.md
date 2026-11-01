# Native exchange across all occupation sectors

The fixed twelve-port carrier has an exact **scalar exchange catalogue for
every binary occupation number**, together with two compatible complete C/E
quantum instruments and a source-executed measurement distinguishing them.
This extends the [two-particle classification](STATISTICS.md) on the same
native seams and declared local deformations. It does not assign a particle
species or construct a spatially protected core.

## Domain and complete catalogue

Let X_N have the N-element subsets of the twelve native ports as vertices.
An edge moves one occupied port along a native seam into a vacancy. Attach
exactly two kinds of cells: a one-particle traversal of a native face with
all other particles outside that face, and a square of two moves on
vertex-disjoint seams. Face fillings are declared path relations; they do
not follow merely from allowing the three edge moves. All occupation states
and all such cells are included, independently of a mass target.

| N | States | Movement edges | Face cells | Square cells | Integral H1(X_N) |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 1 | 0 | 0 | 0 | 0 |
| 1 | 12 | 30 | 20 | 0 | 0 |
| 2 | 66 | 300 | 180 | 315 | Z/2 |
| 3 | 220 | 1350 | 720 | 2520 | Z/4 |
| 4 | 495 | 3600 | 1680 | 8820 | Z/6 |
| 5 | 792 | 6300 | 2520 | 17640 | Z/8 |
| 6 | 924 | 7560 | 2520 | 22050 | Z/10 |
| 7 | 792 | 6300 | 1680 | 17640 | Z/12 |
| 8 | 495 | 3600 | 720 | 8820 | Z/14 |
| 9 | 220 | 1350 | 180 | 2520 | Z/16 |
| 10 | 66 | 300 | 20 | 315 | Z/18 |
| 11 | 12 | 30 | 0 | 0 | Z^19 |
| 12 | 1 | 0 | 0 | 0 | 0 |

For a connected complex, scalar flat unitary transports modulo vertex-frame
changes are Hom(H1,U(1)): fix the frames along a spanning tree and impose the
integral cell relations on the remaining edge phases. Thus at N=2,...,10
the possible scalar exchange phases are precisely the 2(N-1)th roots of
unity. At N=11 there are nineteen independent U(1) phases. These are
transport choices, not nineteen physical charges or species.

[`many_body_topology.py`](many_body_topology.py) proves the table by integral
presentations. Contracting a checked spanning tree gives a basis of the
integral cycle lattice. Every eliminated relation has a coefficient +1 or
-1, so eliminating its generator is an invertible integer presentation
operation. An independent replay checks each substitution against the
original cell boundaries. At N=2,...,10 exactly one generator remains;
the gcd of its residual relations is 2(N-1). Reverse substitution produces
an edge cocycle taking value one on **every** coherently oriented local
triangular exchange with N-2 fixed outside spectators. These exchanges
therefore generate, rather than merely detect, each cyclic group.

At N=11 no cell is admitted. The connected movement graph has first Betti
number 30-12+1=19. Its twenty coherently oriented hole-face loops have sum
zero; deleting any one leaves an integral basis, checked by determinant
equal to +1 or -1. This proves completeness of the face-loop coordinates.

The full integer elimination traces are reproducible through
`reduction_certificate(N)` and checked by `verify_reduction(N, trace)`.
The retained receipt stores their hashes and census rather than thousands
of redundant intermediate rows. No floating rank or sampled loop replaces
an integral relation. Independent tests reconstruct boundaries, use modular
rank controls and reject missing cells and corrupted elimination steps.

For N=2,...,10 the abelian groups match the known spherical-braid
abelianization; see [Goncalves and Guaschi](https://arxiv.org/abs/math/0603701).
That comparison supplies context, not the finite proof. No equivalence to
the continuum configuration space, nonabelian fundamental group for N>=3,
or higher-rank representation classification is claimed. The N=11 result
explicitly prevents extending the continuum comparison to every density.

## Exactly two spectator-compatible scalar families

A single family of identical excitations can require that exchanging two
particles locally have the same holonomy when stationary spectators are
inserted outside the operation. This is a precise compatibility requirement
on the proposed family, not a conclusion from population probabilities.
[`many_body_compatibility.py`](many_body_compatibility.py) proves:

**On the declared finite complexes, exactly two scalar flat families, up to
vertex gauge, satisfy this spectator compatibility: one common exchange
character +1 or one common exchange character -1.**

Each local exchange is a primitive cyclic generator for N=2,...,10.
Spectator insertion therefore fixes its phase to the N=2 value z, which
obeys z^2=1. The twenty generating N=11 face loops must have the same value;
their sole relation imposes z^20=1, satisfied by both signs. The remaining
sectors are trivial. Existence is constructive: follow the actual ordered
labels along each move, and assign either +1 or its permutation sign. Both
are flat on every admitted cell. This proves exhaustiveness and existence
without assuming a creation algebra to select a sign.

For a vacant spectator port p, the insertion chart is

```
J_p^chi |S> = chi^(number of occupied ports below p) |S union {p}>.
```

It intertwines every edge and filling on the subcomplex whose vertices,
edges and fillings avoid p. This domain restriction matters: a filled
one-particle face containing p visits a configuration already occupying p,
so insertion cannot extend over that filling. There is no unrestricted creation map of the entire declared
two-complex. The exact inventory checks 24,576 vertex inclusions, 153,600
edge inclusions and 368,640 cell inclusions. Independent sector choices such
as the legal N=3 phase i fail this common-family contract; they remain valid
on that isolated sector.

## Complete quantum continuation and particle-hole distinction

[`many_body.py`](many_body.py) uses the existing full-domain repair operation
on all 4096 binary occupations and every native seam. For each coin and raw
record, its unweighted map is a partial injection, multiplied by either
one or the CAR sign across occupied spectators. Both retain the fair factor
1/sqrt(2), every C/E raw slot, zero branches and abort slots. Complete effects
sum to I and output projectors sum to I. Total occupation and the pathwise
local continuity equation hold on every column. Classical record and
occupation laws agree for both lifts, including arbitrary untouched
references; coherent controls can distinguish them.

The negative lift is the even CAR instrument: moving blocks are
c_b^dagger c_a and holding blocks are endpoint-number projectors. Every
raw branch therefore commutes with creation and annihilation at an
unaddressed port. Exact matrix-column checks verify the full CAR algebra
and these continuation identities over every number sector. The positive
lift supplies the hard-core occupation control, not unrestricted bosonic
CCR. Odd creation operators are algebraic maps here; a physical
number-changing service needs its reservoir and preparation account.

Particle-hole complement preserves the movement graph but **not** the
declared face fillings. A filled N=1 face maps to a nonzero unfilled N=11
cycle. The full CAR basis map

```
C|S> = (-1)^(sum of occupied port indices) |complement(S)>
```

has C^2=I and C c_p C=c_p^dagger. Complementing the coin and reversing the
recorded current gives
`C K_h = (-1)^abs(j_h) K_complement(h) C`.
The extra minus on each move cannot be removed by site phases on the native
triangular graph: a three-move particle loop has phase +1, its one-hole
partner phase -1. Each indicated complete history has probability 1/8,
not one. The minus is constant within each raw branch and cancels in its
completely positive map. The committed-record instruments are therefore
particle-hole covariant after the stated input/output and record relabeling;
these loop phases concern coherent path amplitudes, not distinguishable
classical history probabilities. This is an exact CAR/native transport distinction, not an
identification of electromagnetic charge conjugation or a revision of the
separate electron/positron field benchmark.

## Executed observer and complete budget

[`many_body_service.py`](many_body_service.py) compiles the endpoint rule on
native seam (0,5), with active native ports (0,1,5,6) and all other ports
vacant. These are local mode indices (0,1,2,3). The circuit prepares a fair
coin, records the full current and E parity, applies the chosen sign
feedback, and moves the occupation. All sixteen input columns are checked
against independent native endpoint arithmetic; the all-4096 proof above
is exact algebra, not a claim that all those source columns were executed.

The observer prepares a fixed-total-two state proportional to
`|{0,1}> + |{0,6}>`, applies complete repair, and recombines the spectator
modes. The moving branch has opposite interference in the two lifts;
the holding branch and every other output remain. An explicit which-path
pointer removes this distinction in the complete classical terminal law,
while retaining the pointer itself. Keeping that pointer coherent for a
different erasure experiment would be a different readout.

| Executed observer | C pulses | E pulses |
| --- | ---: | ---: |
| Positive lift, preparation through final read rotation | 305 | 311 |
| Negative lift, preparation through final read rotation | 325 | 331 |
| Additional retained which-path pointer | 3 | 3 |

Each pulse includes its proper-code load and export. Every first failed
export remains in an immutable prefix archive, including its distinct
complement Kraus outputs. C uses seven blank work bits and 128 terminal
addresses; E uses eight and 256. The which-path control adds one bit and
doubles the address count. Terminal measurements retain every address.

Fresh complete execution admits a fixed Stinespring amplitude error ceiling
of 10^-8: the Frobenius success-map error plus the square root of all failed
export weights must lie below it. This controls arbitrary untouched
references for the sixteen-column map. The prepared observer is separately
checked end to end, including its entire terminal law. These are numerical
admission gates, not outward rounding or hardware error certificates.
Initialization, memory and detection are declared apparatus resources;
their physical costs and an SI time or energy conversion are uncalibrated.

```
python -m model.exchange_sector.many_body_packet --check
python -m pytest -q -W error model/exchange_sector -k many_body
```

Replay owns a fresh isolated interpreter and a fresh bytecode-cache prefix,
checks the complete source manifest before and after execution, and compares
every exact receipt leaf. The current result completes scalar occupation
statistics and compatible finite continuation. Native localized protection,
physical spin, interactions, a matched spectrum and the physical particle
dictionary remain the owning CP-03/CP-05 obligations; neither gate is closed.
