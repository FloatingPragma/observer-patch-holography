# Spatial transport of a protected flux pair

This packet declares a new source, `coherent-mobile-paired-flux-v1`.
It changes the allowed link domain and allows a carrier's boundary data to
move between spatial cells. It does not assert that the previous rule, which
fixed every carrier's regional odd charge, already performed this transport.

## 1. A vacuum and an actual topological object

Let `s = 6553634` be the vertex-5 star. The old ten-word invariant domain is
enlarged by the two flat words `0,s`. All twelve words are literal assignments
to the same thirty icosahedron edges. The face holonomy is the parity of its
three edge bits. Both new words have zero nontrivial face holonomies. Each
old word has exactly two: one inside the chosen five-face region, one outside.
Consequently the regional parity `Q` is zero on the new orbit and one on the
five old orbits, whereas the full sphere defect count is `D=2Q`.

The old interior transport maps preserve both new flat words. The same
link-conditioned carrier wave `W_g`, parity-controlled transports `T_e` and
their inverses therefore define unitary operators on this larger domain.
The local one-wave dimension is 86,400, not 72,000. No old trace, reference
constant or mass inference is transferred to this new source automatically.
The full CAR lift uses the same 7,200 wave modes per cell. Its empty-wave
sector is stationary under the onsite update, for every allowed link word.

## 2. Gauss frame and two valid exchange extensions

Write a link word as `(r,sigma)`, with `r=min(g,g xor s)` and
`sigma=0,1`. The local Gauss involution flips sigma and multiplies the wave
occupation by `chi=(-1)^(number of doublet modes at port 5)`.
The diagonal unitary `F|r,sigma,n> = chi(n)^sigma |r,sigma,n>` obeys
`F Gamma F=X_sigma`. Thus the positive constraint is exactly an arbitrary
physical orbit/wave state tensored with the gauge-pair state `|+>`.

The dressed complete link-word SWAP
`S_full=F_x F_y SWAP_(r,sigma) F_x F_y` conjugates `Gamma_x` into `Gamma_y`.
It preserves the common positive subspace, although it need not commute
with each constraint separately. A stronger extension swaps only `r` while
leaving each sigma in its cell, then dresses by F. This operator `S` commutes
with both constraints separately. The two extensions agree exactly on the
common positive subspace: both gauge-pair qubits are `|+>`.

Both S operators are Hermitian involutions. With
`a=(1-i)/2`, `b=(1+i)/2`, the exchange `E=a I+b S` is unitary,
`E^2=S`, and its inverse is obtained by conjugating a,b. E leaves all CAR
occupations in their cells. It preserves `Q_x+Q_y` but not `Q_x` individually.
It moves the entire physical orbit, including the exterior compensator and
boundary bits. This is transport of a globally neutral flux pair; it is not
an unpaired monopole crossing a boundary that was held fixed.

The exact implementation uses the stronger, individually covariant S.
Tests also check equality with the full-word extension on positive states.
Onsite U_F at both cells followed by E is one declared clock tick. Its inverse
reverses the gates. An empty-wave, flat-orbit product is a stationary vacuum.

## 3. Preparation, complete position records and wave response

A reversible preparation swaps the orbits represented by 0 and 1, leaving
the other four odd orbits fixed. In the literal edge register it flips edge
bit 0 precisely on these two orbits. It commutes with Gauss and creates the
inside/exterior pair from the flat state. This preparation is an intervention;
the idle update conserves the total number of pairs and does not create one
from the empty vacuum.

Preparation from definite words by the two Gauss projectors has four flag
records, each of probability 1/4. Only the `++` branch is the positive physical
preparation. Every other flag is retained as a preparation failure, rather
than discarded from the probability account. Normalized probabilities below
are conditional on success; unconditional successful records acquire a factor
1/4.

For one prepared pair in cell x and a flat cell y, E gives probabilities
1/2 at x and 1/2 at y. The complete charge PVM retains all four records
`00,01,10,11`, including the zero-probability outcomes. A second E transports
the pair with probability one. The state remains coherent before readout.

The receiver experiment additionally prepares one wave at cell y, uniformly
over addresses and over ports 0,1, with internal index 3 and blank repair bit.
It executes E followed by the actual onsite U_F at y. Charge at y and the
wave-basis pointer `(port=0,internal=3,repair=0,address=0)` commute. Their joint
complete four records are computed from that same output state. They test
that moved link holonomy affects the native wave update, rather than merely
renaming a particle position. The receipt contains the exact probabilities.

## 4. Spatial extension and optional energy convention

On a periodic cubic array, use the same onsite layer and E on disjoint
even/odd matching edges of each coordinate. Every layer is unitary, preserves
all positive constraints and total Q, and has finite propagation radius.
The fixed matching clock is covariant under two-cell translations. This
algebraic extension includes all occupation sectors; no tensor factor is
cloned or reset. Intercell CAR hopping can be included as another declared
Gauss-compatible layer; it is not required for flux mobility here.

In the empty-wave positive code each cell has six physical states: one flat
state and five locations of its inside flux. The exchange is exactly the
ordinary partial SWAP of these six states. As a separate continuous clock one
may declare `H_flux=sum_x D_x+(pi/4)sum_edges(I-S_xy)`. Its terms are positive,
and in the one-pair sector its fivefold-degenerate exact band is

`E(k)=2+(pi/2)sum_i(1-cos(k_i))`.

Here D is an actual face-defect count, but its unit coefficient and the
Hamiltonian clock are construction choices. This quadratic lattice band is
not a relativistic rest mass or a GeV prediction. The separate spinor-routing
construction states the extra architecture needed for a controlled free
relativistic one-object limit.

## Boundary of the result

The packet supplies a vacuum, protected objects that really change spatial
cell, complete charge-position readout, and a finite coupling to native waves
on a shared link register. It does not establish a Standard Model identity,
fermionic statistics of a link object, interacting relativistic scaling,
measured masses, or a unique physical vacuum. Those statements require their
own premises; in particular a new routing coin is an input, not an inferred
spin of the old topological register.
