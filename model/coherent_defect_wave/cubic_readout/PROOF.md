# Three-dimensional records of the spatial occupation law

The construction uses the [shared-link CAR law](../species/PROOF.md) at
each cell and its Gauss-dressed one-body root-SWAP hopping. The cubic
placement, ordered matching layers and coherent apparatus are declared
architecture. The source tick is unchanged from that spatial occupation
construction. It is distinct from the tensor-wave exchange circuit and
from any continuous-time Hamiltonian extension.

## A physical encoding with all local link information

Let `s=6553634` be the star mask and choose the smaller representative `g`
of each pair `(g,g xor s)` in the ten-word link orbit. There are five such
representatives, corresponding to five regional face positions. For a wave
mode `w=(port,j,repair,address)`, let `chi(w)=-1` precisely when `port=5`
and `j>=3`, and otherwise let it be `+1`. On one wave the physical encoding is

\[
 J|g,w\rangle={|g,w\rangle+\chi(w)|g\mathbin{\mathrm{xor}}s,w\rangle\over\sqrt2}.
\]

The vectors on the right are an orthonormal basis of the positive local
Gauss sector. The analogous empty-wave vector has `chi=1`. At arbitrary
occupation the product of the occupied chi signs replaces the single sign.
This defines an isometry on the full local occupation algebra, retaining
every one of the five physical link-orbit values.

The parent source proves `[U_F,G_F]=0`, so `J* U_F J` is unitary on the
physical sector. On occupation one it is exactly `J* U J`. In the sparse
implementation, a normalized input pair has common norm denominator two.
After the literal parent tick, the coefficient on the canonical word is
the coefficient in the reduced basis: its partner has coefficient
`chi(w)` times that value, and the square-root factors cancel in `J*`.
`local_column` executes this construction using the owning U and its actual
inverse, rather than rebuilding an approximate local operator.

Conjugating by the parent's diagonal frame changes the constraint into a
pure link-pair flip. In that frame the spatial hopping is the exterior lift
of `a I+b SWAP`, where

\[
 a=(1-i)/2,\qquad b=(1+i)/2.
\]

Consequently the reduced one-wave hopping has these two coefficients and
leaves both local face labels fixed. Tests embed the reduced vectors into
the existing two-cell implementation, including a doublet mode with a
negative chi sign, and verify the exact onsite and hopping intertwiners.

## Finite cubic program and state

Use an even periodic cube of side L, with L at least four. A full tick V
first applies U_F at every site, then the six hopping matchings in order
`(axis,parity)=(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)`. Each matching pairs a
site whose coordinate has the specified parity with its positive-axis
neighbor. Every layer is unitary on the full finite CAR algebra. Reversing
the layers and local gates implements the inverse.

The complete law preserves total occupation. In its exactly invariant
one-wave sector, only the occupied cell has a nonempty onsite update. A
sparse basis key therefore contains the occupied cell, its wave mode and
the face values at all cells that differ from the chosen reference. A
missing face entry means the reference value; it does not trace out or
discard that cell's field register. As the wave visits cells, any number
of changed face values can be retained. The representation is exact in
this sector for arbitrarily many ticks; executing only one tick is a
resource choice for the displayed experiment.

The reference is the empty-wave product with face representative one at
every cell, in the positive Gauss sector. Empty occupation makes every
onsite and hopping gate the identity on this state. There are `5^(L^3)`
orthogonal stationary empty-wave basis states, so this reference is not
claimed to be the unique physical vacuum or a ground state of a selected
energy generator.

The experiment prepares a wave in mode `(7,3,0,0)` with equal amplitudes
on the eight sites `(0 or 1,0 or 1,0 or 1)`, divided by square root eight.
It uses the same reference link values. Pure coherent packet preparation
is an apparatus resource. Starting with definite link word one at every
cell, all local Gauss projections produce `2^(L^3)` retained flag patterns,
each with probability `2^(-L^3)`. The normalized positive-sector experiment
is conditional on the all-positive flag. Its complete original-trial law
multiplies every displayed detector probability by `2^(-L^3)` and retains
every other preparation flag. This small success probability is not omitted
from the resource claim.

Each spatial cell keeps its own odd regional flux throughout this program.
The permitted link repairs change the defect's face within that carrier;
the hopping gates keep both local link registers in their cells. A change
of carrier face therefore differs from transport of regional topological
charge between spatial cells. The wave occupation propagates across the
cubic array and can leave changed local face records. Identifying this
coupled excitation with a mobile protected physical particle requires its
own observable attachment. Carrier-face position and cubic-cell wave
position are distinct readouts.

## Operational translations and momenta

Translate the complete cells, including their link values and CAR modes.
The full Fock translation includes the fermionic mode-reordering sign.
It restricts to an ordinary position/field permutation at occupation one.
Translation by two cells preserves every ordered matching and the identical
onsite layer, so it commutes with V. It permutes the local Gauss constraints
and preserves their common positive sector. Translation of the wave alone
does not have this property when its changed field is left behind. On the
executed L=4 state the squared difference between the complete and wave-only
translated outputs is exactly `1/225`.

Set `M=L/2` and let `T_d` translate by `2d` cells, for `d in Z_M^3`. The
complete commuting translation projectors are

\[
 P_q={1\over M^3}\sum_{d\in\mathbb Z_M^3}e^{iq\cdot d}T_d,
 \qquad q_i={2\pi k_i\over M},\quad k_i=0,\ldots,M-1.
\]

They are orthogonal, sum to the identity and commute with V. The source
packet and its M-cubed translates have disjoint occupation support, hence
every q has conditional probability `1/M^3`. The labels q are dimensionless
translation labels in two-cell units. Attaching a spatial unit a would give
`p_i=q_i/(2a)` modulo the reciprocal supercell period, but neither a physical
momentum scale nor a physical position/energy dictionary is inferred here.

For L=4, the eight translation orbits are evaluated using an additional
exact compression. Each basis configuration has a unique representative
whose occupied coordinates are all zero or one. Its orbit has eight members
because a nonzero two-cell translation changes the occupied position. The
vectors

\[
 |[b];q\rangle={1\over\sqrt8}\sum_d e^{iq\cdot d}|T_d b\rangle
\]

are an orthonormal momentum-orbit basis. `compress` sums the corresponding
signed coefficients and multiplies the norm denominator by eight.
`expand` restores all eight translated configurations. Tests compare this
encoding with the direct eight-term projector on the actual full-tick
output. No momentum, link configuration or detector outcome is sampled away.

## Complete frequency instrument

The apparatus chooses two quadratures with equal probability, retains that
choice and reads both signs. Its four Kraus operators are

\[
 K_{b,\epsilon}={I+\epsilon i^b V\over2\sqrt2},
 \qquad b\in\{0,1\},\quad\epsilon\in\{+1,-1\}.
\]

For fixed b, summing `K* K` over the signs gives `I/2`; all four sum to I.
These are two complementary two-outcome phase instruments with retained
quadrature flags. They are neither four spectral eigenprojectors nor a
length-four Fourier phase estimation circuit. Every translation projector
commutes with these Kraus operators. The complete joint operators
`K_(b,epsilon) P_q` therefore retain all momentum and quadrature outcomes,
also when an untouched quantum reference is present.

Let `c(q)=<psi_q,V psi_q>` for the normalized projected packet. Its joint
probabilities conditional on positive Gauss preparation are exactly

\[
 p(q,b,\epsilon)={1+\epsilon\operatorname{Re}(i^b c(q))\over4M^3}.
\]

They are obtained from the actual source autocorrelations
`<T_d psi,V psi>`, rather than from a guessed dispersion. In the L=4
execution, all 32 branch norms are evaluated in the orthonormal orbit basis.
For L=8, 256 probabilities use the same proved kernel and the executed
autocorrelations; those full branch vectors are not materialized.

## Exact three-dimensional response

The L=4 and L=8 full ticks each execute 9344 nonzero basis amplitudes, have
norm one, and move at least one cell's defect away from its reference face
with probability `1/450`. The L=4 inverse returns the exact input vector.
The diagonal matrix element of the source wave mode is `-29/30`: the
address reflection contributes `2/60-1`, and at address `(0,1)` port seven
is outside both the seam repair and the controlled interior edge endpoints.

Only the same wave mode and unchanged reference links contribute to the
first moment against the prepared packet. For each axis its equal even/odd
sublattice vector is fixed by the even matching; the expectation of the odd
matching at coarse momentum q is `a+b cos(q)`. The three spatial factors
commute at occupation one. Consequently, for every allowed q and even L,

\[
 \boxed{\quad c(q)=-{29\over30}\prod_{i=1}^3
 \left({1-i\over2}+{1+i\over2}\cos q_i\right).\quad}
\]

The receipt checks this formula against the actual tick at all eight L=4
momenta and all 64 L=8 momenta. This is an exact expectation of V, not a
dispersion relation for an eigenphase or a particle pole.

For L=4, let r be the number of pi components of q. Then
`c(q)=(-29/30)(-i)^r`. The conditional phase records at each momentum are
permutations of `(1/120,59/120,1/4,1/4)`, and each momentum has weight `1/8`.
The squared unitary fluctuation is `1-|c(q)|^2=59/900`, proving that none
of these prepared packets is a frequency eigenvector. On L=8 it is at
least `59/900` at every tested momentum.

The L=4 first moment also gives a spectral support bound. If `z=e^(i phi)`
is the V eigenvalue and `phi0=arg c(q)`, then the conditional expectation
of `cos(phi-phi0)` is `29/30`. Outside angular distance pi/3 it is at most
one half, while everywhere it is at most one. Thus at least `14/15` of
the conditional spectral weight lies within pi/3 of phi0. This broad
finite-frequency bound is not an isolated particle pole or a positive
physical mass gap.

Two controls delimit the interpretation. Removing all onsite evolution
makes the tested L=4 projected packets exact hopping eigenstates and gives
zero defect movement. Removing only the five occupation-controlled link
moves, while retaining the native wave and hopping, also gives zero defect
movement but preserves these first-moment phase records exactly. Therefore
this frequency instrument alone cannot certify the interaction or determine
the complete spectrum. The defect readout separates the two programs.

## Controlled spatial comparison and physical boundary

One onsite layer has support in one cell; each of the six hopping layers
crosses at most one nearest-neighbor edge. The backward support of a local
record after t ticks grows by at most `2t` in each coordinate and `6t` in
Manhattan distance. The disjoint gate argument applies to the full even
CAR observable algebra and complete local instruments with arbitrary untouched
references. If preparations and gates agree throughout that cone, their
local record/reference laws agree exactly. Global momentum projectors do
not become local apparatus by this argument.

The receipt also executes one tick in side-eight and side-twelve cubes,
with packet origins `(2,2,2)` and `(4,4,4)`. After the corresponding spatial
identification, the entire sparse output vectors agree exactly; the remote
boundaries are outside the occupied causal cone. This fixed-spacing test
does not supply a zero-spacing continuum limit or uniform particle-pole
estimates.

The construction supplies three-dimensional simulation-generated records
with retained gauge fields, preparation failures and a common tick. It
supplies no chosen physical energy logarithm, Lorentz-spin identification,
chiral nonabelian currents, asymptotic scattering, rest mass or GeV value.
Resolving higher spectral information and the appropriate physical energy
definition is necessary before using these records for a particle mass.
