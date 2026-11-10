# Exact physical neutral odd response of the native quantum cell

The source is the full thirty-mode fermion Fock space tensor the five-state
quantum scalar, of dimension `5*2^30`. Its Hamiltonian is the sum of all
three actual native quantum Yukawa interactions and the scalar excitation
projector, with every coefficient one. The code consumes
`family_cell.column`, whose integer action is compiled from the full
periodic source's local Yukawa words. Dense occupation and scalar-column
controls independently compare this action with the actual native compiler.
No link, intersite hopping or scalar expectation is added to this local
calculation.

## 1. A physical two-component seed

Let Omega be the empty local fermion/scalar preparation and set

`psi = D_up^dagger Omega`,

`D_up=(phi_0 c_up,11 - phi_1 c_up,10)/sqrt(2)`.

The two nonzero components of `sqrt(2) psi` are `(mask=2048,scalar=1)`
with coefficient plus one and `(mask=1024,scalar=2)` with coefficient
minus one. The relative sign is verified from the actual quantum
composite operator. Replacing it with a plus sign fails a native weak
Gauss constraint. Both terms are retained throughout the calculation.

All twelve full normal-ordered gauge currents annihilate psi and commute
with H. All three local spin generators commute with H. The seed has
`S_z psi=psi/2` and `S_plus psi=0`, so `S^2 psi=3 psi/4`. Every vector
in its H-cyclic subspace is therefore a physical gauge singlet of local
spin one-half. This is a statement about the actual local SU(2) action;
it does not identify a relativistic little group.

The [physical composite preparation](../neutral_mobility/PROOF.md)
provides an even, gauge-invariant preparation gate with a retained
fermionic reference. Its two spin branches have the same local energy
distribution by spin rotation. The source Omega is a declared input and
is not a stationary ground of the full quantum Yukawa Hamiltonian.
The present seed is `D^dagger Omega`. It is distinct from applying D to
the interacting vacuum-connected ground; no latter spectral function is
claimed.

## 2. Complete graph and exact finite intertwiner

Closing both seed basis states under every nonzero full Hamiltonian
column gives 113948 basis states and 1111588 matrix entries. Every
column is reconstructed from the same full native source. The resulting
graph is invariant and connected, and its matrix is Hermitian. Its
largest distance from the two seed basis states is fourteen.

A diagonal phase S with entries plus or minus one makes every off-diagonal
matrix entry nonpositive. The physical seed has positive coefficients
on both of its coordinates in this same phase convention. Perron-Frobenius
therefore gives a unique component ground with strictly positive overlap
with the physical seed. Since spectral functions of H preserve the Gauss
and spin constraints, this ground is physical and has spin one-half.
The argument does not require other invariant components to have larger
energies.

The supplied equitable partition has 403 classes. For each class c the
inclusion sends its coordinate to the signed sum of all basis vectors in
that class. The code checks the equality of signed row sums on every
one of the 113948 full basis columns. Writing I for that inclusion,

`H I = I Q`, `I^dagger I = G`, `G Q = Q^transpose G`,

where G is the positive diagonal matrix of class sizes. The two-component
seed is exactly I v with `v^transpose G v=2`. This is an invariant
intertwiner, rather than a compression with an omitted full residual.

There are 104 scalar-blank classes and 299 excited classes. In this
ordering,

`Q = [[0,B],[C,I]]`, `A=C B`.

The physical seed lies in the excited block. A is self-adjoint in the
excited-block metric and has a 25-dimensional cyclic space on this seed.
Finite-field pivot rows establish exact independence of the proposed
columns. Exact integer/rational arithmetic then proves closure and the
full intertwining equation on every quotient coordinate. The finite-field
step alone is not used to infer an exact closure relation.

## 3. Minimal polynomials and all energy probabilities

The squared cyclic minimal polynomial is

`p(x)=x (x-105)(x-77)(x-53)(x-33)(x-17)(x-6)(x-5)`

`     *(x^2-19x+22)`

`     *(x^3-313x^2+31708x-1041468)`

`     *(x^3-229x^2+16756x-392988)`

`     *(x^3-157x^2+7684x-117948)`

`     *(x^3-97x^2+2764x-23388)`

`     *(x^3-49x^2+556x-1308)`.

Every root is real, nonnegative and simple. From the exact normalized
moments `m_k=<psi,A^k psi>`, the code constructs the unique degree-24
resolvent numerator n satisfying

`<psi,(z-A)^(-1)psi> = n(z)/p(z)`.

The receipt gives n and the moments explicitly. The polynomials are
coprime, and the weight of each squared root x is
`w(x)=n(x)/p'(x)>0`. For every positive root,

`E_plus/minus(x)=(1 +/- sqrt(1+4x))/2`,

`probability(E_plus/minus)=w(x)*(1 +/- 1/sqrt(1+4x))/2`.

The sign in the probability is that of an excited-scalar input. Using
the blank-scalar splitting formula would exchange these probabilities.
The root x=0 contributes only the energy-one atom, because the seed is
in the excited block; it contributes no energy-zero atom. There are
therefore exactly 49 visible energy records.

The full H-cyclic minimal polynomial is

`p_H(E)=p(E(E-1))/E`.

The code checks `p_H(Q)v=0` on all 403 quotient coordinates and exact
independence of 49 full-H Krylov columns. This yields the same complete
physical H-cyclic space without interpreting its nonorthogonal coordinates
as Euclidean ones. The source metric is retained in every moment.

Root isolation, rational interval arithmetic and outward dyadic square-root
bounds give every energy and probability interval in the receipt. All
49 probability lower bounds are strictly positive, and their interval
sums enclose one. The numerical eigensolver used for an exploratory
proposal supplies no certificate input.

A complete energy instrument consists of the 49 spectral projectors in
this physical cyclic space and its orthogonal complement on the entire
`5*2^30` carrier. The latter is an explicitly retained zero-probability
outcome for this preparation. Selecting the lowest visible energy gives
a physical ground-preparation filter with its exact success interval and
the retained complementary failure interval. It is not a claim about
successful preparation from every source state.

## 4. Comparison with the vacuum-connected minimum

Let xi be the largest root of

`x^3-313x^2+31708x-1041468`.

The odd component ground is `E_odd=(1-sqrt(1+4xi))/2`, approximately
`-11.23058529424765`. Its certified interval is

`[-1580564367072699/140737488355328, -3161128734145397/281474976710656]`.

The independently recomputed [vacuum-connected source](../family_cell/PROOF.md)
has minimum approximately `-11.14102996`. Exact outward subtraction
places `E_odd-E_vac` in

`[-12603792280123/140737488355328, -25207584559989/281474976710656]`,

strictly below zero. Consequently the vacuum-connected minimum is not
the global local minimum when all conserved source sectors are admitted.
This comparison does not determine the global minimum among every other
sector. It also does not compare two states in a fixed raw-charge sector,
or imply a negative particle mass or a physical instability of the model.

## 5. Centered charge and the actual conjugation partner

Define the raw occupation charge

`L=sum_spin (-n_ec + n_L10 + n_L11)`.

This is a conserved source generator. Its full one-body trace is two,
so its origin must be stated before interpreting the sector number.
The code verifies the antiunitary A defined by

`A c_i A^(-1)=i c_i^dagger`,

and by exchange of the scalar particle and antiparticle states. On the
complete thirty-mode operator algebra,

`A H A^(-1)=H`, `A G_a A^(-1)=-G_a`,

`A S_i A^(-1)=-S_i`, `A L A^(-1)=2-L`.

Thus the centered generator `L_centered=L-I` changes sign under A. Every
basis state in the odd graph has raw L=1, centered L=0, three times raw
baryon number zero, and `2S_z=1`. Raw L=1 is not identified here with
an observed lepton charge relative to a physical vacuum.

For a two-cell reference assigned raw L=2, raw sectors L=1 and L=3 have
relative charge deviations minus one and plus one. This is a consistent
choice of charge origin; it transfers no ground, energy or probe response
between different cell constructions or source couplings.

For an ordered occupation mask m, the full Fock action is

`A|m> = (-i)^popcount(m) (-1)^(sum occupied mode indices) |complement(m)>`,

with coefficient conjugation and the scalar exchange included. It obeys
`A^2=(-1)^N`. On the odd component its square is minus identity. The code
checks every graph label and this square, and verifies that A followed
by a pi spin rotation maps the whole 113948-state component bijectively
to itself. The unique up-spin component ground is therefore mapped by
A to its down-spin partner within raw L=1. There is no additional
observed particle/antiparticle species identification in this statement.

## 6. Requirements for a many-cell reference

A vacuum-connected spectral minimum alone is insufficient to establish
a global local gap for a product reference. A claimed product ground
must specify the admitted conserved sectors and compare their minima.
This calculation supplies one exact physical odd sector and a strict
comparison; it does not exhaust the full local source.

The spacing between visible energy lines does not bound an excitation
orthogonal to this cyclic input. A quantitative gap for the full
113948-state component needs a certificate covering those dark states,
in addition to the comparisons between different conserved sectors.

If a spin-one-half local ground is used, a rotation-invariant reference
can require correlations between cells. For two uncoupled cells, the
singlet of their up/down odd component grounds is a physical even-parity
state with energy `2 E_odd` and centered total charge zero. This is an
available finite reference construction, not the ground of the interacting
periodic Hamiltonian. Intersite interactions, state selection and a
volume-independent control must be calculated on their declared source.
The finite-range composite hopping and the long quantum paths have
different locality estimates; their effects are not imported by matching
onsite coefficient names.

The [receipt](results.json) binds the actual source code, both exact
partitions and the physical operator derivation. Missing declared or
transitive inputs are rejected. Partition source mismatch, failure of
the complete equitable identity and an incorrect relative seed sign are
separate rejected controls.
