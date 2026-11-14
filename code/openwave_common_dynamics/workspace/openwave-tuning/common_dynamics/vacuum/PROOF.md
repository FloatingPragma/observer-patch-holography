# A finite Dirac reference with one shared Gauss field

The source is the explicit fermionic lift of the routed four-component
Dirac proposal, on its eight-cell periodic member. Fermionic statistics,
the mass coin and the principal one-body phase lift are declared choices.
The field has one Z16 register on each of the actual288carrier edges.
The construction below gives native preparations and complete finite
occupation/field instruments. It does not identify an interacting vacuum,
a physical electron or a continuum particle pole.

## The free reference and its excitation records

On the zero-character field fibre, either direction along an axis reaches
the same opposite cell on the size-two torus. The two spin projectors sum
to the identity. Consequently the actual three-axis source is

`U = T_antipodal exp(-i mu beta)`, with `mu=h m`.

For `0<mu<pi`, each antipodal cell pair and each spinor component has two
orthogonal orbitals. Components with beta=+1 have negative principal energy
`-(pi-mu)/h` on the antisymmetric orbital and positive energy `mu/h` on the
symmetric orbital. Components with beta=-1 exchange symmetric and
antisymmetric roles; their negative energy is `-mu/h` and their positive
energy is `(pi-mu)/h`. The executable checks all32columns against the
actual routed source matrix.

Fill all16negative orbitals. There is exactly one occupied orbital in
each of16disjoint two-mode pairs. The reference is their exterior product,
with every fermionic ordering sign retained. Its complete occupation law
has65536outcomes, each with probability1/65536. The signed coefficient of
every outcome is computed and checked under the actual antipodal
permutation. Its mass phases cancel and its full macro eigenphase is one.
Each root port has mean occupation2; this mean does not assert a definite
local occupation. The source retains all local fluctuations.

The chosen energy is

`K=-(1/h) Arg(U)`, `H_F=dGamma(K)-E_sea`.

The sea energy is `-8pi/h`. The normal-ordered finite Fock Hamiltonian is
positive, with particle and hole costs equal to the absolute one-body
energies. The lowest cost is `m`; a neutral particle-hole excitation can
cost `2m`. There are eight low-energy states of each charge in this finite
member. The zero and pi Floquet valleys are retained. These degeneracies
and the selected logarithm branch do not identify one observed species.
The discrete source alone does not select this additive energy lift.

A particle excitation fills the unoccupied orbital of a chosen two-mode
pair. A hole excitation empties its occupied negative orbital. Their
complete occupation laws have32768outcomes each. A particle and hole in
distinct pairs have16384outcomes. All four laws retain every configuration;
no spectral or occupation outcome is discarded or normalized by success.

## A Gauss isometry on the full occupation space

Let B have +1 at the oriented source endpoint and -1 at the target. Reuse
the actual tree particular fields P and primitive cycle matrix C, which
obey

`B P = delta0 1^T-I`, `B C=0`, `C[chords]=I193`.

For the reference set G=2 at each of the eight root ports and zero at the
other ports. For an occupation configuration n with total16, define

`F(n,z)=P(n-G)+C z`, `J_g|n>=|n> sum_z g(z)|F(n,z)>`.

The equation `n+B F=G` holds modulo16 for every basis component. The
identity on chord coordinates makes z recoverable from F; the map is an
isometry, including arbitrary matter-reference entanglement. Every
particle shares this one field. It is not a tensor product of independent
particle fields. The zero-character choice uses the normalized constant
g on all193finite cycle coordinates. The Gaussian choice uses a product
of normalized positive cycle amplitudes.

A single extra particle has total17, and a single hole has total15.
Neither fits the same closed-torus constraint with total G=16. Their
separate preparations use G+delta0 and G-delta0 respectively. A neutral
particle-hole pair has total16 and uses the unchanged reference. The
executable rejects an extra particle with the unchanged G. Comparing
charged-sector energies therefore also requires the stated compensating
reference; its field cost is not removed by a label.

For a native routed hop x to y with current j, the difference

`P_x+j-P_y = C h`

is checked on all48directed size-two routes. Thus the field coordinate
transforms as z to z+h. A constant finite fibre is invariant under every
such shift, so the Gauss isometry intertwines the free dressed Fock
transport with its zero-character CAR lift. Onsite spinor rotations
preserve local number and require no field correction. The comparison
retains the same reference and every complementary occupation sector.

The existing encoded one-particle Dirac compiler is not a32mode Fock
compiler. The preparation here uses its own allocation and native
fermionic Gaussian words. A complete full-Fock evolution and its error
budget are separate from this initialization theorem.

## Finite-energy cycle states and the transfer error

For Z16 labels with centered value z in -8,...,7, use exact integer
amplitudes

`a_b(z)=b^(64-z^2) (b-1)^(z^2)`.

After normalization this is a finite discrete Gaussian. Its translated
overlap is computed as an exact rational number,

`r_b(k)=sum_z a_b(z)a_b(z+k)/sum_z a_b(z)^2`.

The product-fibre overlap for a route correction h is `product_j r_b(h_j)`.
The isometry defect of one axis is bounded by the largest branch norm
`sqrt(2-2 product_j r_b(h_j))`. Orthogonal rail projectors and the route
permutation justify this maximum for a one-particle input, including an
arbitrary untouched reference. Sum the three axis bounds for one macro
and telescope for a finite number of macros. A full Fock comparison must
account for every transported rail current; the one-particle bound is not
reused unchanged at charge16.

The executable reports five widths and exact overlaps. Increasing the
width improves free transport but increases electric energy. For the
base16state, six named occupation outcomes have field-energy expectations
between1266.9and1268.8. These floating expectations evaluate the complete
finite diagonal law by cyclic convolution of independent cycle
coordinates; they are not certified continuum energies. Every finite
state has energy at most9216 under the stated centered-square observable.

The constant fibre makes each edge uniform: every edge belongs to a cycle
with coefficient one. Its exact electric expectation is

`288 (m_field^2+2)/24`,

which equals3096atmodulus16 and diverges quadratically as the field
alphabet increases. It is not a finite-energy infinite-field vacuum.
An electric kick also destroys the constant-fibre code. The fields F=0
and F=A0 are in the same zero-charge fibre but have electric energies0and
3/2. The kick gives different phases for generic alpha. For alpha=1/16
their phase contrast exceeds0.09. Plaquette kicks alone act as the retained
scalar `exp(i sum_f beta_f)` on the zero-character fibre; that scalar is
not dropped from a controlled phase experiment.

These two facts prevent transferring the free Dirac sea to the coupled
electric/plaquette process by its Gauss label. The Gaussian preparation
is an explicit finite state with known free-transfer error, not an
interacting eigenstate. Positive-energy vacuum and field-limit theorems
for a different Wilson-family source are not imported into this one.

`energy_records.py` computes the complete sea's electric expectation,
including all65536occupation outcomes and all193shared cycle coordinates.
For each edge it convolves the independent cycle distributions and the16
independent diagonal Slater choices. Integer arithmetic with160outward
fractional bits encloses every convolution and final energy. These are
certified finite-record expectations, separate from the conditional
configuration examples and from physical rest energy. The compiled base16
target has expectation1264.8057854194997, enclosed in an interval of width
less than3.1times10^-43. The uniform comparison has exact expectation3096.

## A native test of electric leakage from the reference code

Deleting any two carrier edges leaves the actual graph connected. The
executable checks this by deleting each edge in turn and verifying that
the remaining graph is connected and bridgeless. Therefore any two edge
values can be fixed independently on a uniform affine Gauss fibre: the
remaining graph solves the adjusted Gauss equation through an integer
spanning tree. This argument works over the composite alphabet Z16 and
does not divide by a nonunit. Each one-edge marginal is uniform and each
two-edge marginal is a product of uniforms.

For `H_E=sum_e centered(F_e)^2/2`, every occupation-conditioned uniform
fibre consequently has exact mean3096 and variance27234. Let J denote
the uniform-fibre isometry. After an electric kick, the probability of
leaving its range is the occupation-weighted average of

`1-|E exp(-i alpha H_E)|^2`.

For independent copies H and H' with that conditional distribution,
this equals `E[1-cos(alpha(H-H'))]`. The inequality
`1-cos x >= x^2/2-x^4/24`, together with
`|H-H'|<=R=9216`, gives

`p_out >= alpha^2 Var(H) [1-alpha^2 R^2/12]`.

At `alpha=1/4096`, the bracket is37/64, giving a positive lower bound
greater than9/10000. It is uniform over matter states and their untouched
references in the code. This is a finite source intervention; the field
coefficient differs from the larger diagnostic coupling used elsewhere.

The actual native test prepares the Slater sea and uniform cycle state,
applies the electric word, reverses the shared Gauss encoder, applies
Hadamards to the772chord qubits and reads all1152field bits. Any nonzero
field result witnesses the code complement. The precise word, reset and
pointer error is subtracted from the ideal lower bound in ELECTRIC_CODE.json.
The ideal lower bound is0.0009384546428918839; after the complete native
error the conclusive nonzero-field probability is at least0.0009264166969.
The1538828-pulse native word is streamed and hashed, with all failure
outcomes retained. The result tests the chosen reference's response to its own
electric generator; it does not certify a replacement interacting vacuum.

## Actual native preparation and readout

The source allocation is32fermionic occupation qubits,1152field qubits
and3clean scratch qubits. The1187source/work qubits are reset by the actual
bounded native reset instrument. A seed occupation pattern followed by16
native Jordan-Wigner Givens words prepares the negative-band Slater state.
The same procedure admits the stated particle, hole and neutral-pair
seeds. The actual spectator strings and scalar phases are retained.

Prepare the four-bit Gaussian word on every chord register. Set the tree
reference offsets, add the exact cycle matrix C from chords to tree
registers, and add P times every occupation bit. This implements precisely
`F=P(n-G)+Cz`. The complete coefficient matrices and each controlled
increment's clean and dirty truth columns are checked. The controlled
arithmetic returns its scratch. All these blocks are reversible; the
factorization specifies the full unitary action without enumerating a
2^1187array.

The amplitude target uses exact integer weights. The reused outward Gray
certificate checks every rotation and native response, rather than treating
a floating residual as an error bound. The complete native word is streamed
and hashed. Its total error adds all Gaussian, Slater and arithmetic block
bounds. The1026974-pulse preparation has source error at most1.515944times
10^-10. The terminal instrument reads all1184data bits with the actual
native pointer; its auxiliary pair starts at1187/1188. No clock or spectral
probe is included in that allocation.

Every reset timeout, export complement and fine pointer record is an
outcome. The actual reset success Kraus words are scalar `|target><x|`,
so the success effect is diagonal in the full native input basis. Its
minimum eigenvalue bounds each conditional success even for coherent or
correlated inputs. The union bound needs no independence assumption.
The reported complete data-record bound adds these failures, the native
state-preparation error and every pointer error; it is at most1.220644times
10^-5. Private reset labels are
retained by the actual instrument; comparison with a blank tensor the
original reference marginal discards those private labels at that named
comparison boundary. It does not postselect successful runs.
