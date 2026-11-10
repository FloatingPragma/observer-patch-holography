# Contractive quantum links with proper cubic rotations

The finite construction uses the actual native fifteen-state current
representation, thirty fermion modes and five-state scalar memory at every
site. A link has 119 states. The Hamiltonian acts on every state of

`F(C^(30 N^3)) tensor (C^119)^(3 N^3) tensor (C^5)^(N^3)`

for even N at least two. Every source coefficient and the source clock are
specified below. This is a declared gauge-covariant quantum architecture,
with exact finite spatial symmetries and bounded link transport. The theorem
does not identify its excitation energies with measured particles or supply
an interacting relativistic limit.

## 1. Native representation and source ancestry

`physical_attachment.chiral.representation()` supplies twelve exact Hermitian
native currents T on the fifteen actual exterior-representation states.
Their common irreducible blocks, in the native index order, are:

| Native indices | Dimension d | Hypercharge |
| --- | --- | --- |
| 0, 1, 4 | 3 | -2/3 |
| 2, 3, 5, 6, 7, 8 | 6 | 1/6 |
| 9 | 1 | 1 |
| 10, 11 | 2 | -1/2 |
| 12, 13, 14 | 3 | 1/3 |

These are representation labels. The matter carrier contains one declared
native family; doubling the link flux orientations adds no conjugate matter
family. The actual source embedding of the thirty local fermion modes is
the [sitewise native current construction](../../chiral_source/PROOF.md).
The five-state scalar operators and local invariant Yukawa tensors are
consumed from [the quantum Higgs construction](../../chiral_higgs/PROOF.md).
The exact axial momentum coefficients are consumed from
`periodic_family.kinetic_terms`. No prepared state, spectrum or ground from
the sixty-state link Hamiltonian is used.

The link carrier and its operator amplitudes are declared architecture
choices. An implementation by the narrower native idle law is not proved.
The native matter-current ancestry does not turn this link choice into a
deduction from that idle law. The present package defines a distinct
Hamiltonian and clock `exp(-it H_oriented)`. Its energy differences require
their own calculation and physical reference.

## 2. Full link carrier and endpoint generators

For every native block r introduce orthonormal states `|+,a,b>` and
`|-,a,b>`, and one common blank state `|0>`. Thus

`dim K = 1 + 2 sum_r d_r^2 = 119`.

In a forward sector the endpoint generators are `(-T^transpose,T)`;
in a reverse sector they are `(T,-T^transpose)`. The two tensor factors
refer to a and b. They vanish on the blank state. Each is Hermitian,
the endpoints commute, and each endpoint represents the same native
Lie algebra. Define, for a and b in the same block,

`U_ab = (|+,a,b><0| + |0><-,a,b|)/sqrt(d)`.

These equations specify every matrix entry, including annihilation on the
reverse flux sector. The complete endpoint Ward identities are

`[L_0(T),U_ab] = -sum_c T_ac U_cb`,

`[L_1(T),U_ab] = +sum_c T_cb U_ac`.

For example, L_0 acts on `|+,a,b>` as `-sum_c T_ac |+,c,b>`.
Its right action on the reverse bra is
`-<-,a,b|L_0 = -sum_c T_ac <-,c,b|`. Both terms therefore have
the required left transformation. The right endpoint follows in the same
way. This verifies the convention rather than assigning an orientation
label without its actual current matrices.

For `Q_x(T)=c_x^dagger T c_x`, the canonical anticommutation relations give
`[Q_x,c_x,a^dagger]=sum_c c_x,c^dagger T_ca`. The two contributions cancel
in `sum_ab c_x,a^dagger U_ab c_y,b`, and the right endpoint cancels in the
same manner. At an internal path vertex, the right Ward term of one edge
cancels the left Ward term of the next edge. This is a normal-ordered
operator identity on the entire fermion Fock space, with arbitrary
spectator occupations. The local Yukawa and scalar-current identities
are those of the same consumed scalar matrices and invariant tensors.
Thus every local Gauss generator commutes with the full Hamiltonian.

`matrix_checks()` verifies all 1416 endpoint identities using the actual
twelve native 15-by-15 matrices and every column of each 119-by-119 link
matrix. The proof requires no vacuum restriction. The link energy is
`kappa P_excited`, which commutes with all endpoint generators.

## 3. A unitary orientation reversal

Define the linear permutation

`J|0>=|0>`, `J|+,a,b>=|-,b,a>`, `J|-,a,b>=|+,b,a>`.

On the complete link carrier,

`J^dagger=J`, `J^2=I`, `J L_0(T) J=L_1(T)`,

`J U_ab J=U_ba^dagger`.

The index transpose is necessary. A control that swaps sectors while
omitting this transpose fails an off-diagonal matrix identity. The
transformation exchanges geometric endpoints and their flux currents;
it does not conjugate the internal representation of a matter mode.

The raw raising-only sixty-state hardware has
`rank(sum_ab U_ab^dagger U_ab)=1` and
`rank(sum_ab U_ab U_ab^dagger)=59`. This excludes a unitary or antiunitary
that exchanges the complete raising family with its adjoint. For the
declared 119-state operators the corresponding sums are

`15 P_0 + sum_r P_(minus,r)/d_r`,

`15 P_0 + sum_r P_(plus,r)/d_r`.

Both have rank sixty and J exchanges them exactly. The construction
provides a concrete orientation reversal on the full finite basis.

## 4. Contractive transport on every link state

On the internal block tensor the complete link carrier, let

`V_r = sum_ab |a><b| tensor U_ab`.

Define the orthonormal vectors

`eta_b = sum_a |a> tensor |+,a,b>/sqrt(d)`,

`chi_a = sum_b |b> tensor |-,a,b>/sqrt(d)`.

Then

`V_r = sum_b |eta_b><b,0| + sum_a |a,0><chi_a|`,

`V_r^dagger V_r = I_d tensor P_0 + sum_a |chi_a><chi_a|`,

`V_r V_r^dagger = I_d tensor P_0 + sum_b |eta_b><eta_b|`.

Both products are rank-2d orthogonal projections. Consequently
`||V_r||=1` on all `119 d` internal/link states. Other irreducible link
sectors are retained and have zero action under this block transporter.
The normalization is an architectural coefficient. Omitting it gives
`||V_r||^2=d`, as the explicit dimension-six control verifies.

A product along any path is a product of contractions on the common
internal space and the retained link factors, so its norm is at most one.
For a path with distinct edges, a blank-link input saturates this bound:
there are `d^ell` orthogonal output terms for each fixed incoming internal
index, each with squared coefficient `d^(-ell)`. Thus the norm is exactly
one. The code executes every native incoming index at lengths one through
four, absorbing words in every native block, and a mixed absorbing/creating
two-edge word. No reverse flux state is discarded from the operator bound.

The matrices are quantum operators: `U_ab^3=0`, while `U_ab^2` is nonzero.
They are not a substitution of classical identity links.

## 5. Distinction from Peter-Weyl compression

The carrier has the vector-space dimensions of the blank representation,
the five native matrix-coefficient spaces and their conjugates. Those ten
nontrivial representations are inequivalent because their hypercharges
are distinct. This gives the 119-state carrier when those exact sectors
and their conjugates are retained; no global hardware minimality is claimed.

The declared U matrices are not literal compressions of multiplication by
group matrix elements. Such a compression includes every retained
Clebsch-Gordan transition, not only transitions incident on the vacuum.
The general representation-basis multiplication formula is given in
[Zohar and Burrello, *Building projected entangled pair states with a local
gauge symmetry*](https://doi.org/10.1088/1367-2630/18/4/043008).
The distinction here has an exact native witness.

Let e be native index 9, the one-dimensional hypercharge-one character,
and let L be indices 10 and 11, the hypercharge-minus-one-half doublet.
For `epsilon=[[0,1],[-1,0]]`, every actual current satisfies

`epsilon (T_ee I_2 + T_L) epsilon^(-1) = -T_L^transpose`.

Thus `e tensor L` is equivalent to the conjugate L representation.
Multiplication by the one-dimensional e character maps the complete
four-dimensional matrix-coefficient space of L isometrically into the
retained conjugate-L space. Its coefficient intertwiner has rank four,
represented by `epsilon tensor epsilon` up to the choice of dual bases.
The declared `U_ee` annihilates every forward L basis state, so its rank
on that input is zero. It therefore omits an actual retained nonvacuum
transition. The code checks this using all twelve currents, without
assuming the tensor-product identification from particle names.

The construction claims the explicit covariant quantum-link architecture
of sections 2 through 4. A literal Peter-Weyl implementation requires its
additional Clebsch transitions and a separate operator specification.
Quantum links as finite noncommuting link matrices are established in
[Chandrasekharan and Wiese](https://arxiv.org/abs/hep-lat/9609042); that
general construction supplies no continuum transfer for these operators.

## 6. Complete periodic Hamiltonian and proper cubic rotations

Let every spatial direction have N sites, spacing a, and antiperiodic
fermions. Gauge and scalar memories are periodic. With the uniform exact
real Yukawa coefficients `y_up,y_down,y_lepton`, nonnegative scalar and
link penalties, and positive a, define

`H = H_Y + mu sum_x P_scalar,x + kappa sum_e P_excited,e`

`    + sum_(x,i,ell=1..N-1) [c(x,i,ell) c_x^dagger sigma_i W_(x,i,ell)
                            c_(x+ell e_i) + adjoint]`,

where W is the actual ordered product of the 119-state link matrices and

`c(x,i,ell) = q(x,i,ell) [-i pi/(2 N a)] (-1)^ell / sin(pi ell/N)`,

`q(x,i,ell)=(-1)^floor((x_i+ell)/N)`.

The default a is pi and all five source coefficients are one. The one-link
kinetic coefficients are consequently `+i/4,-i/4`, before the declared
block normalization in U. `column` acts on every full key
`(fermion_mask,119_state_links,5_state_scalars)` and retains local quantum
Yukawa interactions. `gauss_column` uses the actual endpoint currents;
`translation_column` supplies the three exact antiperiodic translations.

For a proper signed coordinate permutation R, choose an SU(2) lift S with

`S sigma_i S^dagger=sum_j R_ji sigma_j`.

Write `R x = x_R + N w_R(x)` with `x_R` in the canonical cube. The fermion
one-particle rotation is

`|x,s,a> -> (-1)^sum_j w_R(x)_j sum_t S_ts |x_R,t,a>`.

The full fermion action is its exterior lift, retaining the permutation
signs of every occupied mode. The native index a is unchanged. Scalar
memories are permuted by `x -> x_R` and have spin zero. A positive edge
`(x,i)` maps to the positive edge at `x_R` if `R e_i=+e_j`; for
`R e_i=-e_j` it maps to the positive edge at `x_R-e_j`, with J acting on
its entire link state. This is a unitary on the complete finite carrier.

For an axis with positive image, the transformed hopping term is the
corresponding positive-path term. For an axis with negative image, each
U becomes the transposed adjoint on the reversed geometric edge. After
reversing the path, the result is the adjoint term with its endpoints
exchanged. The spin matrix contributes a minus sign. Since the scalar
path coefficient is purely imaginary, complex conjugation supplies the
same minus sign. The antiperiodic seam phases satisfy

`q(target positive path) = (-1)^(sum w_R(x)+sum w_R(y)) q(original path)`

for the original endpoints x,y, in both sign cases. This proves exact
kinetic invariance for every even N. Omitting the seam phase fails the
side-four control.

The local Yukawa pair is a spin singlet because
`S^transpose epsilon S=epsilon`. Uniform scalar and link penalties are
unchanged; J preserves the excited-link projector. Hence `U_R H U_R^dagger=H`
on the complete carrier. The two choices of S differ by fermion parity,
so the physical unitary action is the double cover of the 24-element
proper cubic group. A two-pi rotation is precisely total fermion parity.
This finite statement does not assert continuous rotation symmetry,
boost covariance or a relativistic little group.

The executable checks all 24 exact spin lifts and every complete path
word for N=2 and N=4, respectively 48 and 1152 words. Expanding a word
uses the full matrix identities of sections 2 and 3, so this check acts
on all internal indices and every link column. An additional direct
execution checks 33 Hamiltonian/pi-rotation commutators, including thirty
occupied, nonblank controls spanning all five native blocks and both flux
orientations. Their coefficients are unequal rational numbers, and their
scalar states are retained. Three actual Gauss commutators and three
translation commutators are also executed. These basis controls complement
the full operator proof and are not a claim to enumerate the enormous
periodic Fock basis.

## 7. Quantitative bounds and scope

On a blank-link one-fermion input, both forward and adjoint hopping terms
create flux, into orthogonal orientation sectors. The two orientations
therefore contribute equally. With local Yukawa and energy terms set to
zero for this control member,

`||H_hop psi_blank||^2 = 6 sum_(ell=1)^(N-1) |c_ell|^2`

`                         = pi^2 (N^2-1)/(2 N^2 a^2)`

by `sum_ell csc^2(pi ell/N)=(N^2-1)/3`. This is bounded by
`pi^2/(2a^2)`, uniformly in the volume at fixed spacing. Direct complete
Hamiltonian executions give `3/8` and `15/32` at a=pi for N=2 and N=4.
These are charged full-carrier operator controls; no physical preparation
or ground is inferred from them. The response differs from the normalized
raising-only architecture because the adjoint can create reverse flux.

Let `P_e=(-1)^N_excited,e`. Every U anticommutes with P_e, including both
of its terms. Thus the local dressed parity

`C_x=(-1)^N_matter,x product_(e incident to x) P_e`

commutes with this Hamiltonian as well. On the all-blank-link subspace,
each C_x reduces to local matter parity. Its exact compressed propagation
therefore cannot move a single odd matter occupation between cells while
returning every link to blank. Orientation reversal and contractive paths
do not remove this independent obstruction to that particular field
comparison.

The construction supplies exact Gauss covariance, three translations,
proper cubic rotations and volume-independent norms for individual path
transporters. The axial kinetic prescription has unbounded path range as
N grows; an individual path bound does not prove a uniform finite-range
Hamiltonian bound or a Lieb-Robinson estimate. There is no explicit magnetic
plaquette term. Interacting ground selection, a physical dressing and
clock, controlled refinement, generation multiplicity and the measured
particle dictionary require their own results on this Hamiltonian. No
sixty-state spectrum or GeV calibration is transferred by the matching
native labels or kinetic coefficients.

The [receipt](results.json) uses exact rational and algebraic values and
hashes its executable import closure. Its check command recomputes all
matrix and operator-word tests. Deliberately omitted transposes,
normalizations and antiperiodic seam phases are rejected by separate
controls.
