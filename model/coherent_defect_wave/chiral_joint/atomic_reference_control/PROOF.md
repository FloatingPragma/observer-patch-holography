# Atomic reference control on the complete Peter-Weyl source

## Source, physical space and the local reference

The source is the [complete Peter-Weyl Hamiltonian](../peter_weyl_links/PROOF.md)
on a periodic cube with even side n, spacing pi, thirty CAR modes and a
five-state scalar at every cell, and a complete 119-state link on every
oriented edge. All native onsite Yukawa coefficients and the scalar penalty
equal one. The source carrier is `2^(30V)*119^(3V)*5^V`, where `V=n^3`.
Every occupation is retained. Physical states satisfy all twelve actual
local Gauss generators. No empirical parameter enters the construction.

The [exact local reference](../local_reference/PROOF.md) is the complete
zero-charge onsite ground doublet, of spin one-half, raw L=1 and B3=0.
Its energy E0 is the smallest real root of

`[E(E-1)]^3-313[E(E-1)]^2+31708E(E-1)-1041468=0`.

It lies near -11.230585. Its next physical onsite energy has a certified
gap greater than 0.089555. The local ground density is `P_g/2`. It is a
quantum many-fermion/scalar state; each native family has exactly half its
available mean occupation, for total mean native occupation fifteen.

That local result imposes zero local Gauss charge. An excited lattice link
can compensate charged local matter. Its ground theorem alone therefore
does not bound the lattice's charged local sectors. The unrestricted bound
below supplies that input.

## Every unconstrained onsite sector

The actual full-CAR particle-hole transformation gives

`H_local = [[0,A^dagger],[A,I]]`

in scalar-blank and scalar-excited blocks. The quark transitions generate
Sp4 on six spinor copies; the lepton transitions act on all sixty-four
lepton Fock states. The complete quark occupation character accounts for
all `2^24` states and all collective Casimirs. The unique C=54 module has
dimension 140 at transformed quark number twelve, with multiplicity one.
Every other module has `C<=95/2`. No color-singlet occupation filter is
used in this statement.

The complete lepton transfer obeys `||A_l||^2<=6`, and the collective
quark inequality is `||A_q||^2<=2C`. Consequently every nonmaximal module
obeys

`||A_q+A_l||^2 <= (sqrt95+sqrt6)^2 = 101+2sqrt570 < 595/4`.

The strict comparison follows by squaring the positive rational
`(595/4-101)/2=191/8`, whose square exceeds 570.

For the maximal module, retain all `140*64*5=44800` states. On the complete
8960-state scalar-blank part, set `S=A^dagger A`. It splits into 1686
unrestricted weight blocks, of dimension at most thirty, labeled by
transformed lepton number, twice spin projection, twice weak weight and
twice hypercharge. All nonzero weak and hypercharge weights are included.
No Gauss kernel or highest-weight projection is taken.

The code forms S from all four actual collective transitions and their
metric adjoints. It rejects any outgoing coordinate outside its declared
weight block. In every block with exact Gram G it checks `GS=S^T G` and
computes all rational LDL pivots of both G and `(595/4)G-GS`. Every pivot
is strictly positive. The complete lists are in the receipt. Thus

`S < (595/4) I`

on the entire maximal module. Its 35840 scalar-excited states and every
scalar-dark state are retained in the exact scalar-star Hamiltonian.
The equivalent 140-dimensional module is justified by the complete native
character and actual Sp4 intertwiners, rather than a variational restriction.

Combining both parts covers all `5*2^30` unconstrained onsite states.
The scalar-star energy formula gives

`H_local > (1-sqrt596)/2 > -1171/100`.

The final rational comparison is checked by
`(-1171/100)(-1271/100)>595/4`. An exact unconstrained ground energy is
not identified; the rigorous lower bound is sufficient here.

## A volume-independent atomic gap at unit link penalty

Define the atomic Hamiltonian

`H0=sum_x H_local,x + kappa sum_e P_excited,e`.

It contains no hopping. Every link excitation projector commutes with the
endpoint currents and with H0. Fix any set of m excited links, and let a
be the number of incident cells. Each edge has two endpoints, so `a<=2m`.
At every other cell all incident links are blank. Its Gauss equations then
reduce to the local physical matter/scalar equations, giving energy at
least E0. At an incident cell the unrestricted bound gives energy above
`L0=-1171/100`, irrespective of the link's compensating charges or
entanglement. Put `delta=E0_upper-L0`, using the exact certified upper
endpoint for E0. On the complete m-link physical sector,

`H0 >= V E0 + m(kappa-2delta)`.

The proof uses an operator lower bound at each cell; it does not assume
product states in charged sectors. Every link count, native flux block and
fermion occupation is retained. The exact graph incidence check uses the
actual periodic edges.

At kappa=1 the quantity `1-2delta` exceeds 0.04117. The all-blank sector
is the tensor product of complete local physical spaces, with first gap
at least the local physical gap. Hence

`g=min(g_local,1-2delta)>0.04117`

is a volume-independent lower bound on the full physical atomic gap.
Its exact ground projector is

`P0=(tensor_x P_g,x) tensor |all links blank><all links blank|`,

with energy `V E0`, rank `2^V` and stationary density `P0/2^V`.
Both side two and side four are executed. The proof applies to every
declared even finite side. This is a theorem for H0; the original hopping
Hamiltonian requires the separate controls below. The atomic gap is not
identified with a particle mass.

## Full hopping norm and zero first-order compression

Let T be the complete actual Peter-Weyl kinetic operator. Every native
internal transporter is a contraction on its full internal/link carrier;
its exact Gram projection identity is replayed on all 119 link states.
Distinct-edge path products are contractions. The direct sum over native
blocks, tensored with the unit-norm Pauli matrix, is a contraction on the
thirty-dimensional spin/internal index.

Write a kinetic monomial without its coefficient as

`C_left^dagger V_path C_right`.

On the full Fock space, each CAR row or column has norm at most sqrt30,
because its squared norm is an endpoint number operator bounded by thirty.
Thus the monomial norm is at most thirty, including every occupied state.
Adding its adjoint and every actual path gives

`||T|| <= tau = 60 sum_paths |coefficient_path|`.

At side two this is `tau=360`. The onsite word norm is at most 57.
These are norm improvements on the same operator, with no source change.

All 59 native forward coordinates and their adjoints have zero blank-to-
blank entry. Every kinetic path is nonempty and has distinct edges.
Therefore a path acting on all links blank either vanishes or leaves an
excited link. It follows exactly that `P0 T P0=0`, independently of the
many-fermion occupations in the actual onsite ground states.

## Controlled finite-volume hopping and complete residuals

Consider the explicitly declared comparison

`H_lambda=H0+lambda T`.

All operators in this section act on the complete physical lattice Gauss
kernel. In particular `Q=I_physical-P0`; the atomic physical gap is not
asserted on the unconstrained lattice. Set `v=|lambda| tau` and require

`v<=g/5`, equivalently `|lambda|<=g/(5tau)`.

At side two this bound is about `2.29*10^-5`. Unit hopping is outside it.
The continued low-energy band has rank `2^V`, energies within v of `V E0`,
and distance at least `g-2v` from the remaining spectrum.

Relative to `V E0`, write the full blocks as

`H_lambda-V E0 = [[0,B^dagger],[B,D]]`,

where `B=lambda QTP0` and `D=Q(H_lambda-V E0)Q>=g-v`. Define

`X0=0`, `X_(j+1)=-D^(-1)B+D^(-1)X_j B^dagger X_j`.

No state of Q is omitted from D or its inverse. On the ball
`||X||<=r=2v/(g-v)<=1/2`, this map has contraction factor
`q=4v^2/(g-v)^2<=1/4`, preserves the ball, and converges to the exact
invariant graph X with `||X-X_j||<=r q^j`. By continuation from lambda=0
and the band separation, that graph spans the complete continued band.

Normalize `J_j=(I;X_j)(I+X_j^dagger X_j)^(-1/2)` and use the same-H
generator `K_j=J_j^dagger H_lambda J_j`. If

`h=V(57+3kappa+12)+v >= ||H_lambda-V E0||`,

then the exact graph isometry J and every iterate satisfy

`||J_j-J|| <= (1+2r) r q^j`,

`||H_lambda J_j-J_j K_j|| <= 4h(1+2r) r q^j = epsilon_j`.

The inverse-square-root integral bound and expansion of the full residual
give these inequalities. Duhamel gives state and complete-record total
variation error at most `min(1,|t|epsilon_j)`. The bound includes arbitrary
retained reference extensions and every outcome of a common complete
instrument. The bare atomic subspace separately has residual at most v
and record error at most `min(1,|t|v)`.

The normalized projector on the whole continued band, `P_lambda/2^V`, is
stationary. It is not asserted to be the ground projector inside a split
band. The approximate density `J_j J_j^dagger/2^V` is within the stated
isometry bound in trace distance. A spectator reference extends P0 and
these maps by its identity; fixing one particular purification is not used
to manufacture a gap in the reference complement. The admitted lambda
range may shrink with volume and is not a thermodynamic stability theorem.

## Actual local ground-spin coefficients

Let S be the spin-one-half operator within the exact local ground doublet.
Gauge covariance, native irreducibility and spin covariance imply

`P_g c_(s,a)^dagger c_(t,b) P_g`
` = delta_ab [delta_st I/2 + (beta_f/d_f) sum_i sigma_i,ts S_i]`

for internal indices in a native family f of dimension d_f. The scalar
coefficient follows from the exact half-occupation result. The quark
collective module has multiplicity one and its commuting quark-spin
representation is trivial. Each quark-family spin vector therefore has
zero compression. Only beta_ec and beta_L contribute, and
`beta_L=1-beta_ec` because the scalar carries no spin.

The code computes beta_ec from the actual native `S_ec,z`, after the exact
particle-hole transformation, orthogonal compression in the complete
physical highest-spin block and its rank-one ground spectral projector.
It uses a square-free characteristic polynomial, its exact projector
numerator and a rational isolating interval for E0. The receipt records
the resulting rational function `beta_ec=N(E0)/D(E0)` and interval. In
particular beta_ec is negative, about -0.1214. No naive transformed lepton
occupation is substituted for the physical spin operator.

Put `chi=beta_ec^2+(1-beta_ec)^2/2`; its certified interval is near 0.6435
and lies strictly above one-third. The actual blank-link Gram matrix for
all 118 forward/adjoint native channels is diagonal with entries `1/d_f`.
Exact normal ordering of 472 complete CAR words for each axis proves

`P0 T_e^2 P0 = [15 I-4chi S_x dot R_axis S_y]/16`,

where `R_axis` is the pi rotation about the edge axis, with diagonal +1
along that axis and -1 on the other two axes. This is a nonconstant
second moment about the actual local quantum ground. It alone is not an
energy-denominator-weighted exchange coefficient.

## Exact energy-denominator-weighted exchange

At side two every path is a single distinct edge. `T_e P0` lies in the
sector with exactly that link excited, and H0 preserves each such sector.
Consequently all cross-edge terms vanish in the complete second-order
operator

`K2=-P0 T Q [Q(H0-V E0)Q]^(-1) Q T P0`.

For one edge the denominator is bounded below by `1-2delta>0`.
Its inverse is the norm-convergent Laplace integral on the complete
physical excited sector. All other cells stay in their ground spaces;
on the two endpoints the integral contains the full unconstrained onsite
semigroup `E(t)=exp[-t(H_local-E0)]` and the link factor `exp(-t)`.
No charged local intermediate state is deleted.

For one native internal index define

`Fminus_st=P_g c_s^dagger E(t)c_t P_g`,

`Fplus_st=P_g c_s E(t)c_t^dagger P_g`.

Gauge covariance supplies the internal delta. Spin covariance gives

`Fminus_st=a_f(t) delta_st I+b_f(t) sum_i sigma_i,ts S_i`.

The two diagonal ground-spin expectations of `Fminus_ss` are real because
E(t) is Hermitian; they fix real a_f and b_f. The actual antiunitary obeys
`A c_s A^(-1)=i c_s^dagger`, preserves H_local and P_g, and sends S to -S.
Its ground action is `i sigma_y` followed by complex conjugation. Applying
it gives, with the index transpose retained,

`Fplus_st=a_f(t) delta_st I-b_f(t) sum_i sigma_i,st S_i`.

All thirty native CAR transforms, three spin transforms and the finite
response-tensor identities are executed. The two orientations contract as

`sum_stuv conjugate(sigma_st) sigma_uv Fplus_su(x) Fminus_tv(y)`

plus

`sum_stuv sigma_st conjugate(sigma_uv) Fminus_su(x) Fplus_tv(y)`.

Exact matrix multiplication gives
`4[a_f(t)^2 I-b_f(t)^2 S_x dot R_axis S_y]`.
The actual blank-link Gram contributes `1/d_f`, and the two independent
internal sums give `d_f^2`, leaving the factor d_f. Therefore the complete
edge term is

`K2_e=-A_edge I+J_edge S_x dot R_axis S_y`,

`A_edge=4|c_e|^2 sum_f d_f integral_0^infinity exp(-t) a_f(t)^2 dt`,

`J_edge=4|c_e|^2 sum_f d_f integral_0^infinity exp(-t) b_f(t)^2 dt`.

Both scalar and spin terms retain all five native families and every time
in the complete charged denominator. Here `|c_e|^2=1/16` on each of the
24 oriented edges. Uniform onsite data make A_edge and J_edge the same
on every edge. The nonnegative square integrand is a consequence of the
actual antiunitary and real tensor coefficients, not a sign assumption.

A rational positive lower bound needs only the L contribution near t=0.
The full local centered norm is at most 69. Since
`b_L(0)=(1-beta_ec)/2`,

`|b_L(t)-b_L(0)| <= 2(exp(69t)-1)`.

For `0<=t<=t0=1/1104`, the series bound
`exp(x)-1<=x/(1-x)` gives a deviation at most `2/15`. If beta_hi is the
certified upper endpoint of beta_ec, put

`b_min=(1-beta_hi)/2-2/15>0`.

Using `exp(-t)>=1-t0` on this interval gives

`J_edge >= t0(1-t0)b_min^2/2 = J_lo >0`.

Its exact value is

`14364888670332257029952361432834743 /`
`173816193821374085199143154070899916800`.

Every other channel and the rest of the integral is retained as a
nonnegative contribution. An upper bound follows from
`||T_e P0||^2<=(15+3chi_upper)/16` and the complete denominator gap.
The receipt supplies J's finite interval and retains the scalar coefficient
with its own bound. This resolves the energy denominator at the actual
atomic ground energy, rather than at a freely selected off-shell value.

The exact site rotations `O_x=R_1^x1 R_2^x2 R_3^x3` obey
`O_x^T R_axis O_y=I` on each edge of the even cube. Their spin lifts give
a unitary change of the retained spin basis. At side two the spin term is
therefore `2J_edge` times the antiferromagnetic Heisenberg operator on the
twelve undirected cube edges, plus the scalar shift `-24A_edge`. This is a
finite effective-spin identity, with no imported phase or particle mass.

## Controlled nonconstant interacting band

Let `D0=Q(H0-V E0)Q`, and let K be the canonical exact graph Hamiltonian.
The graph equation and resolvent identity give

`||X+D0^(-1)B|| <= v^2/[g(g-v)] + r^2 v/(g-v)`.

Graph normalization changes `B^dagger X` by at most `v r^3`.
Consequently, throughout the admitted interval,

`||K-V E0-lambda^2 K2||`
` <= v^3/[g(g-v)] + 12v^4/(g-v)^3 <= 6v^3/g^2`.

The last rational coefficient is at most `95/16<6` for `v/g<=1/5`.
This is a full-Q remainder, not a remainder after projection to selected
virtual states. On the side-two cube, the K2 expectation in an all-up
atomic spin product and that in the same product with the origin spin
flipped differ by exactly J_edge: the six incident oriented edges give
two contributions of each axis.

The explicit positive coefficient

`lambda_split=min(g/(5tau),g^2 J_lo/(24tau^3))`

therefore gives an exact canonical dressed-band expectation contrast at
least `lambda_split^2 J_lo/2>0`. The interacting band is nonconstant at
this controlled coefficient. Its lowest eigenprojector and full internal
spin spectrum are not identified by that contrast alone. Unit hopping is
not transferred into this perturbative interval.

## Complete unit-hopping quench records

A separate experiment uses the original `H=H0+T`, with unit hopping.
Purify each local ground mixture using two gauge-neutral fermionic reference
modes of spin one-half and raw L=-1. Retain their complete four-state Fock
space and zero reference Hamiltonian. The source/reference ground singlet
has even total parity and total raw L=0. A finite controller swaps it with
the full local vacuum and acts as identity on the complement; both states
have the same Gauss, spin, L and parity labels. The source density after
tracing these retained references is exactly `P0/2^8`.

Read the up-mode occupations of the references at cells zero and one.
Each reference effect is `n_up` or `I-n_up`, complete on all four Fock
states. On the prepared one-particle subspace these distinguish its two
spin states. All four records have probability one-quarter and are
retained. Their corresponding source ground spins are opposite. The other
six reference systems remain unmeasured and retained.

The link detector is the selected z edge excited with every other link
blank, and its full complement. It is gauge invariant. For aligned and
antialigned endpoint ground spins the exact first link moments are

`m_aligned=(15-chi)/16`, `m_antialigned=(15+chi)/16`.

Use `tau=360`, `h=936>=||H-8E0||` and
`t=1/(2048 h tau)`. Since `(H-8E0)P0=TP0`, the complete first-order state
remainder is at most `r_t=h tau t^2/2`. The selected edge's first-order
amplitude has norm at most `15t`. Its conditional probability error is
therefore at most

`eta=30t r_t+r_t^2`.

The actual antialigned-minus-aligned contrast is at least
`t^2 chi_lower/8-2eta>0`. Each original joint success interval includes
the reference probability one-quarter; its failure interval is one-quarter
minus that interval with reversed endpoints. All eight reference/detector
outcomes, all scalar states, every reference occupation and the entire
source complement are retained. The total probability enclosure contains
one. This finite-time theorem does not treat the atomic reference as a
stationary state of the unit-hopping source.

## Exact replay and scientific boundary

The receipt binds the complete local-reference and Peter-Weyl inputs,
their transitive code closure, this package's code and proof, and the
actual exact local algebra. Check mode independently recomputes every
record without writing and uses no numerical eigensolver. Tests reject
indefinite metric/bound proposals, insufficient link penalties and hopping
coefficients outside the admitted interval. The atomic ground, continued
band, resolved exchange and unit-hopping quench retain their different
Hamiltonians and domains. No physical mass, relativistic field or
thermodynamic hopping phase is inferred from these finite controls.
