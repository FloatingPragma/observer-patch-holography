# A nonzero finite-window two-packet collision

## Source, input and complete experiment

Use exactly the paired replicated-block Hamiltonian `H_(a,N)` from the
[quadrature construction](../interacting_block_quadrature/PROOF.md).
The [block source and preparation](../interacting_block_propagation/PROOF.md)
provide the immutable native input bounds. Each constituent is the actual interacting
eight-cell source, with the same normalized native triplet matrix units.
Every physical complementary sector is retained. No coefficient or clock
is changed for this collision experiment. Set Delta to one only when
printing dimensionless calculations; actual times satisfy `Delta*T=1`.

The input consists of two distinct native triplet components, alpha and
beta, in finite free normal-mode wave packets. It is projected onto the
actual local-capacity sector and normalized. The output packet is projected
and orthogonalized against that input. These are explicitly prepared,
nonstationary physical source states. They are not substituted for the
actual interacting ground. Their source coefficients are determined by
the finite Hamiltonian, spacing and packet menu alone.

The two outcomes are the exact output-packet projector and its full
complement. Preparation failures carry separate retained records. The
proof bounds all state-preparation and detector differences at the
amplitude scale before giving a nonzero probability. In particular, no
`O(N^-1/2)` interacting-ground approximation enters the collision bound.

## Free channel and the exact interaction vertex

Write `A(k)=2 sum_i cos(k_i)` and
`omega(k)^2=D^2-2Dt A(k)`. In real normal modes,

```
b_k=u_k beta_k+v_k beta_(-k)^dag,
u_k=(sqrt(D/omega_k)+sqrt(omega_k/D))/2,
v_k=(sqrt(D/omega_k)-sqrt(omega_k/D))/2.
```

The free Hamiltonian acts independently on each triplet component and
preserves its total momentum separately. Thus an input with alpha momenta
in I_alpha and beta momenta in I_beta has identically zero amplitude for
an output whose corresponding momentum supports are disjoint. This is
an exact symmetry-forbidden channel for every free time, including the
Gaussian vacuum dressing.

The finite-capacity square root breaks the separate component momentum
symmetries while retaining total momentum. Expanding its exact matrix
elements to first order gives

```
H_N=H_free+V1/N+R_N,
V1=(t/2) sum_<xy>,gamma [Y_xgamma X_ygamma+X_xgamma Y_ygamma],
Y_xgamma=b_xgamma^dag n_x+n_x b_xgamma,
X_xgamma=b_xgamma+b_xgamma^dag.
```

The cross-component part for alpha different from beta is

```
(t/2) sum_<xy> [X_xalpha X_yalpha (n_xbeta+n_ybeta)
               +X_xbeta X_ybeta (n_xalpha+n_yalpha)].
```

Taking its matrix element between one particle of each component gives,
when `k+l=k'+l'` on the finite torus,

```
<k',l'|V1|k,l>=(t/(2V)) {
 s_k s_k' c_l,l' [A(k)+A(k')]
 +s_l s_l' c_k,k' [A(l)+A(l')] },
s_k=sqrt(D/omega_k),
c_k,k'=u_k u_k'+v_k v_k'.
```

It vanishes unless total momentum is conserved. Terms confined to one
component give zero in this channel because the other component's input
and output states are orthogonal. Vacuum contractions likewise cannot
change the separate component momenta. This is the connected quartic
vertex on the identified neutral modes.

For four momenta on the same cubic frequency shell, the coefficient is

```
M1=(t A(k)/V)[(D/omega_k)^2+1].
```

It is strictly positive whenever `t>0` and `A(k)>0`. It is not selected
using a measured mass or a desired collision probability.

The executable uses two finite Fourier packets. Let `u=2pi/L`, and put

```
I_alpha={(2u,u,0),(2u,-u,0)},  I_beta=-I_alpha,
O_alpha={(0,u,2u),(0,-u,2u)},  O_beta=-O_alpha.
```

Each packet has equal positive amplitudes `1/sqrt(2)`. The inputs carry
opposite x-directed carrier momenta with a transverse envelope; the
outputs carry opposite z-directed carrier momenta. Every mode has the
same exact frequency. There are six momentum-conserving combinations
among the sixteen four-mode combinations, each with amplitude weight
`1/4`. Thus the packet vertex is exactly `M_packet=3 M1/2>0`.
These are finite periodic packets overlapping in the observed region.
No asymptotic, spatially separated scattering matrix is claimed.

## Global weighted bounds without a fitted oscillator cutoff

For the proof only, extend the symmetric finite-capacity representation
to full bosonic Fock space by replacing each square root by
`f(u)=sqrt(max(1-u,0))`. Every local sector `n_x>N` is inert under the
quadrature exchange and retains onsite energy `D n_x`. The physical
sector `n_x<=N` is exactly invariant and has precisely the actual finite
source Hamiltonian. The extension introduces no physical oscillator.

Its global ground occurs in the physical sector. In an overflow sector,
replacing each inert site by its empty state removes a strictly positive
onsite cost while leaving all interaction expectations unchanged. The
same finite variational argument excludes the higher native block sectors.
The extension has compact resolvent because its exchange is bounded and
its onsite number energy diverges.

For every nonnegative u,

```
|f(u)-1|<=u,
|f(u)-1+u/2|<=u^2.
```

Consequently
`|f(u)f(v)-1+(u+v)/2|<=(u+v)^2`. On each creation/annihilation monomial,
`u+v<=n/N` and its matrix-element norm is bounded by `n+2`, with
`n=sum_x n_x`. Summing all four terms on every bond and all three
components gives the following conservative constants for a cubic graph:

```
C_delta=C_R=12 z V t,  z=6,
||(H_N_ext-H_free)psi|| <=(C_delta/N)||(n+2)^2 psi||,
||R_N psi|| <=(C_R/N^2)||(n+3)^3 psi||,
||(n+4)^2(H_N_ext-H_free)psi||
 <=(C_delta/N)||(n+6)^4 psi||.
```

These are global weighted bounds, including the auxiliary overflow
sectors. They do not require dropping a square-root boundary or assuming
that a Gaussian has compact occupation support.

Let `M=3V`. Free evolution transforms each canonical b operator into a
linear combination of b and its adjoint whose two coefficient norms sum
to at most `2/a`. Expanding `(n(t)+2)^2` into its degree-four, degree-two
and constant terms therefore gives

```
||(n+2)^2 U_free(s)psi|| <=G||(n+4)^2 psi||,
G=16 M^2/a^4+16 M/a^2+4.
```

This bound is uniform in s. All packet modes lie on one frequency shell,
so their free state changes only by a phase.

## Exact moment and second-order Born error

Let B bound `||(n+6)^4 psi||` for either free two-packet state. It is
computed without numerical tail fitting. In a real normal-mode basis,
the largest squeezed-vacuum occupation mean is at most
`u_max=(a+1/a-2)/4`. For two distinct-component one-particle excitations,
the total canonical number is stochastically bounded by `2+2X`, where X
has negative-binomial generating function

```
E[z^X]=[(1-q)/(1-q z)]^r,
r=(M+4)/2,  q=u_max/(1+u_max).
```

Cross terms between different real one-particle modes vanish for functions
of total canonical number by their distinct mode parities. Hence this
bound applies to the stated complex Fourier wave packets as well.
The domination follows, for example, by writing a negative-binomial
variable as a Poisson variable of mean `u Z`, with Z gamma distributed
with shape r and unit scale. Increasing u admits an increasing Poisson
coupling. The factorial moments are `E[(X)_j]=(r)_j u_max^j`. The executable uses
Stirling numbers to compute the exact rational eighth moment
`E[(8+2X)^8]`, and takes its integer square-root ceiling as B.

Duhamel's formula applied twice, with the global weighted bounds above,
gives

```
|<psi_out,U_N_ext(T)psi_in>
   +i (T M_packet/N) exp(-i E T)| <=R_B/N^2,
R_B=B[T C_R+(T^2/2) C_delta^2 G],
E=E_free+2omega.
```

The first term of R_B bounds the square-root remainder. The second
bounds two insertions of `H_N_ext-H_free`: the intervening free evolution
uses G, and the weighted second insertion uses the same B. All relevant
vectors have finite number moments, so the identities follow on the
finite-occupation core and extend by the displayed dominated bounds.
The finite Hamiltonian evolves exactly; only its comparison estimate is
expanded.

## Admission, preparation and the nonzero complete record

Let P_N project onto `n_x<=N` for every site. The moment bound gives

```
||(I-P_N)psi_in||, ||(I-P_N)psi_out|| <=tau=B/(N+7)^4.
```

Normalize the projected input. Normalize the projected output and then
subtract its component along the admitted input. If `tau<=1/4`, the sum
of the resulting two state-vector errors is at most `8tau`. Their exact
source inner product is zero, so the collision detector has zero
probability at time zero.

Let `M_lower>0` be a rigorous rational lower bound on M_packet. It uses
`cos x>=1-x^2/2` and the rational inequality `pi<22/7`. Choose the finite
integer N so that

```
R_B/N^2 <=T M_lower/(8N),
8B/(N+7)^4 <=T M_lower/(8N),
tau<=1/4.
```

The admitted physical amplitude is then at least
`3 T M_lower/(4N)`. Detector matrix elements are exact finite spectral
services and introduce no additional approximation.

Preparation uses a specified finite source unitary mapping the product
of native block grounds to the admitted input. A rank-two Householder
unitary supplies such a map; its coefficients are the finite projected
Gaussian coefficients and source-defined spectral vectors. It is not an
efficient implementation claim. A retained finite translation and
internal-character reference carries the inverse packet momentum and
native triplet character. Equivalently, translate the packet and its
position reference together and retain the reference. The clock during
the observation is exactly `H_(a,N)` with a zero-Hamiltonian reference.

For a translation-covariant implementation, use the finite position
reference `C[Z_L^3]` and replace each packet phi by
`V^(-1/2) sum_R T_R phi tensor |R>`. On the admitted even sector the
internal translation image is `Z_2^3`; its eight characters extend to
the corresponding even characters of the actual native cover `Z_4^3`.
The eight-character reference starts in its trivial basis state `|0>` and ends
in the distinct basis state `|-X_alpha-X_beta>`. Every reference basis
state has zero clock energy. The complete Householder input is the native
product ground tensor the uniform position state tensor `|0>`; its
target is the relational packet tensor `|-X_alpha-X_beta>`. Both joint
vectors have total character zero. The rank-two unitary acts within that
invariant sector. Translation
invariance of H makes its matrix element between the two relational
packets equal to the source matrix element computed above. The reference
is retained; it is not replaced by a classical symmetry-breaking vector.

The complete finite orientation reference has an equally direct
implementation. Let G be the finite group of the actual source's cube
and macro-graph symmetries, using the double cover if required by its
fermion representation. On its regular reference, with zero Hamiltonian,
define

```
J psi=|G|^(-1/2) sum_g U_g psi tensor |g>,
O_tilde=sum_g U_g O U_g^dag tensor |g><g|.
```

Then `O_tilde J=J O`, and `[H,U_g]=0` gives
`(H tensor I)J=J H`. This applies to the Householder preparation, target
projector and its full complement. Their joint transformations use the
left regular action on the reference. Every probability, vector norm
error and occupation tail is unchanged. It supplies an explicit retained
orientation reference for the directed packet menus. No classical frame
or unrecorded symmetry-breaking state is substituted.

The immutable quadrature parent supplies the actual native cover
`G_native=Z_4^3 semidirect 2O`, of order 3072, and its full native operator
replays. For the five-block cycle use `G_native times C5`, with the
second factor acting by ring translations. Its retained regular reference
has dimension 15360. The additional binary detector register gives
dimension 30720, beyond the separately retained native preparation
references. No three-dimensional macro cubic action is asserted for the
ring. The component projectors lifted by this reference are joint,
reference-relative records. Tracing out the orientation reference instead
gives the source group twirl: each bare component marginal is one third
of the sum of the three relative component marginals. The receipt gives
both tables. The native block-ground projector is invariant.

The native product ground can be supplied by the same complete finite
agreement-filter preparation as the block-propagation theorem. During
this explicitly declared preparation protocol, inter-block exchange is
gated off; the observation clock begins after the state preparation.
For `M_source=NV` constituent blocks, its frozen exact odds bounds give
product-ground infidelity at most `M_source r0 c^d` after depth d.
The initial state is the unfiltered atomic/cube trial u. Its local
vacuum-to-singlet swaps and cube Householder synthesis are deterministic,
with production success exactly one. Its overlap p0 with the interacting
block ground is a fidelity bound, not a preparation-success probability.
The depth d counts every optional agreement-filter round, including the
first. The first round has success at least `p0^2`; this factor is
contained in the j=0 term of the tree bound when d is positive. At d=0
no agreement filter is run, so no such factor occurs, while the full
unfiltered trial error remains in the fidelity bound.
Choose d with

```
2 M_source r0 c^d <=[T M_lower/(8N)]^2.
```

The retained purification then has source-reference vector error at most
`T M_lower/(8N)`. The same unitary maps that error to the admitted input,
and exact source evolution preserves it. Thus, conditional on all
preparation success flags, the final hit probability is at least

```
p_hit_given_success >=[T M_lower/(2N)]^2>0.
```

The unconditional probability includes the strictly positive preparation
success lower bound, which is recorded in factored form and by its
logarithm. Every failed filter, every retained reference and the entire
non-hit detector complement are part of the complete experiment. The
resource sizes and success rates are extreme; no efficient simulation or
experiment is asserted.

## Bounded continuum window and scientific scope

The checked members use physical box side eight, `a=1/4,L=32` and
`a=1/8,L=64`. Their packet components obey `|p_i|<=pi/2`, and
`Delta*T=1`. The frozen exact cosine remainder controls their free
Klein-Gordon dispersion error in this fixed momentum/time window. The
actual finite-N interaction is nonzero and the complete rare collision
record is bounded below after every approximation error.

The admitted state also has a controlled actual energy window. The
extension's variational ground comparison gives
`0<=E_N-E_free<=C_delta B/N`. Since P_N commutes with H_N_ext and the
free packet is an exact energy eigenstate,

```
||(H_N-E_N-2omega)phi_in|| <=3 C_delta B/N.
```

Therefore the exact admitted packet's probability outside an actual energy window of
half-width w about `E_N+2omega` is at most
`[3 C_delta B/(Nw)]^2`. For the physically prepared source with retained
purification error epsilon_p, the complete outside-window probability
is at most

```
min(1,[3 C_delta B/(Nw)+epsilon_p]^2),
epsilon_p<=T M_lower/(8N).
```

This follows from the norm of the actual energy-window projector, so it
does not assume an energy bound on the preparation error. The Hamiltonian
residual itself is asserted only for the exact admitted packet. This is
a separate energy concentration statement; no energy postselection is
used to establish the rare collision lower bound.

The result is a finite-time, nonstationary two-packet collision on the
same declared source as the free relativistic subsector. Its connected
amplitude is of order `T/N` and vanishes in the increasing-capacity limit.
It proves a nonzero controlled finite-window interaction without claiming
a nontrivial interacting relativistic ultraviolet limit, an asymptotic
scattering matrix, physical Lorentz spin, or Standard Model particle
identities.

## Direct full-capacity finite collision

A separate executable diagnostic uses five full eight-cell blocks on a
cycle, N=1, `D/Delta=2`, `t/Delta=1/8` and `Delta*T=1`. The complete
three-component low carrier has dimension `4^5=1024` and is exactly
reducing inside the full native physical source. Its macro graph has
degree two; it carries no three-dimensional continuum claim.

The full native source register has kinematic capacity
`5*(240 log(2)+8 log(5)+24 log(119))` nats. The numerical propagator
executes its exact invariant 1024-state carrier, of capacity `log(1024)`;
the full native matrix is not numerically enumerated. The complement is
controlled analytically by invariance and the full-source preparation
error. Separately retained registers contribute `40 log(4)` nats for
the native trial spin references, `log(15360)` for the regular reference,
and `log(2)` for the binary detector pointer. No filter register is used
at depth zero. The full regular reference includes the translation
characters, so no additional eight-character register is counted.

The input has a free alpha particle at `k=2pi/5` and a free beta particle
at `-k`. The output swaps their momenta. Both are projected onto actual
capacity, normalized, and the output is orthogonalized to the input.
At this small capacity the projection is substantial: its weight is
approximately 0.789394. No small-capacity Gaussian approximation is used.
The ideal unprojected channel is free-symmetry-forbidden, but this does
not hold for the admitted five-block input and detector. Indeed
`H_N1=P_N H_free P_N`, so their matched full-free initial transition
derivative equals the finite source's nonzero derivative. Its real
generator matrix element exceeds 0.0522218. Orthogonalization enforces
only the zero-time overlap. The finite time series is therefore a
momentum-exchange response; it is not labelled a channel that free
evolution cannot open on the same admitted states.

The actual 1024-state Hamiltonian has rational entries and norm at most
`95/8`. Its exponential is evaluated to Taylor order 100 with the explicit
operator remainder

```
exp(95/8) (95/8)^101 /101!.
```

All Gaussian coefficients, admissions, orthogonalization and propagation
are evaluated with outward interval arithmetic. The fifth-root phases
use exact nested square roots. The certified imaginary amplitude exceeds
`0.0263136461392` after the Taylor remainder. The receipt stores a strict
dyadic lower bound, rather than rounding that decimal upward.

Reserving half the amplitude for the complete native preparation gives
hit probability at least `0.000173101993285...`. The same frozen native
trial-overlap bound suffices at filter depth zero for this diagnostic,
so no filter success is suppressed in that lower bound. The initial
source is the explicit prepared packet, with its proven preparation
error, and the two detector outcomes exhaust the physical carrier.

The receipt also gives complete hit/complement intervals at
`Delta*t=0,1/4,1/2,3/4,1`, both for the ideal admitted input and including
the actual native trial preparation error. The five-block filter depth is
zero, so those preparation-inclusive intervals are also unconditional:
the deterministic trial is produced with probability one. At the final time the stronger
complete prepared hit interval is approximately
`[0.00296630042443,0.00308409427253]`. It includes the entire certified
preparation error, rather than only the deliberately loose factor-two
margin above. Ground and each triplet-component occupation marginals are
given per block with outward intervals and the same error propagation.
Translation symmetry makes the ideal one-block marginals uniform around
the ring. They are quantum expectation values and describe no classical
trajectories. Any higher-native-sector population caused by imperfect
preparation has its own retained upper bound.

## Matched full-free comparator for the five-block response

The comparator uses exactly the same admitted input and orthogonalized
target as the finite source. Between them it evolves with the full
quadratic Fock Hamiltonian. This is a mathematical comparison, not an
extra oscillator added to the native physical carrier.

For one component and Fourier mode m, put

```
c_m=cos(2pi*m/5), h_m=D-2t*c_m, g_m=-2t*c_m,
omega_m^2=h_m^2-g_m^2,
A_m(tau)=cos(omega_m*tau)-i h_m sin(omega_m*tau)/omega_m,
B_m(tau)=-i g_m sin(omega_m*tau)/omega_m.
```

The coherent-state propagator has the exact Gaussian kernel

```
K(z,w;tau)=C(tau) exp[(z.Z.z)/2+z.X.w+(w.Y.w)/2],
X_m=1/conj(A_m), Z_m=B_m/conj(A_m),
Y_m=-conj(B_m)/conj(A_m),
C(tau)=exp(i D*5*tau/2) product_m conj(A_m)^(-1/2).
```

The matrices in position coordinates are the finite Fourier transforms
of these mode functions. The square roots are the continuous branches
starting at one. On the reported positive times every `conj(A_m)` lies
in the upper half-plane, so the executable implements that branch
explicitly with positive real and imaginary square roots.

The kernel follows from the Heisenberg transformation
`U^dag b U=A b+B b^dag`: its intertwining equations fix X, Y and Z.
The Schrodinger equation fixes C with `C(0)=1`; the factor shown removes
the specified `-5D/2` quadrature constant. For zero pairing it reduces
exactly to the number-Hamiltonian kernel. Finite-mode Gaussian symmetries
and their coherent-state kernels are treated by
[John and Parthasarathy](https://arxiv.org/abs/1911.06555).

Both admitted vectors have finite Fock support. Their matrix element is
therefore a finite sum of derivatives of K, evaluated as sums over pair
contractions. The three components factor. This evaluates the full free
evolution without any intervening oscillator cutoff. Outward interval
arithmetic bounds every coefficient and contraction. An independent
number-basis matrix evolution also checks representative one-component
kernel entries; its numerical cutoff is a validation control and does
not enter the certified interval calculation.

For any initial purification within vector error epsilon_p of the exact
admitted input, a hit probability p lies in
`[max(0,sqrt(p_lower)-epsilon_p)^2,
  min(1,(sqrt(p_upper)+epsilon_p)^2)]` for either unitary experiment.
To define both experiments on the same imperfect preparation, decompose
the full native physical source as `H_low direct_sum H_rest`, where H_low
is the exactly reducing 1024-state carrier. Embed it isometrically into
`K=Fock(15 modes) direct_sum H_rest` by the canonical capacity embedding
on H_low and the identity on H_rest. Tensor this embedding with every
retained native preparation reference, symmetry reference and record.
The prepared purification has exactly the same vector-error bound after
this isometry.

Define the comparator on K as `U_free direct_sum U_rest`, where U_rest
is the actual native evolution restricted to H_rest. Any other fixed
unitary extension on that complement would obey the same error estimate.
The native experiment has its actual unitary on the embedded physical
carrier and the identity on the unused extra Fock complement. Both use
the target-packet hit projector, extended by zero on H_rest, and its
complete complement. The free comparator's higher oscillator occupations,
every native complementary sector, and all references are retained.
Its symmetry reference is the same regular lift. No native error state
is silently replaced by an oscillator vacuum or discarded.

This construction depends only on the common initial vector error and
the complete detector norm. The receipt subtracts the resulting intervals
to certify the native-minus-matched-free response. At `Delta*t=1` the
ideal matched-free hit probability is approximately 0.00223731693905;
both preparation-inclusive intervals are disjoint. The comparator has a
nonzero response, and the finite source has a larger response on this
specified experiment. No Lorentz-covariant scattering interpretation or
empirical collision prediction is inferred from the five-block graph.
