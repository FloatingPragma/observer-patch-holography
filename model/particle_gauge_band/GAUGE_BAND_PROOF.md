# GAUGE-COVARIANT COMPACT BANDS, CANONICAL OCCUPATION AND COMPLETE RECORD TRANSPORT

1. Source and domain

The spatial operator is the literal compact source on the radius-three
cubic lattice, with 343 sites, 882 undirected edges and two spin coordinates
per site. Its hopping from x to x+e_j is i sigma_j/2, its reverse is the
adjoint, and its onsite term is 64 n(x).sigma. The texture n is the normalized
source direction (x,y,z^2-1) inside its prescribed support and is constant
outside. Its amplitude vanishes only at c_-=(0,0,-1) and c_+=(0,0,1).
The source shell degrees are -1,+1 and zero on the exterior shell.

An internal representation of dimension d is tensored with this source.
The charged source has d=45, consisting of the actual fifteen-state
exterior representation and its three family coordinates. Spatial spin
and the internal representation act on separate factors.

For an oriented edge x->y, supply a unitary internal matrix U_xy, with
U_yx=U_xy^dagger. Its hopping block is

 h[U]_xy = (i sigma_j/2) tensor U_xy.

The onsite texture is unchanged. These classical parallel transports are
additional field data. Their relation to finite nonunitary native links is
addressed in QUANTUM_LINK_PROOF.md. No measured mass selects any U_xy.

2. A uniform gauge-covariant spatial band

Let g_x be a unitary internal transformation at each site, and set

 U_xy^g=g_x U_xy g_y^dagger,
 G=direct_sum_x (I_spin tensor g_x).

Then h[U^g]=G h[U] G^dagger exactly. The texture term commutes with G.
The complete operator-valued row and column norm sums of the hopping are
at most 3, since there are at most six neighbors and every link block has
norm 1/2. The block Schur test gives ||T[U]||<=3 independently of d and
of the link holonomies.

The onsite operator has precisely 4d zero eigenvalues. Every other
eigenvalue is +64 or -64. Ordered Weyl inequalities therefore give exactly
4d eigenvalues of h[U] in [-3,3], with the complete complement outside
(-61,61). Define the Riesz projector

 P_h[U]=(1/(2 pi i)) integral_{|z|=32} (z-h[U])^(-1) dz.

It has rank4d and obeys P_h[U^g]=G P_h[U] G^dagger. The separating gap is
at least58 for every classical unitary link field in this domain. This
proves covariance by changing the background and the projector together.
It does not require a fixed rank-four spatial projector to commute with
every site projector.

The localization estimate also extends directly. With P_C the projector
onto the two texture cores and Q_C=I-P_C, write

 h[U]=[[0,B^dagger],[B,D]].

For half-separation1, B^dagger B has eigenvalues5/4 and7/4, each with
multiplicity2d. Each core has six neighbors, giving3I/2 on each diagonal
block. The cores share one neighbor, whose cross block is one quarter
of a unitary matrix. Thus ||B||^2=7/4 for arbitrary unitary links.
The complementary operator D has |spectrum(D)|>=61. For each low
eigenvector, its Feshbach equation is

 E u = -B^dagger(D-E)^(-1) B u.

Since |E|<=3, this gives |E|<=7/232. A sharper root of
e(61-e)=7/4 is also available. In either form the exterior/core probability
ratio is bounded by (7/4)/(61-e)^2. Gauge transport preserves the texture
degrees but need not preserve the individual energies or their color
degeneracies. A curved background supplies physical data beyond a gauge
choice; the spectrum of a pure-gauge background is unchanged.

3. Pairing must enter the band definition

Let K_r be an actual symmetric all-left source kernel on one retained raw
record r. The source has thirteen distinct kernels: the disabled kernel
and twelve kernels obtained by inserting g P_j in the lepton family block.
The quark coefficients are the fixed six-power construction, with maximum
1/sqrt(2). Each P_j is an exact real rank-one projector. The source's
rational interval calculation gives

 0.39831556363708714589 <= g <= 0.39893988270234903798 < 1.

The neutral-Higgs tensors have disjoint support and unit nonzero singular
values. Consequently ||K_r||<=1 for every raw kernel. The execution checks
the exact projector and tensor identities and recomputes that gain interval.
This bound retains the declared, target-informed source architecture.

The creation pairing is Delta_r=C_spatial tensor conjugate(K_r), where
C_spatial is the source's staggered i sigma_2. For a varying scalar
background the corresponding onsite Delta transforms as described below;
the theorem requires ||Delta_r||<=1. Set

 B_r[U]=[[h[U],Delta_r],[Delta_r^dagger,-conjugate(h[U])]].

Delta_r^T=-Delta_r, so B_r is a valid Hermitian Bogoliubov-de Gennes
matrix. In Nambu coordinates Psi=(c,c^dagger), define

 Xi = tau_x followed by complex conjugation,
 tau_z = diag(I,-I).

Then Xi B_r Xi=-B_r. The onsite Nambu operator has8d zeros, with every
other eigenvalue at+64 or-64. The hopping and pairing perturbation has
norm at most4. Therefore the symmetric low projector

 Pi_r=(1/(2 pi i)) integral_{|z|=32} (z-B_r)^(-1) dz

has rank8d, its energies lie in[-4,4], and its complete complement lies
outside(-60,60). The gap is at least56. Also Xi Pi_r Xi=Pi_r.
For d=45 this is the full360-coordinate low Nambu band, with no discarded
neutral, color or family coordinates. The full field has61,740 Nambu rows.

The scalar spatial projector P_h generally fails to reduce the onsite
pairing on a curved link field. Defining Pi_r from B_r retains that cross
term in the spectral problem. It does not delete it by hand.

4. An intrinsic integer occupation on the complete band

The symmetric Nambu projector alone determines a canonical anticommutation
subalgebra, but an occupation operator also needs an annihilation/creation
splitting. Here that splitting follows from the source particle grading
and an explicit bound.

For the number-conserving matrix B_0=diag(h[U],-conjugate(h[U])), let
Pi_0=diag(P_h,conjugate(P_h)). On |z|=32, the two resolvent norms are
at most1/29 and1/28, respectively. The resolvent identity gives

 ||Pi_r-Pi_0|| <= 32/(29*28) = eta = 8/203.

Since Pi_0 commutes with tau_z,

 ||(I-Pi_r) tau_z Pi_r|| <= ||[tau_z,Pi_r-Pi_0]|| <= 2 eta.

Define C_r=Pi_r tau_z Pi_r. On range(Pi_r),

 C_r^2 = Pi_r-Pi_r tau_z(I-Pi_r)tau_z Pi_r
       >= (1-4 eta^2) Pi_r >0.

Its sign is therefore defined without a chosen eigenbasis:

 J_r=C_r(C_r^2|_{range Pi_r})^(-1/2), extended by zero on the complement.

It satisfies J_r^dagger=J_r, J_r^2=Pi_r and Xi J_r Xi=-J_r. Thus

 Q_r,+ = (Pi_r+J_r)/2,
 Q_r,- = (Pi_r-J_r)/2,
 Q_r,0 = I-Pi_r

are orthogonal projectors. Xi exchanges Q_r,+ and Q_r,-. Each has rank4d.
Choose any orthonormal frame F of Q_r,+ and put b=F^dagger Psi. Then

 {b_i,b_j^dagger}=delta_ij,
 {b_i,b_j}=0,

because F^dagger Xi F=0. The number operator is basis independent:

 N_r=sum_i b_i^dagger b_i = (1/2) Psi^dagger J_r Psi + (4d)/2.

Its spectrum is the integer set0,...,4d, with the usual multiplicities.
For the actual charged source, N_r has180 canonical modes and centre90.
Consequently

 V_r=lambda(N_r-90)^2, lambda>=0

is a positive, even interaction supported on that full band. The
quadratic B_r reduces the band and its complement; its pairing inside the
band is retained, so it need not conserve N_r. The resulting finite Fock
Hamiltonian has a low-factor interaction and a complementary quadratic
factor in a suitable Bogoliubov coordinate system. This is a static
classical-background statement.

5. Exact reduction to the literal source at trivial links

At U_xy=I, the spatial source obeys

 h C_spatial = -C_spatial conjugate(h).

Its symmetric rank-four spectral projector P_W is particle-hole compatible.
Therefore

 Pi_source=diag(P_W tensor I_45,conjugate(P_W) tensor I_45)

reduces every actual raw B_r. On this band the quadratic energies are
source spatial energies plus/minus kernel singular values, all inside
[-4,4]. The complementary quadratic energies have absolute value at least
62-1=61. The rank and spectral separation identify Pi_source with the
Riesz projector Pi_r exactly.

Pi_source commutes with tau_z. Hence C_r=tau_z Pi_source, J_r=C_r, and
the displayed N_r is precisely dGamma(P_W tensor I_45). The interaction
lambda(N_r-90)^2 equals the original compact source interaction. No
small-perturbation approximation is used in this reduction.

6. Local covariance of the complete interaction

For the transformation G above, let

 Delta_r^g=G Delta_r G^T,
 G_N=diag(G,conjugate(G)).

This is the transformation induced by transforming the Higgs coefficient
with its all-left invariant tensor. It is required when changing a local
frame. Holding the scalar coefficient fixed while transforming only the
fermions changes the physical background.

Then B_r^g=G_N B_r G_N^dagger. The same relation holds for Pi_r, C_r,
J_r and all three Q projectors, because G_N commutes with tau_z and
with Xi. The Fock implementation of G gives covariance of N_r and V_r.
Thus the compact interaction and the local gauge frame are compatible
on this classical-link branch.

These equations are exact for spatially local, time-independent gauge
transformations. For time-dependent frames an additional temporal
connection transforms as A_0^g=G_N A_0 G_N^dagger+i dot(G_N)G_N^dagger.
The spatial spectral projector is defined from B_r, excluding that temporal
connection. Gauge-dependent temporal potentials are not rest energies.

7. Complete-record polar transport

Two raw records r,s have the same link field and ||Delta_r-Delta_s||<=2.
The contour identity gives

 ||Pi_r-Pi_s|| <= xi = 64/28^2 = 4/49.

Every eigenvalue of C_r on its band has absolute value between
sqrt(1-4 eta^2) and1. Therefore

 ||J_r-C_r|| <= 1-sqrt(1-4 eta^2) <= beta=4 eta^2.

Also ||C_r-C_s||<=2xi. It follows that

 ||Q_r,+-Q_s,+||, ||Q_r,--Q_s,-|| <= (3/2)xi+beta,
 ||Q_r,0-Q_s,0|| <= xi.

Define the alignment operator and its positive polar part by

 A_rs=sum_{a in {+,-,0}} Q_r,a Q_s,a,
 U_rs=A_rs(A_rs^dagger A_rs)^(-1/2).

Since sum_a Q_s,a=I,

 ||A_rs-I|| <= 4xi+2beta
             =16/49+512/41209 <1.

A_rs is invertible. A_rs Q_s,a=Q_r,a A_rs and A_rs^dagger A_rs commutes
with every Q_s,a, so U_rs Q_s,a U_rs^dagger=Q_r,a. In particular U_rs
transports both the low band and its occupation grading exactly.

Xi exchanges the+ and- summands and fixes the0 summand; hence Xi U_rs Xi=U_rs.
The interpolation of A_rs to I stays invertible, so the induced real
Majorana orthogonal transformation lies in the identity component. It has
an even Gaussian Fock implementation, obtained from its continuous Spin
lift. Gauge covariance of its projectors gives covariance of U_rs.

A controller may apply this transformation when the recorded branch changes,
then evolve the same full band Hamiltonian and its quartic. Every disabled
and enabled record is retained with its original probability. The controller
must retain which preceding branch projector was used; the complete record
contains that information. This transport is an added apparatus operation,
with its own work and duration. It is not an identity of the unmodified
native memory channel, and no native pulse compilation is asserted here.

8. Changing link data and the tracking term

For differentiable link/scalar data, the projectors Q_a(t) are differentiable
under the same gap bound. Set

 K(t)=sum_a dot(Q_a)Q_a = (1/2)sum_a[dot(Q_a),Q_a].

Differentiating Q_a Q_b=delta_ab Q_a proves K^dagger=-K and
[K,Q_a]=dot(Q_a). The additional quadratic Hamiltonian B_track=iK
alone transports all three projectors, including their occupation
splitting. The full B+B_track transports the low band Pi, since [B,Pi]=0.
It generally does not parallel-transport Q_+ and Q_- individually: the
retained pairing can give [B,J]!=0, so the full dynamics need not conserve
the band occupation. B_track is particle-hole admissible because Xi K Xi=K.
With a temporal connection, replace dot(Q) by D_t Q=dot(Q)+i[A_0,Q].
Then B+B_track+A_0 has the required time-dependent covariance.

For the same interacting Hamiltonian, changing the band also changes the
interaction. Its exact derivative is

 dot(V)=lambda {N-90,dot(N)}.

The tracking term is an explicit control law. A physical repair process
would have to generate that control or supply a justified approximation
to it. Instantaneous spectral isolation alone does not prove autonomous
tracking, a stationary vacuum or a physical time calibration.

9. Executed controls

The producer constructs exact Gaussian-rational color, weak and central
unitaries with cosine3/5 and sine4/5. Four edges of one source plaquette
carry a nontrivial Wilson product. These inputs were specified without
empirical particle data. The five matter blocks retain every343-site and
882-edge coordinate. Their low-band ranks are12,8,24,12,4 for u^c,L,Q,d^c,e^c.

The active L/e^c calculation uses4,116 Nambu rows and24 low coordinates,
including the neutral L slot. Exact family tensor identities retain all
three families, giving12,348 Nambu rows and72 low coordinates. The twelve
actual family projectors and the disabled case are retained. All91 unordered
pairs are executed, covering169 ordered pairs by symmetry of distance.
Representative polar maps include a noncommuting family pair.

On the declared curved field, the scalar-band pairing leakage has Frobenius
norm0.012335169221822694. The complete BdG projector differs from the
number-conserving projector by about9.651951e-5. Site-dependent gauge changes
give a full BdG covariance error below6e-13. Polar transport of the occupation
grading has numerical error below4e-14. The proof of the projector ranks,
canonical occupation and polar domains uses the exact bounds above; these
floating residuals are execution checks and are not outward root certificates.

10. Scientific scope and reproduction

The construction supplies a finite, gauge-covariant classical-background
defect band with a canonical integer occupation and its original quartic
at trivial links. It retains the original pairing in the spectral problem
and provides complete-record alignment on the declared additional control.

It does not select a link background, gauge-field state, kinetic energy,
physical Higgs background, interacting continuum limit or particle masses.
The original source texture and charged coefficient architecture are inputs.
Neither a Wilson matrix nor a numerical spectral gap is identified with a
bare C/E repair operation by this theorem. The finite native quantum-link
case has a separate exact analysis in QUANTUM_LINK_PROOF.md.

Run from the repository root:
 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_gauge_band.gauge_band
 python -m model.particle_gauge_band.verify --check

The JSON receipt binds the literal source graph, actual kernel construction,
loaded source dependencies, exact rational bounds and binary witnesses.
The mathematical argument is stated here; no Lean formalization is asserted.
