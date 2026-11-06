# FINITE QUANTUM LINKS: COORDINATE COMPLETION, GAUSS LAW AND THE OCCUPATION BOUNDARY

1. The actual link representation

The declared matter representation has five irreducible blocks

 Q:(3,2)_{q=1}, u^c:(bar3,1)_{q=-4}, d^c:(bar3,1)_{q=2},
 L:(1,2)_{q=-3}, e^c:(1,1)_{q=6}.

The subscripts are the primitive integer q=6Y. Every matrix used here is
obtained from the native source's exact spin/exterior current lift on its
fifteen nonvacuum even masks. This is the same native representation as
model/native_quantum_gauss. Its link on x->y is

 L_xy=C|0> direct-sum sum_r (V_r,x^* tensor V_r,y),
 dim L_xy=1+36+9+9+4+1=60.

The endpoint generators on block r are -T_r^T tensor I and I tensor T_r.
The source's matter-coupled link entry is the nonunitary raising operator

 U^r_ab=|r;bar a,b><0|.

Its combined one-particle operator on V_r tensor L_xy is

 A_r=sum_ab |a><b| tensor U^r_ab.

For any input label b, A_r(|b> tensor|0>) is a sum of d_r orthogonal
terms. These outputs are orthogonal for different b. Thus

 A_r^dagger A_r=d_r I_{V_r} tensor|0><0|,
 ||A_r||=sqrt(d_r).

The link matrix cannot be substituted for a unitary Wilson matrix while
retaining the same numerical hopping bound or canonical-coordinate algebra.

2. A compact-group coordinate attachment on the same60 states

There is a precise conditional completion of the raising rule. Assume a
compact-group rotor realization L^2(G), with normalized Haar measure and
the five source representations. Let

 phi^r_ab(g)=sqrt(d_r) D^r_ab(g)

be normalized matrix coefficients, using the endpoint convention V_r^*
tensor V_r. The constant function is the vacuum. Their selected span is
the same60-dimensional representation as L_xy. Let P_menu project onto
that span. Multiplication by a matrix coefficient is denoted M(D^r_ab).
Then define

 U_tilde^r_ab=sqrt(d_r) P_menu M(D^r_ab) P_menu.

The factor sqrt(d_r) is required to recover the literal native amplitude:
U_tilde^r_ab|0>=|r;bar a,b>. This definition adds a rotor and its coordinate
multiplication. Native local-current covariance alone does not select that
attachment.

For normalized Clebsch-Gordan tensors C^{t}_{i;(a,c)} embedding a target
irrep t in r tensor s, the coordinate product is

 U_tilde^r_ab|s;bar c,d>
 =sum_{t,i,j} sqrt(d_r d_s/d_t)
        conjugate(C^t_{i;(a,c)}) C^t_{j;(b,d)} |t;bar i,j>.

Only target representations in the retained menu survive the projection.
The tensors used here are real exterior tensors, so the conjugation is
numerically trivial but its representation convention is explicit.

3. Complete product classification and exact tensor normalization

Primitive-charge addition restricts the nonvacuum products whose targets
can lie in the menu to

 Q tensor Q -> d^c,
 Q tensor u^c -> L,
 u^c tensor Q -> L,
 u^c tensor e^c -> d^c,
 e^c tensor u^c -> d^c.

These are link-flux representation transitions in coordinate
multiplication. The transported matter block r is unchanged; s and t
label the link's initial and final flux. They are not particle-conversion
reactions or decay channels.

There is no missing vacuum target: none of the five retained primitive
charges is the negative of another retained charge. Every displayed
product has multiplicity one. The corresponding exterior operation is
the wedge of two two-slot states into a four-slot state.

Use the actual native occupation masks. For disjoint two-slot masks a,c,
the coefficient is the sign of the permutation that sorts their
concatenated occupied indices. It vanishes for overlapping masks. Each
target row has respectively4,3,3,1,1 nonzero coefficients for the five
products above. Dividing by the square root of this number gives a tensor
C with C C^dagger=I_target.

The producer verifies, over exact algebraic numbers, all twelve equations

 C(T_r tensor I+I tensor T_s)=T_t C.

It also checks the partial contraction

 sum_{a,i} conjugate(C_{i;(a,c)}) C_{i;(a,d)}
       =(d_t/d_s) delta_cd.

The resulting product factors sqrt(d_r d_s/d_t) are2sqrt(3),3,3,1,1.
Together with the row normalizations they determine every extra link entry.
No fit, continuous search or particle-mass value enters these coefficients.

4. Local Gauss identities and exact norms

The product tensors are endpoint intertwiners, so the completed operators
obey exactly the same covariance equations as the vacuum-raising entries:

 [L_x(T),U_tilde_ab]=-sum_i T_ai U_tilde_ib,
 [L_y(T),U_tilde_ab]= sum_j T_jb U_tilde_aj.

They cancel the fermion endpoint commutators in

 H_xy=J sum_ab c_xa^dagger U_tilde_ab c_yb + adjoint.

Therefore [G_v(T),H_xy]=0 on the full Fock/link space, for every source
current and endpoint. The executable checks1,416 complete matrix identities,
including every matrix element in every native matter block.

On the full rotor Hilbert space, the matrix-valued multiplication operator
D^r(g) is pointwise unitary, hence has norm1. Compression cannot increase
its norm. The displayed sqrt(d_r) normalization gives an upper bound
sqrt(d_r), while its action on vacuum inputs attains that bound exactly.

In this specific menu an additional algebraic identity holds. For the
completed combined operator A_tilde_r,

 (A_tilde_r^dagger A_tilde_r)^2=d_r A_tilde_r^dagger A_tilde_r.

Thus A_tilde_r/sqrt(d_r) is a partial isometry. Its exact ranks are

 block      native rank   completed rank   norm squared
 u^c              3             18               3
 L                2              2               2
 Q                6             30               6
 d^c              3              3               3
 e^c              1             10               1

The extra transitions change the finite source dynamics even though the
vacuum action and local Gauss identities agree. Those extra transitions
are forced by the declared coordinate-multiplication attachment. The
raising-only source and the coordinate completion are distinct dynamics.

Finite compression does not preserve the coordinate multiplication algebra.
For the completed e^c entry, the exact squared Hilbert-Schmidt norm of
[U_tilde,U_tilde^dagger] is20. It maps the vacuum and u^c link blocks to
e^c and d^c blocks and has no unitary inverse on the finite menu. This
noncommutativity matters for projected fermion operators.

5. What local covariance selects about link energies

The independent endpoint group is G_left times G_right. Each summand
V_r^* tensor V_r is irreducible for that product group. The summands are
inequivalent, as their primitive endpoint charges differ; the vacuum has
charge zero. The actual native matrices generate the full matrix algebra
on each V_r. The producer checks dimensions d_r^2 for those associative
algebras, which verifies the required source irreducibility.

Schur's lemma therefore gives the complete Hermitian commutant

 H_link=epsilon_0 P_0+sum_r epsilon_r P_r.

There are six real coefficients. A common energy offset removes one,
leaving five independent energy differences. Positivity or declaring the
vacuum to be the lowest link state imposes inequalities, not numerical
equalities. The source rule epsilon_r=epsilon for all five excited blocks
is an additional design choice. Neither the current representation nor
the coordinate completion determines these energy differences.

6. A positive one-particle spectral statement with quantum links

Put either the native raising operator or its coordinate completion on
every edge of the343-site compact graph. The combined one-particle hopping
block has norm at mostsqrt(6)/2, because6 is the largest matter dimension.
The complete hopping norm is therefore at most3sqrt(6).

The onsite texture acts as identity on the link Hilbert space. Its zero
eigenspace has dimension4d D_link, where d=45 and D_link=60^882 when every
edge has the native finite link register. Ordered Weyl inequalities give

 central band:       |E| <= 3sqrt(6),
 complement:        |E| >= 64-3sqrt(6),
 band gap:          at least64-6sqrt(6)>0.

The corresponding one-particle spectral projector commutes with every
local total matter-plus-link Gauss generator. Its rank is
180 times60^882. This is a finite one-particle result on a joint
matter/link space; it is not a180-mode canonical fermion factor.

This bound excludes the link-energy term. If an included link Hamiltonian
has spectral width W, subtract its spectral midpoint. The sufficient gap
bound becomes64-6sqrt(6)-W. For the sum of882 independent source link
penalties epsilon(I-|0><0|), W=882epsilon. The sufficient domain is

 0<=epsilon<(64-6sqrt(6))/882.

Failure of this inequality does not prove that no useful low-energy
sector exists. It means that this uniform global-band estimate does not
establish one. A selected link state, energy window or different dynamical
analysis would be required.

7. Why an operator-valued projector is not a fermion number

Let P be a projector on a one-particle matter space tensor a quantum link.
Its entries P_ij are link operators. The naive lift

 N(P)=sum_ij c_i^dagger P_ij c_j

is Hermitian and positive. If P commutes with the total one-particle Gauss
action, then N(P) commutes with the full Fock/link Gauss action. Its
spectrum need not be integral, because the P_ij need not commute.

An exact example uses two matter coordinates and a two-state link. Set

 |Phi>=(|0,0>+|1,1>)/sqrt(2),
 P=|Phi><Phi|,
 P_ij=|i><j|/2.

P has rank1 and is invariant under V tensor conjugate(V), for V in SU(2).
On the complete eight-dimensional matter-Fock/link space, N(P) has
eigenvalues

 0 with multiplicity5,
 1 with multiplicity1,
 1/2 with multiplicity2.

On the two-occupied matter sector, N(P)=I_link/2. Add a static reference
transforming in the fundamental SU(2) representation. The link/reference
singlet tensored with the filled two-mode matter state is a total Gauss
singlet, and its N(P) eigenvalue is1/2. The discrepancy is not removed by
restricting to a gauge-invariant state.

The projected fields b_i=sum_j P_ij c_j also fail the projected canonical
anticommutation identity. The producer checks this as an exact matrix
identity and records a nonzero squared error trace.

NATIVE_LINK_CAR_WITNESS.md supplies a second exact example using the
literal native e^c raising operator and a two-site spectral projector.
Its naive number has eigenvalues 1 plus or minus
(1-64/sqrt(4097))/2, including a noninteger Gauss-singlet eigenvalue.
That two-site example is a test of the operator algebra, not the compact
343-site texture.

Neither example asserts that its P is the spectral projector of the
compact343-site source. It disproves the general inference that a
gauge-invariant one-particle projector with quantum-link entries defines
an integer occupation or a canonical fermion factor. A construction for
the actual joint projector must establish that additional algebra.

For comparison, the producer checks an exact four-mode Bogoliubov
transformation with commuting rational coefficients. Its two selected
modes obey the canonical anticommutation relations and their number
satisfies N(N-I)(N-2I)=0. This is the algebra used by the classical-link
construction in GAUGE_BAND_PROOF.md.

8. Scientific scope and reproduction

The finite coordinate completion, all local Gauss identities, operator
norms and link-energy commutant are exact. The rotor realization is a
stated attachment. It supplies neither a selected link Hamiltonian nor
a physical gauge-field state. The one-particle band result does not
identify a canonical180-mode interacting field in the joint quantum-link
space. No quark mass, confinement scale or continuum limit follows from
these finite algebraic identities.

Run from the repository root:
 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_gauge_band.quantum_link_completion

The output includes all completed60-by60 coordinate matrices, exact
Clebsch-Gordan tensors, source hashes and the full small counterexample.
No floating eigenvalue test is used for the algebraic conclusions.
