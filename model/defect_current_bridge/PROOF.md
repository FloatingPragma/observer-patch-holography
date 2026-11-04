# A common defect, current and mass operator

Let h be the retained radial lattice Hamiltonian on ell2(Lambda) tensor C2,
with lattice spacing and hopping velocity one and texture coefficient64.
Lambda is an odd open cube centered at the origin. These numbers use the
chosen common energy unit v(P); identifying that unit with a laboratory clock
is a separate physical premise. The radial kernel theorem supplies an
isometry W:C2 -> ker(h), a complete two-dimensional kernel and a complementary
gap d>=62. Its proof applies to the full lattice operator, including every
lattice valley. The numerical witness uses27 sites and54 rows.

Write J=i sigma2 and D_x=(-1)^(x1+x2+x3), C=D J. The Pauli identities give
C^T=-C, C^dagger C=I and h C=-C conjugate(h). At the origin fix the quaternion
frame of the kernel. Inversion evenness and the source antiunitary symmetry
then give C conjugate(W)=W J. This identity follows throughout the kernel
from its origin evaluation, which is an isomorphism. The normalization Gram
is scalar, so normalization preserves the identity.

## The chosen carrier

Choose the physical CAR algebra on

    V = ell2(Lambda) tensor C2 tensor M15 tensor F3.

M15 is the actual retained native exterior matter representation, transported
by its source intertwiner into the action basis. F3 comprises three declared
gauge-trivial family copies. They are design inputs, not consequences of the
two-dimensional radial index. Let K be any symmetric45-row all-left Higgs mass
kernel. Define the local creation pairing Delta=C tensor conjugate(K) and
the quadratic Hamiltonian with one-body term h tensor I45. Delta is
antisymmetric, as required by fermion creation. The corresponding Nambu
operator is

    B = [[h tensor I, Delta], [Delta^dagger, -conjugate(h) tensor I]].

Nambu rows duplicate annihilation and creation coordinates; no additional
physical mirror fermion is inserted. The complete source has2430 physical
annihilator modes and4860 Nambu rows in the displayed cube.

The isometry E=diag(W tensor I45, conjugate(W) tensor I45) satisfies

    B E = E [[0, J tensor conjugate(K)],
             [(J tensor conjugate(K))^dagger, 0]].

This follows by direct multiplication using hW=0 and C conjugate(W)=WJ.
Thus the same localized90 annihilator modes carry the180-row core operator.
Its eigenvalues are plus and minus every singular value of K, each twice,
including all zero singular values. These are exact identities for the
mathematical source coefficients. Floating-point matrix checks in the
receipt are diagnostic controls, not interval proofs of exact zero.

## Native currents and states

For every actual anti-Hermitian native generator G on M15 tensor F3,
(I tensor G)(W tensor I)=(W tensor I)(I2 tensor G). On Nambu coordinates
the generator is diag(G, conjugate(G)); a Hermitian charge instead uses
diag(Q,-conjugate(Q)). The source M15 intertwiner proves this identity for
the original twelve native port currents, without assigning their names to
unrelated operators. Gauge covariance acts on the Higgs background as well
as on matter. This intertwiner is for the spatially uniform internal
action; no dynamical gauge links or local Gauss-law constraint are supplied.
A fixed neutral Higgs does not make the full broken gauge
algebra commute with its mass term.

The induced CAR inclusion c(f) -> c((W tensor I)f) is injective because W is
an isometry. Quadratic Heisenberg evolution restricts to the displayed core
operator. The reducing orthogonal complement gives a graded CAR product.
For parity-even density matrices, rho -> rho tensor I_complement/dim is an
exact state inclusion with partial-trace recovery, intertwining evolution
and the gauge action on the family of Higgs backgrounds. The normalized
trace on the complement is not a vacuum. No full many-body ground selection
is inferred from this construction.

Every normalized core spinor has the same site probability rho_x because
W_x^dagger W_x=rho_x I2, sum rho_x=1. Inversion gives mean position zero.
Translation of the chosen texture translates this density. Fourier
transformation supplies a momentum distribution, not a momentum eigenstate
or a relativistic energy-momentum relation. Occupying the core does not
change the prescribed texture winding.

## Pixel-dependent coefficients in this operator

Use the source memory construction's real positive matrix M(P), with
eigenvalues the three retained charged-lepton ratios m_i(P)/v(P).
Choose H=(0,v(P)/sqrt(2)) and Ye=sqrt(2) M(P). The charged-lepton term in
K/v(P) is the actual native neutral-Higgs tensor times M(P), hence has
exactly those singular values with the source charge-conjugate multiplicity.
The memory construction supplies M(P) from its unconditioned C-duty
stationary-mean coefficient and declared source gain. A chosen coherent
pairing service implements that coefficient in the Hamiltonian defined
here. This equality does not identify a finite unread memory channel with
a unitary Hamiltonian. The full native repair instrument is not proved to
generate or preserve this coherent service. The two spatial spin rows do
not establish relativistic spin or chirality.

For a complete displayed45-row witness, choose the retained Stage5 quark
design epsilon=1/6 with Yukawa powers (6,3,0) and (6,4,2), as defined in
`code/P_derivation/paper_math.py`. Set the three minimal neutrino couplings
to zero. These choices give rank42 and nullity3 for K, and84 positive,
84 negative and12 zero Nambu core rows. The zero rows are retained; no
massless relativistic neutrino interpretation follows. The quark powers
are declared architecture, not selected by topology. The retained
charged-lepton architecture has calibration provenance; absence of measured
lepton values in this calculation does not establish independence from
that provenance.

## Complete record implementation

There is a controlled implementation of this same Hamiltonian by a chosen
record service. The native current probes act on fresh auxiliary left and
right legs in the invariant input preparations of the memory construction.
They do not measure the spatial fermion field. Thus their transition
probabilities and the twelve-state memory chain are independent of that
field, without deleting any outcome. If j is the updated memory port, let
L(P_j) be the local lepton pairing Hamiltonian with family coefficient
P_j. For every record use the finite Fock Hamiltonian

    H_record = H0 + 1_C * gain * L(P_j),

where H0 contains the full radial kinetic operator and the declared quark
pairing. Linearity in the family matrix makes its stationary mean exactly
H0+L(gain Z_C), the Hamiltonian defined above. Every record Hamiltonian
individually has the same core intertwiner. Its unitary therefore preserves
the embedded CAR algebra and the product state inclusion before averaging.
This is an explicit coherent service driven by the actual complete record
law; it is not an identification with the original site-field repair.

For completeness, on the fixed2430-mode Fock space let M bound every
record Hamiltonian norm. The one-body radial norm is at most66. The norm
bound follows because the kinetic term is sum_i sigma_i D_i, where the
three centered-difference operators commute on the Cartesian open cube
and each has norm at most1. The square of the kinetic term is sum_i D_i^2,
so its norm is at most sqrt(3)<2; the texture multiplication has norm64.
The norm
of dGamma(h) is at most2430*66. A pairing kernel Delta on n modes has
Hamiltonian norm at most n||Delta||/2, by its skew singular-block form.
Both the quark kernel and the unit-gain lepton kernel have norm at most1.
Consequently M=2430*66+1215*(1+gain) is sufficient. This bound does not
materialize the exponentially large Fock matrix.

The complete memory transition is T=(1-p)R+pS, where R resets every input
to port2. On zero-sum operator-valued memory blocks its trace norm
contracts by at most p. Each record unitary of duration h differs from
identity as a trace-norm channel by at most2Mh. The stationary mean step
has generator -i[H0+L(gain Z_C),.] and one-step remainder at most4M^2 h^2.
The resulting deviation of the memory blocks from their stationary
product is at most2p^k+4Mh/(1-p). Summing the step errors at h=T_phys/N
gives, for arbitrary states and finite untouched references,

    error <= 4 M T_phys/[N(1-p)]
             +(8/(1-p)+4) M^2 T_phys^2/N.

For a stationary product memory preparation the first term vanishes.
The service converges to the stated full Fock Hamiltonian on every fixed
finite volume. The finite-step unread service can be nonunitary. The bound
grows with volume and supplies no interacting continuum theorem.

## Complete spectrum and domain

With T=diag(I,C^dagger tensor I), direct multiplication gives

    T^dagger B T = [[h tensor I, I tensor conjugate(K)],
                   [I tensor K, h tensor I]].

The two tensor factors commute, so the complete spectrum consists of
lambda(h) plus or minus the singular values of K. It has additive bands;
no square-root relativistic dispersion is asserted. If ||K||<=M and2M<d,
the core lies in[-M,M] and the complement outside(-d+M,d-M). This is
single-quasiparticle band isolation, not isolation of every sum of core
many-body energies. The native neutral-Higgs tensors split K into disjoint
charged Dirac blocks, so its norm is the maximum of the family norms.
The source rational certificate gives sum_i m_i/v<1, and both declared
quark norms are below1. Thus for the displayed design M=1 is sufficient in v(P)
units, so the complement has absolute energy at least61 and separation60.

Interactions can mix the two factors; their exact spectrum is not given by
this quadratic argument. The original site-resolved repair generally fails
to preserve W. A direct origin-site projection leaks out of the kernel,
and removing the stagger sign from pairing also breaks the intertwiner.
These countercontrols distinguish exact core preservation from accidental
agreement of selected eigenvalues. Dynamical texture selection, physical
pole masses, a common interacting continuum and responsive Einstein
geometry are outside the asserted theorem.
