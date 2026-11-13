# Native correlated baryon propagation

The complete infinite-chain baryon source and its compact-to-complement
intertwiner are those of `openwave_hadron_refinement_flow`. Its calibrated
rest energy does not determine the coefficient of spatial motion. With the
supplied strong length, the unchanged ground band rises by only about
33.1 keV across the Brillouin interval. Exact rational inertia and the full
Schur bound certify this mismatch with relativistic dispersion.

## Complete native construction

A normalized colour-epsilon triquark state of occupied span at most one
has a native creation polynomial C* with three CAR creators. Its adjoint
annihilates it to the vacuum. The sum runs over every spin and flavour
coordinate. For four light spin/flavour coordinates there are20 co-located
and80 edge shapes per translation cell. The exact multiset normalization
is1,3 or6; different multisets have disjoint occupation supports. Therefore
C_alpha C_beta* on the vacuum is delta_alpha,beta. The actual four-cell
native contraction checks all100 shapes, the complete rest matrix, both
spin/isospin Casimirs and its five-coordinate neutron reduction.

The entries used in exact bounds come from a separate occupation-coordinate
CAR enumeration in `exact_native.py`. For a multiset n, the unnormalized
epsilon column has Gram weight g=3!/product(n_i!). A one-body generator
E_ij has coordinate n_i+1 after replacing j by i, and its weighted numerator
is the exact integer identity g_new(n_i+1)=g_old n_j. The diagonal coordinate
is n_i. Spin and isospin squares use4J²=Z²+2Z+4J_minus J_plus. Applying these
identities to every declared source term constructs the complete rational
matrix directly. Its equality to the separately evaluated native contraction
is checked; floating rounding is only a consistency diagnostic and supplies
none of the exact matrix entries.

On the physical tree carrier, parallel transport every colour index to the
leftmost occupied cell before contracting epsilon. This gives gauge-invariant
creation polynomials. The parent zero-exterior-flux Gauss isometry gauges all
tree links to the identity, giving exactly the checked epsilon columns.
The same polynomials act on the full Fock carrier. An antiquark-conjugate
copy is included. Extra families obey the same rule: for d spin/flavour
coordinates the local count is binomial(d+2,3) and the edge count is
binomial(2d+2,3)−2binomial(d+2,3). No neutron projection or occupation-three
projector defines the operator.

Let M be the complete positive compact rest matrix, including all charges
and spin/isospin sectors. Positivity is checked by exact directed Gershgorin
bounds. Declare J=zeta M^−1, where

    zeta=(hbar c/a)^2/(2u^2).

Here u is the inherited electron energy unit and a is the independently
strong-calibrated length. The coefficient uses the native rest operator;
there is no separate fitted proton or neutron kinetic coefficient. The
physical speed in this assignment is a declared input, not a predicted c.
The additional Hamiltonian is

    sum_(x,alpha,beta) (C_(x+1,alpha)−C_(x,alpha))*
         J_(alpha,beta) (C_(x+1,beta)−C_(x,beta)).

It is positive on the full Fock domain because J is positive. Each term has
support in at most three cells. A graph coefficient equal to one when every
cell in that support has degree two and zero otherwise makes the complete
two-cell operator change zero. The infinite path has coefficient one.
Full-Fock positivity is preserved without omitting any populated sector.

In the N=3 fibre, each C maps its own compact shape to the vacuum and kills
all other relative shapes. Consequently the added action is
2(1−cos k)J on the full compact space and zero on Q and on the P–Q blocks.
At k=0 it is identically zero on the complete fibre, including Q. This
cancellation is not asserted for N>3: the annihilators can then leave
spectator particles, and nuclear energies require a separate calculation.
The inherited nuclear finite-box masses are not transferred to this action.

## Exact band and complement bounds

The neutron compact symbol is the native five-coordinate matrix. Scaling
its co-located coordinate by sqrt6 gives rational Gram diag(6,1,1,1,1).
At a rational unit-circle phase, its Hermitian matrix is Gaussian rational;
realification retains both quadratures and duplicates every eigenvalue.
Rational inertia bisection encloses all five bands. In these coordinates
J has numerator zeta G A(0)^−1 G. This is the restriction of the full100-channel
rule, verified through the native intertwiner.

The inherited Q floor C and squared coupling norm5/64 are unchanged. An
upper bound on the compact spectrum is T+4||J||, strictly below C in the
executed neutron/proton sectors. Every low energy is thus enclosed by
[e_j−(5/64)/(C−T−2(1−cos k)||J||),e_j]. Every higher state is retained above C.
The proton uses its own block of the same inverse operator; a constant
neutron-to-proton shift is asserted only at k=0.

The source coupling to each quark neighbour gives ||H'(0)||,||H''(0)||<=6t.
For a simple full ground state with certified gap g, perturbation theory
bounds its unmodified curvature in absolute value by6t+2(6t)^2/g. The
correlated term has first derivative zero and second derivative2J.
Let e bound the compact ground energy, s bound its normalized P-vector
sine error, and q bound full Q probability. The inherited finite theorem
gives these quantities. Hence the expectation of J lies between
(1−q)zeta/e_upper−2||J||s and zeta/e_lower+2||J||s. Combining these inequalities
encloses the full repaired curvature without treating a Ritz eigenvector
as exact. The physical inertial mass follows from
m_inertial=(hbar c/a)^2/(u E''(0)). Its certified interval is approximately
[0.91038,0.97069] GeV and contains the calibrated neutron rest mass.

The ground-band kinetic bandwidth is about0.145065 GeV. The finite-range
cosine law supplies controlled small-momentum motion; it is not exact
relativistic dispersion over the full Brillouin interval. A three-dimensional
Lorentz limit and physical scattering poles are separate requirements.

## Observer domain

A normalized fibre or wave packet is the preparation domain. Spectral
readout retains all five low outcomes and the entire high complement,
commutes with the autonomous Hamiltonian and has zero mean energy work.
The operator estimates hold with arbitrary idle quantum references. No
uniform distribution on the infinite carrier is used. All empirical inputs
and the selected architecture are retained in the receipt.
