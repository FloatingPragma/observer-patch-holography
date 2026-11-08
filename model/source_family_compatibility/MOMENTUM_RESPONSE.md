# Magnetic momentum and the actual composite response

The centered full-CAR source has a three-correlator reduction of its first
composite response on an invariant full-ground average. Its momentum is the
joint momentum of the actual two-step magnetic translations. The reduction
retains the omitted response Gram at every momentum; it does not set the
memory to zero or identify an eigenvalue of a compressed matrix as a particle.

## Source, state and translation convention

Use the signed chart of `COMPOSITE_CAR.md` and
`ALIGNED_HOLE_DISPERSION.md`, on all occupation sectors:

    K=T+U sum_x[n_xup n_xdown-N_x/2],
    h_(x,x+ej)=(c/2) eta_j(x),
    eta_j(x)=(-1)^(sum_(l>j) x_l),     c>0, U>0.          (1)

The original Weyl chart is related by its specified onsite unitary.
Transporting an experiment transports its probes and instruments too.
The independent original-Weyl local-history comparison does not require
invariance of bare repair under that transport.

Start with the full centered-ground projector, tracial exterior and
translation averaging of `centered_source_ground`. The normalized full
projector is invariant under hidden spin, staggered eta spin, particle-hole
conjugation and complex conjugation in this real chart. These properties
survive the exterior extension and the weak-star ground construction.
On one-body vectors the actual magnetic translations are

    (U_j psi)(x)=(-1)^(sum_(l<j) x_l) psi(x-ej).

They preserve (1), obey U_i U_j=-U_j U_i, and have commuting squares.
An original-Weyl ordinary translation transports to U_j times the constant
hidden-spin rotation i sigma_j. Hidden-spin invariance therefore turns the
specified ordinary-translation average into a magnetic-translation invariant
state in the signed chart. This reuses the same ground construction and
clock; it is not translation of an open finite cube as though it had no edge.

For an axis permutation pi and signs s_j, put

    (R x)_(pi(j))=s_j x_j,
    g_R(x)=(-1)^(sum_(i<j,pi(i)>pi(j)) x_i x_j).          (2)

If det R=+1, the 24 resulting signed permutations are the proper cubic
rotations. Multiplying a field at x by g_R(x) and moving it to R x is an
exact symmetry of the real hopping and onsite interaction. Equation (2)
follows by comparing the literal link signs on each transformed bond.
The same formula works when a bond reverses orientation because eta_j
is independent of x_j. These lifts normalize magnetic translations up to
global gauge transformations.

Use the arithmetic average over these 24 automorphisms of the specified
translation-invariant ground. Ground inequalities are linear in the state,
and each automorphism preserves K, so the result is a ground state on the
same action. It retains hidden spin, eta spin, particle-hole symmetry,
reality and magnetic translations. All expectations below use this one
averaged state omega. The finite average keeps its original probabilities;
there is no momentum-dependent conditioning or change of state. A
symmetry-broken pure ground is not silently identified with this average.
The compactness construction supplies no effective volume convergence rate.

The physical-cell Bloch convention is

    (F psi)_s(k)=sum_n exp[-ik.(2n+s)] psi(2n+s),
    s in {0,1}^3, k in [-pi/2,pi/2)^3, measure dk/pi^3.

Thus U_j^2 multiplies by exp(-2ik_j), and the actual source symbol is

    h(k)=c sum_j cos(k_j) Gamma_j,
    h(k)^2=c^2 S(k) I8,   S(k)=sum_j cos(k_j)^2.          (3)

The eight parity components and two hidden-spin components are retained.
They are not eight particle species. Open-cube Walsh preparations do not
have the translation-eigenvector interpretation of this infinite-lattice
Fourier decomposition.

## Three explicit local source correlations

Let b_xs=(2n_x,opposite-1)a_xs. The source derivative is the exact operator

    Y_xs=[b_xs,T]=sum_(y neighbour x) h_xy D_xy,s,
    D_xy,s=(2n_x,opposite-1)a_ys
       +2(a_x,opposite* a_y,opposite
           -a_y,opposite* a_x,opposite)a_xs.             (4)

No interaction term occurs in Y because the full equation is
Lret V1=(U/2)V0+Y. The retarded Gram is
G_(x,s),(y,t)=omega({Y_xs,Y_yt*}); its original addition/removal components
have their distinct physical positive-energy clocks as specified in
`COMPOSITE_CAR.md`.

Write Q_xy=S_x.S_y+eta_x.eta_y with the actual staggered eta signs. At the
origin all six incident real links have sign +1. Define the normalized
neighbour orbital r_s=6^(-1/2)sum_(|y|_1=1) a_ys. It is a canonical
orbital of the original source and anticommutes with a_0s. It adds no
mode. Its eta raising operator has sign -1, while the origin's has sign
+1. Let Q_0r be the same spin-plus-eta dot product on these two orbitals.
Define the actual ground correlations

    q1=omega(Q_0,e1),  q2=omega(Q_0,2e1),
    qstar=omega(Q_0r),
    kappa=4q1,  g0=c^2(9/2-12qstar),  g2=c^2 q2.        (5)

These are respectively two-site, two-site and seven-site observables of
one source state. Cubic and magnetic symmetry equate their choices of
axis and origin. In particular the second compression is exactly

    C_xy=kappa h_xy,                 C(k)=kappa h(k).     (6)

Equation (6) consumes the bond identity proved in `COMPOSITE_CAR.md`.
It requires the stated invariant average. Gauge or dimension counting
alone would not equate different bond expectations.

There are direct operator checks for both remaining coefficients. In the
following identities the sum averages the two actual hidden-spin labels:

    (1/2)sum_s {D_0r,s,D_0r,s*}=3I-8 Q_0r,
    (1/2)sum_s Re{D_0,e1,s,D_2e1,e1,s*}=4 Q_0,2e1.      (7)

Here D_0r is (4) on the central and normalized neighbour orbitals, and
Re A=(A+A*)/2. Consequently G_00=(6c^2/4)(3-8qstar), and the axial
coefficient G_0,2e1=(c^2/4)4q2. The middle site in the second identity
cancels from the anticommutator. `momentum_response.py` verifies (4) and
(7) over exact Gaussian rationals in the full Majorana algebra, before
any state is supplied. This avoids a contraction or factorization of
source four-point functions.

Each Q is the sum of two dot products on mutually exclusive spin and
eta blocks, so its spectrum is contained in {-3/4,0,1/4}. Therefore

    -3 <= kappa <= 1,
    3c^2/2 <= g0 <= 27c^2/2,
    -3c^2/4 <= g2 <= c^2/4.                              (8)

The producer's generic norm error bound permits the larger symmetric
interval |kappa|<=3. No sign of q2 is inferred from symmetry or from a
numerical ground approximation.

## Exact range-two symbol and positivity at every momentum

Two odd operators with disjoint site supports anticommute on the full
CAR algebra. Each term of Y_x is supported on x and one neighbour;
hence G_xy=0 for graph distance greater than two. Hidden-spin symmetry
makes G_xy proportional to delta_st.

Particle-hole conjugation maps a_xs to epsilon_x a_xs*, b_xs to
-epsilon_x b_xs*, and Y_xs to epsilon_x Y_xs*. Thus
G_xy=epsilon_x epsilon_y G_yx. Reality and Hermiticity give G_xy=G_yx,
so every odd-separation coefficient vanishes.

For a mixed-axis displacement delta=+/-ei+/-ej, choose the proper rotation
that exchanges axes i,j with signs chosen to fix delta and reverses the
third axis. It fixes both 0 and delta, while (2) has

    g_R(0)=1,                g_R(delta)=-1.

Invariance therefore gives G_0,delta=-G_0,delta=0. The executable checks
all twelve such displacements with their literal source phases. The
remaining axial coefficients have no magnetic phase under a translation,
because their endpoints have identical parity. Cubic rotations and
Hermiticity equate all six to the real number g2. It follows that

    G(k)=[g0+2g2 sum_j cos(2k_j)] I,
    R*R(k)=[g0-6g2+(4g2-kappa^2 c^2)S(k)] I.            (9)

The identities act on all eight parity and two spin components. They are
Fourier decompositions of the actual source anticommutator kernels, not
an assumed scalar self-energy. In particular C^2=kappa^2 h^2 and
R*R=G-C^2 reuse the exact orthogonal response embeddings.

Since S(k) ranges over [0,3] and the symbols are continuous, positivity
of the full source Gram is equivalent to the two endpoint inequalities

    g0-6g2 >=0,
    g0+6g2-3kappa^2 c^2 >=0.                             (10)

They are source constraints on the three local ground correlations.
Similarly G>=0 gives g0>=6|g2|. Equations (8)-(9) imply the useful
uniform source bound

    0 <= G(k) <= 18 c^2 I.                              (11)

This also bounds ||Y f||^2 by 18c^2||f||^2 for every original normalized
coefficient packet. It does not determine the distribution in energy or
the memory R*(z-Lperp)^(-1)R. A nonzero residual is compatible with
sharp subbranches as well as continuous response.

For reference the first four original signed moments are explicitly

    M1=h, M2=h^2+alpha^2 I,
    M3=h^3+alpha^2(2+kappa)h,
    M4=h^4+alpha^2(3+2kappa)h^2+alpha^4 I+alpha^2 G,
    alpha=U/2.                                         (12)

These are measured response moments. Diagonalizing the two-block
compression without its memory does not turn them into a complete
spectrum. Equation (9) is a scalar Gram; no scalar frequency-dependent
memory is inferred from it.

## Complete local readouts and controlled input errors

The producer expands (5) in Hermitian Majorana words
P_I=i^(|I|(|I|-1)/2) product_(j in I) gamma_j, each with P_I^2=I.
They are even finite Pauli observables in the fixed occupation chart.
The exact score plans contain 6, 6 and 216 terms for q1, q2 and qstar,
with coefficient absolute sums 3/4, 3/4 and 9/2 respectively.
For O=sum_l a_l P_l, select l with probability |a_l|/M, where
M=sum_l|a_l|, measure both signs of P_l and report M sign(a_l)s.
The full effects sum to I on every input, the mean is omega(O), and
no source outcome is divided by a success probability. A zero observable
uses an identity record and score zero. Native coherent compilation uses
the specified complete local services, including their failure records;
the formula does not derive those services from unread repair.

A source/reference output half-trace error eta changes this mean by at
most 2M eta. A local state half-trace enclosure gives the same bound.
If admitted errors in kappa,g0,g2 are delta_k,delta_0,delta_2, with both
bond estimates clipped to [-3,3], the whole Brillouin zone satisfies

    ||delta G|| <= delta_0+6delta_2,
    ||delta(R*R)|| <= delta_0+6delta_2+18c^2 delta_k.      (13)

Thus finitely many actual local readouts control the complete spatial
symbol at this response level. A quantitative growing-volume ground
marginal estimate is needed to turn finite estimates into certified
values for a selected infinite-volume state. The symmetry average,
source clock and complete probabilities are fixed throughout.
