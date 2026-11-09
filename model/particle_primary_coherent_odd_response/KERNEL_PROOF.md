# Four original spectral measures determine the summed spatial response

## 1. The original source kernel and its exact symmetries

Use member zero and the complete instruments of the [spatial response
proof](PROOF.md). Fix one of its actual G1 or H creation multiplets,
with component index j. Let E_i denote the actual source energies, Pi_i
their spectral projections, and rho the unchanged selected polynomial
density. For a Borel set B of source-energy transfers define

    M_xy^lambda(B)
      = sum_(i,k) 1_(E_k-E_i in B)
          sum_(j in lambda) Tr[rho Pi_i F_(x,j)* Pi_k F_(y,j) Pi_i].   (1)

This is the original unnormalized creation kernel. Every finite M(B) is
Hermitian positive semidefinite: for a spatial coefficient vector f,
f* M(B) f is a sum of positive squared amplitudes with the original
state weights. Neither the input state nor the clock is changed.

The [actual spin covariance](../particle_primary_selected_spin_response/PROOF.md)
gives U(g)F_(x,j)U(g)*=sum_l D_lambda(g)_(l,j)F_(g x,l).
The same U(g) commutes with H, rho and every energy projector, and
D_lambda(g) is unitary. Summing all internal components in (1) therefore
gives M_(gx,gy)(B)=M_xy(B). This summation is essential; it does not
reduce the full component-resolved or cross-irrep kernel to scalar rows.

The literal 48 source spin lifts have 24 distinct proper spatial cube
rotations. Acting on ordered vertex pairs, these give four orbits,
with sizes 8,24,24,8, exactly the Hamming-distance shells d=0,1,2,3.
Every ordered pair has a proper rotation that exchanges its endpoints.
Hermiticity then makes its scalar spectral measure real. Thus

    M_xy^lambda(B) = m_d^lambda(B),   d=Hamming(x,y).                   (2)

The diagonal m0 is positive. The off-diagonal m1,m2,m3 are real signed
measures; they are not individual probabilities. Their time Fourier
transforms need not be real. The code verifies the ordered-pair orbits
using the actual source rotations and reuses the proved internal
coefficient covariance.

Global native occupations do not remove these off-diagonal terms.
Every F_x has the same increment Delta, so F_x* Pi_k F_y has total
increment zero for any x,y. Local occupations are changed by the native
hopping and are not conserved charges of the selected state.

## 2. Exact four-measure reconstruction and apparatus normalization

Write m=(m0,m1,m2,m3)^T. The normalized Walsh vectors are
f_q(x)=(-1)^(q dot x)/sqrt(8). A radial matrix (2) has Walsh eigenmeasure
Lambda_s for s=|q|, where

    Lambda = K m,

        [1  3  3  1]
    K = [1  1 -1 -1],          K²=8I,
        [1 -1 -1  1]
        [1 -3  3 -1]

    m = K Lambda/8.                                                   (3)

Each Lambda_s is a positive measure, with multiplicities 1,3,3,1.
The existing source and clock normalization is retained:

    separately selected phase:  p_s^lambda = Lambda_s^lambda/6;
    one row of the full menu:    p_s^lambda = Lambda_s^lambda/48.      (4)

In particular the inverse from the four representative complete-menu
row measures is m=6Kp, and their complete sum is

    sum_s (1,3,3,1)_s p_s^lambda = m0^lambda/6.                        (5)

These factors do not change with the internal dimension of lambda: the
complete apparatus divides every G1 and H success operator by sqrt(6).
The code verifies (3) by the full symbolic eight-by-eight Walsh transform,
checks the inverse, and directly sums all ordered vertex pairs in the
independent test. Four representative row measures determine the summed
creation response of each multiplet, rather than requiring 64 independent
spatial entries. They do not determine its full internal-component kernel.

The same reduction holds for the implemented finite-clock event kernels:
the actual clock uses the same functions of H for every address and
component, so it preserves the covariance used in (1). Its finite event
measures differ from the ideal transfer measures by the declared clock
error; they are not silently identified with ideal spectral projections.

## 3. The calculation that determines contrast

For integer r>=0, stationarity gives the exact ideal transfer moment

    integral omega^r dM_xy^lambda(omega)
       = sum_j Tr[rho F_(x,j)* ad_H^r(F_(y,j))],
    ad_H(F)=[H,F].                                                    (6)

For even order 2r this is the rho-weighted Gram expression

    sum_j Tr[rho (ad_H^r F_(x,j))* ad_H^r F_(y,j)].                     (7)

Equations (6)-(7) reduce the search to representative original-state
matrix elements at each distance, or equivalently four representative
Walsh event measures. The positive finite density floor only transfers
positive-effect lower bounds. It supplies no lower bound on a signed
difference of spatial rows.

Within this component-summed domain, some m_d with d>0 is nonzero if
and only if the four unnormalized Walsh measures are not all identical.
This is a contrast test, not a dispersion theorem. For example,

    m0=mu, m1=mu/8, m2=m3=0

gives positive eigenmeasures (11/8,9/8,7/8,5/8)mu. Their intensities
differ, but every normalized energy shape is the same mu/mu(R).
The exact control uses a point measure at source transfer9, without
assigning it a primary-source eigenvalue or a physical mass.

Conversely, a diagonal positive kernel m_d=0 for d>0 populates all eight
Walsh rows equally. The code uses diagonal weights9/5 and18/5 for the
two multiplets; these exceed all stated spatial-row probability floors
while exhibiting no contrast. This is a linear-algebra countercontrol,
not a computation of the actual source kernel. Populated rows, covariance
and their bounds therefore do not by themselves prove a nonzero m_d.

## 4. Why prepared transport moments cannot replace (6)

The [prepared composite calculation](../particle_primary_composite_transport/PROOF.md)
proves an off-diagonal third matrix element of H on specified vectors.
Those vectors are not the selected polynomial density, and the filled-beta
anchor is not an eigenvector of the interacting H. Its prepared final-energy
moments therefore cannot be substituted for the original transfer moments
(6). The [fourth-moment resolvent](../particle_primary_odd_resolvent/PROOF.md)
retains this distinction. Its source-density-floor bridge is a positive
final-event bound, without a signed inter-site contrast claim.

An exact finite countercontrol shows that stationarity, faithfulness,
conserved occupation labels and nonzero prepared cubic transport alone
cannot justify that substitution. Take three CAR flavors at two sites,
and an independent two-state sign Z. Set

    H = sum_(a=1)^3 (c_(1,a)* c_(2,a)+c_(2,a)* c_(1,a)) tensor Z,
    F_x = c_(x,1)* c_(x,2)* c_(x,3)* tensor I,
    rho = (4I-H)²/2240.                                               (8)

Each flavor number is conserved. The three hopping summands have norm
one, so ||H||<=3 and 4I-H>=I: rho is a faithful stationary polynomial
state. The unitary S=(-1)^(N_2) tensor X commutes with H and rho, fixes
F_1 and negates F_2. It forces M_12(B)=0 for every spectral event B.
On the vacuum with sign Z=+1, however, the normalized three-particle
vectors have vanishing inter-site moments through order two and absolute
third moment six. The code verifies the full128-dimensional identities,
all three conserved numbers, the polynomial normalizer and the first
five original transfer moments. This comparison is not the primary source
and is not a claim that its additional sign symmetry occurs there.
Replacing rho by (4I-H)^(2k)/Tr(4I-H)^(2k), for any positive integer k,
preserves faithfulness and the exact cancellation symmetry. In this fixed
finite dimension, increasing k concentrates on the lowest eigenspace of H,
so the mean ground excess tends to zero. Even arbitrarily small cooling
error does not remove this counterexample to the proposed state substitution.

The source-specific unresolved quantity is a nonzero signed shell measure
or moment in (6), with its actual state and full complement controlled.
If it is nonzero, further energy-shape, isolated-branch and refinement
estimates are required before a relativistic dispersion or GeV rest mass
can be assigned.
