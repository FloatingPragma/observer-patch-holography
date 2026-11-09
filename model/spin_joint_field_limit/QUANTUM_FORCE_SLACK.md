# Quantum correction to the mean-spin force

The cubic spin-one-half source has a strictly positive difference between
its transverse bond coefficient and the squared staggered magnetization.
This excludes replacing each longitudinal spin by its mean in the
unfiltered microscopic field equations. It does not exclude a sharp
infrared component with corrected coefficients.

## Source and state

Use exactly

    H_B=J sum_<xy>(S_x.S_y-1/4)-B sum_x epsilon_x S_x^1,
    p=n^-2, B=J n^-5, L_p=K_B/(Jp), n>=2.

The ground is the same unique infinite-volume positive-field ground used
in `spin_joint_field_limit/PROOF.md`, with its actual positive generator
K_B. The supplied joint trajectory is unchanged. Its internal axis-one
rotations, cubic symmetries, beta translations and real structure give
the common nearest-neighbor correlations

    <S_x^a S_y^b>=diag(c_parallel,c_transverse,c_transverse),
    <S_x^1>=epsilon_x m_B,
    e_B=-3(c_parallel+2c_transverse),
    kappa_B=-(c_parallel+c_transverse).

The spin-dot energy density is -J e_B; the shifted source energy density
is -J(e_B+3/4). The selection term contributes -B m_B. Finite even cubic
tori of side at least four satisfy the same marginal calculation. The
infinite-volume statement follows from the specified fixed-B ground
limit; no zero-field state is substituted at finite p.

## An exact two-site inequality

Rotate both sites so axis one is the diagonal axis, and choose x even.
The two-site density matrix has the form

    rho=[[a,0,0,0],[0,b,z,0],[0,z,c,0],[0,0,0,a]],
    a,b,c>=0, 2a+b+c=1, bc>=z^2, b-c=2m_B.

Internal rotations eliminate entries between different total-spin
sectors. The opposite one-site means give the two equal corner entries.
Reality makes z real. Direct spin matrices give

    c_parallel=a-1/4, c_transverse=z/2,
    d=e_B/3-1/4=-a-z,
    kappa_B-m_B^2=(bc-z^2)+a(2d+1/2)+d^2+d/2.          (1)

Consequently d>=delta>=0 implies

    kappa_B-m_B^2>=delta^2+delta/2.                    (2)

The scalar implication is formalized in
`Lean/SourceLaw/SpinMarginalForce.lean`; the source marginal and
variational attachment here are separate hypotheses of that implication.

The entangled-dimer trial in `spin_charge_dispersion/DIMER_BOUND.md` has
spin-dot energy -31J/40 per site and staggered mean sqrt(6)/5>12/25.
Its field is aligned with the actual selection field. The variational
ground inequality therefore gives

    e_B+(B/J)m_B>=31/40+(B/J)sqrt(6)/5.

Since m_B<=1/2,

    e_B>=31/40-B/(50J),
    d>=delta_n=1/120-1/(150n^5)>=13/1600, n>=2.        (3)

The torus has N/2 dimer bonds and 5N/2 inter-dimer bonds; it is the
literal three-dimensional source. Open-box boundary corrections disappear
in its ground energy density. The trial is used only variationally.
Equations (2)-(3) prove

    kappa_B-m_B^2>=delta_n^2+delta_n/2,
    liminf_(n->infinity)(kappa_B-m_B^2)>=61/14400.      (4)

The zero-field selected ground obeys the last bound as well, using the
unchanged ordered limit. No measured energy or magnetization is used.

## The force approximation that (4) excludes

The actual fields are

    Q_p(F)=p^2 sum_x epsilon_x F(px)S_x^3,
    P_p(G)=-p sum_x G(px)S_x^2.

For a real smooth compactly supported F, define
S_F=p^3 sum F(px)^2 and D_F=sum_<xy>(F(px)-F(py))^2.
There are six bonds at every site. The exact ground double commutator
and canonical cross moment give

    Im<Q_p(F)Omega,(iL_p Q_p(F)-12m_B P_p(F))Omega>
      =6(kappa_B-m_B^2)S_F-(kappa_B/2)p^3D_F
         +(B/J)(m_B/2)S_F.                            (5)

The identity sum_edges(F_x+F_y)^2=12 sum_x F_x^2-D_F
accounts for the gradient correction. It follows directly from the
literal source current equation in `WARD_DYNAMICS.md`. The replacement
S_x^1 by epsilon_x m_B in those currents predicts the coefficient
12m_B in (5). This is an approximation being tested, not an inherited
source theorem.

For G define the sampled positive graph Laplacian

    (-Delta_p G)_x=p^-2 sum_(y neighbor x)(G_x-G_y).

The other exact mismatch is

    Im<P_p(G)Omega,(iL_p P_p(G)+m_B Q_p(-Delta_p G))Omega>
      =(kappa_B-m_B^2)p D_G/2
         +B m_B S_G/(2Jp^2).                          (6)

In the joint limit, S_F tends to integral F^2,
p^3D_F tends to zero, and pD_G tends to integral |grad G|^2.
The source variance envelopes in `spin_joint_field_limit/PROOF.md`
bound the norms of Q_p(F)Omega and P_p(G)Omega uniformly. Cauchy-Schwarz
applied to (5)-(6), with (4), proves that these unfiltered source-vacuum
force residuals cannot converge in norm to zero for nonzero F and
nonconstant G. For a unit-L2 annular test r/2<=|k|<=r the respective
liminf squared norms are at least

    6r(61/14400)^2,       (r^3/16)(61/14400)^2.         (7)

Schwartz annular tests follow by the same summable approximations used
for the source field bounds. The compact-space proof (5)-(6) does not
require momentum projection.

If kappa_B tends to kappa along a selected subsequence, m_B tends to m,
and unfiltered source-vacuum linear force residuals with coefficients
u,b did vanish, the same first-moment test would instead require

    u=12kappa/m,          b=kappa/m.                   (8)

Equation (8) is only necessary for that strong microscopic convergence.
It proves neither that these corrected forces close nor that their
frequency product is a measured propagation speed. Projection to the
bounded-history representation or a fixed energy window can discard
part of an unbounded force vector. Thus (5)-(7) do not exclude a
renormalized infrared shell, a sharp subbranch, or a different limiting
coefficient when energy escapes the chosen representation.

## Exact finite controls and theorem scope

The implementation verifies (1) by rational arithmetic and independently
by full two-spin matrices. At J=1,B=2/3 the actual dimer source has the
rank-one ground marginal with a=0,b=9/10,c=1/10,z=-3/10.
Its m=2/5,kappa=2/5 saturate (2) at delta=3/10. This is a sharp local
matrix control, not a replacement for the cubic ground. Both force-score
identities are checked on this exact source with its degree one in
place of the cubic degree six, preserving the source field and clock.

Koma's `arXiv:2105.04970v2`, Eq. (2.11) and Theorem 2.1, bound an energy
expectation of a positive-energy-filtered wavepacket between two linear
functions of momentum. They supply no variance saturation or pole. Their
volume-then-B limit also does not identify the joint family with the
zero-field history representation. The proof here consumes the actual
finite-B ground symmetries, dimer trial and direct source identities;
it supplies no additional sharp-dispersion premise.
