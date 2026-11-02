# Actual lattice entropy on a common shrinking family

FA101 extends FA87/91/94/96 on precisely the U=0 free source.
It does not supply the area, edge, Newton, or stationarity inputs of the
existing EinsteinRegisterData structure. All eight lattice valleys remain
in every finite regional covariance; a chosen shell prepares one continuum
summand. Units below use x0=v t, so the Cauchy Weyl Hamiltonian is alpha.p.

## 1. Natural shell vectors remove the logarithm regulator

Let P be any orthogonal projection, E a regional projection, A=EPE on EH,
and let g be exterior: Eg=0. Set F=EPg. The off-diagonal block B=EP(1-E)
satisfies BB*=A(1-A). Polar decomposition therefore gives

    F=[A(1-A)]^(1/2) h,             ||h||<=||g||.

This is an operator factorization, not an estimate inferred from a trace.
For k(x)=log((1-x)/x) and
k_delta(x)=log((1-x+delta)/(x+delta)), 0<delta<=1/4,

    sqrt(x(1-x)) |k(x)-k_delta(x)| <= 2 sqrt(delta),
    sqrt(x(1-x)) |k(x)| <= 2.                              (1)

Indeed log(1+y)<=sqrt(y), because 1/(1+y)<=1/(2sqrt(y)),
and integrating starts at zero. Apply this separately to delta/x and
delta/(1-x). The second bound follows from -log x<=1/sqrt(x) and
-log(1-x)<=1/sqrt(1-x), or the same elementary integral bound.
Consequently F is in Dom k(A), and

    ||[k(A)-k_delta(A)]F||<=2sqrt(delta)||g||,
    ||k(A)F||<=2||g||.                                    (2)

The finite lattice A is strictly between zero and one by FA96. The same
paper proves that the continuum regional C has no endpoint eigenvectors.
No uniform gap is required in (1)-(2). The continuous endpoint value of
sqrt(x(1-x))k(x) also makes these assertions valid on the relevant range
if unused endpoint eigenspaces were present.

For unit exterior g put p=<g,Pg> and B_g=c*(g)c(g)-p I. On the regional
algebra the local states omega_theta(X)=omega((1+theta B_g)X), |theta|<1,
are positive and normalized. Wick's determinant formula and its rank-one
update give their exact covariance A_theta=A-theta FF*. They are
gauge-invariant quasifree states: every determinant contains at most one
column from the rank-one update. Equivalently, when p>0 and theta>=0 this
is the mixture (1-theta p)omega+theta p omega_hole; the normalized hole
covariance is A-FF*/p. The p in this conditional covariance is NOT the
normalization of a cell-averaged shell. For negative theta positivity is
proved by the local positive multiplier, not a global affine vacuum family.

For raw lattice g_a=M_b J_a g, put r_a=||g_a||^2. The actual normalized
number mode is g_a/sqrt(r_a), and its local covariance is therefore

    A_a(theta)=A_a-(theta/r_a) F_a F_a*,  F_a=E_a P_a g_a.  (3)

The finite many-body density is faithful for |theta|<1, since the local
state dominates (1-|theta|) times the faithful vacuum density. Its ordinary
von Neumann entropy derivative is exactly

    s_a'(0)=-<F_a,k(A_a)F_a>/r_a.                         (4)

No absolute entropy is assigned to a continuum type-III algebra.

## 2. Explicit estimates on the unit ball

Write eta=a/R, let E_eta select sites eta n in the unit ball, and let J_eta
average over centered cubes. Work in one fixed valley frame; conjugating
P_eta by M_b gives the corresponding multiplier with signed sin(eta p)/eta.
Write Pi for the corresponding continuum negative-band projector.
For 0<s<=1 and 0<eta<=1,

    ||(P_eta J_eta-J_eta Pi)f||<=3 eta^s ||f||H^s,
    ||(J_eta*J_eta-I)f||<=3 eta^s ||f||H^s.                (5)

Here P_eta in (5) denotes the conjugated full lattice projector, not a
one-valley truncation. To prove the first bound, split the input at
|p|=1/eta. On the low part the cell Fourier formula has no aliases;
|sin(eta p)/eta-p|<=eta^2|p|^3/6, and the difference of the rank-one band
projectors is at most eta^2|p|^2/6. The cell sinc factor is a contraction.
The high part costs at most 2 eta^s ||f||H^s. Thus 13/6<3 suffices.
The second estimate follows by the same split from the H1 cube-average
bound sqrt(3) eta ||grad f|| and the norm-one projection bound on the high
part. The decomposition is a proof device, not a bandlimit assumption.

We use the deliberately coarse multiplier estimate

    ||chi_B f||H^(s/12) <= (2048/s)||f||H^s.               (6)

For completeness, set p=2/(1-s/3). Unitary Fourier Hausdorff-Young and
weighted Holder, using integral(1+|xi|^2)^(-3) dxi=pi^2/4, give
||f||Lp<=2||f||H^s for two-component spinors. A ball boundary layer of width
h<=1 has volume at most 64h. Holder bounds its L2 mass by
4 h^(s/6)||f||H^s. Thus translations of chi_B f obey
||tau_h(chi_B f)-chi_B f||<=6 |h|^(s/6)||f||H^s for |h|<=1;
for larger h use 2||f||H^s. In the difference-quotient identity with
t=s/12, the integral is at most (960pi/s)||f||H^s^2<3072||f||H^s^2/s.
The Fourier coefficient c_t=integral |exp(ih1)-1|^2/|h|^(3+2t) dh exceeds
1/1024: restrict to the radius-1/2 ball centered at pi e1, use numerator>=3,
|h|<4, and 3+2t<5. Hence (6) follows from
1+3145728/s <= (2048/s)^2. This also proves all boundary-layer estimates
used below rather than postulating smoothness after a sharp restriction.

Put T_eta=E_eta M_b J_eta, C=chi_B Pi chi_B. For f supported in B,
the intertwining defect obeys

    R_eta(f)=A_eta T_eta f-T_eta Cf,
    ||R_eta(f)||<=12 eta^(s/6)||f||H^s.                  (7)

The two boundary-strip terms are <=4 eta^(s/6)||f||H^s each and (5) supplies
the third term. Centered cubes touching the sphere lie within distance
eta, so this estimate applies to the actual site-center convention.
Also

    ||(T_eta*T_eta-I)f||<=7 eta^(s/6)||f||H^s.             (8)

For unit exterior g with ||g||H1<=B, B>=1, write F=chi_B Pi g and
F_eta=E_eta P_eta M_b J_eta g. Once g_eta is exterior,

    ||F_eta-T_eta F||<=8B eta^(1/6),
    1-r_eta<=3eta^2 B^2.                                (9)

If eta<=1/(6B), then r_eta>=11/12>1/2. For a general exterior g, the separate
condition that its averaged cells remain exterior must also be imposed;
an arbitrary exterior support can approach the boundary. In the particular
shell applications below, eta<=1/(6B) ensures disjoint lattice supports.

Iterating Y=1-2C loses controlled fractional regularity, rather than
assuming any power C preserves H1. Set

    s_j=12^(-j-1),
    K_j=2048 * 8192^j * 12^(j(j+1)/2).

Then ||Y^j F||H^s_j<=K_j B. The initial bound is (6), and each step costs
at most 1+4096/s_j<=8192/s_j. In particular for j<=m,

    ||R_eta(Y^j F)||<=12 K_m B eta^[1/(6*12^(m+1))].       (10)

## 3. An effective natural-vector entropy modulus

Let r=(1+2delta)^(-1), and truncate
k_delta(C)=2 sum_(j>=0) r^(2j+1) Y^(2j+1)/(2j+1) after N terms.
Its norm tail is bounded by

    tau_N=2 r^(2N+1)/[(2N+1)(1-r^2)].

Power telescoping with Y_a T-TY=-2R gives

    ||k_delta(A_a)T F-T k_delta(C)F||
      <=2tau_N||F||+[4r/(1-r^2)] max_(j<2N-1)||R(Y^jF)||.

For delta<=1/4, 4r/(1-r^2)<=2/delta and ||k_delta||<=1/delta.
Combine (2), (8)-(10), both natural-vector replacements, and r_a>=1/2.
With m=2N, a_m=1/(6*12^(m+1)), one obtains

    |<F_a,k(A_a)F_a>/r_a - <F,k(C)F>|
    <= 8sqrt(delta)+4tau_N
         +(256 B^2 K_(2N)/delta) eta^a_(2N).              (11)

For explicit accounting before normalization the scalar bound is
4sqrt(delta)+2tau_N+(64 B K_m/delta)eta^a_m. Normalization doubles it and
adds at most 12eta^2 B^2; these fit within (11). This proves convergence
of the actual unmodified shell vector, not only of a specially matched
graph-core approximation.

Given 0<epsilon<=1, define

    delta=(epsilon/32)^2, N=ceil(16/(epsilon delta)),
    eta <= min(1/(6B),
       [epsilon delta/(1024 B^2 K_(2N))]^[6*12^(2N+1)]).   (12)

Because 1-r^2>=delta, tau_N<=1/(N delta). Each term of (11) is <=epsilon/4,
so the total is <epsilon. The astronomical powers in (12) are an explicit
mathematical modulus, not a feasible simulator resolution. The certificate
retains them symbolically and does not underflow them to numerical zero.

## 4. Continuum entropy/stress identification on the natural vector

FA87 uses the same Weyl vacuum and fixes the modular sign and 2pi factor.
The imported La Piana--Morsella Theorems 3.4/4.1, interpreted with their
wave complex structure, and FA96's CAR modular uniqueness identify k(C)
with the action of q_R=pi{w_R,h} on smooth spinors compactly supported inside
the ball. This uses the selfadjoint modular generator furnished by the
imported theorem; no new assertion of essential selfadjointness is needed.
For the exterior shells, F=chi_B Pi g is smooth up to the boundary by the
separated kernel, even though zero extension need not lie in H1. Multiply
F by smooth interior cutoffs removing a boundary strip of width h. These
vectors belong to Dom k(C) by the just-stated interior action. Their q_R
images converge: w_R=O(h), the cutoff derivative is O(1/h), and the strip
has volume O(h). The remaining smooth terms converge in L2 as well.
Closedness of k(C) therefore gives k(C)F=q_R F, with q_R evaluated on the
interior smooth spinor. Thus

    -<F,k(C)F> = 2pi integral_B w_R(x) delta T00(x) dx.     (13)

The minus sign is the hole covariance tangent. The right side is exactly
the normal-ordered stress response already derived in FA91, not a newly
chosen classical stress or a finite entropy assumed to obey a first law.

## 5. Actual finite entropy differences and the shrinking limit

The positive endpoint states omega_±=omega((1±B_g) .) have mean omega.
Classical-bit to quantum-state data processing gives

    0<=D(rho_theta||rho_0)
      <=[(1+theta)log(1+theta)+(1-theta)log(1-theta)]/2
      <=theta^2,                        |theta|<=1/2.     (14)

At finite lattice spacing the Gibbs identity is exact:
S(rho_theta)-S(rho_0)=theta s_a'(0)-D(rho_theta||rho_0).
It follows that the central entropy difference differs from s_a'(0)
by at most |theta|. There is no large-dimension entropy-continuity estimate
and no limit of the divergent vacuum entropy in this argument.

For FA91's unit normalized shells at scale L, put rho=R/L<=1/8 and
g_scaled(y)=R^(3/2)g_L(Ry). Homogeneity of Pi and exact cell scaling reduce
the problem to the unit ball with eta=a/R. Their gradient norm is <=5/L,
so ||g_scaled||H1<=sqrt(1+25rho^2)<2. A fixed B=6 is sufficient.
The lattice number operator always uses its actual normalization r_a.
Choose epsilon=rho^5, theta=rho^5 and eta by (12). Then the lattice and
finite-entropy errors are together <=2rho^5. Divide by

    N_R=(8pi^2 R^4/15) gamma_L,
    gamma_L=c_gamma/(pi^3 L^4), c_gamma=4249388/150088295.

The normalized extra error is at most

    C_E rho,       C_E=165/(14 c_gamma),                  (15)

using pi<=22/7. The continuum ball Taylor error from FA91 is <=C_B rho^2,
where C_B=(3/14)*195689447424*64/c_gamma. For a rest ball every original
null-direction shell has normalized limit 1. In particular all matter
kernel/refinement and perturbation errors are o(R^4) for fixed L and
rho->0. This is an actual common shrinking family for these matter errors;
it does not assert that geometric chart or area errors have been supplied.

## 6. Observer-frame preparations and the full response matrix

For the nine FA91 frames u_i=(5/4,3n_i/4), boost rapidity is log 2.
Let U_i be the free signed-CAR Cauchy-data Lorentz unitary chosen so that
the rest-ball response of U_i g_j equals the boosted-ball response of g_j.
It preserves the filled projector and L2 norm. On the two null momentum
sheets the momentum magnitude changes by factors in [1/2,2], so
||grad U_i g_j||<=2||grad g_j||. Free finite propagation also gives

    support(U_i g_j) subset {L/2<=|x|<=4L}.

Indeed the transformed initial surface samples |t|+|x|<=2|y| and
|x|-|t|>=|y|/2. The original shell solution vanishes if |x|+|t|<L or
|x|-|t|>2L. Thus these are genuinely exterior unit Cauchy spinors, and B=6
still suffices after scaling with rho<=1/8. The same proof applies to all
81 pairs with one common eta, epsilon, theta.

Their normalized actual finite entropy matrix converges entrywise to

    D_ij=(5/4+(3/4)n_i.n_j)^2,

with error <=C_E rho+C_B rho^2. FA91's exact inverse norm is68/5; hence
the Neumann condition 9*(68/5)*(C_E rho+C_B rho^2)<1/4 gives full rank and
inverse norm <=272/15. For example rho<=2^(-32) satisfies this with exact
rational arithmetic. All finite rows use ordinary rest-ball covariances
in the same source vacuum but distinct observer-frame shell preparations.
We do NOT identify a sharp boosted continuum algebra with a lattice lens,
nor establish a finite lattice boost symmetry or common geometric graph.

The geometric modular input is the one already used in FA87:
[La Piana--Morsella, Theorems 3.4/4.1](https://arxiv.org/html/2405.07888v2).
The CAR modular description is the existing FA96 import of
[Araki (1970/71), Corollary 4.10, Lemma 4.11 and Theorem 3](https://doi.org/10.2977/prims/1195193913).
All new endpoint,
refinement, fractional-boundary and finite-entropy estimates are proved
above. EinsteinRegisterData.generalizedEntropyStationarity,
fixedVolumeAreaExpansion, centralChargeIsEdge, smallBallNewton_pos and
the full common physical tower are not constructed here.

## 7. What is executed and what is proved analytically

The exact calculator verifies a full four-mode CAR reduction, the projection
block identity, the rank-one covariance and fourth-order Wick response, the
nine-direction matrix and the rational error budgets. Finite numerical
controls test the ordinary density-matrix entropy identity, the wrong hole
sign and wrong mode normalization, endpoint behavior and every valley frame.
Those finite controls are not a numerical execution of the enormous lattice
family (12). Existence and its rate are supplied by the analytic proof.
No native entropy-measurement circuit is claimed. Each finite regional
matrix is nevertheless a definite, finite object of the source, not a new
postulated continuum state.

The receipts bind the actual mathematical sources narrowly. They neither
rewrite predecessor receipts nor attest to a new full-repository CI replay.
Run `python -W error -m pytest -q model/shrinking_modular_entropy`, then
`python -m model.shrinking_modular_entropy --write` and `--check`.
