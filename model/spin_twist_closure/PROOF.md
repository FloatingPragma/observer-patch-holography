# Exact source twists and dynamical closure

The source, state and clock are those of the ordered-spin construction:

\[
 H_J=J\sum_{\langle xy\rangle}(S_x\cdot S_y-1/4),\quad J>0,
 \qquad L_p=(H_J-E_0)/(Jp).
\]

On the infinite lattice, `H_J-E0` means the positive generator K of the
infinitesimal-staggered-field-selected ground. Finite local commutators do
not require a value of E0. The selected finite-volume field is removed in
the prescribed order before p tends to zero. The spin is one half and each
unoriented nearest-neighbor bond occurs once. The source derivation of this
interaction and the complete repair comparison are inherited from
`protected_spin_limit`. The fields are exactly

\[
 Q_p(F)=p^2\sum_x\epsilon_x F(px)S_x^3,\qquad
 P_p(G)=-p\sum_xG(px)S_x^2.
\]

## 1. An operator-norm quadratic expansion

For one internal axis a put

\[
 A_\theta=\sum_x\theta_x S_x^a,\quad
 T_{xy}^a=\sum_{b\ne a}S_x^bS_y^b,\quad
 C_{xy}^a=(S_x\times S_y)^a,\quad d_{xy}=\theta_x-\theta_y.
\]

The exact norms of T and C are both 1/2. Direct spin rotations give

\[
 e^{-iA_\theta}H_Je^{iA_\theta}-H_J
 =J\sum_{\langle xy\rangle}
       [(\cos d_{xy}-1)T_{xy}^a-\sin d_{xy}C_{xy}^a].       \tag{1}
\]

The fixed scalar bond shift cancels. The commutators are
`i[H_J,A]=-J sum d C` and `[A,[H_J,A]]=-J sum d^2 T`.
Taylor's integral remainder applied to each rotated bond, whose third
derivative has norm `J|d|^3/2`, proves

\[
 e^{-iA}L_pe^{iA}=L_p+D_p+\tfrac12 Z_p+R_p,
 \quad D_p={i[H_J,A]\over Jp},\quad
 Z_p=-{1\over p}\sum d_{xy}^2T_{xy}^a,
 \quad \|R_p\|\le {1\over12p}\sum|d_{xy}|^3.               \tag{2}
\]

This is uniform over all input states and untouched references. It is an
operator identity on Dom K, and hence a quadratic-form identity on the
form domain. Compactly supported tests have finite sums. Summable tests
with the displayed difference sums finite follow by local truncation.

For smooth tests supported in `[-R,R]^3`, take `R>=1`, `0<p<=1` and
`N0=(2R+1)^3`. At most `6N0 p^-3` bonds meet the support. If
`|partial_i G|<=g`, then the P twist has `|d|<=g p^2`. If `|F|<=f`,
the Q twist has `|d|<=2f p^2`. Thus the respective bounds are

\[
 \|R_{P,p}\|\le N_0g^3p^2/2,\qquad
 \|R_{Q,p}\|\le4N_0f^3p^2.                                \tag{3}
\]

The higher pure-direction BCH terms therefore vanish in norm on the
actual scaling. This conclusion uses neither a ground-state factorization
nor a sharp spectral branch.

At finite positive selection field the same computation also has its
onsite field term. Its third-order remainder is bounded by
`B sum|theta_x|^3/(12Jp)`. It is not dropped at fixed B. The prescribed
field limit removes it at each fixed p; no interchange of limits is used.

## 2. The exact coefficient is a bond observable

For a positively oriented bond `(x,x+e_i)`, (2) gives

\[
 Z_{P,p}=-p^3\sum_{x,i}
  \left({G(px+pe_i)-G(px)\over p}\right)^2T_{x,x+e_i}^2,
 \quad
 Z_{Q,p}=-p^3\sum_{x,i}[F(px)+F(px+pe_i)]^2T_{x,x+e_i}^3. \tag{4}
\]

They are bounded weighted averages of the actual source bond operators.
The selected state's cubic and internal-axis-one symmetries give
`omega(T^2)=omega(T^3)=-kappa`, so their expectations tend to

\[
 \kappa\int|\nabla G|^2,\qquad12\kappa\int F^2.           \tag{5}
\]

These are the coefficients of the established microscopic twist-energy
form. Equations (2)--(5) also identify the full coefficient operators,
whose fluctuations are not replaced by their expectations.

Each finite coefficient has a complete source readout. For
`Z=sum_b z_b T_b^a`, put `s=sum|z_b|` and choose setting `(b,c)` with
probability `|z_b|/(2s)` for each `c!=a`. Measure the two local Pauli
operators, retaining both signs `u,v`. The score
`(s/2) sign(z_b) uv` has expectation exactly `Tr(rho Z)`, on any source
state with any reference. The branch Kraus operator is
`sqrt(|z_b|/(2s)) Pi_(x,c,u) Pi_(y,c,v)`; summing every effect gives I.
This is a score for Z, with its specified measurement disturbance. It
does not claim the spectral instrument of the whole noncommuting sum.
For s=0 the complete instrument is the identity branch with score zero;
no setting normalization or repetitions are required.

The native Pauli compiler supplies the stated finite local projectors
through the admitted coherent services. An unhalved complete-output error
eta changes the mean score by at most `s eta/2`. N independent complete
repetitions have absolute sampling error epsilon except with probability
at most `s^2/(4N epsilon^2)`. No success branch is renormalized. The code
emits the complete rational setting/score plan; it does not report an
executed large detector tape or an infinite-ground preparation rate.

## 3. Bond phases can pass through the same time histories

Let `Z=sum_b z_b T_b^a` have real summable coefficients on a simple graph
of maximum degree Delta. A Hamiltonian bond intersecting one chosen bond
has at most `2Delta-1` choices. A further bond intersecting their union
has at most `3Delta` choices. Every double commutator term has norm at
most J, using `||T||=1/2` and `||S.S-1/4||=1`. Therefore

\[
 \|[Z,[H_J,Z]]\|
 \le3\Delta(2\Delta-1)J\|z\|_\infty\sum_b|z_b|.           \tag{6}
\]

For the cubic source the constant is 198. If `|z_b|<=a p^3` and
`sum|z_b|<=a V0`, the energy-kick lemma from `spin_canonical_field`
applies to `exp(isZ)` with
`delta<=|s|a sqrt(99JV0) p^(3/2)`.
For a normalized history vector of energy at most `E Jp`, its exact
conjugated-flow estimate at source time `T/(Jp)` is

\[
 \|(e^{iTK/(Jp)}-e^{iTK'/(Jp)})\psi\|^2
 \le4|Ts|a\sqrt{99V_0E}\,p+990|T|s^2a^2V_0p^2,
 \quad K'=e^{-isZ}Ke^{isZ}.                                \tag{7}
\]

All fixed finite canonical histories have the needed energy bound from
that same lemma. Equal-time commutation follows independently:
`||[Z,A_theta]||<=||theta||_infinity sum|z_b|`, giving O(p) for P and
O(p^2) for Q. Integrating the commutators gives the corresponding unitary
phase bounds. Equations (6)--(7) therefore carry these actual coefficient
phases through every fixed finite multitime history. They do not use a
uniform operator-norm approximation of the macroscopic time evolution.

## 4. Central limits and their precise representation scope

For each fixed local bond observable, the selected state is invariant
under the beta translations used in `spin_canonical_field`. Beta acts on
each same-axis bond product as ordinary translation. The mean-ergodic
argument for bounded spatial averages gives a strong limit on Omega.
Their norm-small commutators with every fixed local observable then give
a bounded central limit on the microscopic selected-state GNS. Weighted
smooth averages follow by step-function approximation and uniform l1
bounds. Denote these central microscopic bond limits by `Tbar_i^a`.
The coefficient limits on fixed microscopic vectors are consequently

\[
 Z_P=-\sum_i\overline T_i^2\int(\partial_iG)^2,\qquad
 Z_Q=-4\sum_i\overline T_i^3\int F^2.                       \tag{8}
\]

The coefficients of (4) differ from their fixed derivative/test weights
by an operator-norm error tending to zero, by Taylor bounds and Riemann
sum estimates. Central projections preserve the microscopic ground
condition by the positive-generator center argument in the inherited
canonical-field proof. Their time action on that ground representation
is trivial.

One may adjoin these bounded coefficient phases to the common source
history family before taking a subsequence. Equations (6)--(7) and the
equal-time bounds prove that their limits commute with every field
history and are time invariant in this enlarged representation. Their
microscopic joint law, including every central component's original
weight, is retained. The original minimal field-history representation
embeds in this construction. It is not assumed to reduce each added
coefficient or to contain all of their functional calculus.

If the selected microscopic state is translation ergodic, the relevant
bond averages are scalar. The stated selected-state theorems do not
establish that extra hypothesis for this chosen state. Maximality of
magnetization fixes its macroscopic magnetization; it does not identify
every central bond correlation. Equations (5) and (8) are therefore
distinct statements.

### The total bond energy is scalar without ergodicity

One aggregate coefficient is fixed more strongly by the actual ground
condition. Let Lambda be an L-by-L-by-L open cube, with V=L^3 sites,
`M=3L^2(L-1)` internal bonds and `b=6L^2` crossing bonds. Let H_Lambda
contain the internal shifted bonds, and let E_Lambda be its lowest
eigenvalue. Its finite rational spin matrix determines E_Lambda without
a unique-ground premise. For a normalized finite ground g, the complete
reset Kraus family `A_j=|g><j|` on Lambda satisfies
`sum A_j^dagger A_j=I` and
`sum A_j^dagger H_Lambda A_j=E_Lambda I`.

Summing the actual infinite-source ground inequalities
`omega(A_j^dagger[H_J,A_j])>=0` gives

\[
 0\le\omega(H_\Lambda-E_\Lambda)
 \le\omega(\Psi(H_\partial)-H_\partial)\le Jb.             \tag{9}
\]

Here Psi is the complete reset channel in the Heisenberg convention.
Each shifted source bond is a negative projector times J. Consequently
`-Jb I<=H_partial<=0` and the same bounds hold after the unital positive
map Psi, which gives the last inequality with J per boundary bond.
Also `0<=H_Lambda-E_Lambda<=JM`. Squaring this positive operator gives
the quantitative ground-vector bound

\[
 \left\|{H_\Lambda-E_\Lambda\over V}\Omega\right\|^2
 \le {J^2Mb\over V^2}
 =18J^2{L-1\over L^2}.                                   \tag{10}
\]

The thermodynamic ground density e0 has the explicit finite-box bracket

\[
 {E_\Lambda\over V}-{3J\over L}\le e_0
 \le {E_\Lambda\over V}.                                 \tag{11}
\]

To see this, tile a large cube by copies of Lambda. Crossing bonds have
operator bounds `-J<=h_b<=0` and density at most `3/L`. The lower bound
uses `H_Lambda>=E_Lambda`; the upper bound uses a product of finite
ground vectors and the nonpositive crossing-bond expectations. The
outer boundary disappears in the large-cube limit. The same reset
inequality (9) shows that every beta-translation-invariant source ground
has this energy density: its average local energy differs from the
finite minimum density by a surface-to-volume term. Thus the bracket
does not select a different state from the prescribed ordered ground.

Combining (10)--(11) proves strong convergence to e0 I on Omega, with
norm error at most
`J sqrt(18(L-1))/L+3J/L`. Norm-small local commutators extend convergence
to every fixed local vector. Including the forward boundary bonds costs
at most another `3J/L`. Since
`sum_a T_(x,x+ei)^a=2 S_x.S_(x+ei)`, the central limits satisfy

\[
 \sum_{i=1}^3\sum_{a=1}^3\overline T_i^a
   =\left({2e_0\over J}+{3\over2}\right)I=-2eI,           \tag{12}
\]

where e is the unshifted source energy coefficient bounded by the
entangled-dimer calculation. Equations (6)--(7) carry this identity to
the enlarged field-history representation as well. This proves a scalar
aggregate and a quantitative microscopic ground-vector convergence
bound without asserting scalarity of each transverse coefficient. It
does not give a trace-distance rate for every finite-volume ground
marginal; (10) concerns the selected infinite state's local energy.

## 5. Quadratic twists do not supply a closed field equation

Even scalar pure-direction BCH coefficients do not make the observed
canonical algebra dynamically invariant. A positive exact control is

\[
 H={1\over2}(p^2+\pi^2+q^2+x^2+qx),\qquad[q,p]=[x,\pi]=i.
\]

Translation of q or p gives exactly a polynomial of degree two, with
scalar second coefficient one and zero higher pure-direction terms.
Nevertheless `ddot q=-q-x/2`. The two positive normal-frequency squares
are 1/2 and 3/2, and q has nonzero residue at both in the unique Gaussian
ground. The hidden force x is not a function of the observed q,p pair.
This is a logical countercontrol, not an alternate spin-source model.

The actual source has an exact finite nonlinear force as well. On the
three-site path 0--1--2 and endpoint observable `S_0^a`, define

\[
 F^a=J^{-2}[H_J,[H_J,S_0^a]]-(S_0^a-S_1^a)/2.
\]

It is trace orthogonal to every one-site spin and has normalized
trace square 1/16 on all three axes. The code verifies all of these
identities. Thus replacing the finite source equation by a closed linear
spin equation is false. This finite tracial statement does not decide
whether the force vanishes or linearizes on the actual low-energy ordered
histories. The source-specific filtered force estimate in
`spin_shell_diagnostic` is exactly the stronger condition needed there.

Passing first derivatives of time histories or microscopic energy forms
through a weak bounded-history limit also needs its own domain/tail
justification. Vanishing R_p alone gives neither that interchange nor
the reducing property of the equal-time field algebra. The construction
proves exact quadratic source twists, operational coefficient scores and
coefficient-phase transport. It supplies no sharp shell, common Lorentz
speed, surviving scattering amplitude or physical species identification.
