# A canonical continuum field sector from the ordered spin ground

The same ordered protected-spin branch has a nontrivial regular continuum Weyl-field net with positive energy and a finite propagation cone. This constructs a canonical field algebra and its full bounded time histories. It does not establish a sharp particle shell, Lorentz covariance or a surviving connected interaction.

Source: H=J sum_<xy>(Sx.Sy-1/4), spin 1/2 on Z^3, in the same infinitesimal-axis-1-staggered-field-selected ground omega as [the ordered ground](../spin_goldstone_sector/README.md) and [its source scaling](../spin_particle_scaling/README.md). Its actual order parameter m is not numerically identified, but m>1/1000. The statements below use the exact ground, not a product Neel vector.

## 1. Maximal order is a scalar macroscopic observable

Let gamma be the global pi rotation about internal axis 2, and beta_x=tau_x gamma^(x1+x2+x3 mod 2). These commuting automorphisms form Z^3; omega is beta-invariant; beta_x(S0^1)=(-1)^x Sx^1. On the observable algebra gamma^2=id. Beta commutes with the zero-field dynamics.

In the GNS representation, the mean ergodic theorem applied to the beta implementing unitaries makes the full-cube averages M_L converge on Omega. Their commutators with every fixed local A vanish in norm, and ||M_L||<=1/2. Therefore they converge strongly on the dense local vectors A Omega, hence strongly everywhere, to a beta-fixed central selfadjoint Z. In particular omega(Z)=m.

Here is a proof that the central decomposition preserves groundness. The closure of Center(pi(A)'') Omega is a reducing subspace for U(t)=exp(itK), K>=0. On it the map C(Z Omega)=Z* Omega is antiunitary: centrality makes the two norms equal and makes the polarization relation hold. It commutes with U(t) as an antilinear map, so C K C=-K. The restricted positive K is consequently zero. Thus every central spectral projection P has K P Omega=0, and its normalized conditional state is a ground state. Beta-fixed P also preserves beta invariance.

[Koma, arXiv:1712.09018v4, Theorem 2.1](https://arxiv.org/html/1712.09018v4) explicitly identifies m with the supremum of staggered magnetization over ALL infinite-volume zero-field ground states (not just field-selected states). If Z had nonzero spectral support above m, conditioning there would produce a beta-invariant ground state whose staggered mean exceeds m. Hence Z<=mI. Since omega(Z)=m and Omega is separating for the center, Z=mI.

For any real Schwartz h, the same conclusion for weighted Følner sums follows by step-function approximation and uniform l1 tail bounds:

    M_p(h)=p^3 sum_x h(px)(-1)^x Sx^1 -> m integral_R3 h(x) dx I

strongly as p decreases to zero. The integral and sums use no Fourier normalization. No convergence rate for this strong limit is claimed.

The argument also proves zero transverse Bragg mass if needed for compact spatial test functions. Average the full staggered vector over even translations and the two sublattices; all three limits are central commuting operators. The first is mI. If a conditional state had a nonzero transverse mean, a uniform internal rotation would turn its vector (m,z2,z3) into first-axis magnitude greater than m, contradicting the same maximum theorem. Thus both transverse central averages vanish. The spectral theorem for translation averages identifies their squared norms with the transverse structure-measure atom at Q.

## 2. Exact operator-norm energy lemma

Take a real absolutely summable single-axis angle theta and A=sum_x theta_x Sx^a. Finite support first suffices. The exact current identity gives

    C_A := ||[A,[K,A]]|| <= (J/2) sum_<xy> (theta_x-theta_y)^2.

K is the positive ground GNS Hamiltonian; the scalar energy subtraction has no effect. The same bound holds after conjugation by exp(itA). Consequently

    exp(-itA) K exp(itA) = K + t D + R_t,
    D=i[K,A],    ||R_t|| <= t^2 C_A/2.

This is an operator equality on Dom K, extended to quadratic forms. Positivity of the left hand side for every real t gives, for every normalized finite-energy psi,

    |<psi,D psi>| <= sqrt(2 C_A q(psi)),  q(psi)=||K^(1/2)psi||^2.

At t=1 the upper remainder bound therefore yields

    sqrt(q(exp(iA)psi)) <= sqrt(q(psi)) + sqrt(C_A/2).       (E)

The same holds for the inverse and time translates of exp(iA). Finite-support approximation extends this to summable theta with summable edge differences. This is a uniform energy-form estimate, not an expectation in Omega only.

## 3. Fields and equal-time Weyl relation

For real Schwartz F,G define

    Q_p(F)=p^2 sum_x F(px)(-1)^x Sx^3,
    P_p(G)=-p sum_x G(px) Sx^2.

Then [Q_p(F),P_p(G)]=i M_p(FG). The same-axis families commute exactly. The source frame gives Sx^a=(-1)^(sum x-x_a) sx^a. Thus the uniform S2 part has original carrier (pi,0,pi), and the staggered S3 part has original carrier (0,0,pi). Both reduce to (0,0,pi) modulo the actual unbroken physical translation lattice 2Z x Z x Z. They form one canonical pair in that residual valley. The spin indices remain internal; this does not assign Lorentz spin or helicity. A second pair cannot simply be appended: the commutator of its two uniform momenta contains longitudinal spin fluctuations for which the present estimates give no sufficient rate.

A useful exact bounded-matrix inequality, valid for all real a,b, is

    ||exp(iaS2)exp(ibS3)exp(-iaS2)exp(-ibS3)
      -exp(-iabS1)|| <= (|a|^2|b|+|a||b|^2)/4.            (W)

Proof: conjugating B=ibS3 by exp(A), A=iaS2, gives B+[A,B]+R with ||R||<=|a|^2|b|/4 by the second-order integral remainder. Comparing exp(B+C)exp(-B) with exp(C) costs at most ||[B,C]||/2<=|a||b|^2/4. All exponentials are unitary, so no exponential norm factors are inserted.

Apply (W) independently at sites and telescope tensor products. The reordering error for Q_p(F),P_p(G) is at most

    (p^4 sum_x |G(px)|^2 |F(px)|
      +p^5 sum_x |G(px)| |F(px)|^2)/4 = O(p).

Moreover for fixed real h,

    ||[M_p(h),P_p(G)]|| <= (p^4/2) sum_x |h(px)G(px)|=O(p),
    ||[M_p(h),Q_p(F)]|| <= (p^5/2) sum_x |h(px)F(px)|=O(p^2).

Consequently the strong scalar limit from section 1 remains valid on every finite equal-time Weyl word applied to Omega. With scalar symmetric phase convention

    W_p(F,G)=exp(i m integral FG/2) exp(iQ_p(F)) exp(iP_p(G)),

every subsequential all-word limit obeys

    W(F,G)W(F',G')
      =exp[-i m (integral FG'-integral GF')/2]
       W(F+F',G+G').

## 4. Regular vacuum state and source-exact energy cost

First choose test functions with smooth compact Fourier support away from zero. All moment estimates then pass from the retained finite-field reflection-positive inequalities without an infrared endpoint assumption. With normalized Fourier integral d^3q/(2pi)^3, the existing uniform and staggered transverse static bounds imply

    limsup omega(P_p(G)^2) <= (1/4) integral |q| |Ghat(q)|^2,
    limsup omega(Q_p(F)^2) <= 3 integral |Fhat(q)|^2/|q|.

Here e(k)<=|k|^2/2, e(k+Q)>=5, e(k)>=|k|^2/3 for |k|<=1 were used. These are upper bounds, not equalities or a Gaussian assertion. The zero transverse Bragg mass in section 1 allows extension of the pointwise measure domination through Q. Poisson summation and Schwartz tails then give the same limiting bounds for compactly supported smooth spatial functions, and more generally Schwartz functions. Alternatively the equal-time regular state can be extended in these weighted Sobolev norms from the Fourier-cut test class.

The norms of (exp(itQ)-1)Omega and (exp(itP)-1)Omega are at most |t| times these second-moment bounds. The limiting state on the equal-time Weyl algebra is therefore regular and has nondegenerate symplectic coefficient m>1/1000. Dividing both fields by sqrt(m) makes the canonical coefficient one, but does not calculate m. A zero field limit is excluded by the commutation relation.

The exact energy lemma gives, for one twist,

    cost_P,p = (J p^2/4) sum_<xy> [G(px)-G(py)]^2,
    cost_Q,p = (J p^4/4) sum_<xy> [F(px)+F(py)]^2.

Their limits divided by Jp are respectively

    (1/4) integral |grad G|^2,      3 integral F^2.

For a finite product of arbitrary time translates of these twists,

    sqrt(q(word Omega)) <= sum_j sqrt(cost_j).

Thus every fixed finite multitime word has energy O(Jp). All time translates preserve energy. This proves uniform rescaled-time equicontinuity of all word correlation functions: to change the time of one factor use

    ||(U(t)-1)psi||^2 <= 2|t| q(psi)

on the two relevant word vectors and telescope. In physical scaled time t/(Jp), all resulting bounds are uniform in p.

## 5. Propagating the central phase to multitime words

This step is necessary: a strong scalar limit on fixed vectors alone cannot be applied to p-dependent histories.

A general consequence of (E) is as follows. Suppose a unitary V and its inverse change energy square roots by at most delta. Put K'=V* K V. Then for normalized finite-energy psi,

    |<psi,(K'-K)psi>| <= 2 delta sqrt(q(psi))+delta^2
      <= epsilon q(psi)+delta^2(1+1/epsilon).

Thus D=K'-K is bounded in both signs by the positive quadratic form G_epsilon=epsilon K+delta^2(1+1/epsilon). Its polarized form satisfies |D(phi,psi)|<=sqrt(G_epsilon(phi) G_epsilon(psi)).

For q(psi)<=E, both exp(itK)psi and exp(itK')psi have K-energy at most (sqrt(E)+2delta)^2. Differentiate their overlap and choose epsilon=delta/(sqrt(E)+2delta) (delta=0 is trivial). Integration gives

    ||(exp(itK)-exp(itK'))psi||^2
       <=4|t| delta sqrt(E)+10|t| delta^2.                (C)

Domain approximations justify the overlap differentiation for energy-form vectors. For V=exp(is M_p(h)), delta=O(sqrt(J) p^(3/2)) by (E). With E=O(Jp), t=T/(Jp), the right hand side is O(p). Hence this macroscopic phase asymptotically commutes with time evolution on every finite-history vector. Combined with the equal-time O(p) commutator, it can be moved through every fixed multitime word to Omega. Section 1 now proves that the phase is scalar throughout the full limiting history representation. Therefore the same-time Weyl relations are operator relations there, not just vacuum equalities.

## 6. Full-history regularity and spatial continuity

There is a useful second application of (C). Take V_s=exp(isQ_p(F)) or exp(isP_p(G)). Its energy kick is bounded by delta=|s| kappa sqrt(Jp), where kappa depends on the fixed test function. For a finite-history vector of energy at most C Jp, and time T/(Jp), equation (C) gives a squared commutation error at most

    4 |T| |s| kappa sqrt(C) + 10 |T| s^2 kappa^2.

This tends to zero with s, uniformly in p. Move V_s through a fixed multitime word one factor at a time. Passing its time transports costs this bound; passing an equal-time Weyl factor costs the scalar symplectic phase, which tends to one with s, and the already proved vanishing finite-p remainder. Once V_s reaches Omega, its norm difference from identity is O(|s|) by the second-moment estimates. Therefore every limiting Weyl one-parameter group is strongly continuous on the dense set of history vectors, hence on the full GNS space. Regularity has not been inferred merely from the vacuum at one time.

The same argument applies to a test increment instead of an amplitude. Its small parameters are the energy seminorms ||grad dG||_2 and ||dF||_2, together with the vacuum norms ||dG||_(homogeneous H^(1/2)) and ||dF||_(homogeneous H^(-1/2)). The symplectic phases are continuous on compact smooth test families. Consequently the Weyl fields are strongly continuous under smooth translations of their tests on every finite-history vector. Lattice translations with even parity can approximate any fixed continuum shift; the bounded rounding error changes the scaled tests by O(p) in all these norms. The rounded translations therefore have an exact group law in the limit, and give strongly continuous spatial translations. Equivalently, beta translations give the same continuum envelope action after removing the fixed carrier convention. This does not imply spatial rotation or boost covariance.

## 7. Positive-energy local field construction

Choose a countable rational span of smooth compact tests, closed under a countable dense translation set, and rational times. All finite words are products of bounded unitaries. Compactness of their expectations yields a common subsequence; Gram positivity gives the GNS representation. Section 5 proves the same-time Weyl relations in this full representation, and section 6 proves their strong regularity. The identity has expectation one and m>1/1000, so the canonical algebra is nonzero. Neither a Gaussian state nor a quantum central-limit theorem is assumed.

The energy estimate in section 4 gives continuity of every time-translated word: for a changed single factor use ||(U(t)-1)psi||^2<=2|t|q(psi) on the two adjacent history vectors, then telescope. Complete the rational times and tests by this continuity. Stationarity gives a strongly continuous implementing unitary group. Its generator is nonnegative: each matrix coefficient of a finite word has spectral support in [0,infinity) under K/(Jp); testing against a smooth function supported at negative frequency passes this property to the limit. Spatial and time translations commute, as their lattice predecessors do.

For compact spatial supports the retained spin Lieb-Robinson bound is

    ||[alpha_t(A_X),B_Y]|| <= ||A|| ||B|| |X| 2^(-d) (exp(48J|t|)-1).

Here d is graph distance and the field unitaries have norm one. In scaled coordinates X=p*x and T=J*p*t, a positive separation d_1(O1,O2)-72|T1-T2| makes the bound tend to zero: |X|=O(p^-3), whereas the exponential is at most 2^(-[d_1-72|T|+o(1)]/p), using log(2)>=2/3. Thus the limiting Weyl algebras commute outside this finite l1 cone. Isotony follows by taking the algebras generated by tests supported inside each region. The chosen constant is a rigorous upper bound and is not an observed light speed.

This is a regular local bounded Weyl-field net in three continuous spatial coordinates and continuous time, constructed from the actual ordered spin ground and dynamics. Each regular one-parameter Weyl group has its Stone field generator. A common invariant domain carrying every unbounded field product, a Wightman reconstruction, a sharp mass shell, Lorentz covariance and nonzero connected continuum scattering are not claimed. In particular canonical commutation relations alone do not distinguish an interacting state from a Gaussian state.

## 8. Source realization and exact checks

The construction first uses the same fixed-spacing infinite spin interaction proved by [the protected-spin limit](../protected_spin_limit/README.md), then its long-distance/long-time scaling. It does not interchange this with an unproved simultaneous Hubbard ground-state limit. The known finite source comparison and arbitrary finite precision preparation can approximate each fixed finite operation; there is still no effective box-size/field-selection modulus for the ordered state. No finite native run is presented here as the infinite state.

`algebra.py` checks all three exact two-site double commutators, the transverse bond norm, the canonical sign and the nonzero cross-commutator that prevents declaring two independent pairs. `bounds.py` keeps the Weyl and energy constants rational and evaluates the central-phase transport bound. The examples normalize the energy and kick coefficients to one; they are error-bound examples, not measured excitation energies. Tests also check the general energy lemma and Weyl remainder on finite matrices, and retain the countercontrol that positive mean order without maximality permits a nontrivial central mixture.

The receipt binds the current mathematical input texts and the new source bytes. It explicitly does not pretend to rerun every inherited producer. The infinite theorem is proved above; the finite checks audit its algebra and quantitative constants.

```sh
python -W error -m pytest -q model/spin_canonical_field
python -W error -m model.spin_canonical_field --write
python -W error -m model.spin_canonical_field --check
```
