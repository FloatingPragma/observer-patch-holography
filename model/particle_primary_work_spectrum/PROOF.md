# Original work spectra, modular fluctuations and a criticality boundary

The same selected primary source obeys exact modular fluctuation and
affinity identities. They strengthen its downward-energy tail and give
positive finite distinguishing bounds for the critical local phase
apparatus. A full-CAR commutator estimate also improves its selected second
moment cap. That cap grows as the inverse fourth power of spacing. Neither
the modular identities nor the improved cap proves surviving probability
at fixed positive energy.

## 1. Complete original transition law

Use the source, polynomial density, idle reference and source clock of the
[primary contract](../particle_common_contract/PROOF.md). The
[critical unit-tent apparatus](../particle_primary_critical_refinement/PROOF.md)
has unitaries U_(X,sigma), one for each origin and sign, each weighted by
1/(2V). Its adjoint is the opposite sign at the same origin. Every label,
source output and reference is retained. Let E_i be the finite physical
energies, p_i=p(E_i)>0 their original selected populations, and K_i=-log p_i.
The complete ideal two-energy law is

    P(i,X,sigma,j)=p_i |(U_(X,sigma))_ji|²/(2V),
    omega=E_j-E_i,       kappa=K_j-K_i.                         (1)

Initial energy recording preserves the source marginal and generally
changes its reference correlations. The formulas concern that marginal's
energy-record probabilities, not an unchanged joint purification.
Reversing i,j and sigma gives the exact path ratio

    P(reverse path)=exp(-kappa) P(path),
    E exp(-kappa)=1.                                          (2)

The second equality follows directly from both row and column sums of a
unitary modulus-square matrix. There is no Gibbs approximation. On the
actual spectral interval the
[selected polynomial slope](../particle_primary_gravity_premises/PROOF.md)
satisfies 0<beta_low<=K'(E)<=beta, with beta_low>=2beta/3. Thus a positive
transfer has beta_low omega<=kappa<=beta omega. The exact identities imply

    beta_low W <= E kappa <= beta W,    W=E omega.              (3)

Pair each positive transfer with its reverse when proving (3). Each pair
has nonnegative population difference; a pointwise slope inequality cannot
simply be integrated across both signs without this pairing.

The energy-index marginal kernel T_ij=average |U_ji|² is symmetric and
doubly stochastic. With uniform reference pi_i=1/d and q=push(p,T), all
premises of `integral_fluctuation` and `crooks_pointwise` in
[FluctuationTheorems.lean](../../Lean/Thermodynamics/FluctuationTheorems.lean)
are supplied: p,pi,q are positive, both probability laws normalize, and
T preserves and is in detailed balance with pi. Their entropy-production
variable is sigma_EP=log p_i-log q_j. It is different from kappa in (1):

    E sigma_EP=H(q)-H(p),
    E kappa=H(q)-H(p)+D(q||p).                                (4)

Here the source energy-index law is a specified marginal of the retained
experiment. It is not the entropy of all apparatus and reference outputs.
The elementary source attachment supplies the hypotheses of the finite
theorems; it does not formalize the full primary Hamiltonian in Lean or
identify the unitary intervention with conditional repair.

## 2. Sharper actual downward tails

For a positive energy pair put r=p_j/p_i<=exp(-beta_low omega).
The fraction of its total forward-plus-reverse weight carried by the
negative transfer is r/(1+r). Since all pair weights together are at
most one, for every E>0,

    mu(omega<=-E) <= 1/(1+exp(beta_low E)).                     (5)

This is an original complete-probability bound. For exact rational
reporting, e>2 gives the outward bound
1/(1+2^floor(beta_low E)); the potentially enormous exponent is stored
in factored form rather than evaluated or rounded to zero.

The elementary inequality 1+exp(x)>=3x for x>=0 follows from
exp(x)>=1+x+x²/2. Therefore

    integral_(omega<0) |omega| dmu <= 1/(3 beta_low),
    integral |omega| dmu <= W+2/(3 beta_low).                   (6)

The receipt also retains the independently proved ground-excess bound and
takes the smaller upper bound. The first member has beta_low>22271, so
(6) is much smaller than its ground-excess budget1/16. These estimates
control the negative tail; they impose no positive energy-weighted tail
bound or fixed positive-energy residue.

## 3. Affinity and finite distinguishability

On source marginals define

    A=average Tr(sqrt(rho) U sqrt(rho) U*)
      =sum_paths P(path) exp(-kappa/2),       0<=A<=1.         (7)

This is both the quantum affinity between the source classical-quantum
output and the same-label no-operation source, and the classical affinity
between (1) and its reversed path law. It is not the overlap of the
actual retained canonical purification before and after U.

Let Delta be the authenticated upper minus lower source spectral bound.
For a positive energy pair its contribution to 1-A is
(sqrt(p_i)-sqrt(p_j))² times its transition coefficient, whereas its work
is (p_i-p_j)omega times that coefficient. Their ratio is
tanh(kappa/4)/omega. From exp(2x)>=1+2x one gets
tanh(x)>=x/(1+x), x>=0. Hence

    1-A >= beta_low W/(4+beta_low Delta).                      (8)

The trace distance T_CQ, defined as half the trace norm between the two
source classical-quantum marginals, also satisfies

    T_CQ >= W/Delta.                                          (9)

Indeed the Hamiltonian, after subtracting its interval midpoint, has
norm at most Delta/2. Apply its trace duality bound on each retained
origin/sign block and average. Both (8) and (9) have strictly positive
actual source lower bounds from the critical-profile work result.
They deteriorate with capacity because Delta grows. Taking a partial
trace shows that (9) is also a lower bound on the distance of the actual
retained experiments; no reference-inclusive upper bound is inferred.

There is a complete finite-record comparison as well. The inherited
two-clock plan has norm cap R_H and signed-mean error e_mean. Compare the
phase apparatus with a no-operation apparatus using the identical source,
labels and clocks. The latter has zero signed mean because its identical
clock kernels commute. The reported transfer range is at most56R_H, so

    TV(phase records, no-operation records)
        >=(W_lower-e_mean)/(56R_H)>0.                         (10)

The numerator, clock sizes and implementation budgets are all taken from
the frozen primary work consumer. This proves a finite observable
difference without a full histogram or a useful capacity-uniform rate.

## 4. A stronger full-CAR second-moment cap

The exact original second moment is

    M2=average Tr rho ([H,U]*[H,U])
      =average Tr rho (U*HU-H)².                              (11)

The phase commutes with electric, onsite and magnetic terms. Before the
product link and physical Gauss compressions, U*HU-H is multiplication
in the gauge coordinates by the full-CAR second quantization dGamma(t).
The one-body Hermitian matrix t has zero trace. For R=1/a>=2 its nonzero
bonds touch at most

    M=(2R-1)^3+6(2R-1)^2=(2R-1)^2(2R+5)                      (12)

source sites, so rank(t)<=60M. Missing open-boundary edges can only
reduce this count. The periodic coefficient list has the same bound
for wrapped origins. At member zero the source corner star has M=4.

The literal source matrix S_i=(-beta+beta gamma_i)/2 has rank two and
two singular values one. Its full fifteen-state internal representation
is unitary at every gauge coordinate. If delta_b is the phase change
on one bond, the two orientations therefore give exactly

    Tr(t²)=120 a^(-2) sum_b [1-cos(delta_b)].                  (13)

For a finite trace-zero Hermitian t, the full-Fock norm of dGamma(t)
is half the sum of absolute eigenvalues. Cauchy--Schwarz yields

    ||dGamma(t)||² <= rank(t) Tr(t²)/4.

This is an ambient operator bound, not an expectation in a Haar average
or in the product anchor. The phase commutes with both source projections,
so compression preserves the bound on the actual Gauss carrier. Averaging
origins uses the exact same k_n as the selected work identity, giving

    M2 <= 1800 M B k_n/(a²V)
        =1800 M B W/(a V h_n).                               (14)

With the exact discrete gradient sum D2 and pi²<10, k_n<=5a²D2/3, hence

    M2 <=3000 M(B/V)D2.                                     (15)

For R>=2 this improves the inherited324000R^4 bound, but remains
O(a^(-4)). The member-zero estimate keeps the smaller inherited
full-star cap19620. Write w>0 for the selected work lower bound and C
for the smaller second-moment cap. Below w/4 the signed first moment
is at most w/4, including every negative transfer. Above 4C/w it is
at most C/(4C/w)=w/4. Therefore the middle band's first moment is
at least w/2, and Cauchy--Schwarz gives the original probability

    mu([w/4,4C/w]) >= w²/(4C).                               (16)

The identical critical-profile two-clock plan has resolution eta,
one-clock tail tau and complete implementation error e_clock. Its
reported event on [w/4-2eta,4C/w+2eta] has original probability at least
w²/(4C)-2tau-e_clock. The actual receipt checks that the lower endpoint
is positive and the loss is less than half the ideal floor. All origin,
sign and clock records, the source and reference are retained; the
reported probability is not normalized on success. The ideal member-zero
floor exceeds2.58e-5. The floor decays as a^4 while the upper endpoint
grows as a^(-4); this is a finite observable, not a fixed-energy survival
theorem. None of these norm bounds supplies the selected current-correlation
cancellation needed for uniform integrability.

## 5. A faithful polynomial control for the entropy shortcut

The following is a different two-level carrier, not the primary source.
It tests the implication from polynomial cooling, modular identities,
positive work and a uniform phase energy-form bound to fixed-energy
survival. Copy the literal primary polynomial parameters beta_n,m_n,s_n
and upper endpoint U_n. On

    H_control=diag(U_n-1/a,U_n)

apply exactly Q(E)=T_m(beta(U_n-E)/(2s))^s and normalize Q² on this
two-level carrier. Both eigenvalues lie in the original source polynomial
interval. The positive population ratio r obeys

    exp(-beta/a)<=r<=exp(-beta_low/a),
    mean ground excess <=1/beta_low<=epsilon_n.               (17)

This state is faithful and is not replaced by an exact ground state.
Choose F=sqrt(a) sigma_x and U=exp(iF), with both signs retained.
The linear symmetry C=sigma_z commutes with H_control and its polynomial
state and interchanges the two phase signs. Exact calculation gives

    ||[F,[H_control,F]]||=2,
    mu=(1-sin²(sqrt(a))) delta_0
       +p_ground sin²(sqrt(a)) delta_(1/a)
       +p_excited sin²(sqrt(a)) delta_(-1/a),
    W=(p_ground-p_excited) sin²(sqrt(a))/a ->1,
    M2=sin²(sqrt(a))/a² ~1/a.                                (18)

All modular fluctuation, detailed-balance and affinity identities in
Sections1--3 hold exactly. The mean modular work divided by beta has
liminf at least2/3, while A tends to one. Every fixed positive-energy
band has zero weight for sufficiently small a. Since ||U-I||<=sqrt(a),
all fixed finite time-translated unitary histories tend to identity in
operator norm, even with a retained reference. Multiplying F by any
fixed real amplitude gives the quadratic limiting work for that amplitude.

Thus faithfulness, the actual polynomial slope, positive work, and the
common energy-form estimate do not repair the missing implication.
This control proves no primary-source elastic limit. The primary selected
positive energy-weighted tail, or a directly populated fixed-energy
window, is the distinct source estimate required for that decision.
