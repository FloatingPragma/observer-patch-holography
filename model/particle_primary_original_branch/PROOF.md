# Original energy records with conserved source labels

## 1. The selected polynomial fixes the forward/reverse ratio

Use member zero of the [primary contract](../particle_common_contract/PROOF.md):
the physical Gauss carrier of the interacting open Wilson source, its exact
selected polynomial density rho, retained idle reference and source clock.
The [coherent odd apparatus](../particle_primary_coherent_odd_response/PROOF.md)
contains the actual six weak-pair/singlet fields F_(x,j), with two G1
components and four H components. Their native occupation increment is
Delta=(0,2,0,0,1). All gauges, sixty native modes per cell, product link
cutoffs and failures are retained.

On the physical carrier write rho=p(H)>0 and K(E)=-log p(E). The
[exact polynomial slope](../particle_primary_gravity_premises/PROOF.md#2-exact-modular-slopes-of-the-selected-polynomial)
is derived from the actual selected source polynomial:

    Q(E)=T_m(beta(U-E)/(2s))^s,
    K'(E)=beta T_(m-1)(y(E))/T_m(y(E)),
    beta_low=beta(1-delta) <= K'(E) <= beta.                 (1)

Here beta=22272, m=11 and delta is the exact rational endpoint bound
returned by the source preparation. The program has beta_low>22271.
No Gibbs density replaces rho. Beta is the specified cooling parameter;
this calculation does not identify a proper temperature or a modular
clock with the source clock.

For any source operator A define the original creation transfer measure

    mu_A = sum_(i,k) p(E_i)|A_ki|² delta_(E_k-E_i).          (2)

This is the energy-record law of the stated intervention with the
original state weights. Degeneracies may be resolved in any energy basis.
At a positive transfer omega, exchanging i and k relates the adjoint
measure to (2), with coefficient p(E_k)/p(E_i). Integrating (1) gives

    exp(-beta omega) dmu_A(omega)
      <= dmu_(A*)(-omega)
      <= exp(-beta_low omega) dmu_A(omega),  omega>0.       (3)

The inequalities hold as inequalities of positive measures, so they
also hold after summing components with the same apparatus normalization.
This is a paired-adjoint detailed-balance bracket, not an exact Gibbs
identity. It uses source probabilities; it does not assert domination
of the retained source-reference density by a different purification.

A* annihilates the positive-beta composite. It is not the literal
source-C conjugate negative-beta hole creator. The latter has the same
positive-energy law as A by the source C symmetry, whereas (3) suppresses
the reverse annihilation event. No particle/antiparticle abundance
asymmetry follows from (3).

## 2. Controlled positive energy weight at every spatial phase

Let A_(f,j)=sum_x f_x F_(x,j), with |f_x|²=1/8. The complete selected-setting
apparatus has successes A_(f,j)/sqrt(6), A_(f,j)*/sqrt(6), and the positive
square-root failure. Its original signed success weight is at least w,
with w>0.84402. The actual compressed Wilson residual gives signed
second-moment cap110/3 around centers+9 and-9. Therefore the union

    J_+=[1/2,35/2],    J_-=-J_+

has original signed weight at least w-440/867. These statements are
[full-CAR consequences](../particle_primary_coherent_odd_response/PROOF.md)
on the same input state, not a projection to its empty occupations.

Equation (3) supplies mu_-(J_-)<=r mu_+(J_+), where
r=exp(-beta_low/2). Consequently

    mu_+(J_+) >= (w-440/867)/(1+r).                        (4)

For exact arithmetic, e>2 implies e^(-x)<=2^(-floor(x)) for x>=0.
The authenticated slope gives floor(beta_low/2)=11135. The receipt
retains this exponent and uses the weaker rational r<=2^(-16) in (4).
No floating exponential or underflow enters the bound.

The same two complete source clocks have resolution1/8 and total event
loss at most1/4096. Their reported interval [1/4,71/4] has original
selected-setting creation weight greater than1/3. For each individual
G1 or H multiplet, use the respective signed success floor d_lambda w/6
and signed moment38/3 or24. Applying (3) to that same multiplet gives

    G1 reported original weight > 5/48,
    H  reported original weight > 11/48.                  (5)

The complete eight-row Walsh experiment has Kraus operators
A_(f_q,j)/sqrt(48) and their adjoints, plus eight retained address failures.
Its success effects, including the same clock, are exactly1/8 of the
selected-setting effects. Thus every positive spatial row has original
weight greater than1/24, with G1 greater than5/384 and H greater than11/384.
All96 successes, eight failures, clocks, native-number meters and reference
outputs are retained. The finite Walsh label is not conserved momentum.

## 3. Actual total-cube and native-number records

The source's five native numbers N_r commute with H, its Gauss projection,
and every physical cube transformation. The actual source group is the
binary cube group of48 proper spin lifts, together with the separate
[physical inversion](../particle_primary_pair_channels/PROOF.md#3-physical-rotations-and-inversion):
x maps to1-x, oriented links reverse, and the Dirac lift is beta. Write
S_i=(-beta+beta gamma_i)/2 for the Wilson hopping spin matrix and S(g)
for a proper source spin lift. Exact matrix identities beta S_i beta=S_i* and [beta,S(g)]=0 preserve the
Wilson action and commute with the proper lifts. The link-word cutoff
and physical graph are invariant. This inversion is not the central
2pi spin lift (-1)^N. Let Pi_alpha denote its16 central isotypic projectors, with
labels A1,A2,E,T1,T2,G1,G2,H and inversion parity g or u. In particular

    Pi_alpha=(d_alpha/96)sum_g conjugate(chi_alpha(g)) U(g).

Character orthogonality gives a complete orthogonal projective instrument
on the physical source. These projectors commute with H, rho, all source
energy-clock effects and the native numbers. They measure total source
cube sectors; they are not projectors for a continuous Lorentz spin.

Insert a complete Pi_alpha record before the probe and a complete Pi_beta
record after it. Native-number projectors may also be inserted, or the
existing coherent difference meters may record Delta without an initial
absolute-number readout. Summing every added label leaves the original
classical source-energy law unchanged. Before the probe this follows from
sum_alpha Pi_alpha rho Pi_alpha=rho (and the analogous identity for N_r);
after it, trace preservation of the complete final record gives the result.
The initial measurement can change source-reference coherences. They are
retained as part of the enlarged experiment, not asserted unchanged.

For the uniform Walsh row q=0 the positive creation multiplets are G1g
and Hg. The map v_j tensor psi -> F_(q=0,j)psi is an intertwiner from
lambda tensor alpha into the source representation. Hence

    Pi_beta F_(q=0,j) Pi_alpha = 0
         if beta is absent from lambda tensor alpha.      (6)

The exact actual-group character calculation gives28 distinct allowed
initial/final pairs for G1g and44 for Hg. Multiplicities are retained;
counting an allowed pair does not assign it positive weight. The native
labels obey n_final=n_initial+Delta on every creation event. The binary
central spin lift flips sign, while inversion parity is preserved by the
uniform positive-beta operator. Failures have their full separate records.

Equations (5)-(6) show that at least one allowed total-cube transition in
the complete uniform Walsh row has original reported weight greater than
5/10752 for G1 and1/1536 for H. Which transition attains that bound is not
computed. The receipt contains the full fusion tables and the original
probability bounds, rather than assigning a guessed final spin sector.

## 4. Explicit original-state total-cube branches

A separate original-state experiment identifies particular sectors with
a small authenticated weight. It records initial native numbers and total
cube irrep, applies the same uniform creation apparatus, and records final
native numbers, total cube irrep and source energy. There is no initial
energy measurement in this experiment.

The physical filled-negative-beta, constant-link vector Omega has
native occupations n0=(48,32,96,48,16), total cube label A1g and norm one.
Each actual proper spin lift has determinant one on the occupied beta
rows. Each site contains30 occupied modes, so site permutations introduce
no extra sign. Inversion contributes(-1)^240=1. These facts establish its
A1g label without assuming that Omega is an H eigenvector.

The normalized states phi_(lambda,j)=A_(uniform,j)Omega have native
occupations n0+Delta and total cube label G1g or Hg. The
[full-sector bound and uniform vectors](../particle_primary_composite_spin/PROOF.md#5-populated-complete-symmetry-and-energy-records)
give source energy floor L_Delta>=L=-18055527/25000 and mean M=-711.
A complete final clock with resolution1/64 and event loss at most1/512
reports [L-1/64,2M-L+1/64] on every such vector with probability
strictly greater than255/512. This interval is above the authenticated
full-source ground upper; the inferred excitation interval keeps both
ground-energy bounds explicit.

The original selected source satisfies rho>=c P_G with the authenticated
factored c=2^(-582407133). For the joint initial label projector Pi,
Pi Omega=Omega and Pi rho Pi>=c|Omega><Omega|. For one successful
uniform Walsh component K=A_(uniform,j)/sqrt(48),

    K Pi rho Pi K* >= (c/48)|phi_(lambda,j)><phi_(lambda,j)|. (7)

The final number and cube projectors fix phi. Applying the actual positive
final-clock event effect to (7) and summing its d_lambda components gives
explicit original branches

    (n0,A1g) -> (n0+Delta,G1g): weight > (85/4096)c,
    (n0,A1g) -> (n0+Delta,Hg):  weight > (85/2048)c.         (8)

The clock error is included inside the witness expectation before
multiplication by c. Subtracting an unscaled detector error from c times
an ideal probability would not prove (8). Every other initial label,
final label, failure and energy result is kept. Omega is a positive-effect
witness, not the source input or a replacement purification.

Literal source C maps these branches to n0-Delta with G1u or Hu,
respectively, with identical complete classical record probabilities.
It commutes with proper rotations and obeys C I C*=(-1)^N I, so odd
conjugate sectors have opposite inversion parity. The untouched reference
outputs need not be equal as quantum states.

An initial energy record may be added and summed without changing this
original final-energy marginal, since rho commutes with H. This does not
turn (8) into a transfer-window bound. Conditioning on an initial energy
dephases the witness in H's eigenspaces; Omega is not an interacting
eigenvector. The transfer windows in sections2-3 have their own original
two-energy derivation.

## 5. The same labels contain distinct source energies

Actual cube and native-current labels do not make the Hamiltonian a scalar
on the entire sector. Let e be the native charge-six index, and use its two
positive-beta rows p0,p1 and two negative-beta rows m0,m1. Define the local
neutral operator

    Q_y=c*(p0)c*(p1)c(m1)c(m0).

It adds four beta defects, preserves all five native numbers and is a
physical gauge singlet. Proper rotations act on its two exterior products
by determinants one, and inversion has two negative factors. Thus Q_y
is an even A1g scalar. Its action on the filled-beta anchor has unit norm.
The code checks its literal source modes and all48 proper spin lifts.

In addition to the uniform three-defect vectors phi_j, form

    psi_j=(1/sqrt(56))sum_(x!=y) F_(x,j) Q_y Omega.           (9)

The weak-pair site x and the separate four-defect site y distinguish every
ordered pair. The12 vectors phi_j,psi_j are exactly orthonormal on the
actual480-mode matter carrier, with all links constant. The code constructs
them as sparse native CAR vectors and checks every occupation label and
Gram entry. Both copies transform as G1g plus Hg with the same n0+Delta.

They have defect numbers3 and7. Constant-link Haar expectation removes
all hopping and magnetic terms, while electric energy is zero. Their
exact source energy means are consequently

    <phi_j,H phi_j>=-711,   <psi_j,H psi_j>=-699.           (10)

Each same-current G1g and Hg reducing sector therefore has at least two
multiplicity directions and spectral diameter at least12. Equation (10)
does not identify either mean with an eigenvalue. It does not exclude an
isolated low branch inside a sector. It proves that the conserved labels
and their nonzero readout weights alone cannot select one mass or replace
a full-complement spectral analysis.

## 6. Nonzero width of the actual original transfer measure

A positive norm supplies an original-state energy-shape consequence without
substituting a prepared signed moment. For one normalized component write
A_f=sum_x f_x F_(x,j), sum_x|f_x|²=1, and R_f=[H,A_f]-9A_f. On the
filled-beta physical anchor, a negative-beta creation in the commutator
annihilates the vector. The actual surviving positive-beta hopping block
is -R(U)/2 in both edge orientations. Unitarity of the native R(U) and
three distinct neighbors give the one-particle image Gram (3/4)I6.

Replacing one of the three creation factors by its external image has
canonical vacuum Gram3I6 across the six normalized internal components.
The code computes this exact exterior-CAR Gram, without taking a Fock
operator norm as a vacuum equality. It verifies the literal beta blocks.
Different origins have different matter profiles: two positive-beta
particles at the origin and one at a neighbor. They are orthogonal,
including for coherent spatial coefficients f. The created anchor has
three particles at a single origin and is orthogonal to every residual
profile. Therefore

    ||A_f Omega||²=1,
    ||R_f Omega||²=9/4,
    <A_f Omega,R_f Omega>=0.                                (11)

All terms have at most one native matrix-coordinate word on one link.
The actual degree-one per-link cutoff admits this entire first action.
F and [H,F] are gauge invariant, so the resulting vectors are physical;
Gauss restriction removes no part of (11). No closure of arbitrary
products of compressed link operators is used.

For any real center t, stationarity and the physical density floor give

    integral (omega-t)² dmu_(A_f)(omega)
      =Tr rho ([H,A_f]-tA_f)*([H,A_f]-tA_f)
      >=c ||([H,A_f]-tA_f)Omega||²
      =c((9-t)²+9/4).                                     (12)

The inequality acts on a positive source-marginal operator. Omega is a
witness for that inequality; its prepared energy law is not substituted
for the original measure. Since A_f is a contraction, the total measure
is at most one. At its actual mean, (12) implies normalized variance
at least(9/4)c. The apparatus scaling1/6 or1/48 cancels in this normalized
statement. The same bound holds for each component-summed G1/H measure
by summing (12) and its corresponding success cap.

Thus every actual original component transfer measure has at least two
distinct frequencies. The tiny finite density floor controls this result;
no line locations or a useful separation at the reported clock resolution
are computed. A measure with a pole plus other spectral weight satisfies
(12), so this consequence does not exclude an isolated particle pole.

The supported records are finite original-state source energies, actual
native-current increments and total cube sectors. Their relation to an
isolated relativistic particle, a conserved spatial momentum, a controlled
refinement limit and a measured GeV mass is not established by these bounds.
