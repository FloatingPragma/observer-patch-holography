# Composite scalar radial records on the original source

## 1. One actual scalar contraction

Use the unchanged [primary source contract](../particle_common_contract/PROOF.md),
with all native fermions, full Gauss support, product link cutoff, selected
polynomial density rho, retained canonical reference and Hamiltonian H.
The spacing is a and the number of sites is V. At each site use the six
negative-beta modes of the native weak doublet and charge-six singlet:

    Phi_a=(l_a0 e_1-l_a1 e_0)/sqrt(2),
    S=sum_a Phi_a* Phi_a,       A=S/2.                         (1)

The [charged scalar calculation](../particle_candidate_charged/PROOF.md)
constructs Phi on all64 occupations, verifies all twelve native currents
and the actual48 proper spin lifts, and supplies its Higgs-conjugate
comparison weights. S is its gauge-invariant radial contraction. It is
an even one-site source operator, preserves all five native occupations,
has scalar proper-cube coefficients and even inversion parity. Spatial
rotations move its site. No intrinsic Lorentz scalar or total-spin
assignment to an arbitrary output state is inferred.

The exact full local spectrum and multiplicities are

    S: 0(29), 1/2(16), 1(10), 3/2(8), 2(1).                 (2)

The sole eigenvector at2 fills all six modes. Removing one weak mode
gives eigenvalue3/2, while removing one singlet mode gives eigenvalue1.
Let h=6-N_L-N_e be their negative-beta hole count. It commutes with S.
Diagonalizing every h block proves the sharp inequalities

    0 <= 2I-S <= h,
    0 <= I-A² <= (3/4)h.                                    (3)

Both constants are attained by a one-singlet-hole state. The code checks
every occupation block and recomputes the native-current commutators of
S. No filled-state restriction enters (2)-(3).

Phi is charged, so its one-point function vanishes in the invariant
source state. The expectation of S is a different observable. A positive
radial contraction is not a selected Higgs vacuum, and a zero charged
one-point function does not exclude a gauge-invariant Higgs phase.

## 2. Complete records and a capacity-uniform total probability

The complete radial detector has V successful rows and one failure:

    K_x=A_x/sqrt(V),
    K_fail=sqrt(I-sum_x K_x* K_x).                            (4)

Equation(2) gives 0<=A_x<=I. Hence (4) is a complete instrument on the
full physical source. All its rows preserve native occupation, total
fermion parity and Gauss support. The source, address, failure and
reference outputs are retained. Its Stinespring isometry can be extended
within the same finite conserved blocks, with a neutral apparatus register.

Let delta be the authenticated original weak-plus-singlet beta-defect
upper bound from the [selected odd source estimate](../particle_primary_selected_odd_response/PROOF.md).
Actual source charge conjugation exchanges positive-beta occupation and
negative-beta holes in each native block and preserves the selected
marginal. Thus sum_x Tr(rho h_x)<=delta/2. Equation(3) gives

    (1/V)sum_x <A_x> >= b:=1-delta/(4V),
    (1/V)sum_x <A_x²> >= 1-3delta/(8V).                     (5)

Positivity of the density and Jensen's inequality also give

    P_radial=(1/V)sum_x <A_x²> >= b².                       (6)

The source estimate proves b>=59/400 for every selected member. Therefore

    P_radial >= 3481/160000.                                (7)

At member zero the sharper effect inequality in(5) gives P_radial>383/500,
with the receipt retaining the stronger rational value. This is the
unnormalized probability of the complete experiment, with no selection
of sites or postselection on a prepared state. The lower bound(7) is
uniform in capacity, although the number of detector outcomes grows.

Complete initial and final source-energy clocks may accompany(4).
The initial unread clock fixes rho because its Kraus operators are
functions of H; the complete final clock is trace preserving. Summing
their records preserves (5)-(7). The initial clock can change correlations
with the reference, whose outputs are retained. The probability bounds
do not require a joint pure-state operator inequality.

These success probabilities include elastic response and any identity
component of A. They do not prove a uniform centered variance, a
positive-energy window, a particle pole or continuum spectral survival.

## 3. Exact radial dynamics in the full interacting source

Let Omega fill every negative-beta mode and take constant links. It is
a normalized physical vector, used only to evaluate source operators.
At the corner x=0 it obeys S_x Omega=2 Omega. Onsite, electric and
magnetic terms commute with S_x. Other hopping bonds have disjoint
support. Only the three actual corner Wilson bonds contribute to
[H,S_x]Omega.

The native Wilson matrix has beta-changing entries of magnitude1/(2a).
A hop creating a positive particle at the neighbor and a negative hole
at the corner changes S_x by nu_r, where

    nu_weak=1/2,       nu_singlet=1,       nu_other=0.        (8)

Hops with a negative hole at the neighboring site leave S_x unchanged.
For each native representation r of dimension d_r, the summed squared
norm of the former direction is d_r/(2a²). This follows from the two
actual beta spin rows and Haar orthogonality
integral |R_r(U)_ij|²=1/d_r. All five native blocks and all three Wilson
matrices are included in the executed120-mode CAR calculation.

Different links have orthogonal nonconstant words; different native
representations and matter occupations are orthogonal too. Each first
word is admitted by the actual cutoff. Since S_x commutes with that
projection and with Gauss, these are full-source identities:

    ||[H,S_x]Omega||²
       =(3/a²)[(2/2)(1/2)²+(1/2)1²]=9/(4a²),
    ||[S_x,[H,S_x]]Omega||²
       =(3/a²)[(2/2)(1/2)^4+(1/2)1^4]=27/(16a²).          (9)

The code computes the CAR signs and every admitted Haar coefficient,
rather than treating a trial compression as invariant. Its independent
one-representation finite controls include both weak and singlet holes.
For A=S_x/2, the corresponding coefficients are9/16 and27/256.

Let mu_A be the ideal two-energy transfer measure of the single-site
success row A on the original rho. The authenticated physical source
floor rho>=c P_G, c=2^(-N), and the exact selected polynomial's positive
logarithmic slope beta_low give

    integral omega² dmu_A >= c*9/(16a²),
    integral omega dmu_A
       >= (c beta_low/2)||[H,A]||_HS²
       >= c beta_low*9/(32a²).                             (10)

The second line is the stationary spectral double-commutator identity.
The first uses Tr rho [H,A]*[H,A]; Omega is a norm witness, never the
selected input. In an energy basis, p_i-p_j>=c beta_low(E_j-E_i)
for E_j>E_i. Hermiticity of A and monotonicity of p also imply that
its negative-transfer second moment is at most its positive one.
If ||H||<=R, the original ideal positive-transfer weight therefore obeys

    mu_A((0,2R]) >= c*max(9 beta_low/(64a² R),
                          9/(128a² R²)).                  (11)

This is a strictly nonelastic source spectrum. The bound c deteriorates
with capacity; (11) is not a resolved finite clock bin or a uniform
inelastic residue. The comparison between its small floor and (7) is
not omitted: the large total radial detection can be elastic.

## 4. One complete finite-time scalar response

A finite two-record experiment resolves noncommutation using the same
source clock without allocating a denominator of size2^N. Put

    tau=1/(64aR²),
    Q_tau=(I+i tau H/2)^(-1)(I-i tau H/2),
    B_tau=Q_tau* A Q_tau,
    L_tau=[A,B_tau]/2,     L_fail=sqrt(I-L_tau*L_tau).       (12)

Q_tau is the specified rational Cayley approximation to the unchanged H
clock, with the error in(14); its logarithm does not replace H. Both A
and B_tau are contractions, so ||L_tau||<=1. Equation(12) is a complete
physical source instrument; its source operations, quantum output,
failure and reference are retained. It is an intervention with a
specified delay and detector, not a replacement of rho by Omega. At
tau=0 the success row vanishes identically.

This two-record instrument has an explicit finite unitary dilation.
L_tau is anti-Hermitian and commutes with F=L_fail. Thus the block
matrix W=[[L_tau,F],[F,L_tau]] is unitary. Its first ancilla column gives
exactly the two rows in(12), including the retained failure. The matrix
functions preserve the same finite conserved blocks and act trivially
on the idle reference.

For ideal evolution, A(t)=exp(itH)A exp(-itH) satisfies

    ||[A,A(t)]/2-it[A,[H,A]]/2|| <= 2R²t².                (13)

The exact witness(9) gives leading norm sqrt(27)t/(32a)>t/(8a).
Functional calculus bounds the Cayley error by

    ||Q_t-exp(-itH)||<=t³R³/12,
    ||[A,Q_t* A Q_t]/2-[A,A(t)]/2||<=t³R³/6.             (14)

The actual source enclosure has aR>=1. At t=tau, (13) is at most
tau/(32a), and (14) is at most tau/(64a). Thus

    ||L_tau Omega|| > tau/(16a),
    Pr_original(success at tau) > c tau²/(256a²)>0.         (15)

This positive finite-record law is contrasted with the exact zero at
delay zero. The receipt binds the actual rational duration and stores
the probability floor in exact factored form. The source-reference
purification is retained; only its source marginal enters the trace
inequality. Equation(15) is a two-time scalar response, not a measured
energy bin, Higgs mass, or scattering amplitude. Its tiny floor supplies
no practical finite-shot guarantee or continuum survival theorem.

## 5. Reuse and the quantitative boundary

The scalar and weak-current representation premises of the finite
[electroweak composition](../../Lean/Screen/ElectroweakBreakingComposition.lean)
have an actual Higgs-conjugate operator witness in(1). That theorem's
declared vacuum, kinetic normalization, quartic potential and couplings
are not produced by the local contraction. The
[charged-current dictionary](../../Lean/Screen/ChargedCurrentDictionary.lean)
retains those same attachment requirements. The
[scalar regional theorem](../../Lean/QFT/ScalarRegionalTimeSlice.lean)
proves velocity reconstruction for an assumed canonical scalar split;
(1) supplies no canonical field/momentum pair, so its force-readout
premise is not claimed instantiated. The supplied Higgs curvature and
interacting-Higgs controls have distinct source actions and states.
Their mass and vacuum results are not transferred to this source.

One direct attempt to replace the faithfulness factor in(10) by a large
selected-energy bound is quantified. Delete only the weak and charge-six
hopping from the member-zero source. The actual all-occupation electric
Schur estimate on the three retained native representations gives

    inf H_deleted >= -724.469,
    <H>_rho <= -723.7053,
    <H_hop,weak+singlet> <= 0.7634.                        (16)

Exact outward rationals are in the receipt. The positive upper in(16)
does not force a strictly negative restricted hopping expectation and
therefore cannot certify a macroscopic radial susceptibility. It is a
limitation of this estimate, not an exclusion of scalar dynamics, which
(9)-(15) prove directly. A connected spectral law with controlled
positive-energy weight and its common field-vacuum attachment are the
remaining physical requirements. No GeV value or named particle is
assigned by this package.

The native-number obstruction for a single sharp charged scalar does
not prohibit correlations between neutral radial contractions of
different bilinears: every such contraction has zero native increment.
That selection rule permits a common radial-channel calculation; it
does not establish its matrix elements, vacuum or spectral poles.
