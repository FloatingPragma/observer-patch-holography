# A bounded test for closure of a source response space

## 1. Exact reducing defect on the actual positive clock

Let K>=0 be the chosen source ground generator and V:H_r->H an isometry
onto a closed response subspace. V must be constructed from the chosen
source preparations; it is not an assumed particle embedding. For z>0 put

    S=z(K+z)^(-1),  0<S<=I,
    A=V* S V,  D=V* S² V-A²,
    R=(I-VV*) S V.                                         (1)

All operators in (1) are bounded, even when K is unbounded. Direct
multiplication proves

    D=R*R>=0,  ||R||=sqrt(||D||).                           (2)

Consequently D=0 if and only if the response space reduces S. Indeed
R=0 gives invariance, and selfadjointness gives invariance of the
orthogonal complement. The positive resolvent is an injective function
of K; its spectral projections generate those of K. Thus the same
condition is equivalent to the response space reducing K, with its
spectral domain. This criterion is exact for one fixed z>0. It assumes
no first or second energy moment of the prepared vectors.

If D=0 and the compressed resolvent is independently identified as

    A=z(omega(P)+z)^(-1),                                  (3)

for the actual spatial momentum on the response space, then the source
generator restricted to that space is omega(P). This use of momentum
requires that the response space also reduces the actual spatial
translation representation; a compressed momentum is insufficient. Neither (2) nor positive
energy proves (3). A reducing response space can contain several
spectral branches. A nonzero D for one preparation likewise does not
exclude a different reducing subspace or a sharp spectral subbranch.

## 2. A finite error estimate and its energy limitation

For a polynomial f(u)=sum_(j=0)^N c_j u^j, telescoping gives

    S^j V-V A^j=sum_(k=0)^(j-1) S^(j-1-k) R A^k,
    ||f(S)V-Vf(A)||<=sqrt(||D||) sum_(j=1)^N j|c_j|.         (4)

The same bound applies to every vector, with its norm. If a continuous
function g on [0,1] is uniformly approximated by f within eta, (4) gives

    ||g(S)V-Vg(A)||<=2eta+sqrt(||D||) sum j|c_j|.            (5)

This is a useful bounded-observable comparison. It is not an unrestricted
time-evolution bound: exp[-it z(1/u-1)] is not continuous at u=0 for
nonzero t. Small D at one z can hide leakage to arbitrarily high energy.
To compare dynamics with (5), retain a finite energy window or a continuous
cutoff vanishing near u=0 and control the omitted original probabilities.
No conditioning removes those tails.

A concrete control makes this qualification necessary. On two levels let
K_M=diag(M,3M) and V=(1,1)/sqrt(2). Then D_M tends to zero as M tends to
infinity, while at t=pi/(2M) the source unitary sends V into its orthogonal
complement. Along odd integer M, the same escape occurs at the fixed
time t=pi/2. Thus even a fixed-time conclusion needs energy control; a
small bounded-resolvent defect alone is insufficient.

## 3. What a complete source readout must supply

For finitely many actual prepared vectors psi_j, form their Gram matrix
G_ij=<psi_i,psi_j>. If G is strictly positive, set V=Psi G^(-1/2). Then

    A=G^(-1/2) M1 G^(-1/2),
    D=G^(-1/2) M2 G^(-1/2)-A²,
    (M_l)_ij=<psi_i,S^l psi_j>,  l=1,2.                   (6)

Equation(6) identifies the required data. It does not implement the
abstract projection VV* as a local detector. Off-diagonal entries require
coherent access to the specified preparation branches and compatible
same-clock bounded filters or interference measurements. A positive
lower bound for G is required to control inversion errors. A finite
nondegenerate Gram theorem does not supply such a uniform bound for a
growing family. The source preparation probabilities and all filter
complements remain part of the measurement.

For a single source vector of original norm squared w>0, complete energy
records with bounded scores S(E) and S(E)^2 determine

    D=b/w-(a/w)^2,
    a=<psi,S psi>,  b=<psi,S² psi>.                       (7)

If w is exact and score errors are at most epsilon_1 and epsilon_2,
then a bounded estimate clipped to [0,w] gives

    |D_est-D|<=epsilon_2/w+2epsilon_1/w.                    (8)

This error bound explicitly retains the cost of a rare preparation.
If w is also uncertain, first certify a positive interval and propagate
(7) with that interval; no source-independent denominator is asserted.
Existing complete source energy readouts and resolvent filters supply
such scores only within their declared finite preparation, clock,
control-error and access domains. An arbitrary continuum isometry V is
not declared operational by (6).

## 4. Exact native dimer controls

The producer uses the unchanged two-site `native_ordered_ground` source
at J=1,B=3/8. In its Marshall basis the unique ground is

    Omega=(1,3,3,1)/sqrt(20),  E0=-9/8,
    spectrum(H_B)={-9/8,0,0,1/8}.

These values are checked against the actual source matrix. K=H_B+9I/8
is the positive clock, and z=1/8 is shared by all readouts below.
The local Z0-created vector has a single positive energy9/8 and hence
D=0. No mass is assigned to that finite excitation.

For a nonreducing response, prepare the complete local instrument

    A=(X0+Z0)/2,  A_fail=I/sqrt(2),

then apply F=K/(K+z)=I-S to the A branch, retaining its complementary
filter outcome. The selected original vector is psi=F A Omega. There
are three final records: preparation failure, filter failure and success.
Their effects sum to I. The complete positive-energy success weights are

    E=9/8: 81/400,     E=5/4: 16/121,
    w=16201/48400.                                         (9)

The source energy detector retains both levels and every complementary
record. Equation(7) gives a strictly positive exact reducing defect. The
span of these two source spectral components does reduce K, so its
matrix-valued D is exactly zero. This demonstrates both a genuine
nonreducing prepared response and a different reducing response space
on the same source, state and clock. The small operator calculation
specifies a complete instrument and inherited finite readout inputs;
it executes no additional native pulse tape and supplies no continuum
implementation theorem for a Gram projection.

## 5. Why the extra CCR pair and internal circle do not close the force

The two-speed local-field control in `asymptotic.two_speed_control()`
satisfies the stated single-pair CCR, positive energy, source
first-moment values and a finite cone. Take two independent copies and
rotate their real components together. The observed variables have two
regular CCR pairs, and this internal U1 preserves the Hamiltonian and
its vacuum. Scalar pure-direction twist Hessians are unchanged.

In visible and hidden coordinates per real component, the same control is

    H=(1/2) integral [12(pi_q²+pi_x²)+|grad q|²+|grad x|²
                       -grad q.grad x],
    Q=sqrt(m)q,  P=sqrt(m)pi_q,  m=kappa=1/4.               (10)

Its normal-frequency squares are6|k|² and18|k|². The visible source-form
quadratic kick costs are6kappa integral F² and
(kappa/2) integral |grad G|², exactly as in the inherited control.
Nevertheless

    d_t² Q=12 Delta Q-6 Delta X.                            (11)

The hidden force has nonzero vacuum response. No single wave equation
for Q follows from the two observed canonical pairs or their circle
action. This is an implication countercontrol, not a claimed state or
limit of the spin source. Its dynamics is free on the larger field system;
nonzero projected memory does not prove interacting scattering.

The coefficient choice in (10) is an exact illustration, not a constraint
on the source stiffness. For any m,kappa>0, use

    a=12 kappa/m,  b1=kappa/(2m),  b2=3kappa/(2m),
    Q=sqrt(m/2)(q1+q2),  P=sqrt(m/2)(pi1+pi2).             (12)

The two original fields are independent, with kinetic coefficient a and
gradient coefficients b1,b2. The canonical coefficient is m and the exact
vacuum energy forms are ma/2=6kappa and m(b1+b2)/4=kappa/2. Their squared
speeds are6(kappa/m)^2 and18(kappa/m)^2. Doubling the real components gives
the full internal circle as above. In visible/hidden coordinates the force
is12(kappa/m)^2 Delta Q-6(kappa/m)^2 Delta X, so no positive choice of the
source coefficients removes the hidden branch in this control.

The equal-time Q and P norm coefficients are
m(sqrt(6)+sqrt(2))/2 and m(sqrt(6)+3sqrt(2))/48. Their inverse-energy
coefficients are2m^2/(3kappa) and m^2/(24kappa). For the stated joint-source
range1/10<m<=1/2 and m^2/2<=kappa<=1/2, these satisfy the inherited upper
bounds3,1/4,4 and1/12, respectively. The Euclidean maximum speed
sqrt(18)kappa/m fits inside the source l1 cone of speed72, since its l1
speed squared is54(kappa/m)^2<1350<72^2. Thus actual stiffness values
within this range cannot select a single branch from the listed moments,
CCR, cone and circle alone. This family matches those constraints only;
it is not a construction of the actual source's higher correlations.

The finite `Dynamics.PrivateInner` theorem constructs the unitary
implementation after an automorphism of the specified full matrix block
is supplied. It does not prove invariance of a chosen continuum response
algebra. The mode-energy theorems in `FieldSectorEnergyInnerProduct`
consume a supplied mode flow and cannot select that flow from canonical
commutators alone. Applying them here requires the missing source
operator/dynamics identification.

A source proof of D=0, together with (3), would close the corresponding
response-generator identification. A full field equation additionally
needs the same action on the relevant field/history domain. Neither the
finite controls nor the scalar bond-twist coefficients assert these
vanishing conditions. The explicit source force estimates remain the
way to decide them.
