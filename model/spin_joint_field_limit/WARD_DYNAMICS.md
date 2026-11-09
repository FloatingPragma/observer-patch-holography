# Source susceptibility and the two force equations

Use the positive-field ground and the same positive clock K in `PROOF.md`.
K can denote a finite source clock or its specified positive-field ground
representation. Canonical P is the uniform internal-spin field; it is not
the generator of spatial translations. All finite examples below use the
actual source and the joint-family value p=1/4, B=1/32, J=1.

## 1. Both source derivatives

For positively oriented bonds y=x+e_i put

    C_xy^a=(S_x cross S_y)^a,
    Q_p(F)=p^2 sum_x epsilon_x F(px) S_x^3,
    P_p(G)=-p sum_x G(px) S_x^2.

The exact source continuity equation gives

    i[L_B,p,P_p(G)]=-p^-2 Z_p(G)-lambda_p Q_p(G),
    Z_p(G)=p^3 sum_(x,i) D_i^p G(px) C_(x,i)^2,
    lambda_p=B/(Jp^2),                                  (1)

as in `FORCE_READOUT.md`. The other equation is

    i[L_B,p,Q_p(F)]=T_p(F)+(B/J)P_p(F),
    T_p(F)=-p sum_<xy> epsilon_x (F(px)+F(py)) C_xy^3.   (2)

Both identities hold before any approximation. The current operator has
bond norm 1/2. Replacing a spin in a current by its mean is not used.
On the joint family lambda_p=1/n and B/J=n^-5; the selection torques
vanish on the vacuum vectors by the proved uniform second-moment bounds.
Equations (1)-(2) retain the distinct exchange currents that must close.

The ground marginal imposes a strict constraint on a mean-spin
approximation to these currents. With `kappa_B` the actual transverse
bond coefficient and `m_B` the actual staggered order, the variational
trial and two-site positivity in [QUANTUM_FORCE_SLACK.md](QUANTUM_FORCE_SLACK.md)
give

    kappa_B-m_B^2 >= delta_n^2+delta_n/2,
    delta_n=1/120-1/(150n^5), n>=2.

Its limiting lower bound is 61/14400. Exact first-moment tests therefore
exclude vanishing unfiltered microscopic residuals for the replacements
`iLQ=12m_B P` and `iLP=-m_B Q(-Delta_p G)`. This constraint leaves
corrected kinetic coefficients and projected infrared subbranches
undetermined; it is not a no-shell conclusion.

## 2. An exact reduction to static source responses

Let q=Q Omega and a=P Omega for two Hermitian source probes with zero
ground mean. Assume their positive first and inverse energy moments are
finite. At finite positive B this follows from the gap and the finite
source. Write

    h_Q=||K^(1/2)q||^2,    chi_Q=||K^(-1/2)q||^2,
    h_P=||K^(1/2)a||^2,    chi_P=||K^(-1/2)a||^2,
    c=Im<q,a>.

For a vector not in Dom K the expressions below mean the form vectors
i K^(1/2)q-u K^(-1/2)a and i K^(1/2)a+b K^(-1/2)q. No unbounded
operator derivative is asserted merely by writing their norms.
If chi_Q,chi_P>0, set

    u=c/chi_P, b=c/chi_Q,
    delta_Q=h_Q-c^2/chi_P,
    delta_P=h_P-c^2/chi_Q.                              (3)

Expansion of the Hilbert-space norms proves

    delta_Q=||i K^(1/2)q-u K^(-1/2)a||^2>=0,
    delta_P=||i K^(1/2)a+b K^(-1/2)q||^2>=0.           (4)

Indeed Re<i K^(1/2)q,K^(-1/2)a>=c and
Re<i K^(1/2)a,K^(-1/2)q>=-c. The same expansion for any real coefficient
gives delta_Q+chi_P(u_trial-u)^2 or
delta_P+chi_Q(b_trial-b)^2. These are source-determined best coefficients
in the displayed form norm. They are not empirical velocity fits.
The argument does not require Gaussianity, a pure imaginary cross
correlation, a single spectral line or time-reversal symmetry.

For the unfiltered source probes Q_p(F),P_p(G), the canonical spin algebra
determines the cross moment exactly:

    c_p(F,G)=(p^3/2) sum_x epsilon_x F(px)G(px)<S_x^1>.
                                                               (5)

On the infinite positive-field ground this is
(m_B p^3/2) sum FG. Its limit is (m/2) integral FG.
The first moments also are static source expectations:

    h_A=(1/2)<[A,[K,A]]>.                              (6)

The two-bond identities in `PROOF.md` evaluate (6) as actual transverse
bond and selection-field scores. Thus (3)-(4) eliminate the need to
measure separate current-current and current-field correlators for this
form-norm test: the extra dynamical inputs are two ordinary one-field
inverse-energy susceptibilities. Their saturation is an additional source
statement, not a consequence of their proved upper and lower bounds.

## 3. Filtered equations and an original-weight spectral tube

Let R be a bounded spectral multiplier of K, of norm at most one, supported
in [ell,M], with 0<ell<=M. It may also include a compatible actual spatial
spectral filter commuting with K. Set q_R=Rq, a_R=Ra and replace every
quantity in (3) by its filtered value. In particular

    c_R=Im<Q Omega,R*R P Omega>                         (7)

is a separate original-weight cross moment. It cannot be replaced by (5).
The finite-window susceptibility is bounded and involves no inverse-gap
extrapolation. Filtered vectors belong to every power domain of K.
With r_Q=iKq_R-u a_R and r_P=iKa_R+b q_R, (4) implies

    ||r_Q||<=sqrt(M delta_Q),
    ||r_P||<=sqrt(M delta_P).                           (8)

For c_R!=0, omega=sqrt(ub)>0. Direct substitution, with no asymptotic
ansatz, gives

    (K^2-ub)q_R=-iK r_Q-u r_P,
    (K^2-ub)a_R= b r_Q-iK r_P.                          (9)

Hence, for any epsilon>0,

    ||1_(|K-omega|>=epsilon)q_R||^2
      <= M (M sqrt(delta_Q)+|u| sqrt(delta_P))^2
                                      /(omega^2 epsilon^2),
    ||1_(|K-omega|>=epsilon)a_R||^2
      <= M (|b| sqrt(delta_Q)+M sqrt(delta_P))^2
                                      /(omega^2 epsilon^2).     (10)

The scalar inequality used is |E^2-omega^2|>=omega|E-omega| for E>=0.
If the complete Q preparation is Q/Lambda_Q, the original joint
preparation/filter/tube probability is at least

    max(0,(||q_R||^2-right_side_Q_of_(10))/Lambda_Q^2).   (11)

No success probability has been divided out in (11). All other preparation,
filter and energy records remain in the instrument. A common nonvanishing
lower bound and a shrinking tube require a uniform source estimate on
the right of (10), not merely a finite positive denominator.
The direct Pauli-sum normalization Lambda can grow with the source support.
A nonzero unnormalized Stone-vector norm alone therefore does not supply
a nonvanishing original preparation probability in this protocol.
The case c_R=0 supplies no positive frequency by this construction.

These scalar-coefficient tests concern the specified pair of wavepackets.
A dispersion graph requires compatible momentum localization and control
of the coefficients across those packets. A positive full-response slack
cannot exclude a different spectral subbranch. Nor does a small unfiltered
slack justify exchanging source first moments with a projected continuum
limit: high-energy or orthogonal norm may be discarded in that limit.
Filtered (7)-(11) keep those possibilities distinct.

## 4. The spatial Laplacian test

For a real compact smooth G take a=P_p(G)Omega and
q_D=Q_p(-Delta G)Omega. Applying (3)-(4) to this pair gives

    c_D=(m_B p^3/2) sum_x (-Delta G)(px)G(px),
    b_D=c_D/chi_D,
    delta_D=h_P-c_D^2/chi_D.                            (12)

The limit of c_D is (m/2) integral |grad G|^2. The proposed P-field
equation is i[L,P(G)]=a Q(Delta G), or iL a+a q_D=0 on the ground.
The exact form-norm defect for a specified coefficient a is

    delta_D+chi_D(a-b_D)^2.                             (13)

With a fixed compact energy filter, multiply the square root of (13) by
sqrt(M) to bound the filtered ordinary norm. This is a sufficient source
route to the spatial-force equation if the source expressions vanish
and one common coefficient a works for the admitted tests. It does not
replace Q(-Delta G) by a multiple of Q(G) unless that additional relation
has been established for the chosen spatial mode. The unfiltered Ward
moment (12) and filtered cross moment (7) retain their different scopes.

## 5. An exact finite-field Ward specialization

On any finite bipartite source graph let

    A_0=-sum_x S_x^2,  M_3=sum_x epsilon_x S_x^3.

The exchange Hamiltonian commutes with the uniform total spin, so

    i[H_B,A_0]=-B M_3,
    K A_0 Omega=iB M_3 Omega.                            (14)

The unique positive-field ground has zero transverse means. Equations
(14) and the source commutator give the exact susceptibility

    <M_3 Omega,K^-1 M_3 Omega>
                 =(sum_x epsilon_x<S_x^1>)/(2B).        (15)

For any source energy multiplier R, (14) persists after R. It fixes the
entire relation between the two filtered vectors rather than only their
first moments. It does not prove that their energy distribution is sharp:
the first-force slack in (4) measures the independent obstruction.
The mode is uniform on the finite graph. It is not a nonzero-momentum
continuum particle or a replacement for spatially localized preparation.

The remaining obstruction for this uniform pair is an explicit weighted
energy variance. Let q=R M_3 Omega be nonzero, a=R A_0 Omega, where R
commutes with K and is supported in a positive energy window. For the
original measure nu(dE)=||E_K(dE)q||^2 put M_j=int E^j nu(dE).
Equation (14) gives a=iB K^-1 q, hence

    c=B M_-1, chi_P=B^2 M_-3, h_P=B^2 M_-1,
    delta_P=0,
    delta_Q=M_1-M_-1^2/M_-3
           =M_-3 Var_pi(E^2),
    pi(dE)=E^-3 nu(dE)/M_-3,  u b=M_-1/M_-3.          (15a)

All moments exist on the positive window. The same identities hold on a
finite positive-gap source without the window. The auxiliary probability
pi proves the variance identity; it does not renormalize any preparation
or retained detector outcome. Since E>0, delta_Q=0 is equivalent to this
specific response measure being supported at one positive energy.
The exact susceptibility (15) fixes M_-1 only. It supplies neither M_-3
nor the determinant M_1 M_-3-M_-1^2. For a spectral or charge filter the
full susceptibility (15) cannot be substituted for its filtered value.

For example the abstract measure (delta_1+delta_2)/2, with
a=iB K^-1 q, obeys the exact first Ward equation but has
M_-1=3/4, M_-3=9/16, delta_Q=1/2 and Var_pi(E^2)=8/9.
This is an implication control, not a claim about the native continuum
measure. The actual dimer instead has zero weighted variance, checked
from its complete source spectrum. Higher-volume charge blocks need not
be one dimensional, so their conserved charge alone supplies no such
vanishing estimate. A useful continuum estimate must retain original
response weight and use the same scaled clock: replacing K by K/(Jp)
divides the weighted mean of E^2 by (Jp)^2 and its variance by (Jp)^4.
Collapse of an unscaled energy to zero is insufficient for a sharp
nonzero scaled response.

## 6. Reading the source quantities without conditional renormalization

At finite volume Q and P are known Pauli sums. A declared coherent source
service supplies complete block-encoded preparations Q/Lambda_Q,
P/Lambda_P and (Q+-iP)/Lambda_mix, with
Lambda_mix>=Lambda_Q+Lambda_P. Every preparation failure is retained.
Use the same complete energy readout after these branches. Its scores

    |R(E)|^2, E|R(E)|^2, |R(E)|^2/E

give norms, first moments and inverse moments, multiplied by their declared
Lambda^2. The last score is zero outside the positive window and bounded
by 1/ell; no ground-energy division is performed on an unretained event.
For an unfiltered finite source the proved gap supplies that bound.
The two interference preparations determine

    c_R=(||R(Q-iP)Omega||^2-||R(Q+iP)Omega||^2)/4.       (16)

They can be settings of one complete instrument; the setting probabilities
must then be included in the scores. A score f of absolute bound F on a
branch normalized by Lambda has absolute moment error at most
2 Lambda^2 F eta for complete-output half-trace error eta, plus the
certified score/energy-discretization error. The factors in (16) and any
fixed setting probabilities are explicit additional multipliers.
Certified positive susceptibility intervals permit interval inversion in
(3); a merely positive numerical central value does not.
The third inverse moment in (15a) is read by the additional complete
score |R(E)|^2/E^3, bounded by ell^-3 on the same retained window.
Its error is therefore charged to the same original preparation and
setting probability by the complete-score bound, with score bound ell^-3. The kernel
and every unretained outcome have score zero; no division at zero energy
or implicit removal of preparation failure is used.

The existing source inverse-response filter and fixed-scale positive-field
preparation/clock interfaces supply the finite ingredients under their
stated coherent-access and precision premises. This construction specifies
complete finite effects and scores. It does not execute an additional
large native pulse tape or give a continuum convergence rate.

## 7. Exact controls on the complete joint-family dimer

The actual two-site source at n=2 has the full spectrum and middle-level
projector in `FORCE_READOUT.md`. Its transverse Q/P preparations have all
their energy at

    E=(16+sqrt(257))/8

on the scaled clock. With profile F=G=(1,1), the two ground-vector
force slacks vanish exactly, and the source coefficients are

    u=sqrt(257)+513/32, b=1/2, ub=E^2.                  (17)

The exact susceptibility (15) and both current signs are checked against
all four source rows. These equations concern the specified source-created
ground vectors; they do not close the operator field algebra on every
history or define a continuum particle.

With F=G=(1,2), the same complete middle-energy line carries all Q/P
weight, but the two positive form slacks are

    delta_Q=81/4658+20781 sqrt(257)/19153696,
    delta_P=72/2329+2349 sqrt(257)/1197106.              (18)

Both exceed 1/100. The source quadratures have different polarizations
inside the degenerate level. Thus a positive whole-pair force slack does
not even imply absence of a sharp energy line in this finite source.
It cannot be used to exclude an independent continuum subbranch.

The actual conserved internal source charge in the physical spin frame is
C=(X_0+X_1)/2. Its projectors

    R_+=(I+X_0)(I+X_1)/4,
    R_-=(I-X_0)(I-X_1)/4,
    R_dark=I-R_+-R_-

commute with this same H_B. Filtering the profile (1,2) by either R_+ or
R_- gives both force slacks exactly zero, with ub=E^2. Each filtered Q
and P creation has original probability greater than 1/25 under its
declared Pauli-sum normalization. The two cross moments have opposite
signs and sum to the full c. Replacing either filtered c_R by that full
c would give the wrong source coefficients. This is an actual finite
charge/polarization resolution with complete retained records, not an
assumed continuum projection or physical species assignment.

For each profile the four settings Q,P,Q-iP,Q+iP each retain preparation
failure and all three source energy projectors, including zero-probability
labels. The resulting sixteen original effects sum to identity and their
probabilities sum to one. The interference scores reproduce (5) exactly.
No conditioning, empirical fit, omitted source level or added native pulse
execution is used in these controls.
For the profile (1,2), the additional charge-resolved instrument replaces
the three energy records by R_-, R_+ and R_dark in each of the same four
settings, retaining sixteen original records including preparation failure.
