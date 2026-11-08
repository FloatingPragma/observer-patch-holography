# A prepared joint positive-field family

## 1. The state and clock

The physical source is the spin-one-half cubic interaction

    H_B = J sum_<xy> (S_x.S_y-I/4) - B sum_x epsilon_x S_x^1,
    epsilon_x=(-1)^(x1+x2+x3),  J>0.

For each positive B, the uniform finite gap B/8 and the interior marginal
modulus in `spin_fixed_field_locality/PROOF.md` give a unique infinite-volume
ground state obtained from open boxes with these couplings. Write it omega_B,
with positive ground generator K_B. Passing the finite variance-gap inequality
on local vectors gives

    Var_B(A) <= omega_B(A* K_B A)/(B/8).                    (1)

Here A is local; its centered vector is understood on the right. The ground
vector subtraction has zero energy. This also proves the corresponding
spectral gap in this ground representation by density. The finite gap is an
explicit substitution into Rayudu--Takahashi, Definition 1 and Theorem 11,
[arXiv:2607.10765v1](https://arxiv.org/abs/2607.10765v1), as proved in
`spin_twist_closure/FINITE_FIELD_GAP.md`. No zero-field gap is asserted.

Define the rational sequence

    p=n^-2, B=J n^-5, n>=2,  L_n=K_B/(Jp).               (2)

Every operation at this scale uses omega_B and H_B time. Equation (2) defines
a chosen joint preparation and scaling prescription within the source
Hamiltonian family; the unread repair law does not select this power of n.
Equality of its full time histories to the construction
that sends B to zero before p is a separate statement, not used here.

The source has beta translations (a lattice shift followed, for odd parity,
by the global pi rotation about axis 2), rotations about internal axis 1,
and cubic spatial symmetries. The positive-field ground inherits them:
finite-box uniqueness and the common local limit identify the transformed
state. In particular the transverse means vanish and
omega_B(epsilon_x S_x^1)=m_B independently of x. The same selected-state
maximum-order theorem used in `spin_canonical_field/README.md`, Koma,
Theorem 2.1 of [arXiv:1712.09018v4](https://arxiv.org/abs/1712.09018v4), gives
m_B tending to the maximal selected order m as B decreases to zero.
The exact source variational certificate in `spin_charge_dispersion` gives
m>1/10. No numerical value or convergence rate for m_B-m is assumed.

## 2. Quantitative order and bond concentration

For real f supported in [-R,R]^3, R>=1, let N0=(2R+1)^3 and

    M_p(f)=p^3 sum_x epsilon_x f(px) S_x^1.

There are at most N0 p^-3 support sites and at most 6N0 p^-3 incident
bonds. For a single-axis A=sum a_x S_x^a, exact two-site algebra gives

    [A,[H_J,A]]=-J sum_<xy> (a_x-a_y)^2 T_xy^a,
    T_xy^a=sum_(b!=a) S_x^b S_y^b, ||T_xy^a||=1/2.

The selection term commutes with M_p. Since
|a_x-a_y|<=2p^3||f||_infinity, groundness and (1) imply

    q_B(M_p)<=6J N0 ||f||_infinity^2 p^3,
    Var_B(M_p)<=48 N0 ||f||_infinity^2/n.                  (3)

The exact mean is m_B p^3 sum f(px). Thus the limiting phase is the scalar
m integral f, with a quantitative fluctuation bound and a separate
unquantified convergence of its mean.

For a real bond sum Z_p=sum_b z_b T_b, assume
|z_b|<=a p^3, sum_b|z_b|<=a V0, and ||T_b||<=1/2. The local overlap count
from `spin_twist_closure/PROOF.md` gives the exchange double commutator
bound 198J ||z||_infinity ||z||_1. The selection term contributes at most
6B ||z||_infinity ||z||_1: each inner bond has two sites, each site belongs
to at most six outer bonds, and each nested onsite term has norm at most
4(1/2)(B/2)(1/2)=B/2. Consequently

    q_B(Z_p)<=(99J+3B)a^2 V0 p^3,
    Var_B(Z_p)<=(792/n+24/n^6)a^2 V0 <=816a^2 V0/n.       (4)

This centers on the actual source expectation, not on a fitted coefficient.
Weighted bond means have bounded subsequences. On a common such subsequence
the bounded bond phases become scalar. No numerical stiffness follows from
their concentration alone.

To pass (3)-(4) through p-dependent time histories, use the energy-form
transport lemma in `spin_canonical_field/README.md`. A phase exp(isZ_p)
and its inverse change energy square roots by at most
delta=|s|sqrt((99J+3B)a^2V0p^3). On a normalized vector with energy<=EJp,
at physical time T/(Jp), its squared time-commutation error is at most

    4|T| delta sqrt(EJp)/(Jp)+10|T|delta^2/(Jp)=O(p).     (5)

Equal-time commutators with the Q and P twists below are O(p^2) and O(p),
using their single-site angles and the bounded l1 coefficient sum. Move
the phase through each fixed finite word, then apply its vacuum
concentration. This proves scalarity on the complete bounded-history
representation, without assuming state ergodicity. Equation (5) uses the
actual H_B generator, including the selection term. It is not a comparison
between H_B and H_0.

## 3. Transverse fields and their energy

Use the same source normalization

    Q_p(F)=p^2 sum_x epsilon_x F(px) S_x^3,
    P_p(G)=-p sum_x G(px) S_x^2.                          (6)

For transverse A=sum a_x S_x^a, a=2,3, the additional double commutator is
B sum epsilon_x a_x^2 S_x^1. Its norm is at most (B/2)sum a_x^2.
The uniform energy-form lemma therefore gives the following squared
increments, valid on every finite-energy vector and its time translates:

    cost_Q = (Jp^4/4)sum_edges(F_x+F_y)^2
               +(Bp^4/4)sum_x F_x^2,
    cost_P = (Jp^2/4)sum_edges(G_x-G_y)^2
               +(Bp^2/4)sum_x G_x^2.                   (7)

After division by Jp, the extra costs are (B/J)(p^3 sum F^2)/4 and
(B/(Jp^2))(p^3 sum G^2)/4. On (2) these decrease as n^-5 and n^-1.
Every fixed finite time word has energy O(Jp), by adding the square roots
of (7). The same bounds hold for the two transverse axes interchanged.

For completeness, the finite-positive-field reflection-positive bounds
consumed in `spin_goldstone_sector/README.md` are uniform before the limit.
For a unit-normalized transverse mode, with e(k)=sum_i(1-cos k_i), they give

    integral E^-1 dmu_(B,k) <=1/(2J e(k+Q)),
    integral E dmu_(B,k) <=J e(k)+B/2,
    mu_(B,k)(R) <=sqrt((e(k)+B/(2J))/(2e(k+Q))).          (8)

The last inequality is Cauchy--Schwarz. A possible transverse Bragg atom
vanishes at each fixed B: a volume-normalized transverse average has
energy O(J/volume)+O(B/volume), and (1) makes its variance tend to zero.
Thus the density bound passes through its exceptional point without an
atomic contribution. The same finite-field infrared assumptions and
periodic-box limits as the cited source theorem are consumed here; the
positive-field marginal uniqueness identifies their local limit with
omega_B.

For |k|<=1, e(k)<=|k|^2/2, e(k+Q)>=5, e(k)>=|k|^2/3.
In the uniform channel the rescaled density is bounded by
sqrt((|q|^2+B/(Jp^2))/20). In the staggered channel the added B term tends
to zero against the integrable 1/|q| envelope. Since B/(Jp^2)=1/n,
Poisson summation and the Schwartz bounds used in the canonical-field
proof give

    limsup omega_B(P_p(G)^2)<=1/4 integral |q| |Ghat(q)|^2,
    limsup omega_B(Q_p(F)^2)<=3 integral |Fhat(q)|^2/|q|.  (9)

Fourier integrals use d^3q/(2pi)^3. Start with smooth compact Fourier
support separated from zero, then use the same density envelopes to
extend to smooth compact spatial tests. These are second-moment upper
bounds, not Gaussian covariance formulae.

## 4. The direct H_B history limit

The same-site identity [Q_p(F),P_p(G)]=iM_p(FG) and the tensor-product
Weyl remainder O(p) are source operator identities, independent of B.
Equations (3),(5),(7),(9) discharge the four inputs of the full-history
construction in `spin_canonical_field/README.md`: scalar symplectic phase
through all fixed words, uniform energy, vacuum regularity, and strong
regularity on history vectors. Choose one common subsequence for the
countable test/time word family and the bounded bond expectations.
Gram positivity gives a nonzero regular Weyl net with coefficient m>1/10.
Its state is stationary; its strongly continuous time generator is
nonnegative, because every approximating ground generator is K_B/(Jp).

The onsite field can be removed in an interaction picture for the
Lieb--Robinson estimate; this changes neither bond supports nor their
norm J. The same source bound therefore gives a finite l1 cone of speed
72 in X=px, T=Jpt. Even-parity lattice translations approximate all
continuous spatial translations, with test continuity following from
(7),(9). These translations commute with the direct H_B time limit.

This conclusion is existence along a common subsequence. It supplies
neither a convergence rate for the continuum histories nor their
uniqueness. It does not identify them with zero-field-clock histories.
The scalar quadratic coefficients in (4) do not establish a closed
linear field equation: the exact hidden-force control in
`spin_twist_closure` excludes that inference from quadratic twist
identities alone. No mass shell, Lorentz covariance, scattering matrix
or physical particle assignment is inferred here.

`SECOND_PAIR.md` proves the second canonical pair on this same joint
history subsequence. It controls the mixed uniform-momentum commutator
through a transverse auxiliary phase and then a longitudinal phase, each
on all preliminary energy-bounded histories. Its cubic source expansions
also identify rotated momentum exponentials, giving the full internal
circle action on both pairs without assuming convergence of sums of
unbounded generators. That proof requires the positive-B gap and the
family (2); it makes no assertion about the other order of limits.

## 5. Finite preparation and complete local records

At any specified n and rational error, the fixed-positive-field marginal
planner supplies a finite source radius for any observation region. Use
an even outer radius, so the literal source graph's root bipartition is
the same epsilon_x as in (2). Its site and edge counts are

    N=(2R+1)^3, E=3(2R)(2R+1)^2.

The finite adiabatic preparation in `spin_twist_closure/FINITE_FIELD_GAP.md`
starts from the staggered X product. With D=JE, g=B/8 and A=JE+BN/2,
the rank-one Jansen--Ruskai--Seiler bound gives

    C=2D/g^2+7D^2/g^3, T=3C/eta,
    M=ceil(3(DT+A^2T^2)/(2eta)).                          (10)

The adiabatic and piecewise-constant/product errors are at most eta/3
each. The literal source circuit has at most

    16N+2N_odd+M(62E+11N)

native pulses. Choose rational port precision b so 29N_pulses 2^-b<=eta/3.
All native failure records are retained. Equation (10) is a finite
resource bound, not an executed large ramp or an efficient laboratory
protocol. `ground_preparation_counts` is the same literal ramp formula
evaluated from cube counts, without allocating the source Hilbert space.

`local_preparation` divides the requested half-trace error between this
preparation and the fixed-B marginal estimate. Trace-distance contraction
then bounds every complete local instrument and its retained output,
including an initially independent reference. No measurement outcome is
discarded and no target probability is renormalized. Exterior degrees
of freedom are excluded from this local-output comparison. Native
history clock and detector compilation errors add separately; their
budget cannot be replaced by an equality to H_0 dynamics.

The construction supplies a computable finite source for each chosen
positive-field scale and local accuracy. It supplies no effective error
to the limiting continuum state or to nature's measured spectrum.

## 6. A complete local history on that prepared source

`construction.probe_history` joins the preceding bounds for specified local
encoded-Pauli probes and rational clock durations. Each probe is the literal
native pointer instrument in `native_spin_observer/circuits.py`; all pointer
and native failure records are retained. The exact source frame fixes any
physical-spin sign relabeling. The target is the complete local history in
omega_B with H_B time.

Let the probes occupy s sites in radius r, let there be k interventions,
and let t be the sum of absolute physical durations. The inherited complete
history cutoff is bounded in unhalved trace norm by

    C_R=min(2,(2k+1)6Js(2R+1)^2 t
                   2^[-(R-r)+ceil(96Jt)]).              (11)

The selection field leaves this bound unchanged in its onsite interaction
picture. Choose R until C_R<=eta/4. Prepare a larger even-radius box so its
entire R-box marginal has half-trace error at most eta/4 relative to omega_B;
this error includes both the actual finite ground preparation and the
fixed-B marginal error. Run the retained H_B clock on this larger source
with the count-only native error budget eta/4 from `CLOCK_COMPARISON.md`.
Compile the prescribed probes to complete isometry error at most eta/4.

For the state comparison, insert the R-box clock twice. Relative to the
infinite target it costs C_R in unhalved trace norm; relative to the exact
larger-box history it costs another C_R. Together they cost C_R in half
trace norm. The inner marginal comparison contracts through the complete
inner protocol. Adding the finite clock and probe compilation errors gives
a total half-trace error at most eta, including all local retained outputs,
pointer records, failure records and an initially independent reference.
No source exterior is part of this output comparison.

This is an explicit finite protocol for each selected n and error. The
resource planner evaluates integer counts without allocating its enormous
finite source or executing its preparation ramp. It does not turn a
subsequence existence theorem into a computable continuum error bound.
