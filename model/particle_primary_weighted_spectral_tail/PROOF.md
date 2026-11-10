# Selected current moments and finite energy-transfer windows

The actual critical phase apparatus has second energy-transfer moment
at most `4000000/a`, with `a=2^-n`, independently of link degree `ell=n+1`.
This state-weighted estimate uses stationarity, the selected polynomial
populations and the literal local source action. It improves the
full-CAR operator cap of order `a^-4`. A complete original positive band
has a certified weight of order at least a, with controlled clocks.
Its upper energy endpoint grows as `1/a`; no fixed positive-energy
survival or physical mass follows.

## 1. Exact signed third moment on the selected state

Keep the [primary source contract](../particle_common_contract/PROOF.md),
its actual selected polynomial density rho, physical Gauss space,
retained canonical reference and Hamiltonian clock. The
[critical phase apparatus](../particle_primary_critical_refinement/PROOF.md)
retains all V origins and both signs with weight1/(2V). At an origin X,
write `U_sigma=exp(i sigma F_X)`, where
`F_X=pi a sum_x f_R(x-X)(N_x-30)` and `R=1/a`. Define

    D_sigma=U_sigma* H U_sigma-H=C+sigma S.                 (1)

Both C and S are Hermitian. They consist only of hopping terms, since
the phase commutes with the other source terms. On a bond b, with phase
difference delta_b, C has coefficient `cos(delta_b)-1` multiplying the
Hermitian hop, and S has coefficient `sin(delta_b)` multiplying its
Hermitian current. No Taylor replacement of these operators is made.

Let m_r be the signed r-th moment of the complete ideal two-energy law
of this menu, with every initial energy, final energy, origin and sign
retained. The stationary source marginal satisfies

    M2=m2=average Tr rho D_sigma²=average Tr rho(C²+S²),
    m3=average Tr rho(half[D_sigma,[H,D_sigma]]+D_sigma³).
                                                               (2)

Expand `U_sigma* H³ U_sigma=(H+D_sigma)³`, subtract the other terms
in `(E_final-E_initial)³`, and move H through rho to prove (2).
It uses the original density; it is not an expectation in a trial vector
or in a Haar state. Initial energy recording preserves this source
marginal and can change its correlations with the retained reference.

Average the two signs. The exact identity

    average_sign D_sigma³=C³+3SCS+[S,[S,C]]                 (3)

avoids an absolute bound on S³. If `||C||<=c`, then
`C³<=c C²` and `SCS<=c S²` as quadratic forms. Consequently

    m3 <= A+3c M2,
    A >=half average_X Tr rho([C,[H,C]]+[S,[H,S]])
          +average_X Tr rho[S,[S,C]].                       (4)

This is a closed inequality in the actual selected M2. An absolute
operator bound on M2 itself is not substituted for its right-hand side.
Local norms provide one A. The selected electric-energy estimate below
bounds the expectations more sharply than the electric operator norm.

## 2. Local bounds on the actual product-cutoff carrier

Use the owner-site factors of the [source locality proof](../particle_primary_locality/PROOF.md).
Matter and electric terms have one owner, hopping has two and a plaquette
has at most three. Their norms, before restricting to physical states,
give the per-owner sum

    sum_(t touching x) ||h_t|| <= kappa_H/a,
    kappa_H=270+9/7+(1827/2)ell² <1200ell².                (5)

Here each cell meets at most six hopping terms of norm30/a, nine
plaquettes of norm1/(7a), one matter term of norm90/a, and three
electric terms of norm `609ell²/(2a)`. All dangling-link terms may be
retained in this upper bound. The full word projection is a product of
link projections, so it preserves these supports and does not increase
the displayed norms. The physical Gauss restriction is applied after
the local estimates, and is not treated as a tensor-local projector.

The following finite support bound is useful. Let `Q=sum_b Q_b` be
Hermitian even bond terms, with `d_b=||Q_b||`. Suppose each owner meets
at most six such bonds and an interaction H has supports of size at
most three and per-owner norm sum kappa. Then

    ||[Q,[H,Q]]|| <=240 kappa sum_b d_b².                 (6)

To verify the constant, expand the two commutators. A term indexed by
outer bond b, interaction t and inner bond c can be nonzero only if
t touches c and b touches `t union c`. Its norm is at most
`4 d_b ||h_t|| d_c`. For fixed c, t has total norm at most2kappa,
and its union with c touches at most30 outer bonds. The resulting
matrix K_bc has column sums at most60kappa. For fixed b, split the
same incidence into b touching c, or b touching t. The first contributes
at most `12*2kappa`; the second at most `2kappa*18`. Its row sums are
also at most60kappa. Applying `d_b d_c<=(d_b²+d_c²)/2` proves (6).
Extra commutator cancellations are not required.

For the actual tent, `|delta_b|<=pi a²`. Thus an individual cosine
bond has norm at most150a³, so its per-owner norm sum is at most
`kappa_C=900a³`. The complete cosine operator satisfies

    ||C_X||<=c_n,
    c_0=180,    c_n=150a D2 for n>=1,
    c_n<=2025/4.                                          (7)

D2 is the exact sum of squared discrete gradients of the three-dimensional
unit tent on the periodic coefficient list. The actual Hamiltonian is
open; only the apparatus coefficients wrap. With B physical interior
bonds and V cells, averaging origins gives

    average_X sum_b (||C_b||²+||S_b||²)
       <=1800 B k_n/(a²V)
       <=J2_n :=3000(B/V)D2.                              (8)

The first bound uses `(1-cos delta)²+sin² delta=2(1-cos delta)` and
the literal full-Fock hop/current norm30/a. The second uses the exact
critical profile coefficient `k_n<=5a²D2/3`. For every origin the
source bonds are a subset of the coefficient-list bonds, so boundary
loss cannot invalidate (7). No periodic source dynamics is introduced.

Apply (6) to C and S with H, then to S with interaction C. Equations
(4)-(8) imply the actual state-weighted inequality

    m3 <= A_local+B_n M2,
    A_local=(120 kappa_H/a+240 kappa_C)J2_n,
    B_n=3c_n.                                            (9)

Since `a D2<=6` and `B/V<3`, these explicit source coefficients obey

    A_local<=19440000000 ell²/a²,    B_n<=1519.           (10)

The local estimates remain valid after the link cutoff and on every
physical occupation sector. They do not bound the selected electric
covariance by an imported continuum or ground-state assumption.

## 3. The actual electric energy removes link-degree growth

The filled negative-beta vector with constant links is a physical trial
of mean `-90V/a` at every member. The selected cooling bound gives
`Tr rho H<=-90V/a+epsilon`, where `epsilon=a/16`. Onsite matter is
bounded below by `-90V/a`, all hopping by `-90V/a`, and all magnetic
terms by `-3V/(7a)`. Thus the actual positive electric energy satisfies

    Tr rho H_electric <=(90+3/7)V/a+epsilon <91V/a,
    sum_b Tr rho Casimir_b <=182V/105.                    (11)

This uses the unchanged polynomial state, not the trial preparation.
The sum may include all allocated links; dropping nonnegative terms
retains the bound for the interior bonds.

Choose the twelve Hermitian link derivatives L_(b,A) orthonormal in the
actual native trace metric. Their squares sum to Casimir_b. The
[native generator calculation](../particle_gauge_characteristics/PROOF.md)
verifies that their fifteen-state matrices T_A obey
`sum_A Tr(T_A²)=Tr(Casimir_M15)=12`. Hence for `t_A=||T_A||`,
`sum_A t_A²<=12`. These are analysis operators on the product carrier,
not declared gauge-invariant detector operations. Each preserves the
link-word cutoff, so the projection P commutes with L_(b,A).

For Q=C or S, put `d_b=||Q_b||` using its displayed coefficient norm
cap. Literal differentiation of the native link matrix and the full-CAR
bond norm give

    ||[L_(b,A),Q]||=||[L_(b,A),Q_b]||<=t_A d_b.

The same inequality holds after projection because
`[L,PQ_bP]=P[L,Q_b]P`. With Hermitian L and Q, set `J=i[L,Q]` and
`D=[Q,[L,Q]]`. The exact product-carrier identity is

    [Q,[L²,Q]]=2J²+{L,D}.                               (12)

Only the at most twelve Q bonds touching b enter D, so
`||D||<=24 t_A d_b d_max`. Cauchy--Schwarz in the original rho gives
`|Tr rho{L,D}|<=2||D|| sqrt(Tr rho L²)`. This remains valid when rho
is extended by zero off the physical carrier. It does not assume that
an individual L preserves Gauss support.

For the unit tent, `d_max,S<=30pi a` and `d_max,C<=150a³`.
For every fixed physical bond, averaging the V origins yields

    average_X d_(b,S)<=60pi/(V a²),
    average_X d_(b,C)<=50a D2/V.

The first bound uses the exact axis absolute-gradient sum `2R²`;
the second uses the axis square-gradient sum `D2/3`. By (11),
`sum_b sqrt(Tr rho Casimir_b)<=V sqrt(182/35)`.
Summing (12) for both C and S, using `pi²<10` and
`sqrt(12*182/35)<8`, therefore proves

    average_X Tr rho([C,[H_electric,C]]+[S,[H_electric,S]])
       <=E_n :=105/(2a) *
             (24J2_n+6912000/a+2880000a^4 D2).             (13)

The finite link cutoff contributes no factor ell in (13). The derivative
and commutator identities were proved on the actual compressed product
carrier; no commutation of P with multiplication was assumed.

Apply (6) to the nonelectric action, with per-owner norm sum
`kappa_0/a`, `kappa_0=270+9/7`. Equation (4) holds with

    A_selected=120 kappa_0 J2_n/a+E_n/2+240 kappa_C J2_n,
    A_n=min(A_local,A_selected),
    A_selected<=14100000000/a².                          (14)

For the last envelope use `J2_n<=54000/a`, `a D2<=6`, and a<=1.
The original electric expectation, rather than a uniform bound on every
electric eigenvalue, is what removes the cutoff-degree dependence.

## 4. The selected spectral inequality closes the second moment

The [actual polynomial adjoint law](../particle_primary_work_spectrum/PROOF.md)
pairs the two signs with population ratio at most
`exp(-beta_low omega)` for omega>0. For each pair, its negative fraction
is at most `1/(1+exp(beta_low omega))`. The same pairing that proves
the negative first-moment bound gives

    N1=integral_(omega<0)|omega| dmu <=1/(3 beta_low),
    N2=integral_(omega<0)omega² dmu <=2/beta_low²,
    N3=integral_(omega<0)|omega|³ dmu <=6/beta_low³.        (15)

For the latter two inequalities use `exp(x)>=x^r/r!`, r=2,3,
and the fact that all paired masses sum to at most one. No thermal
density is substituted for the selected polynomial. Write Wbar for
the authenticated signed work upper and `P1=Wbar+1/(3 beta_low)`.
For compact exact arithmetic the receipt weakens beta_low to its positive
integer floor and records that bound. This only enlarges the certified
negative moments and changes neither the state nor its temperature.
Cauchy--Schwarz on the positive part of the original measure, together
with (9) and (14), yields

    (M2-N2)_+² <=P1(A_n+B_n M2+N3).                      (16)

If M2<=N2 the bound below is immediate. Otherwise solving (16) gives

    M2 <=N2+B_n P1/2
          +sqrt((B_n P1/2)²+P1(A_n+N3+B_n N2)).          (17)

Every coefficient is regenerated from the original source member. The
square root and final bound are enclosed upward by exact dyadic rationals. The finite
code takes the smaller of (17) and the full-CAR operator cap from the
work-spectrum proof. This combines compatible upper bounds on the same
observable, state and clock.

The source bounds give P1<444 and N2,N3<1. The weaker expression
`N2+B_n P1+sqrt(P1(A_n+N3+B_n N2))` also bounds (17).
Use (10), (14), `1/a>=1`, `1+1519*444<700000` and
`444*(14100000000+1520)<3000000²`. The dyadic rounding error is less
than one, so a uniform family statement is

    M2 <=4000000/a.                                      (18)

This is not a uniform second-moment bound: 1/a diverges. It is a
stronger actual source estimate than order a^-4, obtained without
assuming unproved selected-current decorrelation.

## 5. Complete original finite windows and the remaining limit question

Let w_n>0 be the authenticated selected signed work lower. With the
smaller cap C_n supplied by (17)-(18), the elementary original-window
argument gives

    mu([w_n/4,4C_n/w_n]) >=w_n²/(4C_n).                  (19)

The lower and upper discarded signed first moments are each at most
w_n/4; Cauchy--Schwarz on the retained middle band proves (19). No
success conditioning is used. The same critical-profile two-clock plan
has resolution eta, one-clock tail tau and complete implementation
error e. The reported window expands each endpoint by2eta and retains
original probability at least `w_n²/(4C_n)-2tau-e`. The receipt checks
that its lower endpoint is positive and its probability exceeds half
the ideal floor. Every clock label, origin, sign, source output and
reference remains in the complete experiment.

Because the source work floor is uniform and positive, (18) supplies
a window floor of order at least a with an upper endpoint of order
1/a. The finite first-moment tail also satisfies
`integral_(|omega|>R)|omega| dmu <=C_n/R`. Neither estimate is uniform
in capacity at fixed R. A nonzero fixed-energy response or uniform
integrability needs a sharper selected correlation estimate. In this
derivation the specific unbounded term is the expectation of the two
local current double commutators in (4), bounded by O(a^-2) in (14).
Controlling that selected quantity uniformly would close (16) with a
uniform M2 cap. The present local norm bound does not supply that result,
and neither a gap nor a source elastic limit is asserted. Rare cutoff-scale
weight carrying finite work is compatible with these actual bounds;
excluding or proving that behavior needs further source correlations.
