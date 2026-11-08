# Strictly local charged preparation on the same field source

## 1. A compact source profile with retained overlap

Fix a real unit-L2 Schwartz test F with smooth Fourier support in
1/2<=|k|<=1. All Fourier integrals use d^3k/(2pi)^3. Let chi_R be a
smooth cutoff equal to1 for |X|<=R, equal to0 for |X|>=2R, and between
zero and one. Choose R>=1 so that the full tails of F outside R have
L1 and L2 norms at most1/4. Set F_c=chi_R F without renormalizing it.
Then

    ||F_c||2<=1,    integral F_c F>=15/16.                  (1)

Such an R is computable from supplied analytic Schwartz seminorms.
For example, |F(X)|<=C(1+|X|)^-6 gives the explicit sufficient bounds

    4pi C/(3R³)<=1/4,    4pi C²/(9R^9)<=1/16.

The producer uses pi<22/7 and exact integer comparisons. These are
analytic bounds on the chosen test, not estimates from sampled tails.
The source state is not an input to this cutoff choice.

For h=F-F_c, the existing complete charged-preparation cutoff theorem
in `spin_energy_polarization` proves

    integral |hhat(k)|²/|k| <= ||h||1²/(4pi²)+||h||2²
                            <=37/576.                    (2)

Since the parent annular norm is at most2, the weighted triangle
inequality gives a source variance bound for each Q component:

    3 integral |Fchat(k)|²/|k|
       <=3(sqrt(2)+sqrt(37)/24)²
       <=5547/576<12.                                    (3)

The complete three-record vacuum-preparation truncation also has its
inherited state-specific isometry error2|s|sqrt(V(h)). No operator-norm
approximation of an unbounded smeared field is inferred. In particular,
the fixed tail allowance1/4 is not used to transfer a tiny probability
floor by continuity; the local creator's positive floor is proved directly.

## 2. A local positive energy-momentum band

For any r>0 define

    F_r(X)=r^(3/2) F(rX),
    F_(c,r)(X)=r^(3/2) (chi_R F)(rX),
    x=sqrt(r)Q3(F_(c,r)),    y=sqrt(r)Q2(F_(c,r)),
    p=P2(F_r)/sqrt(r),    ell=L/r.

The configuration test has support in |X|<=2R/r. The original probe
p retains its annular spatial spectral support, norm squared<=1/4,
first ell moment<=1/4 and inverse ell moment<=1/12. Equation(3) bounds
each of ||xΩ||²,||yΩ||² by12. The actual Weyl shift coefficient is
m c with c=integral F_c F>=15/16; p fixes y. No extra momentum operator
or rotation of p is supplied.

Apply the existing bounded angular charge projection to this local pair:

    B_(c,+)=(x-i y)/[sqrt(1+x²+y²)(sqrt(1+x²+y²)+1)],
    B_(c,-)=B_(c,+)*.

The complete third Kraus operator is the positive square root of
I-B_(c,+)*B_(c,+)-B_(c,-)*B_(c,-). All three lie in the local
configuration algebra of the stated compact region. Each charged
operator has norm<=1/2 and exact internal charge±1.

With A=Re B_(c,+), the same weak canonical differentiation and
configuration-vacuum rotation invariance give

    |< [p,A] >| >= m c/[2(1+24)^(3/2)]
                 >= d0=3m0/1000,       m0=1/10.           (4)

This uses only the source p-shift and the configuration distribution;
it does not extend the circle symmetry to the full canonical net.
Set a=3d0²/2, b=2/d0² and let Pi select
r/2<=|P|<=r and ar<=L<=br. The source p-tail bound is

    ||(I-Pi)pΩ||²<=a/12+1/(4b)=d0²/4.

The weak commutator and ||A||<=1/2 imply ||Pi AΩ||>=d0/2.
Orthogonal opposite charges and the source conjugation then prove

    ||Pi B_(c,±)Ω||² >= d0²/2
       =9m0²/2000000=9/200000000.                         (5)

This is an original complete-instrument probability, uniform in r.
The nonlinear creator can have other momenta and energies; all those
outcomes and the zero preparation record are retained. The inherited
all-input form estimate applies with ||F_c||2<=1, giving the full
unnormalized first energy moment<=3r. Thus the actual field has strictly
local bounded charged interpolating operators with nonzero positive
energy-momentum response. No pole, velocity or physical particle label
is inferred from this band.

## 3. Energy-selective, differentiable source operators

Choose a smooth compactly supported real frequency function h equal to1
on [a,b] and supported in (a/2,2b). With Fourier convention
h(E)=integral g(t)exp(itE)dt, its inverse g is Schwartz. Put
h_r(E)=h(E/r), g_r(t)=r g(rt), and choose a computable
M>=max(1,||g||1). Define actual source operators

    C_(r,±)=M^-1 integral g_r(s) alpha_s(B_(c,±)) ds.      (6)

Strong operator integration suffices. Their norms are at most1/2,
charges are±1, and the positive complement completes the three-record
instrument. Charge conjugation interchanges the two raising operators.
They are not asserted to be each other's adjoints: the adjoint of a
positive-energy-transfer operator lowers energy.

On the vacuum, C_(r,±)Ω=M^-1 h_r(L)B_(c,±)Ω. Consequently the original
probability in the inner joint band is at least(5)/M². The source clock
and momentum detector are unchanged.

The spectral transfer of (6) lies in [ar/2,2br]. Explicitly, if bounded
energy intervals I,J have (J-I) disjoint from that interval, then
E_L(J) C_(r,±) E_L(I)=0. This follows by Fourier smearing the left/right
spectral evolution. For a direct justification, test on vectors with
compact energy support: their matrix coefficient is entire of the
corresponding difference exponential type. The inverse transform of h_r
has rapid horizontal decay and exponential type fixed by its compact
frequency support. Shifting the integration contour upward or downward
through a positive separating energy gap makes the smeared coefficient
zero. Approximation by compact intervals gives the stated projection
identity. No positive-frequency assumption is imposed on B itself.

Time differentiation is bounded:

    i[L,C_(r,±)]=-M^-1 integral g_r'(s)alpha_s(B_(c,±))ds,
    ||i[L,C_(r,±)]|| <=||g_r'||1/(2M).                    (7)

Changing variables under a time translate and differentiating g_r in L1
proves(7); it does not assume norm continuity of the unsmeared field.
Higher time derivatives follow from higher L1 derivatives of g_r.

## 4. Source locality gives controlled tails

The actual field net commutes outside the inherited l1 cone of speed72
in the same clock. Let two compact preparation supports have l1 distance
d>0, put v=72, and T=d/(4v). Truncate both time integrals in (6) to
|s|<=T. Every pair of retained source factors then commutes, since their
time separation is at most2T<d/v. No time-slice equality of local algebras
is needed for this argument.

Write epsilon_g(T)=integral_(|s|>T)|g_r(s)|ds. The truncated and full
operator norms are at most1/2, and their difference is at most
epsilon_g(T)/(2M). Therefore, for either charge and also their adjoints,

    ||[C_(r,sigma)(x),C_(r,tau)(y)]||
        <=epsilon_g(d/(4v))/M.                            (8)

The same expansion with g_r' bounds a derivative commutator by

    epsilon_(g')(T)/(2M)
       +||g_r'||1 epsilon_g(T)/(2M²).                     (9)

Every tail decreases faster than any inverse power at fixed r. For
translated support centres in max-norm separation D, the actual supports
have l1 separation at least D-4R/r when this is positive. These are
operator bounds on every input and untouched reference. They are
stronger than a vacuum-only spatial-tail estimate.

Equations(5)-(9) provide source-derived local positive charged inputs,
energy-transfer filters, bounded time derivatives and almost-local
commutators. They do not select a one-particle dispersion. A packet
comparison must additionally match the actual L/P spectral response and
bound its one-packet dynamical defect; locality alone controls only the
inter-packet commutator terms.

## 5. Analytic filter normalization

For the Fourier convention in section3, integration by parts gives

    ||g||1 <= (||h||1+||h''||1)/pi.

The producer bounds pi from below by3 and uses the resulting rational
upper bound for M. It constructs a compact piecewise-polynomial filter
with any integer smoothness K>=2. Its rising transition is

    s_K(u)=[(2K+1)!/(K!)²] integral_0^u t^K(1-t)^K dt.

Use s_K((E-a/2)/(a/2)) on [a/2,a], the constant1 on [a,b],
s_K((2b-E)/b) on [b,2b], and zero outside. Every endpoint derivative
through orderK vanishes. The weak derivatives through orderK+1 are
integrable and their L1 norms have exact coefficient bounds. Therefore

    integral_(|s|>T)|g_r(s)|ds
       <= ||h^(k)||1/[pi(k-1)(rT)^(k-1)],
       2<=k<=K+1.                                         (10)

For g_r', replace h by E h and multiply the bound by r. The same
construction bounds ||g_r'||1 through the first and third expressions
(Eh) and (Eh)''. These are analytic bounds, with no sampled Fourier
integral used as a norm estimate. Choosing K as large as needed supplies
any specified inverse-power locality estimate. The single C-infinity
filter in section3 supplies all powers simultaneously; the finite
polynomial examples claim only their certified powers.
