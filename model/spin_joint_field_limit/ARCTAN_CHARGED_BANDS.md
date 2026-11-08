# Charged energy–momentum bands on the positive-field source

Use the same source family, selected ground and clock as `PROOF.md`:

    p=n^-2, B=J n^-5, L_n=(H_B-E_B)/(Jp), n>=2.

The conclusion concerns the two source-even collective charge channels in
that joint history limit. It does not transfer the zero-field ordered state
to this family, identify the original odd CAR fields with these channels,
or assert a sharp dispersion, Lorentz covariance or a physical species.

## 1. Energy bounds in this state and clock

For either transverse spin component let
`kappa_B=-omega_B(S_x^1 S_y^1+S_x^2 S_y^2)` on a nearest-neighbour bond.
The internal circle, cubic symmetries and real positive-field ground give
the same diagonal transverse covariance on every bond. On the actual
positive-field ground the exact first moment is

    <A_f(H_B-E_B)A_f>
       =(J kappa_B/2) sum_edges(f_x-f_y)^2
          +(B m_B/2) sum_x f_x^2,                     (1)

where `A_f=sum f_x S_x^a`, `a=2,3`, and `m_B=omega_B(epsilon_x S_x^1)`.
For staggered coefficients the first difference is the corresponding sum.
The operator bound gives `kappa_B<=1/2`. The finite-field susceptibility
bound consumed in `PROOF.md`, the canonical commutator and Cauchy–Schwarz
also give `m_B^2<=2 kappa_B+B m_B/(6J)` by taking a uniform momentum packet.
Thus a common subsequence has `m>=m0=1/10` and
`m^2/2<=kappa<=1/2`. No stiffness value is fitted.

The selection terms in (1), after scaling, are proportional to
`B/J=n^-5` for Q and `B/(Jp^2)=n^-1` for P. They vanish. Applying the
Stone-vector enlargement and reducing projection in
`spin_canonical_energy/README.md` to this joint family gives

    int E dnu_P,G <= (kappa/2) int |k|^2 |Ghat|^2,
    int E^-1 dnu_P,G <= (1/12) int |Ghat|^2,
    int E dnu_Q,F <= 6 kappa int |Fhat|^2,
    int E^-1 dnu_Q,F <= int |Fhat|^2/|k|^2.             (2)

All integrals use `d^3k/(2pi)^3`. The inverse score takes value infinity
at zero, so (2) excludes a zero-energy part in these field vectors. This
argument uses the same positive L_n; no H_B/H_0 clock comparison occurs.
Projection can discard source-vector norm and first energy moment, which
is why (2) consists of inequalities. The field vectors are in the original
bounded-history representation, rather than only an auxiliary Gram space.

Fix a real unit-L2 Schwartz F whose smooth Fourier support lies in
`r/2<=|k|<=r`, and put

    X=sqrt(r) Q3(F), Y=sqrt(r) Q2(F),
    P=P2(F)/sqrt(r), ell=L/r.

Here P is a canonical internal-spin field, not spatial momentum. The same
state obeys

    ||X Omega||^2, ||Y Omega||^2 <=6,
    ||P Omega||^2 <=1/4,
    int E dnu_P^ell <=1/4, int E^-1 dnu_P^ell <=1/12.  (3)

`SECOND_PAIR.md` supplies the full circle on both pairs. The canonical
P group shifts X by its scalar CCR coefficient and fixes Y. Even microscopic
translations give the common spatial translations; their spectral measure
places P Omega in the stated momentum annulus. An odd microscopic
translation combined with the internal pi rotation about axis two preserves
H_B and its unique ground. Its continuum internal action sends
`(X,Y,P)->(X,-Y,P)` and reverses circle charge. It commutes with the limiting
clock and even-translation spatial group. The infinitesimal spatial shift
vanishes by the proved test continuity. This supplies equal opposite-charge
weights on this joint family without replacing its state.

## 2. A bounded monotone creator

Define

    f(t)=atan(t)/pi
        =(1/pi) int_0^infinity e^-u sin(ut) du/u,
    X_theta=X cos(theta)+Y sin(theta),
    B_+=(1/2pi) int_0^(2pi) e^-i theta f(X_theta) dtheta,
    B_-=B_+*, A=(B_++B_-)/2.                            (4)

The bounded functions have norm at most 1/2. Circle covariance makes B_+
and B_- exact charge +1 and -1 operators. Their complete preparation is

    B_+, B_-, B_0=sqrt(I-B_+*B_+-B_-*B_-).             (5)

The expression under the root is at least I/2, including on the finite
source where X and Y need not commute. No complementary record is removed.

The canonical shift and the bounded derivative of f give the weak identity

    |< [P,A] >|=(m/(2pi)) <(1+X^2)^-1>
                >=m/[2pi(1+<X^2>)] >=m/(14pi).        (6)

For the equality, differentiate the bounded function under the canonical
shift, then average cos^2(theta). Circle invariance makes the derivative's
expectation independent of theta. The last step is Jensen's inequality
for the convex function `(1+t)^-1`, t>=0. Only P Omega and finite second
moments are needed for this weak commutator. No Gaussian law or radial
configuration density enters the estimate. The positive derivative is
useful here: it supplies a bound without finding a favourable sine amplitude.

Use the rational lower bound `d=m0/44`, from `pi<22/7`, and set

    a=3d^2/2, b=2/d^2,
    Pi=1_(r/2<=|spatial momentum|<=r) 1_(a<=ell<=b).    (7)

Equation (3) gives `||(I-Pi)P Omega||^2<=a/12+1/(4b)=d^2/4`.
Since `||A||<=1/2`, the commutator contribution outside Pi is at most d/2.
The inside part is at most `2||Pi P Omega|| ||Pi A Omega||<=||Pi A Omega||`.
Thus `||Pi A Omega||>=d/2`. Charge orthogonality and the conjugation from
section 1 imply

    ||Pi B_+ Omega||^2=||Pi B_- Omega||^2
                      =2||Pi A Omega||^2 >=d^2/2.     (8)

For the actual certified `m0=1/10`, each original charge record has
probability at least `1/387200` jointly with momentum in `[r/2,r]` and
energy in `[3r/387200,387200r]`. All preparation, energy and momentum
complements remain allocated. These are broad bands on one state and clock,
not a pole or a rest-mass determination. The lower probability is uniform
as r decreases. Thus this joint family has positive collective response at
arbitrarily small energy and momentum scales, despite the separate odd-CAR
escape bound on a different source family.

The same source sine energy bound and (4) give
`||ell^(1/2) B_+ Omega||,||ell^(1/2) B_- Omega||<=sqrt(3)/pi`.
The integral is finite because its factor 1/u cancels the linear sine
energy bound. This bounds the entire created response, not only Pi.

## 3. The finite-source commutator error

The monotone construction admits a source estimate before the CCR limit.
For summable real F define the finite/source operators with the same
normalization as section 1, and put

    M=p^3 sum_x epsilon_x F(px)^2 S_x^1,
    mu=<M>, sigma^2=<(M-mu)^2>,
    C=(sqrt(r) p^5/2) sum_x |F(px)|^3.                 (9)

Exact spin algebra gives `[P,X]=-iM`, `[P,Y]=0` and
`||[M,X_theta]||<=C`. Duhamel's formula, applied before taking an expectation,
therefore gives

    |< [P,e^(iu X_theta)] >
         -u cos(theta) mu <e^(iu X_theta)>|
       <=|u cos(theta)| (sigma+|u|C/2).                (10)

To check the error, move `M-mu` to the right of the second exponential.
Its action on Omega has norm sigma; its commutator with that exponential
has norm at most the time parameter times C. Integrating the intermediate
time gives C/2. Integrate the sine form of (10) against `e^-u/(pi u)`, then
multiply by cos(theta) and average. The error is `(sigma+C/2)/(2pi)`.
Circle invariance and Jensen consequently give the rigorous lower bound

    |< [P,A] >|
      >=max(0,[mu/(1+<X^2>)-sigma-C/2]/(2pi)).         (11)

For the actual positive-field ground the gap B/8 and the exact local
double commutator give

    sigma^2 <= (2Jp^6/B) sum_edges(F(px)^2+F(py)^2)^2. (12)

The selection term commutes with M. For a Schwartz F the sum in (12) is
O(p^-3), hence sigma^2=O(1/n); C=O(p^2). The expectation is exactly
`mu=m_B p^3 sum F(px)^2` on the infinite homogeneous ground. Equations
(9)-(12) prove the finite-source approach to (6), on this family's state
and clock. Finite-box marginals and truncation errors add under the existing
positive-field locality bound. A finite source can also use its actual
mu, variance and covariance scores in (11); no continuum order is inserted
into a small box without its marginal error.

At finite N the circle Fourier integral is exactly a discrete sum with
`2N+1` equally spaced phases: all operator charge differences lie in
`[-N,N]`, so the target charge has no modular alias. The onsite spin sum X
has a known finite diagonalization. Its bounded arctangent can therefore
be evaluated with certified scalar intervals and compiled as a complete
finite isometry by the stated native coherent-control service. This is
not a positive-weight LCU claim for the divergent absolute coefficient
integral at u=0. The Laplace integral in (4) proves the analytic bounds;
the finite spectral table specifies the instrument. No large native tape
is executed by the bound producer.

## 4. A bounded original-probability Ward test

Let Pi be any actual joint spectral projection commuting with the source
circle and charge conjugation. Put `v=||P Omega||^2`,
`eta=||(I-Pi)P Omega||`, and let `d<=|<[P,A]>|` be a certified lower
bound, from (6) or its finite-source version (11). The same decomposition
used in section 2 proves

    ||Pi B_+ Omega||^2=||Pi B_- Omega||^2
       >=max(0,d-eta)^2/(2v),   v>0.                  (13)

The outside contribution is at most eta because `||A||<=1/2`. The inside
contribution is at most `2 sqrt(v)||Pi A Omega||`. Charge orthogonality
then gives (13). Upper bounds on v and eta and a lower bound on d preserve
the inequality. These are original probabilities of the bounded instrument
(5); there is no Pauli-sum preparation normalization growing with volume.

For a concrete Ward input choose a compact original source energy window
R supported in `[ell_0,M]`, within the actual momentum annulus. Use
`WARD_DYNAMICS.md` to measure the two filtered inverse susceptibilities,
first moments and original cross moment c_R for the pair X,P. They determine
`u=c_R/chi_P`, `b=c_R/chi_Q`, `omega^2=ub>0` and the positive force slacks
`delta_Q,delta_P`. Let Pi retain `|E-omega|<epsilon` inside R. For sharp
commuting window projections the existing force estimate gives

    eta^2 <= ||(I-R)P Omega||^2
          +M (|b|sqrt(delta_Q)+M sqrt(delta_P))^2
                                     /(omega^2 epsilon^2).    (14)

All quantities in (14) use the same clock and retained filtered cross
moment. Smooth detector effects require the existing boundary-weight and
resolution errors before a spectral projection is inferred. Equations
(13)-(14) turn a small measured or proved Ward-tube leakage into an original
charged preparation probability. The broad-band theorem establishes no
vanishing force slack or shrinking relative spectral width. In particular
the nonlinear bounded creator can have multiparticle spectral weight even
in a free field. Only the selected projection's original weight enters
(13); the complete created state is not required to be a single pole.

The actual uniform dimer in section 6 saturates (13) with Pi its complete
middle-energy projector. It has `eta=0`, `v=p^2(1-s)` and
`d=p k/sqrt(257)`, and

    d^2/(2v)=2k^2 s.

This is exactly each original charged record's probability, with both
complementary energy outcomes retained. No continuum saturation follows
from this finite equality.

## 5. Spatial cutoff and finite resolution

For `F_R=chi_R F` the configurations at a fixed angle commute. The scalar
Lipschitz constant of f is 1/pi, so

    ||(B_+(F)-B_+(F_R))Omega||
       <=sqrt(r)/pi ||Q3(F-F_R)Omega||.                (15)

The same estimate holds for B_-. In the limiting commuting configuration
algebra, (5) and its cutoff have complete isometry difference at most twice
the right side of (15): the two charged amplitudes each cost delta and the
complement costs at most sqrt(2)delta, because its square is bounded below
by 1/2. The inherited covariance bound makes delta tend to zero for Schwartz
tails. This supplies strictly local approximations retaining all three
records. It is a vacuum-input statement with an initially independent
reference, not a channel norm bound on arbitrary correlated source states.

A continuous energy–momentum effect equal to one on the compact band (7)
and supported in a slightly larger positive band has at least the weight
(8). Its translation-average approximation gives finite spatial and clock
resolution with the existing complete detector services. Boundary weights
and service errors must be charged before inferring a sharp spectral event.
The analytic joint-history limit gives no effective n-to-continuum rate or
unique limiting state; the band theorem does not supply those stronger claims.

## 6. Complete source control

The executable reconstructs the actual two-site `native_ordered_ground`
Hamiltonian at n=2, J=1, B=1/32 and p=1/4. A global rotation places the
selection axis along Z; the clock, state and charge are rotated together.
Set `k=atan(1/16)/(2pi)`, `s=(1+16/sqrt(257))/2`. Then

    B_+=k(S_0^+-S_1^+), B_-=B_+*,
    H_B spectrum={-(16+sqrt(257))/32,0,0,(-16+sqrt(257))/32}.

Each charge branch has original probability `2k^2 s` at scaled energy
`(16+sqrt(257))/8`. Put `c=sqrt(1-4k^2)`. The retained complementary branch
has ground probability `[1-s(1-c)]^2` and upper-energy probability
`s(1-s)(1-c)^2`; its middle-energy probability is zero. Every one of the
nine preparation/energy labels is retained. The probabilities sum to one
as an identity, and the full effect matrices sum to identity on all inputs.

Rational Machin-series bounds for pi, alternating arctangent bounds and
integer square-root enclosures certify the emitted intervals. The complete
dimer control tests source signs, charge grading and complement accounting;
it is not an annular continuum packet or a demonstration of a physical mass.
