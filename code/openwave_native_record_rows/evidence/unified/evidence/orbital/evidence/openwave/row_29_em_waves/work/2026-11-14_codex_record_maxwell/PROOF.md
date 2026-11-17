# A transverse Maxwell reconstruction from native occupation records

Written conditional source theorem; the finite companion experiment is
separate, and no owner pin or accepted row is asserted.
The construction uses four copies of a specified native-record experiment and
a local difference decoder. It derives the Gaussian electromagnetic field
strength and its two-state reflection quotient. It does not run a Maxwell update
on the source, identify native attempt order with physical time, or derive a
charged interaction.

## 1. Existing results used, and the new statement

The finite source theorem in
`openwave/shared/native_occupation_fluctuation_limit/PROOF.md` proves that the
normalized occupation block of one native unit converges jointly, including
every fixed moment, to a Gaussian vector with covariance

    C = L^+ - Pi/(4W).

Here L is the weighted graph Laplacian, W is total seam weight and Pi removes
the constant mode. This is a consequence of the native transition
P=I-L/(2W), rather than an installed Gaussian dynamics or a stipulated Wick
rule. The reviewed row30 slice theorem gives its four-axis routed version,
native internal eigenvalues, exact slice pole and uniform refinement bounds.

The old row29 conditional field theorem and E2 free calculation already give
the transverse projector and lattice pole for supplied field covariances. The
new composition here begins with the occupation-record Gaussian limit, defines
the complete readout, and proves reflection positivity on its field-strength
algebra. Four positive scalar channels do not automatically imply this last
claim: a vector reflection makes their time component negative. The explicit
quotient calculation below is needed.

## 2. Source, reference, channels and readout

Use the declared row30 family: twelve native ports at each point z of a
four-axis periodic lattice, thirty weight-one internal seams per cell and
forty-eight external seams per cell, each of weight a^(-2). The one-unit
transition is

    P_a=I-L_a/(2W),
    L_a=L_int+a^(-2) sum_mu(2-S_mu-S_mu^(-1)).

For a=1/m the weights are m² repeated external addresses. Each endpoint event
is the pinned native event. Its extension to this routed graph, refinement,
scheduler and observation registers is a proposed N+ apparatus, not the pinned
two-event controller A.

Run four copies, indexed mu=0,1,2,3, with independent initial uniform locations,
seam draws and placement coins. This is one explicitly declared product source.
It can be executed as four disjoint routed components containing one unit each,
with one native event in each component per round, using independent choices.
Each round consists of four ordinary endpoint events; it is a scheduler block,
not a physical-time increment. No interactions between components are claimed.

Each observer retains its disclosed initial location and the current tape.
Continuity reconstructs the occupation indicator. For n local events define

    Y_n^mu = sqrt(M/(4Wn)) sum_(r=1)^n [e_(X_r^mu)-1/M],

where M=12 times the number of lattice cells, and W is for one component.
The four block vectors have the independent joint Gaussian limit proved by the
occupation theorem. Their independence is declared apparatus, not a theorem
that four channels are required or uniquely selected by OPH.

Project each channel onto the normalized constant internal vector u=1/sqrt(12).
This selects the exact native internal eigenvalue zero. Embed channel mu at
the midpoint of the oriented mu-link and call its projected fluctuation A_mu.
The link-centre assignment is a readout convention. It is specified before any
data are inspected. At physical spacing a use field normalization a^(-2)u^T Y;
with the volume element a^4 its Fourier covariance has the symbol below.

The observer compares neighboring records and reports only

    F_mu,nu = D_mu A_nu - D_nu A_mu,

on plaquette centres. D is the oriented forward difference divided by a. All
these numbers are local linear functions of retained occupation records. The
source never reads F, implements an oscillator, or feeds a Maxwell equation
back into an event. The antisymmetric comparison, four channel labels, routing,
normalization and Euclidean slice assignment are disclosed readout inputs.

For nonzero lattice four-momentum let

    kappa_mu = (2/a) sin(a p_mu/2),
    Q_4² = sum_mu kappa_mu².

Midpoint Fourier phases make each difference multiplier i kappa_mu. The exact
Gaussian-limit covariance, with the usual reciprocal-volume normalization, is

    <A_mu(p) A_nu(p)*> = delta_mu,nu [1/Q_4² - 1/(4W)].      (1)

At the constant momentum the source mode is removed. The curl annihilates it,
so this removal introduces no constant subtraction in the F covariance. The
second term in (1) is the actual occupation correction and is retained. The
first term in the F covariance is

    [k_mu k_rho delta_nu,sigma - k_mu k_sigma delta_nu,rho
     -k_nu k_rho delta_mu,sigma + k_nu k_sigma delta_mu,rho]
       / Q_4²,                                               (2)

where k=kappa. Its second term is the same numerator times -1/(4W), a local
differentiated contact. Thus the nonlocal Maxwell denominator comes from the
native occupation correlations. The antisymmetric decoder supplies its tensor
interpretation, but no inverse Laplacian or wave pole is put into that decoder.

## 3. Exact slice pole and the reflected Gram

First take the infinite slice-length limit at fixed a and fixed nonzero spatial
momentum. It is the ordinary limit of the finite periodic Green sum; there is
no zero-mode ambiguity in this sector. Finite-volume constant modes are not
silently used to prove positivity. Write

    k=(kappa_1,kappa_2,kappa_3), Q=|k|>0,
    P_T=I-kk^T/Q², C=i[k cross],
    E_a=(2/a) asinh(aQ/2),
    g_a(t)=a exp(-E_a |t|)/(2 sinh(aE_a)).

Then C*=C, C²=Q²P_T, CP_T=C, and rank P_T=2. The scalar slice covariance is
g_a. Set E_i=F_0i and B_i=(1/2)epsilon_ijk F_jk. E lives at half-integer slice
centres; B lives at integer slice centres. On separated slices the centred
time difference of g_a has multiplier

    -(2/a)sinh(aE_a/2)=-Q                                    (t>0).

The product of two time differences gives -Q²g_a away from coincidence.
Consequently the separated positive-separation covariance of the column (E,B)
is

    g_a(t) [ -Q² P_T    -Q C  ]
           [  Q C        Q² P_T ].                           (3)

The signs follow from differentiating with respect to the first or the second
time argument. They also show why positivity of four scalar copies alone is
not the required vector-reflection argument.

Reflect the slice coordinate and the two-form indices: E is odd and B is even.
For fields on strictly positive slices, the reflected two-point kernel is

    g_a(t+s) [ Q² P_T    -Q C  ]
             [ -Q C       Q² P_T ].                          (4)

Restrict the test support to slice centres at least 2a from the reflection
plane. Then every term of the differentiated contact in (1) has zero
cross-reflection pairing. This support restriction states exactly where (4)
holds for the native occupation covariance; it is not a claim of positivity
for tests crossing the reflection plane or touching its contact layer.

For finitely supported complex test sequences e(t),b(t), with the appropriate
integer and half-integer centres, put

    e_tilde=sum_t exp(-E_a t)e(t),
    b_tilde=sum_t exp(-E_a t)b(t).

Multiplying out (4) gives the reflection quadratic form

    a Q²/(2 sinh(aE_a))
       ||P_T e_tilde - (C/Q)b_tilde||² >= 0.                 (5)

Both null longitudinal directions and two further redundant transverse
combinations disappear in this quotient. Its nonzero one-particle fibre is
exactly two-dimensional for every Q>0: the expression in the norm is an
arbitrary transverse vector, and P_T has rank two. Translating the positive
test support by one slice multiplies its quotient representative by
exp(-aE_a). Thus the reconstructed positive energy is E_a, the same source
pole as the massless band in row30.

The spatial zero mode of F has only contact terms and has no positive-energy
fibre in this calculation. Harmonic constant electromagnetic backgrounds are
not produced by a curl and are not included in this vacuum sector.

## 4. Gaussian reflection positivity and free Maxwell dynamics

The occupation theorem supplies convergence of all fixed moments, so the
limiting F variables are jointly Gaussian. Reflection preserves their law when
channel labels transform as a vector and F as a two-form. Gaussian reflection
positivity on polynomial functionals follows from (5): after Wick ordering,
cross-reflection expectations of degree r are sums of products of the positive
kernel (4). They are the symmetric tensor powers of its positive one-particle
form. Different Wick degrees are orthogonal. This proves positivity on the
finite polynomial algebra of the stated positive-support tests.

The full vector potential with vector reflection is not positive. Its A_0
sector has a negative reflected scalar norm. The physical positive algebra
used here is generated by F and quotiented by the null vectors of (5). No
claim of a four-polarization physical field is made.

Along a->0, kappa->p and

    E_a -> |p|,
    g_a(t) -> exp(-|p||t|)/(2|p|),
    Q² a/(2sinh(aE_a)) -> |p|/2.

The continuum reflected form is therefore

    (|p|/2) ||P_T e_tilde - (i[p cross]/|p|)b_tilde||².       (6)

The two transverse positive-energy states and the Gaussian tensor powers give
the free electromagnetic field-strength representation. The Bianchi identity
dF=0 is an exact decoder identity. To check the second vacuum Maxwell identity
without substituting a pole into its own denominator, choose a basis
epsilon_1,epsilon_2 perpendicular to p. The two quotient states have energy
q=|p|. After the usual continuation of the odd Euclidean electric component,
their mode amplitudes can be written

    E_phys=i q epsilon_r exp(-i q t+i p.x),
    B_phys=i(p cross epsilon_r) exp(-i q t+i p.x).

These exhaust the two-dimensional quotient and satisfy div E=div B=0,
partial_t B=-curl E and partial_t E=curl B by direct differentiation. Their
relative normalization is the one fixed by (6). Thus continuing the positive
slice generator yields the two transverse vacuum waves with dispersion
omega²=|p|². This is a consequence of the reconstructed correlations; no
real-time Maxwell step was executed.

## 5. Common slope and a controlled limiting order

The exact pole is the lambda=0 case of the row30 native-band formula. Hence

    0 <= |p|²-E_a²
      <= (a²/12)[sum_j p_j^4+|p|^4],
    |grad_p E_a|<=1.

The same routed source has the massive internal bands of row30. Their common
infrared slope tends to one with deficit at most a²lambda/6. No independent
Maxwell speed coefficient is selected to match them. This statement compares
the same coordinate rulers and the same slice generator; it is not an
attachment to laboratory seconds or metres.

A safe limiting order is: fix a finite source, take the record-block Gaussian
limit n->infinity, take volume to infinity on finitely many smooth smeared
observables with nonzero spatial momentum, then refine a->0. A diagonal choice
can approximate every finite list of moments, but no uniform finite-resource
central-limit rate is asserted. The contact coefficient 1/(4W) vanishes as
volume grows. Even before that limit, its cross-reflection term vanishes for
the strict-support domain used in (5). Standard smooth test smearing controls
the ultraviolet derivative distributions; pointwise continuum F variance is
not finite and is not claimed.

## 6. Law dependence and the remaining operational boundary

The readout takes only nearest-neighbor differences. It cannot supply a pole
absent from its input record covariance. For example, an independently
resampled-location control has P_control=11^T/M and occupation covariance
Pi/(4W) in the same normalization. Four such record channels give only local
contact terms after the same curl decoder, hence zero cross-reflection kernel
on the separated-support domain and no nonzero transverse propagating fibre.
This exact changed-law prediction is evaluated as a population calculation
in the companion experiment; no actual control tape is claimed.

This theorem removes the need to assume a Gaussian higher-function rule or to
install an independent oscillator generator for this free-field reconstruction.
It retains declared four-dimensional routing, four channels, their vector
interpretation, the curl readout, uniform preparations, long-block limits and
the Euclidean slice interface. It does not establish a same-history
intervention/readout correspondence between native events and the reconstructed
real-time field, a physical clock, charged matter coupling, alpha, or atomic
levels. It does not prove the source architecture is unique, and uniqueness is
not imposed as a new acceptance requirement.

For OpenWave this is a conditional positive mathematical component of row29.
The companion `experiment/` now retains a preregistered finite native branch
census, four actual current tapes and their finite-block population covariance
calculation. The changed law is checked at the population level; no changed-law
native tape is asserted. No owner pin choice has been inferred, and no
row-level gate or official cell is promoted by this theorem and finite check.

## 7. An optional native contact completion

The shared occupation theorem sections6–8 gives a stronger observation
protocol: independently record native coins on a disclosed zero-load pair,
form xi=2coin-1 and add Pi sum xi/sqrt(4Wn) to each occupation channel. Its
covariance is exactly Pi/(4W) for every n; its Gaussian limit cancels the
contact term in (1). These are native event records, not supplied Gaussian
draws. Their scheduler, normalization and extra event cost are explicit
additional apparatus. The same Maxwell quotient then uses an exact limiting
Green covariance. The present companion experiment uses occupation records
without this addition and does not borrow an execution of the auxiliary coins.
