# Native leptonic widths in GeV

## Source and normalization

The native charged analyzer and the electron anchor define the three charged
energies. The calibrated Dirac tensor of `openwave_neutrino_calibration`
defines the three neutral masses and its unitary family frame U. The six
oscillation inputs, chosen lightest neutral mass, charged analyzer phase and
electron anchor keep their stated ancestry.

This member has the same positive-frequency current construction and smooth
ultraviolet amplitude as `openwave_inclusive_decay`. Each heavier charged
parent couples to every allowed lighter charged daughter, all neutral mass
pairs, spin/helicity labels and momenta, with the Hermitian reverse vertex.
For a parent of family beta and daughter alpha, the current coefficient is
U_(beta,i)* U_(alpha,j). These are actual coefficients of the native weak
raising matrix in the declared mass frame. The complete angular dark space
and all free spectator sectors are retained. For the tau-like parent both
charged daughter continua occur with the same coupling and parent amplitude.
Each prepared-parent member is a declared source. The tau member omits a
subsequent muon-daughter decay vertex; the separate muon member calculates
that parent process. No simultaneous hierarchy or full cascade theorem is
transferred from these reducing sectors.

The conserved number combinations prove an exact parent-plus-three-daughter
sector for each admitted parent input. The positive-frequency vertex is an
explicit architecture; it is not asserted equivalent to the unfiltered
native contact source. No selected two-state approximation is used.

The native internal factor is one. Direct Pauli traces give the unpolarized
amplitude square

    8 G² (P.p_antineutral)(p_charged.p_neutral).

The trace uses the massive left spin sums D(p)=p0 I-sigma.p. It sums both
massive-neutral helicities. Neutrino and antineutrino are distinct Dirac
outputs, so there is no identical-particle factor. The coefficient is
recomputed by the imported exact phase-space and native-current checks.

Each finite-box field has u(p)/sqrt(2E V). Integrating the quartic vertex
leaves 1/V and a momentum Kronecker delta. The two independent momentum sums
have density V²/(2pi)^6. The parent flux is 1/(2M). With
s=(p_charged+p_neutral)², masses m,a,b for daughter, neutrino and antineutrino,
the complete three-body width is

    G²/(128 pi³ M³) integral ds
      (M²+b²-s)(s-m²-a²)
      sqrt(lambda(M²,s,b²)) sqrt(lambda(s,m²,a²))/s,

with (m+a)²<=s<=(M-b)². The scalar integral is sufficient because the spin
sum depends only on these invariants. All nine mass pairs have the weights
|U_(beta,i)|² |U_(alpha,j)|². Unitarity makes their total one, but their
separate massive phase spaces are retained. In the massless-neutral limit,
the exact integral gives

    Gamma = G² M⁵ f(m/M)/(1536 pi³),
    f(r)=1-8r²+8r⁶-r⁸-12r⁴ log(r²).

The native spin, Fourier and phase-space factors reconstruct 1536; it is
not inferred from a measured lifetime.

The conventional four-component current with (1-gamma5) equals twice the
left-Weyl current. Thus an empirical coefficient GF/sqrt(2) multiplying
two such currents corresponds to G=2sqrt(2) GF in this source. G has units
GeV^-2; Gamma has units GeV. In parent-energy units the dimensionless
coupling is G M². The exact SI Planck constant and elementary charge give
hbar in GeV seconds. This is a physical unit dictionary, not an emergent
observer-clock theorem.

The named PDG 2025 GF value is a weak normalization input with muon-lifetime
ancestry. The calculated muon characteristic time is consequently a
normalization check. A separate control uses alpha(0), W and Z masses in
the on-shell tree relation; it fails the precision muon comparison. No
radiative correction or loop-matched weak coupling is installed.

## Uniform massive-neutral correction

Scale the parent mass to one. Every admitted charged ratio is less than
1/16, and every neutral ratio is less than mu=1/(2*10^9). Let I0 be the
massless-neutral integral; f(r)>=1-8r²-r⁸ proves I0>2/25. The massive
integrand is nonnegative and at most its massless counterpart. For the
outer factor this follows by squaring

    (1+b²-s) sqrt[(1-s-b²)²-4sb²] <= (1-s)²;

the kinematic domain has 1-s>=2b-b². The inner factor is bounded by
(s-r²)² in the same way.

Set epsilon=1/20000 and retain the central interval
r²+epsilon<=s<=1-epsilon. The two omitted massless integrals are at most
epsilon²/2 and epsilon³/3: use (s-r²)²/s<=s-r² near the lower endpoint.
The shifted physical thresholds are inside those strips. On the central
interval the outer factor is at least its massless value times
1-4mu²/epsilon². The inner factor is at least its massless value times
1-5mu²/epsilon², using sqrt(1-x)>=1-x. Therefore

    0 <= (I0-I_massive)/I0
       <= (epsilon²/2+epsilon³/3)/(2/25)+9mu²/epsilon²
       < 1.7e-8.

The receipt checks a slightly sharper rational floor. Complete normalized
family weights preserve this bound. The displayed numerical massive
integrals are diagnostics inside this analytic rate enclosure; the
quadrature tolerance is not substituted for the proof.

## Finite survival comparison and its domain

For the massless-neutral comparison member, rho(E)=
E⁵ f(r/E)(2-E)^4/(3072pi^4), r<=E<=2. The parent survival obeys the exact
Volterra equation derived in `openwave_inclusive_decay`. For all three
charged channels, 1/4000<=r<=1/16. Termwise differentiation with
log(2/r)<10 gives the rational derivative envelopes in the producer.
Integration by parts proves

    integral |C(t)| dt < 1/1000,
    integral t |C(t)| dt < 1/125.

The survival probability error against exp(-Gamma t) is at most
2G² B1(1+G² B0 t), where G is in parent units. Rounding each reservoir energy
to a bin midpoint of width h gives the independent probability error t h.
The receipt evaluates one and four lifetimes with the supplied physical
weak normalization and h=2^-110. This large finite-bin existence statement
is analytic; it is not an executed enumeration of all bins. For the common
tau parent the two continuum correlation bounds add and its rate is the
sum of the two partial rates. All angular dark states remain in the source.

The explicit 64-bin computations retain the complete bright spectrum,
all emission-bin probabilities and survival, with actual physical weak
normalization, autonomous energy and readout work. They exhibit finite
unitary dynamics and do not establish irreversible exponential behavior
at arbitrarily long times. The Volterra bound applies directly to the stated massless-neutral
comparison source. The weighted spectral argument below supplies an explicit
transfer for massive parent survival and complete inclusive energy bins.

## Comparison roles

The ratio of the two tau-like leptonic widths is independent of GF and
the electron energy anchor. The native charged analyzer supplies the mass
ratio; the neutral family sum and its bounded mass correction are explicit.
The PDG 2025 HFLAV ratio is a comparison value that enters no source or
normalization parameter. The tree result agrees with that ratio at its
quoted precision. The muon width does not match the measured lifetime at
precision level; its electromagnetic radiative terms are absent.
Hadronic tau vertices are outside this declared source, so the summed
leptonic width is not a prediction for the physical total tau lifetime.

## Exact massive spectral transfer at weak lifetimes

A stronger weighted estimate transfers the survival theorem to the actual
calibrated massive channels. Put D=1-r² and mu=max(a,b), where all masses
are scaled by the variable daughter total energy. On the physical interval,
write A=1-s and B=s-r². The massive outer factor X and inner factor Y obey

    0<=A²-X<=b^4+8sb²,
    0<=B²-Y<=2a²(s+r²).

For the first bound, use sqrt(x²-y)>=x-y/x and
(A+b²)/(A-b²)<=2, since b<1/10 and A>=2b-b².
The same square-root inequality proves the second bound. Integration of
the integrand difference over the physical interval is bounded by

    b^4 D²/2 +(8/3)b²D³ +(4/3)a²D³.

The upper missing interval contributes at most (16/3)b³D², since
s>1/2 there. The lower interval has length 2ra+a²<=3ra and contributes
at most 9ra³D². The admitted masses satisfy a<=r; when phase space is
nonempty, b<=D and ra<=D. These inequalities imply the uniform estimate

    0<=I0-Iab<=21 mu²D³.

If the massive threshold exceeds the total energy, D<4mu and r>1/2;
I0<=D⁵/(30r²) gives the same estimate. Thus it covers the threshold strip
as well. Independently, s<=1 gives I0>=D⁵/30.

Let m denote the fixed charged daughter mass in parent units and let
mu_max denote the largest fixed neutral mass in those units. For the
common ultraviolet amplitude (2-E)², the exact coarea densities are
rho(E)=E⁵ I(E)(2-E)^4/(256pi^4). Therefore

    (sqrt(rho0)-sqrt(rhom))² <= (rho0-rhom)²/rho0
      <= 13230 mu_max^4 E(2-E)^4 D/(256pi^4).

Integrating, using D<=1, integral_0^2 E(2-E)^4 dE=32/15 and pi>3,
gives ||sqrt(rho0)-sqrt(rhom)||_2<=2 mu_max². Summing the nine normalized
family weights preserves this bound. Each parent-spin block has the same
bound by rotational covariance.

The exact angular coarea theorem represents the parent and all bright
energy channels on a common direct-integral carrier. On the interval
[m,m+a+b), the massive source has no daughter states; its comparison
representation adjoins uncoupled continuum channels of zero amplitude.
This is a survival-preserving dilation, not an onto identification of the
original free spectral fibres. Every original angular and helicity dark
complement is retained. This comparison does not identify
individual massive and massless angular records. The two full Friedrichs
operators differ by an offdiagonal operator of norm at most 2G mu_max².
Duhamel consequently bounds the complete unitary difference by
2Gt mu_max², and any common inclusive record probability by 4Gt mu_max².
This estimate uses a squared neutral-mass suppression, rather than a
bare mass-times-time estimate. It remains useful at the actual calibrated
weak lifetime. The receipt evaluates it with exact mass squares and rational
bounds for G and t. For a common tau parent with both daughter continua,
the squared coupling norms add; the norm is at most 2sqrt(2)G mu_max².

Combining this transfer, the Volterra bound and the energy-bin bound proves
the stated finite-time massive parent survival interval. The massive and
massless on-shell rates also differ relatively by at most
630mu_max²/(1-r²)². A finite-bin source approximates these inclusive massive
records with its full complement and readout work. The finite construction
has recurrence and supplies no exact all-time exponential law.
