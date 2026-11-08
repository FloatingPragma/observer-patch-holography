# Complete spectral-tube measurements on the selected spin source

The selected cubic spin source has a positive-energy continuum field
representation and nonzero low-energy spectral weight. Those statements do
not identify a sharp particle shell. This construction makes the missing
spectral-component test an explicit finite experiment with a resolution,
success probability and complete error budget. The source is the unchanged
spin Hamiltonian, ground preparation and controlled clock used by
`native_joint_spin_spectroscopy`. No measured mass supplies a coefficient.

## Complete energy and spatial observations

On a finite periodic even cube use the actual commuting combined translations
and their spectral projectors Q_k. On the open eight-site cube the same
construction uses combined graph reflections; its labels are not translation
momenta. For the selected finite ground let K=H_B-E_0 be its nonnegative
excitation generator. A source packet success map C, with every complementary
outcome retained, gives the positive subprobability joint measure

    mu(k,dE)=||Q_k dP_K(E) C Omega||^2,   s=mu(total)<=1.

Examples are a unitary local Pauli, or the existing energy-filtered transverse
packet with its explicit spatial and temporal linear-combination normalization.
The latter normalization must not be omitted when s is small.

For M clock values and tau>0 define the existing finite instrument

    B_j(K)=M^(-1) sum_(r=0..M-1) exp(2 pi i jr/M) exp(-i r tau K).

Every successful joint probability is

    p(k,j)=integral F_j(tau E/(2pi)) mu(k,dE),
    F_j(phi)=|M^(-1) sum_r exp(2pi i r(j/M-phi))|^2.

All packet failures, native processor failures and unused clock bins remain
outcomes. Summing every j gives one exactly, by finite Fourier orthogonality;
summing all successful joint outcomes gives s. No detector branch is
success-renormalized. The native controlled evolutions, group operations,
Fourier register and their complete dilation errors are the existing source
construction. This package evaluates an additional finite spectral certificate;
it does not assert execution of a larger native pulse tape.

The complete finite source spectrum must satisfy

    0<=E<=Emax,    tau Emax<2pi.                              (1)

An independent certified origin error enters the complete clock error. A
clock phase alone is periodic; deleting high energies or interpreting aliases
as low masses would violate (1). When the box grows, its full Emax may grow,
so tau must decrease or an independently justified spectral restriction must
retain its rejected outcomes. A small-momentum preparation alone is not a
bound on the entire energy spectrum.

## A summed clock-tail estimate

Use turns rather than radians. Let d be circular distance in [0,1/2]. Away
from a clock center,

    F_j(phi) <= 1/[M^2 sin(pi d)^2] <= 1/(4 M^2 d^2).

For any margin delta in (0,1/2], sum over all centers at distance at least
delta. On each side their distances have spacing 1/M. The decreasing-square
sum is bounded by delta^(-2)+M/delta. Including both sides therefore gives

    sum_(d>=delta) F_j(phi)
       <= ell(M,delta)
       := min(1, 1/[2(M delta)^2]+1/[2M delta]).              (2)

This bound includes every distant bin. It has no factor proportional to a
manually selected list of observed outcomes. For a requested ell<=epsilon<1,
a power-of-two M at least 1/(epsilon delta) suffices. The finite register
then uses log2(M) clock qubits. The maximum controlled time and the sum of
the binary controlled durations both equal (M-1)tau.

## Tube weights from complete counts

For each admitted k choose one proposed shell center E_c(k), with a certified
turn interval [a_k,b_k] containing tau E_c(k)/(2pi). The proposal may come
from a source calculation. An exploratory scan over speeds remains a scan;
it is not a source derivation of the successful speed. Choose a tube turn
width w and collar delta<w. The center, tube and collar lie strictly within
one phase period and w+delta<1/2.

Define inner clock bins by

    b_k-w+delta <= j/M <= a_k+w-delta,

and outer clock bins by

    a_k-w-delta <= j/M <= b_k+w+delta.

Let v_- and v_+ be their original, unnormalized joint probabilities, summed
over the specified momentum window. Use s for the success weight in that
same window. Every inner center lies at least delta from energies outside
the true tube. Every energy inside the true tube lies at least delta from
the complement of the outer centers. Applying (2) pointwise gives

    max(0,v_- - ell s) <= mu(tube)
                       <= min(s,v_+ + ell s).               (3)

This statement permits a continuum outside the tube and within its energy
window. It does not require the entire response to concentrate.

Suppose certified intervals for the complete measured probabilities are
[s_l,s_u], [v_l^-,v_u^-], [v_l^+,v_u^+]. A total output half-trace error
epsilon, valid for normalized complete instruments, gives

    S_u=min(1,s_u+epsilon),
    L=max(0,v_l^- - epsilon - ell S_u),
    U=min(S_u,v_u^+ + epsilon + ell S_u).                    (4)

Then L<=mu(tube)<=U. If S_u>0, L/S_u is a lower bound on the fraction of the
original successful response in the tube. A strictly positive L certifies
positive success without a conditional-state assumption. When rounded
numerical maps only have a gain bound, first supply their stated complete
event-error bound; an exactly trace-preserving shortcut is not available.

Ground preparation, field-to-zero-field clock comparison, native pulse
precision, splitting, Fourier synthesis and origin errors all contribute
to epsilon. The comparison must use the same prepared input and the actual
observed clock. A finite H_B ground cannot be replaced by the infinite H_J
ground because their clocks have nearby finite matrices.

For N independent normalized repetitions the three event means have
variance at most 1/(4N). A predeclared scan of R candidate shells therefore
has simultaneous sampling error eta, except with probability at most

    3R/(4N eta^2).                                        (5)

The implementation retains N-success nonselected or failed trials in the
counts. It rejects insufficient sample size. Equation (5) is a conservative
finite-shot bound, not a claim that repetitions were executed or that a
source-derived continuum theorem follows from a successful finite scan.

## What a shrinking-tube theorem would establish

Use one common source family, state prescription, translation dictionary,
clock normalization and speed c>0. Let p tend to zero on one subsequence
on which the normalized original broad-window joint spectral measures
converge weakly. Their momentum support lies in a fixed compact annulus
after scaling k=pq. The broad-window energy support is also compact after
scaling E=Jp lambda. The proposed center is c J s(k), with
s(k)=sqrt(2 sum_j(1-cos k_j)); s(pq)/p tends uniformly to |q| there.

Suppose the original response fraction in tubes of energy half-width
Jp delta_p is at least zeta>0, with delta_p tending to zero. For every
fixed closed tube about lambda=c|q|, these shrinking tubes are eventually
inside it. The closed-set inequality for weak convergence gives limiting
weight at least zeta in each fixed closed tube. Continuity from above
then gives at least zeta on the cone itself. Thus a positive sharp component
exists in that limiting response; a continuous complement is allowed.

Equations (2)--(5) provide a sufficient operational way to certify the
finite inequalities only if their bounds are uniform along the sequence.
For a collar of physical energy eta_p=Jp delta_p/2, its turn margin is
tau eta_p/(2pi). Vanishing clock leakage requires

    M tau Jp delta_p -> infinity.                          (6)

Consequently the required observation time grows faster than
1/(Jp delta_p). Fixed-time or fixed-history convergence alone supplies no
such comparison. Preparation and complete source-output errors must be
small relative to the original success probability. The inherited compact
packet has a conservative success lower bound proportional to p^4, with
its stated envelope and filter constants. An absolute error of fixed size
cannot resolve a fixed fraction of that response as p tends to zero.
Finite-shot errors need the same relative control. An effective
thermodynamic state approximation at these growing observation times is a
separate input; no such modulus is derived by this readout construction.

Neither one finite tube nor (6) proves the uniform positive fraction.
The selected-ground current/force estimate, or another source proof of a
positive sharp component, is the substantive missing input. Even a sharp
cone supplies neither Lorentz covariance of the full field algebra nor
interacting scattering, physical species or a mass calibration.

## A certified finite source readout

The positive control uses exactly the selected open eight-site source with
J=1, B=1/4, its complete 256-vector integer witness, and the local encoded-Z
probe. It retains the combined-reflection character 111 and all other
characters as distinct outcomes. Its energy clock has M=128, tau=1/2.
The source-chosen test band [4/5,7/5] surrounds its first excitation gap,
whose certified interval is supplied by the existing full matrix proof.
The band choice uses source spectral information; it is not an independent
experimental prediction.

The calculation bounds every clock kernel by exact rational cosine
intervals, combines integer overlap squares, and adds the inherited polar,
origin and ground-vector errors. Equation (4) certifies more than 1/5 of
the original probe probability in this band and character, without dividing
by the character probability. Every actual finite level in the band lies
in the first isolated excitation cluster; that follows from the full
256-level interval certificate, not from a few selected eigenvectors.

This is a finite spectral excitation of the selected source. The open-cube
character is not a small momentum. The 128-bin source probabilities are
certified spectral computations and a finite constructive instrument, not
a reported native-tape execution or a measurement of a continuum mass.
