# Spectral thresholds on complete source records

A stable particle can have an isolated energy--momentum branch. That is
not the only physical identification route. For a charged sector with a
long-range electromagnetic field, a sharp mass eigenvector can fail under
the hypotheses of [Buchholz's Gauss-law theorem](https://doi.org/10.1016/0370-2693(86)91110-X).
A controlled infrared threshold and suitable inclusive or dressed
asymptotic detectors must therefore remain possible targets. This paper is
motivation, not a claim that the OPH spin charge is electromagnetic charge
or satisfies the theorem's premises.

## An original-weight edge theorem

Let nu be the subprobability spectral measure of an actual prepared source
response in a declared coordinate, theta=inf supp(nu). The missing weight
is retained in the complete preparation/filter/measurement record. Assume
the source supplies constants a>0, alpha>0, epsilon_0>0 such that

    nu([theta,theta+epsilon]) >= a epsilon^alpha
                 for every 0<epsilon<=epsilon_0.                 (1)

For x_-<x_+, suppose certified original CDF intervals give

    nu((-infinity,x_-]) <= u_- < a epsilon^alpha,
    nu((-infinity,x_+]) >= l_+ > 0.                               (2)

Then

    x_- - epsilon < theta <= x_+.                                (3)

If theta<=x_--epsilon, the whole interval in (1) lies below x_-,
contradicting (2). Positive original weight below x_+ proves the other
inequality. No atom or conditional renormalization is used. Conversely,
an empty observed bin alone gives no lower bound on the edge: arbitrarily
small unseen weight can lie below it. The modulus (1) is an independent
source premise, not something inferred from a finite count of zero.

The code uses positive integer alpha for exact rational arithmetic. The
proof also holds for positive real alpha. Detector resolution, source
preparation, finite-clock alias control and complete-channel error must all
be included in the CDF intervals. The bounded resolvent hierarchy in
`../spin_dynamical_memory/SPECTRAL_BRANCH.md` supplies another conditional
interface: Lipschitz inner/outer ramps bound interval probabilities, with
their boundary strips retained. It does not evaluate (1) for a continuum
source.

For N independent complete trials and M predeclared CDF events, the sample
means have variance at most 1/(4N). Chebyshev's inequality and the union
bound give simultaneous additive error eta with failure probability at
most M/(4N eta^2). A complete-output half-trace error d adds d to every
event interval. `sampled_cdf` retains failed preparations and successes
above the last threshold, and rejects insufficient N. Selecting thresholds
after seeing the same data needs a separate uniform bound or fresh data.

The analytic control nu((-infinity,x])=(1/5)clamp(x-9/16,0,1)^2 has no atom.
Taking delta=2^-k, x_-=9/16-delta, x_+=9/16+delta and eta=delta^2/100 gives
the certified bracket (9/16-2delta,9/16+delta]. Its detected near-edge
weight tends to zero while its edge converges. The complete failure weight
is 4/5. This is an implication control, not a computed OPH field state.

## Actual finite-source control

`threshold.py` replays the actual two-site spin source and complete
preparation/filter instrument in `../spin_dynamical_memory/RESOLVENT_CLOSURE.md`.
Its positive-clock response has energies 9/8 and 5/4 and original weights
81/400 and 16/121. All preparation and filter failures remain. The exact
first atom verifies (1) with a=81/400, alpha=1, epsilon_0=1 for this finite
response. With epsilon=1/1000, eta=1/10000 and probes 281/250 and 563/500,
(3) gives (1123/1000,563/500], containing the actual edge 9/8. This is a
source energy threshold; no continuum spatial momentum, rest mass or
additional large native experiment is claimed.

For a physical mass readout, a source-derived relativistic dictionary must
first identify the joint energy--momentum coordinate and common units. A
Lorentz-invariant lower boundary would then have E^2-c^2|P|^2=m^2 c^4.
Writing that equation does not establish it for the source. Spin, gauge
charge, statistics, infrared dressing and appropriate scattering or
inclusive detector predictions remain part of the identification.

## A restriction on identifying the whole spin-charge sector

The existing ordered-zero-field source theorem constructs compact charged
creators with energy transfer in [ar/2,2br], fixed positive original weight
p_* and arbitrarily small r>0. Here a,b and p_* are the actual analytic
constants from `spin_multipacket_comparison.bounds`, including its filter
L1 normalization. The charges are the source's global internal U(1) labels.

Suppose a proposed target has strictly positive energy threshold Delta in
the same charge sector, preserves the vacuum zero and uses one fixed
positive source-to-target energy conversion kappa. Choose r so that
2 kappa b r<Delta. The source event then has original probability at least
p_* and the corresponding target event has probability zero. Every
complete-output identification has half-trace error at least p_*.

Thus this entire gapless spin-charge sector cannot converge to that gapped
target sector under a charge- and energy-preserving dictionary. The code
computes r and this explicit lower bound from the existing source producer.
It does not identify global spin charge with electromagnetic charge, rule
out a massive subbranch, or exclude another source family. A catalogue must
retain those distinctions instead of attaching all integer source charges
to massive measured particles.
