# The collar survival coefficient under the presence reading: λ = 1 − P/24

**Status:** theorem-grade correction to the L5 occupancy clause.
(six theorem statements, built against the v10 tree).

**Date:** 2026-07-14.

## Finding

The canonical value χ_ν^can = e^(−P/24) = 0.9343006394893864… rests on a type
substitution between two different quantities that share the symbol ε(y):

1. The screen microphysics paper defines ε(y) as a ratio of quotient traces of
   projections (`screen_microphysics_and_observer_synchronization.tex`,
   lines 1102–1123). Z₆ is a sum of sector projections (lines 622–649), so
   τ_q(Z₆) = P/24 is the **probability of the reserve-presence event**, and
   ε(y) is its slice-conditional probability. The receipt
   `Z6_NORMALIZED_TRACE_P_OVER_24` computes exactly this.
2. Assumption *Local Poisson reserve survival* (lines 1260–1271) consumes the
   same number as the **mean of an ℕ-valued occupancy count**.

These readings are mutually exclusive at the stated precision. If the
reserve-presence event is Z₆. The *No separate scalar carrier*
(lines 742–789) forces, since an occupancy variable not reducible to the
edge-center sectors would be a second scalar carrier. In that case survival of a
scalar opportunity is the complement event, and the per-slice survival factor
is 1 − ε(y) for **every** occupancy law. No distribution placed between the
presence event and its complement can change the identity
Pr[absent] = 1 − Pr[present].

## Consequences

**1. `presence_gate`.** Under the presence reading, the scalar-weighted mean
receipt Σ w(y)ε(y) = P/24 alone forces

    λ_collar = 1 − P/24 = 0.9320429912748350…

exactly, by linearity of the trace. The uniformity clause (L7, slice-wise
scalar-reserve unbiasedness) is not consumed: linear functions have no Jensen
gap. The corrected exact value therefore needs strictly fewer gate clauses
than `UNIFORM_PRODUCT_THICKENING_EXACT`; L5 and L7 both drop out.

**2. `presence_below_poisson_floor`.** 1 − P/24 < e^(−P/24). The claimed
theorem-level band [e^(−P/24), 1] **excludes** the value implied by the
corpus's own receipt semantics. The band is an artifact of the Poisson
choice, not a floor.

**3. `markov_band`.** Under the mean-count reading (which severs the link to
the projection trace and requires a carrier beyond Z₆), Markov's inequality
alone gives the assumption-minimal band

    1 − P/24 ≤ λ_collar ≤ 1

for any ℕ-valued occupancy. The qualitative Tier-C claim states that χ_ν^can is order
one and bounded away from zero. It survives with strictly weaker hypotheses
than L5. Only the 16-digit exact value is occupancy-law-dependent.

**4. `finite_refinement_below_poisson`.** Every finite m-fold sub-slot
refinement (independent Bernoulli(ε/m) occupancy, survival = all-empty) gives
(1 − ε/m)^m < e^(−ε). The Poisson value is the m → ∞ supremum of the
finite-carrier family and is attained at no finite regulator. A fixed-cutoff
finite carrier cannot realize e^(−P/24); reaching it requires leaving the
declared finite framework.

**5. MaxEnt does not rescue L5.** The L5 clause names "a finite
MaxEnt/occupancy theorem" as the missing derivation. MaxEnt on ℕ with a mean
constraint and counting base measure yields the geometric law, with survival
1/(1+ε) = 0.9363672805…, a third value distinct from e^(−ε). Poisson is MaxEnt only
under the 1/k! base measure, which encodes the independent-rare-placements
structure it was supposed to derive. The three candidate laws give three
distinct coefficients:

| occupancy law | survival at ε = P/24 | status |
|---|---|---|
| Bernoulli / presence (literal register) | 1 − P/24 = 0.9320429913 | forced by the projection semantics |
| Poisson (corpus L5) | e^(−P/24) = 0.9343006395 | declared; unreachable at finite regulator |
| geometric (MaxEnt on ℕ) | 1/(1+P/24) = 0.9363672805 | what the named repair route actually yields |

## Downstream impact

- Force coefficient for the 80×60 mm PoC: 5.130×10⁸ ΔS N (was 5.142×10⁸),
  a −0.24 % shift. No experimental target moves at a precision that matters;
  ΔS remains unknown across orders of magnitude.
- The 16-digit statement of χ_ν^can should be demoted from theorem-grade to
  chart convention, with the receipt-consistent exact value 1 − P/24.
- The numerical proximity of e^(−P/24) = 0.9343 to the binned-RAR profile
  value 0.9367 (used as motivation in the profile-envelope remark of
  `chi_nu_susceptibility_bounds.tex`) loses force under the correction; note
  that the geometric value 0.9364 sits closer to the RAR number than the
  Poisson value does, which illustrates how little selective proximity
  licenses.
- The galaxy-side dark-sector law is untouched: there the Poisson structure
  arises from a rare-events binomial limit over many repair opportunities
  with unbounded mean μ(x) = λ_c√x, a genuinely different object from the
  lab collar coefficient, whose mean is pinned to a projection trace.

## What decides between the readings

The two readings are separated by one question: is P/24 the probability of a
quotient-visible event (the trace of the Z₆ projection, as every receipt in
the corpus computes it), or the mean of a hidden multiplicity? The first is
the corpus's own formal semantics and yields 1 − P/24. The second requires a
new carrier that *No separate scalar carrier* forbids, and even then yields
only the Markov band, inside which Poisson is one unforced point. Either the
receipt or L5 must change; they cannot both stand.

## Scope

Nothing here settles whether any real substrate produces a nonzero ΔS
receipt, the G9 record-to-gravity calibration, or the nature-instantiation
obligation. The correction is internal to the theory side: it replaces the
declared L5 clause with the value its own register forces, strengthens the
nonzero-coefficient claim to a Markov-grade theorem, and retires the
16-digit Poisson value as unlicensed at any finite regulator.
