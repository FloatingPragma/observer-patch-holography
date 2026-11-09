# One-phase charged-lepton sanity test

This test chooses the coherent family phase using the measured muon/electron
ratio. It tests the tau ratio without using the tau measurement in that
choice. The source event weights, balanced modulus, Higgs tensor and kinetic
convention are unchanged. This is an explicitly calibrated test, not a
source derivation of the phase or an identification of laboratory particles.

The comparison values are the
[2022 CODATA extensive listing](https://physics.nist.gov/cuu/pdf/all.pdf):
`m_mu/m_e=206.7682827(46)`, `m_tau/m_e=3477.23(23)` and
`m_e=0.51099895069(16) MeV`. The reference file records the role of each
input. Its quoted standard uncertainties are comparison intervals; no joint
statistical significance or theory-error model is claimed.

## Phase choice and exact source calculation

Write `R=m_mu/m_e`, `u=(sqrt(R)-1)/(sqrt(R)+1)` and

\[
 t=\tan(\theta/2)
   ={u(\sqrt2-1)\over\sqrt3+\sqrt{3-u^2}}.
\]

This solves the ratio of the two smaller positive circulant roots. The
code chooses a rational approximation to `t` with denominator `10^18`,
then defines the exact unit phase
`z=(1-t^2+2it)/(1+t^2)`. Independent rational radical enclosures verify
the resulting coefficients and ratios. The selected specification is

\[
 t={55785269584552807\over500000000000000000},\qquad
 \theta\simeq0.22222204707308785.
\]

For `z=x+iy`, the three ordered roots of the same source `C_ex` are

\[
 \lambda_e={\sqrt5\over6}-{\sqrt{10}\over12}x-{\sqrt{30}\over12}y,
 \quad
 \lambda_\mu={\sqrt5\over6}-{\sqrt{10}\over12}x+{\sqrt{30}\over12}y,
 \quad
 \lambda_\tau={\sqrt5\over6}+{\sqrt{10}\over6}x.
\]

The subscripts describe the calibrated comparison roles. The exact action
coefficients are `lambda_j^2/3`, with positive roots and `Q=2/3`.
The actual native charged-Higgs kernel and complete 180-dimensional free
field symbol are evaluated with this same `Y`; no eigenvalue is replaced
after the phase choice.

| Quantity | Calculated action ratio | Measured reference | Role |
| --- | ---: | ---: | --- |
| muon/electron | 206.768282700000011 | 206.7682827 ± 0.0000046 | Phase calibration |
| tau/electron | 3477.441636283848796 | 3477.23 ± 0.23 | Conditional comparison |

The held-out tau ratio lies inside the reported one-standard-uncertainty
interval. The values do not coincide at their printed central coordinates.
Using the electron mass as one additional common energy calibration gives
`m_tau=1776.9690272267634 MeV`, compared with CODATA's `1776.86±0.12 MeV`.
The scale is supplied by that electron input. The code also propagates the
muon/electron input interval through the conditional Koide tau formula.

## What the test establishes

The chosen phase is compatible with the observed ratios at the stated
comparison precision inside the balanced quadratic action. The
[full-CAR construction](CAR_KOIDE.md) realizes its family coefficient, and
the [free field theorem](RELATIVISTIC_MASS.md) identifies the corresponding
family-resolved Dirac branch masses in the declared reference and clock.
The phase score in [PHASE_AND_PREMISES.md](PHASE_AND_PREMISES.md) gives an
independent operational inverse for a readable coherent phase.

The calibrated match is a consequence of the Koide relation plus one mass
ratio input. It is not independent evidence that the repair law selects
the measured phase. The uncalibrated zero-phase program has two equal lower
coefficients. Changing the phase is explicit configuration, and the four
stationary-feedback policies in the parent package retain their different
coefficients.

The test does not establish an interacting chiral field limit, the selected
Higgs vacuum, measured weak decay widths, physical charged-state dressing,
or control of corrections between action and measured masses. In particular,
the [mass-response bound](MASS_RESPONSE_BRIDGE.md) requires actual bounds
on kinetic and self-energy corrections from that same interacting theory.
No such correction is chosen here to improve the comparison. The physical
particle map is empty.
