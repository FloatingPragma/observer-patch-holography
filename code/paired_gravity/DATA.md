# Measurement provenance and comparison contract

All results were inspected on 2026-10-03. This is a retrospective comparison
of published compressed measurements, with no untouched dataset or new
prediction freeze. The negative source-selection theorem is independent of
these observations. `measurements.json` transcribes numerical facts; copyrighted
article text and PDFs are not redistributed. File hashes identify downloaded
primary versions, not independently acquired instrument records.

The static redshift input is [Pound and Snider, Physical Review 140 B788
(1965)](https://doi.org/10.1103/PhysRev.140.B788). Their nuclear-resonance
reversal observable is normalized by measured `2gh/c^2=4.905e-15`.
The published ratio is `0.9990`, count-statistical standard deviation `0.0076`,
and estimated systematic limit `0.010`. We retain the last as a bounded
nuisance, not another independent Gaussian variance. The supplied same local
metric coupling predicts unit ratio for every constitutive member to the
stated weak-field/laboratory order. Neither g nor the nuclear transition is
derived from OPH. Original count records and their likelihood are not replayed.

The deflection input is [Fomalont et al., Astrophysical Journal 699 1395
(2009), arXiv v1](https://arxiv.org/pdf/0904.3992v1). Table 5 contains four
overlapping reductions of the October 2005 VLBA differential-angle
measurements. All four are retained. The corona-corrected 43 GHz row is the
designated comparison; alternative frequencies and exclusions are sensitivity
controls, not independent experiments. Section 4 uses residual template
`P=(gamma-1)D/2`, where D is the GR angular prediction. Thus the measured total
bending amplitude is `A=1+(gamma-1)/2`, with standard deviation `sigma_gamma/2`.
The original reduction imports source positions, solar ephemerides, station
geometry and propagation corrections. We retain these as empirical/reduction
inputs. The compressed fits provide no cross-row covariance or complete raw
likelihood, so we form no combined chi-square, discovery significance or
probability for core OPH.

The unchanged universal parameter is gamma. The terrestrial and solar channels
use independently measured g,h and solar strength/geometry respectively, with
the same G,c and matter/photon coupling. There is no redshift-specific or
deflection-specific rescaling. The GR reference has unit redshift ratio and
unit total bending amplitude under those same measured inputs.

The interval comparison uses the published one-standard-deviation range,
with no newly optimized numerical decision threshold. A bounded systematic
can shift the redshift residual within its stated limit. All four source
members have the same redshift result, while their leading bending amplitudes
are 0, 1/2, 1 and 3/2. Report membership in the published interval, not a
Gaussian tail probability extrapolated thousands of standard deviations.
This is a first-order-template comparison. The separate analytic nonlinear
ray bound applies to infinity-to-infinity rays; it is not a precision error
budget for the unreplayed finite-observer radio reduction. Accordingly the
packet claims the constructive negative exit, not a positive end-to-end
physical postdiction of the microscopic OPH source.

The Galileo redshift paper was also inspected and its version hash retained.
It is not scored: its clock reduction uses orbital and signal-propagation
models whose changes cannot be ignored when changing spatial response. This
exclusion is based on the observable map, not on the reported agreement.

Source numerical rounding and discretization are separately bounded in the
derivation. The coarse source run is not used to assign astronomical theory
uncertainties. Physical source-to-mass, clock and detector identification stay
with their existing parents. Neither nominal agreement nor this comparison
selects the source law independently of the observations.
