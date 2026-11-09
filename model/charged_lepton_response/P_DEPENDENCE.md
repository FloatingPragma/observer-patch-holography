# Capacity dependence of the declared charged response

The MCPR source has an explicit capacity-dependent charged coefficient.
Its complete native C-duty inverse realizes that coefficient in a specified
correlated weak-pulse limit. This is a positive construction with readable
parameters. The finite repair assumptions alone do not select its response
law, scale ansatz or source coordinates. The balanced coherent analyzer
and the MCPR construction are different configured actions.

`p_dependence.py` authenticates the complete source-binding maps of
[`mcpr_native_memory_bridge`](../mcpr_native_memory_bridge/PROOF.md) and
[`pixel_mass_closure`](../pixel_mass_closure/PROOF.md), then evaluates their
retained source coordinates with rational outward enclosures. It reruns no
root search or native-memory producer and fits no parameter. Experimental
values enter only a separate comparison, after computing each source shape.
The reference is [calibration_reference.json](calibration_reference.json).

## The source law and its capacity variable

Write `u=alpha_U` and let P be the dimensionless pixel capacity of the
declared response. With phi=(1+sqrt(5))/2, the law is

```math
\chi={u^2/512\over1-u/77},\qquad
\zeta={(P-\phi)^2/(21\pi)\over1-u/27-u^2/135},
\quad\delta=2/9+\zeta,\quad\rho=\sqrt2e^{-\chi},
\quad r_j=1+\rho\cos(\delta+2\pi j/3),\quad m_j=v A(P,u)r_j^2.
```

Here A(P,u)>0 is the common determinant normalization in the owning law;
it cancels from all ratios. The ordered positive roots are explicitly
certified. The trigonometric sum identities give the exact quotient

```math
Q={\sum r_j^2\over(\sum r_j)^2}
  ={1+e^{-2\chi}\over3}<2/3\qquad(u>0).
```

Thus using only the phase delta in an exactly balanced coherent analyzer
does not realize this mass law: its amplitude rho changes too. The full
native inverse instead chooses two admissible pairs of probe rates and a
gain g satisfying `spec(g Z_C)={m_j/v}` with all disabled outcomes retained.
The family change of basis is gauge-trivial. Its existence does not select
the rates, gain, clock or capacity law from the unconfigured repair process.

P is neither the number of weak kicks nor the native FIFO/attempt limit.
The closure design imposes `P=phi+sqrt(pi) alpha` as one of its response
relations. No identity identifies that real coordinate with a total number
of physical cells or a hardware memory length. A law connecting those
quantities must specify the same source, reference, clock and readout.

## Four authenticated branches

The first two coordinate pairs are exact declared decimals from the owning
forward and measured-endpoint artifacts. The forward coordinate lacks a
complete physical endpoint-selection proof. The measured-endpoint branch
inherits measured-alpha input. Neither branch reads charged-lepton masses
when this certificate executes. The MCPR architecture itself has the
target-informed provenance stated in its owner.

The other two rows are interval-certified roots of declared response
designs, each unique within its stated interval. The supersymmetric proxy
uses beta=4 with additional fields absent from the native quadratic
carrier. The minimal census uses a chosen color-scale coefficient beta=7.
Neither choice is forced by finite current-representation arithmetic.

| Branch | P | u | delta | rho |
|---|---:|---:|---:|---:|
| Forward coordinates | 1.630972095858897 | 0.041124247441817 | 0.222224763427524 | 1.414208888550144 |
| Measured-alpha endpoint | 1.630968209403959 | 0.041124336195630 | 0.222224761901066 | 1.414208888529965 |
| Supersymmetric proxy, beta=4 | 1.630971514159316 | 0.041124260725897 | 0.222224763199025 | 1.414208888547124 |
| Minimal census, beta=7 | 1.628533369611833 | 0.019531387695090 | 0.222223894366381 | 1.414212508420078 |

The table displays rounded centers. The certificate retains exact rational
enclosures for every coordinate, response and ratio, including the full
closure-root boxes rather than only their centers.

| Branch | m2/m1 | m3/m1 | Q | Individual one-standard-uncertainty intervals containing the whole source interval |
|---|---:|---:|---:|---|
| Forward coordinates | 206.768300595181 | 3477.365536596023 | 0.666664463409039 | tau/e only |
| Measured-alpha endpoint | 206.768282739123 | 3477.365262268937 | 0.666664463399526 | mu/e and tau/e |
| Supersymmetric proxy, beta=4 | 206.768297922262 | 3477.365495531309 | 0.666664463407615 | tau/e only |
| Minimal census, beta=7 | 206.782621310399 | 3477.644439943646 | 0.666666169828635 | neither |

The frozen CODATA2022 comparison coordinates are
`mu/e=206.7682827±0.0000046` and `tau/e=3477.23±0.23`.
The mu/e residuals divided by the quoted uncertainties are approximately
3.8903, 0.00851, 3.3092 and 3117.09, respectively. The tau/e values are
0.5893, 0.5881, 0.5891 and 1.8019. These are descriptive comparisons of
conditional action ratios, without a covariance-aware significance test or
a bound on physical theory errors. The second row's agreement is not an
independent prospective prediction of a historically target-informed model.

The closure receipts retain `v/E_*` approximately `2.01983668924e-17` and
`8.61718526720e-21`. E_* is a declared energy unit with no laboratory value.
Supplying the measured electron value as a common unit would give tau
values approximately 1776.93014, 1776.93000, 1776.93012 and 1777.07266 MeV;
that additional calibration is not a source prediction. These optional
numbers are direct rescalings of the table, not inputs to its calculation.

## Scope and independent checks

The evaluator uses the owning rational pi, square-root and Taylor bounds
for exponential and cosine. Every denominator, positive root and strict
ordering is admitted explicitly. The two endpoint shapes are checked
against the owner's independent rational shape helper. Both closure shapes
are checked against their retained mass-coefficient intervals. The tests
also evaluate the circulant Fourier spectrum independently at high precision,
reject incomplete provenance, and alter the experimental reference while
checking that source coordinates and ratios cannot change.

This calculation supplies the actual conditional spectra of the four
specified C-duty realizations. It identifies no natural particle, no
physical defect, and no interacting mass correction. The same state,
physical clock, observable and renormalized charged response are required
to turn a configured family spectrum into an experimentally identified
electron threshold or muon/tau resonance. A capacity-dependent action law
is available; its selection by the primitive dynamics is a separate claim.
