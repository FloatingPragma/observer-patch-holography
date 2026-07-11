# Ultraluminous X-Ray Sources as Accretion-Record Normal Forms

## Motivating Result

On 8 October 2014, Matteo Bachetti and the NuSTAR team reported a clock inside
M82 X-2. It ticked every 1.37 seconds, and its arrival time moved through a
2.5-day orbit. The pulsed component alone radiated an isotropic-equivalent
$4.9\times10^{39}\ \mathrm{erg\,s^{-1}}$ between 3 and 30 keV; the associated
source had reached $1.8\times10^{40}\ \mathrm{erg\,s^{-1}}$ in the 0.3 to
10 keV band. The clock was a rotating neutron star
[Bachetti et al. (2014)](https://doi.org/10.1038/nature13791).

That observation broke the convenient habit of reading an X-ray luminosity as
a mass label. A source shining at an apparent hundred times the classical
neutron-star Eddington luminosity could, in fact, contain a neutron star. NGC
5907 ULX-1 later made the point at an apparent luminosity near
$10^{41}\ \mathrm{erg\,s^{-1}}$, while its spin changed from 1.43 to 1.13
seconds [Israel et al. (2017)](https://doi.org/10.1126/science.aai8635).

Issue
[#454](https://github.com/FloatingPragma/observer-patch-holography/issues/454)
asked OPH to classify ULXs as beamed neutron stars, super-Eddington
stellar-mass black holes, intermediate-mass black holes, or another branch.
The observations force a sharper answer. Those labels cut across different
physical coordinates and cannot serve as one categorical partition.

**Status:** conditional OPH record-classifier and evidence-audit specification.
The algebraic luminosity degeneracy is exact under the stated definitions, but
individual branch assignments remain conditional on the frozen conventional
forward models and each source’s observations. OPH supplies no native accretion
or ULX population source law.

**Date:** 2026-07-11.

## Standard physics, unresolved target, and OPH contribution

**Standard physics.** Compact-object accretion theory supplies the Eddington
scale, supercritical disks and winds, magnetic-column models, relativistic
timing scales, binary mass functions, radiative transfer, and detector response.
Pulsations, orbital dynamics, spectra, winds, nebulae, and survey selection are
ordinary astronomical evidence channels. Together these channels show that
ULXs are a luminosity-selected mixture rather than one engine class.

**Unresolved target.** Many individual sources lack a dynamical mass,
viewing geometry, bolometric luminosity, magnetic-field measurement, or unique
accretion-mode identification. Population fractions additionally require
exposure, contamination, censoring, orientation, and binary-population models.
This note does not solve those physical inverse problems or derive a ULX source
law.

**OPH contribution.** OPH specifies a candidate bounded self-reading software
classifier over those conventional models. Its local state records candidate
families and nuisance parameters; its ports ingest calibrated observations; its
readback exposes residuals, missing records, and viable-set size; and its
record-conditioned repair moves are restricted by frozen checkpoints. This is
an OPH-specific evidence architecture, not a unique physical explanation of
ULXs.

## Abstract

Ultraluminous X-ray sources form a luminosity-selected mixture. Their correct
classifier output is set-valued and factorized by accretor identity, mass,
intrinsic Eddington ratio, emission anisotropy, viewing angle, and source mode.
The isotropic-equivalent bolometric luminosity measures only $M\lambda/b$
for a fixed declared effective opacity, where $M$ is
compact-object mass, $\lambda=L_{\rm true}/L_{\rm Edd}$, and
$b=L_{\rm true}/L_{\rm iso}$ for a beam directed toward the observer. A
luminosity by itself therefore identifies none of the three quantities. At
$10^{40}\ \mathrm{erg\,s^{-1}}$, a beamed supercritical neutron star, a
supercritical stellar black hole, a sub-Eddington $10^3M_\odot$ black hole,
and an unresolved blend all survive.

The candidate OPH carrier specified here is a bounded, self-reading software
classifier over conventional accretion forward models. Its ports ingest
calibrated spectra, event times, distance and counterpart information, orbital
constraints, nebular calorimetry, and survey selection. It preserves raw records, reports
cross-model residuals and ambiguity, and accepts only predeclared reconciliation
moves. Coherent pulsations can identify a neutron star; a donor mass function
can constrain mass; fast winds favor supercritical-flow models; nebulae
constrain angle-integrated power. No one of these measurements fills every
coordinate.
The conditional audit states which coordinate each observation reaches and
refuses a unique label when several frozen models remain viable.

## 1. The classification question

A ULX is conventionally an off-nuclear, point-like extragalactic source with
an apparent X-ray luminosity above about
$10^{39}\ \mathrm{erg\,s^{-1}}$. The exact band varies among catalogs. The
threshold is a catalog convention rather than a phase boundary. “Extreme ULX”
and “hyperluminous X-ray source” introduce higher conventional thresholds.
Sutton and collaborators selected extreme ULXs above $5\times10^{40}$ and
used the HLX label above $10^{41}\ \mathrm{erg\,s^{-1}}$
[Sutton et al. (2012)](https://doi.org/10.1111/j.1365-2966.2012.20944.x).

The public claim boundary is:

```text
apparent luminosity is not compact-object mass
an Eddington ratio cannot be inferred from apparent luminosity without mass, effective opacity, and angular transfer
a curved spectrum is not an accretor identity
a missing pulse is not a black-hole detection
a QPO frequency is not a mass without a mode and radius
a young environment is not an individual engine measurement
a normal-form label is not an OPH-derived accretion law
```

Published observations support a population-level conclusion while leaving many
individual sources conditional. Multiple ULXs contain neutron stars. Spectra,
winds, bubbles, and nebular ionization support supercritical accretion across
much of the ordinary population. M101 ULX-1 has reported dynamics favoring a
stellar-mass black hole, conditional on the Wolf–Rayet line tracer, donor mass,
and inclination. ESO 243-49 HLX-1 and M82 X-1 retain serious IMBH evidence,
although neither has a dynamical IMBH mass. Supersoft sources, transients,
blends, supernova remnants,
foreground objects, and background active galactic nuclei require explicit
alternatives. A single branch for all ULXs is ruled out by observation.

## 2. Source OPH branch and bounded carrier

The source branch is the finite observer-patch semantics in
[Screen Microphysics and Observer Synchronization](../paper/screen_microphysics_and_observer_synchronization.tex)
and the protected-record, accepted-repair, and normal-form conditions in
[Reality as Consensus Protocol](../paper/reality_as_consensus_protocol.tex).
The compact-object forward-model discipline follows the local
[JWST compact-object audit](jwst_compact_object_source_release.md) and the
exposure and censoring rules in
[Compact Record Transients](compact_record_transients.md).

One classifier patch is declared for each triple

$$
(\text{source identity},\ \text{frozen data generation},
\text{frozen model generation}).
\tag{1}
$$

Its observer-like structure is operational and testable:

| OPH component | ULX classifier realization |
|---|---|
| local state | finite candidate families, nuisance parameters, response and calibration identifiers, uncertainty law, residuals, ambiguity flags, and claim tier |
| ports and boundaries | X-ray counts and backgrounds; event times, gaps, exposure, and dead time; distance and flux; counterparts and environment; orbital, dynamical, magnetic, wind, and nebular constraints; survey selection and non-detections |
| readback | per-family likelihood and port residuals, posterior-predictive checks, viable-set size, missing-record mask, and classification status |
| durable records | immutable observation identifiers and hashes, event selections, extraction regions, response matrices, background products, timing windows, distance posterior, model equations, priors, thresholds, and every candidate score |
| repair moves | canonicalization, nuisance fitting inside frozen families, ingestion of a newly versioned port, and reconciliation of cross-epoch or cross-instrument mismatch |
| checkpoints | preregistration; raw-data and calibration freeze; extraction freeze; fit and residual freeze; held-out epoch or instrument check; selection/censoring audit; public manifest |
| evidence bundle | data locators, raw and reduced hashes, calibrated responses, code and dependency hashes, exact likelihood, all alternatives and upper limits, controls, held-out results, and a machine-readable source table |

Repair applies to the classifier state. It may never rewrite an event list,
change a counterpart, move a threshold, or replace a prior after seeing which
branch wins. A changed family or threshold begins a new model generation. The
accretion flow itself is not called an OPH repair process here; such language
would require a separate physical self-reading and feedback theorem for the
source.

The classifier is therefore a specified candidate self-reading patch. It reads its
own residual and missing-record state, uses those records to determine the next
allowed fit or data-ingestion move, and exposes a reproducible checkpoint. An
astronomical source may also exhibit feedback through torque, radiation
pressure, advection, and winds, but conventional accretion physics supplies
those mechanisms in this article.

## 3. Replace the four labels with physical coordinates

For a **single compact-accretor hypothesis**, let the latent source description
be

$$
z=(a,M,\kappa_{\rm eff},\lambda,b,i,s),
\tag{2}
$$

where

- $a$ is accretor identity: neutron star, black hole, white dwarf, another
  declared object, or unknown;
- $M$ is mass, with “stellar” and “intermediate” treated as declared ranges
  rather than laws of nature;
- $\kappa_{\rm eff}$ is the composition-, field-, and geometry-dependent
  opacity entering the declared Eddington reference;
- $\lambda=L_{\rm true}/L_{\rm Edd}$ is the intrinsic accretion-luminosity
  regime;
- $b$ and inclination $i$ describe angular transfer and viewing geometry;
- $s$ distinguishes persistent binary, transient binary, and other declared
  single-engine modes.

An unresolved blend instead has a finite collection or point process
$\{z_j\}_{j=1}^N$ plus the instrument confusion model. A foreground star,
background AGN, supernova, or remnant is a separate contaminant hypothesis
whose accretor coordinates $(a,M,\kappa_{\rm eff},\lambda,b,i)$ may be
undefined. The classifier must not force those alternatives into the scalar
tuple (2).

The legacy single-engine labels are projections of (2):

$$
\begin{aligned}
\mathcal B_{\rm NS,beam}
  &=\{z:a={\rm NS},\ b<b_{\rm cut}\},\\
\mathcal B_{\rm sBH,super}
  &=\{z:a={\rm BH},\ M\in\mathcal M_{\rm stellar},\ \lambda>1\},\\
\mathcal B_{\rm IMBH}
  &=\{z:a={\rm BH},\ M\in\mathcal M_{\rm IMBH}\}.
\end{aligned}
\tag{3}
$$

Equation (3) exposes the defect in the original ansatz. The first family fixes
identity and geometry, the second fixes identity, mass scale, and accretion
regime, while the third fixes identity and mass scale alone. A neutron-star ULX
can radiate supercritically with weak beaming. An IMBH can also accrete above
Eddington or emit anisotropically. A sub-Eddington stellar black hole can cross
the catalog threshold near the high end of the stellar mass distribution.

The output must retain six statuses that a forced four-way label would erase:

| Status | Meaning |
|---|---|
| identified coordinate | a hard receipt fixes accretor, mass range, regime, or geometry |
| ambiguous | two or more frozen families fit every mandatory port |
| insufficient record | a required port or sensitivity bound is absent |
| out of distribution | the record lies outside the model, calibration, or survey-support domain, so adequacy has not been tested there |
| model-set failure | the record is complete enough to reject every declared family |
| other physical branch | a separately declared alternative forward model passes |

A background quasar is not an exotic ULX engine. It fails host association and
leaves the physical ULX catalog. Ambiguity is likewise not an “other” object.

## 4. The public record and detector forward model

A minimum record packet is

$$
\begin{aligned}
R_{\rm ULX}=(&D,F_X,N_H,C_{\rm bol},S(E,t),P,\dot P,
\operatorname{PSD},P_{\rm orb},K_2,e_{\rm orb},\\
&\mathcal N_{\rm neb},\mathcal E_{\rm host},
\mathcal W_{\rm obs},\mathcal S_{\rm survey}).
\end{aligned}
\tag{4}
$$

Here $D$ is a distance posterior rather than a point copied from a catalog,
$C_{\rm bol}$ records the band-to-bolometric transformation, and
$\mathcal W_{\rm obs}$ contains good-time intervals, cadence, gaps, energy
selection, pile-up treatment, and dead time. The survey term records sky
coverage, sensitivity, source confusion, and non-detections.

For nearby sources,

$$
L_{\rm iso,bol}=4\pi D^2F_{X,{\rm unabs}}C_{\rm bol}
\tag{5}
$$

or numerically

$$
L_{\rm iso,bol}=1.1965\times10^{40}
\left(\frac{D}{10\ {\rm Mpc}}\right)^2
\left(\frac{F_{\rm bol}}{10^{-12}\,
{\rm erg\,cm^{-2}\,s^{-1}}}\right)
\ {\rm erg\,s^{-1}}.
\tag{6}
$$

Ignoring covariance, the leading propagated uncertainty is

$$
\sigma_{\ln L}^2\simeq
4\left(\frac{\sigma_D}{D}\right)^2
+\left(\frac{\sigma_{F_{X,{\rm unabs}}}}{F_{X,{\rm unabs}}}\right)^2
+\sigma_{\ln C_{\rm bol}}^2.
\tag{7}
$$

A spectral model is compared to detector counts rather than to an unfolded
plot. For epoch $e$ and detector channel $k$, a conventional forward model
has expected counts

$$
\mu_{ek}(h,\theta_h)=t_e\int
R_{ek}(E)A_e(E)S_h(E;\theta_h)\,dE+B_{ek},
\tag{8}
$$

where $R$ is redistribution, $A$ effective area, $B$ background, and
$h$ a candidate family. Timing models are multiplied by the actual window
and dead-time operator before comparison with event arrivals. This is why
“luminosity plus timing plus environment” is too small a public record:
spectrum, instrument response, source association, and selection are
load-bearing.

## 5. Luminosity identifiability theorem

For fully ionized material with hydrogen mass fraction $X$, the
electron-scattering opacity and classical Eddington luminosity are

$$
\kappa_{\rm es}\simeq0.20(1+X)\ \mathrm{cm^2\,g^{-1}},
\qquad
L_{\rm Edd}=\frac{4\pi GMc}{\kappa_{\rm es}}
\simeq\frac{2.50\times10^{38}}{1+X}
\frac{M}{M_\odot}\ \mathrm{erg\,s^{-1}}.
\tag{9}
$$

For pure hydrogen, retaining the unrounded
$\kappa_{\rm es}=\sigma_T/m_p=0.39773\ \mathrm{cm^2\,g^{-1}}$, this is
$1.2571\times10^{38}(M/M_\odot)\ \mathrm{erg\,s^{-1}}$; the rounded
$0.40\ \mathrm{cm^2\,g^{-1}}$ opacity gives $1.2500\times10^{38}$.
Composition changes the coefficient. Strong magnetic fields change the
scattering cross-section and column geometry, so (9) is a reference scale
rather than a universal neutron-star ceiling.

In equations (10) through (17), $L_{\rm iso}$ and $L_{\rm true}$ are
bolometric luminosities. A band-limited measurement enters only after the
declared absorption and bolometric transformation in equation (5). For
emission collimated toward the observer, define

$$
b=\frac{L_{\rm true}}{L_{\rm iso}},\quad 0<b\le1,
\qquad
\lambda=\frac{L_{\rm true}}{L_{\rm Edd}}.
\tag{10}
$$

Some papers use the reciprocal beaming convention. A general angular transfer
function can also describe off-axis views; equation (10) is the convention for
the toward-us beaming comparator.

The luminosity factor $b=L_{\rm true}/L_{\rm iso}$ is not, in general, the
fraction of randomly oriented systems visible to a survey. Define that separate
orientation-selection probability as $f_\Omega$, including the full angular
emission pattern, flux threshold, and obscuration. Only for an ideal top-hat
beam with uniform intensity inside a total solid angle $4\pi f_\Omega$, no
off-axis detection, and the same angle-integrated luminosity does
$b=f_\Omega$. Population corrections below use $f_\Omega$, not $b$.

Let $L_{\rm Edd,\odot}$ be the one-solar-mass coefficient for the declared
composition. Then

$$
q:=\frac{L_{\rm iso}}{L_{\rm Edd,\odot}}
=\frac{\lambda}{b}\frac{M}{M_\odot}.
\tag{11}
$$

### Identifiability result 1: luminosity cannot classify the accretor

For any observed $q>0$ and any proposed $M>0$, positive values of
$\lambda$ and $b$ can satisfy (11). In logarithmic coordinates,

$$
\ln q=\ln(M/M_\odot)+\ln\lambda-\ln b
\tag{12}
$$

has a one-row Jacobian with respect to three unknowns. Perturbations satisfying

$$
\delta\ln M+\delta\ln\lambda-\delta\ln b=0
\tag{13}
$$

are invisible to luminosity. Two parameter directions remain unidentified.
This proves that luminosity alone cannot distinguish a neutron star, a stellar
black hole, or an IMBH.

For a comparator family with

$$
M\in[M_-,M_+],\quad
\lambda\in[\lambda_-,\lambda_+],\quad
b\in[b_-,b_+],
\tag{14}
$$

the luminosity-feasible interval is

$$
q\in\left[
\frac{\lambda_-M_-}{b_+M_\odot},
\frac{\lambda_+M_+}{b_-M_\odot}
\right].
\tag{15}
$$

A branch survives whenever its interval overlaps the uncertainty interval for
$q$. If external evidence establishes $b\ge b_{\min}$ and
$\lambda\le\lambda_{\max}$, luminosity gives the conditional bound

$$
\frac{M}{M_\odot}\ge
\frac{L_{\rm iso,-}}{L_{\rm Edd,\odot}}
\frac{b_{\min}}{\lambda_{\max}}.
\tag{16}
$$

With no positive lower bound on $b$, or no finite upper bound on
$\lambda$, equation (16) supplies no positive mass bound.

### A reproducible $10^{40}\ \mathrm{erg\,s^{-1}}$ example

Take $D=10$ Mpc and
$F_{\rm bol}=8.358\times10^{-13}\,\mathrm{erg\,cm^{-2}\,s^{-1}}$.
Equation (5) gives exactly
$L_{\rm iso}=10^{40}\ \mathrm{erg\,s^{-1}}$ to the shown precision. Using
the pure-hydrogen reference in (9), $q=79.55$.

| Candidate | Assumed mass | Required $\lambda/b$ | One feasible realization |
|---|---:|---:|---|
| neutron star | $1.4M_\odot$ | 56.82 | $\lambda=4.42,\ b=0.0778$ |
| stellar black hole | $10M_\odot$ | 7.955 | $\lambda=3.55,\ b=0.446$ |
| heavy stellar black hole | $30M_\odot$ | 2.652 | $\lambda=2.65,\ b=1$ |
| IMBH | $1000M_\odot$ | 0.0795 | $\lambda=0.0795,\ b=1$ |
| unresolved blend | ten $10M_\odot$ binaries | 0.795 per source | ten nearly Eddington components |

The first two realizations come from the illustrative supercritical
prescription. Here King’s convention is
$\dot m:=0.1\dot M c^2/L_{\rm Edd}$, equivalently
$\dot M/[L_{\rm Edd}/(0.1c^2)]$; changing the reference efficiency changes
the numerical $\dot m$. With that convention,

$$
\lambda\simeq1+\ln\dot m,
\qquad
b\simeq\min\left(1,\frac{73}{\dot m^2}\right),
\tag{17}
$$

which gives $\dot m=30.63$ for the neutron star and 12.79 for the
$10M_\odot$ black hole. King proposed the beaming law for
$\dot m\gtrsim8.5$
[King (2009)](https://doi.org/10.1111/j.1745-3933.2008.00594.x).
Poutanen and collaborators derived a related supercritical model with
$L/L_{\rm Edd}\simeq1+0.6\ln\dot m$ and face-on amplification of roughly
2 to 7
[Poutanen et al. (2007)](https://doi.org/10.1111/j.1365-2966.2007.11668.x).
These prescriptions demonstrate feasibility. Neither is a universal inversion
law for an individual source.

## 6. What each additional record can establish

### Pulsations and neutron stars

Coherent accretion-powered pulsations with orbital Doppler motion identify a
rotating magnetic accretor. They identify a neutron star when the period and
other records exclude the declared magnetic-white-dwarf family. For the cited
PULXs that conclusion uses the full fast-period, hard-X-ray luminosity,
spin-up or torque, and orbital record; the period alone is not a universal
white-dwarf exclusion. They do not determine $b$. Pulse absence only
supplies an upper limit over the searched period, acceleration, energy band,
pulsed fraction, and window; scattering, dilution, geometry, and intermittency
can hide a neutron star.

Spin-up gives an accretion-torque consistency check:

$$
N_{\rm obs}=2\pi I\dot\nu
=-2\pi I\frac{\dot P}{P^2},
\tag{18}
$$

$$
R_{\rm co}=\left(\frac{GMP^2}{4\pi^2}\right)^{1/3},
\qquad
R_m\simeq\xi\left(\frac{\mu^4}{2GM\dot M^2}\right)^{1/7}.
\tag{19}
$$

Accretion requires approximately $R_m\lesssim R_{\rm co}$. Torque, a
propeller transition, and luminosity can jointly constrain magnetic moment,
accretion rate, and beaming, although the result inherits the torque model.
A secure cyclotron feature can constrain the local field. Particle
identification, harmonics, phase dependence, and line significance have to
survive the audit.

### Dynamics and black-hole mass

For donor radial-velocity amplitude $K_2$, orbital period $P_{\rm orb}$,
and eccentricity $e$,

$$
f(M_X)=\frac{P_{\rm orb}K_2^3}{2\pi G}(1-e^2)^{3/2}
=\frac{M_X^3\sin^3 i}{(M_X+M_2)^2}\le M_X.
\tag{20}
$$

In convenient units,

$$
\frac{f(M_X)}{M_\odot}=1.036\times10^{-7}
\left(\frac{P_{\rm orb}}{\rm day}\right)
\left(\frac{K_2}{\rm km\,s^{-1}}\right)^3
(1-e^2)^{3/2}.
\tag{21}
$$

A verified donor line and mass function provide the most direct route from the
record to mass. Irradiated winds and disk emission can move the line away from
the donor’s center of mass, so tracer fidelity belongs in the receipt. A mass
lower bound above the declared maximum neutron-star mass selects a black hole
within the declared neutron-star/black-hole families only conditional on the
orbital solution and tracer fidelity. An IMBH branch is forced within the
declared black-hole families only when that lower bound exceeds the
preregistered maximum stellar-remnant mass. That inequality is sufficient to
exclude stellar masses within those declared families, subject to the orbit
and tracer receipt. Donor-mass and inclination constraints are additionally
needed to bound the actual mass from above or confine it to a narrower
interval. None of these results fixes $\lambda$ or $b$.

### QPOs and variability

The gravitational clock and a Keplerian frequency are

$$
t_g=\frac{GM}{c^3}=4.9256\ \mu{\rm s}\frac{M}{M_\odot},
\tag{22}
$$

$$
f_K(R)=3.231\times10^4
\left(\frac{M}{M_\odot}\right)^{-1}R^{-3/2}\ {\rm Hz},
\qquad R=r/r_g.
\tag{23}
$$

A 0.1 Hz feature can correspond to $R\simeq1014r_g$ for a
$10M_\odot$ object or $R\simeq47r_g$ for a $1000M_\odot$ object.
Frequency identifies $M^{-1}f(\text{mode},R,a_*)$, not mass alone. The
3:2 peaks at $3.32\pm0.06$ and $5.07\pm0.06$ Hz in M82 X-1 gave
$428\pm105M_\odot$ under stellar-black-hole inverse-frequency scaling
[Pasham, Strohmayer, and Mushotzky (2014)](https://doi.org/10.1038/nature13710).
That is conditional IMBH evidence rather than a dynamical mass.

### Spectra and winds

High-count ULX spectra often show a soft excess and curvature above about 3
keV. Gladstone, Roberts, and Done fitted cool, optically thick Comptonizing
components with optical depths around 5 to 30 and argued for an ultraluminous
accretion state
[Gladstone et al. (2009)](https://doi.org/10.1111/j.1365-2966.2009.15123.x).
The phenotype supports thick or supercritical flow models. It cannot by itself
separate a neutron star from a black hole, and a fitted cool component is not
an IMBH thermometer when a wind photosphere obscures the inner disk.

Resolved lines in NGC 1313 X-1 and NGC 5408 X-1 indicate winds near $0.2c$
[Pinto, Middleton, and Fabian (2016)](https://doi.org/10.1038/nature17417).
This is direct evidence for fast outflowing gas and strong evidence favoring a
supercritical-flow model. Wind power depends on ionization, density, filling
factor, solid angle, and launch radius. The wind does not announce the compact
object’s identity.

### Nebular calorimetry and environment

A surrounding nebula reads the source in directions other than our line of
sight and over a longer clock. Under case-B recombination, an ionization-bounded
He III region gives approximately

$$
Q({\rm He^+})=
\frac{\alpha_B({\rm He^{++}})}{\alpha^{\rm eff}_{4686}}
\frac{L_{4686}}{h\nu_{4686}}.
\tag{24}
$$

The covering fraction, density, abundance, unobserved extreme-UV spectrum,
and recombination history enter the conversion. For Holmberg II X-1,
$L_{4686}=2.7\times10^{36}\ \mathrm{erg\,s^{-1}}$ and photoionization
models inferred a lower bound of about $4$ to
$6\times10^{39}\ \mathrm{erg\,s^{-1}}$ on the X-ray luminosity, conditional
on extrapolating the fitted spectrum from 300 eV down to the 54 eV ionization
edge. The morphology was inconsistent with narrow beaming
[Kaaret, Ward, and Zezas (2004)](https://doi.org/10.1111/j.1365-2966.2004.08020.x).
That receipt constrains the angular emission pattern and, conditional on the
extreme-UV and nebular model, $b$ for this source. It does not establish
exact isotropy at every epoch.

Young, metal-poor star-forming regions raise the prior probability for
high-mass X-ray binaries and massive stellar remnants. Old clusters change the
prior. Neither environment fixes the individual accretor. In a 40 Mpc census,
Kovlakas and collaborators reported 629 candidates in 309 galaxies and an
estimated foreground/background contamination near 20 percent
[Kovlakas et al. (2020)](https://doi.org/10.1093/mnras/staa2481). Population
classification must carry exposure, sensitivity, source confusion, and the
non-detections. For an ideal two-sided top-hat beam covering total solid-angle
fraction $f_\Omega$ under random orientation, the parent population scales as

$$
N_{\rm parent}\sim\frac{N_{\rm detected}}{f_\Omega}
$$

before other selection corrections. For a structured beam or a flux-limited
survey, $f_\Omega$ must instead be computed by integrating the detection
probability over orientation, distance, luminosity, absorption, and exposure.
It must not be replaced by the line-of-sight luminosity factor $b$ without the
top-hat assumptions above. Extreme collimation can therefore fail a
population-count check even when it fits one source’s luminosity.

## 7. Frozen classifier and identifiability rule

For hypothesis family $h$, parameters $\theta_h$, port data $D_e$, and
forward prediction $F_{h,e}$, define a predeclared mismatch

$$
\Phi_h(\theta_h;R)=
\sum_e w_e\,d_e\!\left(D_e,F_{h,e}(\theta_h)\right)
+\Phi_{\rm audit}.
\tag{25}
$$

The distances may be Poisson likelihood terms for counts, event likelihoods
for timing, astrometric association probabilities, or calibrated residuals for
nebular and dynamical ports. Audit penalties cover a missing mandatory port,
failed convergence, unmodeled pile-up, posterior-predictive failure, or a
held-out threshold violation. Priors derived from environment enter explicitly
rather than masquerading as measurements of identity.

The viable fiber is

$$
\mathcal C_R=
\left\{(h,\theta_h):
\text{every mandatory adequacy, calibration, and held-out gate passes}
\right\}.
\tag{26}
$$

Project this parameter fiber to the declared factorized normal-form cells:

$$
\mathcal V_R=\pi_{\rm NF}(\mathcal C_R).
\tag{27}
$$

The projection removes silent implementation labels and groups parameter
points only when they occupy the same published accretor, mass-range, regime,
geometry, and source-mode cell. It retains the parameter interval inside that
cell and never identifies physically different accretors or declared mass
ranges.

The specified canonical report is

$$
\operatorname{NF}_{\rm spec}(R):=
\operatorname{Canon}(\mathcal V_R,\ \text{missing mask},
\text{domain status},\ \text{claim tier}),
$$

with a frozen sort order and serialization. This is a normal-form contract,
not an earned confluence theorem. An implementation may claim the OPH normal
form only after it publishes protected-boundary preservation, a strict descent
measure for every accepted reconciliation, termination, schedule independence,
and a reproducible evidence bundle. This article supplies the contract rather
than those implementation receipts.

### Conditional classifier rule 2: retain every viable cell

The output has three mathematical cases:

$$
\begin{array}{c|l}
\mathcal V_R=\varnothing & \text{no declared normal-form cell closes}\\
|\mathcal V_R|=1 & \text{one identifiable cell within the frozen model set}\\
|\mathcal V_R|>1 & \text{ambiguity among the surviving physical cells.}
\end{array}
\tag{28}
$$

A highest posterior probability does not change the cardinality in (28).
Within the frozen family list, a conditionally unique output requires
preregistered adequacy, calibration, and separation, not merely a winner. Empty
fibers are split by the missing-record mask: absent mandatory data yield
`INSUFFICIENT_RECORD`; a record outside the declared validation or
survey-support domain yields `OUT_OF_DISTRIBUTION`; complete in-domain data
rejecting every family yield `MODEL_SET_FAILURE`. A declared alternative that
passes becomes `OTHER_PHYSICAL_BRANCH`.

The decision order is:

1. establish host association, off-nuclear position, distance, point-source
   identity, and deblending;
2. apply hard identity receipts such as coherent pulsations or a reliable mass
   function;
3. test the full $(M,\lambda,b,X,C_{\rm bol})$ feasibility intervals;
4. evaluate spectra, timing, winds, nebulae, orbital data, and environment by
   their own forward models;
5. run held-out epochs or instruments, alternatives, and selection/censoring
   checks;
6. emit every surviving coordinate interval and the status in (28).

## 8. Observational branch audit

| Source or evidence | Receipt | Valid normal-form statement | Unresolved coordinate |
|---|---|---|---|
| M82 X-2 | 1.37 s pulsations, 2.5-day orbital modulation, pulsed $L_{3-30}=4.9\times10^{39}\ \mathrm{erg\,s^{-1}}$ | neutron-star accretor identified | true luminosity, beaming, magnetic geometry |
| NGC 5907 ULX-1 | spin changed 1.43 to 1.13 s; apparent peak near $10^{41}$ to $2\times10^{41}\ \mathrm{erg\,s^{-1}}$ | neutron star identified; hyperluminous appearance does not require an IMBH | split among beaming, magnetic-opacity effects, and intrinsic supercritical power |
| NGC 7793 P13 | 0.42 s pulsations and persistent spin-up | neutron star identified | torque-derived field and beaming remain model-dependent |
| NGC 1313 X-1 and NGC 5408 X-1 | resolved $\sim0.2c$ winds | fast outflow identified; supercritical-flow model strongly favored | accretor identity and intrinsic Eddington ratio |
| Holmberg II X-1 | He II nebular photon counter and non-narrow morphology | substantial angle-integrated power; narrow beam disfavored | mass and accretor identity |
| M101 ULX-1 | reported 8.2-day orbit, mass function $0.18\pm0.03M_\odot$, likely $20$ to $30M_\odot$ solution | stellar-black-hole candidate; supersoft temperature alone did not require an IMBH | line origin, donor mass, inclination, and exact regime |
| ESO 243-49 HLX-1 | host-consistent Hα redshift, peak near $1.1$ to $1.3\times10^{42}\ \mathrm{erg\,s^{-1}}$, state changes and transient radio emission | strong IMBH candidate under modest anisotropy and sub- or near-Eddington scaling | no dynamical mass; beamed supercritical alternatives are not mathematically excluded |
| M82 X-1 | stable 3:2 QPO pair and conditional $428\pm105M_\odot$ scaling | model-dependent IMBH candidate | mode identity, scaling transfer, crowding, dynamical mass |
| IC 4320 projected HLX candidate | optical spectrum found a background $z\simeq2.84$ quasar | host association fails; remove from physical ULX class | none within the ULX engine classifier |

The M101 result is from
[Liu et al. (2013)](https://doi.org/10.1038/nature12762); the HLX-1 discovery
and peak luminosity are from
[Farrell et al. (2009)](https://doi.org/10.1038/nature08083); the Hα redshift
and radio/state records are from
[Wiersema et al. (2010)](https://doi.org/10.1088/2041-8205/721/2/L102) and
[Webb et al. (2012)](https://doi.org/10.1126/science.1222779). The background
quasar classification is documented by
[Sutton et al. (2015)](https://doi.org/10.1093/mnras/stv505).

Ultraluminous supersoft sources deserve a phenotype coordinate rather than an
automatic IMBH label. Temperatures near 0.1 keV and large fitted photospheres
can arise from an optically thick supercritical wind, especially away from the
funnel; nuclear-burning white dwarfs and thin IMBH disks occupy parts of the
same spectral space. Absorption and bolometric correction are severe. The
classifier keeps these families alive until timing, radius evolution,
dynamics, or another port separates them
[Urquhart and Soria (2016)](https://doi.org/10.1093/mnras/stv2293).

## 9. Claim ladder and public evidence

| Tier | Claim earned |
|---|---|
| `ULX0_CATALOG_RECORD` | off-nuclear candidate with stated band, flux, distance, astrometric uncertainty, and threshold |
| `ULX1_RESPONSE_CORRECTED_RECORD` | deblended multi-epoch counts, responses, backgrounds, timing windows, absorption, bolometric convention, counterpart probabilities, and upper limits |
| `ULX2_FROZEN_FAMILY_FIT` | all declared accretion and contaminant families fitted under frozen priors and thresholds |
| `ULX3_DEGENERACY_AUDITED_CANDIDATE` | viable fiber, missing-record mask, hard identity receipts, and rejected alternatives published |
| `ULX4_HELDOUT_STABLE_CANDIDATE` | branch coordinates survive an independent epoch, instrument, or preregistered predictive test |
| `ULX5_SELECTION_CORRECTED_POPULATION` | exposure, sensitivity, background, censoring, orientation, and non-detections enter the population likelihood |
| `ULX6_SOURCE_ONLY_OPH_PREDICTION` | OPH-derived accretion source law, radiation and timing transfer, population law, no-target-leak ledger, and frozen likelihood |

`ULX0` through `ULX5` are observational and analysis-maturity tiers for a
classifier built from conventional astrophysics. They are not steps toward
`ULX6`. This note specifies the classifier and receipt contract and summarizes
published evidence, but a source earns a tier only through its own public
bundle; the summary table is not a substitute for those raw and reduced
receipts. `ULX6` is an orthogonal source-origin claim requiring an OPH-derived
accretion and population law. The local OPH corpus contains no such law, so no
ULX abundance, mass function, luminosity function, or branch frequency is
promoted as an OPH prediction.

A public `ULX3` or higher bundle must contain raw observation identifiers,
event and extraction hashes, response and background products, exact timing
windows, distance and association posteriors, model equations and priors, all
candidate scores, residuals and posterior-predictive checks, pulse-search trial
accounting, upper limits, alternative counterparts, nebular assumptions,
selection functions, and the deterministic code needed to regenerate the
classification. Publishing only the winning fit fails the record test.

## 10. Falsifiers and conditional audit outcome

The classifier fails on a source if repeated execution from the same frozen
record changes the viable fiber, if a candidate survives only after moving a
threshold or prior, if held-out data reject its forward prediction, or if an
omitted contaminant closes the record better. A population claim fails when
the inferred parent count, including the $1/f_\Omega$ orientation correction
under the declared angular selection model,
conflicts with the available binary population or when the result disappears
after exposure and background correction.

Each physical coordinate has its own falsifier. A neutron-star claim falls if
the supposed pulse loses coherence, fails source localization, or is explained
by the observing window. An IMBH claim falls if reliable dynamics give a
stellar mass. A narrow-beam claim falls when surrounding gas requires large
angle-integrated ionizing power. A supercritical-flow claim falls if a
sub-Eddington model with independently measured mass and geometry predicts the
full spectrum, timing, and outflow record more accurately.

The conditional classifier and audit specification supports the following
conclusions; it is not a complete physical account of ULXs:

- ULXs are a luminosity-selected mixture rather than one engine class.
- A single-compact-accretor hypothesis factors into accretor, mass, opacity,
  intrinsic regime, geometry, viewing angle, and source mode; blends and
  contaminants use separate composite or alternative state spaces.
- Luminosity measures $M\lambda/b$; equations (11) to (16) prove the
  two-dimensional degeneracy and state the assumptions needed for a mass bound.
- Neutron-star-specific pulsation records identify neutron stars, dynamics
  constrain mass, winds provide strong evidence favoring supercritical-flow
  models, nebulae constrain the angular emission pattern, and environment
  supplies a prior. These receipts are not interchangeable.
- A frozen, set-valued classifier returns an identified candidate, ambiguity,
  insufficient record, out-of-distribution status, model-set failure, or a
  declared other branch.
- The source-only OPH accretion and population certificate is absent and no
  empirical fit is mislabeled as one.

The remaining empirical work is object-specific and population-specific. Within
the frozen candidate set, a source reaches a conditionally unique classifier
output when its multi-port record makes the projected viable set in (27) a
singleton and that cell survives a held-out checkpoint. This does not prove
that the candidate set is exhaustive or that the surviving engine is an
OPH-derived source. Sources without those receipts remain ambiguous, and
reporting that ambiguity is a valid audit outcome rather than a completed
physical classification.

## References

- P. Kaaret, H. Feng, and T. P. Roberts, “Ultraluminous X-Ray Sources,”
  *Annual Review of Astronomy and Astrophysics* **55**, 303–341 (2017),
  [doi:10.1146/annurev-astro-091916-055259](https://doi.org/10.1146/annurev-astro-091916-055259).
- A. King, J.-P. Lasota, and M. Middleton, “Ultraluminous X-ray sources,”
  *New Astronomy Reviews* **96**, 101672 (2023),
  [doi:10.1016/j.newar.2022.101672](https://doi.org/10.1016/j.newar.2022.101672).
- M. Bachetti et al., “An ultraluminous X-ray source powered by an accreting
  neutron star,” *Nature* **514**, 202–204 (2014),
  [doi:10.1038/nature13791](https://doi.org/10.1038/nature13791).
- G. L. Israel et al., “An accreting pulsar with extreme properties drives an
  ultraluminous x-ray source in NGC 5907,” *Science* **355**, 817–819 (2017),
  [doi:10.1126/science.aai8635](https://doi.org/10.1126/science.aai8635).
- F. Fürst et al., “Discovery of coherent pulsations from the ultraluminous
  X-ray source NGC 7793 P13,” *Astrophysical Journal Letters* **831**, L14
  (2016),
  [doi:10.3847/2041-8205/831/2/L14](https://doi.org/10.3847/2041-8205/831/2/L14).
- J. Poutanen, G. Lipunova, S. Fabrika, A. G. Butkevich, and P. Abolmasov,
  “Supercritically accreting stellar mass black holes as ultraluminous X-ray
  sources,” *MNRAS* **377**, 1187–1194 (2007),
  [doi:10.1111/j.1365-2966.2007.11668.x](https://doi.org/10.1111/j.1365-2966.2007.11668.x).
- A. R. King, “Masses, beaming and Eddington ratios in ultraluminous X-ray
  sources,” *MNRAS Letters* **393**, L41–L44 (2009),
  [doi:10.1111/j.1745-3933.2008.00594.x](https://doi.org/10.1111/j.1745-3933.2008.00594.x).
- J. C. Gladstone, T. P. Roberts, and C. Done, “The ultraluminous state,”
  *MNRAS* **397**, 1836–1851 (2009),
  [doi:10.1111/j.1365-2966.2009.15123.x](https://doi.org/10.1111/j.1365-2966.2009.15123.x).
- C. Pinto, M. J. Middleton, and A. C. Fabian, “Resolved atomic lines reveal
  outflows in two ultraluminous X-ray sources,” *Nature* **533**, 64–67
  (2016), [doi:10.1038/nature17417](https://doi.org/10.1038/nature17417).
- P. Kaaret, M. J. Ward, and A. Zezas, “High-resolution imaging of the He II
  4686 emission-line nebula associated with the ultraluminous X-ray source in
  Holmberg II,” *MNRAS* **351**, L83–L88 (2004),
  [doi:10.1111/j.1365-2966.2004.08020.x](https://doi.org/10.1111/j.1365-2966.2004.08020.x).
- S. A. Farrell et al., “An intermediate-mass black hole of over 500 solar
  masses in the galaxy ESO 243-49,” *Nature* **460**, 73–75 (2009),
  [doi:10.1038/nature08083](https://doi.org/10.1038/nature08083).
- K. Wiersema et al., “The nature of the brightest ultraluminous X-ray source:
  host galaxy and optical counterpart of ESO 243-49 HLX-1,” *Astrophysical
  Journal Letters* **721**, L102–L106 (2010),
  [doi:10.1088/2041-8205/721/2/L102](https://doi.org/10.1088/2041-8205/721/2/L102).
- N. A. Webb et al., “Radio detections during two state transitions of the
  intermediate-mass black hole HLX-1,” *Science* **337**, 554–556 (2012),
  [doi:10.1126/science.1222779](https://doi.org/10.1126/science.1222779).
- D. R. Pasham, T. E. Strohmayer, and R. F. Mushotzky, “A 400-solar-mass black
  hole in the galaxy M82,” *Nature* **513**, 74–76 (2014),
  [doi:10.1038/nature13710](https://doi.org/10.1038/nature13710).
- J.-F. Liu, J. N. Bregman, Y. Bai, S. Justham, and P. Crowther, “Puzzling
  accretion onto a black hole in the ultraluminous X-ray source M101 ULX-1,”
  *Nature* **503**, 500–503 (2013),
  [doi:10.1038/nature12762](https://doi.org/10.1038/nature12762).
- K. Kovlakas et al.,
  “A census of ultraluminous X-ray sources in the local Universe,” *MNRAS*
  **498**, 4790–4810 (2020),
  [doi:10.1093/mnras/staa2481](https://doi.org/10.1093/mnras/staa2481).
- A. D. Sutton, T. P. Roberts, D. J. Walton, J. C. Gladstone, and A. E. Scott,
  “The most extreme ultraluminous X-ray sources: evidence for
  intermediate-mass black holes?” *MNRAS* **423**, 1154–1177 (2012),
  [doi:10.1111/j.1365-2966.2012.20944.x](https://doi.org/10.1111/j.1365-2966.2012.20944.x).
- A. D. Sutton, T. P. Roberts, J. C. Gladstone, and D. J. Walton, “The
  hyperluminous X-ray source candidate in IC 4320: another HLX bites the dust,”
  *MNRAS* **450**, 787–793 (2015),
  [doi:10.1093/mnras/stv505](https://doi.org/10.1093/mnras/stv505).
- R. Urquhart and R. Soria, “Optically thick outflows in ultraluminous
  supersoft sources,” *MNRAS* **456**, 1859–1880 (2016),
  [doi:10.1093/mnras/stv2293](https://doi.org/10.1093/mnras/stv2293).
