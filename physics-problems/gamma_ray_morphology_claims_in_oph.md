# Gamma-Ray Morphology Claims In OPH

## Motivating Result

This note entered the queue after Fermi-LAT analyses reported a diffuse
gamma-ray dipole-like feature larger than expected, with an axis different from
the CMB dipole
([NASA Goddard, 2024](https://svs.gsfc.nasa.gov/14476/)). The older Galactic
Center excess debate supplied the second pressure point: the same gamma-ray
residual can be read as dark matter, unresolved pulsars, diffuse-emission
modeling, or foreground structure depending on the chosen morphology template.
The OPH question is therefore not whether there are "extra photons." It is
whether a source-derived, frozen OPH morphology survives ordinary source,
foreground, instrument, held-out, cross-tracer, and null tests.

**Status:** solved as a morphology-first audit theorem package and simulator
receipt contract. The default claim is diagnostic; candidate promotion requires
frozen source artifacts, count-space response, controls, held-out validation,
cross-tracer tests, and null tests.

Date: 2026-07-08

## Origin

This note records the OPH gamma-morphology audit prompted by the motivating
observations. It asks how OPH should treat gamma-ray dipoles, Galactic Center
or halo residuals, and large-scale gamma morphology tests without fitting
arbitrary residual power.

The answer is a paper-side closure of admissible objects, not a detection
claim. OPH gamma signatures are source-derived morphology claims. They are not
excess-power claims.

## Claim Boundary

The closed paper-side objects are

```math
\mathfrak G_{\gamma,r},
\quad
\mathcal C_{\gamma,T},
\quad
\mathcal R_{\gamma,r},
\quad
\mathcal I_{\gamma,r},
\quad
\mathfrak B_{\partial,r},
\quad
R_\gamma .
```

These objects define admissible OPH gamma morphologies. They do not guarantee a
detected signal. The simulator answers the separate empirical question: whether
the frozen consequence survives foregrounds, source catalogs, masks,
instrument response, likelihood comparison, held-out validation, cross-tracer
tests, and null tests.

The observed model is

```math
I_\gamma(\hat n,E)
=
I_{\rm iso}(E)
+
I_{\rm gal}(\hat n,E)
+
I_{\rm src}(\hat n,E)
+
I_{\rm OPH}(\hat n,E)
+
\epsilon .
```

The OPH component is allowed only when it is a forward projection of a frozen
source artifact:

```math
I_{\rm OPH}
=
\Pi_{\gamma,r}(\mathfrak G_{\gamma,r}),
\qquad
\Pi_{\gamma,r}
=
\mathcal I_{\gamma,r}
\circ
\mathcal L_{\gamma,r}
\circ
\mathcal R_{\gamma,r}
\circ
\mathcal C_{\gamma,r}.
```

A component chosen by fitting gamma residuals is not an OPH prediction.

## Gamma Source Artifact

The proof-carrying gamma artifact is

```math
\mathfrak G_{\gamma,r}
:=
(
id_{\gamma,r},
r,
Q_r,
\mu_r,
\mathrm{route}_r,
\mathfrak P_{A,r},
\mathfrak B_{\partial,r},
\mathfrak A_{\rm astro,r},
\mathcal C_{\gamma,T,r},
\mathcal R_{\gamma,r},
\mathcal L_{\gamma,r},
\mathcal I_{\gamma,r},
\mathcal F_{\rm base,r},
\mathcal F_{\rm alt,r},
\mathcal X_{\rm cross,r},
\mathcal M_{\rm freeze,r},
\mathcal E_{\rm err,r}
).
```

The route is one of:

```text
TRANSPORTED_STRESS
BOUNDARY_RECORD_DIPOLE
BOTH
```

The quotient and source law are mandatory:

```math
Q_r=\Sigma_r/\Gamma_r,
\qquad
\mu_r(q)=Z_r^{-1}m_r(q)e^{-S_r(q)} .
```

Normal form is not a probability law. If the source law, parent artifact,
boundary artifact, response, instrument, foreground registry, likelihood, and
freeze manifest are not declared before comparison, the gamma row is only a
diagnostic map.

## Transported-Stress Route

The transported-stress route starts from a finite covariant parent

```math
\mathfrak P_{A,r}
=
(\mathcal X_r,g_r,e_r,\Pi_r,Z_r,f_r,\mathsf T_r,\mathsf C_r,\mathsf M_r,\mathsf R_r).
```

This parent must emit the anomaly stress tensor \(T_A^{\mu\nu}\). A scalar row
is insufficient. Given the baryonic rest-frame four-velocity

```math
u_b^\mu
=
\frac{J_b^\mu}
{\sqrt{-g_{\alpha\beta}J_b^\alpha J_b^\beta}},
```

the default OPH gamma morphology scalar is

```math
\sigma_{T,r}
=
T_A^{\mu\nu}u_{b,\mu}u_{b,\nu}.
```

On the positive kinetic parent branch this contraction is nonnegative. The
simulator must not replace this theorem path by `abs(rho_A)`. Any use of flux,
anisotropic stress, gradients, smoothing axes, or extended contractions is an
extended branch that needs an `EXTENDED_STRESS_CONTRACTION_RECEIPT`.

## No-Direct-Gamma Rule

Gamma rays are electromagnetic Standard Model radiation. On the neutral anomaly
branch, anomaly packets carry no Standard Model electromagnetic, color, or weak
current unless a separate theorem proves otherwise:

```math
q_z^{\rm EM}=q_z^{\rm color}=q_z^{\rm weak}=0.
```

Therefore

```math
j_{\gamma,\rm direct}[T_A]=0
```

unless an `ANOMALY_EM_CURRENT_RECEIPT` passes. Transported OPH stress can only
modulate ordinary photon-production channels or enter through a boundary-record
projection.

## Photon Response

The response decomposes as

```math
\mathcal R_{\gamma,r}
=
\mathcal R_{\gamma,T,r}
\oplus
\mathcal R_{\gamma,D,r}.
```

For transported stress, the ordinary frozen channel basis is

```text
pi0
IC
brem
Phi
```

where `pi0` is gas/cosmic-ray hadronic emission, `IC` is inverse-Compton
emission, `brem` is bremsstrahlung, and `Phi` is a declared baryonic-potential
or stress-lensing-correlated response. For a channel \(c\),

```math
j_{\gamma,T,r}^{(c)}(x,E)
=
\kappa_{c,r}(E)\,
\bar\sigma_{T,r}(x)\,
E_{c,r}(x,E).
```

The spectral weight \(\kappa_{c,r}(E)\), units, normalization, smoothing, and
ordinary emissivity fields must be frozen before comparison. No arbitrary
residual-shaped spectrum is allowed.

## Boundary-Record Dipole

The boundary-record route lives directly on the observer-facing screen:

```math
\mathfrak B_{\partial,r}
=
(
Q_{\partial,r},
\mu_{\partial,r},
b_{\partial,r},
\Pi_{\ell=1},
\mathbf a_{\partial,r},
D_{\partial,r},
\mathbf d_{\rm OPH,r},
s_{D,r}(E),
\mathcal H_{\partial,r},
\mathcal E_{\partial,r}
).
```

Remove the monopole from the source-side boundary scalar \(b_{\partial,r}\),
then compute the dipole vector

```math
\mathbf a_{\partial,r}
=
\frac{3}{4\pi}
\int_{S^2}
b'_{\partial,r}(\hat n)\hat n\,d\Omega .
```

The dipole axis is

```math
\mathbf d_{\rm OPH,r}
=
\frac{\mathbf a_{\partial,r}}{|\mathbf a_{\partial,r}|}.
```

If the dipole amplitude is zero, the boundary-dipole route is disabled. If the
axis is chosen after inspecting gamma maps, the claim fails closed.

The signed dipole template is

```math
\tau_{D,r}(\hat n,E)
=
s_{D,r}(E)\,\mathbf a_{\partial,r}\cdot\hat n .
```

Because this is signed, predicted counts must remain positive in every
pixel-energy bin.

## Instrument And Likelihood

The simulator must compare binned predicted counts, not unconvolved intensity
maps. The instrument operator contains exposure, PSF, energy dispersion, event
class/type, zenith cuts, time interval, mask, pixelization, energy bins, and
invalid-bin policy:

```math
\mathcal I_{\gamma,r}
=
(
\mathcal E_r,
\mathcal P_r,
\mathcal D_{E,r},
\mathcal M_r,
\mathcal B_r,
\mathcal C_{{\rm event},r},
\mathcal Z_r
).
```

The count model is

```math
\mu_{pe}
=
\mathcal I_{\gamma,r}
\left[
I_{\rm iso}
+
I_{\rm gal}
+
I_{\rm src}
+
\sum_c A_{T,c}\tau_{T,c,r}
+
A_D\tau_{D,r}
\right]_{pe}.
```

The likelihood is binned Poisson:

```math
\ln\mathcal L
=
\sum_{p,e}
\left[
k_{pe}\ln\mu_{pe}
-
\mu_{pe}
-
\ln(k_{pe}!)
\right].
```

For a public Fermi implementation, the data/source release must be pinned. The
Fermi 4FGL-DR4 page describes the 14-year LAT catalog over 50 MeV to 1 TeV
([Fermi LAT 14-year catalog](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/)).
The Fermi background-model page also warns that the Pass 8 Galactic
interstellar model is intended for point-source analysis, so morphology tests
must include foreground alternatives rather than treating one diffuse file as
truth
([Fermi LAT background models](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html)).

## Identifiability

Let \(\tilde\tau_{\rm OPH}\) be the instrument-convolved OPH template and
\(\tilde B_j\) the nuisance templates. With

```math
\langle X,Y\rangle_W
=
\sum_{p,e}\frac{X_{pe}Y_{pe}}{\mu_{0,pe}},
```

project away the nuisance span:

```math
\tilde\tau_\perp
=
(I-P_B)\tilde\tau_{\rm OPH}.
```

The identifiability metric is

```math
\eta_{\rm id}
=
\frac{|\tilde\tau_\perp|_W}
{|\tilde\tau_{\rm OPH}|_W}.
```

If \(\eta_{\rm id}\) is below the declared threshold, the OPH amplitude is not
identifiable under the foreground model. A likelihood improvement alone cannot
promote the claim.

## Promotion Receipt

The gamma promotion receipt is

```math
R_\gamma
=
R_Q
\wedge R_\mu
\wedge R_{\rm route}
\wedge R_{\rm source}
\wedge R_{\rm parent/boundary}
\wedge R_{\mathcal C}
\wedge R_{\rm neutral}
\wedge R_{\mathcal R}
\wedge R_{\mathcal L}
\wedge R_{\mathcal I}
\wedge R_{\rm freeze}
\wedge R_{\rm no\ data}
\wedge R_{\rm positive}
\wedge R_{\rm foreground}
\wedge R_{\rm ident}
\wedge R_{\rm like}
\wedge R_{\rm heldout}
\wedge R_{\rm cross}
\wedge R_{\rm null}.
```

The simulator promotion ladder is:

```text
DIAGNOSTIC_GAMMA_MAP
SOURCE_DERIVED_GAMMA_TEMPLATE
INSTRUMENT_CONVOLVED_GAMMA_TEMPLATE
IDENTIFIABLE_GAMMA_TEMPLATE
LIKELIHOOD_EVALUATED_GAMMA_MORPHOLOGY
CROSS_TRACER_VALIDATED_GAMMA_MORPHOLOGY
OPH_GAMMA_MORPHOLOGY_CANDIDATE
```

No receipt may skip a boundary. The default output is diagnostic.

## Cross-Tracer Tests

For the transported-stress route, held-out correlations must favor the frozen
stress template over ordinary gas, dust, source, diffuse, bulge, bubble,
millisecond-pulsar, and dark-matter alternatives. Relevant tracers include
lensing maps, baryonic-potential maps, gas maps, dust maps, source catalogs,
cluster gas/lensing maps, and large-scale-structure maps.

For the boundary-dipole route, the axis must be stable across energy, time,
masks, event classes, and sky regions, while exposure, ecliptic, solar, lunar,
Earth-limb, scanning-pattern, shuffled-axis, and random-rotation nulls must not
explain the signal.

## Falsifiers

The gamma branch fails closed if any of these occur:

- the source DAG reads gamma residuals, likelihoods, posterior summaries, or
  nuisance fits;
- the route, axis, smoothing, response, or spectrum is chosen after inspecting
  residuals;
- direct anomaly gamma is nonzero without an electromagnetic-current theorem;
- a transported template is generated from scalar rows alone;
- signed templates make any predicted count \(\mu_{pe}\le0\);
- the OPH template is collinear with foreground, source, or dark-matter
  templates;
- random rotations perform comparably;
- held-out energy, time, event-class, mask, or sky-region splits fail;
- a single diffuse model is treated as the only foreground alternative;
- immutable data, mask, source, solver, nuisance, and likelihood hashes are
  missing;
- cross-tracer validation fails but the report still emits
  `OPH_GAMMA_MORPHOLOGY_CANDIDATE`.

## Simulator Contract

Both simulator surfaces use the same fail-closed behavior:

- `oph-physics-sim` compiles the typed gamma morphology receipt surface,
  exports it into measurement packs, checks no-data-use tokens, and refuses
  promotion by default.
- `reverse-engineering-reality/code/particles` mirrors the receipt scaffold so
  paper-stack runs can emit the same source, route, instrument, foreground,
  likelihood, validation, and null-test manifest.

The simulator may consume source artifacts, hash them, compile templates, and
run likelihood/validation tests. It must not discover the OPH gamma template
from gamma residuals.

## Outcome

This note closes at the OPH theorem-package level for admissible gamma
morphology claims. It remains empirically open at the detection level. A future
claim must pass \(R_\gamma\), including frozen source provenance, count-space
instrument response, foreground alternatives, identifiability, held-out
validation, cross-tracer tests, and null tests.
