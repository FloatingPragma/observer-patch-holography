# Gamma-Ray Morphology Claims In OPH

## Motivating Result

This note entered the queue after a 13-year Fermi-LAT analysis reported a
diffuse gamma-ray-background dipole of roughly $6.5$--$7\%$ over
$2.74<E\le115.5\,\mathrm{GeV}$, with an axis close to the Pierre Auger UHECR
dipole rather than the CMB dipole
([Kashlinsky, Atrio-Barandela, and Korotkov, 2024](https://arxiv.org/abs/2401.04564)).
The older Galactic Center excess debate supplies a separate pressure point:
there a residual can be read as dark matter, unresolved pulsars,
diffuse-emission modeling, or foreground structure depending on the chosen
morphology template. The diffuse-background dipole and Galactic Center excess
are different datasets and must never be merged into one claimed signal.
The OPH question is therefore not whether there are "extra photons." It is
whether a source-derived, frozen OPH morphology survives ordinary source,
foreground, instrument, held-out, cross-tracer, and null tests.

**Status:** conditional morphology-audit specification. The note defines a
fail-closed route from an assumed OPH source artifact to a testable count-space
template; it does not derive a gamma-ray morphology from the recovered OPH core
and does not report a detection. Candidate promotion requires frozen source
artifacts, count-space response, controls, held-out validation, cross-tracer
tests, and null tests.

Date: 2026-07-08

## Origin

This note records the OPH gamma-morphology audit prompted by the motivating
observations. It asks how OPH should treat gamma-ray dipoles, Galactic Center
or halo residuals, and large-scale gamma morphology tests without fitting
arbitrary residual power.

The answer is a paper-side closure of admissible objects, not a detection
claim. OPH gamma signatures are source-derived morphology claims. They are not
excess-power claims.

## What Is Standard, What Is Open, And What OPH Adds

Standard high-energy astrophysics already supplies photon-production channels,
radiative transport, catalog and diffuse-source models, instrument response,
Poisson likelihoods, and template/null-test methodology. The open physical
question is whether OPH supplies a source law that fixes a nonzero morphology
and spectrum before gamma-ray data are inspected.

OPH's present contribution is narrower and useful: it makes the candidate an
observer-like self-reading patch with bounded state, declared ports, source
readback, immutable records, repair/refinement checks, and a public promotion
receipt. This architecture makes source-to-map provenance unusually explicit;
it is not a unique solution to gamma-ray source identification.

Throughout this note:

- **assumed** means a provisional source or response ansatz;
- **derived within the declared model** means algebra following from those
  assumptions;
- **software-checked** means a finite receipt or invariant was recomputed; and
- **empirically supported** requires held-out photon data and alternatives.

## Claim Boundary

The closed paper-side objects are

$$
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
$$

These objects define admissible OPH gamma morphologies. They do not guarantee a
detected signal. The simulator answers the separate empirical question: whether
the frozen consequence survives foregrounds, source catalogs, masks,
instrument response, likelihood comparison, held-out validation, cross-tracer
tests, and null tests.

The pre-instrument expected sky intensity is

$$
I_\gamma^{\rm sky}(\hat n,E)
=
I_{\rm iso}(E)
+
I_{\rm gal}(\hat n,E)
+
I_{\rm src}(\hat n,E)
+
I_{\rm OPH}^{\rm sky}(\hat n,E).
$$

Foreground-model discrepancy may be represented by separately declared
nuisance fields or priors. It is not an additive random count after a Poisson
count model has already been specified.

The OPH component is allowed only when it is a forward projection of a frozen
source artifact:

$$
I_{\rm OPH}^{\rm sky}
=
\Pi_{\gamma,r}(\mathfrak G_{\gamma,r}),
\qquad
\Pi_{\gamma,r}
=
\mathcal L_{\gamma,r}
\circ
\mathcal R_{\gamma,r}
\circ
\mathcal C_{\gamma,r}.
$$

A component chosen by fitting gamma residuals is not an OPH prediction.

## Gamma Source Artifact

The proof-carrying gamma artifact is

$$
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
$$

The route is one of:

```text
TRANSPORTED_STRESS
BOUNDARY_RECORD_DIPOLE
BOTH
```

The quotient and source law are mandatory:

$$
Q_r=\Sigma_r/\Gamma_r,
\qquad
\mu_r(q)=Z_r^{-1}m_r(q)e^{-S_r(q)} .
$$

Here $\Sigma_r$ and the action of $\Gamma_r$ must be declared, and only
genuine gauge or presentation redundancies may be quotiented out; physical
ports and currents remain visible. On a finite state space, $m_r\ge0$,
$S_r$ is dimensionless, and
$0<Z_r=\sum_qm_r(q)e^{-S_r(q)}<\infty$. On a continuous space the base
measure, domain, and integrability condition must be stated. A quantum source
branch needs a positive trace-one density operator rather than this classical
Gibbs notation.

Normal form is not a probability law. If the source law, parent artifact,
boundary artifact, response, instrument, foreground registry, likelihood, and
freeze manifest are not declared before comparison, the gamma row is only a
diagnostic map.

## Transported-Stress Route

The transported-stress route starts from a finite covariant parent

$$
\mathfrak P_{A,r}
=
(\mathcal X_r,g_r,e_r,\Pi_r,Z_r,f_r,\mathsf T_r,\mathsf C_r,\mathsf M_r,\mathsf R_r).
$$

This parent must emit the anomaly stress tensor $T_A^{\mu\nu}$. A scalar row
is insufficient. Given the baryonic rest-frame four-velocity

$$
u_b^\mu
=
\frac{J_b^\mu}
{\sqrt{-g_{\alpha\beta}J_b^\alpha J_b^\beta}},
$$

the default OPH gamma morphology scalar is

$$
\sigma_{T,r}
=
T_A^{\mu\nu}u_{b,\mu}u_{b,\nu}.
$$

This contraction is nonnegative only on a branch that supplies a
future-timelike $u_b$ and an energy condition such as
$T_A^{\mu\nu}u_\mu u_\nu\ge0$. That condition is an assumption to be
checked, not a consequence of the contraction notation. The simulator must
not replace it by `abs(rho_A)`. Any use of flux,
anisotropic stress, gradients, smoothing axes, or extended contractions is an
extended branch that needs an `EXTENDED_STRESS_CONTRACTION_RECEIPT`.

## No-Direct-Gamma Rule

Gamma rays are electromagnetic Standard Model radiation. On the declared
neutral, electromagnetically decoupled anomaly branch, anomaly packets carry no
Standard Model electromagnetic, color, or weak current:

$$
q_z^{\rm EM}=q_z^{\rm color}=q_z^{\rm weak}=0.
$$

Neutral charge alone is insufficient: a neutral composite can still couple
electromagnetically through moments or polarizability. Therefore the following
zero-current statement additionally requires a decoupling/current theorem:

$$
j_{\gamma,\rm direct}[T_A]=0
$$

unless an `ANOMALY_EM_CURRENT_RECEIPT` passes. Transported OPH stress can only
modulate ordinary photon-production channels or enter through a boundary-record
projection.

## Photon Response

The response decomposes as

$$
\mathcal R_{\gamma,r}
=
\mathcal R_{\gamma,T,r}
\oplus
\mathcal R_{\gamma,D,r}.
$$

For transported stress, the ordinary frozen channel basis is

```text
pi0
IC
brem
Phi
```

where `pi0` is gas/cosmic-ray hadronic emission, `IC` is inverse-Compton
emission, `brem` is bremsstrahlung, and `Phi` is a declared baryonic-potential
or stress-lensing-correlated response. For a channel $c$, a dimensionally
explicit **illustrative linear-response ansatz** is

$$
j_{\gamma,T,r}^{(c)}(x,E)
=
j_{0,c,r}(x,E)
\left[
1+\kappa_{c,r}(E)
\frac{\bar\sigma_{T,r}(x)}{\sigma_{0,r}}
\right].
$$

Here $\kappa$ is dimensionless, $\sigma_0$ is a positive reference
stress in the same units as $\bar\sigma_T$, and positivity of the bracket
is required. A different response is admissible if its units, normalization,
smoothing, and ordinary emissivity fields are frozen before comparison. This
response is assumed until derived from a microscopic OPH-matter coupling; no
arbitrary residual-shaped spectrum is allowed.

## Boundary-Record Dipole

The boundary-record route lives directly on the observer-facing screen:

$$
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
$$

Remove the monopole from the source-side boundary scalar $b_{\partial,r}$,
then compute the dipole vector

$$
\mathbf a_{\partial,r}
=
\frac{3}{4\pi}
\int_{S^2}
b'_{\partial,r}(\hat n)\hat n\,d\Omega .
$$

The dipole axis is

$$
\mathbf d_{\rm OPH,r}
=
\frac{\mathbf a_{\partial,r}}{|\mathbf a_{\partial,r}|}.
$$

If the dipole amplitude is zero or fails a predeclared source-side numerical
threshold, the axis is undefined and the boundary-dipole route is disabled.
The threshold, mask, mode-coupling correction, and axis uncertainty must be
frozen before comparison. If the axis is chosen after inspecting gamma maps,
the claim fails closed.

The signed dipole template is

$$
\tau_{D,r}(\hat n,E)
=
s_{D,r}(E)\,\mathbf a_{\partial,r}\cdot\hat n .
$$

Because this is signed, predicted counts must remain positive in every
pixel-energy bin.

## Instrument And Likelihood

The simulator must compare binned predicted counts, not unconvolved intensity
maps. The instrument operator contains exposure, PSF, energy dispersion, event
class/type, zenith cuts, time interval, mask, pixelization, energy bins, and
invalid-bin policy:

$$
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
$$

The count model is

$$
\mu_{pe}
=
\mathcal I_{\gamma,r}
\left[
I_{\rm iso}
+
I_{\rm gal}
+
I_{\rm src}+I_{\rm OPH}^{\rm sky}
\right]_{pe}.
$$

Two amplitude conventions are allowed, and they must not be mixed:

1. **Source-amplitude prediction.** The frozen source projection supplies
   $I_{\rm OPH}^{\rm sky}$ with physical intensity units, including the
   normalization of $\bar\sigma_T$, $\mathbf a_\partial$, and
   $s_D(E)$. No extra fitted amplitude is applied; equivalently every
   multiplier is fixed to one.
2. **Morphology-only test.** Each raw template is divided by a preregistered
   nonzero norm $N_j$, giving $\widehat\tau_j=\tau_j/N_j$, and
   $I_{\rm OPH}^{\rm sky}=\sum_jA_j\widehat\tau_j$. The fitted $A_j$
   carries the intensity units. This can test a frozen shape, axis, and spectrum
   family, but it is not a source-derived amplitude or absolute-intensity
   prediction.

The template norm, sign convention, energy dependence, and positivity range are
frozen before data comparison. Allowing $A_D$, $s_D$, and
$|\mathbf a_\partial|$ all to float would be an unidentifiable
rescaling, not an OPH prediction.

The likelihood is binned Poisson:

$$
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
$$

For a public Fermi implementation, the data/source release must be pinned. The
Fermi 4FGL-DR4 page describes the 14-year LAT catalog over 50 MeV to 1 TeV
([Fermi LAT 14-year catalog](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/)).
The Fermi background-model page also warns that the Pass 8 Galactic
interstellar model is intended for point-source analysis, so morphology tests
must include foreground alternatives rather than treating one diffuse file as
truth
([Fermi LAT background models](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html)).

## Identifiability

Let $\tilde\tau_{\rm OPH}$ be the instrument-convolved OPH template and
let the columns of $B$ span the local nuisance tangent space at a strictly
positive baseline $\mu_0$. With

$$
\langle X,Y\rangle_W
=
\sum_{p,e}\frac{X_{pe}Y_{pe}}{\mu_{0,pe}},
$$

let $P_B$ be the $W$-orthogonal projector onto that span and project it
away:

$$
\tilde\tau_\perp
=
(I-P_B)\tilde\tau_{\rm OPH}.
$$

The identifiability metric is

$$
\eta_{\rm id}
=
\frac{\|\tilde\tau_\perp\|_W}
{\|\tilde\tau_{\rm OPH}\|_W}.
$$

This ratio is defined only when $\|\tilde\tau_{\rm OPH}\|_W>0$; a zero
template disables the route rather than earning an identifiability score.

If $\eta_{\rm id}$ is below the declared threshold, the OPH amplitude is not
identifiable under the foreground model. A likelihood improvement alone cannot
promote the claim.

## Promotion Receipt

Let $U_T$ and $U_D$ indicate whether the declared route uses transported
stress and boundary dipole, respectively. Define

$$
R_T=(\neg U_T)\vee
\left(R_{\rm parent}\wedge R_{\mathcal C}\wedge R_{\rm neutral}\right),
\qquad
R_D=(\neg U_D)\vee R_{\rm boundary}.
$$

Thus a boundary-only model does not need a neutral-anomaly or stress-parent
receipt, while a transported-stress model does; `BOTH` must pass both branches.
The route declaration itself must be one of the three allowed values and agree
with $(U_T,U_D)$.

The gamma promotion receipt is

$$
R_\gamma
=
R_Q
\wedge R_\mu
\wedge R_{\rm route}
\wedge R_{\rm source}
\wedge R_T
\wedge R_D
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
$$

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
- signed templates make any predicted count $\mu_{pe}\le0$;
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

This note supplies a conditional audit specification for admissible gamma
morphology claims. It remains open both at the OPH source-derivation level and
at the empirical detection level. A future claim must first derive or
explicitly assume a nonzero source morphology and then pass $R_\gamma$,
including frozen source provenance, count-space instrument response,
foreground alternatives, identifiability, held-out validation, cross-tracer
tests, and null tests.
