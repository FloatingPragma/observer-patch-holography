# High-Energy Messenger Coefficient Emission

## Motivating Result

This note entered the queue because high-energy messenger astronomy has started
to isolate source classes that look source-side and hidden rather than simply
bright in ordinary photons. IceCube reported evidence for high-energy neutrinos
from NGC 1068, a gamma-obscured active galaxy, and in 2025 reported evidence
that X-ray-bright active galaxies form a broader neutrino-emitting population
([IceCube, 2025](https://icecube.wisc.edu/news/research/2025/10/evidence-for-neutrino-emission-from-x-ray-emitting-galaxies/)).
At the same time, ultra-high-energy cosmic-ray source inference remains
unsettled: the Telescope Array Amaterasu event arrived with about
$2.4\times 10^{20}$ eV and pointed back toward the Local Void
([University of Utah, 2023](https://attheu.utah.edu/facultystaff/cosmic-ray-2023/)),
while Auger combined analyses favor nearby source populations such as
starburst galaxies over several alternatives
([Pierre Auger Observatory](https://www.auger.org/news/scientific-highlights/321-where-do-ultra-high-energy-cosmic-rays-come-from)).

The OPH question is therefore not whether OPH can fit high-energy neutrino,
cosmic-ray, or gamma-ray maps. The question is whether a finite source-release
artifact can emit the coefficients of a common hidden compact-engine source law
before event coordinates, energies, associations, likelihood residuals, or
diagnostic overlays are read.

**Status:** conditional finite-MaxEnt result plus an open physical source
problem. Existence and uniqueness of coefficients follow under the stated
finite-support assumptions. OPH does not derive the compact-engine loads,
feature map, or astrophysical baseline needed to turn that result into a
source-only physical prediction. The result is not a detected-source claim and
not a UHE event fit.

Date: 2026-07-08

## Claim Boundary

The forbidden shortcut is

```text
UHE events -> fitted eta_A, eta_C, eta_H, eta_AC, eta_AH -> OPH source law
```

The allowed route is

```text
finite OPH source release
-> source-side compact-engine load moments
-> unique MaxEnt coefficients
-> shared source emissivity for neutrinos, cosmic rays, and gammas
-> frozen propagation/detector likelihood
```

OPH does not get a free coefficient for every messenger channel. The source
coefficient must be common at the compact-engine layer; neutrino, cosmic-ray,
and gamma maps are downstream species kernels.

The common-coefficient lock is a declared joint-source hypothesis, not a
consequence of OPH alone. It is useful precisely because the different
messenger data can falsify it after it is frozen.

## What Is Standard, What Is Open, And What OPH Adds

Exponential families, maximum-entropy duality, multimessenger transport, and
Poisson detector likelihoods are standard machinery. The unresolved physics is
the actual source population, acceleration and interaction laws, composition,
propagation, and the derivation of source-side moment constraints.

OPH adds a disciplined observer-like source patch: bounded engine state,
physical ports, source readback, repair/load records, a frozen dependency DAG,
and public receipts proving that target events did not leak upstream. That
makes the source-versus-fit distinction unusually explicit; it does not make
the MaxEnt mathematics or source-identification problem uniquely OPH.

## Source Artifact

For regulator $r$, the finite UHE coefficient artifact is

```math
C_r^{\rm UHE}
=
(
Q_r^{\rm rel},
\mu_r^{\rm rel},
G_r,
U_r,
m_{0,r},
F_r,
L_r^{\rm CE},
N_r,
{\rm Solve}_r,
H_r
).
```

The components are:

- $Q_r^{\rm rel}=\Sigma_r^{\rm rel}/\Gamma_r^{\rm rel}$, with hidden
  gauge-equivalent representatives and demonstrably inert computation
  metadata quotiented out. Physical ports, source-region labels, duty state,
  and repair histories that alter a load remain visible.
- $\mu_r^{\rm rel}(q)=Z_r^{-1}m_r(q)\exp[-S_r(q)]$, or a declared finite
  source MaxEnt law.
- $G_r$, a source-cell partition. A cell has labels such as sky/source region,
  redshift shell, and source class. If those labels come from an external
  survey, the result is source-separated from the target UHE events but is not
  an ab initio OPH population prediction.
- $U_r$, the astrophysical support. A point $u=(g,\theta)$ includes compactness,
  luminosity, obscuration, baryon loading, photon field, magnetic environment,
  maximum rigidity, duty state, and opacity.
- $m_{0,r}(u)$, a frozen baseline probability mass with full support, fixed
  before UHE comparison and normalized by
  $\sum_{u\in U_{r,g}}m_{0,r}(u)=1$.
- $F_r(u)$, a finite feature map. The first UHE lane uses
  $\widehat A_r(g)$, $C_{\rm compact}(\theta)$, $H_{\rm hidden}(\theta)$,
  $\widehat A_r C_{\rm compact}$, and $\widehat A_r H_{\rm hidden}$.
- $L_r^{\rm CE}$, a compact-engine source-load observable.
- $N_r$, a unit and normalization ledger.
- ${\rm Solve}_r$, a deterministic convex-dual solver.
- $H_r$, hashes and receipts for the source DAG, no-data-use manifest, feature
  map, baseline, load observable, solver, coefficients, refinement, and
  nonclaim ledger.

## Compact-Engine Source Loads

The compact-engine source-load observable is not a fitted UHE coefficient. It is
source-side scalar-slot bookkeeping:

```math
L_{\alpha,r}^{\rm CE}(q;g)
=
\Pi_{g,r}
\sum_{C\subset g}
V_{C,r}^{\rm phys}(q)\,
\omega_{\alpha,r}(C;q)\,
S_{\nu,r,C}(q).
```

The gate $\omega_{\alpha,r}$ must be quotient-visible. Allowed channel gates
include compact-record activation, transported stress, compact-engine support,
hidden or obscured-source state, the interaction gates $AC$ and $AH$, pion or
baryon loading, gamma opacity, and maximum rigidity. If any gate reads event
coordinates, energies, associations, residual maps, or likelihood values, the
coefficient is not source-only.

The source-emitted target moment is

```math
c_{\alpha,r}(g)
=
N_{\alpha,r}\,
\mathbb E_{\mu_r^{\rm rel}}
\left[
L_{\alpha,r}^{\rm CE}(Q;g)
\right].
```

This is a modeling bridge and carries no OPH derivation status. The normalization
$N_{\alpha,r}$ must state units and map the load vector to the same
dimensionless feature-moment coordinates as $F_r$; the resulting
$c_r(g)$ must lie in $\operatorname{relint}\operatorname{conv}F_r(U_{r,g})$.
Changing $N$, the gates, or the source-cell definition changes the physical
hypothesis and requires a separately frozen artifact.

No simulator output, likelihood value, observed event target, tuned diagnostic,
or human-selected event pattern may enter this source law or constraint ledger.

## Coefficient-Emission Theorem

**Statement.** Fix a source cell $g$ and finite support $U_{r,g}$. Assume:

```text
full support: m0_r(u) > 0 on U_{r,g}
minimal features: no nonzero v makes v.F_r constant on U_{r,g}
interior moment: c_r(g) lies in the relative interior of conv F_r(U_{r,g})
source-only loads: c_r(g) is built only from quotient-visible source loads
no target leakage: no UHE target data path reaches the source ledger or solver
```

Define

```math
A_{r,g}(\eta)
=
\log
\sum_{u\in U_{r,g}}
m_{0,r}(u)\exp[\eta\cdot F_r(u)].
```

Then there is a unique coefficient vector $\eta_r(g)$ satisfying

```math
\nabla_\eta A_{r,g}(\eta_r(g))=c_r(g),
```

equivalently minimizing

```math
A_{r,g}(\eta)-\eta\cdot c_r(g).
```

It defines the source probability mass function

```math
p_{\eta,r,g}(u)
=
m_{0,r}(u)
\exp[
\eta_r(g)\cdot F_r(u)-A_{r,g}(\eta_r(g))
],
\qquad u\in U_{r,g}.
```

and the OPH source weight

```math
W_{\rm OPH}(u|g)
=
\exp[
\eta_r(g)\cdot F_r(u)-A_{r,g}(\eta_r(g))
].
```

**Proof.** The gradient of $A$ is the expectation of $F$ under the exponential
family. The Hessian is the covariance of $F$. Minimality makes the covariance
strictly positive definite on coefficient directions. Interior target moments
give existence of the dual optimum. Strict convexity gives uniqueness. The
same optimum is the relative-entropy MaxEnt law relative to $m_{0,r}$. $\square$

## Closed-Form Corollaries

For a binary source gate $B$,

```math
\eta_B
=
\operatorname{logit}(p_{\rm OPH})
-
\operatorname{logit}(p_0)
=
\log
\frac{p_{\rm OPH}(1-p_0)}
{p_0(1-p_{\rm OPH})}.
```

For a Poisson opportunity count $K$,

```math
\eta_K
=
\log(\lambda_{\rm OPH}/\lambda_0).
```

This Poisson relation is a separate count-family example on infinite support;
it is not a corollary of the finite-support theorem without an explicit
truncation-and-limit argument.

For small moment displacement $\delta c$ around the baseline,

```math
\eta
=
\Sigma_0^{-1}\delta c
+
O(|\delta c|^2),
```

where $\Sigma_0$ is the baseline covariance of the feature map.

## Quotient Theorem

**Statement.** If $S_{\nu,r,C}$, $\omega_{\alpha,r}$,
$V_{C,r}^{\rm phys}$, and $\Pi_{g,r}$ are quotient-visible and fixed before UHE
comparison, then $L_{\alpha,r}^{\rm CE}$ descends to $Q_r^{\rm rel}$ and is
invariant under hidden representatives.

**Proof.** Equivalent presentations in $\Sigma_r^{\rm rel}$ share the same
quotient state. Each factor in $L_{\alpha,r}^{\rm CE}$ is a function of that
quotient or a predeclared physical normalization. The product and finite sum
therefore assign the same load to equivalent representatives. $\square$

## Source-Only Classifier

Let $G_{\rm dep}$ be the dependency graph for the coefficient artifact.

```text
SOURCE_ONLY_OPH_COEFFICIENT:
  no path from target UHE data into Q, mu, m0, F, L, c, Solve, or eta

FITTED_OPH_COEFFICIENT:
  a declared target-data path exists

INVALIDATED_COEFFICIENT_DAG:
  a hidden target-data path exists while the artifact claims source-only
```

Forbidden upstream objects include event coordinates, energies, association
failures, catalog matches chosen after seeing events, likelihood values,
posterior summaries, diagnostic overlays, residual maps, and human-selected
event patterns.

## Coefficients To Messenger Emissivity

Let $N_g$ be the declared total compact-engine number density in source cell
$g$, and identify the physical baseline density with the MaxEnt base mass by

```math
\rho_{\rm HCE}^{0}(u|g)=N_gm_{0,r}(u|g),
\qquad \sum_{u\in U_{r,g}}m_{0,r}(u|g)=1.
```

The shared compact-engine density is then

```math
\rho_{\rm HCE}^{\rm OPH}(u|g)
=
\rho_{\rm HCE}^{0}(u|g)
\exp[\eta_r(g)\cdot F_r(u)-A_{r,g}(\eta_r(g))].
```

It sums to $N_g$. If a different $\rho^0$ is used, its normalized shape
must replace $m_0$ in the partition function; multiplying an unrelated
baseline by the MaxEnt weight would not preserve normalization.

Species emissivities are downstream:

```math
S_\nu = \rho_{\rm HCE}^{\rm OPH} Q_\nu,
\qquad
S_{\rm CR} = \rho_{\rm HCE}^{\rm OPH}\pi(Z|\theta) Q_{\rm CR},
\qquad
S_\gamma = \rho_{\rm HCE}^{\rm OPH} Q_\gamma.
```

An observed intensity has the form

```math
\lambda_a(y)
=
b_a(y)
+
\sum_Z\int
S_a(E_s,Z,u)
P_a^{\rm prop}(y_{\rm true}|E_s,Z,u)
R_a(y|y_{\rm true})
\,dE_s\,du.
```

For neutrinos or photons the discrete composition sum is replaced by the
appropriate species sum. Every source yield and kernel must state its measure
and units so that $\lambda_a$ has the detector intensity units used by the
likelihood.

The coefficient is common. Channel-specific differences enter through species
yield, propagation, and detector kernels.

## Refinement Stability

On a common coefficient chart, suppose two adjacent regulators have pushed
source laws, loads, and features with normed defects $\epsilon_\mu$,
$\epsilon_L$, and $\epsilon_F$. If both target moments remain in a
common compact subset of the moment-polytope interior and the smallest Hessian
eigenvalue there is bounded below by $\kappa>0$, the inverse-function
estimate has the conditional form

```math
\|\eta_r-\eta_s\|
\le
\kappa^{-1}
O(\epsilon_\mu+\epsilon_L+\epsilon_F).
```

Vanishing adjacent defects alone does not imply convergence: increments can go
to zero while their sum diverges. A stable refinement limit additionally needs
summable defects, a direct comparison to a declared common limit, or another
Cauchy/compactness argument. The norms, constants hidden by $O(\cdot)$, and
coarse-map identifications must be part of the receipt.

## Simulator Contract

The simulator command is

```bash
oph-fpe uhe-emit-coefficients \
  --features-json @features.json \
  --baseline-json @baseline.json \
  --target-json @target.json \
  --source-dag-json @source_dag.json \
  --out runs/uhe/coefficient_emission
```

`--target-json` is the current CLI's legacy name for the source-side moment
constraint $c$; it must not contain target messenger events or fitted
summaries. A future interface should call this `--source-moments-json`. If the
command uses its planted default coefficients or generated moments, the output
is a synthetic algebra/visualization fixture and may not claim a source receipt.

It must emit `uhe_coefficient_emission_report.json` and
`uhe_coefficient_emission_report.md` with:

```text
artifact_type: UHE_COEFFICIENT_EMISSION_RECEIPT
claim_tier: SOURCE_ONLY | CONDITIONAL | FITTED | INVALIDATED
input_hashes
moment_targets
solver convergence and residual
coefficients
source-law probabilities
readiness_gates
nonclaim ledger
```

Required receipts:

```text
BASELINE_FULL_SUPPORT
FEATURE_MINIMALITY
MOMENT_INTERIOR
SOURCE_LOAD_QUOTIENT_VISIBLE
NO_UHE_DATA_USE
REFINEMENT_COMPATIBILITY
COEFFICIENT_SOLVE_CONVERGED
COMMON_SOURCE_LOCK
```

Hard failures:

```text
moment outside convex hull
feature non-minimal unless redundant features are removed
source-only claim with target-data dependency
separate learned coefficient maps for neutrinos, cosmic rays, and gammas
hidden event coordinates, energies, residuals, likelihoods, or overlays upstream
```

The numerical solver is a damped Newton solve on the finite log partition. Its
regression tests must include the binary logit corollary, Poisson log-ratio
corollary, small-signal covariance limit, no-data-use trap, feature-rank trap,
moment-polytope trap, refinement-stability report, and common-source lock.

## What OPH Answers

OPH does not uniquely solve the ordinary astrophysical source-identification
problem. Regular physics has to determine source classes, acceleration,
composition, magnetic deflection, propagation losses, detector response, and
likelihood support.

What OPH adds is a strict source-only coefficient rule for observer-like source
patches. A high-energy messenger source claim is admissible only when the
finite source release has local state, visible ports or boundaries, readback
records, repair/load observables, feedback or source-law selection, and public
receipts. The coefficient is emitted from that source artifact by finite MaxEnt
duality; it is not learned from the events it later tries to predict.

The claim labels are:

```text
SOURCE_ONLY_COEFFICIENT_EMITTED:
  all coefficient-emission receipts pass

CONDITIONAL_SOURCE_MODEL:
  the finite solver works but one or more source, refinement, or no-leak
  receipts are missing

FITTED_OPH_COEFFICIENT:
  target data are deliberately used, so the row is phenomenology only

INVALIDATED_COEFFICIENT_DAG:
  hidden target leakage exists under a source-only label
```
