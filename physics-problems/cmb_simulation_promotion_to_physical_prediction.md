# CMB Simulation Promotion to a Physical Prediction

## Motivating Result

This note entered the queue after the public OPH mini-universe simulator began
showing an all-sky observer-boundary record view that was visually reminiscent
of a CMB map:
[simulation.floatingpragma.io](https://simulation.floatingpragma.io). The
relevant external precision surface is the real CMB likelihood problem, for
example Planck's final cosmological-parameter analysis:
[Planck 2018 VI](https://www.aanda.org/articles/aa/full_html/2020/09/aa33910-18/aa33910-18.html).

**Status:** promotion contract specified; physical CMB prediction unresolved.
The public simulator is an end-to-end visual representation of OPH. Its current
CMB-facing result is a visual/spectrum diagnostic, not an OPH-derived physical
CMB prediction. Provisional geometry, scale, or transfer assumptions are useful
for the demonstration only when they are labeled as assumptions.

Date: 2026-07-08
Audit revision: 2026-07-11

## Problem, Standard Physics, and OPH Contribution

**Standard physics.** A physical CMB calculation starts with a statistically
defined primordial source, evolves perturbations on a specified cosmological
background through a Boltzmann/recombination model, predicts temperature and
polarization observables, applies beams, masks, foreground and noise models,
and evaluates a frozen likelihood. A visually similar map is not evidence that
this chain has been reproduced.

**Unresolved OPH target.** The unresolved question is whether a source law,
physical scale, common primordial surface, and geometry can be derived from the
OPH quotient carrier and then produce competitive frozen TT, TE, EE, lensing,
and cross-probe predictions. This note supplies none of those missing physical
derivations by itself.

**OPH contribution.** OPH contributes a precise claim ledger for its bounded
self-reading software patches: local state, overlap ports, readback, records,
repair moves, and public evidence receipts. The ledger makes the origin of each
visual or numerical field explicit and blocks a diagnostic screen record from
being relabeled as a physical sky prediction. Provenance ledgers, blinding, and
frozen likelihoods are not unique to OPH; the OPH-specific feature is their
attachment to the observer-patch quotient and repair/readback structure.

For the present visual demonstration, the following are assumptions rather
than recovered results:

- a rendered boundary scalar may be displayed on an all-sky projection;
- an imported, frozen FLRW/CAMB packet may be used for a labeled comparison;
- a screen multipole or color scale may be chosen for visualization, but it is
  not thereby the observed CMB multipole or thermodynamic temperature;
- no visual match supplies a primordial-source, geometry, recombination,
  foreground, or likelihood receipt.

## Claim Ladder

The OPH implementation is a promotion ladder:

| Tier | Claim | Required Meaning |
|---|---|---|
| `VISUAL_DIAGNOSTIC` | the simulator emits an all-sky boundary-scalar field | visual evidence about the simulator only |
| `SPECTRUM_DIAGNOSTIC` | reproducible proxy spectrum or residual diagnostics exist | measurement-facing comparison, not a source derivation |
| `SOURCE_ONLY_FINITE_ARTIFACT` | the field is generated from frozen OPH source artifacts | no CMB, BAO, SN, WL, RSD, cluster, posterior, residual, or likelihood path enters the source DAG |
| `CONDITIONAL_PHYSICAL_CMB_SOURCE` | the source is promoted through a physical scale/freezeout bridge using imported frozen FLRW geometry | an externally conditioned physical pipeline, not an OPH-native geometry derivation |
| `OPH_NATIVE_PHYSICAL_CMB_SOURCE` | the same source is promoted through quotient-derived `CosmoGeomRead_r` | native OPH cosmological geometry |
| `LIKELIHOOD_EVALUATED_PHYSICAL_CMB_PREDICTION` | all required parent gates and matching source/solver/likelihood hashes pass | a likelihood-evaluated physical CMB prediction under the declared model |

The simulator can supply finite witnesses, values, residuals, hashes, and
failure receipts. The paper stack supplies the semantic types that decide what
those values are allowed to mean. Passing the final tier would test a declared
OPH cosmological branch; it would not make that branch uniquely correct without
comparison to standard cosmology and alternative models.

For a scalar diagnostic $X_r(\hat n)$, the reproducible angular diagnostic is

```math
a^{X,r}_{\ell m}=\int_{S^2}X_r(\hat n)Y^*_{\ell m}(\hat n)\,d\Omega,
\qquad
C^{XX,r}_\ell=\frac{1}{2\ell+1}\sum_m|a^{X,r}_{\ell m}|^2.
```

At diagnostic tier, $X_r$ is a boundary scalar with declared units and is not
silently identified with $\Delta T/T$. A physical CMB tier must separately define
thermodynamic temperature and Stokes $Q/U$, the TT/TE/EE spectra, lensing,
beam and mask conventions, foregrounds, noise, and the evaluated multipole
range.

## No-Visual-Shortcut Validation Invariant

**Invariant 1** (no visual shortcut to physical CMB prediction). A simulator
output that visually resembles the CMB, or even fits a CMB-like spectrum, is not
a physical CMB prediction unless its terminal prediction artifact is reached
through a typed dependency DAG whose source, geometry, transfer, no-data-use,
freeze, and likelihood predicates all verify from their own artifacts and
parent hashes.

*Conditional verification argument.* Let every artifact node carry its value,
units, semantic type, parent list, data ledger, code/build hash, theorem hash,
witness hash, and error envelope. A terminal
`PHYSICAL_CMB_PREDICTION_RECEIPT` verifies by
recomputing its predicate from those fields, not by trusting a caller-supplied
Boolean. Induct over the dependency DAG. Each passed node implies its parent
predicates. Therefore a passed terminal prediction implies that the source
artifact was frozen before evaluation, that forbidden observational data did
not enter the source path, that scale and freezeout were physically typed, that
the Boltzmann transfer was fixed, and that the likelihood rule was declared in
advance. Visual similarity has no edge in this DAG, so it cannot promote the
claim.

This is a software-validation result, not a CMB-physics derivation. It assumes
an acyclic and complete dependency graph, a trusted validator, authenticated
artifacts, a hermetic or otherwise audited build, and no undeclared human,
network, cache, or environment side channel. Hashes protect artifact identity;
they do not prove that a model is physically correct.

## Finite Source Object

At regulator $r$, the CMB branch object is:

```math
\mathfrak C_r=
(\mathcal F_r,\Gamma_r,Q_r,\mu_r,\mathfrak S^{\rm scr}_r,
\mathfrak S^{\rm prim}_r,\mathcal S_{{\rm OPH},r},
\mathfrak I_{\star,r},\mathfrak B_r,\mathcal L_D,\mathsf{DAG}_r).
```

Here $\mathcal F_r$ is the early OPH screen federation: bounded
observer-like patches with local state, ports or overlap boundaries, readback,
records, feedback or repair moves, and public evidence receipts. The quotient
$Q_r=\Sigma_r/\Gamma_r$ removes only transformations demonstrated to
preserve the source law, admissible transitions, public readouts, and receipt
validators. These may include gauge representatives and inert mesh, worker,
retry, or port *names*. A physical boundary, port placement, carrier variable,
or schedule that changes a future readout is not a quotient artifact.

The source law is not the settled state. For a finite quotient it may be written
as

```math
\mu_r(q)=Z_r^{-1}m_r(q)e^{-S_r(q)}.
```

Here $m_r(q)\ge0$, $S_r(q)$ is dimensionless, and

```math
Z_r=\sum_{q\in Q_r}m_r(q)e^{-S_r(q)}\in(0,\infty).
```

For a continuous quotient, the sum must be replaced by a declared base measure
and integral. Normal forms classify representatives; they do not select the
probability law. Indeed, almost any strictly positive finite distribution can
be written in this Gibbs form. The expression becomes predictive only when
$m_r$ and $S_r$ are independently derived or preregistered rather than
chosen to reproduce a target spectrum.

The simulator therefore needs `source/quotient_ensemble.json` with quotient
schema hashes, base measure, action terms, zero-mode policy, sampler hash,
detailed-balance or stationarity checks, and partition-invariance receipts.
Sampler stationarity shows that the declared law was sampled; it does not show
that the declared law is the physical primordial source.

## Gate Package

Define the route-dependent gates

```math
R_{\rm geom}:=R_{\rm FLRW\ import}\vee R_{\rm CosmoGeomRead},
\qquad
R_{\rm dark}:=(\mathrm{dark\ continuation}=\mathrm{OFF})
\vee(R_{\rm parent}\wedge R_{\rm kernel}).
```

The compact readiness predicate is then

```math
\begin{aligned}
\mathsf{CMB\_READY}_r={}&
R_Q \wedge R_\mu \wedge R_{\rm DAG} \wedge R_{\rm no\ data}
\wedge R_q \wedge R_K \wedge R_A \wedge R_{\rm screen}\\
&\wedge R_{\rm coh} \wedge R_{\rm radial} \wedge R_{\rm prim}
\wedge R_{\rm scale} \wedge R_{\Sigma_\star} \wedge R_{\rm init}\\
&\wedge R_{\rm geom}\wedge R_{\rm dark}
\wedge R_{\rm transfer} \wedge R_{\rm freeze}
\wedge R_{\rm like} \wedge R_{\rm fals}.
\end{aligned}
```

The symbols are gate names, not additional derived equations. Their detailed
artifact schemas and thresholds live in
[`oph_cosmology_finite_source_cmb_program.md`](../cosmology/oph_cosmology_finite_source_cmb_program.md)
and
[`oph_cosmology_data_likelihood_contracts.md`](../cosmology/oph_cosmology_data_likelihood_contracts.md).
The geometry route is a typed alternative: the conditional route requires an
imported-FLRW receipt, while the native route requires `CosmoGeomRead_r`; a run
does not need to satisfy both mutually exclusive origins. An
`OPH_NATIVE_PHYSICAL` label specifically requires $R_{\rm CosmoGeomRead}$,
not merely the disjunction.

The definition of $R_{\rm dark}$ makes the optionality explicit: parent and
kernel receipts are bypassed only when the declared run has
`dark_continuation: OFF`. If the claim includes OPH dark/anomaly cosmology,
both become mandatory.

The important semantic split is:

- `CONDITIONAL_PHYSICAL` may use an imported, frozen FLRW geometry packet while
  still testing OPH-native source artifacts through physical scale, clock,
  freezeout, and transfer gates.
- `OPH_NATIVE_PHYSICAL` additionally requires quotient-derived
  `CosmoGeomRead_r`, source embedding, and geometry-origin receipts.

Imported FLRW geometry may exercise the physical pipeline. It must not be
renamed as OPH-native geometry.

## No-Data-Use Firewall

A source artifact is certified source-only within the declared build boundary
iff no transitive parent, cache, config,
environment variable, dependency version, random stream, network read, file
read, model-selection rule, or human branch choice has a path from forbidden
observational data. Forbidden data include CMB maps, spectra, residuals,
likelihoods, and posteriors; BAO; supernovae; weak lensing; RSD; SPARC;
clusters; $H_0/S_8$ target values; diagnostic overlays; and residual plots.

This forbids common quiet shortcuts:

- choosing the scale bridge by CMB peak alignment;
- fitting $A_\zeta$, $q_{\rm IR}$, $\ell_{\rm IR}$, or kernels from a
  Planck residual;
- using a cached posterior or trained emulator that saw the target spectra;
- exposing a human branch choice to CMB overlays before source freeze;
- silently importing LambdaCDM best-fit parameters as OPH source parameters.

The firewall must also name allowed inputs: mathematical constants, declared
laboratory calibrations, external background assumptions, training data, and a
final held-out likelihood set. Imported best-fit FLRW parameters make a result
externally conditioned; freezing them prevents post hoc tuning but does not turn
them into OPH-derived parameters. A provenance receipt certifies the audited
process and trust boundary, not the unknowable absence of every off-ledger human
influence.

## Simulation Contract

The simulator should emit a promotion ledger in addition to plots:

```bash
python3 -m oph_fpe.cli cmb-promotion-ledger \
  --run-dir runs/<run_id> \
  --out runs/<run_id>/cmb_promotion_ledger
```

The ledger reports:

- current claim tier;
- conditional imported-FLRW scale-bridge status;
- OPH-native `CosmoGeomRead_r` status;
- source-only finite artifact status;
- geometric screen scalar, radial-lift, freezeout, transfer, likelihood, and
  prediction receipts;
- forbidden input classes and semantic type boundaries;
- blockers and the fail-closed outcome.

The measurement pack carries `cmb_promotion_ledger_report.json`,
`cmb_promotion_ledger_report.md`, and summary flags in `claims.json`, so OPH
Sage and public reviewers ingest the claim tier rather than inferring it from
images.

For the current visual demonstration, the expected outcome is
`VISUAL_DIAGNOSTIC` or `SPECTRUM_DIAGNOSTIC`. Missing physical parents are an
expected incomplete state, not a numerical falsification of OPH. A terminal
prediction Boolean supplied by an output file is only a producer assertion; the
promotion ledger must recompute the terminal receipt from the required parent
gates and matching frozen hashes.

## Outcome

The claim-governance problem closes at contract level; the physical CMB problem
does not. The public simulator result is useful because it visualizes how an
observer-boundary record could enter a finite-source CMB branch, not because
visual resemblance supports that branch empirically. The OPH implementation
requires every CMB-facing output to retain its source, geometry, scale,
freezeout, transfer, no-data-use, and likelihood status. Until those receipts
pass in one hash-linked bundle and the frozen prediction is compared with
standard baselines, the correct label remains visual, spectrum-diagnostic, or
externally conditioned—not an OPH-native physical CMB prediction.
