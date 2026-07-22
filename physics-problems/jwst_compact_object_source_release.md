# JWST Compact Objects As Source-Release Record Surfaces

## Motivating Result

This note entered the queue after JWST made the question unavoidable: early red
galaxy **candidates**, interpreted through stellar-population models, appeared
massive only 500 to 700 million years after the Big Bang
([Labbe et al., Nature 616, 266-269, 2023](https://www.nature.com/articles/s41586-023-05786-2)).
The later "little red dot" literature sharpened the ambiguity: compact red
objects with broad lines can be read as starbursts, dust, dense gas, selection,
or accreting black holes
([Rusakov et al., Nature 649, 574--579, 2026](https://www.nature.com/articles/s41586-025-09900-4)).
That analysis found electron-scattering-broadened lines in its high-quality
sample and inferred black-hole masses about two orders of magnitude below prior
line-width estimates, illustrating why the forward model matters. The OPH
question is how to audit compactness, color, luminosity, AGN contribution,
source release, and mass-age interpretation without treating a catalog row as a
physical closure.

**Status:** conditional audit specification and receipt schema. It does not
derive a compact-object source law, execute a physical JWST forward model, or
establish an OPH mass-age prediction. The current lane is source-open and may
support visual or synthetic demonstrations when their assumptions are labeled.

Date: 2026-07-08

## Origin

This note records the OPH compact-object audit prompted by the motivating JWST
observations. It asks how OPH should treat JWST high-redshift compact objects,
early massive galaxy candidates, little red dots, and apparently mature
black-hole candidates without claiming that JWST has confirmed OPH.

JWST sees bright early objects. The audit question is whether luminosity, color,
compactness, apparent maturity, AGN contribution, and selection effects have
been collapsed into one physical claim too quickly.

The OPH answer is:

```text
JWST compact objects are an audit surface for source release, compact record
surfaces, finite object parents, forward modeling, and degeneracy receipts.
They are not confirmation of OPH by themselves.
```

## What Is Standard, What Is Open, And What OPH Adds

Standard astronomy supplies photometric and spectroscopic inference,
SED and line-profile modeling, lensing, radiative transfer, population
selection, PSF convolution, and catalog likelihoods. The unresolved OPH
physics is a source-derived object population and emission history that makes a
distinct prediction after those conventional effects are included.

OPH's useful addition is the observer-like self-reading object patch: bounded
local state, release boundary, readback records, repair/refinement lineage, and
a public evidence bundle linking source assumptions to detector-level rows.
This can expose category errors such as luminosity-to-mass or compactness-to-age
promotion. It is not a unique explanation of early compact sources.

## Claim Boundary

The lane must keep these non-identifications explicit:

```text
luminosity is not stellar mass
red color is not old age
compactness is not completed assembly
AGN contribution is not stellar mass
selection is not physical abundance
source release is not galaxy formation
record density is not stellar-mass density
```

OPH may predict observer-facing record-density or compactness structure only
after it declares a source law, release state, object parent, forward operator,
degeneracy audit, and frozen likelihood. A normal form, attractive image, or
compact surface does not select a probability law.

## End-To-End Object Chain

The theorem-gate chain is:

```text
compact-object quotient
-> compact-object source law
-> source-release state
-> compact record surface and record load
-> finite covariant object parent
-> JWST forward operator
-> degeneracy audit
-> frozen catalog likelihood
```

In symbols, the required paper-side object family is:

```math
Q^{obj}_r
\rightarrow
\mu^{obj}_r
\rightarrow
\mathfrak R^{rel}_{obj,r}
\rightarrow
\mathsf L_{obj,r}
\rightarrow
\mathcal P^{obj}_r
\rightarrow
\mathsf F^{JWST}_r
\rightarrow
\mathsf{Audit}_{deg}
\rightarrow
\mathcal L_{frozen}.
```

The simulator can instantiate and test this chain. It cannot replace the
source law, quotient, finite parent, or promotion theorem.

## Compact-Object Quotient

At regulator `r`, the compact-object presentation space is:

```math
\Sigma^{obj}_r.
```

An element contains finite causal cells, metric/tetrad/parallel transport data,
baryon/gas/star/dust/black-hole/accretion/radiation sectors, object-facing
records, environment data, repair/update maps, checkpoint lineage, emission,
lensing, and radiative-transfer maps.

The redundancy groupoid removes gauge representatives and demonstrably inert
presentation metadata such as mesh names, worker IDs, queue order, retry
counters, and random-number bookkeeping. Physical ports, release surfaces,
packet species, repair histories, and any label that changes emission, lensing,
selection, or dynamics remain quotient-visible:

```math
Q^{obj}_r=\Sigma^{obj}_r/\Gamma^{obj}_r.
```

Only quotient-visible functions may enter source laws, abundance selectors,
physical parents, forward mocks, or claim-tier receipts.

## Source Law

A source-only compact-object claim requires:

```math
\mu^{obj}_r(q)
=
Z^{-1}_{obj,r}m^{obj}_r(q)
\exp[-S^{obj}_r(q)].
```

For a finite classical branch, $m_r^{obj}\ge0$, the action is dimensionless,
and $0<Z_{obj,r}=\sum_qm_r^{obj}(q)e^{-S_r^{obj}(q)}<\infty$. A
continuous branch must declare its base measure, domain, and integrability; a
quantum branch requires a positive trace-one density operator. Writing a Gibbs
symbol does not establish any of these properties.

Every constraint in the action must declare its name, units, support, target,
source derivation, zero-mode convention, sector convention, refinement map, and
no-data-use proof.

The hard rule is:

```text
NORMAL_FORM_MAP != SOURCE_LAW
```

A normal-form map can classify object records. It cannot select object
abundances, formation rates, or physical maturity distributions.

## Source Release

A compact-object source release is the finite event by which a quotient-visible
object packet becomes available to the observer-facing light-cone/readout
pipeline.

It records the release surface, normal/tetrad, release scale factor, release
clock, comoving volume, object projector, branch data, generation hash,
no-target-leak ledger, coordinate chart, and light-cone map.

Source release is a record-availability event. It is not, by itself, physical
galaxy formation.

## Compact Record Surface

A compact record surface records:

```math
\mathcal C^{rec}_r
=
(C_r,E_{rec,r},p_{e,rec},a_e,x_e,A_{rec},R_{rec},
\mathcal E_{env},\Phi_{sync},\mathsf{Emit}_r).
```

Useful compactness/readout indices must compare quantities in one plane. Under
a locally constant scalar magnification, define
$A^{src}=A^{img}/\mu_{lens}$ for both the object and PSF. Then

```math
\kappa_R = R_{rec}^{src}/R_{halo}^{src},
\qquad
\kappa_{PSF}
=
\frac{A_{rec}^{img}}{A_{PSF}^{img}}
=
\frac{A_{rec}^{src}}{A_{PSF}^{src}},
\qquad
\Sigma_{rec}=\frac{\sum_e a_ep_{e,rec}}{A_{rec}^{src}}.
```

For anisotropic lensing, the full Jacobian and source-plane PSF are required;
a scalar $\mu_{lens}$ is insufficient. The record activation
$p_{e,rec}$ is dimensionless, while $a_e$ must declare its record-weight
units, so $\Sigma_{rec}$ has record weight per physical area.

Low synchronization residual `Phi_sync` means observer-facing records align
coherently across the compact surface. It does not imply that the system is
old, massive, or fully assembled.

## Record Density

For an observer bin `B`, define compact-object record load:

```math
\mathsf L_{obj,r}(q;B)
=
\sum_C\sum_{e\in E_{rec,r}(C)}
a_ep_{e,rec}(q)
\mathbf 1_B[\mathsf{Read}_{obj,r}(e,q)].
```

With $a_e$ in record-weight units, $\mathsf L$ is a record weight and
$\rho_{rec}$ below is record weight per comoving volume. A volume-weighted
field is a different construction: it must introduce a proper-volume density
$`w^{vol}_{e,C}`$, use $`V_C^{phys}w^{vol}_{e,C}`$, and apply the explicit
proper-to-comoving factor. The surface weight $a_e$ must not be multiplied
by cell volume without that new definition.

The source-side record density is:

```math
\rho_{rec,r}(B)
=
V^{-1}_{com,r}
\mathbb E_{\mu^{obj}_r}[\mathsf L_{obj,r}(q;B)].
```

This is an object-lane abundance analogue. It is not automatically stellar-mass
density, halo-mass density, black-hole mass density, or galaxy number density.

## Theorem Gates

1. **Normal forms do not select compact-object probability laws.**  
   Let $c:Q_r^{obj}\to Q_r^{obj}$ satisfy $c^2=c$, with normal-form
   image $N=\operatorname{im}c$. Then $`c_\#(c_\#\mu)=c_\#\mu`$, and
   every law supported on $N$ is fixed. Idempotence therefore leaves many
   laws possible; a source law must be supplied separately.

2. **Record density is quotient-invariant.**  
   If load, volume, record activation, readout, and bin membership are defined
   on `Q^{obj}_r`, then hidden representatives cannot change
   `rho_rec,r(B)`.

3. **Compactness is readout, not assembly.**  
   Small effective radius, high record density, or low PSF compactness does not
   imply old stellar age, high stellar mass, high dynamical mass, or completed
   assembly unless the finite object parent and degeneracy receipts close.

4. **Finite total-stress closure requires a parent.**  
   A physical object parent must carry packet sectors, mass shells, momenta,
   weights, exchange channels, stress readouts, radiative transfer, lensing,
   chemistry, and star-formation history. Scalar rows cannot be promoted to
   stress, mass, or maturity.

5. **Source-only abundance needs a no-target source DAG.**  
   The expected object load is source-only only when the quotient, source law,
   release state, object parent, load map, and no-target ledger are frozen
   without paths from JWST catalog counts, anomaly labels, posterior summaries,
   likelihood residuals, or hand-selected interesting objects.

6. **Exact readout degeneracy gives a minimax lower bound.**
   Let maturity values lie in a metric space $(M,d_M)$. If two admissible
   states have the same frozen JWST readout but
   $d_M(m_1,m_2)=\Delta$, then every estimator from that readout has
   worst-case error at least $\Delta/2$, by the triangle inequality. A
   noisy or approximate version needs an explicit readout metric, noise ball,
   and stability modulus; "large degeneracy support" alone is not a theorem.

7. **Mass/age tension promotion requires full degeneracy closure.**  
   A compact object may be labeled `PHYSICAL_MASS_AGE_TENSION` only after
   redshift, SED, dynamical, lensing, AGN, dust, nebular, population,
   morphology, PSF, and selection receipts close and the remaining allowed
   physical support lies outside the conventional baseline by a predeclared
   threshold.

8. **Source-release effects are a provisional hypothesis.**
   A declared OPH response model may shift compactness, surface brightness,
   line coherence, or record density. No such shift follows from the abstract
   release definition alone. Even when assumed in a visual model, it does not
   imply old age, high stellar mass, large black-hole seed mass, or completed
   assembly.

9. **Forward counts are pushforwards.**  
   JWST-facing expected counts are the pushforward histogram of
   `mu^{obj}_r` through the frozen JWST forward operator.

10. **Refinement compatibility bounds bounded loads.**
    For a coarse map $c_{sr}:Q_s\to Q_r$, assume
    $`\|L_r\|_\infty\le M`$,
    $`\|\mu_r-c_{sr\#}\mu_s\|_{TV}\le\epsilon_\mu`$, and
    $`\sup_{q\in Q_s}|L_s(q)-L_r(c_{sr}q)|\le\epsilon_L`$. Then
    $`|\mathbb E_{\mu_s}L_s-\mathbb E_{\mu_r}L_r|\le\epsilon_L+M\epsilon_\mu`$,
    using the convention
    $`|E_\mu f-E_\nu f|\le\|f\|_\infty\|\mu-\nu\|_{TV}`$.
    The convention and finite bound must be recorded. A refinement failure
    blocks source-only promotion.

11. **Simulator soundness is conditional bounded-observable soundness.**
    For bounded catalog statistics, one may add declared sampling, refinement,
    and forward-operator error bounds. Those terms require a concentration
    inequality and operator norm; their mere names are not a bound. Continuum
    or physical-maturity claims need the stronger parent and likelihood
    receipts.

## Two-Axis Claim Status

```text
OBSERVATION AXIS
O0_DIAGNOSTIC_PROXY
O1_CATALOG_RECORD
O2_SPECTROSCOPIC_OR_PHOTOMETRIC_OBJECT
O3_DEGENERACY_AUDITED_OBJECT

OPH MODEL AXIS
M0_NO_PHYSICAL_SOURCE_MODEL
M1_CONDITIONAL_PHYSICAL_PARENT
M2_FROZEN_SOURCE_RELEASE_CANDIDATE
M3_SOURCE_ONLY_OBJECT_ABUNDANCE
M4_FORWARD_MOCK_PHYSICAL_SPECTRUM
M5_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION
```

Every result receives a pair $(O_i,M_j)$. Better observations can raise the
first coordinate without providing an OPH source model; a sophisticated OPH
mock can raise the second only if its parent receipts pass. This prevents a
catalog-quality improvement from being mistaken for OPH confirmation.

Observation promotion requires provenance for O1; object ID, aperture,
redshift posterior, photometry/spectrum, and morphology for O2; and redshift,
dust, AGN, nebular, stellar-population, lensing, PSF, morphology, and selection
alternatives for O3.

Model promotion requires the physical parent and stress/transfer receipts for
M1; a source-release state and response model frozen without object tuning for
M2; source
abundance, no-target leakage, and refinement for M3; a frozen JWST forward
operator for M4; and a frozen likelihood plus predeclared falsification
thresholds for M5.

## Catalog-Level Checklist

Every catalog row should carry:

- object ID, survey, coordinates, aperture definition, and provenance hashes;
- photometry, spectrum, morphology, PSF, and selection hashes;
- redshift posterior and line-identification posterior hashes;
- SED, AGN, dust, nebular, stellar-population, and lensing model-family hashes;
- dynamical-mass receipt status;
- compact-record-surface receipt status;
- degeneracy-audit receipt status;
- observation tier `O0` through `O3` and OPH-model tier `M0` through `M5`;
- blocking receipts;
- printed nonclaims.

## Required Receipts

```text
OBJECT_QUOTIENT_ENSEMBLE_RECEIPT
OBJECT_RELEASE_STATE_RECEIPT
OBJECT_SOURCE_LAW_RECEIPT
OBJECT_PARENT_RECEIPT
PACKET_MASS_SHELL_RECEIPT
FINITE_PACKET_STRESS_READOUT_RECEIPT
TOTAL_STRESS_CLOSURE_RECEIPT
RADIATIVE_TRANSFER_RECEIPT
LENSING_SOURCE_PLANE_RECEIPT
CHEMICAL_SFH_RECEIPT
COMPACT_RECORD_SURFACE_RECEIPT
SYNC_RESIDUAL_RECEIPT
RECORD_DENSITY_RECEIPT
JWST_FORWARD_OPERATOR_RECEIPT
JWST_SELECTION_RECEIPT
DEGENERACY_AUDIT_RECEIPT
OBJECT_ABUNDANCE_SOURCE_RECEIPT
NO_TARGET_LEAKAGE_RECEIPT
LOAD_REFINEMENT_COMPATIBILITY_RECEIPT
FROZEN_CATALOG_LIKELIHOOD_RECEIPT
JWST_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION_RECEIPT
```

## Simulator Contract

A future physical simulator must implement this chain. A schema-only run may
display it, but must remain $(O0,M0)$ and mark every uninstantiated parent
receipt false:

```text
quotient -> source law -> release state -> record density -> object parent
-> JWST forward operator -> degeneracy audit -> frozen likelihood
```

Any implementation must fail promotion when:

- a normal-form map is used as a probability law;
- representative sampling changes orbit-size weights without declaring that
  measure;
- source artifacts read JWST catalog counts, residuals, posteriors, anomaly
  labels, or likelihood values;
- compact record surfaces are promoted to physical assembly;
- red SEDs are promoted to old age without degeneracy receipts;
- luminosity is promoted to high stellar mass without parent, dynamical, SED,
  and lensing receipts;
- broad lines are promoted to large black-hole mass without AGN/cocoon/scatter
  alternatives;
- scalar rows are promoted to stress, gas mass, dynamical mass, or maturity;
- lensing magnification is omitted from luminosity, size, or compactness
  claims;
- nonlinear shard-local reductions are used for global abundance or likelihood;
- the forward operator changes after catalog comparison;
- falsification thresholds are chosen after seeing likelihood residuals;
- source abundance is normalized to observed JWST counts;
- producer-supplied booleans are accepted without recomputation;
- exact source/config/solver/likelihood hashes are missing from a physical
  promotion receipt.

## Falsification Rule

For a frozen generation

```text
G_n = (H_source, H_release, H_parent, H_forward, H_catalog,
       H_likelihood, H_priors, A_n),
```

the lane retracts or fails under any of these conditions:

- `NO_TARGET_LEAKAGE_RECEIPT` fails. This is invalidation, not physical
  falsification.
- a run promotes physical maturity without `OBJECT_PARENT_RECEIPT` or
  `TOTAL_STRESS_CLOSURE_RECEIPT`;
- an object with open redshift, dust, AGN, lensing, nebular, morphology, or
  selection receipts is labeled `PHYSICAL_MASS_AGE_TENSION`;
- refinement violates the declared abundance-defect bound;
- frozen likelihood is worse than the baseline by the predeclared threshold;
- predeclared OPH correlations between record density, compactness, release
  variables, synchronization residuals, and environment vanish within tolerance;
- ordinary astrophysical models close the object and population ledgers after
  dust, AGN, lensing, interlopers, stellar-population effects, and selection are
  marginalized.

## Workbench

The first useful workbench is not a claim that JWST has confirmed OPH. It is a
receipt-producing pipeline that ingests catalog rows, freezes the forward model,
constructs degeneracy pairs, verifies no-target source artifacts, and reports
the strongest allowed observation/model status pair for every object and
population bin.

The current receipt-scaffold surface is `oph-physics-sim/oph_fpe/jwst`, with
commands
such as:

```bash
python3 -m oph_fpe.cli jwst-object-source-artifact --out runs/jwst/source
python3 -m oph_fpe.cli jwst-compact-object-simulation-plan \
  --run-dir runs/jwst \
  --out runs/jwst/plan
```

The paper-stack code mirror is
`reverse-engineering-reality/code/particles/jwst/build_compact_object_source_release_receipts.py`.
It emits the same fail-closed schema into the OPH content repository:

```bash
cd reverse-engineering-reality/code/particles
python3 jwst/build_compact_object_source_release_receipts.py \
  --output runs/jwst/compact_object_source_release
```

These commands do not execute the source-to-spectrum chain or validate a
physical object parent. The remaining hard paper-side object is the actual
numerical compact-object source action `S^{obj}_r(q)`. Until OPH derives it
upstream or declares it as a
conditional branch with no-data receipts, the JWST lane is source-open and
simulation-buildable. It carries no physically closed status.
