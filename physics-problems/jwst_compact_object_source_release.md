# JWST Compact Objects As Source-Release Record Surfaces

## Motivating Result

This note entered the queue after JWST made the question unavoidable: early red
galaxy candidates appeared massive only 500 to 700 million years after the Big
Bang
([Labbe et al., Nature 616, 266-269, 2023](https://www.nature.com/articles/s41586-023-05786-2)).
The later "little red dot" literature sharpened the ambiguity: compact red
objects with broad lines can be read as starbursts, dust, dense gas, selection,
or accreting black holes
([Nature, 2025](https://www.nature.com/articles/s41586-025-09900-4)). The OPH
question is how to audit compactness, color, luminosity, AGN contribution,
source release, and mass-age interpretation without treating a catalog row as a
physical closure.

**Status:** solved as a source-release audit theorem package and simulator
receipt contract. The JWST lane remains source-open and simulation-buildable,
not a physical mass-age or OPH-confirmation claim.

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

The redundancy groupoid removes gauge representatives, mesh labels, port
labels, packet labels, basis choices, hidden carrier coordinates, worker IDs,
queue order, retry counters, repair schedules, and inert ancillas:

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

Useful compactness/readout indices are:

```math
\kappa_R = R_{rec}/R_{halo},
\qquad
\kappa_{PSF}=A_{rec}/(A_{PSF}/\mu_{lens}),
\qquad
\Sigma_{rec}=\frac{\sum_e a_ep_{e,rec}}{A_{rec}}.
```

Low synchronization residual `Phi_sync` means observer-facing records align
coherently across the compact surface. It does not imply that the system is
old, massive, or fully assembled.

## Record Density

For an observer bin `B`, define compact-object record load:

```math
\mathsf L_{obj,r}(q;B)
=
\sum_C V^{phys}_{C,r}(q)
\sum_{e\in E_{rec,r}(C)}
a_ep_{e,rec}(q)
\mathbf 1_B[\mathsf{Read}_{obj,r}(e,q)].
```

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
   The pushforward of a law by a normal-form map is idempotent, so every law
   supported on the normal-form image is fixed. A source law must be
   supplied separately.

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

6. **JWST maturity is non-identifiable when degeneracy support is large.**  
   If two quotient states have indistinguishable JWST catalog readouts but
   separated physical maturity tuples, no estimator using only the catalog
   record can uniformly recover physical maturity below half that separation,
   up to observational stability.

7. **Mass/age tension promotion requires full degeneracy closure.**  
   A compact object may be labeled `PHYSICAL_MASS_AGE_TENSION` only after
   redshift, SED, dynamical, lensing, AGN, dust, nebular, population,
   morphology, PSF, and selection receipts close and the remaining allowed
   physical support lies outside the conventional baseline by a predeclared
   threshold.

8. **Source release and synchronization can change apparent maturity.**  
   OPH source release may shift compactness, surface brightness, line
   coherence, or record density. It does not by itself imply old age, high
   stellar mass, large black-hole seed mass, or completed assembly.

9. **Forward counts are pushforwards.**  
   JWST-facing expected counts are the pushforward histogram of
   `mu^{obj}_r` through the frozen JWST forward operator.

10. **Refinement compatibility bounds finite evidence.**  
    If source-law and load defects are bounded under coarse maps, then object
    abundance changes across regulators are bounded by the total-variation and
    load-defect terms. A refinement failure blocks source-only promotion.

11. **Simulator soundness is bounded-observable soundness.**  
    For bounded catalog statistics, simulation error is bounded by sampling,
    refinement, and forward-operator defects. Continuum or physical maturity
    claims need the stronger parent and likelihood receipts.

## Claim-Tier Ladder

```text
J0_DIAGNOSTIC_PROXY
J1_CATALOG_RECORD
J2_SPECTROSCOPIC_OR_PHOTOMETRIC_OBJECT
J3_DEGENERACY_AUDITED_OBJECT
J4_CONDITIONAL_PHYSICAL_OBJECT
J5_SOURCE_RELEASE_CANDIDATE
J6_SOURCE_ONLY_OBJECT_ABUNDANCE
J7_FORWARD_MOCK_PHYSICAL_SPECTRUM
J8_LIKELIHOOD_EVALUATED_PHYSICAL_PREDICTION
```

Promotion requires:

- J1: catalog ingestion and provenance.
- J2: object ID, aperture, redshift posterior, photometry/spectrum, and
  morphology records.
- J3: redshift, dust, AGN, nebular, stellar-population, lensing, morphology,
  PSF, and selection receipts.
- J4: object parent, packet mass-shell, packet stress, total-stress closure,
  radiative transfer, lensing source-plane, and chemical/SFH receipts.
- J5: OPH release state and source-release residual not tuned to the object.
- J6: object abundance source receipt, no-target leakage receipt, and load
  refinement compatibility.
- J7: frozen JWST forward operator.
- J8: frozen catalog likelihood and predeclared falsification thresholds.

## Catalog-Level Checklist

Every catalog row should carry:

- object ID, survey, coordinates, aperture definition, and provenance hashes;
- photometry, spectrum, morphology, PSF, and selection hashes;
- redshift posterior and line-identification posterior hashes;
- SED, AGN, dust, nebular, stellar-population, and lensing model-family hashes;
- dynamical-mass receipt status;
- compact-record-surface receipt status;
- degeneracy-audit receipt status;
- claim tier `J0` through `J8`;
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

The simulator must implement this chain:

```text
quotient -> source law -> release state -> record density -> object parent
-> JWST forward operator -> degeneracy audit -> frozen likelihood
```

It must fail promotion when:

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
the strongest allowed claim tier for every object and population bin.

The active simulator surface is `oph-physics-sim/oph_fpe/jwst`, with commands
such as:

```bash
python3 -m oph_fpe.cli jwst-object-source-artifact --out runs/jwst/source
python3 -m oph_fpe.cli jwst-compact-object-simulation-plan \
  --run-dir runs/jwst \
  --out runs/jwst/plan
```

The paper-stack code mirror is
`reverse-engineering-reality/code/particles/jwst/build_compact_object_source_release_receipts.py`.
It emits the same fail-closed claim ladder into the OPH content repository:

```bash
cd reverse-engineering-reality/code/particles
python3 jwst/build_compact_object_source_release_receipts.py \
  --output runs/jwst/compact_object_source_release
```

The remaining hard paper-side object is the actual numerical compact-object
source action `S^{obj}_r(q)`. Until OPH derives it upstream or declares it as a
conditional branch with no-data receipts, the JWST lane remains source-open and
simulation-buildable rather than physically closed.
