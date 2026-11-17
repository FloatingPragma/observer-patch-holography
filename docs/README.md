# OPH Reference Documents, Registers And Research Notes

This directory contains the canonical reference documents, the generated
scientific registers, repository policies, and the research notes that
accompany the executable lanes. Research planning and task state live in this
repository's [GitHub issues](https://github.com/FloatingPragma/observer-patch-holography/issues).
Scientific results live in the papers, Lean library, executable code, and
evidence artifacts. The registers index those results, and each research note
states a result that its code package verifies.

If you are new to OPH, the strongest starting material lives outside this
directory. The [technical paper](https://philpapers.org/rec/MUEFOC)
states the primary technical account. The [OPH textbooks](https://learn.floatingpragma.ai/) work through
the basic derivations with the math taught along the way. The
[interactive simulation](https://simulation.floatingpragma.io/) renders real
run data from declared candidate models; it does not demonstrate that a
physical universe has been reconstructed.

## Layout

| Path | Contents |
| --- | --- |
| This folder | The basis, the paper set, canonical references and the common-objections guide |
| [`registers/`](registers/) | Generated V3 registers and ledgers |
| [`research/`](research/) | Research notes, each verified by a code lane |
| [`policies/`](policies/) | Data, hardware-evidence and writing policies |
| [`COMPLETION_PROGRAM.md`](COMPLETION_PROGRAM.md) | Common-source finite checks and physical acceptance conditions |

## Rules for this directory

- The root holds only the selected high-level programs, references and
  ledgers listed on this page. Add a root document only when it is a new
  program or ledger of that kind.
- No file under `docs/` is a daily note, progress report, audit result,
  integration record, PR handoff or status summary. File names and headings
  carry no dates. Integration and review history belongs in the GitHub pull
  request or issue.
- Research results live in the papers, Lean library and code. A note in
  `research/` is allowed only when it states a current result that its code
  package verifies, without history, dated findings or progress narrative.
  Prefer extending the paper or the package README, because separate notes
  go stale.
- Release papers do not cite files in `docs/`. They cite papers, Lean
  modules, code packages, machine registers and public preprints. The one
  exception is `research/PUBLIC_DATA_PREDICTION_AUDIT.md`, which collects
  retrospective public-data comparisons that no single code package
  reproduces.
- Code, receipts, reviews and registers bind code, Lean and receipts as
  evidence, never a note in `docs/`.
- Implementation-specific mathematics belongs with its owning code. Preserve
  valid links and source commitments when consolidating documents.

## Three Reading Routes

- **First encounter:** the [technical paper](https://philpapers.org/rec/MUEFOC),
  [textbooks](https://learn.floatingpragma.ai/), and the
  [simulation](https://simulation.floatingpragma.io/) above, then the
  repository [README](../README.md), the [basis](BASIS.md) and the
  [paper set](PAPER_SET.md).
- **Technical verification:** use the [claim registry](../claims/claim_registry.yaml),
  the [observation ledger](registers/OBSERVATION_LEDGER_V3.md), the [premise register](registers/PREMISE_REGISTER_V3.md),
  the [falsification program](OPH_FALSIFICATION_PROGRAM.md), and the [paper index](../paper/).
- **Build and test:** begin with the repository [reproduction
  guide](../REPRODUCE.md), the [runnable observer model](../model/README.md),
  the [bounded completion campaign](../model/completion/README.md),
  [executable evidence](../code/), and [Lean
  formalization](../Lean/). Follow the [numerical audit
  procedure](policies/NUMERICAL_AUDIT.md) when changing finite numerical
  evidence or addressing review findings.

## Basis And Paper Set

- [The OPH Basis](BASIS.md) distinguishes the closure hypothesis,
  internal-realization policy, agreement, adopted state-selection and
  response-endogeneity principles, and executable specification, and explains
  how the theorem-library premises A1, A2 and A3 relate to them.
- [The Paper Set](PAPER_SET.md) records the model milestone, the role of each
  paper and the archived papers.
- [Executable Hypothesis](OPH_EXECUTABLE_HYPOTHESIS.md) states the common
  clock, common-energy, interacting-field and Einstein-premise target.
  [Completion Program](COMPLETION_PROGRAM.md) assigns its construction gates,
  retained results, dependencies and next work.
- [Current field, particle and entropy results](../model/completion/CONTINUUM_FIELD_AND_ENTROPY.md)
  joins the component proofs and their remaining premises. The
  [mass and energy table](../model/completion/PARTICLE_MASS_TABLE.md) distinguishes
  finite source energies, supplied mass coefficients and physical predictions.

## Canonical References

- [Interpretation and Working Terms](INTERPRETATION_AND_TERMS.md) distinguishes
  complete-history closure, intrinsic computation, emergent-physics proposals,
  mathematical fixed points and the different uses of time.
- [Premise Reference](AXIOM_REFERENCE.md) gives the formal statements of the
  theorem-library premises A1, A2 and A3. It is retained byte-identical for
  pinned receipts; its "three core axioms" wording is superseded by the
  [basis](BASIS.md).
- [Source and Parameter Contract](../model/SOURCE_SELECTION.md) separates
  the executable candidate inputs, derived consequences and physical identification.
- [Frozen-Prediction Ladder](FROZEN_PREDICTION_LADDER.md) records tests whose
  conditions are fixed before the comparison data are examined.
- [Postdiction Ledger](POSTDICTION_LEDGER.md) records comparisons with measured
  values and their input ancestry.
- [Empirical Status](EMPIRICAL_STATUS.md) inventories laboratory comparisons,
  reference-value and hardware tests, and what they can decide about the
  executable repair rules.
- [OPH Falsification Program](OPH_FALSIFICATION_PROGRAM.md) lists only mature
  mathematical and realized-branch falsifiers.
- [Programme-wide obstruction search](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
  seeks inconsistency, joint-realization impossibility and limits on explanatory
  unification, with explicit exclusion and stopping criteria.

The [physics issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
assign bounded work on the common repair model, its tests and existing-evidence
audits. Historical closure issues retain their scientific boundaries;
administrative closure does not discharge a physical premise.

## Executable model

[CCF-LT-2](../model/SPEC.md) is the official finite matter, field and response
construction. Its [implementation](../model/ccf_lt2/README.md) regenerates the
native operator, complete finite diagnostics and local energy response.
The [construction contract](../model/CONSTRUCTION.md) specifies the joint law,
preparations and instruments. The [dependency map](../model/particle_joint_reconstruction/DEPENDENCIES.md)
connects mathematical proofs to the actual source and distinguishes physical
attachment obligations.

[Source inputs](../model/SOURCE_SELECTION.md), [family and controls](../model/CANDIDATES.md),
[observer records](../model/OBSERVER_RELATIVE_RECOVERY.md) and
[physical scope](../model/PHYSICAL_AUDIT.md) state the model's contracts.
The complete operator replay is `python -m model`; finite matrix experiments
have their stated carriers, and full-native column action is not arbitrary
many-body state evolution.

Other model and code packages contain conditional mathematical components.
Their equations and proof domains remain explicit. Source-dependent spectra,
geometry, clocks and numerical records are not predictions of the official
law without a checked connection. [Established Mathematics for Observer
Agreement](AGREEMENT_RESEARCH_FOUNDATIONS.md) describes the general theorem
reuse requirements.

## Registers

Each register is generated from a JSON source by a tool in [`tools/`](../tools/).
The first line of every page names its source and generator; edit the source,
then regenerate. The [five-ledger map](../claims/README.md#the-five-ledgers)
groups the sources into the claims, premises, predictions-and-tests, evidence
and papers ledgers.

- [Premise Register](registers/PREMISE_REGISTER_V3.md) names each V3 input and
  its consumers and classifies each evidence path by its scientific role.
- [Observation Ledger](registers/OBSERVATION_LEDGER_V3.md) records the physical
  targets, adequacy rungs, premise ancestry, frozen-target links, and evidence paths.
- [Emergent-Instrument Register](registers/OBSERVATION_LEDGER_V3.md#emergent-instrument-register) lists the
  simulation-instrument designs and frozen instruments of the adequacy program.
- [Constants Ancestry](registers/OBSERVATION_LEDGER_V3.md#constants-ancestry) declares the complete
  input ancestry of each constants family against the premise register.
- [Selection Ledger](registers/SELECTION_LEDGER.md) records the selections that
  consensus forces, with the exposed premises and proofs behind each row.
- [Gravity Premise Ladder](registers/PREMISE_REGISTER_V3.md#gravity-rungs) types each
  gravity interface with a premise packet, a theorem or countermodel, and a
  deletion test.
- [Standard Model Correspondence](registers/PREMISE_REGISTER_V3.md#sm-coverage)
  classifies each Standard Model term and sector by what OPH supplies.

## Research Notes

Each note states a maintained construction or theorem tied to its verifying
code package. The eight clock and source-realization notes are retained
conditional theorems of the earlier source lane; the live construction is
[`model/`](../model/README.md), which cites them only where stated. Source-pinned notes require matching receipt updates when their
content changes. Consolidate overlapping notes and keep implementation
contracts beside the code; `extra/` holds TeX papers only.

**Clocks and source realization.** These notes build on one another in this
order:

1. [Massive flights, resolved reads and observable clocks](research/MASSIVE_OPERATIONAL_CLOCKS.md),
   verified by [`code/m1_operational_clocks`](../code/m1_operational_clocks/).
2. [Coherent source codes, entropy selection and charged clocks](research/COHERENT_SOURCE_CLOCKS.md),
   verified by [`code/m1_source_realization`](../code/m1_source_realization/).
3. [Fermionic clocks and reads from the proper-code source](research/FERMIONIC_SOURCE_CLOCKS.md),
   verified by [`code/m1_fermionic_source`](../code/m1_fermionic_source/).
4. [Fixed-rate noise and protected source reads](research/FIXED_RATE_SOURCE_READS.md),
   verified by [`code/m1_fixed_noise`](../code/m1_fixed_noise/).
5. [Source histories with noisy records and control](research/NOISY_SOURCE_RECORDS.md),
   verified by [`code/m1_noisy_records`](../code/m1_noisy_records/).
6. [Fixed-strength quantum computation and live public records](research/FIXED_STRENGTH_PUBLIC_RECORDS.md),
   verified by [`code/m1_expander_archive`](../code/m1_expander_archive/).
7. [Full-interface leakage, dissipation and source refinement](research/LEAKAGE_SOURCE_REFINEMENT.md),
   verified by [`code/m1_leakage`](../code/m1_leakage/).
8. [A filled source vacuum, positive excitation energy and readable clocks](research/SOURCE_SEA_ENERGY.md),
   verified by [`code/m1_sea_energy`](../code/m1_sea_energy/).

**Gravity and the dark sector.**

- [Finite geometry: support, refinement and oriented caps](../code/geometry/GEOMETRY_RECEIPT_AUDIT.md)
  audits the geometry readout, preserves higher-dimensional joint support,
  certifies midpoint subdivisions, and reconstructs oriented caps from
  explicit conformal data and a side witness.
- [What a calibrated source clock does and does not fix about light bending](research/PAIRED_SOURCE_GRAVITY.md)
  proves that the calibrated clock-field response does not select the
  light-deflection coefficient, verified by [`code/paired_gravity`](../code/paired_gravity/).
- [Clock agreement does not select a spatial response in the coherent source](research/SOURCE_GRAVITY_ADMISSIBILITY.md)
  carries that result to programs of the constructed coherent source, verified by
  [`code/source_gravity_admissibility`](../code/source_gravity_admissibility/).
- [What the anomalous density must supply to predict lensing](research/DARK_SOURCE_DYNAMICS_LENSING.md)
  shows that the released density interface does not fix dynamics and lensing
  together, verified by [`code/dark_source_lensing`](../code/dark_source_lensing/).

**Independent comparison.**

- [OPH predictions and public-data comparisons](research/PUBLIC_DATA_PREDICTION_AUDIT.md)
  places the fixed predictions against public data: the primordial tilt at
  the lower 0.612 percent tail of ACT DR6 P-ACT-LB, the coupling, tau mass,
  dispersion, Bell and optical-clock comparisons.
- [OPH particle mass response and hierarchy audit](research/OPH_PARTICLE_MASS_RESPONSE.md)
  records that the tested mass prescriptions fail the charged-lepton
  hierarchy and bounds the whole affine response class. Its text is pinned
  by retained receipts. The [finite campaign](../model/completion/README.md)
  carries the subsequent source-charge and clock-access decisions.
- [Observer spacetime and predictions without a laboratory](research/NATIVE_CLOCK_EXPERIMENT_PROTOCOL.md)
  states the internal clock laws and the acceptance contract for a native
  clock experiment.
- [Can the fixed edge dispersion support the proposed photon test?](research/DISPERSION_FEASIBILITY.md)
  bounds the exact pair-production thresholds and distinguishes vacuum emission
  in two specified lepton controls, verified by
  [`code/dispersion_feasibility`](../code/dispersion_feasibility/).
- [Independent alpha/tau reproduction and the first bounded natural comparison](research/INDEPENDENT_POSTDICTION_COMPARISON.md)
  reproduces both certified alpha roots and the conditional tau window, verified by
  [`code/independent_postdictions`](../code/independent_postdictions/).

**Defect register and particle classes.** The
[defect-class certificate](../model/defect_particle_classes/README.md)
identifies each protected defect class with the centre class `6Y mod 6` of a
Standard Model multiplet, and the [integer lift](../model/defect_integer_lift/README.md)
checks the register on the stored runs. The [lepton register](../model/lepton_register_round/README.md)
gives the charged-lepton ratios `206.770316` and `3477.4728` with no measured
mass as input. On the lattice, the declared coupling chain gives the
[electroweak crossover](../model/electroweak_crossover/README.md), the colour
links give [deconfinement](../model/colour_deconfinement/README.md) and
[hadron binding](../model/quenched_hadron_binding/README.md), and the
[seam class balance](../model/seam_class_balance/README.md) measures the sector
law behind the electroweak scale. The pinned notes below state the boundaries
of their own source constructions.

**Field and particle construction.** The source-pinned
[field reconstruction guide](research/OPH_FIELD_RECONSTRUCTION.md) records
the transient-collision test, conditional response constructions and their
physical identification requirements. Its statement that other quantum or
history sectors are untested belongs to that receipt's domain. The
[branch-algebra result](../model/completion/defect_algebra.md) additionally
excludes commuting bare-source charges in every quantum basis; localized
history-dependent sectors, controlled lifetimes and refinements require their
own tests. Both pinned guides retain their original evidence bytes; the
[completion program](COMPLETION_PROGRAM.md) owns the work ordering.

**Lean companion.**

- [The declared repair law on the Euclidean history register](research/SOURCE_LAW_HISTORY_REGISTER.md)
  is the companion note of the Lean `SourceLaw` library: inertia, the golden
  refinement family and the reference state.

**Existing-theory audit and simplification.**

- [MaxEnt inference independent of observable coordinates](../code/maxent/PROJECTION_COORDINATES.md)
  repairs unit, basis and energy-offset dependence, proves a global
  optimization-error bound, and distinguishes that error from closure defect.
- [Finite-state entropy and sector labels](../code/quantum_information/README.md)
  records the entropy and collar audit, gives the direct-sum form of the finite
  A3 objective, and documents the shared implementation and corrected
  Markov-state, alignment and MaxEnt checks.
- [When finite entropy completion is a repair channel](../code/quantum_information/ALGEBRAS_AND_REPAIR.md)
  proves the equivalence between affine full-marginal completion, modular
  compatibility and reference-preserving conditional expectation. It corrects
  regional standardness checks and treats noncommuting primitive repairs.
- [Whitney numerical methods](../code/electromagnetism/WHITNEY_NUMERICS.md)
  gives the positive kinetic reduction, exact neutral-packet projection,
  numerical contracts and independent controls beside their implementation.
- [Finite null tomography: reconstruction, consistency and resolution](../code/geometry/NULL_TOMOGRAPHY.md)
  repairs incomplete and non-null tensor reconstructions, reuses the proved
  nine-direction inverse, and connects the dependent-family test to angular
  harmonics, an exact cubic obstruction and sharp noise amplification.

## Instrument Specifications

- [INS-03 source-bound phase-sensitive readout](../code/instrument_specs/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md)
- [OL-A1 signature replication](../code/instrument_specs/OL_A1_SIGNATURE_REPLICATION_SPEC.md)
  (status `FAILED`; see the [physical audit](../model/PHYSICAL_AUDIT.md))
- [OL-A1 carrier count by absolute support size](../code/instrument_specs/OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md)

## Evidence, Data And Writing Policies

- [Hadron Data Policy](policies/HADRON.md) defines provenance and promotion rules for the hadronic pipeline.
- [Hardware Evidence Bundle H](policies/HARDWARE_EVIDENCE_BUNDLE_H.md) defines the
  evidence and attestation contract for physical hardware claims.
- [Writing Style Guide](policies/STYLE_GUIDE.md) defines the project’s reader-facing prose conventions.
- The [IBM Quantum Cloud code and data archive](../code/ibm_quantum_cloud/README.md)
  contains the reproducible engineering benchmark, receipts, and interpretation
  boundary beside their executable producers.
- Simulator outputs have no theorem or empirical-promotion authority in the
  paper stack. A simulator count, replay, fit, or diagnostic closes no
  analytic or physical receipt; theorem status needs analytic proof packets,
  and empirical promotion needs public physical evidence bundles.

## Reader Support

- [Common Objections](COMMON_OBJECTIONS.md) answers recurring technical and conceptual criticisms.
- [Contributed background notes](../contributions/README.md) record standard mathematics that OPH
  work can reference; they make no OPH claim.

The [screened recovery service](../model/screened_recovery/README.md),
[same-code transport](../model/screened_transport/README.md),
[four shared clock events](../model/clock_defect_energy/shared_event_clock/README.md)
and [conserved-sector action theorem](../model/conserved_sector_action/README.md)
record the finite frontier results. The [SI dictionary](../model/clock_defect_energy/SI_UNITS.md)
provides conditional unit conversion without a particle-mass prediction.
