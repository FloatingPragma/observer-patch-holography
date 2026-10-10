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
states the primary technical account. The [OPH textbooks](https://learn.floatingpragma.io/) work through
the basic derivations with the math taught along the way. The
[interactive simulation](https://simulation.floatingpragma.io/) renders real
run data so you can watch the universe assemble itself.

## Layout

| Path | Contents |
| --- | --- |
| This folder | Canonical references and the common-objections guide |
| [`registers/`](registers/) | Generated V3 registers and ledgers |
| [`research/`](research/) | Research notes, each verified by a code lane |
| [`instrument_specs/`](instrument_specs/) | Instrument designs listed in the instrument register |
| [`policies/`](policies/) | Data, hardware-evidence and writing policies |

The five canonical references keep fixed paths because verifier receipts, the
Lean core and frozen custody records cite them by path. Published essays link
to the common-objections guide at its path in this folder. Records written before
2026-10-04 cite the other documents at `docs/<NAME>.md` or `extra/<NAME>.md`;
each file keeps its name in the folders above.

## Three Reading Routes

- **First encounter:** the [technical paper](https://philpapers.org/rec/MUEFOC),
  [textbooks](https://learn.floatingpragma.io/), and the
  [simulation](https://simulation.floatingpragma.io/) above, then the
  repository [README](../README.md) from the three axioms through the twist.
- **Technical verification:** use the [claim registry](../claims/claim_registry.yaml),
  the [observation ledger](registers/OBSERVATION_LEDGER_V3.md), the [premise register](registers/PREMISE_REGISTER_V3.md),
  the [falsification program](OPH_FALSIFICATION_PROGRAM.md), and the [paper index](../paper/).
- **Build and test:** begin with the repository [reproduction
  guide](../REPRODUCE.md), [executable evidence](../code/), and [Lean
  formalization](../Lean/). Follow the [numerical audit
  procedure](policies/NUMERICAL_AUDIT.md) when changing finite numerical
  evidence or addressing review findings.

## Canonical References

- [Axiom Reference](AXIOM_REFERENCE.md) states the three core axioms.
- [Canonical Repair-Law RFC](CANONICAL_REPAIR_LAW_RFC.md) specifies a proposed
  strengthening of A1 and A2.
- [Frozen-Prediction Ladder](FROZEN_PREDICTION_LADDER.md) records tests whose
  conditions are fixed before the comparison data are examined.
- [Postdiction Ledger](POSTDICTION_LEDGER.md) records comparisons with measured
  values and their input ancestry.
- [OPH Falsification Program](OPH_FALSIFICATION_PROGRAM.md) lists only mature
  mathematical and realized-branch falsifiers.

Each quantitative closure condition is tracked as a
[GitHub issue labeled `closure`](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure),
with its evaluation boundary and required completion stated on the issue.

## Model Construction Status

[Model Search](MODEL_SEARCH.md) records what is constructed, which physical
connections remain missing, and the rule to start from the simplest existing
OPH carrier and repair implementation. It includes the model campaign's
scoped results, source-reread corrections and completion criteria for the
physical-model issues. It is a research status map, not a replacement for the
claim registry or a new verified-result package.

## Registers

Each register is generated from a JSON source by a tool in [`tools/`](../tools/).
The first line of every page names its source and generator; edit the source,
then regenerate.

- [Premise Register](registers/PREMISE_REGISTER_V3.md) names each V3 input and
  its consumers and classifies each evidence path by its scientific role.
- [Observation Ledger](registers/OBSERVATION_LEDGER_V3.md) records the physical
  targets, adequacy rungs, premise ancestry, frozen-target links, and evidence paths.
- [Emergent-Instrument Register](registers/INSTRUMENT_REGISTER_V3.md) lists the
  simulation-instrument designs and frozen instruments of the adequacy program.
- [Constants Ancestry](registers/CONSTANTS_ANCESTRY_V3.md) declares the complete
  input ancestry of each constants family against the premise register.
- [Selection Ledger](registers/SELECTION_LEDGER.md) records the selections that
  consensus forces, with the exposed premises and proofs behind each row.
- [Gravity Premise Ladder](registers/GRAVITY_PREMISE_LADDER.md) types each
  gravity interface with a premise packet, a theorem or countermodel, and a
  deletion test.
- [Standard Model Correspondence](registers/SM_LAGRANGIAN_CORRESPONDENCE.md)
  classifies each Standard Model term and sector by what OPH supplies.

## Research Notes

Each note states a construction, theorem or audit in Markdown; the code
package named with it verifies the note and pins it by path and SHA-256 in its
receipt. A change to a note therefore rebuilds that receipt in the same commit.
New notes go in this folder; `extra/` holds TeX papers only, and the claim
registry check rejects Markdown notes there.

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

- [Can the fixed edge dispersion support the proposed photon test?](research/DISPERSION_FEASIBILITY.md)
  bounds the exact pair-production thresholds and distinguishes vacuum emission
  in two specified lepton controls, verified by
  [`code/dispersion_feasibility`](../code/dispersion_feasibility/).
- [Independent alpha/tau reproduction and the first bounded natural comparison](research/INDEPENDENT_POSTDICTION_COMPARISON.md)
  reproduces both certified alpha roots and the conditional tau window, verified by
  [`code/independent_postdictions`](../code/independent_postdictions/).

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
- [Finite null tomography: reconstruction, consistency and resolution](../code/geometry/NULL_TOMOGRAPHY.md)
  repairs incomplete and non-null tensor reconstructions, reuses the proved
  nine-direction inverse, and connects the dependent-family test to angular
  harmonics, an exact cubic obstruction and sharp noise amplification.

## Instrument Specifications

- [INS-03 source-bound phase-sensitive readout](instrument_specs/INS03_SOURCE_BOUND_PHASE_INSTRUMENT_DESIGN.md)
- [OL-A1 signature replication](instrument_specs/OL_A1_SIGNATURE_REPLICATION_SPEC.md)
- [OL-A1 carrier count by absolute support size](instrument_specs/OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md)

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
