# OPH Reference Documents, Registers And Research Notes

This directory contains the canonical reference documents, the observables
register, repository policies, and the public-data comparison that the papers
cite. Completion of the theory is tracked in the
[observables register](registers/OBSERVABLES.json): every observable the
theory must reproduce and every premise of its finite theorems, each with its
status. Scientific results live in the papers, Lean library, executable code,
and evidence artifacts, and the register indexes those results.

If you are new to OPH, the strongest starting material lives outside this
directory. The [technical paper](https://philpapers.org/rec/MUEFOC)
states the primary technical account. The [OPH textbooks](https://learn.floatingpragma.ai/) work through
the basic derivations with the math taught along the way. The
[interactive simulation](https://simulation.floatingpragma.io/) renders real
run data from declared candidate models.

## Layout

| Path | Contents |
| --- | --- |
| This folder | The basis, the paper set, canonical references, the common-objections guide and the OpenWave submission guidelines |
| [`registers/`](registers/) | The observables register: observables and premises with their status |
| [`research/`](research/) | The public-data comparison that the papers cite |
| [`policies/`](policies/) | Data and writing policies |

## Rules for this directory

- The root holds only the selected high-level programs, references and
  ledgers listed on this page. Add a root document only when it is a new
  program or ledger of that kind.
- No file under `docs/` is a daily note, progress report, audit result,
  integration record, PR handoff or status summary. File names and headings
  carry no dates. Integration and review history belongs in the GitHub pull
  request or issue.
- Research results live in the papers, Lean library and code. No research
  note repeats what the code, the Lean library or the papers contain; extend
  the paper or the package README instead.
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
  the [observables register](registers/OBSERVABLES.json),
  the [falsification program](OPH_FALSIFICATION_PROGRAM.md), and the [paper index](../paper/).
- **Build and test:** begin with the repository [reproduction
  guide](../REPRODUCE.md), the [executable model](../model/candidate_laws/contrast_memory_spherical_v1/README.md),
  [executable evidence](../code/), and [Lean
  formalization](../Lean/). Follow the [numerical audit
  procedure](policies/NUMERICAL_AUDIT.md) when changing finite numerical
  evidence or addressing review findings.

## Basis And Paper Set

- [The OPH Basis](BASIS.md) distinguishes the closure hypothesis,
  internal-realization policy, agreement, adopted state-selection and
  response-endogeneity principles, and executable specification, and explains
  how the theorem-library premises A1, A2 and A3 relate to them.
- [The Paper Set](PAPER_SET.md) records the role of each paper.
- [Executable Hypothesis](OPH_EXECUTABLE_HYPOTHESIS.md) states what the
  executable model has to supply and how it is tested on observer records.

## Canonical References

- [Interpretation and Working Terms](INTERPRETATION_AND_TERMS.md) distinguishes
  complete-history closure, intrinsic computation, emergent-physics proposals,
  mathematical fixed points and the different uses of time.
- [Premise Reference](AXIOM_REFERENCE.md) gives the formal statements of the
  theorem-library premises A1, A2 and A3. The [basis](BASIS.md) places them
  among the principles of the theory.
- [Frozen-Prediction Ladder](FROZEN_PREDICTION_LADDER.md) records tests whose
  conditions are fixed before the comparison data are examined.
- [Empirical Status](EMPIRICAL_STATUS.md) states the observable the
  executable model reproduces and the theory's comparisons with public data.
- [OPH Falsification Program](OPH_FALSIFICATION_PROGRAM.md) lists only mature
  mathematical and realized-branch falsifiers.
- [Programme-wide obstruction search](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
  seeks inconsistency, joint-realization impossibility and limits on explanatory
  unification, with explicit exclusion and stopping criteria.

Open work on the executable model, its tests and existing-evidence audits is
tracked in the [observables register](registers/OBSERVABLES.json). Only
evidence changes a status; administrative closure of an issue does not
discharge a physical premise.

## Executable model

The [executable model](../model/README.md) holds the fixed carrier, the
spherical echosahedral federation, and the candidate laws that run on it. The
official [scorecard](../model/evaluation/SCORECARD.md) states the criterion
for the correct native repair law: it runs on the fixed federation and fulfils
the premises of the finite theorems together. [CANDIDATE.md](../model/candidate_laws/CANDIDATE.md)
lists the candidates with their status, and `python -m model evaluate` scores
each of them against the scorecard. The main candidate is
[Interval Echo](../model/candidate_laws/interval_echo/README.md), with its
[specification](../model/candidate_laws/interval_echo/SPEC.md); its seam interaction is the reference candidate
[A-contrast-memory-spherical-v1](../model/candidate_laws/contrast_memory_spherical_v1/README.md), with its
[specification](../model/candidate_laws/contrast_memory_spherical_v1/SPEC.md)
and [binding theorem](../model/candidate_laws/contrast_memory_spherical_v1/PROOF.md);
the [executable-model paper](../paper/executable_observer_construction.pdf)
gives the full account, and the [observables register](registers/OBSERVABLES.json)
records what the candidates reproduce. The state of the simulation is no
physical space: every test compares distributions of admissible observer
records with what observers in our world record.

Other code packages contain conditional mathematical components with explicit
equations and proof domains. [Established Mathematics for Observer
Agreement](AGREEMENT_RESEARCH_FOUNDATIONS.md) describes the general theorem
reuse requirements.

## Registers

[`registers/OBSERVABLES.json`](registers/OBSERVABLES.json) is the one
completion register. Its rows of kind `observable` list what the theory must
reproduce; its rows of kind `premise`, with ids `PR-nn`, list the premises of
the finite theorems. Every row carries one of four statuses, whose definitions
the file states, with the claims, papers and code behind it.
`python tools/check_observables.py` validates the file, and the public status
page reads it directly. The [five-ledger map](../claims/README.md#the-five-ledgers)
groups the machine sources into the claims, premises, predictions-and-tests,
evidence and papers ledgers.

## Research Notes

Results live in the papers, the Lean library and the code. The one note kept
here is cited by the released papers:

- [OPH predictions and public-data comparisons](research/PUBLIC_DATA_PREDICTION_AUDIT.md)
  places the fixed predictions against public data: the primordial tilt at
  the lower 0.612 percent tail of ACT DR6 P-ACT-LB, the coupling and the
  photon dispersion.

The audit notes of the code packages sit beside their code, for example
[finite geometry: support, refinement and oriented caps](../code/geometry/GEOMETRY_RECEIPT_AUDIT.md).

**Defect classes.** The
[defect-class certificate](../code/defect_particle_classes/README.md)
identifies each protected defect class with the centre class `6Y mod 6` of a
Standard Model multiplet.

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

## Evidence, Data And Writing Policies

- [Hadron Data Policy](policies/HADRON.md) defines provenance and promotion rules for the hadronic pipeline.
- [Writing Style Guide](policies/STYLE_GUIDE.md) defines the project’s reader-facing prose conventions.
- [OpenWave Submission Guidelines](OPENWAVE_SUBMISSION_GUIDELINES.md) bind research,
  audits and filings for the OpenWave comparison table: the platform's bar, the
  form of every test, the research cautions, the filing gates and the mandatory
  submission checklist.
- Simulator outputs have no theorem or empirical-promotion authority in the
  paper stack. A simulator count, replay, fit, or diagnostic closes no
  analytic or physical receipt; theorem status needs analytic proof packets,
  and empirical promotion needs public physical evidence bundles.

## Reader Support

- [Common Objections](COMMON_OBJECTIONS.md) answers recurring technical and conceptual criticisms.
- [Contributed background notes](../contributions/README.md) record standard mathematics that OPH
  work can reference; they make no OPH claim.
