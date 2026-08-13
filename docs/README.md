# OPH Project Documentation Index

This directory contains project policies and indexes, architecture and custody
material, reproduction guidance, and canonical or generated status ledgers.
Scientific results live in the papers, Lean library, executable code, and
evidence artifacts; this directory is not a parallel result-publication layer.

If you are new to OPH, the strongest starting material lives outside this
directory. The [technical paper](../flagship/from_observer_consensus_to_standard_physics.pdf)
states the primary technical account. The [OPH textbooks](https://learn.floatingpragma.io/) work through
the basic derivations with the math taught along the way. The
[interactive simulation](https://simulation.floatingpragma.io/) renders real
run data so you can watch the universe assemble itself.

## Three Reading Routes

- **First encounter:** the [technical paper](../flagship/from_observer_consensus_to_standard_physics.pdf),
  [textbooks](https://learn.floatingpragma.io/), and the
  [simulation](https://simulation.floatingpragma.io/) above, then the
  [compact case](../extra/compact_proof_of_oph.pdf) and the repository
  [README](../README.md) from the three axioms through the twist.
- **Scientific audit:** use the [proof spine](PROOF_SPINE.md), the
  [closure issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure),
  the [falsification program](OPH_FALSIFICATION_PROGRAM.md), and the [paper
  index](../paper/).
- **Build and test:** begin with the repository [reproduction
  guide](../REPRODUCE.md), [executable evidence](../code/), and [Lean
  formalization](../Lean/).

## Reader Support

- [Common Objections](COMMON_OBJECTIONS.md) answers recurring technical and conceptual criticisms.

## Proof And Verification Maps

- Each quantitative closure condition is tracked as a
  [GitHub issue labeled `closure`](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure),
  with its evaluation boundary and required completion stated on the issue.
- [Proof Spine](PROOF_SPINE.md) gives the dependency route through the principal claims.
- [Premise Register](PREMISE_REGISTER_V3.md) names each V3 input and its
  consumers; the generated [Premise Discharge Queue](PREMISE_DISCHARGE_QUEUE_V3.md)
  turns every `remove` or `axiomatize` row into an auditable action item.
- [Architecture Replay Index](ARCHITECTURE_REPLAY_INDEX.md) shows which
  promotions and frozen targets survive the current architecture version;
  [Audit Custody](AUDIT_CUSTODY.md) pins the independent review history.
- [OPH Falsification Program](OPH_FALSIFICATION_PROGRAM.md) lists only mature mathematical and realized-branch falsifiers.

## Evidence And Data Policies

- [Hadron Data Policy](HADRON.md) defines provenance and promotion rules for the hadronic pipeline.
- [IBM Quantum Cloud Archive](IBM_QUANTUM_CLOUD.md) documents the reproducible hardware benchmark and its interpretation boundary.
- Simulator outputs have no theorem or empirical-promotion authority in the
  paper stack. A simulator count, replay, fit, or diagnostic closes no
  analytic or physical receipt; theorem status needs analytic proof packets,
  and empirical promotion needs public physical evidence bundles.

## Editorial Policy

- [Writing Style Guide](STYLE_GUIDE.md) defines the project’s reader-facing prose conventions.
