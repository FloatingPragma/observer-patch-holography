# Candidate laws

Each candidate for the native repair law has a folder here. Every candidate runs on the fixed carrier of
[executable_model/](../executable_model/README.md) and is scored against the [requirements](../REQUIREMENTS.md)
by `python -m model evaluate`. The table is written by that command; the [scorecard](../evaluation/SCORECARD.md)
has every check, and each candidate's `EVALUATION.json` has every premise with its evidence.

<!-- CANDIDATES:BEGIN -->
| Candidate | Role | Carrier checks | Law premises | Joint requirements |
| --- | --- | --- | --- | --- |
| [A-contrast-memory-spherical-v1](contrast_memory_spherical_v1/README.md) | main | 7 of 7 pass | 1 fulfilled, 0 partial, 58 open, 0 failed | 0 fulfilled, 0 partial, 4 open, 1 failed |
| [A-scalar-CAR-pair-transit-v2 (the scalar channel)](scalar_channel/README.md) | second | 6 of 7 pass | 1 fulfilled, 0 partial, 58 open, 0 failed | 0 fulfilled, 0 partial, 5 open, 0 failed |
<!-- CANDIDATES:END -->

## The candidates

**[A-contrast-memory-spherical-v1](contrast_memory_spherical_v1/README.md), the main candidate under test.** Loads
0, 1 or 2 sit on the ports and an eleven-state memory on every seam. A seam event writes the contrast of its two
loads into the memory, gives it a phase and erases it again; on the loads this is an exchange plus an agreement
energy. It runs on one cell and on the federations of 20, 80 and 320 cells, a theorem binds a pair of loads on all
of them, and the papers describe it. [Specification](contrast_memory_spherical_v1/SPEC.md) ·
[proof](contrast_memory_spherical_v1/PROOF.md)

**[The scalar channel](scalar_channel/README.md), the second candidate.** One object or none sits on every seam and
on every wedge, a wedge being two seams at one port. An object hops between a seam and a wedge, which sit at the
two opposite levels `+m` and `-m`. It is written on one cell; its cage and the experiments built on it are in
[experiments/](../experiments/README.md). [Proof](scalar_channel/PROOF.md)

## Adding a candidate

A candidate is a package `model/candidate_laws/<name>` with a `candidate.py` that declares:

| Name | What it is |
| --- | --- |
| `NAME`, `TITLE`, `ROLE`, `SUMMARY` | the identifier, the display name, `main` or another role, one sentence |
| `LABELLING` | the 30 seams of one cell in the candidate's port labelling |
| `PARAMETERS` | every setting with its value, its provenance (`designed`, `calibrated` or `derived`) and its quantity |
| `basis(car)`, `generator(car)` | the registers of the one-excitation sector and the law on them, for one cell and, where the law is written for them, the federations |
| `response(car)` | the operator whose slowest nonconstant level is the carrier response |
| `represent(rotation, car)` | the action of a port permutation on the registers |
| `PREMISES`, `JOINT` | the status and evidence for rows of the premise register and for the joint requirements |
| `run(argv)` | the candidate's own replay, called by `python -m model run --law <name>` |

Exactly one candidate has the role `main`; `python -m model` runs it.
