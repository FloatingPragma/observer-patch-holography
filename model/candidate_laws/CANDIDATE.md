# Candidate laws

The table gives a reproducible measurement and the scientific reason for each candidate’s role. Interval Echo is active for the current experiment on clocked patches and bound matter. The common [scorecard](../evaluation/SCORECARD.md) compares all laws on the same requirements, and the [evidence ledger](../evaluation/EVIDENCE.md) links every assessment to its source.

The electron search uses an experimental extension with written seam fields. It supplies a candidate charged constituent whose flux, distant-probe response and orbital transport are measured. The scorecard lists its evidence separately from results on the registered Interval Echo law. [Electron candidate and its exchange controls](interval_echo/openwave/electron_field/README.md)

<!-- CANDIDATES:BEGIN -->
| Candidate | Role | Measured example | Reason for its role |
| --- | --- | --- | --- |
| [Interval Echo](interval_echo/README.md) | main | Traveller/home elapsed ratio on alternating 2- and 6-service handoffs **0.955354 ± 0.0019 (SE)**. [Data](interval_echo/RESULTS.json) | Full federation generator, a transported and jointly read patch clock, and a current instrument on the bound pair. The 19-port ring comparator contradicts the common rest-energy and inertia assignment for its free clusters. |
| [Contrast memory](contrast_memory_spherical_v1/README.md) | reference | Largest probability of losing pair adjacency across nine sampled waits on 20 cells at κ = 8 **0.00962164**. [Data](contrast_memory_spherical_v1/RESULTS.json) | Reference for the shared seam dynamics, with an all-size binding theorem and local reader controls. The 19-port ring resolves different rest-energy and kinetic ratios for free clusters. |
| [Private record walk](private_record_walk/README.md) | experimental | Home-pair probability with eight sharp comparison signs at exposure 37 **0.978566**. [Data](private_record_walk/results.json) | Each load exchange and each enabled private sign reversal is one edge of the same unweighted graph. |
| [Scalar channel](scalar_channel/README.md) | second | Cage clock frequency in its supplied reference units, equal-step schedule **4.89898 ± 2.1e-06 (SE)**. [Data](../experiments/cage/results.json) | Independent test of a stable cage clock and finite mode dispersion. Its registered generator covers one cell. The larger cage construction has its own apparatus. |
<!-- CANDIDATES:END -->

## The candidates

**[Interval Echo](interval_echo/README.md) is the main candidate under test.** It equips the contrast-memory seam with a finite observer clock. Its hand, elapsed count and interval memory move with the patch. At a confirmed handoff the patch uses its own preceding repair events to decide whether to post the closing mark. The guide, gate selection, prepared hand and cadence are supplied. [Specification](interval_echo/SPEC.md) · [Reunion data](interval_echo/RESULTS.json)

Joint native gateway records resolve a stationary orbital phase. On monitored regions, independently readable exterior word capacity cancels held interior entropy variations. A physical area reading at fixed proper volume is required for the Einstein-tower interpretation. [Interval Echo gravity measurements](interval_echo/README.md#can-observers-measure-gravity)

**[A-contrast-memory-spherical-v1](contrast_memory_spherical_v1/README.md) is the reference candidate.** Loads
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
| `PARTICLE_IDENTIFICATIONS` (optional) | field-type keys from the [particle inventory](../evaluation/particles.json), each with its grade, discriminating observation and evidence |
| `run(argv)` | the candidate's own replay, called by `python -m model run --law <name>` |

Exactly one candidate has the role `main`; `python -m model` runs it.
