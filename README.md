# Observer Patch Holography

**OPH proposes that the Universe is what its observers agree on, and builds that idea as a program you can run.**

<p align="center">
  <img src="assets/oph-schematic.png" alt="Observer Patch Holography schematic. Primitive observers on an icosahedral carrier read, record and act. Overlapping observer patches on a holographic screen agree on shared records, which become shared facts. From inside, one patch encodes a whole history: an emergent spacetime with three space dimensions and one time, where gravity appears as geometry." width="100%">
</p>
<p align="center"><b><a href="https://philpapers.org/rec/MUEFOC">Synthesis paper</a> · <a href="https://oph-book.floatingpragma.ai/">Book</a> · <a href="https://floatingpragma.ai/faq/">Reality FAQ</a> · <a href="https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a">Notebook</a> · <a href="https://particles.floatingpragma.ai/">Simulation</a></b></p>

<p align="center">
<a href="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/ci.yml"><img src="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
<a href="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean.yml"><img src="https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean.yml/badge.svg" alt="Lean"></a>
</p>

## The idea

The core intuition behind OPH is that reality is the self-answering question:
“What would a universe that creates itself consistently be like?” The key word
is *consistently*: mathematical and logical consistency, free of contradiction.
OPH aims to show that the universe it constructs on this basis is precisely
the one we inhabit.

Small observers keep local records and compare them with their neighbours.
The executable construction specifies their matter, fields, preparation and
measurement law before a run. Whatever survives every
comparison is public fact. OPH asks whether such a process, fully specified
and executed, produces observers whose shared records look like our world:
three-dimensional space, quantum measurement, the Standard Model's forces and
particles, and gravity.

The Universe reads itself into existence. Observers recover the rules behind
their records and build the process that realizes that same history, so
nothing lies outside the loop. The [basis](docs/BASIS.md) states this closure
hypothesis in full.

## What comes out

The search for the executable model's specification is ongoing, and it is
the project's current main priority.

The [current executable model](model/SPEC.md) is a calibrated test
candidate, not yet a production-level model. It joins complete native matter,
charged and neutral lepton transport, weak currents, electric and magnetic
fields, and local registers that respond to total source energy. Each vertex
has its own volume, capacity and memory register. Exact checks preserve
charge and Gauss constraints and verify source-sensitive coherent readout.
Disconnected components evolve independently. The full native operator is
checked through exact polynomials and basis columns; small propagated
diagnostics have separately stated carriers.

Its particle masses are inputs, not outputs. The quark and hadron
coefficients are calibrated to the measured proton, neutron and charged-pion
masses, the electron mass sets the energy unit, and the charged-lepton,
neutrino and quark-family values are supplied. OPH derives no particle mass
so far.

The theorem library also establishes conditional gauge and geometry results,
quantum-record identities and finite spectra on their specified sources.
A common interacting continuum, physical particle identification and the
Einstein theorem's source premises require their own demonstrations.

## How it works

One fixed Hamiltonian evolves matter, fields and response registers. A
complete detector retains every outcome, and its measurement work is
separate from energy-conserving free evolution. The creator supplies the
law and initial preparation; simulated observers need not install their
own interactions. The [executable hypothesis](docs/OPH_EXECUTABLE_HYPOTHESIS.md)
states the physical test: a fully specified construction must reproduce
what observers measure.

## Mathematics you can check

The [Lean library](Lean/) holds more than 12900 public theorems and lemmas
with no admitted proofs, covering consensus, quantum records, the gauge
structure, thermodynamics and the Koide identity. Every computational result
ships with a receipt that binds its inputs and replays.

## Tests fixed in advance

Predictions and their rejection rules are frozen before the data, with
OpenTimestamps proofs anchored in Bitcoin, in the
[frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md). The
[postdiction ledger](docs/POSTDICTION_LEDGER.md) and the
[falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) state how the
construction can fail.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
python3 -m model
```

This checks the common source and its complete finite diagnostics. The
[model guide](model/README.md) explains the full source replay and the
separate component proofs.

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main paper. Cite its released version, the
  [PhilPapers record](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer](paper/executable_observer_construction.pdf)**
  specifies the program and its executed components.
- **The theorem library** holds the long proofs:
  [consensus](paper/reality_as_consensus_protocol.pdf),
  [observers and records](paper/observers_are_all_you_need.pdf),
  [spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf),
  [gauge structure](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf)
  and [particles](paper/deriving_the_particle_zoo_from_observer_consistency.pdf).
  The [paper index](paper/README.md) lists every paper.

## Where things are

| Path | What it holds |
| --- | --- |
| [model/](model/README.md) | The executable observer: implementation, tests and receipts |
| [flagship/](flagship/) | The main paper |
| [paper/](paper/README.md) | The construction paper and the theorem library |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Companion papers |
| [Lean/](Lean/) | The machine-checked library |
| [code/](code/README.md) | Exact calculations and verifiers behind the papers |
| [evidence/](evidence/README.md) | Retained evidence packets |
| [claims/](claims/README.md) | Claims and their dependencies |
| [tracking/](tracking/) | Premise, observation, prediction and instrument registers |
| [docs/](docs/README.md) | The basis, working terms, policies and ledgers |
| [book/](book/) | The book |

## Reproduce and verify

The [reproduction guide](REPRODUCE.md) runs the full verification suite from
a clean clone:

```sh
python -m pip install -r model/requirements.txt
python -m model
python -m pytest -q -W error model
```

The [book](book/), the textbooks at [learn.floatingpragma.ai](https://learn.floatingpragma.ai/)
and the [interactive simulation](https://simulation.floatingpragma.io/) are
the gentlest introductions.

## Contributing

Contributions should connect one more piece of observed physics to the
running observer, or expose a failure in it. The
[construction guide](model/CONSTRUCTION.md), the
[model specification](model/SPEC.md) and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out the work. Main-branch PDFs are source previews; deposited releases
keep their original bytes.

Ce README existe aussi [en français](README_FR.md).

## License

Software, including Lean and the model code, uses [Apache-2.0](code/LICENSE).
Papers, documentation, figures, data and ledgers use [CC BY-NC-SA 4.0](LICENSE),
with upstream Particle Data Group terms for the packaged particle values.
The [license map](LICENSE) gives the directory-specific terms.

<details>
<summary>Related projects and funding</summary>

[Cadence](https://github.com/muellerberndt/cadence) is the related learning-machine
project. Its capability evidence is separate from OPH's physics evidence.
[Pragma Research](https://floatingpragma.ai/) brings the projects together.

OPH's official TEO (Theory Of Everything) memecoin was deployed by the Base
community and lives on Base at
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Its trading fees are routed to Bernhard Mueller to fund the project.
The token has no role in the theory. Other tokens using the OPH or TEO name
are unrelated. Nothing here is financial advice.

</details>

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
