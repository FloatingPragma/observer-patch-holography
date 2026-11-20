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
The executable construction specifies their carrier, law, preparation and
readers before a run. The Universe evolves through
local disagreement and repair, and the records that survive every comparison
are its public facts. OPH asks whether such a process, fully specified
and executed, produces observers whose shared records look like our world:
three-dimensional space, quantum measurement, the Standard Model's forces and
particles, and gravity.

The Universe reads itself into existence. Observers recover the rules behind
their records and build the process that realizes that same history, so
nothing lies outside the loop. The event where observers have reverse-engineered
and built that process is a mathematically guaranteed fixed point of the
history, and it exists in our subjective future just as our past exists. The
[basis](docs/BASIS.md) states this closure hypothesis in full.

## Research focus

The theory around the executable model is finished: the closure principle,
the theorem library and the finite constructions below are in place. Under
their premises, observers inside such a computation describe their agreed
records as a relativistic spacetime with three dimensions of space and one of
time, with quantum measurement, the Standard Model's forces and gravity.

What is missing is the native repair law, and the search has narrowed it to
three candidate laws, which are being tested. All three move each constituent
by reflected routing, bind two constituents that meet into a pair whose two
rails form a clock, and form that pair through a capture that sends away a
marker. They differ only in how a newborn pair opens its motion: on a cosine
profile over 32 updates, on the same profile over 128 updates, or linearly
over 128 updates. The [candidate laws](model/README.md) page states what has
been proven about them, and the [simulator](simulator/README.md) runs them
and replays every experiment.

## Finding the repair and memory law

The main focus of the research is to identify the correct repair and memory
law. This is the native repair law, the simple law that governs the
interactions of observers in the early Universe: the one local rule that
every observer follows when it compares its records with a neighbour, which
says how a disagreement is repaired and what the repair leaves behind in
memory. Repair decides what becomes a shared fact. Memory matters just as
much: in the candidate laws a capture leaves a marker that records where and
when a pair was born, and the age a newborn pair carries decides which way
it travels.

The couplings of the candidate laws are supplied, and the other observables
in the papers, from charge and the particle spectrum to gravity, are work in
progress on the same process. Any part of the executable observer may change
as the research goes on.

**Why it matters.** The theorem library states what follows once a law
supplies certain outputs: three directions of space, the rules of quantum
measurement, the forces of the Standard Model, Einstein's gravity. A law has
to supply all of them at once, on one process, with every check and every
repair local. With the correct law, the physics we measure becomes something
you compute: run the law, let its observers keep their records, and read off
what they find. A wrong law gives itself away the same way, because its
observers record something we do not. That makes the search for the law the
sharpest test of the whole idea, and the rest of the project serves it. The
[model overview](model/README.md) lists what the candidates reproduce and what
is open.

## What comes out

Under stated premises, the theorem library and the finite constructions give:

- **Spacetime and fields as effective descriptions.** Coarse-grained, the
  observers' shared records are described by relativistic field theories in
  three space dimensions and one time, together with gravity. The carrier's
  slow response is exactly three-dimensional, and a Lorentzian (1+3)-dimensional
  field theory reproduces the observers' joint readouts.
- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach one shared public result whatever their order.
  [Consensus paper](paper/reality_as_consensus_protocol.pdf).
- **The rules of quantum measurement.** Completed records obey the Born rule,
  repeatable measurement and the Tsirelson bound on Bell correlations, with a
  state that attains it.
  [Finite quantum event algebras](extra/machine_checked_finite_event_algebras.pdf).
- **The Standard Model's forces.** Under stated premises on complete reversible
  response and observer transport, the twelve-port response algebra has the
  Standard Model gauge Lie type, and port phases realize the global group.
  [Gauge derivation](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **Quarks and leptons.** The defects that repair leaves behind carry exactly
  the six centre classes of the Standard Model particles, and each generation
  comes out with the observed charges: three colours of up and down quarks, the
  electron and its neutrino. A specified matter representation gives an
  anomaly-free generation.
- **Thermodynamics.** Finite repair theorems give a Clausius inequality,
  fluctuation relations and an information Landauer inequality. Fair-coin seam
  repair is exact conditional resampling, so entropy never falls on binary
  populations while the complete tape keeps every bit.
  [Thermodynamics account](paper/observers_are_all_you_need.pdf).
- **Gravity.** A realized balance of records gives the Einstein equation.
  [Spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).

## How it works

**Observers come first. What OPH compares with nature is what observers read
from their own records, and what they agree on where their records overlap.
Measurement, consensus and the synchronizing of clocks come after the records.
No objective field and no global clock stand behind them: fields, particles,
spacetime and a common time are effective descriptions that observers agree to
use. The state of the program and its time parameter calculate those readouts.
They are no second universe beside them.**

A fixed law acts on the registers of neighbouring addresses. Under the
candidate laws it reflects and routes amplitude between them, binds two
constituents that meet into a pair, and records each binding in a marker
that leaves on its own wire. Readers keep every outcome, misses included.
Each observer counts its own events, and no common clock is used. The creator supplies the law and the preparation; simulated
observers need not install their own interactions. The [executable hypothesis](docs/OPH_EXECUTABLE_HYPOTHESIS.md)
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
[falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) states how
the construction can fail.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
python3 -m model                  # exactness checks of the candidate laws
python3 -m simulator run all      # every experiment, and the offline replay app
```

The first command checks each of the simulator's exact methods against the
literal update. The second recomputes every experiment of the search and
writes `simulator/runs/replay.html`, which replays them. The
[candidate laws](model/README.md) page explains both.

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main paper. Cite its released version, the
  [PhilPapers record](https://philpapers.org/rec/MUEFOC).
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
| [model/](model/README.md) | The candidate laws and what is proven about them |
| [simulator/](simulator/README.md) | The simulator of the candidate laws and its replay app |
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
python -m pytest -q -W error simulator model
```

The [book](book/), the textbooks at [learn.floatingpragma.ai](https://learn.floatingpragma.ai/)
and the [interactive simulation](https://simulation.floatingpragma.io/) are
the gentlest introductions.

## Contributing

Contributions should connect one more piece of observed physics to the
running observer, or expose a failure in it. The
[candidate laws](model/README.md), the
[simulator specification](simulator/SPEC.md) and the
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
