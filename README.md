# Observer Patch Holography

**OPH proposes that the Universe is what its observers agree on, and builds that idea as a program you can run.**

<p align="center">
  <img src="assets/oph-schematic.png" alt="Observer Patch Holography schematic. Primitive observers on an icosahedral carrier read, record and act. Overlapping observer patches on a holographic screen agree on shared records, which become shared facts. From inside, one patch encodes a whole history: an emergent spacetime with three space dimensions and one time, where gravity appears as geometry." width="100%">
</p>
<p align="center"><b><a href="https://status.floatingpragma.ai/">Research status</a> · <a href="https://philpapers.org/rec/MUEFOC">Synthesis paper</a> · <a href="https://doi.org/10.5281/zenodo.24191882">Paper stack</a> · <a href="https://oph-book.floatingpragma.ai/">Book</a> · <a href="https://floatingpragma.ai/faq/">Reality FAQ</a> · <a href="https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a">Notebook</a> · <a href="https://floatingpragma.ai/oph/simulation">Simulation</a></b></p>

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

## The search for the repair law

The open problem of the project is the native repair law: the one local rule
that every observer follows when it compares its records with a neighbour,
which says how a disagreement is repaired and what the repair leaves behind in
memory. Everything else is in place. The finite theorems say how three
directions of space, the rules of quantum measurement, the forces of the
Standard Model and Einstein's gravity come about once their premises hold. So
the law has one job: fulfil those premises, all of them together, on one
process.

The official [scorecard](model/evaluation/SCORECARD.md) makes this a precise target:

- **Fixed:** the carrier, a sphere of echosahedral cells in which each cell is
  one observer and its twelve ports are its interfaces, and the premises of the
  finite theorems, listed in the [observables register](docs/registers/OBSERVABLES.json).
- **Free:** anything a designer can build in the simulating universe may be
  calibrated so that the observers' records match ours: couplings,
  preparations, feedback, readouts.
- **Forbidden:** putting in by hand what the theorems derive or what observers
  measure, such as particle masses, the number of space dimensions or a
  spatial lattice.
- **Together:** where two theorems talk about the same thing, the law supplies
  one record, one clock and one response for both. The stress that observers
  read from their clocks is the stress carried by the excitations that move,
  and the gauge currents act on those same excitations.

Once a law and its calibrations are frozen, held observer records test the predicted field, clocks, OpenWave phenomena and the identity of every Standard Model particle. The inventory requires measured distinctions among species, including their spins, charges and interactions.

Two simulated observers can send a charge carrying three coherent winding states. Changing how they name and phase their local arms reproduces all 90 measured count banks. A separate probe at the receiving port distinguishes like from unlike charges. From measurements on one sign, a two-state quantum model predicts 48 fresh detector readings on both signs. A second apparatus using the charge itself as an interferometer recovers the same winding state. Its lowest propagating band approaches 16 in guide energy units, and a full turn of the detected winding gives positive phase. In a separate finite quantized Maxwell model, charged source and receiver circuits distinguish the two windings of a one-quantum packet through local measurements. On the Einstein branch, tidal arm readings distinguish two classical gravitational-wave polarizations that leave all ten Einstein readings unchanged. [Three-state charge readout](model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/README.md) · [Charged families in two gauge frames](model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/face_frame/README.md) · [Photon](model/candidate_laws/interval_echo/openwave/photon_identity/README.md) · [Graviton](model/candidate_laws/interval_echo/openwave/graviton_identity/README.md)

**Can observers predict heat flow before heating either end?** We connected four small energy registers and watched energy move back and forth at equal temperatures. Those fluctuations predicted how much heat would flow after warming one end. Reversing the temperature difference gave the predicted response. Disconnecting one end kept the fluctuations and stopped the steady flow. [Measurements, error bounds and replay](model/experiments/interval_echo_consensus/refinement_cooling/dimer_band/heat_transport/README.md)

**Can one quantum field tell us how another one gravitates?** In a simulation with Einstein gravity built in, we calibrated on one field and predicted what radar and a boundary-area meter would read for a second. With two coupled fields, changing their shared quantum state changed both readings, even though each field looked the same on its own. [Experiment, controls and replay](code/matched_einstein_limit/quantum_scalar_radar/README.md)

**Does the repair network produce that gravity?** A flying probe can tell when two nearby seams share a memory. An ordinary matter source and a seam memory can store the same extra energy and give very different probe readings. We can measure the scattering. The next gravity test is to measure directional stress and curvature independently for those same sources. [Native source and probe experiment](model/experiments/interval_echo_einstein_station/shared_geometry_probe/README.md)

Every candidate law lives in [model/candidate_laws](model/candidate_laws/CANDIDATE.md)
and runs on the same fixed carrier. One command scores all of them against
the official [scorecard](model/evaluation/SCORECARD.md). To display it, run

```sh
python -m model evaluate
```

The main candidate is [Interval Echo](model/candidate_laws/interval_echo/README.md). It equips the contrast-memory interaction with a clock carried and read by a finite patch. The contrast-memory reference and the scalar channel appear alongside it in the common scorecard. New candidates are welcome: a candidate is a folder with a
`candidate.py`, and [CANDIDATE.md](model/candidate_laws/CANDIDATE.md#adding-a-candidate)
says what it declares.

## Research focus

Our current focus is the [executable model](model/README.md) and its repair
law. The hardware is the fixed sphere of echosahedral cells. The open question
is which local repair law, run on that hardware, gives observers whose shared
records reproduce the observed Universe. The main candidate, Interval Echo, runs the contrast-memory seam repair and carries a finite observer’s clock through local handoffs. The theory around the model is in place: the closure
principle, the theorem library and the finite constructions below.

Progress is tracked in two places. The [scorecard](model/evaluation/SCORECARD.md)
grades each candidate law against the requirements. The
[research status page](https://status.floatingpragma.ai/) lists every
observable OPH must reproduce and every premise of its finite theorems the
model has to confirm, and links each one to its Lean proofs, model evidence
and papers.

**What the observers inside experience.** The state of the simulation does
not contain the physical Universe. Its ports, seams and cells are the hardware
of a computation, and they are no space. Space, time, fields and particles are
how the observers inside describe the records they agree on. The native cell has a rank-three slow response. On Interval Echo’s programmed scalar instrument, observers recovered three binary spatial directions from eight anonymously labeled stations. The chart predicted 24 held clock-stopped detector distributions with a largest absolute discrepancy of 0.001022. Refining that routing program gives a flat 3+1D free-wave limit on bounded times and momentum bands. On the freely evolving federation, a shared three-dimensional metric and an interacting four-dimensional conformal field theory are open. [Scalar chart and continuum theorem](model/experiments/seam_program/CLOCKED_FIELD.md) · [Geometry evidence](model/evaluation/EVIDENCE.md#patch-05) · [Gravity requirements](model/evaluation/EVIDENCE.md#g1)

A change in the alternating orientation of loaded port pairs travels between cells when the connecting memories stay coherent. Observers measured its wavelength using one set of detector readings and predicted another. On three reserved boards, every measured spatial profile passed while the wavelengths used to launch the waves failed. Two local coupling measurements predicted 22 frequencies within 0.000575 in dimer-clock units. [Experiment and independent spatial test](model/candidate_laws/interval_echo/openwave/submission/row_30_confirmation/README.md)

A history can be a common terminal record of local repairs. If a nonnegative disagreement measure decreases by a fixed amount at every accepted repair, its initial value bounds the number of accepted repairs. A polynomial running-time bound also requires bounds on that initial value and the cost of finding and serving the repairs. Uniqueness requires overlapping repairs to reach a common endpoint on the public records. Interval Echo tests that condition on a two-seam fork. A polynomial bound for its unique federation-wide history is open. [Reference specification](model/candidate_laws/contrast_memory_spherical_v1/SPEC.md) · [Interval Echo repair-order test](model/candidate_laws/interval_echo/CONSENSUS.md)

**How close is the electron?** In the written-field extension, a moving defect carries one unit of signed charge, and a separate probe feels a confined source 22 cells away. In registered Interval Echo, a defect carries a qubit across two cell boundaries. Turning it all the way around gives positive phase. The qubit records an orbital orientation. The electron search needs a travelling spin-½ on the same state as the charge, followed by an independent rest-energy measurement. [Charge and field](model/candidate_laws/interval_echo/openwave/electron_field/README.md) · [Rotation and exchange](model/experiments/interval_echo_angular/charged_orientation/README.md)

**Where to pick up the search.** For matter, start with the preparations that carry a qubit and pass the local colour filter at both ends. The first signal to enter a neighbouring cell has one-quarter colour-triplet weight, so a candidate that stays colour neutral throughout flight needs a different preparation or representation. For gravity, reuse the ordinary-source calibration and its equal-energy seam control. Measure the curvature with a separate distance or area instrument before assigning the probe response to gravity. Each experiment links its protocol, raw counts and replay commands. [Replay the experiments](REPRODUCE.md#replaying-the-matter-and-gravity-experiments) · [Colour and transport](model/experiments/interval_echo_gauge/singlet_flight/README.md) · [Matter criterion](model/evaluation/EVIDENCE.md#ft-matter) · [Universal coupling criterion](model/evaluation/EVIDENCE.md#g6)

The contrast-memory reference and the scalar channel supplied the following observer records. Interval Echo uses the same contrast-memory interaction and a clock filter calibrated to a guided reunion comparison.

- **Particles that hold together.** A pair of loads written on one seam is
  read as a pair at every later reading, with a probability that a theorem
  bounds on carriers of any size, and it moves from seam to seam. The same
  law without the agreement energy loses the pair.
- **A clock inside an object.** An object held on four seams swings at one
  fixed pitch after a pulse. An observer reads that pitch from its own
  counts, in units of its own reference, and reads the same number a second
  way from the two levels of the object.
- **A record that persists on empty ports.** A note left in the memory of
  a seam with nothing around it is never erased by the law. It holds no
  material load, its observers read its level, and matter that reaches it sticks to
  it. Its physical electric neutrality is a separate open question.
- **Looking changes what comes later.** Asking "is a pair there?" and asking
  "where is it?" give the same first answer and different later records.
- **Twelve force currents.** With one declared port service, the response of
  a cell's ports has twelve generators with the algebra of the Standard
  Model's forces.

These outputs plug into the theorem library wherever an output has the form
of a premise: the cell is the carrier of the three-direction theorem, the
readings are projective instruments of the measurement theorems, and the
observers compare clocks through shared records. The
[scorecard](model/evaluation/SCORECARD.md) records which premises each
candidate fulfils, and the [research status page](https://status.floatingpragma.ai/)
shows what the program reproduces.

## What comes out

Under stated premises, the theorem library and the finite constructions give:

- **Spacetime and fields as effective descriptions.** Coarse-grained, the
  observers' shared records are described by relativistic field theories in
  three space dimensions and one time, together with gravity. The carrier's
  slow response is exactly three-dimensional, and a Lorentzian (1+3)-dimensional
  field theory reproduces the observers' joint readouts.
- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach one shared public result whatever their order.
  In an Interval Echo run with 240 ports, two opposite orders each made 74
  repairs and reached the same minimum imbalance. The two runs could write
  different seams. [Observer counts](model/experiments/interval_echo_consensus/routed_counts/)
  · [Consensus paper](paper/reality_as_consensus_protocol.pdf).
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

A fixed law acts on the registers of neighbouring ports. Under the
contrast-memory law it moves loads between them and lowers the energy of a
seam whose two observers agree. Readers keep every outcome, the outcome
"no pair" included. Each observer counts its own events, and no common clock
is used. The creator supplies the law and the preparation; simulated
observers need not install their own interactions. The [executable hypothesis](docs/OPH_EXECUTABLE_HYPOTHESIS.md)
states the physical test: a fully specified construction must reproduce
what observers measure.

## Mathematics you can check

The [Lean library](Lean/) holds more than 13400 public theorems and lemmas
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
python3 -m model                          # the main candidate against its stored results
python3 -m model run --law scalar_channel # another candidate on the same carrier
python3 -m model evaluate                 # every candidate against the requirements
```

The [model guide](model/README.md) explains the layout: the fixed carrier, the
candidate laws, their evaluation and the experiments that observers run.

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main paper. Cite its released version, the
  [PhilPapers record](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer](paper/executable_observer_construction.pdf)**
  specifies A-contrast-memory-spherical-v1, proves its binding theorem and reports the
  records its observers agree on.
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
| [model/](model/README.md) | The executable model: the fixed carrier, the candidate laws, their evaluation and the experiments |
| [flagship/](flagship/) | The main paper |
| [paper/](paper/README.md) | The construction paper and the theorem library |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Companion papers |
| [Lean/](Lean/) | The machine-checked library |
| [code/](code/README.md) | Exact calculations and verifiers behind the papers |
| [evidence/](evidence/README.md) | Retained evidence packets |
| [claims/](claims/README.md) | Claims and their dependencies |
| [docs/](docs/README.md) | The basis, working terms, policies, ledgers and the observables register |
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

The most valuable contribution is a candidate law that fulfils more premises
of the finite theorems on the fixed carrier, or a failure of a candidate shown
on its observers' records. The [scorecard](model/evaluation/SCORECARD.md)
and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out the work.
Main-branch PDFs are source previews; deposited releases
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
