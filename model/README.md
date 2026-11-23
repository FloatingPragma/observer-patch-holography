# The executable model

The theory says the Universe is a very large network of very simple observers. Each one holds a few records and
compares them with its neighbours, and whenever two neighbours disagree they repair the difference between the two
of them. No observer sees the whole, and nothing stands outside. What the network settles on is its history, and
that history contains the observers who work out how the network is built.

This folder holds that network as a program. Its hardware is fixed: the spherical echosahedral federation. Its
software is a candidate for the native repair law, the simple local rule by which neighbouring observers settle
their disagreements. Finding the correct law is the open problem of the project, and the
official [scorecard](evaluation/SCORECARD.md) covers the axioms and finite theorem premises on one process, the patch-clock and field experiments, the 31 OpenWave phenomena, and identification of all Standard Model particle types, and it gives each law’s evidence and the scope of every measurement.

The written-field extension of Interval Echo supplies a strong candidate for the electron’s charge and field sector. An independent boundary reading follows the signed material load, and a confined source changes the signal at a probe 22 cells away. The orbital information that can accompany the load during transport can also be left behind during exchange. Spin and a physically calibrated rest energy require independent measurements on the same excitation. [Evidence, controls and replay](candidate_laws/interval_echo/openwave/electron_field/README.md)

## What is where

| Folder | What it holds |
| --- | --- |
| [executable_model/](executable_model/README.md) | everything fixed for every law: the echosahedral cell, the federations of 20, 80 and 320 cells, their rotations and incidence, and the registry of laws |
| [candidate_laws/](candidate_laws/CANDIDATE.md) | one folder per candidate law with its specification, implementation, tests and evaluation; [CANDIDATE.md](candidate_laws/CANDIDATE.md) lists them with their status |
| [evaluation/](evaluation/README.md) | the systematic evaluation of every candidate and its [scorecard](evaluation/SCORECARD.md) |
| [experiments/](experiments/README.md) | programs that observers run with a candidate law: preparations, instruments and extended constructions |
| [benchmarks/](benchmarks/README.md) | phenomenological benchmarks on supplied architectures; they are no candidates |

Evidence prepared for an OpenWave row, and every filing to OpenWave, follows the
[OpenWave submission guidelines](../docs/OPENWAVE_SUBMISSION_GUIDELINES.md).

## Run it

```sh
pip install -r model/requirements.txt
python -m model                                # replay Interval Echo on 20 cells and compare the reunion cohorts
python -m model --full                         # include 80/320 cells, patch records, the repair fork and port currents
python -m model --full --seed 42 --output /tmp/echo-42  # generate an independent data set
python -m model run --law scalar_channel       # another candidate
python -m model run --law interval_echo --clock-trip 8 8   # five-bit clock, local hands and shared reunion
python -m model list                           # the candidates and their roles
python -m model evaluate                       # every candidate against the requirements
python -m model evaluate --criterion G2        # evidence for one criterion across the candidates
python -m model carrier                        # the checks of the fixed carrier
python -m pytest -q model                      # every test
```

## What the candidates show

| Observable | What the observers read | Candidate | Ledger row |
| --- | --- | --- | --- |
| A clock carried by a patch | Both owners read a travelling nine-bit elapsed word and a home word on the same signed seam at reunion. The patch derives each handoff interval from its own repair record. | [Interval Echo](candidate_laws/interval_echo/SPEC.md) | finite guided comparison |
| Particle stability | Two loads written on one seam are read as a pair at every later reading, on every number of cells, by a theorem; the same law without the agreement energy loses them. A second object, the excitation of a [cage](experiments/cage/README.md), is kept for every interval and read by a local detector. | contrast-memory reference and scalar channel | OL-K5, attained |
| A clock inside an object | After a pulse the seam record of a cage swings at one line with no drive. An observer reads the line at 4.8990 of its own reference and reads the level difference of the two modes separately at the same value. | scalar channel | OL-K4, attained as a prepared clock |
| A stable seam record with zero material load | A seam memory that holds a record while its ports are empty stays that way under every event. Its observers read its level. Physical electric neutrality remains open. | contrast-memory reference | OL-K9, attained for the stored record |
| Waves of an object | Every frequency of the scalar channel obeys `ω² = m² + λ`, with `λ` a level of the seams. One point of the relation is read from records. | scalar channel | OL-K8, partial |
| What a reading leaves | Asking "is a pair there?" and asking "where is it?" give the same first answer and different later records. | contrast-memory reference | |

The last column names the row of the [observation ledger](../docs/registers/OBSERVATION_LEDGER_V3.md). Each result
is a statement about one construction, with everything its makers supplied listed beside it. Interval Echo shares the reference candidate’s seam operator and adds a patch clock. The scalar channel acts on a different state space. Each candidate has its own premise evaluation.

On the 20-cell Interval Echo federation, a patch in cell zero records whether
two initially separated loads reach the same seam and whether they return
after free release. Conditional on a positive first reading, 17 bookings
give simultaneous 99% lower bounds of eight positive gap frequencies and
36.72% of the prepared state's weight outside any two energy levels. A
single particle species can occupy many such levels.
[Public return records and all-size operator bound](experiments/interval_echo_encounter/README.md#how-many-frequencies-can-one-patch-measure)

## Supplied, and open

Supplied by the makers of the construction, and listed with every result: the couplings `κ` of the main candidate
and `m` of the scalar channel; emitters, pulses and their strengths; the reference registers, guides, gate selection, service cadence, readers and
apparatus of the experiments. Open, and work in progress: the premises the scorecard lists as open, a law that
fulfils them together, and the particle census of a frozen law.
