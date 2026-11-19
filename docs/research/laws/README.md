# The search for the repair and memory law

OPH needs one native repair law. It is the local rule every observer follows
when it compares its records with a neighbour: it says how a disagreement is
repaired and what the repair leaves behind in memory. The theorem library
states what follows once a law supplies certain outputs, from three
directions of space to quantum measurement and gravity. This folder documents
the search for a law that supplies them on one executable process.

The [executable model](../../../model/polynomial_sv/README.md) runs two
candidate laws, the contrast-memory law and the
[scalar channel](../../../model/scalar_channel/README.md). The search also
tests laws outside the model, such as the reflected-routing law of the first
step. A result in this folder concerns the law and the preparation it names.
It changes no mark in the
[observation ledger](../../registers/OBSERVATION_LEDGER_V3.md), whose marks
come from the model's own code and receipts.

## Contents

| Document | What it holds |
| --- | --- |
| [Plan](PLAN.md) | The search program: what counts as a candidate, the experiments in dependency order, how candidates are enumerated, and when a family of laws is retired. |
| [Journal](JOURNAL.md) | One entry per step: what the step measured, which classes of laws it excludes, favours or leaves undecided, and what the following step tests. |
| [Clocks and spacetime from observer records](relativistic-observer-records/report.md) | Step 1. A moving clock carried by a freely propagating wave, compared with its rest calibration, and a three-dimensional scalar field measured by equipped observers, with its continuum limit. [PDF](relativistic-observer-records/report.pdf) |
| [A race between the addressed relay and the free wave](signal-race/report.md) | Step 2. Can an inhabitant send a freshly chosen bit ahead of the smooth scalar cone? The addressed relay carries every bit early; the free wave carries them late. |
| [Two clocks on the same journey](second-clock/report.md) | Step 3. One calibrated metric predicts a private-mode beat clock, an optical clock on a transverse ring and a fresh signal on five routes. |
| [A clock held together by an interaction](interaction-clock/report.md) | Step 4. A compact bound pair keeps proper time, the finite law predicts its deviation from the smooth metric ([correction](interaction-clock/correction.md)), and separated constituents never form the pair. |

## The simulator

The [simulator](../../../simulator/README.md) executes the reflected-routing
family with its members and controls, records what observers read, and
replays the records in an offline web app, including the scalar continuum as
an observer at the centre of the chart receives it.

## How a step is added

Each step has a folder named after its subject. The folder holds the report,
its figures and, where one exists, the PDF and the script that typesets it. The step gets one entry
in the journal, which states how the step changes the standing of each class
of laws it tested. The plan changes when a step changes the order of the
experiments or the definition of a candidate. Each report names the
executable, the records and the revision it used.

A claim in these documents is a claim about observers' records: who holds
the record, which reading is computed from it, and what the law predicts for
that reading. Statements about the simulator state are steps of a derivation
or an audit and are labelled as such.
