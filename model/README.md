# The executable model

The theory says the Universe is a very large network of very simple observers. Each one holds a few records and
compares them with its neighbours, and whenever two neighbours disagree they repair the difference between the two
of them. No observer sees the whole, and nothing stands outside. What the network settles on is its history, and
that history contains the observers who work out how the network is built.

This folder holds that network as a program. Its hardware is fixed: the spherical echosahedral federation. Its
software is a candidate for the native repair law, the simple local rule by which neighbouring observers settle
their disagreements. Finding the correct law is the open problem of the project, and the
official [scorecard](evaluation/SCORECARD.md) covers the axioms and finite theorem premises on one process, the patch-clock and field experiments, the 31 OpenWave phenomena, and identification of all Standard Model particle types, and it gives each law’s evidence and the scope of every measurement.

In the written-field extension a boundary reader follows a signed material load, and a probe detects a confined source 22 cells away. Orbital information can follow the load during transport and stay behind during exchange. On a separate eight-site interacting ring, two neighboring receivers record either signed material unit and recover spin entangled with the sender’s reference on the same trial. An electron assignment requires electric charge and independently calibrated rest mass on the lightest charged-lepton branch. [Written field](candidate_laws/interval_echo/openwave/electron_field/README.md) · [Contact-ring records](candidate_laws/interval_echo/openwave/electron_field/lepton_family/theorem_interfaces/interacting_chain/independent/README.md)

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
pip install -r model/requirements-tests.txt
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

| Observable | What the observers read | Candidate | Register row |
| --- | --- | --- | --- |
| A clock carried by a patch | Both owners read a travelling nine-bit elapsed word and a home word on the same signed seam at reunion. The patch derives each handoff interval from its own repair record. | [Interval Echo](candidate_laws/interval_echo/SPEC.md) | finite guided comparison |
| Particle stability | A signed reading from five local queries on a two-tip spectroscopy guide persists at every fresh first read. All 32 simulated launches have projector weight above 0.9999994989 after the numerical allowance. The signed detector clicks with one eighth of that weight. Changing one physical tip during flight lowers every signed return below 0.515 at exposure 37. | [Interval Echo spectroscopy preparation](../openwave/interval_echo/row_04/spectroscopy_release/METHOD_NOTE.md) | `ow-04-particle-stability`, Partial on the submitted pair |
| A clock inside an object | After a pulse the seam record of a cage swings at one line with no drive. An observer reads the line at 4.8990 of its own reference and reads the level difference of the two modes separately at the same value. | scalar channel | `ow-03-de-broglie-clock`, partial on a prepared clock |
| Emission of a light daughter | A lower parent and remote absorber tag the same low mode on powered 23- and 45-dimer guides. A separate signed read detects both charges before absorption. The direct absorber accepts only balanced charge. | [Interval Echo decay guide](candidate_laws/interval_echo/openwave/row_09_neutral/critical_emission/cycle_energy/README.md) | `ow-09-neutrinos`, light tagged mode |
| Neutral massive material | A released polar reversal has native rest cost 39.91317. At time 64 another owner detects the upper dimer in 81.03–81.12% of starts versus 0.47–0.60% on the unreversed host. | [Interval Echo polar rings](experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/README.md) | `ow-19-dark-matter-candidate`, internal host charge |
| Waves of an object | Two measured couplings and independently read wavelengths predict 22 finite-lattice frequencies on three board sizes within 0.000575 clock units. | [Interval Echo prepared strands](candidate_laws/interval_echo/openwave/submission/row_30_confirmation/README.md) | `ow-30-klein-gordon`, one-dimensional dispersion |
| What a reading leaves | Asking "is a pair there?" and asking "where is it?" give the same first answer and different later records. | contrast-memory reference | |

The last column names the row of the [observables register](../docs/registers/OBSERVABLES.json). Each result
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
