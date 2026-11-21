# Experiments

An experiment is a program that observers run with a candidate law on the fixed carrier: a preparation, the
instruments that read it, and sometimes an extended construction with registers or couplings of its own. Each
experiment declares what it supplies. The results are statements about that construction; they enter the
evaluation of a candidate only through the evidence the candidate cites.

| Experiment | Law | What it asks |
| --- | --- | --- |
| [material_inertia](material_inertia/README.md) | contrast-memory reference | what one measured binding energy predicts about the motion of clusters, and whether their rest energy and inertia share one dictionary |
| [interval_echo_clock](interval_echo_clock/README.md) | Interval Echo | two material hands at a local reunion, with raw and accepted repair-event controls |
| [interval_echo_optics](interval_echo_optics/README.md) | Interval Echo | held journeys along three physical home links per board and the paired interval-capacity control |
| [interval_echo_space](interval_echo_space/README.md) | Interval Echo and its shared contrast sector | three spatial response components transmitted to neighbouring patches and a forecast at untrained stations |
| [interval_echo_encounter](interval_echo_encounter/README.md) | Interval Echo and its shared contrast sector | pair formation, skipped and unread measurements, and changes in the stable bound band |
| [interval_echo_consensus](interval_echo_consensus/README.md) | Interval Echo | strict intercell arbitration, complete two-unit endpoint channels and a local phase reader |
| [contrast_clock](contrast_clock/README.md) | contrast-memory reference | a fixed program of seam events and two readings under one time-independent law; its clock register is a sequencer of the computation |
| [seam_program](seam_program/README.md) | contrast-memory reference | a supplied eight-slot gate on a 3 × 2 × 2 × 2 index grid compiled into seam events on the 20-cell federation, with a fuel source and a phase-sensitive receiver |
| [cage](cage/README.md) | scalar channel | an object held on a ring of four seams: its theorems, its two operators as words on registers, and three experiments read from counts |
| [contact](contact/README.md) | contrast-memory and scalar channel | a coupling of the contrast square of the contrast-memory reference to the scalar occupation with a supplied coefficient; collar and fuel readings, a mobile pair, six-station ranging |
| [mobile_memory](mobile_memory/README.md) | contrast-memory and scalar channel | mobile charges of two signs, the bound pair and a written cage link on one carrier; speed and clock tests |
| [radiative_repair](radiative_repair/README.md) | carrier only | two opposite charges that convert into a pair of neutral quanta, with a law of its own on the ports of the carrier |
| [refinement_spectrum](refinement_spectrum/README.md) | scalar channel hop | a four-child register attached to the scalar hop; rest frequencies and the ratios of successive ones |

Every experiment runs from the repository root with `python -B -O -m model.experiments.<name>` and is tested by
`python -m pytest -q model/experiments/<name>`.
