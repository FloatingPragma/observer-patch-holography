# A-contrast-memory-spherical-v1

A-contrast-memory-spherical-v1 is the candidate executable model of OPH. Its law is the contrast-memory law,
a candidate for the native repair law of the theory: the simple local rule
by which neighbouring observers settle their disagreements. Each port of a twelve-port cell (the echosahedron)
holds a load of 0, 1 or 2. Each seam between two ports keeps a small private memory of how much its two ends
disagree. The memory writes the disagreement, gives it a phase, and erases it again, and what is left on the
loads is an exchange of loads plus an energy for disagreeing.

The state of the simulation does not contain the physical Universe or its spacetime. Ports, seams, cells and
federations are the hardware of the computation, and no physical space. Space, time, fields and particles are
effective descriptions that the observers inside the computation use for the records they agree on, so every
test of the model compares distributions of observer records, never state encodings.

## The law in two formulas

On one seam with loads `(a, b)` the private memory reduces exactly to the seam term

$$h_e = \tfrac12 (I - \mathrm{SWAP}_e) + \kappa\,(n_a - n_b)^2 ,$$

and on the whole carrier, with a degree term that makes every port count six seams,

$$H_\kappa = \tfrac12 \sum_e (I - \mathrm{SWAP}_e) + \kappa\Big(6 \sum_v n_v^2 - 2 \sum_e n_a n_b\Big).$$

The coupling `kappa` is the agreement energy in units of the exchange. A seam whose two ends both hold a load
lowers the energy by `2 kappa`, and that is enough to bind loads together. On loads 0 and 1 the seam term is a
ferromagnetic XXZ interaction, so the binding is of the droplet kind known from that literature (Nachtergaele and
Starr, arXiv:math-ph/0009002; Fischbacher and Stolz, arXiv:1712.10276). What the model adds is the memory that
realizes it, the constants on this carrier, and the tests on records.

## What comes out

Two loads written on one seam by a local emitter stay together. The port observers of every seam report one
agreed record, "both ends hold one load", and the pair is read on some seam with probability at least one minus a
proven bound, at every reading, on one cell and on federations of 20, 80 and 320 cells. It is read on other seams
than the first, and nothing holds it but the law. Remove the agreement energy
and the pair falls apart (the largest "no pair" probability goes from 0.0096 to 0.80 on 20 cells). At `kappa = 8`:

| Carrier | Bound with the constants of the carrier | Largest P(no pair) | Largest P(pair on another seam) |
| --- | ---: | ---: | ---: |
| one cell | 0.1488 | 0.00957 | 0.986 |
| 20 cells | 0.1830 | 0.00962 | 0.977 |
| 80 cells | 0.1828 | 0.00964 | 0.977 |
| 320 cells | 0.1830 | 0.00963 | 0.977 |

The bounds of the table use the measured constants of each carrier and include the emitter. Behind them stands a
theorem that needs no constant of the carrier: on every federation of complete cells and for `kappa > 3` the
probability of "no pair" stays below `40 / (2 kappa - 6)^2`, whatever the number of cells ([PROOF.md](PROOF.md)).
Three loads on a face, and bigger clusters, bind as well. Run as local events with no common clock, the pair holds
when every observer serves each of its seams once per round, in any order.

Two more results concern what observers do with the pair and with the memory of a seam. Asking "is a pair
there?" and asking "where is it?" give the same first answer and different later records: the pair is one object
spread over seams, and a reading that keeps its seam places it ([presence.py](presence.py)). And a seam memory
that holds a record while its ports are empty stays that way under every event: it holds no load, its observers
read its level from their own records, and a load that reaches it binds to it ([collar.py](collar.py)). All
numbers are in [RESULTS.json](RESULTS.json) and [presence_results.json](presence_results.json).

## Run it

From the repository root, with NumPy, SciPy and SymPy:

```sh
python -m model.candidate_laws.contrast_memory_spherical_v1.verify            # headline subset, under a minute
python -m model.candidate_laws.contrast_memory_spherical_v1.verify --full     # everything, 80 and 320 cells included (an hour or more)
python -B -O -m model.candidate_laws.contrast_memory_spherical_v1.presence     # replay the two-reading public laws
python -m pytest -q model/candidate_laws/contrast_memory_spherical_v1          # the tests of this candidate
```

`verify` recomputes and compares with `RESULTS.json` within tolerances; `--output PATH` writes a fresh result
file and `--verify PATH` compares another one. The last line it prints is a JSON line with `"verified": true`.

## The modules

| Module | What it computes |
| --- | --- |
| [law.py](law.py) | the local law on the 99 states of one seam and its exact SymPy identities for symbolic `kappa` |
| [carrier.py](../../executable_model/carrier.py) | one cell, the federations of 20, 80 and 320 cells ([carriers/](../../executable_model/carriers/)), and the hypotheses of the theorem |
| [sectors.py](sectors.py) | the law on sectors of fixed load number, the changed laws, the emitter, the reader's outcome law, seam events |
| [theorem.py](theorem.py) | the constants of the binding theorem, its operator inequalities, the emitter bounds |
| [stability.py](stability.py) | particle stability: the pair and the face object on records, the controls, single readings, layered pulses |
| [reader.py](reader.py) | the eleven-bit seam reader and its finite readout with the law acting (161,792 states) |
| [presence.py](presence.py), [proof](PROOF.md#8-pair-presence-and-the-state-after-a-reading) | Two ways to record that a pair is present give different later observer records. The finite-duration difference exceeds 0.0817 on one public event. |
| [events.py](events.py) | the law as local events with no common clock, and the member with committed contrast |
| [catalogue.py](catalogue.py) | the bound objects of one cell, their bindings, the pair on the federation |
| [port_response.py](port_response.py) | the port response for one load and for the bound pair |
| [collar.py](collar.py) | the unrecovered seam (the papers say collar): its kept label, its levels, a load beside it, its reading with four phase settings |
| [extras.py](extras.py) | one-load spectra, three and four loads, the pair band, average adjacency, the branch on a ring |
| [verify.py](verify.py), [test_model.py](test_model.py) | recomputation, comparison and tests |

The full specification is [SPEC.md](SPEC.md); the binding theorem with its proof is [PROOF.md](PROOF.md).

## Supplied

- the coupling `kappa` (tested at 0 to 128; headline values 8 and 128);
- the coherent generator on the private seam memory;
- the reference registers of the degree term on ports with five seams;
- the emitter with its fuel, coupling and duration;
- the evolution parameter of the runs with a continuous parameter (as local events the pair needs each seam to
  act once per round, or small event angles);
- the routing of the external seams of the federation;
- the reader with its stage durations, and the tree of the presence reading;
- the unrecovered seam as a preparation, with its fuel and converter.

## Work in progress

Several parts of the theory are work in progress on this model: a selection of `kappa`, a second sign of load,
a record of currents that later events read, the class of an object as a record, the twelve port currents from
the law's own services on the bound pair, one limiting speed for all objects (on a chain the pair's branch is
`E(K) = Delta - cos^2(K/2)/Delta` with `Delta = 1 + 2 kappa`), a refinement map that preserves seams, an observer
clock normalized to proper time, and the regional response that the Einstein theorem of the theory consumes.
The model is no solver for general constraint problems and claims no polynomial bound for one.
