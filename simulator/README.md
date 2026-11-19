# OPH simulator

This simulator runs the repair laws that are in the search for the native
repair law, records what their observers read, and replays the records in a
web app. It executes the reflected-routing family and its two-constituent
extensions. The [specification](SPEC.md) states every law exactly, and the
[journal](../docs/research/laws/JOURNAL.md) records why each member has its
standing.

| Member | Standing |
| --- | --- |
| Reflected routing, field-only access | candidate: steps 1 to 3 |
| Compact pair with motion-coupled rails | candidate as a prepared clock (step 4); cannot form from separated constituents |
| Compact pair with a capture vertex and outgoing marker | forms a pair (step 5); the pair leaves in both directions and carries no single-traveller clock |
| Reflected routing with an encounter phase | undecided: equal to the candidate for one excitation |
| Compact pair with spectator rails | excluded (control): its clock counts services |
| Reflected routing with an inhabitant-accessible relay | excluded (control) |
| Equal-slot phase changed by 0.05 rad | disfavoured (control) |
| Equal-slot phase i | control |

## Run it

With `model/requirements.txt` installed, from the repository root:

```sh
python -m simulator laws          # members and their standing
python -m simulator check         # exactness and regression checks, about ten seconds
python -m simulator run all       # every successful experiment of steps 1 to 5, and the app
```

`run all` writes one run file per experiment to `simulator/runs/` and the app
`simulator/runs/replay.html`, one offline file with the runs embedded. Open it
in a browser. The app also opens any run file through its file button or by
drag and drop, so a run computed on a larger machine replays the same way.
`python -m simulator run continuum --size 48` and
`python -m simulator app RUN ...` run and package single experiments. Run
files and the built app stay out of the repository.

## What the app shows

Each successful experiment of the journal is one run.

| Run | What it shows |
| --- | --- |
| Step 1 · The 3D scalar continuum, seen from an observer | A field pattern released at rest on $32^3$ addresses, evolved with the exact law. Observer A sits at the centre; a reading reaches A at $c = 1/2$, so A's picture at its clock time $t$ shows each address as it was at $t - r/c$, and the record horizon grows until it covers the carrier. The view switches between the law's readings, the Klein–Gordon field it converges to and their difference, and between A's records and the chart description. |
| Step 1 · A clock carried by the travelling wave | The moving/rest beat ratio against the speed for launch momenta from long to short wavelengths, with the curve of special relativity. |
| Step 2 · A fresh message on the free wave | A compact launch at A read by a detector at B, on a log scale: zero probability before the strict front, precursors below $10^{-38}$ before the early deadline, the message with the smooth cone. |
| Step 3 · Two clocks on one journey | A beat clock and an optical clock travel the diagonal and return; dials and advance charts for the traveller and its stationary twin, and the reunion readings against the forecast. |
| Step 4 · A clock held together by an interaction | The same journey for the bound pair's rail clock. |
| Step 4 · The finite-wavelength correction | How far the smooth metric, the finite law at the central wavelength and the finite law over the whole packet lie from the rates recorded in step 4. |
| Step 5 · A pair formed in an encounter | Two constituents meet head-on, the vertex forms a pair and emits a marker up the junction's wire, and the formed pair leaves its birthplace in both directions. |

## Checks against the research records

The port reproduces the published records where it computes the same
instrument. The axial two-clock journey reunites at 1459 home ticks with a
forecast of 0.71127, the beat clock reading 0.71082 and the optical clock
0.70950. The bound pair on the diagonal reunites at 1787 ticks, its rest clock
advances 1.36699 rad and its local home reading is 0.577424; the finite law
at the central wavelength gives 0.685670, 0.576056 and 0.454292 at the three
step-4 orders. The spectator rails read 1.000000. On $32^3$ addresses the
scalar field stays within 1.8% of the Klein–Gordon solution over three chart
time units, and the deviation falls with the square of the resolution.

## Files

| File | Content |
| --- | --- |
| [law.py](law.py) | Carrier, coin, routing, encounter term, contact law, the members |
| [sector.py](sector.py) | Exact private-sector evolution, the field reading, the continuum solution |
| [spectral.py](spectral.py) | Exact plane-wave evolution of prepared packets and bound pairs, envelope tracking, the clock readings |
| [experiments.py](experiments.py) | The continuum run, the journeys, the signal race and the scans |
| [formation.py](formation.py) | The capture vertex with marker wires, evolved branch by branch |
| [runs.py](runs.py) | The run-file format |
| [app/index.html](app/index.html) | The replay app |
| [test_simulator.py](test_simulator.py) | Checks of every method against the literal update |
