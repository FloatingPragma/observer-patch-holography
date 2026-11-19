# OPH simulator

This simulator runs the repair laws that are in the search for the native
repair law, records what their observers read, and replays the records in a
web app. It executes one family of laws, reflected routing, with the
parameters that distinguish its members: the equal-slot phase of the coin,
the phase two excitations gain when they meet, and the access rule that says
which operations inhabitants can trigger. The
[specification](SPEC.md) states the law exactly. The
[journal](../docs/research/laws/JOURNAL.md) records why each member has its
standing.

| Member | Standing |
| --- | --- |
| Reflected routing, field-only access | candidate |
| Reflected routing with an encounter phase | undecided |
| Reflected routing with an inhabitant-accessible relay | excluded (control) |
| Equal-slot phase changed by 0.05 rad | disfavoured (control) |
| Equal-slot phase i | control |

## Run it

With `model/requirements.txt` installed, from the repository root:

```sh
python -m simulator laws                      # members and their standing
python -m simulator check                     # exactness and regression checks, a few seconds
python -m simulator run continuum --size 32   # the scalar field against its continuum limit
python -m simulator run journeys              # two clocks on the step-3 routes
python -m simulator app simulator/runs/*.run.json.gz
```

The last command writes `simulator/runs/replay.html`, one offline file with
the runs embedded. Open it in a browser. The app also opens any run file
through its file button or by drag and drop, so a run computed on a larger
machine replays the same way. Run files and the built app stay out of the
repository.

## What the app shows

**The continuum from an observer's place.** The continuum run releases a
field pattern at rest on a carrier of $N^3$ addresses and evolves it with the
exact law. Observer A sits at the centre. A reading travels to A at the signal
speed $c = 1/2$, so A's picture at its clock time $t$ shows each address as
it was at $t - r/c$, and addresses that no record has reached stay grey.
The record horizon grows at $c$ until it covers the whole carrier. The view
switches between the law's own readings, the Klein–Gordon field that the law
converges to, and their difference, and between A's records and the chart
description that all observers agree on once every record is collected. A
plane through A and the time series of three stations' readings, as they
arrive in A's records, accompany the 3D view.

**Two clocks.** The journeys run carries a beat clock and an optical clock
on the axial and diagonal routes of step 3 and compares their reunion
readings with the flat-metric forecast.

## A first test of the candidate

On a carrier of $32^3$ addresses the candidate's field stays within 1.8% of
the Klein–Gordon solution, relative to the largest initial value, over three
chart time units; at $16^3$ and $24^3$ the deviation is 7.3% and 3.2%. The
deviation falls with the square of the resolution, as the continuum theorem
states. On the step-3 routes the exact expectations of the two clocks are
0.71082 and 0.70950 against a forecast of 0.71127 on the axial route, and
0.58166 and 0.58043 against 0.58057 on the diagonal route.

## Files

| File | Content |
| --- | --- |
| [law.py](law.py) | Carrier, coin, routing, encounter term, the members |
| [sector.py](sector.py) | Exact private-sector evolution, the field reading, the continuum solution |
| [spectral.py](spectral.py) | Exact plane-wave evolution of prepared packets, envelope tracking, the four-aperture clock reading |
| [experiments.py](experiments.py) | The continuum run and the two-clock journey |
| [runs.py](runs.py) | The run-file format |
| [app/index.html](app/index.html) | The replay app |
| [test_simulator.py](test_simulator.py) | Checks of every method against the literal update |
