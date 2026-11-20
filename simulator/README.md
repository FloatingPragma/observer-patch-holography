# OPH simulator

This simulator runs the OPH law, records what its observers read, and
replays the records in a web app. The [law](../model/README.md) page states
the law with what is proven about it, and the [specification](SPEC.md)
states every part exactly. The experiments also run controls, changed laws
that alter one part of the OPH law, to show which part of the law produces
each result.

## Run it

With `model/requirements.txt` installed, from the repository root:

```sh
python -m simulator laws          # the law and the controls of its experiments
python -m simulator check         # exactness and regression checks, about ten seconds
python -m simulator run all       # every experiment, and the app
```

`run all` writes one run file per experiment to `simulator/runs/` and the app
to `simulator/runs/site/`. Open `simulator/runs/site/index.html` in a browser;
the pages open from disk and request nothing outside the folder. Every
experiment and every OpenWave row has its own page, so each has its own
address once the folder is served:

| Address | Page |
| --- | --- |
| `index.html` | The first experiment |
| `experiments/` | Every experiment |
| `experiments/NAME/` | One experiment; `NAME` is its run file's name without the number, as in `experiments/continuum/` |
| `openwave/` | The 31 OpenWave rows |
| `openwave/row-N/` | One row: the problem, an animated diagram and OPH's result |
| `openwave/row-N/replay/` | The replay of a row's run |

The app holds its code once (`app.js`, `app.css`), the list of pages with the
OpenWave rows (`site.js`) and each run once (`data/NAME.js`). It also opens
any run file through its file button or by drag and drop, so a run computed
on a larger machine replays the same way. A page address also takes the view
after `#`, for example `experiments/continuum/#frame=chart&t=1.5&play=0`.
`python -m simulator run continuum --size 48` and
`python -m simulator app RUN ...` run and package single experiments. Run
files and the built app stay out of the repository.

## What the app shows

| Run | What it shows |
| --- | --- |
| Spacetime · The 3D scalar continuum, seen from an observer | A field pattern released at rest on $32^3$ addresses, evolved with the exact law. Observer A sits at the centre; a reading reaches A at $c = 1/2$, so A's picture at its clock time $t$ shows each address as it was at $t - r/c$, and the record horizon grows until it covers the carrier. The view switches between the law's readings, the Klein–Gordon field it converges to and their difference, and between A's records and the chart description. |
| Clocks · A clock carried by the travelling wave | The moving/rest beat ratio against the speed for launch momenta from long to short wavelengths, with the curve of special relativity. |
| Signals · A fresh message on the free wave | A compact launch at A read by a detector at B, on a log scale: zero probability before the strict front, precursors below $10^{-38}$ before the early deadline, the message with the smooth cone. |
| Clocks · Two clocks on one journey | A beat clock and an optical clock travel the diagonal and return; dials and advance charts for the traveller and its stationary twin, and the reunion readings against the forecast. |
| Clocks · A clock held together by an interaction | The same journey for the bound pair's rail clock. |
| Clocks · The finite-wavelength correction | How far the smooth metric, the finite law at the central wavelength and the finite law over the whole packet lie from the recorded rates of the bound clock. |
| Matter · A pair formed in an encounter | Two constituents meet head-on, the capture forms a pair and sends a marker up the junction's wire, and the pair opens its hop and travels on. |
| OpenWave Row 30 · The matter-wave relation | Resting and moving matter packets coloured by phase, and the frequencies of fourteen directions at three periods against $\omega^2 = m^2 + c^2\vert k\vert^2$. |

The OpenWave menu lists the 31 rows of the OpenWave comparison with an
animated diagram of each problem; a row with a passing result on the OPH law
links to its replay.

## Checks against the recorded values

The simulator reproduces the recorded values where it computes the same
instrument. The axial two-clock journey reunites at 1459 home ticks with a
forecast of 0.71127, the beat clock reading 0.71082 and the optical clock
0.70950. The bound pair on the diagonal reunites at 1787 ticks, its rest clock
advances 1.36699 rad and its local home reading is 0.577424; the finite law
at the central wavelength gives 0.685670, 0.576056 and 0.454292 at central
orders 12, 16 and 22. The spectator control reads 1.000000. On $32^3$
addresses the scalar field stays within 1.8% of the Klein–Gordon solution
over three chart time units, and the deviation falls with the square of the
resolution. The formed matter's frequencies match the recorded rest line
0.098175 and the long-wave coefficient 0.24919629.

## Transverse pulse clocks

A circulating free wave supplies a clock that can be read by one station. Five
of six long-wave winding peaks meet the five-percent timing tolerance. The
coarsest reversed wave misses it, and at size 384 the angled receiver accepts
about one tagged click per 1.75 billion preparations at its peak. The frozen
replay checks every enrolled count, the exact free-wave tails and the costs
of an equipped relay. [Apparatus and results](transverse_clock/README.md)

## Files

| File | Content |
| --- | --- |
| [law.py](law.py) | Carrier, coin, routing, contact, capture and opening; the law and its controls |
| [sector.py](sector.py) | Exact private-sector evolution, the field reading, the continuum solution |
| [spectral.py](spectral.py) | Exact plane-wave evolution of prepared packets and bound pairs, envelope tracking, the clock readings |
| [experiments.py](experiments.py) | The continuum run, the journeys, the signal race and the scans |
| [formation.py](formation.py) | The capture vertex with marker wires, evolved branch by branch |
| [row30.py](row30.py) | The matter-wave relation of a formed pair |
| [runs.py](runs.py) | The run-file format |
| [webapp.py](webapp.py) | Writes the replay app, one page per experiment and OpenWave row |
| [app/index.html](app/index.html), [app/app.js](app/app.js), [app/app.css](app/app.css) | The replay app: the page template, its code and its style |
| [test_simulator.py](test_simulator.py) | Checks of every method against the literal update |
