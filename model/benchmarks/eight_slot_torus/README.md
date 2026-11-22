# The eight-slot torus benchmark

This benchmark runs a walk on a supplied four-dimensional torus: four periodic coordinates, two directions on each,
eight coherent slots per address. One update reflects in the equal-slot ray, shifts one address and reverses the
slot; two constituents that meet in matched slots step together under a contact coin, and a head-on encounter forms
a pair. The code records what the observers of the torus read and runs controls, changed laws that alter one part
of the torus law, to show which part produces each result. The [specification](SPEC.md) states every part exactly.

The torus is an architecture its makers supply, so this is no candidate for the native repair law: the
[scorecard](../../evaluation/SCORECARD.md) fixes the spherical echosahedral federation as the carrier and counts a supplied
spatial architecture as hardcoding. Its clocks, waves and bound pairs are a phenomenological benchmark that a
candidate on the federation can be compared with.

## Run it

With `model/requirements.txt` installed, from the repository root:

```sh
python -m model.benchmarks.eight_slot_torus laws          # the law and the controls of its experiments
python -m model.benchmarks.eight_slot_torus check         # exactness and regression checks, about ten seconds
python -m model.benchmarks.eight_slot_torus run all       # every experiment as a run file
```

`run all` writes one run file per experiment to `model/benchmarks/eight_slot_torus/runs/`, which stays out of the
repository. `python -m model.benchmarks.eight_slot_torus run continuum --size 48` runs a single experiment.

## What the runs show

| Run | What it shows |
| --- | --- |
| Spacetime · The 3D scalar continuum, seen from an observer | A field pattern released at rest on $32^3$ addresses, evolved with the exact law. Observer A sits at the centre; a reading reaches A at $c = 1/2$, so A's picture at its clock time $t$ shows each address as it was at $t - r/c$, and the record horizon grows until it covers the carrier. The view switches between the law's readings, the Klein–Gordon field it converges to and their difference, and between A's records and the chart description. |
| Clocks · A clock carried by the travelling wave | The moving/rest beat ratio against the speed for launch momenta from long to short wavelengths, with the curve of special relativity. |
| Signals · A fresh message on the free wave | A compact launch at A read by a detector at B, on a log scale: zero probability before the strict front, precursors below $10^{-38}$ before the early deadline, the message with the smooth cone. |
| Clocks · Two clocks on one journey | A beat clock and an optical clock travel the diagonal and return; dials and advance charts for the traveller and its stationary twin, and the reunion readings against the forecast. |
| Clocks · A clock held together by an interaction | The same journey for the bound pair's rail clock. |
| Clocks · The finite-wavelength correction | How far the smooth metric, the finite law at the central wavelength and the finite law over the whole packet lie from the recorded rates of the bound clock. |
| Matter · A pair formed in an encounter | Two constituents meet head-on, the capture forms a pair and sends a marker up the junction's wire, and the pair opens its hop and travels on. |
| Matter · The matter-wave relation | Resting and moving matter packets coloured by phase, and the frequencies of fourteen directions at three periods against $\omega^2 = m^2 + c^2\vert k\vert^2$. |

The [single-aperture grating measurements](formed_clocks/README.md#fringes-after-a-marked-birth)
compare a normally formed pair, a coherently loaded pair, and a free scalar
at an optical clock. The outgoing marker determines which phase comparisons
are available to the material detector.

## Ranging with six local ports

A station launches one constituent in its local equal-slot state. Five
addresses along each spatial axis form a periodic carrier with one private
position. Six stations surround a seventh address that applies a phase of
π to its eight slots at ticks 3 and 15. A sham uses the same schedule.
Each local query swaps the equal-slot mode into an empty latch and records
its occupation. A miss leaves the orthogonal wave on the carrier. These
launches, latches, phase drives, synchronized counters and the classical
collection of the two stations' records are supplied services.
[Source, controls and acquisition](six_station.py)

The receiving station can hold a capture until tick 12 and release it.
Both stations record the first click in a twelve-update window. Given a
receiver click at tick $j$ and a sender click at tick $12+k$, subtracting
the recorded dwell $12-j$ from the elapsed time gives the range $(j+k)/2$
in ticks. A pair's range is the equal-weight mean of the two directed
conditional means. The archive retains all 169 joint outcomes in each
directed menu. Out of $30\cdot2^{34}$ launches per pulse setting,
171,311,326 trips completed with the sham and 109,477,123 with the phase
pulse. A separate two-slot challenge keeps both local letters, either
timeout and every wrong reply. A parked receiver yields zero confirmed
returns and sender clicks after a receiver miss. Disabling its detector
allows the unmeasured wave to propagate through the receiver.
[Complete counts and sampling weights](six_station_records.npz)

| Quantity | Sham | Phase pulse |
| --- | ---: | ---: |
| Mean opposite range in ticks, exact source rounded | 8.86958 | 8.22838 |
| Mean adjacent range in ticks, exact source rounded | 4.80869 | 4.96324 |
| Least worst-pair relative error of one 3D placement | 13.2039% | 7.9312% |
| Simultaneous record bounds on that error, rounded outwards | 12.910% to 13.547% | 7.612% to 8.291% |

The count calculation uses only the recorded times. Conditional on the
number of completions, a confirmed half-return lies in $[1,12]$. The reader
applies Hoeffding bounds to the directed success frequencies and the
conditional means, accounts for rounding the acquisition probabilities,
then averages directions. The joint failure probability is below
$10^{-8}$ under the independent-trial sampling law. The five-point Gram
test alone excludes zero subsets at 5%. A different Euclidean identity
excludes a placement even in higher dimension. For any six points with
squared distance matrix $D^2$ and any coefficients with zero sum,
$x^T D^2 x=-2\Vert\sum_i x_i p_i\Vert^2\le0$. The reader checks this
inequality on the simultaneous distance intervals. One octahedron gives
a common upper bound on distortion for all fifteen pairs.
[Geometry certificate](ranging_geometry.py)

The exact law also gives the optimum. Signed spatial axis permutations
preserve the launch, phase pulse and detectors. Write $a$ for the adjacent
range and $b$ for the opposite range. Two opposite pairs with coefficients
$(1,1,-1,-1)$ require $(1-t)b\le\sqrt2(1+t)a$ for any Euclidean fit of
relative error $t$. Both preparations have $b>\sqrt2 a$. A regular
octahedron with diameter $(1-t)b$ and edge $(1+t)a$ attains equality on
all fifteen pairs, so $t=(b-\sqrt2 a)/(b+\sqrt2 a)$. This tests the
conditional half-return ruler and its finite observation schedule. The
free field's long-wavelength continuum uses a different observable.
[Exact source comparison and record reader](station_ranging.py)

```sh
python -B -O -W error -m model.benchmarks.eight_slot_torus.station_ranging
python -B -O -W error -m model.benchmarks.eight_slot_torus.station_ranging --details
```

The default replay authenticates the retained counts and recomputes each
sampling distribution from the literal torus update and independent integer
evolution. The optional `--seed-replay` also repeats the original NumPy
multinomial draws.

## Checks against the recorded values

The benchmark reproduces the recorded values where it computes the same
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

The [formed-clock experiments](formed_clocks/README.md) connect ordinary
capture to a local observer's elapsed-time reading, with separate continuous
source, phase-reference and transport controls. `python -m model.benchmarks.eight_slot_torus run
loaded-clock` and `python -m model.benchmarks.eight_slot_torus run driven-clock` regenerate the run files of the
stationary receiver and the driven reunion.

The [electron tests](electron_wall/README.md) probe spatial rotations of
the eight-slot torus on both mature rails. A separate fifteen-port walk
supplies the private wall. Its charged excitation has a calculable rest line and
spin and orbital frequencies measured with the same clock.

| File | Content |
| --- | --- |
| [law.py](law.py) | Carrier, coin, routing, contact, capture and opening; the law and its controls |
| [sector.py](sector.py) | Exact private-sector evolution, the field reading, the continuum solution |
| [spectral.py](spectral.py) | Exact plane-wave evolution of prepared packets and bound pairs, envelope tracking, the clock readings |
| [experiments.py](experiments.py) | The continuum run, the journeys, the signal race and the scans |
| [formation.py](formation.py) | The capture vertex with marker wires, evolved branch by branch |
| [matter_wave.py](matter_wave.py) | The matter-wave relation of a formed pair |
| [runs.py](runs.py) | The run-file format |
| [test_simulator.py](test_simulator.py) | Checks of every method against the literal update |
