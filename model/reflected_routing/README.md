# What a clock and a message measure on reflected routing

Bernhard Müller's reflected-routing experiment gives a long-wave scalar signal speed of one half and a moving clock made from a beat between two private modes. Can a different clock read the same slowdown? Can observers use the gate to deliver a message faster than that signal? Here a circulating transverse pulse supplies the second clock. Stations that measure and regenerate a pulse deliver a fresh authenticated message ahead of the scalar prediction on three of four routes. The exact passive propagation rate quantifies the cost of doing so without intermediate stations. [Clock reunion](https://github.com/muellerberndt/openwave/blob/aa71bbcb335503eef35cbe14da5b0ad4b7dd083e/simulation/reports/CLOCK_REUNION.md) · [Reflected-routing clock and field](https://github.com/muellerberndt/openwave/blob/260de06352eac5425d78fbbae636625544d439ab/simulation/reports/relativistic-observer-records/report.md)

This experimental member uses the eight-slot wave gate from the [pinned source](https://github.com/muellerberndt/openwave/blob/260de06352eac5425d78fbbae636625544d439ab/simulation/experiments/observer_geometry/reunion/registered/ingredient_sources/radar_shell/replay.py), with supplied sources, a shared scheduler, calibrated local clocks, readers and relay stations. Three periodic coordinates label the stations' spatial dials and one is private. The spatial comparison uses Euclidean distances reconstructed from reciprocal addressed messages on the three axes. Each clock or communication result below is a readout of this member. A transfer to polynomial-sv's twelve-port carrier and repair law remains open. No measured particle property enters the construction.

## A clock with a readable winding

A packet in private order zero travels around a public transverse circuit. Independent preparations measure its longitudinal speed, a reference signal speed and the rest rate of that circuit. These fix the relativistic forecast before observers open the moving transverse counts. Each position distribution includes all eight slots. At each of four sizes from 128 through 1024, all three long-wave clock ratios agree within the registered 0.035 tolerance. At size 1024 the two slower clocks read 0.80066 and the faster one 0.60099. The short-wave clock reads 0.51022 against 0.73270 and misses at every size. Erasing the transverse phase makes the clock unreadable. The exact finite-gate prediction reproduces every coherent bin. [Forecasts](records/forecast.json) · [Readouts](records/result.json)

Can **one local station** detect a completed winding? A second registration uses three fresh grids and three reversible bits that count signed crossings of the transverse cut. A single receiver counts all slots at one address in all four coordinates. The rest clock and the longitudinal speed select that receiver and its observation times before the moving counts are opened. The table gives the peak time for exactly one positive winding, followed by the relativistic forecast. [Local apparatus](local_clock/protocol.json) · [Sealed predictions](local_clock/records/forecast.json) · [Complete local results](local_clock/records/result.json)

| Grid size | Angled | Reversed | Short wave |
| --- | ---: | ---: | ---: |
| 96 | 245 / 253.0 | 270 / 288.0 | 343 / 313.8 |
| 192 | 453 / 453.2 | 502 / 501.9 | 641 / 569.7 |
| 384 | 872 / 872.5 | 963 / 963.3 | 1294 / 1119.1 |

Five of the six long-wave cases meet the 5% tolerance. The coarse reversed clock and all short-wave clocks miss it. The independent finite-gate solution predicts **every count** in all nine experiments, including other windings and misses. The local price of these records matters. At size 384 only 815,120 of the $2^{16}384^4$ prepared exposures at the angled peak reach the receiver with winding one, a fraction of $5.72\times10^{-10}$. Counts are deterministic midpoint quantizations of terminal probabilities. Each wait uses a fresh ensemble. The unwritten-tag control restores interference between windings. This experiment tests a prepared wave clock with a specified aperture, bandwidth and time grid. A clock made of interacting matter on the same law remains open.

## Which signals respect the cone?

For passive free propagation, at every fixed axial speed $1/2<u<1$ we derive the exact exponential decay rate of the probability of a click beyond $ut$. The upper bound covers every polarization, a spectator quantum reference and all directions. A suitable axial preparation attains the rate. Exact integer recurrence gives finite detector counts up to 1024 updates. At $u=3/4$ the finite decay rate is 0.82553 per update and its limit is 0.81753. Independent passive copies growing polynomially with distance still have exponentially small early detection probability. [Derivation and resource bounds](PROOF.md)

An equipped observer can perform a different experiment. Every third axial step has a detector with parallel fresh rails. Detection of a complete frame enables a fresh local preparation after one service tick. The sender chooses a nonce after forecasts are sealed. Only the remote endpoint produces the authenticated reply. The following elapsed times come from complete round trips. Counts include the outward and returning legs, in both directions, with 4096 occupied preparations per transmitted bit at each hop. [Protocol and controls](protocol.json) · [Complete channel counts](records/channel.npz)

| Displacement in the calibrated chart | Relay time | Scalar round-trip forecast | Early / returned / attempted | Occupied preparations per completed trip |
| --- | ---: | ---: | ---: | ---: |
| (24, 0, 0) | 63 | 96.19 | 12 / 12 / 12 | 8,388,608 |
| (18, 12, 0) | 79 | 86.71 | 12 / 12 / 12 | 10,485,760 |
| (−12, 0, 15) | 71 | 76.99 | 12 / 12 / 12 | 9,437,184 |
| (6, 6, 6) | 47 | 41.65 | 0 / 12 / 12 | 6,291,456 |

The addressed courier beats the deadline on every route. The passive, changed-coin, blocked, silent, tampered and stale-request controls produce zero authenticated returns. The passive experiment targets the earliest arrival. The proof gives a separate bound over the whole interval ahead of the scalar cone. The faster relay uses finite fresh banks of parallel preparations and intermediate writes. Its speed depends on the available service latency. These results distinguish a long-wave field cone from a universal causal limit for this equipped observer. They establish no laboratory signal speed, universal proper time, particle mass or OpenWave row.

## Offline reproduction

Install the [model requirements](../requirements.txt) and run from the repository root. The fast tests check a whole four-coordinate upstream gate, literal cut crossing, both flat bands, exact finite tails and rejection of corrupt records. The two replays recompute calibration, forecast, controls and acquisition against the stored counts. A single BLAS thread is sufficient.

```sh
python -B -m pytest -q -W error model/reflected_routing
python -B -O -W error -m model.reflected_routing
python -B -O -W error -m model.reflected_routing.local_clock
```

Both programs accept `--stage seal` for the forecast alone. To repeat a prospective run use a fresh path outside the package with `--record PATH --stage seal`, then `--record PATH --stage open`. The communication run draws new sender entropy at opening. The stored replay uses its published entropy. [Source hashes](sources.json) · [Local source hashes](local_clock/sources.json)
