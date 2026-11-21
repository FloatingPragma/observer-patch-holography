# A free-wave clock at a local station

A pulse circulating around a transverse ring gives a local station a clock reading. Its observable is the time at which the station receives the most pulses carrying one positive winding. Five of six long-wave peaks fall within five percent of the forecast from the stationary twin and the longitudinal speed. A narrow aperture makes these clicks rare. The assembled position readings from a whole ring give a different, distributed clock. [Local readings](local_clock/records/result.json) · [Distributed readings](records/result.json)

## Preparation, update and readout

Both experiments prepare one free constituent of the eight-slot torus benchmark. The five-by-five Fourier amplitudes have binomial weights, a common phase convention and private order zero. The replay uses the simulator's equal-slot coin, slot partners and signed axes. The replay compares all eight columns with the simulator's plane-wave update at 93,576 numerical wave-vector inputs. The maximum measured difference is zero. Translation of any point source supplies the remaining real-space columns. Coin and route act once per update. The pulse preparation, chart calibration, winding register and detector are supplied apparatus. [Executable specification](../SPEC.md) · [Source and decision rules](freeze.json) · [Equivalence checks](equivalence.py)

The distributed detector counts every position on one axis, summing all slots and the other coordinates. Independent preparations supply each wait, the reference speed, the longitudinal speed and the stationary transverse rate. The moving pulse's circular position moment gives its transverse rate. All twelve long-wave ratios on rings of size 128, 256, 512 and 1024 meet the absolute tolerance 0.035. The four short-wave ratios miss it, and translating the source uniformly in the transverse direction erases the hand on all four rings. [Distributed protocol](protocol.json) · [Forecasts](records/forecast.json) · [Moving counts](records/moving.npz)

The local monitor starts a three-bit register at zero on every amplitude. Use transverse addresses from −L/2 to L/2 − 1. After the coin, a forward hop from L/2 − 1 to −L/2 increments the register modulo eight. Its reverse decrements it. Every other hop preserves it. This coherent permutation entangles the register with the walk without measuring a crossing. Routing on a ring of length 8L implements the joint update. At the receiver, eight disjoint addresses of the lift resolve the eight register values. The wait is shorter than 6L, which excludes an aliased winding of −7 or +9 at the receiver. Suppressing the write gives a separate control in which paths with different windings interfere. [Monitor and aperture](local_clock/source.py) · [Literal cut test](test_records.py)

The receiver sums all eight slots at one address in four coordinates. Each registered wait uses a fresh ensemble of 65,536 L⁴ preparations. Midpoint quantization of the terminal probabilities gives eight winding counts and one miss count. The unwritten control has one hit bin and one miss bin. The sampled time of the largest tag-one count is compared with the prediction. An endpoint peak, a tie or fewer than 1,000 counts makes the peak unresolved. The stationary calibration fixes the station and the waiting times before the moving read. [Local protocol](local_clock/protocol.json) · [Raw detector bins](local_clock/records/opened.npz)

| Ring size | Angled, read / forecast | Reversed, read / forecast | Short wave, read / forecast |
| ---: | ---: | ---: | ---: |
| 96 | 245 / 253.0 | 270 / 288.0 | 343 / 313.8 |
| 192 | 453 / 453.2 | 502 / 501.9 | 641 / 569.7 |
| 384 | 872 / 872.5 | 963 / 963.3 | 1,294 / 1,119.1 |

The size-96 reversal and every short-wave peak miss the five-percent tolerance. At the size-384 angled peak, 815,120 preparations yield a tag-one click, a fraction of 5.72027 × 10⁻¹⁰. The finite-walk forecast reproduces every tagged bin and miss at each of the nine moving-panel schedules. These are finite readouts of a prepared free wave. [Local forecasts](local_clock/records/forecast.json) · [Complete readings](local_clock/records/result.json)

## Passive propagation

An equal-slot source occupies one axial address and is uniform on the other coordinates. A terminal detector resolves the equal-slot polarization on the half-line at and beyond a specified speed times the wait. The real-space replay carries arbitrary-precision integer amplitudes through the pinned route on a strip of length 2,049. Its scaled coin is four times the simulator's equal-slot reflection. At each of six waits from 32 through 1,024 updates and each of five speeds, the squared coefficients reproduce the exact archived numerator. The replay also computes the complementary polarization on each side of each cut, so the four integer weights sum to the common denominator. The all-polarization exponential bound holds before periodic images enter the experiment. At axial speed 3/4 the finite decay rate at 1,024 updates is 0.82553 per update, with limit 0.81753. [Derivation](PROOF.md) · [Finite numerators](records/result.json)

## Equipped ranging and relay controls

Each chart range uses a supplied addressed courier that moves a payload one neighbor per service. Three exchanges calibrate each station's affine timestamps. The scalar round-trip deadline uses the Euclidean norm of the three reciprocal ranges and the independently measured reference-wave speed. Returning a station's terminal report uses the courier as well. [Ranging and report delivery](equipped/chart.py)

The authenticated relay places equipped receivers every three axial addresses. A frame has 128 bits. Each occupied bit rail launches 4,096 fresh copies per segment. Detection of a complete frame enables local regeneration, with one service tick for each regeneration or remote reply. Initial launch and final read are boundary services. For m segments outward and m return, the transaction costs twice the total outward hop distance plus 2m − 1 service ticks. Times include both legs. The sender and recipient share an enrolled key, and the replay uses the published trial entropy. Each addressed courier in the ranging experiment transports its payload directly. The regenerating wave relay transmits its bits by detector counts at successive stations. [Equipped controls](equipped/race.py) · [Per-bit detector counts](records/channel.npz)

| Displacement | Relay ticks | Scalar deadline | Early / returned / attempted | Fresh occupied preparations per complete relay trip |
| --- | ---: | ---: | ---: | ---: |
| (24, 0, 0) | 63 | 96.19 | 12 / 12 / 12 | 8,388,608 |
| (18, 12, 0) | 79 | 86.71 | 12 / 12 / 12 | 10,485,760 |
| (−12, 0, 15) | 71 | 76.99 | 12 / 12 / 12 | 9,437,184 |
| (6, 6, 6) | 47 | 41.65 | 0 / 12 / 12 | 6,291,456 |

All 12 diagonal transactions return late. The passive, silent, blocked, altered-coin, tampered and stale-request controls each give zero authenticated returns. The addressed courier returns ahead of the deadline on all four routes. A faster equipped relay depends on its bank of fresh preparations and its service delay. [Per-transaction outcomes and resource totals](records/result.json)

## Reproduction and source provenance

From the repository root, install the simulator requirements and run the frozen replay. The second command prints all four exact passive detector weights for every registered threshold.

```sh
python -B -m pytest -q -W error model/benchmarks/eight_slot_torus/transverse_clock
python -B -O -W error -m simulator.transverse_clock --passive-outcomes
```

The replay checks source and fixture hashes, all opened arrays, forecasts and decisions. The input freeze identifies the source of the original protocols, counter arrays and entropy, and the archived protocols keep their original construction labels for provenance. [Frozen inputs](freeze.json)
