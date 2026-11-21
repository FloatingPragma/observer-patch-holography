# The eight-slot gate as a seam program

This package compiles a supplied eight-slot gate into seam events of law 1 on the
20-cell federation. The gate acts on a 3 × 2 × 2 × 2 index grid with eight slots
per index. One port is a fixed marker, and the first 192 other ports in numerical
order hold the slots. The grid and its four index directions belong to the program;
the carrier supplies only ports and seams. A native fuel source prepares a packet,
and a phase-sensitive receiver reads it at the home seam.

## From native pair fuel to an observer's phase reading

The eight-slot gate can be built from the recovered contrast-memory law. One of two loads marks a reference port. A seam exchange moves the other load. Seven commuting matchings of partial exchanges produce the equal-slot reflection exactly. Full exchanges route the reflected slots. On a 20-cell federation the compiler checks the whole 192-column operator. One routed tick uses 4,302 physical seam events with dependency depth 2,301. The federation has 240 ports and 630 native seam records plus no pair. With the actual finite fuel pulse the largest total variation from the whole-grid forecast is 0.000547 across five waits and six aperture settings. Dropping a matching raises it to 0.4015. Transport alone gives 0.6832. [Program and controls](protocol.json) · [Full-seam results and resources](records/result.json)

The routing pulses preserve the uniform amplitude of a sharp source. At other durations, two native pulses on one cell can take that source to a state with zero uniform amplitude. To prepare a specified packet, two absence reads herald an antisymmetric pair of rails. Exchanges on those rails set each requested amplitude and phase. Reading the unused rails returns the state to the eight-slot encoding. On a sufficiently large connected component the ideal unconditional probability of producing **any specified finite complex packet** by this protocol is exactly one quarter. Every rejected rail has a recorded outcome. The 72-slot four-mode execution used the finite emitter, 4,976 seam events and 74 reader queries on its accepted preparation path. The source accepted 262,399 of 1,048,576 exposures. The largest discrepancy across the complete preparation and subsequent seam records was 0.002721, within the emitter bound of 0.023438. [Proof and error budget](PROOF.md) · [Preparation records](preparation_records/result.json)

The spatial apertures failed their registered phase-control criterion. With phase pulses omitted, the largest total variation between the observed counts and sealed ideal counts was 0.018617, below the emitter bound, even though the conditional target fidelity dropped from 0.9981 to 0.0201. A matched receiver tests this difference in an agreed record. It reverses the successful preparation, including every absence read, before querying the emitter seam. The ideal absolute return probability is one sixteenth of the squared mode overlap. Both source and receiver use the original blank-memory generator at every pulse angle. The table includes preparation and read failures in its denominator. [Forecasts](receiver_records/forecast.json) · [Readings and resource counts](receiver_records/result.json)

| Source phases | Receiver phases | Forecast | Home-seam counts / 1,048,576 | Observed frequency |
| --- | --- | ---: | ---: | ---: |
| Programmed | Programmed | 0.062500 | 65,474 | 0.062441 |
| Omitted | Programmed | 0.001264 | 1,322 | 0.001261 |
| Programmed | Omitted | 0.001264 | 1,338 | 0.001276 |
| Omitted | Omitted | 0.062500 | 65,505 | 0.062470 |

The phase read is an absolute frequency available to the two native port observers. The four-mode run is finite.

## Clocked scalar field and free release

The scalar field result uses a three-dimensional spatial chart, a compact private register and a stationary material clock. A 20-cell execution supplies 144 complete distributions. A second instrument releases the same preparations under free exchange and reads every possible marker partner. At clock count three it detects an excess marker loss of 0.08355 to 0.08832 across all three private orders. The leakage theorem explains why this encoded field requires protection. [Continuum theorem, rest spectrum and free release](CLOCKED_FIELD.md)

## Reproduce

```
python -B -O -W error -m model.experiments.seam_program
python -B -O -W error -m model.experiments.seam_program.preparation
python -B -O -W error -m model.experiments.seam_program.receiver
python -B -O -W error -m model.experiments.seam_program.clocked_run
python -B -O -W error -m model.experiments.seam_program.free_marker
python -B -m pytest -q -W error model/experiments/seam_program
```

Each command checks the source pins, recomputes the sealed forecasts and compares
every opened record with its archive.
