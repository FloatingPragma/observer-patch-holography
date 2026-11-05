# A measured execution attached to a host clock

The executable native observer has been run against a real monotonic host
clock. One local clock calibrates the first block; a different local clock
predicts the durations of three untouched blocks. Their measured relative
errors are **10.74%, 5.68% and 12.27%**, within the declared **25% engineering
tolerance**. The tolerance is fixed in the capture policy; it is not a
statistical confidence interval or a derived physical law.

The [raw records](host_clock_records.json) retain two actual executions of
264,532 native attempts and 40,000 movements each. The second adds host-only
arithmetic while leaving every native choice and current identical. All four
of its duration comparisons reject the original calibration. These are
measured results for an engineered software apparatus. They do not establish
an identification of native event time with physical proper time, modular
flow, an atomic clock or the dynamics of our Universe.

## What is measured

Each execution calls the existing `native()` rule on the original twelve-port,
two-load binary carrier. A seeded software service chooses one of thirty seams
and a fair-coin address. The source starts at a declared nonstationary basis
configuration and continues through four blocks of 10,000 moving currents,
without resetting at a boundary. Every zero-current attempt remains recorded.
The record tape is finite, and an attempt cap bounds the acquisition.

A local clock ticks whenever a recorded current changes its port's occupation.
The acquisition samples `time.perf_counter_ns()` at the initial endpoint and
after each completed block. It retains elapsed readings relative to that
initial endpoint, not absolute wall dates or machine identifiers. Python
[documents this counter](https://docs.python.org/3/library/time.html#time.perf_counter)
as an elapsed-duration clock; its integer-nanosecond API avoids floating
conversion of the readings. The retained platform reports `mach_absolute_time()`
and a resolution of approximately 41.7 ns. The reported resolution is **not an
accuracy or traceability bound**. System scheduling, frequency control,
instrument overhead and service latency have no certified bounds here.

Native execution, pseudorandom choice generation, record collection and the
declared workload all occur inside the timed intervals. Import and source
initialization occur outside them. Endpoint sampling has its own overhead;
the capture does not claim to measure the isolated native function's latency.
The repeated choices come from `random.Random(1729)`, not a physically
authenticated source of independent random bits. Exact replay verifies this
seeded service; it does not prove the physical IID scheduler premise.

## Calibration and held-out measurements

Port 2 records 1,767 ticks in the first baseline block's measured 31.1850 ms.
This fixes the single conversion

```text
host duration per local tick = 547106/31 ns.
```

Port 3 and all later blocks are excluded from the fit. For each validation
block, its count multiplied by that fixed conversion predicts the duration.
The [analysis receipt](host_clock_results.json) retains measured durations,
exact rational predictions and errors, and each individual verdict. A failed
held-out comparison remains a rejection; the analyzer never changes the
calibration or tolerance to make it pass.

The [single-carrier count theorem](COMPARISON.md) separately gives an exact
finite rate-error envelope, independent of stationary preparation. For the
four chosen count comparisons, a deviation tolerance `1/12` from the
expected count rate `1/6` has a combined false-rejection bound below `1/20`
under the ideal independent uniform source. The observed software counts lie
inside that envelope. This source-law check does not turn the seeded sample
into certified random experimental data, nor turn the 25% host-duration gate
into a confidence interval. These two tests have different meanings.

## A measured control that keeps the observer records unchanged

The workload control repeats exactly the same seeded source and continues the
same native state sequence. After every native call it performs 64 iterations
of a specified integer recurrence that cannot affect the native state or its
random choices. The raw tape retains the final workload value as well.

Every recorded native current, local clock tick and block boundary in repair
coordinates agrees between the two executions. The measured duration ratios
of control to baseline are approximately **8.08, 8.14, 8.16 and 8.13**.
The baseline local-clock calibration fails the duration test in every control
block. No uniformity of the added host delay is assumed; the observed variation
is retained.

This is a concrete operational boundary: an unchanged native record history
can run at different rates against its implementing host clock. A hardware
and service contract can establish a useful measured correspondence over a
tested range, but the native record law alone does not select it. This result
does not impose the host clock as a fundamental clock within OPH.

## Replay and scope

```sh
python -m model.native_clock.host_clock --check model/native_clock/host_clock_results.json
python -m pytest -q -W error model/native_clock/test_host_clock.py
```

The replay decodes every packed byte into a seam, coin and observed current,
executes the original rule again, verifies the seeded choice sequence, checks
the workload computation and every block endpoint, then reconstructs the
calibration and all held-out verdicts. Raw bytes encode
`3*(2*seam_index+coin)+current+1`; hexadecimal storage keeps the complete trace
near one megabyte. The raw capture and the analysis both bind their public
source files. Neither file contains user names, machine names or private paths.

`--check` validates the retained data and analysis; it does not generate the
same physical timestamps again. `--capture PATH` performs a new measurement
and refuses to overwrite an existing capture. A new capture can fail the
fixed tolerance and must retain that outcome. Timestamp files are not signed
hardware attestations: replay detects source-record inconsistency but cannot
authenticate fabricated timestamps. The receipt therefore uses
`replay_integrity_passed` separately from each timing verdict.

The achieved attachment is between executed native records and an independently
sampled host counter, with a measured held-out check and a rate-change control.
A traceably calibrated laboratory instrument, identified detector/error model,
modular-phase dynamics and an empirical identification with observed physics
remain separate obligations. This experiment supplies none of those by
renaming host elapsed time.
