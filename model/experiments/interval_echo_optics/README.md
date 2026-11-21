# Optical comparison of a carried interval hand

At fixed-home reunion, the three-bit Interval Echo hand agrees with the prediction from its home observer’s optical records for three held out-and-back schedules. At eight services per handoff, the three-bit interval counter saturates and the clock runs close to its raw primary rate. A five-bit clock on the same optical itinerary brings its elapsed ratio back into the frozen equivalence band. [Primary protocol](PROTOCOL.md) · [Separate capacity protocol](CAPACITY_PROTOCOL.md) · [Material witnesses and estimates](results.json)

| Out/back services | Moving interval width | Out/back speeds from the home optics, in relay units | Traveller/home ratio | Optical prediction | Difference ± one paired SE | Frozen decision |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| 3 / 5 | 3 bits | 0.33348 / 0.20015 | 0.966543 | 0.965885 | +0.000658 ± 0.000382 | Agrees |
| 5 / 3 | 3 bits | 0.20003 / 0.33341 | 0.966049 | 0.965913 | +0.000136 ± 0.000334 | Agrees |
| 7 / 7 | 3 bits | 0.14273 / 0.14282 | 0.989576 | 0.989755 | −0.000180 ± 0.000268 | Agrees |
| 8 / 8 | 3 bits | 0.12503 / 0.12503 | 0.999759 | 0.992153 | +0.007606 ± 0.000248 | Disagrees |
| 8 / 8 capacity comparison | 5 bits | 0.12505 / 0.12501 | 0.992118 | 0.992153 | −0.000035 ± 0.000261 | Agrees |

Agreement requires the entire residual interval of four paired block standard errors to lie inside ±0.004. Its disjoint complement gives disagreement. Every raw-primary clock control disagrees. In the independent 8/8 capacity sample, the narrow control receives exactly the same primary word as the wide hand. At this interval the narrow rule keeps every primary. The narrow residual is +0.007802 ± 0.000264, and the wide-minus-narrow change is −0.007836 ± 0.000020. The temporal sequence, home clock and optical calibration are common to the two readings in this paired comparison. [Complete joint packet histograms](packets.npz)

## Local hardware and optical events

On the 80- and 320-cell spherical federations the home station is installed at cell zero. On each board its three different physical intercell contacts each carry an outward and a return journey. The geometry and register layout are fixed independently of the clock outcomes. The 23-cell home patch includes its own hand, elapsed register, six nine-bit recording banks and the shared reunion cell. The far station belongs to a distinct cell outside home. The travelling patch initially contains 12 cells with three-bit interval fields, or 14 cells with five-bit fields. It physically moves its signed-contrast source, hand, elapsed word and interval memory. The two owners read the same signed seam word at return. Each material witness stores the actual contact endpoints, custody path, home banks, all local handoff and primary records, and both final readings. [Native material replay](experiment.py) · [Registered operations](../../candidate_laws/interval_echo/clock.py)

For ranging, a prepared double makes a round trip along a reciprocal guide. At the far station an occupancy read triggers its reflection. The home patch keeps emission and receipt times. The trip uses a second prepared response double on a remote reservoir. The far station releases it on physical arrival of the moving marker. Home reads it on the next optical opportunity and keeps that local time. It also keeps trip start, reunion and the published elapsed word. Turning off the far trigger removes the response and leaves the completed handoffs and final hands intact. [Control and independent estimator tests](test_optics.py)

Half the calibration round-trip time gives the one-way light duration under reciprocity. Subtracting that duration from the home interval between trip start and report receipt gives the outward leg time. Reunion gives the remaining leg time. Each speed is the measured light duration divided by its leg time. The predicted proper-time ratio is the sum of the two leg chords divided by the elapsed home time. All quantities in this estimator are averages of retained words. Individual preparations can remain quiet. The nonnegative integer count histogram retains every realized word tuple, with each block's complete enrolment and relative multiplicities.

The four narrow schedules each have six material witnesses. The capacity control has six more. The guide paths and contact endpoints distinguish physical installations. Their handed-off, received and reunited event sequences coincide for each schedule. Each such sequence supplies one statistical cohort of 2,097,152 simulated preparations, analyzed in 64 independent paired blocks. Cohort replication resamples both local fair sources and the registered finite closing operation conditional on that chronology. The six copies contribute hardware checks. Their equality of responses is built into the common service and signal timing. They provide no separate measurement of a rest-frame anisotropy. The two optically compared widths in the 8/8 control share all primary marks in the statistical experiment.

The supplied apparatus selects the route, source preparation, optical sense, adjacent gate sequence, dwell opportunities and readout protection. This tests the clock against finite optical records under that reciprocal service. Physical speed for a freely propagating signal or a moving bound excitation and agreement of independently coupled clock constructions require their own measurements. The unfiltered-primary control and the saturation failure delimit the measured result.

## Reproduction

From the RER repository root the first command verifies all material witnesses, the five statistical cohorts and the retained block histograms. The second checks the receiver trigger, ownership and independent optical inversion.

```sh
python -m model.experiments.interval_echo_optics
python -m pytest -q model/experiments/interval_echo_optics
```

The run accepts `--output DIRECTORY` for independent deposition. The integer histogram arrays use five bits per observed field in ascending order of significance. The JSON lists the seven field names and the block pointer arrays partition each schedule’s packed reads and multiplicities. All readings in these experiments fit their material nine-bit counters and the compact serialization.
