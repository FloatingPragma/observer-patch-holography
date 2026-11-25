# How far can a clocked receiver decode one bit?

On the intact 20-cell screen, the closest sender cleared the specified 52% accuracy at the receiver's first tick. The passive receiver missed the forecast at each longer range. Its address-only decisions gave 132,275 correct answers in 262,144 trials at the middle sender's forecast of two ticks, and 131,072 at the far sender's forecast of three ticks. A single addressed relay on the middle route raised the first-tick result to 138,291 correct. [Counts and control arms](records.json)

| Sender port | Addressed round trip, receiver ticks | First-tick passive accuracy | Second-tick passive accuracy | Third-tick passive accuracy | Marked passive round trips at the three deadlines |
| --- | ---: | ---: | ---: | ---: | ---: |
| 82 | 2 | 52.441% | 55.098% | 55.435% | 70, 389, 619 |
| 72 | 4 | 50.000% | 50.459% | 50.173% | 0, 0, 0 |
| 182 | 6 | 50.000% | 50.000% | 50.000% | 0, 0, 0 |

Each accuracy uses 262,144 simulated launches, including quiet outcomes and timeouts. Each passive round-trip count uses 131,072 separate launches. The simultaneous 99% half-width for trained decoders is 0.42736 percentage points (rounded up). The closest first-tick lower limit is 52.014%. Both longer routes have upper limits below 52% at every deadline through their forecast. Even allowing any of the 8,192 possible fixed address decision tables with the same timeout default, a separate simultaneous 99% upper bound at the far sender is below 50.63% at each tick. This last analysis uses the evaluation data, bounding all tables together to account for selection. A zero passive round-trip count retains a nonzero sampling uncertainty. [Count and uncertainty calculation](check.py)

Ports 82 and 72 share cell 6. Port 182 lies in cell 15. The route includes an internal cell-6 link, so the three addressed distances characterize service between port addresses. They provide two distinct sender cells.

## What counts as a distance or a deadline?

The receiver holds a double load. An echo and its inverse on its two local clock ports register a tick when the seam's total load is two. A single carrying the bit and the clock can interact on the same 20-cell screen. The sender chooses its balanced bit after setup and uses a sender-local exchange to prepare one of two sharp releases. The receiver reads the single's local port, or quiet, when its first, second or third positive tick occurs. A failure to reach the requested tick within four interrogation opportunities leaves its default answer in force. Timeouts count in the success rate. The reader receives neither the bit nor the simulator's free-flight exposure. [Established generator and clock instrument](../contact_deadline/README.md) · [Complete design](PROTOCOL.md)

Independent calibration launches transport a single out from port 1, mark its arrival at the selected station and return it, while the receiver reads the same clock seam. One addressed native exchange per clock interrogation moves the single along the route. The observed round-trip words have two, four and six ticks. Removing the last outward contact prevents the remote mark at all three stations. These are exact outcomes of the addressed service with the hand on its protected port. The controller supplies the duration and routing. Half the return word length gives a ranging unit for this service. A separate passive bank also records each remote mark and subsequent local return during free evolution. These two instruments measure different transmission processes. [Calibration and complete count banks](records.json)

For the signal test the release and hand evolve freely on every seam. The relay arm adds one addressed exchange toward the receiver after encoding. The sham uses the same release for both input labels. A fourth arm cuts every free exchange rate in half at the native contrast. The near station and its one-hop pilot were known at design time. The 52% practical threshold was set before measurement on the two longer routes. Each table is trained on 131,072 labelled launches per bit, then scored on a separate bank of the same size. A single familywise 99% Hoeffding interval covers the clock-only and address-only rates on every positive deadline and flight arm. [Acquisition](measure.py) · [Count-only check](check.py)

## Scope of the timing test

The proportional-travel-time hypothesis uses zero intercept and time in receiver ticks per unit of independently measured addressed range. Tick quantization makes an upper-bound prediction: if one unit clears the threshold by tick one, stations at two and three units should clear it by ticks two and three. Each deadline denotes a separate bank. The clock can miss the deadline, and its interrogation cadence is supplied. The test concerns usable one-shot communication with this detector at the named ticks. A physical limiting speed would require a common signal and ranging standard on freely moving observers, with a justified intervention menu. [Observer criterion](../../../evaluation/EVIDENCE.md#obs-03) · [Lorentz row](../../../evaluation/EVIDENCE.md#ow-21)

For an independent numerical check, represent the two loads as the row and column of a 240 by 240 array and evolve the first interrogation with a Chebyshev series. The implementation uses the public board seam list and shares no signal propagator with the acquisition. The contrast has spectrum in [−32,0], the exchange Laplacian lies in [0,12], and the affine transformation (H+10I)/22 gives spectrum in [−1,1]. At exposure 1/2 the expansion has Bessel argument 11. Its degree-60 tail is bounded by the absolute Bessel power series and a geometric tail. This checks the terminal address probabilities on the first positive interrogation before any subsequent intervention. The binomial prediction intervals compare the first-slot histograms on all station and label combinations.

From the RER repository root run:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.ranged_deadlines.measure
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.ranged_deadlines.check
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.ranged_deadlines.first_service
```
