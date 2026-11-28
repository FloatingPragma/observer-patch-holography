# Do two material probes see the same change in transverse area?

On the same held change of joint head correlation, a vacancy probe’s image shrank while an excess probe’s image expanded relative to its reference. Both clocks resolved the same appointment, and each observer measured the area of the triangle that received its image. Calibrating both instruments on the same ordinary change also gave opposite signed image responses. Direct radar surveys of the receiving triangle gave different responses in a precision repeat. These simulated records reject a common response for this choice of probes, clock branches, ranging protocol and camera.

## Independent clocks and receiver areas

The 80-cell preparation has two separate written heads. The held intervention varied their real exchange coherence in the even reversal sector. Each head’s own readings and the full parked-energy distribution remain fixed. The endpoints have positive preparation weights. The ordinary comparison used the previously specified energy-changing intervention. All settings, source weights and clocks were fixed before measuring the held area. Earlier data for a different source sector and optical observable motivated the design. [Protocol and input ancestry](PROTOCOL.md) · [Native source and acquisition](measure.py)

The vacancy observer read a blank two-site exchange clock. The excess observer read a parked reversal sign on an independent preparation. The sign read left the flight sample untouched. The two clocks both gave 4.375 to within one part in ten billion of a tick at the equilibrium camera appointment. Each observer independently surveyed ports 39, 47 and 42 using timed round trips, including timeouts. On equilibrium, the inferred triangle areas were 3.892167 and 3.894323 in the respective squared radar units.

Three senders at ports 0, 9 and 7 produced separate images on that triangle. Each observer recorded all three receivers and the outside outcome. The accepted coordinates of a sender’s centroid are its three relative hit frequencies. The determinant of the three centroid columns gives the signed image-area fraction. Multiplying by the separately surveyed triangle area gives the image area in that probe’s own radar units. The equilibrium throughput was 0.244639% for vacancy and 0.215605% for excess. [All starts, clock counts, return histograms, centroid counts and nominal probability laws](results.json)

## Held comparison after ordinary calibration

Each number below divides the held area change per unit preparation step by the ordinary change per unit step for the same instrument. This permits each species its own constant gain and reference. The radar values are multiplied by a million for readability. The first bank and precision bank use the same camera data.

| Area and acquisition | Vacancy | Excess |
| --- | ---: | ---: |
| Image in independently measured radar units, first bank | [0.0056075, 0.0056080] | [−0.0091708, −0.0091568] |
| Receiver triangle alone, first bank, ×10⁶ | [96.69, 101.50] | [−252.75, 273.02] |
| Receiver triangle alone, precision repeat, ×10⁶ | [98.99, 99.60] | [−33.88, 31.61] |

The nominal finite-model ratios for the image are 0.00560778 and −0.00916355. For radar they are 99.2853 and −0.7168 in the scaled units. The primary radar comparison leaves their difference unresolved. After a fresh acquisition of all radar and clock settings the excess-minus-vacancy radar difference lies in [−133.48, −67.38] in those units. The excess radar change by itself remains compatible with zero. The finite-count intervals are rounded outward and have at least 99.1% simultaneous coverage within each acquisition under the simulator's sampling laws. The repeat followed inspection of the initial intervals and nominal results. [Precision design](RADAR_PRECISION.md) · [Precision counts and reuse of primary camera](radar_precision.json)

Removing squared-contrast coupling during flight made both held image changes compatible with zero. Exchanging the two clock records in the active run preserved the difference between the image responses. Providing both observers the same complete probe records left their difference compatible with zero. An independent synthetic triangle of known side lengths recovered its geometric area and reversed its signed image on exchanging two senders. [Instrument and control checks](test_measure.py)

The experiment measures how these particular materials transport and detect the two charges. A clock agreeing on parked copies fixes the measurement appointment. An arrival-event sequence, a transported screen, and a physical relation between capped round-trip distance and proper transverse separation remain to be established. A next test could use independently prepared neutral probes at successive receiving stations, enroll their clocks by coincident arrivals, and compare one predicted change in transverse area at matched arrival distance with the ordinary coupling held fixed.

Run from the repository root with `PYTHONPATH=code:.` and one BLAS thread.

~~~sh
python -Bm openwave.interval_echo.row_19.screen_rulers.measure
python -Bm openwave.interval_echo.row_19.screen_rulers.radar_precision
python -Bm pytest -q openwave/interval_echo/row_19/screen_rulers
~~~

Both experiment modules accept `--counts-only` for a replay of the observers’ readings from the stored counts.
