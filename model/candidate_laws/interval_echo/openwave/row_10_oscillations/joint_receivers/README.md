# Two observers recover a winding from their joint records

The charged emitter prepares one of two opposite windings. After conversion and free neutral flight, the receiving owner reads three local loads. When a vacancy and an excess arrive together, the conversion owner can use a borrowed vacancy to interrogate two original direction records. In nine of twelve enrolled comparisons the conditional frequency of a successful local read distinguishes the winding signs. For either sign, the distant pair alone has the same arrival probability. [Frozen design](PROTOCOL.md) · [Count banks](counts.json) · [Count decoder](readback.py)

At receiver 12, appointment 256 and positive analyzer setting, the positive winding produced 2,618 successful conversion-site reads among 17,144 distant pairs. The negative winding gave 1,549 among 17,088. The difference is 6.21 percentage points, with simultaneous 99% limits of 3.87 to 8.54 points. The pair is rare. The two independent cohorts each contain 2⁵² original simulated starts, sampled in aggregate from the computed outcome law. Refused launches, missed charged arrivals, donor failures, missed pairs, wrong local flags and local probe misses each have a place in the record. [Acquisition](measure.py)

## What the conversion owner does

The distant owner halts and reads its three physical trits at the specified appointment. The conversion site is held during the exchange of classical records. A remote pair signals that the two original source leaves are unit-filled. The conversion owner can open known unit shutters and bring an existing vacancy from stock port 77 to hub 160 along a route avoiding both distal owners and all six original direction memories. The compensating stock stays parked. The following seven-port reader is entirely within owner 13. These addressed holds, pulses, shutter controls and borrowed stock are supplied apparatus. The appointments refer to scheduled exchange time.

The vacancy explores triangles (165,164,159) and (158,162,166). In each triangle one edge carries the original written direction, with a two-edge blank return path. A split, a calibrated phase on the second branch, one marked traversal per branch and a recombination prepare a superposition of the two one-written-head words. Reversing these native pulses and checking for a hub vacancy with two unwritten flags measures that superposition. The flag tests use sharp vacancy positions and their inverse motion. Either value of the unread third head gives the same effect. [Native pulses and access checks](apparatus.py)

For reference phase φ = ±2π/3 the success effect is |v⟩⟨v| where v has coefficients 1/2 and exp(iφ)/2 on numeric head masks 1 and 2. Its eigenvalue on the matching normalized superposition is 1/2. This efficiency follows from the material splitter and the accepted hub aperture. Both phases and the two marked crossings were selected to test the forecast of a coherent source reference. Separate local standards give at the positive panel 1,047,999 successes in 4,194,304 original matched-source starts and 262,466 on the opposite reference. Roughly half the forward preparations refuse. The sharp stocks gave one-quarter success within counting error for either single-written word and zero for blank and double-written words. [Independent ideal comparison](test_joint.py) · [Source-reference forecast](../ERASER.md)

The park-and-release audit replays the actual movements of the controlled stock before and after conversion, covering 79 distinct blank edges. None addresses an original marked memory. The registered full donor pulse acts on separate contacts and returns the selected charged/target/donor word with unit modulus. All shutter operations used while staging touch known stocks and distinct memories. On this addressed schedule each operation therefore either has the same scalar on the retained directions or acts entirely on the other factors. A bath dephasing the held source would alter the experiment.

## Effect of the two controls and changed contrast

The entries give differences in local success frequency, positive minus negative winding, in percentage points conditional on a remote pair. Each pair of columns gives the negative and positive phase panels. The third and fourth columns use the native law and intact instrument. The last two use contrast four during neutral flight. Every displayed value comes from distinct count banks. [Complete intervals and controls](readback.py)

| Owner | Appointment | Native − | Native + | First head read − / + | Second crossing skipped − / + | Changed − | Changed + |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 256 | −5.49 | +6.21 | +0.17 / +0.07 | +0.36 / −0.05 | −9.72 | +10.03 |
| 12 | 512 | −6.10 | +2.49 | +0.20 / −0.19 | +0.20 / −0.08 | −13.26 | +14.11 |
| 12 | 2048 | −3.80 | +4.71 | −0.13 / +0.25 | −0.01 / +0.03 | −14.37 | +15.23 |
| 17 | 256 | −0.10 | −1.18 | +0.06 / +0.08 | +0.09 / −0.01 | −9.19 | +10.21 |
| 17 | 512 | +1.02 | +1.33 | +0.20 / −0.10 | +0.05 / −0.03 | −9.90 | +10.43 |
| 17 | 2048 | −2.51 | +5.26 | +0.30 / +0.57 | +0.01 / −0.11 | −15.10 | +15.38 |

The simultaneous binomial intervals resolve zero of twelve winding differences after a sharp first-head read and zero of twelve for the selected success click after suppressing the second marked crossing. In the latter circuit that click reads a single sharp word. The changed contrast resolves twelve of twelve. All 48 far-only comparisons include zero. The deciding observable is a correlation with an accessible source reference and depends on the neutral flight. Identifying an oscillating neutrino family requires three weak production/detection tags and held oscillation parameters. The original neutral receiver by itself cannot distinguish these two windings. [Receiver theorem](../theory/PROOF.md)

Run from the repository root with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.joint_receivers
python -Bm "$p.measure" --acquire /tmp/oph-joint-receivers.json
cmp /tmp/oph-joint-receivers.json model/candidate_laws/interval_echo/openwave/row_10_oscillations/joint_receivers/counts.json
python -Bm "$p.readback"
python -Bm pytest -q --pyargs "$p.test_joint"
~~~
