# Angular messages carried by a scalar Interval Echo load

A receiver on the southern pentagon of one intact cell can read the orientation and relative phase of a load released from the northern pentagon. The registered one-load generator has exactly one amplitude per port. This gives a concrete scalar control for interpreting signs of mediator spin in observer records. [Protocol](DESIGN.md) · [Executable](measure.py) · [Receipt](results.json)

### What the southern observer sees

Ports 1, 5, 11, 10 and 7 form the source ring about port 0. Receiving ports 2, 6, 8, 9 and 4 lie on the opposite ring about port 3. The sender prepares each load with four native seam splitters and degree-reference phases confined to the northern ring. At the scheduled terminal read the receiver reports one of its five ports or quiet. For each label, 8192 independently simulated copies were read at exposure 1. A fixed likelihood decoder scores every outcome, including quiet. [Source preparation](../../../../experiments/interval_echo_space/shell_source.py) · [Pointer counts and likelihoods](results.json)

| Binary alternatives given to sender | Receiver | Orbital order 1, correct / 16,384 | Orbital order 2, correct / 16,384 |
| --- | --- | ---: | ---: |
| Two standing orientations one rotation step apart | Position | 10,687 (ideal 65.33%) | 8,475 (ideal 51.67%) |
| Opposite circular signs | Position | 8,192 (ideal 50%) | 8,192 (ideal 50%) |
| Opposite circular signs | π/2 seam exchange, then position | 9,858 (ideal 59.96%) | 8,386 (ideal 51.09%) |

The probability of a southern arrival at exposure 1 was 52.3610% for order 1 and 9.2490% for order 2. Rotating either standing preparation by one step cyclically shifts the entire receiving pattern by one port. The two orientations carry respectively 0.202164 and 0.014424 bits per terminal record for equal input probabilities. The amplitude order and the position pattern have different Fourier indices. Squaring the standing amplitude gives density indices ±2 for order 1 and ±4 for order 2, with ±4 equivalent to ∓1 on five sites. Suppressing all hopping while retaining the native diagonal gives zero southern detections at every tested exposure. [Distributions and Fourier magnitudes](results.json)

A rotation eigenstate has uniform intensity around the receiving pentagon. A reflection fixing the poles and reversing both rings commutes with the intact generator and exchanges the circular sources of opposite sign. These two symmetries make the ordinary terminal position distributions equal at every exposure. Adding one addressed exchange between receiving ports 2 and 6 before reading position converts their relative orbital phase into a population difference. Distinct circular signs, a rotating standing pattern and nonzero message information can all occur in this scalar sector. A physical spin or helicity assignment requires a receiver that distinguishes the proposed mediator from this orbital explanation. The run establishes an intracell angular signal on a supplied schedule. [Native generator](../../candidate.py) · [Symmetry and seam checks](measure.py)

### Information across a cell boundary

A separate registered experiment used two sharp-port encodings in a neighboring cell and a thirteen-outcome position receiver in cell 0. Replaying its ten station/windows recovered all retained decisions and every reported probability to within 1.6 × 10⁻¹⁵. At exposure 4 the decoded counts were 45,559 of 65,536 on 20 cells and 47,252 of 65,536 on 320 cells. The no-encoding arm gave 32,909 and 32,698 correct. Both studies have scheduled releases and reads. The intercell trial used a sharp-port bit, so transport of either signed angular mode across a seam remains a separate measurement. [Intercell experiment](../../../../experiments/interval_echo_message/README.md) · [Independent checker](../../../../experiments/interval_echo_message/headline.py)

From the repository root, reproduce the receipt with

```sh
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.mediator_controls.measure --output /tmp/interval_echo_mediator_controls.json
cmp /tmp/interval_echo_mediator_controls.json model/candidate_laws/interval_echo/openwave/mediator_controls/results.json
```
