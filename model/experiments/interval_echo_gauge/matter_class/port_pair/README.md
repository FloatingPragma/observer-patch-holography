# Reading the angle between two marked half-turns

Cell 0 programs two half-turns of the even doublet, sends a single unit of load along the existing 127-contact route, and cell 13 reads its interference with a third address serving as the odd reference. In the first receiver bank the bright, unused, dark and empty exits recorded **4,955, 26,208, 34,373 and 0** starts. Reversing the second half-turn's polarity gave **34,526, 25,994, 5,016 and 0**. The second doublet address and a rotated receiver gave the same sign changes. Each bank contains 65,536 simulated starts, including every empty reading. [Frozen appointments and controls](PROTOCOL.md) · [All 24 banks](records.json)

The physical preparations use the two marked generators of the specified P1 bracket, a private phase reference, and addressed Interval Echo pulses. The source adds one selected doublet address coherently to the odd reference. At the receiver, undoing that splitter exposes the real part of the doublet's return amplitude. Subtract dark counts from bright counts and divide by **all** starts. The other doublet output carries about 40% of the pair launches. Rejecting those readings would bias the inferred angle. The two addresses, the port labels, their polarities, and the code's even/odd labels are supplied by the apparatus. [Native pulse implementation](measure.py) · [Origin of the two signs](../MARKED_PORTS.md) · [Lean proof for either doublet address](Ramsey.lean)

| Program | Minimum to maximum bright-minus-dark across the four source and frame combinations | Forecast |
| --- | ---: | ---: |
| No turn | 1 | 1 |
| One positive half-turn at port 0 | −0.004425 to +0.000183 | 0 |
| Two positive half-turns at port 0 | −1 | −1 |
| Positive half-turn at 0, then at 1 | −0.448883 to −0.444839 | −1/√5 |
| Positive at 0, reversed at 1 | +0.445938 to +0.450287 | +1/√5 |
| All inter-owner transfers parked | 0; every remote detector reads empty | 0 |

The simultaneous 99% sampling half-width is 0.016084 after rounding outward. Counts are seeded draws from the native pulse predictions on a fixed integer probability grid. A prospective source identification needs the signed physical port couplings calibrated independently of the prepared representation. The chosen P1 line, a compatible P2 assignment and G2 give the same positive-axis geometry on this module. F2 predicts the mirror under the stated marking; deliberate reversal in the control demonstrates the polarity ambiguity. A freely arriving lepton and its spatial spin require separate measurements.

Run the replay from the repository root.

~~~sh
python -B -m model.experiments.interval_echo_gauge.matter_class.port_pair.measure
~~~
