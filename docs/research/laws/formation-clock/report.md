# The clock of a pair formed in an encounter

**OPH · Pragma Research · 20 November 2026.** Two independently launched constituents met, produced an outgoing record and were subsequently detected together. The internal dial of the newly formed pair gave the forecast local interference signal. **Zero of three moving preparations passed the single-traveler test.** The signed detectors recorded substantial traffic in both directions, and the independently surveyed center did not make the enrolled lap. Keeping the object bound and carrying a clock along one measurable path are distinct requirements.

## What changed in the law

The earlier contact law could transport a prepared compact pair. Its protected subspace was invariant in both directions, so an initially separated state had zero probability to enter it. The successor adds a neutral output wire at each potential junction. The local rule mixes the incoming two-body channel with a direction-matched pair plus one outgoing marker. It treats both clock rails identically. The marker moves away on its own wire. An incoming marker allows the transition to be reversed. Once the emitted marker has left the junction, the two constituents continue under the original contact law. During the enrolled interval no marker can circle its wire and return to the emitter. This proves persistence for the guarded branch of this finite construction.

The marker is an information carrier supplied by the successor. The experiment assigns it no measured energy. The beams, phase reference, empty wires, ideal counters, coherent marker collector and addressed reporting service are supplied. The four-site approach cycle is a feature of this apparatus. The comparison uses timestamps, two-body coincidences, signed arrivals and phase-analyzer counts recorded by the equipped observers. Port distances between sites of the simulator do not serve as physical distances.

## What the observers recorded

**At one aimed crossing.** There were 4,096 independent trials at each terminal age. A marker in the first encounter bin and two coincident constituent detections occurred in 1,994, 2,034, 2,063 and 2,054 trials at ages 0, 1, 4 and 8 after the completed capture update. That gives **8,145 of 16,384 launches**, and zero first-marked separations. Each age used a new cohort because terminal detection consumes the pair. Disabling the vertex or offsetting one incoming beam gave zero first-bin marks. The offset beams did produce later marks, which stayed among the first-bin misses. When the receiver deliberately displaced one constituent at readout it reported all 8,230 first-marked pairs as separated. The detector could therefore report either outcome.

**On the flight apparatus.** One constituent was a smooth transverse packet, the other an independent broad beam. They started on disjoint approach planes. The marker alone went through a coherent mixer. A click at its common output left the birthplace unresolved. Each pair detector operated at a booked appointment, even when the corresponding marker outcome was still being collected. Trial identity joins the two readings. The analysis includes every other marker port and every missed first encounter. Across the flight and geometry instruments, a first marker occurred in about one in 524,288 incident trials. About one in 262,144 such markers reached the common port, giving **one qualified marker per roughly 137 billion incident trials**. The simulation drew aggregate counts representing 1.9854 × 10²⁴ independent attempts. It did not simulate that many interacting bodies.

| Moving preparation | Home appointment | Target x-winding / all home hits | Sealed endpoint phase (rad) | Recorded endpoint phase (rad, 1 SE) |
| --- | ---: | ---: | ---: | ---: |
| Axial | 736 | 511,310 / 629,082 (81.28%) | 0.39875 | 0.39950 ± 0.00201 |
| Diagonal | 896 | 496,041 / 642,693 (77.18%) | 0.39635 | 0.39475 ± 0.00205 |
| Reflected diagonal | 896 | 493,910 / 640,127 (77.16%) | 0.39622 | 0.39667 ± 0.00205 |

The table's phase is the change read by the selected signed aperture between the initial and final appointment. Each of the five coherent phase comparisons, including the two stationary references, agreed with its sealed prediction within one counting standard error. Randomizing the initial rail phase erased the resolved dial. These are tests of the implemented simulated instrument.

**Why the clock test failed.** The registered criterion required 90% of home hits in the intended signed sector. On the forward diagonal, observers counted 496,041 forward arrivals, 132,953 reverse arrivals and 13,699 with zero net winding. The phase-blind survey put the center at just 0.034 and 0.035 laps along the two enrolled axes at the final appointment. There was no observed one-lap crossing of the center, and the phase trackers also disagreed with the signed endpoint. Dividing that endpoint phase by a rest reading would attach a duration to a journey the observations did not establish. Even the central stationary aperture showed a changing average rate across appointments. Its advance was 0.07440 ± 0.00128 rad at 736 and 0.04692 ± 0.00074 rad at 896, in independent trials with the same stationary preparation. A resolved local phase is insufficient evidence for a uniform clock.

## Which components went backwards?

A separately registered receiver combined all four approach-layer inputs coherently at the terminal station. Port zero selected the component whose phase is the same on every layer. At birth each of the four output ports received approximately 25% of the marked pairs. At the diagonal home appointment, port zero recorded **1,977,622 forward and 530,047 backward**. Only 78.86% of its home arrivals had the intended winding. The corresponding axial fraction was 85.47%. On both routes more than 99.99% of backwards arrivals came through port zero. Selecting the smooth longitudinal component alone would still fail the 90% criterion. Phase components interfere at the primary single-layer aperture. The two receiver totals answer different questions.

The finite gate explains how this can occur. For a central transverse mode on a rail of angle φ, let $a$ be the mean cosine over its four routing axes, with $\cos E=a\cos(\varphi/2)$. The equal-slot output of the capture vertex has mean one-step amplitude $a e^{-i\varphi/2}$. Resolving that amplitude into its two unit phases $e^{\mp iE}$ gives weights $(1\pm a\sin(\varphi/2)/\sin E)/2$. The two phases have opposite group velocities. For the central diagonal mode at port zero the forward weights are 78.73% and 78.88% for the two rails. This central-mode calculation is a diagnostic of the specified vertex. The recorded fractions also depend on the packet and detector.

## Implication for the search

The obstruction to forming this particular protected pair can be removed by exporting a distinct record. Its capture vertex supplies the binding channel and its π/4 conversion angle. This establishes an executable possibility with prepared incoming beams. The separate question of what an ordinary unprepared repair history produces is open. The earlier prepared-pair finite-wavelength result retains its original scope. This experiment establishes neither a common metric across inhabitants nor a particle identity. [Previous prepared-pair evidence](../interaction-clock/correction.md) · [One-law plan](../PLAN.md)

The vertex emits an equal sum over routing slots. At a given nonzero transverse momentum that state generally contains components traveling in opposite directions. Birth on a single approach plane also populates all four longitudinal Fourier modes. These are deductions from the source operator. A discriminating next formation test keeps this law fixed and uses longer, initially disjoint wave trains. Admit a finite interval of encounters across multiple planes. Track earlier emission and depletion explicitly, and compare separately recorded output directions, stationary phase advances at multiple times, and the original trajectory criteria on the entire enrolled output. This tests whether the encounter itself selects the motion needed for a clock.

## Reproduce

The protocol and generator were registered at `a903c9af`. Predictions were published at `d48b8119` before opening the six flight panels. From `simulation/experiments/observer_geometry/formation_clock` in the research checkout at revision `bf4a5b1f`, use a Python environment with NumPy. The archived replay used Python 3.13.7 and NumPy 2.5.4.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -B reader.py --self-test
python -B reader.py evidence
out=$(mktemp -d)
for stage in prospect reveal compare
do
    python -B run.py "$stage" "$out"
done
diff -qr evidence "$out"
```

To replay the supplementary four-port measurement, use the same pinned environment and a fresh directory.

```sh
probe=$(mktemp -d)
for stage in register prospect reveal verify
do
    python -B z_probe.py "$stage" "$probe"
done
diff -qr supplement_evidence "$probe"
```
