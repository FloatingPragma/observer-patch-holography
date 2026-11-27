# A source record that preserves the active Interval Echo instruments

The traveller and home observers in the active law already post a sharp bit at each clock service. Each post is independently read by the existing mechanical and coherent hands. The comparison sector can listen to those published values without touching the original rails or inserting a measurement before the original read. Its interface takes the resulting pair of bits only. It sends no gate, timing change, entropy value, random number or controller choice back to the native system. [Original publication instrument](../../../experiments/interval_echo_clock/PUBLICATION.md) · [Published 80-cell records](../../../experiments/interval_echo_clock/publication.json)

## Why every old experiment is preserved

Let $h$ be the old public transcript to date, including any settings chosen by an old observer from earlier outcomes. Let $\rho_h$ be its unnormalized conditional state on *all* the old quantum systems and an arbitrary reference. Its trace is the probability of $h$. Once a genuine classical outcome $x$ has been emitted, the old instrument has produced the classical–quantum blocks $\rho_{hx}$. For each block an extra archive can store a copy of $x$ or a random tag with any normalized kernel $K(a\mid hx)$. The extended blocks are $K(a\mid hx)\rho_{hx}$. Summing over the added labels gives back the entire old operator $\rho_{hx}$ for each original transcript. This holds for entangled inputs and reference systems and preserves conditional states even on rare branches. A private quantum memory driven only by the already published transcript also has this property when its conditional updates preserve trace.

At the next step the old observer uses the original instrument chosen from its original transcript. Linearity carries the equality through each conditional outcome map. Induction gives equality of the complete old transcript law and its conditional reference states for every finite adaptive policy, including future coherent operations on old registers. Equivalently one can first run any old experiment to completion and then apply the enrichment to its public copy. The result concerns the entire original system when the newly added degrees of freedom are ignored. An experiment that feeds the new sector back into native dynamics defines a different instrument. A new sector requiring hardware to maintain this isolation would need a physical implementation of the one-way interface. [Conditional record channel](archive.py) · [Adaptive measurement with entangled references](test_archive.py)

Timing matters. Write one bit into orthogonal private states *before* an uncommitted path superposition $|+\rangle$ is reversed. Ignoring the copied path now gives the dephased old density $I/2$. A subsequent Hadamard and zero read returns probability $1/2$ in place of the original probability one. This is a reversible binary-path control on a quantum carrier. The separate 80-cell halfway marker gives a measured instance of premature monitoring changing a later native image. Both illustrate why the enrichment listens only after the old observer's actual publication. [Coherent-reversal test](test_archive.py) · [Native path-marker comparison](../../../../code/matched_einstein_limit/source_implication/full_geometry/station_balls/timelike_areas/exterior_screens/README.md)

## A held change with both ordinary singles fixed

Freeze the 80-cell, width-three, seed-31 journey and its 33 common service addresses. Use the traveller bit $T$ and the home bit $H$ **at the same address**. Uniform reading of the archive gives the observed integer counts
\[
(N_{00},N_{01},N_{10},N_{11})=(12,10,3,8).
\]
Two event observables drive the comparison sector, $T$ and $T H$. The latter asks whether both observers posted in the same service. Counting the two singles separately would lose this information. The three settings below are new choices of sampling measure over the stored addresses. The ordinary setting chooses the unique zero-coincidence joint distribution with the held single rates. Within each selected class all eligible original addresses have equal probability. Both nonblank settings sample only realized, same-service words. Their requested supports differ at 11 by design. No synthetic pairing of bits at different times is used.

| Comparison preparation | Probability weights on 00, 01, 10, 11 | Mean $T$ | Mean $H$ | Mean $TH$ |
|---|---|---:|---:|---:|
| Blank | $(33,0,0,0)/33$ | 0 | 0 | 0 |
| Ordinary | $(4,18,11,0)/33$ | $11/33$ | $18/33$ | 0 |
| Held | $(12,10,3,8)/33$ | $11/33$ | $18/33$ | $8/33$ |

The source API returns the first and last moment in each row as $(\rho_{\rm ordinary},\rho_{\rm joint})$. These are frozen empirical input frequencies for the comparison sheet. Preparing its ordinary and held sectors changes the archive lottery applied downstream. The physical native trajectory, both old single-bit streams and every old observer distribution retain their original law. There is no inference of a universal intensity from 33 temporally related services. The proposed new response to $TH$ requires an independently calibrated material and geometric measurement before it can be identified with collar stress or any astronomical observation.

## Independent live check

The tests execute the native clock journey at the archived seed and at a fresh seed 2026112607, holding the schedule and board fixed. The terminal analyzer is selected from the first five published traveller outcomes. A second replay taps the returned bit only *after* the original coherent hand performs its read. The receiver sees just that returned integer. The complete tapped and untapped native returns compare equal, including both clock words, oscillator outputs, native gate counts, individual service records and random branch draws. All 33 archived same-service address pairs are reproduced in time order. The fresh journey gives $(9,11,7,6)$ in the same four classes, with 6 actual double posts. The extra journey tests the interface independently of the frozen input used for the comparison geometry. [Live native replay and all four checks](test_archive.py)

~~~sh
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/published_strings/test_archive.py
~~~
