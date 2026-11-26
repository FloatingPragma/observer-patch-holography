# A formation-frequency clock search on the neutral emission material

The material that delivers a light neutral excitation to a different observer provides an independent candidate for OW-03. Its geometry and all nine seam signs were specified or sampled for the decay experiment before this clock comparison. A clock carried by an excitation can be read after removing coherence with the material ground. We tested whether the remaining energies support the required oscillation for **any** fixed terminal measurement. [Frozen protocol](PROTOCOL.md) · [Light neutral emission](../../row_09_neutral/extended_decay/README.md)

| Native material | Lowest formation frequency, exchange units | Closest an excited-to-excited frequency comes |
| --- | ---: | ---: |
| Two receiving dimers, connecting collar live (256 histories) | 0.854638 or 1.048574 | At least 0.19393656 |
| Same two dimers, connecting collar cut (256 histories) | 1 | Exact equality |
| Entire joined source and remote receiver, 32 orientation words | 1.263631 | At least 0.00644841 |

Of the 256 histories with a conducting receiver, 128 also conduct at the source. All 512 histories enter the census. For each history the two-dimer energy standard is exactly the isolated receiver used in the decay calorimeter. On the independent owner-crossing trial we instead use the **entire joined assembly** and its own ground, with every contact and the source intact. The bounds apply to arbitrary initial preparations after ground coherence is removed. A comparison involving a further external phase reference or a time-dependent probe would change the clock test. [Complete frequency certificates](frequencies.json)

**Derivation.** In a fixed terminal probability each time-dependent term has angular frequency equal to a difference of two occupied energies. Form twice the native generator on each stock, using the integer diagonal from the physical seams and minus one on every legal orientation flip. Isolate every real root of its integer characteristic polynomial in rational intervals of width at most 10⁻¹⁴ and divide both endpoints by two. Subtract intervals for each distinct excited pair and then subtract the interval for the lowest formation gap. The lower bounds above round the resulting absolute distances down. Degeneracies can only repeat a positive difference or contribute a constant. [Exact calculation](measure.py)

When the receiver coupling vanishes, its diagonal splits as the sum of two one-dimer potentials and the off-diagonal flips are independent. Relative to the ground its four energies are 0, a, b and a+b. If a ≤ b, a transition from b to a+b has exactly the one-dimer formation frequency a, irrespective of the other dimer's occupation. The surviving equality therefore has an independent local oscillator explanation. With the connecting collar live the frequency exclusion covers the complete receiver, including the 128 source-connected histories. It leaves the separate OW-09 observation of light neutral transport intact.

Run from the repository root with one numerical thread.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_03_clock.emitted_material.measure
~~~
