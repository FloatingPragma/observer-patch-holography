# A polarized quantum sent and read by independent charged circuits

**A guided charge prepares an almost pure one-photon packet in the finite quantized Maxwell field. A second charge distinguishes its two circular preparations at a distant station using only its paths and its clock.** On the side-four box, conditioning on a local source pointer gives one field quantum with probability 0.9999999822. At a receiver wait of four field units the two preparations give opposite differences in local pointer probabilities. The measured difference of differences is (1.839 ± 0.084) × 10⁻⁶, with one nominal sampling standard error over independently simulated trials. The much larger quiet and background populations are included. [Design](DESIGN.md) · [Executable](measure.py) · [Complete outcomes](results.json)

## What the charge prepares

The open boxes use the archived Wilson–Maxwell vertex–edge gradient, edge–face curl and oriented Peierls phase. Their incidence arrays match direct replays of the original constructor. The measured field is the transverse occupation space through total number three. A directed hop of a unit charge writes the unitary exp(igQₑ) on the field in the hop's unit singular channel. At every corner the apparatus aligns the spin with the next addressed bond. The π/2 gate swaps the charge between its two endpoints and carries the field phase. Each branch takes six steps around one rectangle and returns to the same material and reference state. The relative phases of four quarter-rotated circuits specify the circular setting. The local analyzer recombines the paths and reports one of four outcomes. [Archived field and ordered link operators](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/native_wilson_maxwell/field.py) · [Archived Wilson bond and its spin channels](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/native_wilson_maxwell/matter.py)

Write Uⱼ for the *ordered* field phase of path j. The field state associated with local path result k is the sum of Uⱼ|vac⟩ weighted by exp[i(m−k)jπ/2]/4, where m is the chosen sign. The paths return the charge to its common starting vertex. The outcome k = 0 prepares opposite quarter-turn characters for the two settings. The one-quantum parts are orthogonal. In a separately enrolled number calibration the smaller box gives a factorial correlation ⟨N(N−1)⟩/⟨N⟩² of 1.0664 × 10⁻⁷ on this local herald. A mean-matched calibration mixture gives one. The remote experiment can operate with both field-number counters omitted. The conditioned emission probabilities include all four number outcomes, including the small three-quantum component. [Source probabilities, number correlations and controls](results.json)

## Flight and number-resolved absorption

The end stations are separated along the z axis and share no vertices or edges. The charges park while the free field evolves. In the one-quantum sector its native clock gives the positive frequencies of the square root of the transverse curl-squared matrix. For each mode, an independent circuit at the far end traverses its six links in reverse. The field-wide number meter identifies arrival of vacuum after the receiver pulse in these trials. Across 2⁴² side-four simulated releases of the positive preparation at wait four, 221 trials returned vacuum with the matching charged-receiver pointer. All opposite pointer counts were zero. This coincident event is rare. It uses a number read after emission and another after the driven absorber.

| Side in lattice units | Transverse modes | Lowest observed generator frequency with nonzero emitted weight | One quantum with source pointer 0 per launch | Largest sampled joint matched absorption | Wait of maximum in [0, 2n] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 176 | 0.874032 | 1.63213 × 10⁻⁴ | 2.15252 × 10⁻¹⁰ | 7.25 |
| 6 | 540 | 0.629384 | 1.63039 × 10⁻⁴ | 7.26182 × 10⁻¹¹ | 12.00 |
| 8 | 1,216 | 0.491151 | 1.63020 × 10⁻⁴ | 5.87033 × 10⁻¹¹ | 16.00 |

All three lowest frequencies agree with 2√2 sin[π/(2(n+1))] in lattice units. The finite-box gap closes in that prescribed Maxwell family at fixed spacing. The positive-frequency packet has a nonzero equal-time overlap with the receiver. On the side-four box, free evolution raises matched absorption at wait 7.25 to 629 times its equal-time value. Suppressing propagation preserves the equal-time response. Turning off either local charge coupling removes the joint event. Longitudinal gradients have zero overlap with the transverse coordinates. The largest two boxes reach their largest surveyed absorption at the final window. These three maxima do not determine a signal speed. [Dispersions, time survey and controls](results.json) · [Archived fixed-spacing spatial family](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/field_record_addenda/native_wilson_maxwell_spatial_family/PROOF.md)

## Reception using two local pointers

An instantaneous four-path receiver gives almost identical records for the two circular signs when the final field count is unread. The largest difference in total variation over the surveyed single-pulse stations and waits is 6.1 × 10⁻¹⁵. At side four and wait 7.25 the same receiver distinguishes the prepared packet from its vacuum calibration by 2.64 × 10⁻⁶. Reading only the charge at one time combines amplitudes for removing a quantum and adding a quantum.

Give the receiver a clock qubit that selects between the same circuit at an early or a late time, with relative phase i. Recombine the two times locally and retain both clock outcomes and all four path outcomes. The clock separation π/(2ω₀) is independently calibrated from the lowest box frequency. Closing paths does work through external controls. The clock read accesses relative field phase through the two interaction times. All field occupations are traced over. On comparing source and receiver records after the run, use the source path result k = 0 to identify the packet. Averaging all four local source outcomes makes the remote marginal exactly independent of the sign. The intact joint table exhibits the difference.

| Side | Early wait | Ideal probability contrast for + and − at the receiver | Measured difference in units of one simulated standard error |
| ---: | ---: | ---: | ---: |
| 4 | 4 | +8.8681 × 10⁻⁷ and −8.8681 × 10⁻⁷ | +21.8 |
| 6 | 6 | −9.0819 × 10⁻⁸ and +9.0819 × 10⁻⁸ | −1.8 |
| 8 | 8 | −5.3075 × 10⁻⁸ and +5.3075 × 10⁻⁸ | −2.5 |
| 4 | 7.25 | −1.4468 × 10⁻⁶ and +1.4468 × 10⁻⁶ | −33.9 |
| 6 | 12 | −2.5418 × 10⁻⁶ and +2.5418 × 10⁻⁶ | −60.9 |
| 8 | 16 | −1.4694 × 10⁻⁶ and +1.4694 × 10⁻⁶ | −34.7 |

The observable subtracts path-3 clicks from path-1 clicks for one clock outcome and reverses that subtraction for the other. Each setting enrolled 2⁴⁹ simulated launches per sign using exact multinomial aggregation. The lower three rows reuse the maxima of the independently calculated earlier number-resolved survey. The fixed-side waits and their weak larger-box count contrasts are also reported. With both circuits at the same field phase the two-time contrast is zero. A direct calculation on the span of the evolved source links and both sets of detector links keeps every reachable Fock occupation through the chosen cap. It needs no factorization of the field into independent spatial detectors. [Local pointer probabilities and raw count multiplicities](results.json)

The object identified here is a finite Maxwell-model photon with a measured quarter-turn polarization in a specified geometry. The cochain metric, quantization, occupation cutoff, step clock, coherent local path controls, and compensation of the unused Hamiltonian terms during addressed hops are supplied. During flight both charges are held fixed. The geometrically disjoint charged circuits act on one transverse field subject to a common occupation cap. Strict local commutation of all projected operations and a light cone are additional questions. The calculation uses the same field throughout preparation and reception. It establishes a controlled comparison additional to the archived autonomous matter-pair response. Relativistic helicity in an unconfined continuum and an identification of this quantum in registered Interval Echo require further measurements on their respective sources. [Earlier current-sourced Interval Echo extension and model comparison](../README.md)

From the repository root, reproduce the probabilities and count tables with a fresh output path.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/photon_identity/quantum_transfer/measure.py --output /tmp/finite-photon-transfer.json
```
