# Binding read on the written four-face fixture

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

Observers can compare energy-filtered arrival records against a separated charge and doublet on the same two-cell fixture. For either charge sign, the released spectrum has exactly **one** level below the independently calibrated parked threshold. A selected energy outcome of multiplicity one transmits at most one of the two independent source amplitudes. Some allowed logical preparation therefore has zero chance of passing the filter, regardless of when the receiver opens. The two-input transfer measured with the unrestricted release depends on higher energies. [Protocol](PROTOCOL.md) · [Receipt](results.json)

## Separated energy reference

The fixture uses the archived 14→24 gateway and four written faces on each side. Parking the sole intercell material link removes its entire rate-4 Laplacian. Exact matrix comparison shows that each occupied-cell block is the sum of the charged fan and the neutral opposite fan acting on independent factors. The charged fans have ground energy about 3.5012284. The opposite fan's independently measured doublet adds 1, giving a separated threshold of about **4.5012284**. When the same source is launched against this parked reference, the opposite fan stays in its vacuum and the charge stays in its source cell. Both joint arrival reads and both opposite-field detectors have exactly zero probability at all times. The executable checks these zeros by exact projector identities.

An integer symmetric elimination counts eigenvalues below rational energy cuts by the signs of consecutive leading principal determinant ratios. Both parked charged fan grounds lie strictly between 3.50 and 3.51. Both released charge sectors have exactly one eigenvalue below 4.50 and exactly one below 4.51, which brackets the true parked separated threshold. Numerically the two lowest released energies are 3.5480114 and 4.6541089. Exact factorization and inertia establish the rank restriction. The probability estimates use numerical eigensolvers. [Executable](measure.py)

## Arrival records

The sole retained level passes between **0% and 1.10716%** of launches, depending on the source amplitudes. After a pass, the receiving gateway records charge with its local doublet on **1.12267%** of trials. The opposite-cell doublet registers on **1.85932%** of passes. With an equal mixture of the two archived input preparations, the subthreshold level alone gives a receiver charge-and-doublet click on **0.0062149%** of all launches.

| Release time | Receiver record | Full release | Below threshold alone | Above threshold alone | Interference in full release |
| --- | --- | ---: | ---: | ---: | ---: |
| 27.0625 | Gateway charge and local doublet | 20.14017% | 0.0062149% | 20.33933% | −0.20538 percentage points |
| 30.8125 | Gateway charge, local doublet, and source vacuum | 13.15899% | 0.0061290% | 12.79328% | +0.35958 percentage points |

These percentages keep all launches in the denominators. The below- and above-threshold columns describe ideal projections with declined outcomes counted as misses. The unfiltered release lets the two amplitude contributions interfere. At 30.8125 its opposite-cell detector registers the doublet on **37.53468%** of the same launches.

The result applies to the two-cell written-field candidate and its parked spectral reference. The [three-arm winding](../../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/echo_turn/intact_exchange/README.md) has different source and receiver records. Testing a longer written-field guide for a bound doublet would require an independently calibrated separated channel on that guide and two source-overlapping energy degrees of freedom whose opposite-field detector weight can be tracked with distance. The present records are simulated.

Run from the repository root to verify the saved receipt.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/binding_audit/measure.py
```
