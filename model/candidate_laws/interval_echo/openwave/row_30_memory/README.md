# Retained seam memory during free matter propagation

An ordinary unit of load can change a written gateway coin during free flight. The mobile double used in the Klein–Gordon guide cannot couple to this coin. Both results use the 99-state Interval Echo seam interaction on the spherical federation. All contacts remain active between the source preparation and the terminal readings. [Protocol](PROTOCOL.md) · [Detector counts](results.json) · [Uniform-background certificate](algebra.json)

## What a patch can read

One standard repair of a single unit across an external seam writes its signed contrast and a coin. A sharp coin reading leaves the unit at one of the two endpoints. At fixed service exposure a terminal patch reads the coin again and tests for a load on an internal seam at the far endpoint. The detector accepts the absolute written contrast after a native echo. Both original coins and all four terminal outcomes enter the counts.

| Screen | Largest probability of a changed coin in seven exposures | After moving the unit to a distant port | With endpoint–coin coupling removed |
| --- | ---: | ---: | ---: |
| 20 cells | 74.4948% | below 0.000001% | 0 |
| 80 cells | 74.4842% | below 0.000001% | 0 |

For example, on the 20-cell sphere at exposure one, the initial coin put the load at the near endpoint in 8,238 of 16,384 simulated trials. In 5,975 of those trials the coin changed and the receiver on the far side reported a load. The initial coin put the load at the other endpoint in the remaining 8,146 trials. The distant controls use the same written coin and start the unit 13–14 or 30–31 blank-seam steps away. The receiver measures occupation of its two ports. These local records establish a physical interaction with the retained coin. A dispersion measurement would also require the momentum and frequency of the same propagating excitation. [All preparations and four-outcome counts](results.json)

The mobile double sees a marked gateway as a scatterer. On the 20-cell sphere at exposure four the far-side detector fired in 223 of 16,384 trials for a blank gateway. A native write on the empty gateway, made before release of the same sharp double, gave zero clicks in 16,384 trials. The corresponding calculated probabilities are 1.3909% and about 0.0000013%. [Source and receiver implementation](experiment.py)

## Why this memory cannot supply the double's internal gap

On a nonblank collar the operator is QGQ, where Q = I − WW†, W is the native contrast write, and G changes the record sign with a contrast-dependent phase. The range of W reduces G. Both blank/nonblank status and absolute written contrast commute with every free term. Other seams can act on shared endpoint loads while leaving those collar projectors invariant. The ten marked collar states with endpoint occupancies (0,0), (2,0) and (0,2) give 30 native columns. On every column the nonblank Hamiltonian is exactly G on the memory, with identity on the load. The test uses integer matrices and compares the result with the registered numerical implementation. [Exact reduction and replay](experiment.py)

Consequently, for any fixed pattern of marked and blank collars the double's Hamiltonian splits into a scalar hopping operator on its position and a commuting sum of collar operators. The blank edges give exchange and contrast energy. Marking an edge removes these terms there, and the degree-reference energy retains the original carrier degree. Changing the pattern can change the spatial spectrum. A precessing sign on a marked edge cannot furnish a second, co-moving component for this excitation.

## What a uniform background of written single-unit records allows

Put absolute contrast one on each seam. The sign-reflection parity σ is conserved on each collar. For a single unit let a and b be the two endpoint positions and c its gateway's coin. In a fixed sign-parity sector the native forward transfer from a to b is −g |0⟩⟨1| / 2 on the coin, and its reverse is the adjoint, with g = 8 for even parity and 9 for odd parity. The corresponding amplitudes are −4 and −4.5. One such edge changes the all-blank single-particle Hamiltonian with a spectator collar by an operator supported at two ports and four collar states. Its perturbation rank is eight. [Block matrices and two-edge check](algebra.json)

Suppose every allowed position shares the same internal space K, so the proposed matter-wave sector contains every position tensored with K. Project free evolution between adjacent positions. The two directed transfers preserve K. Their products include both diagonal projectors on that edge coin, so the transfers generate all four coin matrix units. With both sign parities present, the positive sign factor has eigenvalues 8 and 9 and is invertible. The polar part of each directed transfer belongs to the generated matrix algebra and acts as the same coin matrix unit on both parities. On distinct edges the coin units act on different factors. On M marked contrast-one seams their products therefore generate all matrices on M coins. Every nonzero common invariant K has dimension at least 2ᴹ. Two seams suffice to exclude a common two-component internal space. The full boards have 630, 2,520 and 10,080 seams.

This result concerns an exact position-independent truncation. A position-dependent pattern of correlated records carried with a disturbance, or a collective mode defined from local correlation functions, could have a smaller effective description. The single-gateway experiment measures neither its existence nor its dispersion.

From the repository root, verify both saved receipts with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_memory.experiment
```
