# Does a freely moving pair have an atomic spectrum?

Two distinguishable port loads bind during free evolution of registered Interval Echo. An observer can resolve a discrete standing wave in their relative separation, using the amount-two port as the moving centre. On the complete 20-, 80- and 320-cell boards the entire bound band falls below the expected second and third hydrogenic shells. This rules out the simple two-load atom on the blank, empty host. OpenWave lists M10 as untested for orbital quantization. [Public criterion and M10 decision](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L468) · [Preparation, tests and controls](PROTOCOL.md) · [Full timing and separation counts](results.json)

## What the two receivers measure

Take one amount-one and one amount-two load on different ports with empty surroundings and blank memories. No centre is clamped. A signed echo against an empty auxiliary reads the amount at each port after free evolution. The terminal spatial record is the graph distance between the two nonempty ports. A supplied 128-position clock entangles its addresses with successive free evolution times, subtracts the phase of two independently prepared uniform single-load states and recombines the addresses. A separate clock comparison with stored loads reads 48 and 192 exchange units for these two free states. Their sum, 240, sets zero binding. The timing interferometer, conditional evolutions, joint addressing and isolated storage ports are apparatus. The relative energies come from the complete native interaction. The terms “amount one” and “amount two” refer to distinguishable port occupations; their identification with laboratory species would require other measurements. [Clock and detector program](measure.py)

On one isolated cell the three rotationally invariant, equal-phase waves at separations 1, 2 and 3 form an exact invariant subspace. In this normalized basis the relative generator in exchange units is

$$
\begin{pmatrix}
 -30 & -2 & 0\\
 -2 & 3 & -\sqrt5\\
 0 & -\sqrt5 & 5
\end{pmatrix}.
$$

Its levels are **−30.12129, +1.63884 and +6.48245**. The first standing wave has 99.632% of its probability at separation one. The other two are above the dissociation reference. For an adjacent-shell launch, the bound-window clock tag retained 1,029,639 out of 1,048,576 copies and 1,025,871 of those read adjacent. Removing the contrast cost lifts the lowest radial line to zero.

On the complete 20-cell board, a sharp oriented emitter seam gave **1,034,018** bound-window clock outcomes in 1,048,576 launches. Of those, **1,032,974** measured separation one. The exchange-only control gave 656 bound-window outcomes. The separate window from −9 to −2, containing both prospective excited shells, recorded 752 launches. Every clock has finite frequency resolution: for spectra in the two rigorously allowed regions, summing the squared Dirichlet sidelobe bounds over that window permits at most 2.36% leakage for any source. The native count is 0.072%. An assigned three-line timing control on the isolated cell put 697,832 of 1,048,576 launches in the same window. These instrument counts preserve the different distance records for every clock bin.

## Why another detector cannot find the two excited shells in this sector

Let $P$ select adjacent addresses and $Q=I-P$. For amounts $a,b$ on the otherwise empty board, the measured uniform single-load energies sum to $6\kappa(a^2+b^2)$. Subtracting that independent reference gives the exact pair generator $H=K-2\kappa ab P$, where $K$ is half the Laplacian of legal swaps. In particular, $K$ is positive and $QHQ$ is positive. For an occupied seam with endpoint degrees $d_1,d_2$ and $c$ common neighbors, the diagonal plus absolute off-diagonal row sum of $PKP$ is $(d_1+d_2)/2+c$ for distinguishable amounts and one less for equal amounts. This counts the permitted swap of the two unlike amounts and both ways to form another triangle. On every federation $d_i\le6$ and $c\le2$.

Writing $G=2\kappa ab$, every normalized vector has energy at least $-G$. All of $P$ has energy at most $-G+8$ for unlike amounts or $-G+7$ for equal amounts. The min–max principle on $P$ and $Q$ therefore gives exactly $\operatorname{rank}(P)$ negative levels (with multiplicity) and excludes every level between the top bound shown and zero. This is a full-sector result, including arbitrary centre motion and detector selection. Many negative levels inside a narrow contact band do not by themselves constitute multiple radial shells.

| Amounts | All negative energies | Negative levels on 20, 80 and 320 cells | Conservative range for hydrogenic $n=2$ on 20 cells |
| --- | --- | --- | --- |
| 1 and 1 | [−16, −9] | 630, 2,520, 10,080 | [−4.000, −2.209] |
| 1 and 2 | [−32, −24] | 1,260, 5,040, 20,160 | [−8.000, −5.959] |
| 2 and 2 | [−64, −57] | 630, 2,520, 10,080 | [−16.000, −14.209] |

The table even allows the dissociation reference $D$ to range from zero up to the Rayleigh quotient of a uniform wave restricted to nonadjacent configurations. This upper bound on the bottom of $QHQ$ is half the number of swaps from $Q$ to $P$ divided by the number of $Q$ configurations. It is 0.05348, 0.01311 and 0.003261 on the three boards. For the mixed pair the second and third predictions $D+(E_1-D)/n^2$ lie in [−8.000, −5.959] and [−3.556, −2.619] even with this finite-volume allowance and without fitting the lowest line. Both intervals lie entirely inside the forbidden gap. The equal-amount results and all three sizes obey the same exclusion for both shells. [Full sector bounds and spatial counts](results.json)

A vacancy and an excess on the unit-filled host provide a separate signed case. Compress the sum of the two independent defect generators to distinct addresses. The native pair differs from this compression only at contact. In an oriented pair and its reverse the difference has diagonal 15.5 and off-diagonal −0.5. Its eigenvalues are 15 and 16. Therefore the actual pair energy is bounded below by the sum of the separate one-defect ground energies on the same board. Repulsive binding above a band is compatible with this inequality. [Native opposite-sign sector](../../../../experiments/material_inertia/neutral_defect.py) · [Independent full-cell matrix check](test_orbitals.py)

## Can a third load mediate the attraction?

Two empty writes at disjoint seams admit an indirect interaction. A third, mobile amount-one load can lower its ground energy by visiting both cores. Each write removes the blank exchange term at its seam and the local eight-unit contrast energy at each endpoint. A wider source distribution could therefore affect the interaction even where the original pair had no contact. The two writes are stationary and supplied. The source compiler prepares the positive ground vector computed for each surveyed Hamiltonian. Receiver phases test these conditioned states. [Protocol](MEDIATOR_PROTOCOL.md) · [Two-core measurement](mediator.py) · [Every placement, phase calibration and count](mediator.json)

Use as reference the lower of the two single-write ground energies measured separately. The two-write generator adds a negative semidefinite term to each single-write generator. The variational principle therefore fixes the sign of its ground-energy shift. On 20 and 80 cells the energy gain is 0.34006 at a two-link midpoint separation and 0.07219 at three. Independent two-quadrature Ramsey banks measure the three-link gain as 0.07220 and 0.07221, with simultaneous 99% calibration allowances of ±0.00067. These readings fix the coefficient of a candidate inverse-distance potential. At four links it predicts a gain of about 0.05415. The computed gains are 0.001893 and 0.000204.

At each held separation choose the wait at which the measured inverse-distance forecast predicts a dark-port probability of at least 0.99979 under the calibration allowance. At four links the fresh 20- and 80-cell clocks recorded **3,119 and 34 dark outcomes per 1,048,576 launches**. At five links they recorded one and six. At all seven separations from six to twelve they recorded zero on both boards. The zero-contrast generator has zero uniform-mode energy for every placement and gives a bright return at every wait. The full scan retains the choice of orientations and all outcomes, including the initial numerically calibrated acquisition. Each compared arm can use a distinct prewritten replica board, with one mobile load shared coherently between arms. The empty board contributes zero vacuum energy in either arm. Very small computed gains at distant placements are limited by numerical precision. The finite scan establishes indirect binding on these supplied cores and rejects the calibrated inverse-distance law at all 18 held settings. It gives no atomic mass to an immobile source.

A common mobile field or a polarizable many-body medium could change the interaction between separated constituents. The single static write and the free two-body test cover different sources. Neither supplies the long-distance attraction or the hierarchy of radial levels needed for an atomic identification. [Static written core](../row_23_coupling/README.md) · [Field-extension static potential](../electron_field/README.md)

## Replay

From the RER root with the scientific Python environment, run

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.measure
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.mediator
python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_31_orbitals
```

The largest propagation uses 57,360 ordered positions. The tests independently assemble the complete 132-dimensional ordered generator from native seam columns on one cell, compare its radial spectral measure with the deposited results, check the filled-host contact correction, and decode the full-federation clock counts.
