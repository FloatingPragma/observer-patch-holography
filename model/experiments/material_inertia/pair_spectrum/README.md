# Multi-cell flight of the lowest neutral-pair excitations

On the registered 20-cell Interval Echo sphere, an excess and a vacancy travel together for at least two cell steps while occupying adjacent ports. The source installs one load per port and prepares the same three sharp emitters used for the earlier recoil measurement. An ideal energy filter selects the lowest isolated cluster of their bound band. All 630 seams evolve freely during flight. The carrier has 30 internal seams per cell and 30 bridges. [Protocol and detector amendment](PROTOCOL.md) · [Earlier preparations](../pair_recoil/README.md)

## What the position records measure

The strict terminal meter accepts a regular-port excess adjacent to a gateway vacancy in a cell at least two cell-owner steps from launch. It retains every other outcome in the sample denominator. A second meter resolves shared ownership including nonadjacent positions. The spectral preparation succeeded on 13.269% to 13.275% of sharp sources on the three listed cells. Launches at the opposite ends of those same bridges were also enrolled. Their admissions ranged from 3.52 × 10⁻⁸ to 9.07 × 10⁻⁷, below the declared flight threshold. The table reports probabilities conditional on filter admission. [Address, filter and terminal records](results.json)

| Launch cell | Initial adjacent pair at ≥2 steps | After 100,000 | After 300,000 | After 1,000,000 | Final qualifying clicks / 67,108,864 |
|---|---:|---:|---:|---:|---:|
| 0 | 1.53 × 10⁻¹² | 0.01966% | 1.1553% | 6.8287% | 4,579,453 |
| 5 | 1.51 × 10⁻¹² | 0.01987% | 1.2864% | 12.6035% | 8,456,695 |
| 10 | 1.59 × 10⁻¹² | 0.01995% | 1.3290% | 18.8158% | 12,628,476 |

Times use the native exchange unit. All nine post-launch increases clear five sampling standard errors plus the numerical allowances. The worst allowance is 4.12 × 10⁻⁶ in absolute probability. The synthetic terminal counts use a fresh random stream for the stricter adjacency test. The initial implementation of the owner meter included nonadjacent configurations in a shared cell. Both sets of readings are archived. The gain in the strict meter above the time-zero tails certifies multi-cell flight within the selected bound sector.

## Which branch moves

The shifted energies of the lowest 24 lines span 36.11562618 to 36.11563856, a bandwidth of 0.0000123782. Adding the filled-reference offset of 96 gives the corresponding formation energies. The next eight lines lie between 36.11715 and 36.11720. Every normalized superposition in the lowest cluster has at least 99.36999% probability on adjacent regular-excess and gateway-vacancy addresses at every time.

The three gateway distances within a cell predict where the low cluster can be prepared. Summing the 24 spectral overlaps over its sharp high addresses gives 1.98739 to 1.98741 per cell on the twelve (2,2,2) cells and 5.11 × 10⁻⁶ to 5.13 × 10⁻⁶ per cell on the eight (1,2,2) cells. The entire cluster therefore cannot be a homogeneous single orbital per owner on this carrier. The two six-cell connected regions of the first type each have diameter three. On the 80-cell refinement, 56 cells of that type form one connected region of diameter 16. The 320-cell refinement has two connected 106-cell regions of diameter 19 and two two-cell regions. These larger regions provide concrete sites for a long-wavelength phase and dispersion measurement. Only their gateway geometry has been surveyed here.

OW-19 gains a source-resolved, multi-cell bound-state flight on the native law. Its remaining measurements are inertial dispersion with a calibrated spatial scale and electrical neutrality with a physical field probe on the same preparation. For OW-26 the flights provide a source-free baseline for subsequent force experiments. They provide no acceleration under a second source, force strength or clock reading.

## Do independent charged probes see the moving pair?

After a million exchange times, flag a common arrival cell at least two owner steps away. Test a randomly addressed port there or elsewhere on the board with a separately supplied excess or vacancy on a two-rail probe. One probe rail contacts the selected port and then interferes with the reference rail. The five signed echo bins include probes that changed occupation during contact. Each shot enrolls the original launch, including spectral-filter and arrival failures. Three disjoint choices of probe rails suffice to find two host loads, using local zero-contrast reads. The source of the fresh charged probe, the terminal shutters, and the phase reference are supplied apparatus. [Probe protocol](PROBE_PROTOCOL.md) · [Three-port circuit](probe.py) · [Complete probe counts](probe.json)

At a contact duration of 1/32 the two-probe score, calibrated on sharp charges, gives +1.0005 for an independent single excess and −0.9991 for a single vacancy when the sampled ports are summed. On the three remote arrivals, regular ports contribute +0.9957 to +0.9997 and gateways −0.9969 to −0.9968. The net is −0.0012 to +0.0029, with sampling standard errors of 0.0025 to 0.0033. At half the contact duration the three net readings range from +0.0057 to +0.0079. Removing contrast from the probe–sample contact changes its unit-charge sensitivity from 0.47197 to 0.00782 at the original dwell, with the splitter, reader and incoming pair held fixed.

The cancellation has a simple origin. Reflection about unit filling exchanges the two probes. Their difference vanishes on a host port and has opposite values on an excess and a vacancy. Once calibrated, the full-screen sum therefore vanishes on **any** one-excess, one-vacancy state. Where the opposite signs accumulate after free flight is a measurement of the bound branch. The probes also have a **common** change relative to an empty host, from +0.2191 to +0.2239 in the same units at 1/32 and +0.1169 to +0.1178 at 1/64. Thus the net sign cancels while the pair scatters both probes. Even in the short-contact limit the common response survives relative to the calibrated charge response. The native diagonal seam energies give slopes of 1.5 for the combined common change and 16.25 for one odd charge unit, a ratio of 6/65. These are local contact measurements. Identifying the interaction with electromagnetism requires an independent physical field and coupling.

## Does a wavelength predict the pair’s kinetic energy?

The cell-owner graph supplies a wavelength standard independent of the pair generator. Its first three nonconstant Laplacian shells have eigenvalues 0.763932, 2 and 3 in inverse squared owner steps. For each of the same three sharp launches, select the lowest unresolved doublet in the 24-line bound cluster. Probe the centre density with an owner-shell amplitude grating. Interfering 128 possible application times gives energy differences on the same bound sector. The rest and final band filters, centre-density grating and timing reference are specified measurement apparatus. The grating can excite multiple lines in an inhomogeneous medium. [Spatial and timing protocol](DISPERSION_PROTOCOL.md) · [Source, reader and controls](dispersion.py) · [Timing counts](dispersion.json)

| Launch owner | Lowest-doublet admission | Measured peak bins, shells 1, 2, 3 | Shells 2 and 3 forecast from shell 1 |
|---|---:|---:|---:|
| 0 | 0.01224% | 19, 19, 23 | 49.74, 74.61 |
| 5 | 0.02617% | 9, 9, 18 | 23.56, 35.34 |
| 10 | 0.30687% | 9, 9, 18 | 23.56, 35.34 |

One bin spans 0.000000374507 native energy units. Every held peak misses the prediction from a single positive kinetic coefficient in ΔE = Aλ. The maximum fraction in any three neighboring bins of the second shell is bounded above by 0.413, 0.457 and 0.409 for launch owners 0, 5 and 10. These upper bounds include the numerical allowance and a simultaneous 99% sampling allowance. The frequency profiles therefore also fail the declared majority criterion. A scalar reference particle on the same owner graph recovers its two held frequencies with the same timer. Halving the material generator's rates after preparation shifts the simulated timing spectrum while preserving the spatial profiles. These results exclude the specified single-line scalar reading of the pair’s centre-density response on this board. A mass from an internally resolved dispersing branch across refinements remains open. The study uses native exchange time and graph steps.

## Numerical checks and replay

Let P be the 284 high sharp addresses and Q the other 57,076 ordered configurations. For shifted energy E, the Schur matrix is H_PP + H_PQ (E − H_QQ)⁻¹ H_QP. A row-sum bound gives H_QQ ≤ 27 on this board. Above 30, the number of eigenvalues of the Schur matrix below E equals the number of native upper-band lines below E. The linear solve residual bounds its operator error by 1.23 × 10⁻¹³ at the five reported counting thresholds. Two thresholds enclose exactly the 24 lines and two others bracket the next eight. The independent count at the filter setting also returns 24.

Resolvent amplitudes at the calibration center, corrected through their first energy derivative, lift a 24-dimensional trial space to the full generator. Orthonormalizing and diagonalizing there gives a Frobenius residual of 4.09 × 10⁻¹². Dividing by the distance to the excluded spectrum bounds the projector error by 5.37 × 10⁻⁹. For sharp-source admission a, the bound on normalized preparation error is twice that projector error divided by √a. At time t, Duhamel's formula adds at most t times the generator residual. These bounds apply to each detector event, including coarse distance bins. Residual norms and spectral margins are evaluated in double precision. The replay is a numerical check of the analytic error estimates. The terminal samples are generated from this law and independently drawn at each appointment. The supplied filter, ideal address meter and exchange clock are explicit parts of the experiment. [Full-space solver and sampling code](measure.py)

From the repository root, replay with single-thread BLAS.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.material_inertia.pair_spectrum.measure --output /tmp/oph-pair-spectrum.json
```

To reproduce the two follow-up receipts, run `python -m model.experiments.material_inertia.pair_spectrum.probe` and `python -m model.experiments.material_inertia.pair_spectrum.dispersion` with the same single-thread BLAS settings.
