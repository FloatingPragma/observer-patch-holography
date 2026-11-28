# Does the regional modular charge weigh the same as ordinary energy?

Two stationary preparations on a native binary contact have equal mean energy and identical endpoint load distributions. The first occupies the symmetric one-load eigenstate. The second uses the antisymmetric eigenstate on eight ninths of trials and an equal mixture of empty and doubly loaded ports on the remaining ninth. The eigenenergies are 8, 9 and 0 respectively. Both banks therefore release ordinary energy tickets with probability 1/64. With fresh contacts on every trial, the downstream material-probe distributions coincide, including all drains and missed tickets. Energy spectroscopy can distinguish the two preparation distributions. [Protocol and order of the tests](PROTOCOL.md) · [Material readouts](reservoir.json)

## A same-energy test of the original reference

For a masked sharp head pair with measured mismatch probability $w$, the two exchange parities have identical reference weight $w/2$. Either equal-load word has weight $(1-w)/2$. The modular difference between the two banks is $\log[(1-w)/w]/9$. The two original source references give 0.23548 and 0.22746 nats. Calibrating on the ordinary change from the balanced zero-energy bank to the symmetric state would assign the equal-energy banks a difference of **8/9 native energy units** for either reference. The actual difference is exactly zero. At the measured face, the assigned change would give a signed endpoint excess of −0.00011574 on each probe. Both actual banks give the same −0.00104167 per preparation. These numbers are predictions of the specified source and terminal laws. [Ordinary calibration and finite proof](../../../../../../../../Lean/ObserverPatchHolography/EinsteinBranch/ClosedOverlap.lean) · [Executable comparison](reservoir.py)

An ordinary equilibrium reference provides a positive control with different parity weights. Give $M$ identical binary contacts exactly $M$ total loads and $2M$ total native energy units. Assign equal weight to the admitted eigenwords. The local weights for energies 0, 8 and 9 follow from exact counts of the remaining $M-1$ contacts. No temperature or modular slope is fitted to the held correlation. First calibrate a slope $\beta_8$ with the zero-to-eight ordinary transition. Independently compare the slope $\beta_1$ from the eight-to-nine transition. The fractional error on a sharp decorrelation energy equals $(\beta_1/\beta_8-1)/17$.

| Contacts | Ordinary slope from 0 to 8 | Ordinary slope from 8 to 9 | Fractional error on the correlation energy |
| ---: | ---: | ---: | ---: |
| 16, 32, 64 | | | Missing parity state makes the local log singular |
| 128 | 0.156595 | −0.133531 | −10.898% |
| 256 | 0.156054 | −0.125938 | −10.629% |
| 512 | 0.139147 | 0.125230 | −0.5883% |
| 1,024 | 0.138126 | 0.141396 | +0.1393% |
| 2,048 | 0.138302 | 0.138603 | +0.01282% |
| 4,096 | 0.138314 | 0.138418 | +0.004449% |

For 4,096 contacts the assigned energies are 0.228651 and 0.295837, compared with 0.228641 and 0.295824. Exact rational bounds on the log ratios enclose the fractional error in [0.000044487, 0.000044488]. Both inputs vary the same sharp mismatch observable. The independent equal-energy challenge leaves a residual of −0.00009298 nats in this bath. Reaching equilibrium is an additional dynamical premise. Each redistribution that preserves contact count, load and additive energy changes the populations of (empty, symmetric, antisymmetric, double) by an integer multiple of (1,−18,16,1). Permuting contact states or resonantly mixing fewer than eighteen at once cannot connect the different parity populations. This restriction concerns the additive isolated-contact energies; interactions during equilibration can store energy themselves. [Exact multiplicities and rational bounds](reservoir.py)

## Preparing an exact reference from ordinary material records

There is a finite way to implement the reference without assuming a thermally mixed bath. Prepare uniformly chosen empty, double, symmetric and antisymmetric contact states. Produce each pilot from the same classical preparation instruction as an untouched target. The ordinary meter has eight units of load-mismatch energy and one unit of odd exchange energy. For each active unit, an independent three-bit shutter rejects the target with probability 1/8. The two pilot observables are measured on separate copies. The readout reconstructs both standing-wave states using the inherited exchange and phase instruments.

A target in energy level $e$ survives with probability $(7/8)^e$. The resulting normalized reference has $K=-\log\rho=\log Z\,I+\log(8/7)H$ on the full four-state contact. This gives **exact equality of calibrated modular and ordinary energy variations for every source on that contact**, including coherent sources. The original paired heads give 0.228641 and 0.295824 energy units, and the equal-energy challenge gives exactly zero. Removing the odd shutter or using one mismatch shutter spoils the identity. On 2,097,152 simulated preparation attempts 1,385,685 targets survive the registered gate sequence. The ordinary energy law fixes the weights of the preparation. This shows how to build the reference; it gives no independent evidence for a gravitational response beyond the energy already used to prepare it. [Material shutter implementation](herald.py) · [Surviving and rejected preparations](herald.json) · [Exact transfer theorem](../../../../../../../../Lean/ObserverPatchHolography/EinsteinBranch/ClosedOverlap.lean)

Conditioning independently prepared survivors on total energy and load gives exactly equal weights to every admitted eigenword, since the product of all survival factors depends only on total energy. At 4,096 contacts the probability of obtaining the declared two totals per independently filled bank is $6.67\times10^{-6}$. This is a selective preparation; all rejected preparations remain in its reported cost. It bypasses the need for an eighteen-contact resonance during preparation. Whether a self-evolving physical region generates this reference is open.

## Where does the energy go in a closed system?

Prepare the original sharp heads at two nonadjacent ports one hop from an anchor. Keep every internal edge active after preparation on a native cell and on two cells joined by their original seam. The two banks differ only in the initial pair covariance. No bond sees an energy difference initially. The time-independent Hamiltonian conserves the complete energy. For every subsequent linear read with outputs $F_{ij}$ on the four sharp starting words, its independent-minus-retained response is $\operatorname{Cov}(X,Y)(F_{10}+F_{01}-F_{00}-F_{11})$. The sign of the regional response therefore needs the actual dynamics. [Source propagation and all apertures](measure.py)

| Enclosed graph radius | Single-cell long-time average | Two-cell mean through 65,536 | Two-cell numerical infinite-time projection |
| ---: | ---: | ---: | ---: |
| 0 | −1.02605 | −0.53681 | −0.27901 |
| 1 | 0 | −2.92240 | −0.81232 |
| 2 | +1.02605 | −1.67644 | +0.17517 |
| 3 | 0 | −1.06073 | −0.81215 |
| 4 | | 0 | 0 |

Entries are native energies per unit of initial covariance, averaged over every enrolled pair. The two source covariances on the original 80-cell placement are 0.01344946 and 0.01740140. The one-cell long-time signs follow from exact rational projections. Diagonalize in spectral blocks algebraically using the square-free characteristic polynomial of the integer native generator. The residues of its resolvent give the projectors on each distinct energy. Lagrange interpolation sums the equal-energy components over the full spectrum without discarding any degeneracy. The two-cell long-time entry is a numerical projection. Its radius-two sign differs from the finite window even at time 65,536. Neither long-time average implies that the instantaneous energy settles. Graph radii count native seams. Their physical distance requires a ruler. [Exact native projection](exact.py) · [All positive and negative shell records](closed.json)

Independently compute the four-dimensional density on each of the thirty single-cell edges from its true stationary ensemble, take its full logarithm, and repeat the two ordinary energy pilots. The zero-to-eight and eight-to-nine slopes differ on every edge both for the formed and the independent source references. The smallest absolute discrepancy is 0.01687 per energy unit. This checks physical regional references at the same locations used to weigh the energy. A supplied Gibbs reference gives equality in the same code. [Thirty-contact comparison](modular.json)

## Astronomical interpretation

The closed region can move energy to an aperture only from another region that initially had energy to supply. Under equal mass conversion and spherical weak lensing the KiDS measurements require between 0.402 and 3.423 trillion solar masses outside a 291.07 kpc aperture across four weighted galaxy samples. The comparison uses favorable observational errors and twice the catalogue baryon mass ceilings. An absolute background of already gravitating overlap stress is a different physical hypothesis and requires its own abundance prediction. For the adopted M31 transport source, changing among the tested faithful heat-bath references permits at most a 1.124-fold outer-to-inner gain ratio even with local ordinary calibration. The steep Chemin profile needs at least 1.888 across the tested source variations and two-error boxes. The Corbelli uncertainties admit the model. [Compensation calculation](../../../../../../../../code/cosmology/rar_deep_regime/compensation/README.md) · [M31 gain comparison](../../../../../../../../code/cosmology/rar_deep_regime/m31_transport/reference_gain.json)

From the repository root, set `PYTHONPATH=code:.` and `OPENBLAS_NUM_THREADS=1`. Run each module below with `python -m`.

```text
model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.measure
model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.exact
model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.modular
model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.reservoir
model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.herald
cosmology.rar_deep_regime.m31_transport.reference_gain
```
