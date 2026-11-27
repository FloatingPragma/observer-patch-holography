# A matter wave on material grown from sharp repairs

The first connector reports one of two contact strengths with equal probability. An observer can use that record to grow either of two uniform materials. Independently timed contacts predicted wall gaps of **15.00062** and **15.49947** in exchange-clock units. After free release, observers measured 34 frequencies and spatial wavelengths across the two materials on three boards. Every lattice forecast agreed within **0.000708**. Substituting the other material's calibration missed by about **0.499** on each board. Single repairs without feedback produced temporal lines too. Their spatial profiles failed the same wavelength test. [Preparation terms](PROTOCOL.md) · [Complete contact, calibration and frequency receipts](records.json.gz)

## How a sharp repair supplies the material

An ordinary write on the sharp load pair 20 or 02 makes both loads one. The head stores their old signed contrast. In either case its two coin-even sign-parity components have equal weights. At the native contrast strength eight, contact with a fixed unit neighbour costs zero for an occupied facing port and costs 32 or 33 for its empty alternative, depending on sign parity. A blank seam splits the probe load between that port and a borrowed vacancy. An exposure of π on the written head followed by recombining the probe reads the parity as a returned or transferred load. A parked wait always returns the initial load. The terminal load read and a conditional swap restore the two incident units, preserve the selected parity and reproduce the same outcome on repetition. These checks act on the 99-state seam for both sharp input directions. The π appointment relies on the declared strength eight and the calibrated exchange clock. [Native read and routing](growth.py)

On a mismatch the observer applies the inverse echo to the restored 11 word. Reading the endpoints again gives 20 or 02 with equal probabilities. A new sharp repair has a fresh chance of matching the first connector. The controller uses a recorded sign, a bound of 24 attempts and the given dimer layout. The first-pointer sampler gave 524,814 returned and 523,762 transferred loads on 1,048,576 simulated starts. All later detector settings retain both first outcomes. The six assembly histories exhausted zero budgets. For a board with N connectors the exact all-start completion probability is (1−2⁻²⁴) raised to the power N−1. It is 0.99999487, 0.99997813 and 0.99991143 for the respective 87, 368 and 1,487 connectors. All sign and reset histories are in the receipt. A physically routed spare vacancy returns along blank seams after each measurement. The complete assembly conserves the stock and leaves the same endpoint word needed for ordinary strand orientation and spectator writes.

## The frequency and wavelength comparison

For a resolved sign, the connector gives an interaction coefficient J of 8 or 8.25. A blank dimer exchanges with h=1/2. Both coefficients are independently measurable through the local contact meter. On each open uniform strand, a load reversal across adjacent dimer orientations defines a wall. The native time evolution gives its squared angular frequency as 4(J−h)² + 16Jh sin²(k/2), where k is the spatial phase in radians per dimer. This follows by applying the physical Ising evolution to even wall products. The first term measures a bulk rest intercept. In the long-wave limit the coefficient of k² is 4Jh. Across all 22 waves with measured k ≤ 0.5 the ordinary continuum frequency forecast agreed within 0.003272. This range was used in the earlier protocol. [Operator derivation](../README.md) · [Acquisition and instrument](measure.py)

| Board | Dimers | Modes per first-pointer outcome | Largest lattice error, returned | Largest lattice error, transferred |
| ---: | ---: | ---: | ---: | ---: |
| 20 cells | 14 | 4 | 0.000708 | 0.000545 |
| 80 cells | 20 | 5 | 0.000633 | 0.000564 |
| 320 cells | 36 | 8 | 0.000659 | 0.000568 |

The receiver combines measured frequencies of overlapping pairs and estimates the spatial waves on even appointments. Both first-outcome materials passed the independent odd-appointment profile test at each size. The control using the nominal launch wavelengths failed on all six materials. The readouts are simulated draws from the computed native outcome probabilities. The layouts were exposed, the sharp material tests were registered locally and the v3 spatial decoder was used without changing its thresholds. The records and seeds regenerate all 212 frequency, clock and spatial count arrays. [Frozen spatial decoder](../confirmation/spatial.py)

## Do single repairs support the same waves?

The first six-dimer open component lies on the 20-cell board and includes ports from two owners. Each of its five ordinary sharp connector repairs was read once. Fresh runs kept every one of the 32 patterns, each of probability 1/32. The complete 64-orientation Hamiltonian retains the actual contact and the two end terms for each pattern. On the two constant patterns it reproduces the physical even-product transfer computed independently from the uniform wall generator. On a mixed pattern the connector contributions leave spatially varying diagonal fields. [Single-write control](SINGLE_WRITE.md) · [Spectral test registered after the impulse](UNSORTED_SPECTRUM.md)

At four clock units the squared transfer length from the central dimer flip into the 91 even wall products is **0.108573** on the pooled probability law. Averaging the lengths conditional on each pattern gives **0.322910**. Both uniform controls have unit length. Every one of the 30 mixed patterns lost resolved length in the simulated terminal counts. The complete simultaneous 99% count allowance is 0.00451. This tests the quadratic evolution of that local observable. [Transfer counts](single_write.npz) · [Decoded probabilities and lengths](single_write.json)

The direct time-and-space measurement used the same ordinary observer on all three ensembles. The two uniform controls scored 79.30 and 55.05 against the odd-appointment spatial cutoff of 121.01. Their four resolved frequencies meet their own independently metered lattice forecasts within 0.000707 and 0.000543. Pooling all single-write histories gave **83,987.74** against the identical cutoff. Temporal pair lines remained measurable. The spatial rejection prevents assignment of the four frequencies to the proposed waves in that ensemble. Contributions above the hand's Nyquist frequency change any normalized signed mean by at most 2 × 10⁻⁹ in the pooled control. Other spatially disordered modes and other instruments remain possible. [Frequency and spatial decisions](unsorted.json) · [Full orientation-space acquisition](unsorted.py)

OW-30 therefore has a measured massive wall dispersion on a strand whose uniform contact response is actively commissioned from sharp writes. Unfiltered repairs fail the enrolled spatial test on the six-dimer component. The public M10.4 pin specifies a different pair and detector. At official revision c531265eb the task directory ends at M10.4 and the coherent-strand confirmation still awaits admission and its public draw. [Public canonical specification](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/openwave/xperiments/m10_oph/research/m10_theory_canonical.md) · [Public tasks](https://github.com/openwave-labs/openwave/tree/c531265eb6085725c6268e156804499c810ba5e4/openwave/xperiments/m10_oph/research/tasks)

From the repository root with one numerical thread run

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
P=model.candidate_laws.interval_echo.openwave.row_30_collective.sharp_growth
.venv/bin/python -m "$P.measure"
.venv/bin/python -m "$P.single_write"
.venv/bin/python -m "$P.unsorted"
.venv/bin/python -m pytest -q --pyargs "$P.test_growth"
~~~

The two spectrum acquisitions accept an optional `--save-counts` directory that retains each full raw array.
