# Can individual source counts and galaxy signals determine the mean dark density?

Two admissible preparations of the same finite record meter have identical total mark injection on every trial, identical one-cell injection probabilities and identical resolved galaxy signals. On the 1,280-cell carrier their **homogeneous record charges differ by a factor of eight** at the same known gain. Including the common, localized galaxy gives a **7.836-fold difference in the volume mean**. The ambiguity comes from the joint occupation of adjacent source cells, which the individual-cell census leaves open. [Protocol](PROTOCOL.md) · [Full measurements](results.json)

## What the native record meter measures

In the registered regional array, each seam reads its two owning cells. If their placement probabilities are $u_i$ and $u_j$, the two-outcome comparator has signed bias $y_e=z(u_i-u_j)/3$ before the pairing fraction $s$. Its probabilities are $(1\pm s y_e)/2$. At $s=1/2$ the modular response per unit weakening of pairing is exactly $y_e\operatorname{atanh}(s y_e)$. Independent seam counters add in the logarithm of their product distribution. Direct propagation on the formed 80-cell collar gives the same material contrast $z=-0.01437036554$ in both preparations. The return readout, reference state and source conversion are held fixed. [Material collar and return calibration](../../matched_einstein_limit/collar_identifiability/refinement/admission.py) · [Regional array and its exact output law](../../matched_einstein_limit/collar_identifiability/refinement/composition/README.md)

Take two nonadjacent ancestral triangular faces. At subdivision length $L$ each face has $L^2$ native cells and $3L$ exterior seams. A source cell inside their union $U$ either supplies a uniform placement $u$ or remains inactive. Each is active with probability one half. A fair preparation bit selects one of two complementary spatial patterns and stays in the record. Conditional seam records are independent as in the registered array. Since the bit's probability is independent of pairing, the modular response is the mean of the two conditional responses. Erasing that bit before taking the logarithm would change the experiment.

Let $B$ be the cut size of the active subset. Every seam with exactly one active owner contributes the same positive value $q_u=(zu/3)\operatorname{atanh}(szu/3)$. All other seams contribute zero. Therefore $Q_{\rm bg}=q_u\,\mathbb E B$ **at finite placement probability**, without expanding the logarithm. The $6L$ exterior seams of $U$ have one active owner with probability one half. Every internal seam contributes its two-cell mismatch probability. Thus

$$
\frac{Q_{\rm bg}}{q_u}=3L+\sum_{(i,j)\ \mathrm{internal}}\Pr(A_i\ne A_j),\qquad
3L\le\frac{Q_{\rm bg}}{q_u}\le3L^2.
$$

Both bounds are attained with exactly $L^2$ active cells on every preparation. Selecting a whole face gives no internal disagreements. Each face's triangles admit a two-coloring across shared edges. Selecting one color on the first face and the opposite color on the second, then complementing both when the fair bit changes, makes every internal edge disagree. The one-cell laws averaged over the fair preparation bit coincide. For each bit value separately, the two policies also inject exactly the same total. The adjacency between active cells is different and is independently measurable. A retained random selection between the two policies fills the interval in expectation. The 80-cell check also exhausts all 70 half-filled subsets and their complements. [Graph enumeration and exact probability calculation](measure.py) · [Independent controls](test_record_abundance.py)

| Carrier cells | Active background cells on each trial | Mean cuts, whole face | Mean cuts, staggered faces | Background charge ratio | Total mean ratio including galaxy |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 80 | 4 | 6 | 12 | 2 | 1.91227 |
| 320 | 16 | 12 | 48 | 4 | 3.86237 |
| 1,280 | 64 | 24 | 192 | 8 | 7.83566 |

These are within-carrier comparisons. The three refinements do not establish a relation between native cell count and physical cosmological volume. The table uses the labelled native 80- and 320-cell boards and the deposited 1,280-cell incidence file, verified against the [model’s registered SHA-256](../../../model/executable_model/refinement.py). [Deposited carrier bytes](FEDERATION_LEVEL3.npz)

## What matched galaxy signals can measure

The transport chart is a periodic $25^3$ box with the existing walk that stays with probability $1/4$ and visits each of six neighbouring sites with probability $1/8$. At survival $7/8$, each active background cell has uniform injection per chart site $1/1024$ and uniform stationary placement $1/128$. One other native cell supplies a compact galaxy by injecting strength $1/4$ at the chart center. Its three detector seams are disjoint from the switched regions. Its full, nonuniform stationary charge is therefore identical under both switching policies, including at finite rates. Direct calculations from both binary probabilities on every seam check this decomposition at radii 0, 1, 2 and 4. [Transport law](../../matched_einstein_limit/collar_identifiability/refinement/transport.py) · [Field calculation](measure.py)

For the gravitational comparison, assign the same positive conversion to both record densities and hold the distance chart fixed. In the periodic Newtonian peculiar-field equation $\nabla^2\Phi=4\pi G(\rho-\bar\rho)$ a constant difference cancels exactly. Under ordinary thin-lens projection it adds a uniform sheet, which cancels from $\Delta\Sigma=\bar\Sigma(<R)-\langle\Sigma(R)\rangle$. On the largest board the two mean total charges were $1.72096\times10^{-8}$ and $1.34849\times10^{-7}$ nats per observer site per unit weakening. At two chart links the inward peculiar acceleration was $2.21914\times10^{-6}$ and the surface-density contrast was $4.70816\times10^{-7}$ in both preparations, with $G$ and the mass conversion set to one in chart units. The maximum difference across every component of the acceleration field and all enrolled projected radii was below $4\times10^{-21}$. Redistributing the same total background injection into a long-wavelength spatial modulation changes both diagnostics. [All settings](results.json) · [Redistribution control](test_record_abundance.py)

A different homogeneous density can affect expansion, angular-diameter distances, growth and absolute magnification. This calculation compares local fields and excess surface densities at a common epoch and fixed geometry. It derives neither a Friedmann background nor a gravitating OPH stress. The address chart, signed comparator, retained preparation bit, and prescribed spatial preparation of the marks are all part of the conditional construction.

## A native moving source whose neighbour statistics are fixed

The stationary screen already has a physical binary occupation $A_i$. A repeatable local read distinguishes a mobile double from a vacancy. During each scheduled service a disjoint set of contacts offers exchange. The 99-state native half-pulse swaps unequal sharp inputs with probability $1/2$, with routing that restores borrowed ports. In a separate, declared readout we retain a sharp snapshot of these flags and use each double as the active owner of the finite collar meter above. The placement $u=1/128$, pairing $1/2$, and formed-collar gain are inherited without fitting; each owner mismatch then supplies $q_u=7.00234\times10^{-10}$ nats per unit decrease in pairing. This is a static read of a frozen or copied configuration. A continuously driven galaxy would need to produce and retain its marks as its matter moves. [Original mobile-double instrument](../../../model/candidate_laws/interval_echo/observations.py) · [Previously specified transport and ordinary preparations](../screen_coherence/PROTOCOL.md) · [New experiment registration](STATIONARY.md)

On a connected three-regular board of $N$ cells with exactly $K$ doubles, the stationary distribution is uniform over $K$-subsets. Randomly selected neighbouring endpoints then differ with probability $2K(N-K)/[N(N-1)]$. Summing all $3N/2$ seams gives the absolute stationary mean $\mathbb E Q=q_u\,3K(N-K)/(N-1)$. At the instrument's original independent Bernoulli preparation density $p=0.16$, the unconditional ratio of expected charge to expected number of doubles is $3(1-p)q_u=1.76459\times10^{-9}$ nats per double per unit weakening. The density $0.16$ belongs to the original screen preparation. Equilibrium at fixed total also follows because each matching can perform a single prescribed exchange or idle with positive probability, and edge exchanges on a connected board generate all placements. [Original stationary-screen analysis](../../../cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.tex) · [Stationary cut mean and variance](stationary.py)

For the independent sharp preparation filling the half-board nearest cell zero, the eleven-service measurement probes the transient. Write $P_{ij}=\mathbb E[A_iA_j]$. On distinct serviced pairs transport each index independently; on the two ends $i,j$ of a single serviced edge their product is invariant, and $P_{ii}=\mathbb E A_i$. These updates predict the complete edge sum $\sum_{ij\in E}(\mathbb E A_i+\mathbb E A_j-2P_{ij})$. Fresh independent material trajectories under the same sharp exchange law give the following cut counts. Uncertainties are one standard error over 4,096 simulated preparations per board. [Predictions, count histograms and charge gains](stationary.json)

| Cells | Forecast cut at service 11 | Held mean cut ± SE | Using products of marginal means | Forecast charge (nats per unit weakening) |
| ---: | ---: | ---: | ---: | ---: |
| 80 | 43.028 | 43.120 ± 0.098 | 40.873 | $3.0130\times10^{-8}$ |
| 320 | 83.223 | 83.188 ± 0.143 | 78.443 | $5.8275\times10^{-8}$ |
| 1,280 | 177.096 | 177.067 ± 0.208 | 167.599 | $1.2401\times10^{-7}$ |

An earlier 320-cell archived preparation supplies a separate check. Given its 59 doubles, the stationary prediction is $1593/11=144.818$ mismatches with standard deviation 6.551. Directly reading all 480 neighbouring pairs in its packed terminal word at service 128 gives **133**, or $9.3131\times10^{-8}$ nats per unit weakening at the fixed gain. This is one realization, 1.80 standard deviations below the conditional mean. The original preparation drew independent Bernoulli flags; conditioning on their conserved total makes the location ensemble stationary from the start. Its acquisition manifest records a modified source checkout, and the replay verifies the archive's SHA-256. No galaxy or CMB observation entered the calculation. Converting these measured binary flags into cosmological sources still requires identifying the physical occupation density, mark routing and lifetime, proper volume, and gravitational mass per unit response. [Original archive summary and provenance](../screen_coherence/measurements.json) · [Direct replay](stationary.py)

## Cosmological normalization

In the stationary 0.16 preparation, a planted double adds 2.04 mismatching seams on average. Attributing each seam locally to its two owners yields measurable delayed correlations and a hard density limit of 3.721 times its cosmic mean. With uniform proper volume per cell this limit is too low for the nonbaryonic mass of galaxy clusters. [Dynamical prediction and cluster-mass measurements](CLUSTERING_RESULTS.md)

Allowing the local number of cells per proper volume to vary still leaves at most three mismatches per occupied double. Against a cosmic preparation with 2.52 per double this permits at most 19.05% enhancement. At the externally supplied Planck normalization, the measured stars and neutral gas in NGC 3741 could supply at most 30.11% of the inferred nonbaryonic mass. Closing the enclosed mass budget would require at least another $3.67\times10^8\,M_\odot$ of baryons. The high-mass cluster mean passes the analogous bound. [Aperture proof, gas uncertainties and cluster comparison](MASS_RATIO_RESULTS.md)

Draco’s enclosed stellar light alone requires a proper site density at least 27,400 times cosmic for a stellar mass-to-light ratio of one. Its occupation remains unmeasured. A separate comparison with its kinematic mass and independent neutral- and warm-gas measurements permits the two-step recovery route with a half-light occupation ceiling of 0.0487 under the stated inventory. [Draco photometry, radio spectra and stellar kinematics](SPHEROIDAL.md)

For scale, the Planck Collaboration's base-$\Lambda$CDM fit to temperature, polarization, low-multipole polarization, lensing and BAO gives $\omega_c=0.11933\pm0.00091$ and $\omega_b=0.02242\pm0.00014$ (individual 68% limits). Their central-value ratio is **5.32248**. This is an observational inference in the fitted cosmology. No source parameter was set from it. A predicted absolute abundance would require the physical density of source cells, a measured stationary joint injection law including nearest-neighbour disagreements, a common proper-volume assignment, and a gravitational mass conversion. Even with the same conversion and ruler, the exact example above leaves the cell-correlation contribution free. [Planck Collaboration VI, Table 2](https://arxiv.org/pdf/1807.06209#page=15)

Run from the repository root with one BLAS thread and `PYTHONPATH=code:.`.

~~~sh
python -Bm cosmology.record_abundance.measure
python -Bm cosmology.record_abundance.stationary
python -Bm cosmology.record_abundance.stationary --archive ../oph-model-sim/runs/L2/universe.npz
python -Bm pytest -q code/cosmology/record_abundance
~~~
