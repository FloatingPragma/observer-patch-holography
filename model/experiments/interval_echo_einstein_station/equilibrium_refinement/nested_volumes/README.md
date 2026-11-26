# Neighboring volumes from reversible tracer flights

The same five-rim source admits independent volume measurements at radii one and two. One free marker measures available sites. Two freely moving markers detect a change in how those sites are distributed across source preparations, even when their mean number stays fixed. Both instruments recover their respective volumes at two independently clocked exposures. **Their ordinary-calibrated shell responses miss the held entropy change.** A separate reversible aperture experiment reproduces the missing information in nats. [Native flight and error bounds](flight.py) · [All model counts](results.json)

## Bidirectional volume measurement

Use the authenticated one-cell apparatus on the 20-cell host. The center is port zero. The timed balls contain six and eleven sites respectively. The monitored uniform source has 55 admitted two-resident words on the larger ball. Each resident has load two. Add one or two identical load-one tags at distinct vacant sites. All residents and tags evolve together under the active blank-memory Hamiltonian during each flight. The original exterior contacts remain parked. The finite sectors have 660 and 2,970 states. One tag always occupies one slot. Two tags use separate slots and can exchange with residents. [Source, timing and native rights of access](../CLOCK_BALL.md) · [Original one-tag apparatus](../tracer_flux/README.md)

For each radius, launch one bank from the full twelve-site bath and one from inside the timed ball. Pick a uniform subset for the tags and count an insertion failure if a resident already occupies it. On a successful launch the observer opens the shutter after 1/8 or 1/4 generator-time units. The receiver reads the outgoing two-resident word against a freshly and independently sampled word from the same incoming stock. Local comparisons against blank and load-one standards distinguish a resident, a tag and an empty site. The completed tags must all be in the destination. The five retained outcomes classify arrival or miss with or without a source match, plus failed insertion. Source preparation, the independent reference draw and the shutter are supplied. A resident may move and an authenticated arrival can have a different resident word from the one launched.

Let $A_w$ be the set of allowed tag placements in the reference bath for incoming resident word $w$. Let $B_w$ give the placements in the requested ball. For any incoming distribution $p$, put $Q_A=\sum_w p_w|A_w|$ and define $Q_B$ similarly. Conditional on admitted insertion, a particular permitted start has weight $p_w/Q_A$ or $p_w/Q_B$. The free transition table $K_t$ is symmetric because the native Hamiltonian is real symmetric. Thus the measured forward and reverse authenticated probabilities obey

$$Q_A f=\sum_{w,w'}p_wp_{w'}\sum_{x\in A_w,\,y\in B_{w'}} K_t((w',y),(w,x))=Q_B b.$$

The independently drawn word supplies the factor $p_{w'}$ at the receiver. No knowledge of the weights enters the experiment. Where both fluxes are positive the ratio recovers $Q_B$ for every stock, including stocks that change in flight. The observed admission frequency in the twelve-site bath calibrates $Q_A$ independently. With $j$ tags and $v_r$ vacant slots in the timed ball, the measured quantity is $Q_j(r)=E\binom{v_r}{j}$. Its separate forward, reverse and reference-admission counts are deposited. [Acquisition](measure.py)

## What the neighboring volumes predict

Calibration moves one resident from the inner ball to the five-site rim, taking $k$ from zero to one and the exterior completion score $C$ from zero to $\log5$. The held contrast changes a split stock with $k=1$ to equal shares of $k=0,2$. These held stocks have identical individual occupancies, incident energy distributions and outer volumes. Their independently read change in $C$ is $\tfrac12\log(2/5)=-0.458145\ldots$ nats. This is minus the interior entropy derivative at the uniform monitored source.

| Measured quantity | Ordinary before → after | Held split → spread |
| --- | ---: | ---: |
| $Q_1(2)$ outer slots | 9 → 9 | 9 → 9 |
| $Q_1(1)$ inner slots | 4 → 5 | 5 → 5 |
| $Q_2(2)$ outer unordered placements | 36 → 36 | 36 → 36 |
| $Q_2(1)$ inner unordered placements | 6 → 10 | 10 → 10.5 |

The two radii differ by one scheduled exchange. For a sharp stock the tag volume is $v=Q_1$ or equivalently $v=(1+\sqrt{1+8Q_2})/2$. The latter is a different mean on mixtures. Use the difference of neighboring volumes per radial step as the shell response. At uniform stock the derivative of either smooth inner scale cancels in the ordinary-to-held ratio, since the measured outer volume is fixed. One tag predicts zero. Two tags predict $\log5\,(10.5-10)/(10-6)=+0.20118\ldots$ nats. On the two acquired exposures the 99% intervals are contained in $[-0.0029,0.0027]$ and $[0.1980,0.2045]$ respectively. A separate interior-word acquisition gives $[-0.45920,-0.45717]$ nats. Directly differencing the four **finite**, independently calibrated square-root volumes instead gives $[0.1713,0.1821]$ nats across exposures. This empirical sign is also opposite to the entropy response.

The failure has a combinatorial cause. At radius one $v_1=4+k$. For any $j$, the ordinary change in $Q_j(1)$ is $\binom4{j-1}$ and the held change is $\tfrac12\binom4{j-2}$, with invalid binomials zero. Both are nonnegative. Every differentiable function nondecreasing in all these independently insertable marker counts has the same sign for both changes at the reference, whenever its ordinary derivative is nonzero. This includes positive mixtures of marker sizes. Source authentication matters for the volumetry itself. Without it, two tags on the split stock give an apparent outer count in $[36.0220,36.0247]$ at the shorter dwell and $[36.0529,36.0555]$ at the longer one.

These are operational volumes for tag insertion under a measured finite clock. Treating their difference as proper area would require a common spatial metric, proper radial thickness, convergence of the tag measure and control of the variation across the shell. The existing five-rim timing menu includes branching equalities that obstruct an exact smooth geodesic realization. [Metric tests](../rim_geometry/README.md)

## A positive test with material shutters

Place five double-load plugs on spare rails 12–16 outside the cell. While the source is parked, reversible local comparisons on the six inner stations determine $k$. Choose a random order of rim sites. At each step test the next site using gateway 1 and blank 82. If vacant, move one plug there by native endpoint exchange. Stop at the first resident or once all $5-k$ vacant rim stations have been filled. Undo the insertions. Every one of the 55 resident words and its 120 possible orders was checked, including restoration of the complete physical stock. The deposited counters include 21,690 echo comparisons. [Material shutter operations](shutters.py)

When $m$ sites remain after successful closures, exactly $k$ of them hold original residents. The next closure succeeds with probability $(m-k)/m$. The chance of filling every vacancy is therefore $1/\binom5k$. Independent probabilities for sequential success multiply. A continuous additive score on multiplicative admission probabilities is proportional to their logarithm. In nats the measured insertion score is $-\log P(\text{complete}\mid k)=\log\binom5k$. Its mean over incoming trials is precisely $C$ for arbitrary sharp-word stocks. The measured ordinary change lies in $[1.60923,1.60960]$ and the held change in $[-0.45846,-0.45786]$ nats. Taking the logarithm only after pooling classes gives 0.59784 nats on the spread stock, while the mean conditional insertion score is 1.15129 nats. This assay makes the entropy-balancing information physically accessible through pluggable apertures. A mechanical work or stress interpretation would require an independent energy, force or temperature measurement.

## Verification

The flight code assembles all nine exact native seam columns, the original host diagonal and every allowed exchange. It sums a projected Chebyshev expansion with exact rational Bessel-series proposals and a tail bound. The symmetric scaled generator has absolute row sum at most one. The bounded Chebyshev recurrence controls rounding for all propagated sharp starts. The maximum column probability error is below $3\times10^{-9}$. Each volumetry, bulk and shutter bank uses $2^{34}$ independent simulated preparations. Independent native dimer hands use $2^{26}$ starts per quadrature. The reported 99% intervals use a common binomial Hoeffding bound and union bound over at most 2,000 margins and sums, with a further $10^{-7}$ probability allowance in the flight analysis. [Regression checks](test_nested_volumes.py)

From the repository root, set a single numerical thread and run

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.nested_volumes
python -m "$pkg.measure"
python -m pytest -q --pyargs "$pkg"
~~~
