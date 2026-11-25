# A common cell speed and what its cone can measure

A vacancy and an excess can each transfer occupation across a native external bridge at the exact instantaneous rate **one half cell per exchange-time unit**. The registered head state, the two marked internal contacts and all cell-internal matter energies drop out of the cell-velocity operator. The same construction gives a uniform exponential upper bound on propagation outside that speed and a spatial distance on the cell labels. With the positions and geometric prescription held fixed, an area or volume computed solely from these distances has zero response to the ordinary and collar perturbations. The populations received at a finite time do respond to these perturbations. [Propagation and clock measurements](results.json) · [Source bounds](certificate.json)

## A speed that observers can attain

Let $P_c$ project onto all ports and both coin degrees of freedom in cell $c$. For real numbers $f_c$ define $F=\sum_c f_cP_c$. An occupancy record supplies the cell label and hence $f_c$; for a single-cell indicator, only its receiver's click and complement are needed. The initial derivative of the mean is the expectation of $i[H,F]$. Every term whose endpoints are in one cell commutes with $F$, including both marked contacts. Each external bridge has hopping $-1/2$ times the identity on the heads. Each of its ports meets only this one external bridge. Thus $i[H,F]$ is a direct sum of Hermitian two-port blocks. For the bridge from $x$ in $c$ to $y$ in $d$ the block has eigenvalues $\pm(f_d-f_c)/2$. Consequently

$$
L(F)=\|i[H,F]\|=\tfrac12\max_{c\sim d}|f_d-f_c|.
$$

The coherent incident state $(|x\rangle\pm i|y\rangle)/\sqrt2$ on the maximizing bridge realizes both signs with any head preparation. A blank seam isolated during preparation supplies the two spatial phases. Restore the entire free source for the read. Distinct phase preparations have the same initial cell populations. The receiving-cell click difference divided by twice the elapsed dwell estimates the signed velocity. An independent blank dimer measures the dwell at its short positive branch by its opposite-port probability $\sin^2(t/2)$. Free-flight timesteps serve the calculation; these local clock and occupancy counts supply the estimator. A finite-difference error bound follows from $\|[H,[H,P_c]]\|\le2\|H\|\|[H,P_c]\|$. [Test of the native generators on all three stored boards](measure.py)

| Charge at the stationary 80-cell source | Measured speed | Simultaneous 99% interval | Measured after halving the external hopping |
| --- | ---: | ---: | ---: |
| Vacancy | 0.500034 | [0.49378, 0.50629] | 0.250080 |
| Excess | 0.499764 | [0.49351, 0.50602] | 0.250072 |

The reported intervals include the calibrated clock, eight signed click settings, finite differences and a common union bound. Each phase setting uses $2^{44}$ virtual acquisitions and the clock uses $2^{52}$. Halving only the external hopping in free flight keeps the original emitter, clock and readout. The exact speed for that changed generator is one quarter. Other directed states can move more slowly. In particular a sharp localized release has zero initial occupation current. [Complete synthetic records and changed-law control](results.json)

The distance inferred by admitting every cell observable with speed at most one is

$$
d_{\rm vel}(c,d)=\sup_{L(F)\le1}|f_d-f_c|=2d_{\rm cell}(c,d).
$$

The upper bound follows along a shortest path. The witness $f(z)=2d_{\rm cell}(c,z)$ attains it. The dual diameters of the 20-, 80- and 320-cell boards are respectively 10, 22 and 46 exchange-time units. More generally, a fixed scalar hopping magnitude $g_e>0$ on each disjoint external bridge gives path lengths $1/g_e$ by the same proof. This measurement uses cell labels as positions. The associated distance need not equal a remote packet's arrival time. The commutator distance construction in smooth geometry has a related variational form. [Connes, Proposition 2](https://alainconnes.org/wp-content/uploads/reality.pdf)

## Exponential propagation bound

Start with a packet entirely inside a set $S$ of cells. Write $D(c)=d_{\rm cell}(S,c)$ and $H_\mu=e^{\mu D}He^{-\mu D}$. Only external bridges have a nonzero gradient of $D$, each of absolute value at most one. On each two-port block the norm of the anti-Hermitian part of $H_\mu$ is at most $\tfrac12\sinh\mu$. As the bridges form a matching, the entire norm has the same bound. Differentiating the squared norm of the weighted evolved vector and applying Grönwall's inequality gives

$$
\|P_{D\ge r}e^{-itH}P_S\|\le
\exp[-\mu r+|t|\sinh(\mu)/2].
$$

For $v=r/|t|>1/2$, optimizing at $\mu=\operatorname{arcosh}(2v)$ bounds the probability by $\exp[-2|t|I(v)]$, with

$$
 I(v)=v\operatorname{arcosh}(2v)-\frac12\sqrt{4v^2-1}>0 .
$$

Cell size, internal energies, marked internal couplings, native charge and head density do not enter this bound. The argument applies unchanged to any continuation with the same disjoint external bridges. It proves an **outer** ballistic cone. Establishing a long-distance wave actually travelling on its boundary would also require a lower bound. Absolute first support, measured at arbitrary sensitivity, starts at time zero on the connected finite board. At dwell 4 the actual sharp-launch probabilities at least three cell steps away are $2.65\times10^{-8}$ for vacancy and $1.93\times10^{-6}$ for excess. Both are below the rigorous outer bound 0.272. The different populations permit different finite-threshold arrivals. [Earliest-support calculation](../README.md) · [Shell probabilities and derivatives](results.json) · [Related general estimates for translation-invariant quantum walks](https://doi.org/10.1007/s00023-025-01598-4)

## Which source changes can move a front?

For a fixed protocol take any event $E$ in the complete record. With the two prepared heads as the only input changed, its probability is $p_\chi(E)=\operatorname{Tr}(M_E\chi)$ for a common positive effect. The same applies to a first-detection event, a stopped multi-time record and an adaptive measurement on one copy when its gates are fixed functions of the observed outcomes. A calibrated positive squared-commutator measurement also has this property. An unnormalized signed difference of two outcomes need not obey this likelihood bound.

On the certified equilibrium, each of the two head density blocks satisfies $\chi_b\ge m I_4$ with $m=19999/1000000$. Mixing the blocks equally preserves this lower bound. In the unpairing preparation each left partial trace is at least $2m I_2$ and the common right partial trace is at least $4m I_2$, so each product is at least $8m^2 I_4$. Convex interpolation and the compensating unitary rotation preserve the appropriate lower eigenvalue bound. All normalized source blocks are also bounded above by the identity. For either finite preparation $\tau$ on these paths,

$$
 8m^2\,p_\chi(E)\ \le\ p_\tau(E)\ \le\ m^{-1}p_\chi(E).
$$

The factors are about 0.00319968 and 50.0025. This comparison holds within each charge experiment or a fixed charge mixture. It does not equate the two charges’ finite-population laws or their actual asymptotic speeds. They apply uniformly over distances, appointments, positive detectors and all source-independent ancillary states. When an event probability tends to zero under one source it tends to zero under the other. Apply this to each event $|Z_n-z|>\epsilon$. If a common single-copy readout $Z_n$ concentrates on a deterministic velocity, flight time or normalized correction $z$ in any refinement under either preparation, it concentrates on the **same** $z$ under the other. Any deterministic centering and scaling common to the experiment is allowed. Applying the probability bound to $L_n^{-1}\log p(E_n)$ for any diverging scale $L_n$ also makes upper and lower exponential rates source independent, whether or not a rate limit exists. This conclusion is conditional on the same finite head changing and the measurement protocol having a common physical definition at each level. It does not assume that the 80-cell station is an equilibrium of larger boards.

At the reference, the outward rational certificate proves $-B\chi\preceq\dot\chi\preceq B\chi$ for each matched tangent. The integer bounds are **479 for ordinary matter** and **144 for the collar**. The certificate encloses the true volume-control angle and Bloch axis, checks all corners of a matrix interval by rational Sylvester tests, and includes the certified stationary-source radius. It follows that $|\partial_s\log p(E)|\le B$ wherever $p(E)>0$. In particular the derivative of any log probability divided by $L_n$ vanishes uniformly as $L_n\to\infty$. Finite amplitudes and their source derivatives can carry information even in arbitrarily distant cells. [Exact matrix checks](certificate.py) · [Certificate](certificate.json)

For $M$ independent copies the likelihood factors become $(8m^2)^M$ and $m^{-M}$. The log-rate statement at scale $L_n$ persists for $M=o(L_n)$. If the number of independently prepared heads grows, the single-copy concentration argument need not hold. For example, a first parity record is random, but many independent records distinguish the stationary odd fraction 0.20244 from the ordinary equalized value 0.5. Repeated tomography at each radius can resolve different source means; a distributed material state or source-dependent physical route also changes the hypotheses. Those distinctions matter when proposing an observer-defined macroscopic metric.

Repeated native rail coincidences discriminate collar preparations, even with arbitrary compensating right rotations. Two charge channels also measure shorter time-four mean flight for the same collar intervention, with incompatible ordinary-calibrated response fractions. The distance from the intensive speed operator stays fixed for every independent-copy count. [Detectors, acquisition budgets and certified flight ratios](copies/README.md)

## Relation to the Einstein tower

The exact velocity distance gives the same balls, edges and any fixed geometric weights on them for all the tested head states. At fixed measured volume their area has zero derivative even though the ordinary and collar entropy responses on the selected equilibrium are strictly positive. Its midpoint-subdivision continuation has quadratic mesoscopic packing: the cell graph is the dual graph of the refined two-dimensional triangles. A uniformly bounded-distortion dense parametrization of a proper three-dimensional ball by these cell positions would require cubic packing. Thus this physically measured speed geometry supplies a common bound on signal transport and fails to supply the tower's source-responsive proper ball. Relational positions built from many independent measurements or additional measured conformal scales fall outside this argument. [Quadratic packing proof](../../../../../model/experiments/interval_echo_einstein_station/collar_scaling/README.md#which-metric-do-these-clocks-determine-at-large-scales) · [Signed entropy measurements](../../full_algebra/certificate.json)

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.velocity.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.velocity.certificate --write
python -m pytest -q code/matched_einstein_limit/source_implication/common_front/velocity
```
