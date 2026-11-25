# A common material entropy reading across three regional instruments

On the marked 80-cell Interval Echo source, observers can recover the signed entropy response from spatial measurements of momentum transfer alone. Three different partial translations each distinguish all 31 trace-one source coordinates. Moving the material instrument therefore preserves the measured source information. Pooling the readings across mounts substantially improves the precision. At the calibrated four-station equilibrium, the common modular response gives the ordinary-calibrated signed area change on all thirty directions of fixed survey volume. [Native effects](effects.json) · [Rank, transfer bounds and reader prices](transfer.json) · [Survey stationarity and ordinary gain](../../../PROOF.md)

## Physical experiment

The two pairs of marked ports, the equal lottery over their ten launch ports, both probe charge preparations, both reversal sectors, and the four-tick material flight are identical for each mount. A mount changes the partial translation of the decorated port graph and the physical edges interfered at readout. Its readings include the supported shift, the momentum flux leaving each side, the force on the marks and the individual oriented local fluxes. The appointment includes misses with zero scored transfer. Observers use separate copies of the source for incompatible readings. [Elementary momentum detectors](../../../../directional_stress/translation/README.md)

| Mount | Centre to target | Orientation | Cells in domain |
| --- | --- | ---: | ---: |
| A | 0 to 19 | 0 | 12 |
| B | 3 to 4 | 1 | 10 |
| C | 3 to 16 | 2 | 9 |

For example, A puts cell 0 in the domain while B puts it in the image. B includes cells 64 and 65 absent from A's domain, while C includes cells 12 and 15 absent from B's. Neither additional translation is A or its reversal. Each mount lies on the same board with the same marked Hamiltonian. The observer does not translate the source together with the detector. [Chart geometry and native apparatus](mounts.py)

## Recovering the source from momentum records

Let $D$ be a trace-zero Hermitian source change on $M_4\oplus M_4$. In each parity block use nine two-head Pauli products and six single-head Pauli products, each divided by two. The final coordinate uses $(I_4,-I_4)/\sqrt8$. This basis is orthonormal. Its coordinates $d$ obey $\|d\|_2=\|D\|_{\rm HS}$. The difference between two material means on chart $i$ is exactly $B_i d$ because launches and readout are affine in the input. The 31 rows of each $B_i$ are actual momentum channels selected solely from the native material response. No parity or head-coordinate reading enters the resulting record vector.

Each real coefficient comes from the native integer Hamiltonian for the corresponding probe charge and head parity. The degree-640 Chebyshev flight uses 225 fractional bits. A selected four-by-four port propagator has error below $10^{-56}$. At $n$ selected ports with detector norm bounded by $M$, each normalized Pauli coefficient has error below $5n\max(1,M)10^{-56}$. The trace contrast sums both block trace bounds. The archive rounds outwards on a $10^{-58}$ grid before normalization of the trace contrast. [Kernel and rational compression](native.py) · [Directed propagation bound](../radial_profiles/README.md)

Let $R_i$ be the displayed algorithm's 55-significant-digit inverse of the native matrix centre. Interval arithmetic verifies $\|I-R_i B_i\|_F\le h_i$ for the exact matrix. The largest $h_i$ across individual mounts is $1.05\times10^{-32}$. Thus $B_i$ is invertible and $\|B_i^{-1}\|_2\le\|R_i\|_F/(1-h_i)$. For example, with any finite source difference, another chart's means can be inferred by $B_j B_i^{-1}(B_i d)=B_jd$. Transfers obey the composition rule $T_{kj}T_{ji}=T_{ki}$ on the entire trace-one source. The decimal implementation uses a rounded centre $C_j$ and $T=C_jR_i$. Its deterministic error per unit source change is at most $\|B_j-C_j\|_F+\|C_j\|_Fh_i$. For all six directed single-mount transfers this is below $3.49\times10^{-36}$ in native detector units. These identities apply to finite differences between arbitrary physical source preparations, independently of stationarity. [Interval inverse and transfer proof](transfer.py)

Combining mounts allows the observer to use stronger channels. Each pooled instrument selects 31 material readings by pivoting on the raw native responses across the named mounts. The table gives rounded upper bounds in native units. The left column of bounds multiplies the Euclidean error of the 31-component detector difference. The right column multiplies the maximum absolute error of any one detector difference to bound the induced error in the modular reading.

| Mounts | Source reconstruction error gain | Modular reading error gain |
| --- | ---: | ---: |
| A | $1.42\times10^{21}$ | $3.19\times10^{17}$ |
| B | $9.16\times10^{21}$ | $1.43\times10^{20}$ |
| C | $2.42\times10^{20}$ | $4.41\times10^{17}$ |
| A, B | $1.09\times10^{11}$ | $2.29\times10^9$ |
| A, C | $1.11\times10^{11}$ | $4.65\times10^8$ |
| B, C | $1.89\times10^{10}$ | $7.51\times10^7$ |
| A, B, C | $2.00\times10^9$ | $5.01\times10^6$ |

The three-mount instrument selects seven readings from A, twelve from B and twelve from C. Its smallest singular value lies between $5.0218\times10^{-10}$ and $7.7394\times10^{-10}$. Its $L^1$ price bound permits an absolute error of $1.99\times10^{-10}$ in each input to contribute at most $10^{-3}$ nats to the reader. This is a requirement on expectations. Sampling costs depend on the detector distribution and the number of preparations.

## Two independent determinations of the modular response

One construction transports the fixed material-and-marginal operator of the reference chart. All seven instruments use the same target obtained from that reader's published prices and response intervals. No area, volume, ordinary coupling or spectrum enters this transport. The native source spectrum provides an independent check. All seven decimal readers differ from $-\log\rho_*$, modulo scalar identity, by at most $1.01\times10^{-5}$ nats per unit source variation. The discrepancy from the material-and-marginal operator on the reference chart is at most $7.83\times10^{-6}$. The comparison with the true entropy derivative includes the possible displacement of the stationary state inside its certified enclosure. The largest disagreement between any two individual mount readers is at most $1.62\times10^{-34}$ nats per unit source variation. At the survey equilibrium the exact transported entropy derivative obeys $dS=-\eta\,dA$ whenever $dV=0$. The reported decimal readers approximate that exact derivative with the stated uniform bound. [Reference instrument](../README.md) · [Full signed operator comparison](transfer.json)

The other construction reads the reference itself from material measurements. Prepare a fair mixture of the four computational head products with the unread left reversal write, so that each parity block is $I_4/8$. This reference stock is reachable by the native product lottery. Subtract its 31 material expectations from those of the equilibrium stock. Reconstruction gives $\widehat\rho=I_8/8+R(y_*-y_{\rm fair})$ in the displayed trace-zero coordinates. The observer diagonalizes the two reconstructed blocks to obtain $-\log\widehat\rho$ directly. For eigenvalue floor $m$ and reconstruction error $r<m$, the resolvent identity gives $\|\log\widehat\rho-\log\rho_*\|_{\rm HS}\le r/(m-r)$. For example, a Euclidean contrast error of $10^{-16}$ at the three-mount instrument gives $r\le2.01\times10^{-7}$, including the deterministic inverse defect, and log error below $1.01\times10^{-5}$ for $m=0.02$. The three-mount material measurements reconstruct the numerical stationary preparation to $2.05\times10^{-10}$ in Hilbert–Schmidt norm and predict the twelve directed radar means to $2.76\times10^{-14}$. Reading the logarithm of this reconstruction gives slopes $-0.06948170435$ and $+0.12554097720$ on the two held variations. The single-mount material logarithms are more sensitive to finite numerical precision. On mount A the second slope differs from direct entropy by $9.44\times10^{-5}$. [Reachable reference lotteries](../../../PROOF.md) · [Independent spectrum reconstruction](held.py)

## Whole-region readings at two nested apertures

The observer can dispense with all the separately read bonds. Restrict each partial translation first to its central cell and then to its domain cells within one graph step. The second domains are $\{0,3,16\}$ for A, $\{1,3\}$ for B and $\{1,2,3\}$ for C. The image and each local port permutation follow the original shift on these cells. The domain and its image therefore remain isomorphic. Recomputing the two cut operators on each restriction makes transfers across its boundary part of the regional flux. The observer measures the whole-region shift, domain flux, image flux and marked force for vacancy and excess, eight numbers per aperture. The experiment retains the same source and detector time. These are radii in the actual translation graph.

Any two mounts at the two apertures give 32 candidate whole-region means. In each case 31 of them have an invertible native response on the complete trace-one source. The worst inverse defect across the three pairs is below $2.83\times10^{-43}$. Pooling all three mounts gives minimum singular value above $3.48\times10^{-9}$ and modular noise gain below $2.61\times10^6$. The pooled regional reference reconstructed from the absolute material means differs from the prepared density by $3.51\times10^{-11}$ in Hilbert–Schmidt norm. It gives signed entropy slopes $-0.06948170477$ and $+0.12554097844$ on the two equal-volume preparations. The omitted full mount, including the local fluxes throughout its larger domain, gives a separate test of reconstructions using only the two small regions on each of the other mounts. All six such held comparisons have relative detector-vector discrepancy below $6.05\times10^{-11}$. [Whole-region Ward operators and native coefficients](regional_effects.json) · [Reference reconstruction and held predictions](regional.json)

## Equal-volume preparations and held charts

Two pairs of native preparations change different mixtures of parity, local polarization and joint correlations. A second source coordinate matches the measured four-station volume separately in each arm. Both signs are within the $10^{-4}$ reachable neighborhood. An independent sparse-exponential flight measures the actual momentum effects for each chart. The readout of the omitted mount supplies a test of the other two mounts' inferred state change. Across the three mount pairs and the two preparations the relative Euclidean error of the 31 predicted held readings is below $5.72\times10^{-11}$. The largest absolute vector discrepancy is below $6.55\times10^{-15}$ in native material units. The ordinary gain is the one calibrated on a separate ordinary source. The equal-volume slopes below are central differences with arm step $4\times10^{-5}$. [Absolute controls and independent flights](held.json)

| Source difference | Three-mount material | Direct entropy | Ordinary-calibrated area |
| --- | ---: | ---: | ---: |
| First | $-0.069481705$ | $-0.069481706$ | $-0.069481705$ |
| Second | $+0.125540978$ | $+0.125540990$ | $+0.125540979$ |

Switching off the marked-head coupling retains the graph, launch ports, charge signs, clock and detector alignment. The blank Hamiltonian then acts identically on all heads in both parities. Each spatial detector numerator factors as a mobile operator times $I_4$ and the marked force vanishes. Every trace-zero source change has zero material response in this control, verified by exact sparse identities. [Changed-interaction check](test_transport.py)

These results identify one common entropy observable on a selected finite source across distinct material supports. The measured four-station volume and area refer to a single four-station survey. Identifying small proper balls, selecting their local boost weights and relating their stress to null focusing across physical scales require additional geometry and dynamics.

Run from the repository root with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.native
python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.transfer
python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.held
python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.apertures
python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.regional
python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/correlation_rank/spatial_cuts/chart_transport
~~~
