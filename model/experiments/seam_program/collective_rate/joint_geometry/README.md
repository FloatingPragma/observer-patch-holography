# Material timing on the thermal cell

A load signal crosses the coherent 12-port cell on the same preparations used to measure regional modular energy. On seam (0,5), two changes with the same one-port marginals raise the energy by one and eight units. At a material clock reading corresponding to exposure 1/8, the first change delays a receiver signal level by 0.0130424 per unit of preparation change. The second advances it by 0.0199979. Calibrating on the eight-unit change would predict an advance of 0.00249974 for the one-unit change. [Regional modular measurements](../README.md) · [Finite signal probabilities](results.json)

## Two independent records on the same preparation

The reference is the all-blank cell at inverse temperature 1/16 and chemical potential 8. All 4,096 binary words are included with the registered exchange, squared contrast and degree references. Let $D_e$ transfer probability from symmetric to antisymmetric one-load words and let $D_a$ transfer it from the equal mixture of aligned 00 and 11 words into the symmetric one-load word. Tensor both with the normalized identity on the ten remaining ports. The resulting tangents have zero one-port marginals and native energy derivatives 1 and 8. The calibrated eight-port modular derivatives are 0.06253325 and 0.50000345. The prepared densities are the faithful reference plus or minus a sufficiently small multiple of the tangent. The comparison supplies these preparations. The executable gives exact-ensemble probabilities by numerical diagonalization. The common admissible amplitude for both signs is 0.00139872. Finite sample counts would require a separate error analysis.

On fresh copies, the sender at port 1 swaps its load with a sharp empty or full rail. The old load is discarded inclusively. This trace-preserving preparation leaves the source-port marginal unchanged. All cell contacts then act freely. The receiver at port 11 reads its sharp occupation, including zero. The two terminals share two-edge paths through both source ports. Write $S$ for the full-launch click probability minus the empty-launch click probability. A separate isolated blank dimer released with one load has terminal occupation probability $\sin^2(t/2)$, giving an independent material phase reading on these appointments. Its common release and read windows are supplied. At a reference appointment with nonzero $S'_0(t)$ the change in clock time of the crossing of that same signal level is $-\delta S/S'_0(t)$ by implicit differentiation. The signed signal is a density response of this correlated medium. Its threshold defines a material timing observable. It has no established identification with the travel time of a null ray.

| Reference exposure | $S_0$ | $\delta S(D_e)$ | $\delta S(D_a)$ | Clock shift $D_e$ | Clock shift $D_a$ |
|---:|---:|---:|---:|---:|---:|
| 1/32 | −0.00005497 | +0.0000005781 | −0.0000009212 | +0.00016846 | −0.00026841 |
| 1/16 | −0.00020469 | +0.000008435 | −0.00001329 | +0.00141982 | −0.00223683 |
| 1/8 | −0.00063642 | +0.00009263 | −0.00014204 | +0.01304242 | −0.01999793 |
| 1/4 | −0.00131786 | +0.00038609 | −0.00063234 | +0.10422211 | −0.17069789 |

## An exact rank test

The local native Hamiltonian on this sector, after removal of a constant and a multiple of the conserved total number, is $4H'=-\sum_e(XX+YY+17ZZ)_e$. All stages between preparation and occupation measurement preserve load number. For an observable $O$ form $L(O)=2i[H',O]$. Its action on Hermitian Pauli words has integer coefficients. Expand the terminal number operator and take the trace against the changed initial state. For either endpoint tangent the change in the two launches is $-Z_1\otimes D_{05}\otimes I/2^9$. Its pairing with the unpaired Pauli operators vanishes. Odd powers vanish by reality and time reversal. Applying $L$ on both sides of the Hilbert–Schmidt pairing gives the even coefficients with exact rationals. [Exact commutators](short_flight.py) · [Full finite propagation](measure.py)

The two leading responses on the selected route are $5t^4/8$ and $-t^4$. Removing the squared contrast throughout the cell preserves the alphabet and routing and changes those coefficients to $-1/24$ and zero, respectively, in units of $t^4$. Averaging over every ordered pair of terminals disjoint from the endpoints also gives a test without choosing an orientation.

| Port graph distance | Ordered terminal pairs | Quartic coefficient $D_e$ | Quartic coefficient $D_a$ |
|---:|---:|---:|---:|
| 1 | 42 | −2/9 | −11/24 |
| 2 | 40 | 27/160 | 43/96 |
| 3 | 8  | 17/48 | −17/96 |

The first two rows give determinant $-767t^8/34560+O(t^{10})$ for the matrix of finite source derivatives. Hence the two readings have rank two for sufficiently small positive exposure. Every differentiable common scalar source parameter gives a response matrix of rank at most one by the chain rule. This excludes that factorization for these material channels, even if each has its own calibration. A tensor source, direct material scattering, or a separate geometric response requires independent measurements. The port graph, time appointments and source tangents in this test supply no proper-volume match or curvature.

## A four-port change invisible to the energy meter

Use the four-port tangent $P/4096=Z_0Z_1Z_7Z_{10}/4096$, with identity on the other eight ports. Every three-port or smaller marginal has zero derivative. Its one-port and two-port energy response is exactly zero and the eight-port regional modular response is $+3.72076\times10^{-6}$. Resetting port 1 traces out a nonidentity factor and erases the entire tangent. This explains its zero signal on the first route. [Source selection and independent shell measurement](../shell_response/measure.py)

For the second route, require both terminals outside the union of the two source supports and separated by two graph edges. Maximize the number of incidences from the pair to the four parity ports, choosing the lexically first tie. This selects 6 and 11. Resetting 6 preserves the tangent. Its leading signed receiver response is $-9137t^6/5760$; the coefficient becomes $-1/5760$ with exchange alone. At exposure 1/8 the full native response is $-3.10070\times10^{-6}$ and the crossing shifts by $-0.00043656$. Calibrating the same route's response on $D_a$ with its eight-port modular charge would predict $-1.05533\times10^{-11}$ for $P/4096$. The four-body correlation affects load propagation directly. That sensitivity does not establish a geometric response to modular energy.

Replay from the repository root with the following commands.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.seam_program.collective_rate.joint_geometry.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m pytest -q model/experiments/seam_program/collective_rate/joint_geometry
~~~
