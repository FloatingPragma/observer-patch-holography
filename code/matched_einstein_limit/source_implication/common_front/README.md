# What the earliest arrivals measure on the 80-cell collar

Both probe charges define the same reciprocal, source-sensitive distance when two stations compare their earliest **occupation probabilities** in both directions. The three-hop clock sees both the ordinary and held pairing interventions on the same preparation. The four-hop clock gives the original four stations a nondegenerate tetrahedron. Its independently measured area response fails the ordinary-calibrated entropy balance on the held pairing intervention, even at its **own** fixed tetrahedral volume. This separates discovering a common clock from establishing the metric and stationarity required by the Einstein tower. [Certified intervals](certificate.json) · [Measurements](results.json) · [Equilibrium preparation](../full_algebra/README.md)

## An absolute arrival time and an arrival coefficient give different information

Take the certified full-algebra equilibrium, with faithful four-dimensional head blocks $X_b$ for even and odd parity. Leave both heads in place. Launch one native probe at port $i$ and record whether it occupies $j$ at time $t$ on a fresh copy. All misses count. The 80 cells have 960 ports. Write $d$ for their shortest unweighted seam distance and $B_{ji,b}$ for the $ji$ head block of the $d$-th power of the off-position hopping matrix $-2H_{\rm off}$. The first nonzero term in the propagation amplitude is $(it/2)^d B_{ji,b}/d!$. Every shortest route is simple. Blank hops have identity head matrix in these units, and each of the two marked contacts acts on a separate coin and can occur at most once on a shortest path. The marked hop matrices are nonzero and entrywise nonnegative after extracting the common sign. Each route therefore contributes a nonzero matrix, with addition across routes preserving a positive entry.

It follows that the single-appointment probability is strictly positive for all sufficiently small positive times for every connected pair and either charge. Its first temporal support time is zero. The first nonzero power of the probability is $2d$. Any smooth clock change with a finite positive rate preserves this exponent. A ball defined solely by these hop orders has exactly zero area or volume response to either physical intervention. This conclusion covers all graph distances on this board.

For a more informative timing comparison, choose the sending end uniformly at random and pool orientations for **each charge separately**. For a given parity and positive charge let $B=B_{ji,b}$. Reversing the link transposes $B$. Flipping the charge transposes the individual directed marked hops, and their separate coin factors commute on each shortest path, giving $B_{ji,b}^{(-)}=B^T$. Consequently the pooled leading probability has exactly the same coefficient for both charges. Set $O_b=B^T B+BB^T$. For endpoints $k$ hops apart,

$$
 P^{\leftrightarrow}_{ij,q}(t)
 = {c_{ij}\over 2^{2k}(k!)^2}\,t^{2k}+O(t^{2k+1}),\qquad
 c_{ij}=\tfrac12\sum_b {\rm Tr}(X_b O_b)>0 .
$$

Integer matrix powers check the charge identity across every enrolled link. As an independent timing check, direct continuous-time evolution of the two charges on port pair 0–39 gives a rescaled sixth-order coefficient within 0.044% of the formula at $t=0.005$. Finite-time occupation probabilities differ by charge. [Implementation](measure.py) · [Direct evolution test](test_common_front.py)

Compare the first small time at which the ensemble probability reaches a common threshold $\epsilon$. Analyticity and the positive leading coefficient give a unique local crossing. A unique all-blank $k$-hop control path on the same board has $c=1$. Dividing the two crossing times and taking $\epsilon\downarrow0$ gives the positive reciprocal link length $w_{ij}=c_{ij}^{-1/(2k)}$. The control paths for orders 1 through 4 are respectively 2–3, 12–44, 13–792 and 21–792. Their blank status and uniqueness are checked on the unweighted graph. A finite-copy experiment needs the mean number of threshold clicks $N\epsilon$ to grow for a vanishing relative counting error.

For each order separately use all pairs at that order as links, and take least sums of their lengths. Each of the four resulting graphs is connected. This builds in a common order, bidirectional appointments, relay composition, and repeated preparation. Comparing absolute threshold times across hop orders without relays would violate the triangle inequality at small thresholds. For example, a two-hop threshold time scales as $\epsilon^{1/4}$ and the sum of two one-hop times scales as $\epsilon^{1/2}$.

## Same-source geometry and held transfer

The four survey ports are 0, 1, 38 and 39. Port 5 is held out. The ordinary direction redistributes the two parity blocks equally while preserving their sum. The pairing intervention takes each left parity marginal times the shared right marginal and subtracts the original block. A rotation about the right marginal's Bloch axis supplies a volume control with zero entropy derivative. [Physical intervention and rotation](../signed_transfers.py)

| Order | Links | Geometry and response |
|---|---:|---|
| 1 | 2,520 | The ordinary direction changes the left marked link. Pairing preserves every leading length. The four original stations lie on a line in these distance units. |
| 2 | 3,600 | The four-station centered Gram matrix has a negative eigenvalue at the population center. |
| 3 | 4,704 | The direct 0–39 length is in [0.29311137, 0.29311143], while every path of two or more links costs at least 0.58622182. On the previously registered cap-six survey's fixed-volume directions its slopes are [0.06500, 0.06607] ordinary and [0.03576, 0.03621] pairing. |
| 4 | 6,672 | The four-station distance Gram is positive definite after anchoring one point. The tetrahedron has area [6.2337402, 6.2337431] and volume [0.6342285, 0.6342298]. Adding port 5 gives a certified negative Gram direction, excluding a common Euclidean embedding of these five distances in any dimension. |

At order three every shortest 0–1 path passes through 39 and every shortest 39–38 path passes through 0. If those triangle equalities arose from unique smooth geodesics on a proper ball, they would place the points on one continuing geodesic in order 38, 0, 39, 1. When this entire segment is within the minimizing scale, its endpoint distance must be its length. The measured 38–1 distance falls short by at least 0.56813. This excludes that local minimizing interpretation for this four-station menu. A negative Gram eigenvalue by itself only excludes a flat Euclidean embedding.

At order four measure area and volume of its own tetrahedron from its six ranges. For each raw direction $D$, use $D+\alpha R$ with $\alpha=-V'_D/V'_R$, where $R$ is the same marginal-axis rotation and $V'_R<0$. The two resulting directions therefore keep precisely this tetrahedron's volume fixed to first order. The certified signed responses are

| Intervention | Entropy derivative | Own-volume area derivative |
|---|---:|---:|
| Ordinary | [0.4094266, 0.4094802] | [−0.00492282, −0.00486967] |
| Pairing | [0.2310794, 0.2311330] | [+0.0000159354, +0.0000161264] |

In population units the ordinary standard gives $\eta=-S'_o/A'_o=83.62566$. The ratio of the predicted pairing entropy change $-\eta A'_p$ to the independently measured $S'_p$ is in **[−0.005871, −0.005732]**. The required value is $+1$. Direct finite changes with numerical volume matching at step $10^{-6}$ reproduce the area slopes. Larger steps at $10^{-4}$ exceed a local bracketing rule on one sign for each intervention and are recorded. This sign obstruction concerns this clock, four-station tetrahedron, rotation, source and entropy. A finite curved embedding, another observer instrument and a macroscopic refinement remain distinct physical questions.

The certificate propagates the equilibrium enclosure in Hilbert–Schmidt norm through integer hopping matrices, whole-graph path selection, tetrahedron derivatives and the signed transfer. It bounds competing shortest routes with directed integer rounding, cancels shared link types to certify differentiability, and evaluates Gram inequalities with rational arithmetic. For pairing, the fixed marginals reduce the variable link effects to two types. Their two-by-two Jacobian determinant retains the sign under interval evaluation. [Source and verification](certificate.py)

A separate cell-velocity measurement gives exactly one half cell per exchange-time unit for both charges, with an attainable local current and an exponential outer cone. Its velocity-dual distance is twice the unweighted cell distance. On the faithful heads, a uniform likelihood bound proves that any positive single-copy arrival instrument concentrating on a definite rescaled speed or travel-time correction gives the same limit along the ordinary and collar interventions. Finite arrival populations do change. [Velocity, source dependence and scope](velocity/README.md)

A second remote test measures the phase of interference between two native launch ports. The ordinary and collar interventions change that phase, but different frequency settings of the same receiving records exclude a common delay response even with a freely recalibrated affine clock. [Fringe data, occupied spectral lines and the clock test](fringes/README.md)

Eliminating the stationary heads also gives a source-dependent optical self-energy supported on four marked ports throughout the retarded frequency plane. A connected pilot detects both a charge-common delayed hopping and a charge-common loss of interference. Their signed responses disagree after ordinary calibration. [Connected retarded propagation and native pilot counts](influence/README.md)

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python3 -m matched_einstein_limit.source_implication.common_front.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python3 -m matched_einstein_limit.source_implication.common_front.certificate --write
OPENBLAS_NUM_THREADS=1 python3 -m pytest -q code/matched_einstein_limit/source_implication/common_front
```
