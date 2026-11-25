# Reading the source from spatially resolved material flux

At the calibrated 80-cell source, the four-tick material experiment can
recover every change in two-head correlations that preserves the two
marginals and weight within each parity. The observer records ten oriented
cell-to-cell momentum transfers as well as the eight existing independent
regional translation means. These eighteen means distinguish every
direction in the eighteen-dimensional correlation space. This includes
all nine directions invisible to the regional means and the radar volume
taken together. [Regional experiment](../../../directional_stress/translation/README.md) · [Equal-volume comparison](../README.md)

The translation chart, vacancy and excess launches, read time, two
contrast parities, and detector operations are those of the regional
experiment. Its operator identity expresses momentum transfer as a sum
of directed edge terms. Opposite internal contributions cancel in the
regional sum. Here the observer retains the reading at each interior
or exterior cut. There are 25 oriented edges after discarding the
redundant reverse orientations. The two mobile species give 50 distinct
means, each with terminal misses carrying zero. The ten selected
species/edge pairs were chosen using just the material responses. Each
measures interference between native ports with the usual outward sign.
[Oriented operators and interval propagation](measure.py) · [All 50 measured effects](effects.json)

## Rank, source and precision

Use the nine orthonormal joint Pauli products in each parity as source
coordinates. Let $B$ map the 18 coordinates to the eight regional
means and ten selected spatial means. The fixed-point flight supplies
rational interval entries for the true material response. For each edge,
the measured operator has norm at most one and uses at most four ports.
Every selected four-by-four propagation block has operator error below
$10^{-19}$. Stacking the four blocks gives error at most $2\times10^{-19}$.
Compression against the propagated state, followed by contraction with
a normalized Pauli product of trace norm two, gives entry error below
$10^{-18}$. The original regional channels retain their published error
bounds. A proposed right inverse $L$ has $\|I-BL\|_F<0.03636$.
Thus $B$ is invertible. Its smallest singular value lies between
$7.4211\times10^{-13}$ and $8.3337\times10^{-13}$ in the recorded native
units, and its condition number exceeds $4.14\times10^7$.

Head tomography supplies the reference density matrix, including its
eigenprojectors and eigenvalues, and therefore the entropy derivative $K=-\log X_*$. If $k$ consists of its 18 joint Pauli
coefficients, $w=kB^{-1}$ is the unique row giving $dS=w\,dm$ on this
source slice. The published decimal $w$ uses the calibrated preparation’s
density matrix and the material effects. It uses no area measurement,
volume reading or entropy-to-area gain. At the true equilibrium, the
read weights have worst-case bias below $1.587\times10^{-5}$ nats per
unit Hilbert–Schmidt variation. The bound includes the rational spectral
calculation, rounding the weights, the material effect uncertainties, and
the state-localization error in the matrix logarithm. [Reader weights and conditioning](results.json)

The sum of absolute weights is below $2.268\times10^9$. An absolute
uncertainty $\epsilon$ in each material change can therefore add up to
$2.268\times10^9\epsilon$ nats to the inferred change. Limiting this
contribution to $10^{-3}$ nats, for example, takes individual absolute
errors below $4.4\times10^{-13}$ in the worst case. The held check uses
computed expectation values. A separate sparse-exponential flight
compares two reachable preparations displaced by $\pm5\times10^{-5}$
along the previously published central correlation direction invisible
to the regional means. It reads a signed entropy slope of 0.0015940293
from the material channels. Direct differencing of the two state entropies
gives 0.0015940511. These numerical arms use the published central
direction, while the rank and response bounds cover the exact equilibrium.
[Held preparations and independent flight](held.py)

## The complete signed response on one source

On the 80-cell source the same material effects, augmented with thirteen
local coordinates that each commute with the parity measurement, determine the entropy derivative for
**all 31 trace-one variations**. Use the six one-head Pauli means in each
parity, in the order $XI,YI,ZI,IX,IY,IZ$, each divided by two. A mean
in one parity is a joint expectation with that parity's projector, so it
includes its preparation weight. The thirteenth mean is the parity
contrast $(I_4,-I_4)/\sqrt8$. Each setting can use a fresh copy of the
preparation. [Head and parity analyzers](../../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/README.md)

Write $dp$ for the 18 joint Pauli coordinates and $dm$ for these 13 local
coordinates. They are an orthonormal basis of trace-zero block variations.
The same material Hamiltonians and edge detectors give
$dJ=B\,dp+C\,dm$. The fixed-point flight encloses each entry of the
$18\times13$ matrix $C$ independently. Let $k_p,k_m$ be the corresponding
coordinates of the independently reconstructed $-\log X_*$. Define
$$
 w=k_pB^{-1},\qquad v=k_m-wC .
$$
Substitution gives $w\,dJ+v\,dm=k_p\,dp+k_m\,dm=dS$ for every trace-one
variation. The reference is the calibrated maximum of $S+\eta A$ at
the measured cap-six volume, and $dV$ has rank one there. Therefore
$$
 w\,dJ+v\,dm=dS=-\eta\,dA
 \qquad\text{on all 30 directions satisfying }dV=0.
$$
All quantities refer to fresh instances of the same two-head preparation.
The material and local weights use the reference density matrix and the
material operators. The ordinary exchange gain and six-unit radar enter the
last equality through the finite stationarity theorem.
[Fixed-volume theorem](../../PROOF.md) · [Certified material effects](local_effects.json)

The archived decimal weights have worst-case error below
$1.783\times10^{-5}$ nats per unit Hilbert–Schmidt variation on the entire
trace-one space, compared with $1.587\times10^{-5}$ on the joint sector.
The full worst-case read-noise multiplier is below $2{,}267{,}315{,}366$.
The bound includes localization at the true constrained maximum. A
symmetric two-arm check changes the left head in the even sector, the
right head in the odd sector, and the parity probabilities. It uses the
independent radar to match both arms to the same survey volume. The check
compares the weighted four-tick material and head measurements to the
direct spectral entropy and the ordinary-calibrated area difference.
Both arm volumes are 2.763471836788309 to the reported numerical
precision. The material derivative of 17.6617753854 and local derivative
of −17.3273213219 sum to 0.3344540635 nats per unit tangent. Independently,
the spectrum gives 0.3344540613 and the calibrated area gives 0.3344540641.
[All coefficients and held results](full_source.json) · [Native reconstruction](full_source.py)

## A fixed spatial weighting

The chart was registered with source cell 0. Before reading the response,
take $R=1+\max_{u\in P}d(0,u)=6$ graph steps for the chart $P$. Give its
cells the weight $f(u)=(R^2-d(0,u)^2)/(2R)$ and exterior cells weight zero.
Sum the local outward momentum transfer at each cell with this weight.
Equivalently, the multiplier at the oriented edge $u\to v$ is $f(u)-f(v)$.
Keep the vacancy and excess readings separate. This uses the quadratic
profile of a ball with distance measured on the native cell graph.

Even when both weighted aperture readings, all eight regional channels and
the surveyed volume are retained, their differential has rank 11 on the
18 correlation coordinates. An exact tangent in the remaining seven
dimensions has entropy slope between 0.00072217 and 0.00073132 nats per
unit Hilbert–Schmidt norm. The interval right-inverse calculation projects
the area gradient off all eleven rows and applies the independently fixed
ordinary gain. Full row rank and the interior preparation neighborhood
integrate the tangent to a local curve on which the eleven readings are
exactly constant. [Aperture certificate](results.json)

Material tomography recovers the signed source using weights calibrated
to this particular head spectrum. A physical null-balance theorem would
also have to select the stress weighting and normalization from the
material dynamics on a common family of proper balls. This particular
graph-distance weighting loses information needed for that identification.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.full_source --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/correlation_rank/spatial_cuts
```
