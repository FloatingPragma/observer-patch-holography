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

Independently of the radar, the reference head spectrum supplies the
entropy derivative $K=-\log X_*$. If $k$ consists of its 18 joint Pauli
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
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/correlation_rank/spatial_cuts
```
