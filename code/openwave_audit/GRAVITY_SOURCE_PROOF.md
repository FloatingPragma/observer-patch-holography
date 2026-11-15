# Native clocks and a supplied gravitational lapse

For a graph with M undirected seams, uniform seam attempts and fair native
coins give

\[
\mathbb E[n'\mid n]=(I-L/(2M))n.
\]

The two coin currents on an oriented seam sum to the endpoint difference.
Summing the endpoint changes over seams proves this identity for arbitrary
integer loads. Accumulated current records do not occur in the native
transition rule. A field reconstructed from those records therefore has no
reciprocal effect under this scheduler.

For a tagged hard-core unit among K other units on the twelve-port carrier,
uniform preparation on labelled configurations gives

\[
P(\text{tag moves})=\frac{1}{12}\left(1-\frac K{11}\right),
\qquad P(\text{attempt touches tag})=\frac16.
\]

Each vertex has five incident seams among thirty. Conditional on a contact,
the neighbour is empty with probability (11−K)/11, and one of the two fair
coins moves the tag. The contact clock retains every hold, so its rate is
independent of K. A reduced hop rate is consequently not a universal dilation
of these two record clocks. The executable compares the actual native function
on 578 signed endpoint branches and independently enumerates the existing
development preparations K=0,1,2. It does not execute the reserved total-four
preparation.

## Stationary states on a pinned graph

The supplied pin intervention forbids motion of the background loads. The
tag's transition matrix on available ports is

\[
T=I-L_{\rm free}/(2M).
\]

It is symmetric. On each connected component, its stationary probability is
uniform. The weights assigned to different components are conserved and
depend on the initial preparation; there is no unique global stationary
state without that specification.

The executed flat(3,3) preparation has 216 ports, 23 pinned ports and free
component sizes 1, 5, 7 and 180. One free port is isolated. Thirteen free
ports lie outside the largest component. The stationary probability 1/180
belongs to a preparation supported on the largest component. Independent
adjacency traversal establishes these counts without a transition-matrix
eigensolver. Symmetry gives zero stationary current along every free bond.
This statement concerns the supplied pin intervention and does not establish
attraction or a gravitational equilibrium.

## Lapse and frequency

On the specified compact two-level subspace, the supplied coherent generator
is N times a matrix with eigenvalues ±E, where E²=m²+6. A port preparation
has

\[
p_{\rm port}(t)=1-\frac6{E^2}\sin^2(NEt).
\]

The sine-argument rate is NE; the population angular frequency is 2NE.
The ratio between two such clock frequencies is N_a/N_b. The executed
four-clock receipt at m=.4 agrees with that ratio. It is a comparison under
a supplied lapse and coherent source, with no self-starting clock mechanism
or attachment to the native attempt clock.

For a static metric with flat spatial sections,

\[
ds^2=-N(x)^2dt^2+\delta_{ij}dx^idx^j,
\quad R_{tt}=N\Delta N,
\quad R_{ij}=-N^{-1}\partial_i\partial_jN,
\]

so

\[
G_{tt}=0,\qquad
G_{ij}=N^{-1}(\delta_{ij}\Delta N-\partial_i\partial_jN).
\]

A Newtonian-looking harmonic lapse therefore does not supply a vacuum
Einstein solution. An independent linearized calculation uses
g_tt=−1+2μ/r and g_ij=(1+2γμ/r)δ_ij. Away from the origin,
G_ij=(1−γ)∂_i∂_j(μ/r). At (3,0,0), G_xx/μ=2(1−γ)/27. The γ=0 supplied
rate metric has nonzero exterior spatial stress; γ=1 cancels this linear
residual. This calculation tests the stated metric interpretation, without
asserting a physical spatial metric for the native source.

## A common dimensional attachment

For the registered geometry comparison, let u be its electron-anchored
energy unit, a its independently calibrated length, and g its declared
dimensionless Newton response. Its own formulas give

\[
G_{\rm model}=\frac{ga}{u\hbar c},\qquad
\frac{c_{\rm model}}c=\frac{ua}{\hbar c},\qquad
\frac{G_{\rm model}}{c_{\rm model}/c}=\frac g{u^2}.
\]

The retained inputs give G_model=79.2936249871 GeV⁻² and
c_model/c=.00195750508412. Fixing the light speed by changing a leaves
G_model/G_observed=6.03793713114×10⁴². Thus this length parameter cannot
repair both observables at fixed u and g. This excludes that dimensional
attachment, not every source or normalization in the repository.

`gravity_source_audit.py` reconstructs the native rates, component counts,
linearized tensor and common-unit invariant. `gravity_source_results.json`
retains all source hashes and checks their stability. The reproduction command
from the workspace root is:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 reverse-engineering-reality/code/openwave_audit/gravity_source_audit.py --workspace "$PWD" --output /tmp/gravity_source_results.json
```

**Scope**: finite native event clocks, a supplied pinned graph, compact
coherent lapse clocks and one dimensional comparison; the calculation
establishes no source-derived gravitational field, physical Newton strength
or OpenWave cell award.
