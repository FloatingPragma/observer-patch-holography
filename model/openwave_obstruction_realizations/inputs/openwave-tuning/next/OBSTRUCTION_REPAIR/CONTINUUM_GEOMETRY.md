# Continuous geometry and the remaining gravity construction

Gravity need not be formulated on a discrete spatial lattice. In general
relativity the metric is a continuous field, and proper distances can vary
continuously. A numerical mesh, a discrete event program and a continuous
effective metric are three different objects.

The [repository construction contract](../../../reverse-engineering-reality/AGENTS.md)
requires executable events and observer comparisons. It explicitly permits
smooth continuum limits with stated domains and error bounds. Therefore a
continuum metric may be the effective description of finite records; no
requirement here makes spatial contraction occur by removing whole cells.
A particular microscopic geometry still needs an implementation and a
controlled relation between those records and the metric.

For example, a spatial metric g_ij(x)=exp(2 sigma(x)) delta_ij gives a curve
the proper length integral exp(sigma(x)) |dx|. Varying sigma continuously
varies that length continuously. A fixed-connectivity numerical mesh can
represent this through its metric-dependent weights. A model with every
edge weight and regional algebra fixed is a more restrictive construction;
its failure to respond is not a theorem about every continuum metric.

There is also a distinction between physical curvature and a coordinate
rescaling. For spatially constant sigma on infinite flat space, adjusting
the coordinate radius to keep a ball's proper radius fixed describes the
same geometry. Its entropy should be unchanged. The exact identity in
`gravity/continuum_source/PROOF.md` establishes that invariance for the
declared source. A nonconstant sigma can produce curvature and a genuine
fixed-proper-volume area response; `GRAVITY_KERNEL_PROOF.md` computes the
corresponding smooth geometric variation. Neither construction removes
cells to obtain it.

The continuum source H_g=D_g+beta R^3 D_g^4 has a proper length R and
metric-dependent Dirac operator D_g. Its regional entropy is finite and
its spatial coordinate covariance is proved. The existing finite lattice
calculations approximate that kind of continuous coframe response.
The original time generator fails a separate temporal-observer test.
A distinct flattened generator now has a complete Poincare representation,
and its boost-orbit CAR apparatus preserves the finite entropy exactly.
A continuum operator bound controls its nonlocal tails at large geometric
separation. These are observer-component results; continuous geometry by
itself does not supply curved observer transport or matter-selected gravity.

The Einstein argument requires an area response, the appropriate matter
entropy/stress relation, and a stationarity condition for compatible
regions at every point and for all timelike observers. It also needs one
physical clock and normalization. These are explicit hypotheses in
[Jacobson's entanglement-equilibrium derivation](https://arxiv.org/html/1505.04753v4),
particularly equations 23–27. Writing down a smooth metric or solving an
assumed Einstein equation does not establish that this OPH source produces
those hypotheses. A successful construction must derive the response to
matter rather than prescribe the desired gravitational field.

Accordingly, this lane is pursuing both continuum source proofs and finite
executable approximants. The fixed-lattice negative results remain in the
failure register with their restricted scopes. They neither require
discrete contraction nor prove a general impossibility of OPH gravity.

There is now an explicit observer-refinement theorem in
`cell_entropy_limit/PROOF.md`: nested scalar cell observations, dense in
a fixed region and including all spin components, recover that region's
actual continuum entropy. The proof uses the declared source's localized
Schatten estimate and includes the unobserved subcell modes. It does not
assume that coarse-cell entropy is an area law, and entropy along the
refinement need not be monotone. Rates for a particular mesh, shape
derivatives and the large-region limit remain separate calculations.
Thus a finite observer approximation has a defined continuum target;
discrete cell deletion is not a necessary gravity mechanism.
