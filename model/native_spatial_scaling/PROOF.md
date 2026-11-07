# Spatial displacement and local clocks from native records

## The selected source and the result

Use the original fair binary `port_pair` federation, independent uniform
seam selection, repair probability one and the native C/E instruments.
No endpoint weight, kinetic term, Hamiltonian, spatial dimension or clock
service is inserted. For the displacement result the prepared state is
the uniform diagonal **one-unit** sector. This stationary preparation is
declared; passive repair does not select its total or prepare it from every
input. The geometry is read from the simulator specification, which the
observer is allowed to inspect.

The spatial dictionary is the existing auxiliary affine cell-centroid
dictionary of [`uniform_green.py`](../physical_space_time/uniform_green.py).
The icosahedron's circumsphere has radius one. Every port of a cell receives
the same coordinate. These are numerical comparison coordinates, not
physical positions, a laboratory metric, or the floating vectors used to
route intercell seams.

Write s=4^ell. There are V=240s ports and M=630s seams. For the actual
one-unit path X_T at T global attempts, comparison coordinate x and r>0,

\[
 \mathbb E\|x(X_T)-x(X_0)\|^2\le {T\over12600\,16^\ell},\qquad
 \Pr(\|x(X_T)-x(X_0)\|\ge r)
 \le\min\left(1,{T\over12600r^2 16^\ell}\right).                 \tag{1}
\]

These are finite record statistics, without a Gaussian limit or
Poissonization. The all-sector native E-tree in
[`native_sector_readout`](../native_sector_readout/PROOF.md), appended
after the body, identifies its terminal occupation. Reversing the complete
body currents then identifies every preceding occupation. The detector's
own displacement, records, failures and host limits are retained. Initial
occupation disclosure is unnecessary for this decoder; the specified
initial probability law remains an input to (1).

## The graph attachment at every combinatorial level

The original source starts with the twenty distinct triangular faces of
the icosahedron. `_refine_level` caches a fresh midpoint identifier for
each unordered edge and replaces (a,b,c) by

```
(a,ab,ca), (b,bc,ab), (c,ca,bc), (ab,bc,ca).
```

On each side of an original edge the same two subdivided edge identifiers
are used. Interior child edges have exactly two incident children; boundary
child edges have one on each side. Distinct children are distinct triangles.
Induction therefore preserves a simple closed triangular complex and its
connected cell adjacency, independently of midpoint coordinate rounding.
Each face has three distinct adjacent faces: two distinct simplices sharing
two edges would be the same triangle. Thus the dual has 20s vertices and
30s edges.

The source keeps all thirty original icosahedron seams in every cell.
`assign_echosahedral_ports` assigns each of the three cell-neighbour
endpoints a distinct port. If its initial choices collide, its rectangular
three-by-twelve assignment uses distinct columns. The returned federation
also explicitly checks uniqueness of every used global port slot. These
are combinatorial postconditions, independent of the numerical choice of
which three ports to use. Hence there are 600s internal and 30s crossing
seams. Exactly 60s ports have degree six; the other 180s have degree five.
Connectivity follows from connected internal carriers and the connected
dual.

For any whole cell, its three external halo ports belong to three distinct
neighbours. Each halo port has one crossing seam back to the root cell,
five internal seams leading outside the fifteen-port window and no seam
to either other halo port. It has no second crossing seam because routing
is injective at each endpoint. This proves the geometry premise of the
four-state halo calculation in
[`native_local_response`](../native_local_response/PROOF.md) at every
combinatorial level and for every such routing. Its complete-channel bounds
continue to require the stated record roots, independent uniform selection
and repair probability one; observing crossing records enlarges the roots.

The code checks the actual routed graphs and every cell halo at levels
zero through three. Its all-level conclusion follows from the induction
and source postconditions, not extrapolation of those finite checks.
It does not assert that a floating-point constructor returns at arbitrary
levels or that sufficient memory exists to allocate them. Native admission
accepts a formal federation level, while its expanded evaluator supports
only levels zero through two and retains its existing step/branch limits.

## Exact all-level comparison energy

Put q=sqrt(5). The raw base vertices in `uniform_green.py` are
permutations of (0,+/-1,+/-(1+q)/2). Their squared circumradius is
R0²=(5+q)/2, and the raw side has squared length four. Division by R0
therefore gives side squared

\[
 a^2={4\over R0^2}=2-{2q\over5}.
\]

For adjacent original faces, their raw centroids differ by one third of
the difference of their two nonshared vertices. The literal thirty pairs
give squared difference (6+2q)/9 before radius normalization. For example,
faces (0,11,5) and (0,5,1) differ through vertices11 and1; substitution of
the listed coordinates gives exactly 6+2q. The finite base checks verify
the same value for every pair. Thus the normalized base crossing energy is

\[
 E_0=30{6+2q\over9R0^2}={20+4q\over3}.                     \tag{2}
\]

Let n=2^ell. An original triangular face has n² small triangles and 3n
boundary edges. Counting its triangle-edge incidences shows it has
3(n²-n)/2 internal dual edges. Across all original faces there are
30(n²-n) such edges and 30n original-face hinge edges. An interior centroid
pair has squared distance a²/(3n²). A hinge pair's nonshared vertices
differ by 1/n times the base difference, so its squared distance is
E0/(30n²). This yields the exact source comparison energy

\[
 E_\ell:=\sum_{\{u,v\}\in\mathcal E}\|x(u)-x(v)\|^2
 =10a^2+{E_0-10a^2\over n}
 =20-4q+{-40+16q\over3\,2^\ell}<12.                        \tag{3}
\]

Internal carrier edges contribute zero. The coefficient of 1/n is negative
because q<5/2, and 20-4q<12 because q>2. No floating routing affects (3):
every routed port in a cell has that cell's coordinate and every dual edge
occurs exactly once. The implementation evaluates (2)–(3) in Q(sqrt5)
and compares them with the literal edge sums at levels zero through three.

## From complete native attempts to a finite displacement bound

In total number one, a selected seam moves an occupied endpoint across
with probability one half and otherwise holds it; an unoccupied seam
holds. This is a statement about the sum of the original raw outcomes,
including zero addresses, and is the existing exact source identity

\[
 P=I-{L\over2M}.
\]

Each fixed-seam operator is (I+Swap)/2, so P is symmetric, bistochastic
and positive. Its uniform law is stationary. For any real port function f,
stationarity gives

\[
 \mathbb E(f(X_T)-f(X_0))^2
 ={2\over V}\langle f,(I-P^T)f\rangle
 \le {2T\over V}\langle f,(I-P)f\rangle
 ={T\over MV}\sum_{\{u,v\}}(f(u)-f(v))^2.                  \tag{4}
\]

The inequality follows eigenvalue by eigenvalue from
1-p^T<=T(1-p), 0<=p<=1. Sum (4) over the three coordinate components
and use (3), MV=151200*16^ell, and Markov's inequality to obtain (1).
Equivalently E_l/(MV) is the exact stationary mean one-attempt squared
coordinate increment. Correlated increments are not treated as independent.

At the count b_l=1260*8^ell and T=floor(t*b_l), (1) is at most
t/(10*2^ell) for the mean square. The displacement probability is at most
t/(10*r²*2^ell). This finite one-unit statistic freezes at that proposed
ballistic count scale. It complements the existing density-field path
result in [`observer_time_fields/REFINEMENT.md`](../observer_time_fields/REFINEMENT.md);
it neither proves a new continuum limit nor excludes a different sector
or observation dictionary. In particular it does not call one source unit
a measured particle.

## A stopped local clock without an invalid time substitution

Fix a nonempty seam set I used for clock ticks, of size d, independently
of source state and all records. Let a=d/M. A receiver's incident seams
give d=5 or6; a whole cell's internal seams give d=30. Let U and V be
the normalized remote and incident one-unit transitions. They are
symmetric bistochastic averages of the same seam maps. For d<M the
completed one-tick transition is

\[
 K=VC,\qquad C=a\sum_{j\ge0}(1-a)^jU^j.                    \tag{5}
\]

This sums every remote waiting word followed by one incident attempt.
Elapsed global count is a retained input-independent clock record; the
formula is its marginal, not conditioning on a duration. The full instrument
keeps all raw waiting/incident records. For d=M use C=I.

Let epsilon_U and epsilon_V denote their stationary mean-square coordinate
increments. Equation(4) for U gives the C increment at most
(1-a)epsilon_U/a. The squared triangle inequality at the intermediate
state yields

\[
 \epsilon_K\le2\left({1-a\over a}\epsilon_U+\epsilon_V\right)
 ={2\epsilon_P\over a}.
\]

K is bistochastic even when U and V do not commute. At k ticks the pathwise
inequality ||sum_j Delta x_j||²<=k sum_j||Delta x_j||² and stationarity
give

\[
 \mathbb E\|x(X_{\tau_k})-x(X_0)\|^2
 \le{2k^2\epsilon_P\over a}
 \le{k^2\over10d\,4^\ell}.                                \tag{6}
\]

This proof does not substitute E tau_k=k/a into (1): the stopping time
and source path are correlated. With a finite global cap H, the timeout
probability is Pr(Binomial(H,a)<k). Report timeout as its own outcome and
assign no successful-event score to it. Completion-event probabilities
remain bounded by (6), without normalization; comparing capped and completed
record laws incurs at most this timeout mass. General adaptive FIFO
selectors, state-dependent ticks and reading unretained external controller
information require their actual dependencies and are outside (1),(5),(6).

## First arrival and the support of the local clock

For any distinct sender u and receiver v choose a shortest native port path
u=x0,...,xD=v. Its first D-1 edges avoid v. Choose those seams in that order
and their native coins placing the unit at the next endpoint. With pure E
readout this is one complete positive raw word of probability (2M)^(-D).
All its earlier attempts are remote to v, and the last is its first
incident tick. The receiver's incoming contrast and current record that
arrival. Hence first-tick arrival probability is at least (2M)^(-D),
not zero at distances D>1. Other words and timeouts are retained; this
word is a lower bound, not an assumed unique history.

The port diameter is unbounded in the mathematical graph family. The
opposite base-face centroids6 and8 have z-coordinate separation greater
than4/3 after division by R0, verified by exact squaring. Their descendants'
centroid averages preserve those values, so some descendants have at least
that separation. Each dual-edge displacement has length at most
2/(sqrt3*2^ell), by the existing affine hinge inequality; internal edges
have zero displacement. The triangle inequality forces their port distance
to exceed2^ell. Thus a completed incident tick has no capacity-independent
strict port-distance support bound. Its far-arrival probability can be
extremely small; this support statement is consistent with the stationary
mean-square bound (6). It supplies no violation of locality in global
attempt count and no relativistic causal interpretation.

Concrete native paths at levels zero through three are checked against the
original endpoint rule. Formal arbitrary-level paths and completed waiting
laws are mathematical continuations. A single admitted program retains its
attempt ceiling; the expanded interpreter retains its geometry, step and
branch ceilings. The package does not claim unbounded actual execution.

## Record resources and reuse boundaries

The combined body/detector program uses the actual pure E uniform body,
the deterministic native spanning tree, FIFO capacity one and oldest-record
export. Every selected raw word is checked with the original `_outcomes`
kernel, with its unnormalized state and complete archive. `decode` reads
the detector records and reverses the body currents, never consulting the
initial working mask. The verifier's known basis input is a control, not
an input to that decoder. The detector changes the final source; its
post-detector state is kept separately from the inferred body snapshots.
The record dictionary does not replace the complete quantum instrument
by a classical one.

On the level-zero241-command example the original expanded evaluator
returns `RESOURCE_STOP` after one body attempt. The complete raw-word query
is a separate exact factorization of its admitted kernel; it does not
relabel that host stop. The full raw alphabets, including impossible words,
remain part of the instrument. Access to every body record plus the final
tree differs from a receiver's incident tape alone.

Existing [`federated_observer/SCALING.md`](../federated_observer/SCALING.md)
supplies density closure, full-sector positive transfer, logarithm bounds
and Poissonized propagation. Existing physical-space-time and observer-time
field proofs supply the controlled subsequential diffusive limit and its
two-dimensional spectral result on their dictionary. These results are reused under their stated hypotheses. The additions are the all-level halo attachment, exact
crossing-energy formula, finite native one-unit displacement and stopped-clock
bounds, and complete body-path readback. They do not select physical density,
spacetime, a common interacting action or any measured particle spectrum.
