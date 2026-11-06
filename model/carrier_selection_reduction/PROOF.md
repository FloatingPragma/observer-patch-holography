# Minimal quadratic readout and stable carrier incidence

This is an alternative conditional carrier selector. Its input is an abstract
real three-dimensional response space with an inner product. No icosahedral
incidence, spherical triangulation, group action, vertex degree, Euler
characteristic, equal channel weights or port count is assumed. The native
rank-three example in the existing source stack uses a declared icosahedral
carrier; it therefore cannot independently discharge this theorem's
three-dimensional input.

## Operational premises

An unknown response is any real symmetric operator `X` on this space. A
directional channel with unit vector `u_i` reads `tr(P_i X)`, where
`P_i=u_i u_i^T`. Channels have positive resource weights `w_i`, with
`sum_i w_i=1`. Opposite directions give the same quadratic channel, and
duplicate channels are merged by adding their weights. The weighted readout
map and its information operator are

```
A(X)_i=sqrt(w_i) tr(P_i X),
M=A* A=sum_i w_i |P_i><P_i|.
```

The norm on symmetric operators is Frobenius norm. Thus
`||A(X)||^2=<X,M X>` and the worst-direction information is `lambda_min(M)`.
This defines the criterion without a probability model. Independent readout
errors whose precision is proportional to resource weight give its usual
least-squares interpretation, but that noise law is not derived here.

The selector first minimizes the number of distinct unoriented channels
subject to full quadratic informational completeness (`M>0`), then maximizes
`lambda_min(M)` at that count. The minimum and the optimum are attained by
the explicit construction below. A separate inverse-port convention exposes
both signs of each selected channel as distinct ports. Without that
convention the readout has six channels, not twelve oriented ports.

## Sharp weighted theorem

Let `n` be the number of positive-weight channels. Then:

1. Informational completeness requires `n>=6`.
2. For every number of channels and every choice of weights,
   `lambda_min(M)<=2/15`.
3. At `n=6`, equality holds exactly for equal weights `w_i=1/6` and the six
   unoriented icosahedral axes, up to an orthogonal change of frame,
   channel permutation and reversal of channel representatives.

Consequently the selector gives six axes and, under the inverse-port
convention, twelve icosahedral directions. Joining two distinct directions
when their scalar product is positive gives degree five and thirty edges;
the three-cliques give twenty faces, with two faces per edge and a five-cycle
link at every port. These counts and incidence properties are conclusions,
not premises. A physical interpretation of this graph as native seams is a
separate attachment.

### Universal bound and equality

The space `Sym(3)` has dimension six, while each term in `M` has rank one,
proving the first assertion. Write

```
e0=I/sqrt(3),    Q_i=P_i-I/3,    gamma=2/15.
```

Every channel has `<e0,P_i>=1/sqrt(3)` and `||Q_i||^2=2/3`.
The compression of `M` to the five-dimensional trace-free space therefore
has trace `2/3`. Its smallest eigenvalue is at most `2/15`; the smallest
eigenvalue of the full operator is no larger.

If the full operator attains this bound, its trace-free compression is
exactly `gamma I`. In the block matrix of `M-gamma I`, that compression is
zero. Positivity forces the off-diagonal block to vanish: a nonzero cross
term would make the quadratic form negative on some vector `(t,x)` for
sufficiently small `t` of the opposite sign. Hence

```
M=M0,    M0(X)=(2 X+tr(X) I)/15.
```

In particular, both weighted second-moment and quadrupole isotropy follow
from the one information objective; they are not separate assumptions.

At `n=6`, let `F` have columns `sqrt(w_i) P_i` in any orthonormal basis of
`Sym(3)`. It is square and invertible, so

```
F^T M^(-1) F=I_6.
```

The eigenvalues of `M0^(-1)` are `3` on the scalar line and `15/2` on its
orthogonal complement. Thus `<P_i,M0^(-1)P_i>=1+5=6`, and the diagonal
identities give `w_i=1/6`. For `i!=j`, the off-diagonal identities give

```
0=<P_i,M0^(-1)P_j>=(15/2)(u_i.u_j)^2-3/2,
```

so every squared overlap is `1/5`. Conversely the explicit axes with these
weights have `M=M0`, which proves attainability.

### Geometry and incidence without a carrier template

Let `G_ij=u_i.u_j` and write `G=I+S/sqrt(5)`. Then `S` is symmetric,
has zero diagonal and off-diagonal entries in `{+1,-1}`. The weighted
second-moment identity gives `sum_i P_i=2I`, hence `G^2=2G` and `S^2=5I`.
Switch signs of the representatives so `S_0i=1` for all `i>0`.
The entries `(S^2)_0i=0` say that each of the remaining five vertices has
exactly two negative edges among its other four vertices. Their negative
graph is a simple two-regular graph on five vertices, necessarily one
five-cycle. This determines the switching class uniquely.

Put `phi=(1+sqrt(5))/2`. An attaining set of representatives is

```
(0,1,phi), (0,-1,phi), (1,phi,0), (-1,phi,0),
(phi,0,1), (phi,0,-1),
```

each divided by `sqrt(phi+2)`. Its Gram matrix has the switching class just
proved. Two spanning vector families with the same Gram matrix differ by
an orthogonal map: the map sending one family to the other is well-defined
on linear combinations because the Gram matrix gives the squared norm of
each combination, and it preserves inner products. Including opposite
representatives gives the regular icosahedron.

The exact certificate independently enumerates the `2^10` normalized sign
matrices. Exactly twelve satisfy `S^2=5I`, and they are exactly the twelve
labelled five-cycles. For the induced twelve-direction graph it verifies
all degrees, edges, triangular faces, two-face edge incidences and links.
Its Euler value two is computed from these derived counts.

## Stability of the discrete incidence

Exact metric symmetry need not be demanded to recover this combinatorics.
With six channels and the same positive normalization, assume

```
lambda_min(M)>=2/15-10^(-7).
```

Then every weight lies in `[50/301,50/299]`, every squared overlap differs
from `1/5` by at most `1/375`, and the positive-overlap graph on the paired
directions has exactly the icosahedral incidence above. This does **not**
assert exact icosahedral angles, equal weights, an isometric `A5` action, or
an approximate continuum or physical gauge theorem.

Here are conservative bounds; no claim of a sharp tolerance is made. Write
the block matrix of `M` on the scalar and trace-free spaces as

```
M = [ 1/3  b^T ; b  gamma I+Delta ],    epsilon=10^(-7).
```

The lower eigenvalue bound and fixed trace imply
`||Delta||<=4 epsilon`. Positivity of `M-(gamma-epsilon)I` gives
`||b||^2<=(1/5+epsilon)5 epsilon=epsilon+5 epsilon^2`.
Thus `||b||<1/3000` and

```
||M-M0|| <= 1/3000+4 epsilon,
||M^(-1)-M0^(-1)||
 <= (1/3000+4 epsilon)/(gamma(gamma-epsilon)) < 1/50.
```

These are exact rational inequalities except for the displayed norm bound,
whose square is bounded rationally. Applying `F^T M^(-1)F=I` as above yields
`|1/w_i-6|<1/50` and `|(u_i.u_j)^2-1/5|<1/375`.

Let `T=sum_i P_i`. Since `sum_i w_i Q_i=sqrt(3)b`,

```
||T-2I|| <= 6 sqrt(3)||b|| + sum_i |1-6w_i|
           < 7/2000+6/299 < 3/125 = tau.
```

All off-diagonal overlaps are nonzero. Set `S_ij=sign(u_i.u_j)` and
`A=I+S/sqrt(5)`. For `i!=j`,

```
|A_ij-G_ij|
 = |(u_i.u_j)^2-1/5|/(|u_i.u_j|+1/sqrt(5))
 < sqrt(5)/375 < 3/500.
```

There are five such entries per row, so `||A-G||<3/100=sigma`. Also
`||G||=||T||<=2+tau` and `||G^2-2G||<=(2+tau)tau`. Consequently

```
||S^2-5I|| = 5||A^2-2A||
 <= 5[(2+tau)tau+sigma(6+2tau+sigma)]
 = 57729/50000 < 2.
```

Every off-diagonal entry of `S^2` is an even integer (a sum of four signs),
and every diagonal entry is five. An absolute bound strictly below two
therefore forces `S^2=5I` exactly. The five-cycle argument identifies the
same discrete graph. This integer gap is the step that makes a finite
information tolerance sufficient for exact combinatorial incidence.

The perturbed witness rotates one axis by a rational Cayley rotation with
parameter `t=10^(-8)`, keeping the weights equal. Its information operator
changes in norm by at most `4t/3<10^(-7)`: a unit vector changes by at most
`2t`, its rank-one projector by at most twice that in Frobenius norm, and
the one weighted information term changes by at most one third of the
projector change. Its metric is not equiangular, while its discrete graph
satisfies the stable conclusion. This is a finite exact bound, not a
floating eigenvalue diagnostic.

## Controls and comparison with the existing selector

The tetrahedral four-channel and octahedral three-channel readouts are
rank deficient on all six symmetric response coordinates. Euler curvature
alone places no restriction on their information rank. A different six-axis
configuration is informationally complete but fails the optimal bound, so
channel minimality alone does not select geometry. Unequal icosahedral
weights likewise fail the optimum. A union of two rotated optimal frames
is optimal with twelve axes and twenty-four paired ports, so optimal
information alone does not select minimal capacity. The six unpaired
representatives give the same quadratic readout without twelve ports.
The analogous two-dimensional response has three minimal channels and six
paired directions, making the supplied response dimension load-bearing.

The variational sieve in
`paper/screen_microphysics_and_observer_synchronization.tex` assumes a
spherical triangulation, integer defect charges, a strict additive defect
cost, an inverse-port convention, feasible isotropic six-axis data, and a
vector/quadrupole determinant objective. This alternative replaces the
topology/defect-cost route, the preassigned channel count, equal channel
weights and the two-factor objective by full quadratic readout, minimum
channel count and optimal worst-direction information. Exact icosahedral
uniqueness itself is shared finite mathematics, not claimed as a new result.
The weighted reduction and the stated stable-incidence implication are the
additional results supplied here.

Neither minimum resource selection nor the information objective follows
from arbitrary record consensus or deterministic execution. This theorem
does not derive the response dimension, the quadratic readout law, inverse
ports, physical seams or their native selection. The marked-sphere degree
floor in `code/a5_closure/continuous_equivariant_degree_proof.md` concerns
the topological degree of equivariant maps between two already specified
`A5` sphere actions. It is distinct from vertex valence and supplies none of
these operational premises. No frozen source or prediction is changed.
