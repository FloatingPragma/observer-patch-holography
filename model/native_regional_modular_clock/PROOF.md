# Regional modular clocks of the native transfer reference

This construction uses the actual fair binary exclusion transfer on the
source graph. It attaches regional modular data to the same reference
`rho=P/Tr(P)` used by the existing native phase preparation and energy
measurement on the 66-state cell. Its growing-family reference is a
specified mathematical state. The existing 66-state preparation is not a
scalable preparation of that global state, and no such producer is inferred.
All regions keep every supported local occupation sector.

## Exact restriction, without enumerating the global carrier

Let the connected graph have V vertices and M edges, with binary occupation
q. On its q-subset Hilbert space put

```
G = sum_edges (I - Swap_e),   P = I - G/(2M),
u = q(V-q)/(V(V-1)),
Z = Tr(P) = binom(V,q)(1-u),  rho = P/Z.
```

These are the literal fair endpoint population counts. The existing
federated-observer triangle argument proves `10/21 <= P <= I` on every
sector of the growing fair federation. Thus rho is faithful on its fixed-q
carrier. H=`-2M log(P)` has K=`-log(rho)=H/(2M)+log(Z)I`.
This identifies mathematical generators; it is not a conversion from
coherent service pulses to physical or native-repair time.

Partition vertices into A and B, with v=|A|, b=|B|>=2, internal edge counts
m_A,m_B and crossing count m_C. In a supported local block k write
l=q-k, d_l=binom(b,l), r=l/b and

```
c_B = 2 m_B l(b-l)/(b(b-1)),
N_boundary = sum_(i in A) crossing_degree(i) n_i,
beta = c_B + m_C r,
D_k = G_(A,k) + (1-2r) N_boundary,
a_k = 1 - beta/(2M).
```

The actual unnormalized regional density is exactly

```
R_(A,k) = (d_l/Z) [a_k I - D_k/(2M)].                 (1)
```

Proof: internal-A exchanges survive the exterior trace d_l times. Internal-B
exchanges contribute the scalar trace of their exclusion Laplacian. An
exterior vertex is occupied with probability l/b under its uniform l-sector,
so a crossing edge contributes diagonal `r+(1-2r)n_i`; its off-diagonal
entries change the local number and disappear under partial trace. Summing
these terms proves (1), including its normalization. Empty or impossible
blocks have zero mass and are retained; no local number postselection is
hidden. The formula needs only local adjacency, V,M,q and small binomial
ratios, not a matrix of dimension binom(V,q).

## Controlled rescaled regional modular generator

The algebra is the ordinary number-preserving regional algebra. Blockwise
scalar shifts generate no action on it. Therefore its modular action at
parameter `s=2M t` is generated exactly, up to such scalars, by

```
L_(A,k) = -2M log[I - D_k/(2M a_k)].                  (2)
```

The sign convention is explicit: `rho_A^(i 2M t) O rho_A^(-i 2M t)`
equals `exp(-i t L_A) O exp(i t L_A)`. For the Heisenberg convention
`alpha_H^t(O)=exp(i t H) O exp(-i t H)`, this is time `-t`.
The flow-error bounds use |t| and hold with either consistent orientation.

Let L=2m_A+m_C. Since `||G_A||<=2m_A` and
`||(1-2r)N_boundary||<=m_C`, `||D_k||<=L`. Also beta<=m_B+m_C<=M,
so a_k>=1/2. When z=L/(2M a_k)<1, the convergent logarithm series gives

```
||L_(A,k) - D_k/a_k||
    <= L^2 / [4M a_k^2(1-z)].                       (3)
```

This bound is valid for negative as well as positive eigenvalues of D_k.
It follows by bounding each n>=2 series term in operator norm and using
1/n<=1/2. No commuting assumption between source edge terms is used.

For a declared limiting density p in [0,1], set

```
a_p = 1-p(1-p),
D_p = G_A + (1-2p)N_boundary.
```

The additional error is at most

```
2m_C |r-p|/a_k + L |a_k-a_p|/(a_k a_p).              (4)
```

For fixed A in a bounded-degree growing graph and q/V -> p, r->p and
a_k->a_p. Explicitly,

```
|r-p| <= (V/b)|q/V-p| + v/b,
|a_k-a_p| <= |r-p| + (m_A+3m_C)/(4M) + 1/[4(b-1)].
```

Equations (3)-(4) prove the regional modular/exchange limit in operator norm.
For an observable O in a supported regional block, its flow error through
|t| is at most `2|t| ||O||` times that generator bound. Maximizing over all
supported k gives the complete local algebra bound. At half density the
boundary potential vanishes in the limit and the limiting generator is
`(4/3)G_A`. This is a rescaled high-temperature exchange limit, not a vacuum
ball modular boost. Unscaled local states approach the tracial reference;
the growing factor 2M is essential and must be charged in any implementation.

## Fixed-region modular overlaps do not become compatible

This local limit does not produce a compatible modular flow on fixed nested
regions. On a native triangle with A={0,1}, B={0,1,5}, choose O=n_0. At
half density the two limiting generators differ, in their action on O,
by the added crossing exchange on (0,5). Its commutator norm is one.
Consequently

```
||[(4/3)G_B, n_0] - [(4/3)G_A, n_0]|| = 4/3.        (5)
```

The (1,5) exchange commutes with n_0. Equation (5) follows on the
three-site one-occupation block and is an exact, capacity-independent lower
bound. Boundary potentials also commute with n_0, so omitting them cannot
repair this example. Compatible interior dynamics requires a separate collar
limit and a locality estimate with the observer support held away from both
boundaries. A claim of fixed-region modular covariance fails (5).

The full common state also fails the stronger state-preserving conditional
expectation premise. For any nonempty crossing set,

```
||[rho,N_A]||_HS^2
 = 2m_C binom(V-2,q-1)/(2M Z)^2 > 0                 (6)
```

when 0<q<V. Every crossing move changes N_A by one, and each crossing edge
has binom(V-2,q-1) undirected allowed occupation pairs; no entries are counted
twice. A unital CP retraction onto the number-preserving regional algebra
fixes its central local-number projections. Their multiplicative-domain
identity kills cross-block operators under that map. If it preserved rho,
all cross-block entries of rho would vanish, contradicting (6). Ordinary
restriction remains valid. The 66-state first-seam example gives exactly
`1/70560`. It cannot inherit the sigma_tree reference's conditional
expectations or compatible intrinsic flows.

## Geometry cannot be inferred from support response alone

For any fixed complete finite native history instrument {K_w}, a faithful
reference rho and a perturbation obeying
`(1-eta)rho <= rho' <= (1+eta)rho`, with eta<1, positivity gives

```
(1-eta) Pr_rho(E) <= Pr_rho'(E) <= (1+eta) Pr_rho(E)
```

for every event E, including arbitrary unions of histories and all failure
records. Thus exactly the same events have zero probability. A geometric
readout using only this history-support relation has zero response throughout
that open state neighbourhood. It cannot supply nonzero matter backreaction.
Responsive weighted arrival/readout data remain possible; their observer
rule, sampling errors, volume, clock and geometric interpretation must be
constructed. This is not a no-go for every emergent geometry.

## Quantitative completion criteria

A complete same-reference construction must meet the following checks together.

1. Produce the global rho_V or a joint local approximation with certified
   source/reference trace error; a preparation archive may be marginalized
   only by a stated observation map. Fixed-q boundary sectors and rejects stay.
2. Bound (3)-(4) uniformly on the chosen growing region and local density
   regime. Choosing region size requires L^2/M -> 0, not merely M -> infinity.
3. Supply an interior-collar locality estimate and make the overlap generator
   error tend to zero on the same observable support. Fixed-region (5) fails.
4. Attach the resulting controlled exchange clock to the actual complete
   unread repair process, or retain their different dissipative generators.
5. Derive responsive weighted geometry from specified complete observations;
   support alone fails. Supply count/volume calibration, dimension, Lorentz
   propagation, independent stress/Ward and geometry/Bianchi variations.
6. Establish geometric modular normalization, area/edge identification,
   fixed-proper-volume generalized entropy stationarity, a common small-ball
   limit and universal positive Newton coupling before using G1--G6.

The existing source-causal and celestial/null-ray finite theorems remain
conditional on their actual coverings, read laws, Lorentz coordinate data
and null-balance premises. They supply no omitted premise in the list above.
No G1--G6 physical premise is discharged by this regional clock result.

## A concrete preparation target and its exact repair-only obstruction

There is a positive decomposition of this same transfer reference. For a
configuration x let d_x be the number of source seams with unequal endpoint
occupations, and let x~y denote an unordered allowed occupation exchange.
Then

```
rho = [sum_x (M-d_x)|x><x|
       + sum_(unordered x~y) |+xy><+xy|] / (M Z),
|+xy> = (|x>+|y>)/sqrt(2).                          (7)
```

All coefficients are nonnegative because d_x<=M. Each off-diagonal entry is
1/(2M Z), while each diagonal is (M-d_x+d_x/2)/(M Z), proving (7).
This supplies a simple finite protocol conditional on two explicit services:
uniform fixed-q configuration preparation/readout, and an addressed coherent
mixer that prepares |+xy> from its known basis endpoint. Select one seam
uniformly. Equal endpoint occupations accept and retain |x>. Unequal ones
accept only the orientation with the smaller port occupied, then execute
the coherent mixer; the other orientation rejects. Each unordered exchange
pair is counted once. Acceptance probability is exactly

```
c = Z/binom(V,q) = 1-u,
```

and the accepted source marginal is rho. The rejection probability is u;
all records and rejected source states remain in the complete protocol.
For H66 this acceptance is 28/33. It is a different declared sampler from
the existing binary-slot selector, whose acceptance is 105/128. Raw records
are correlated with the accepted state; no factorization of their archive
or of an initially correlated external reference is asserted.

This decomposition is not itself an admitted growing-source implementation.
It identifies the exact services and their positive probability cost, instead
of postulating a global arbitrary-state preparation. The supplied coherent
service must be attached to the source, the complete flags and its resource
budget before this protocol discharges a preparation premise.

In particular the mixer cannot be obtained by concatenating only the fair
C/E repair instruments with classical feedback. Each native Kraus branch is
a scalar partial injection in the occupation basis. It sends a diagonal
input to a diagonal output. Products, classical mixtures, record-adaptive
choices and retained-flag postselection preserve this property, by induction
on the finite history tree. A nonzero native exchange gives a nonzero
off-diagonal entry of rho in (7). Therefore a repair-only protocol starting
from the canonical diagonal source cannot prepare rho, even probabilistically.
The existing H66 preparation uses its declared coherent source services;
its source-readability does not remove this distinction. This is an exact
obstruction for the stated operation class, not a no-go for admitted coherent
extensions or different initially coherent states.
