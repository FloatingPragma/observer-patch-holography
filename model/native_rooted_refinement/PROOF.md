# Rooted source refinement

The source is the fair binary `port_pair` federation. Its graph grammar,
global conserved occupation, scheduler and arithmetic are readable parts
of its specification. This construction evaluates a rooted neighborhood
from that specification without allocating its remote sphere. It proves
an all-radius limit for the shared-midpoint **cell topology**, certifies
finite port-routed windows, and attaches the canonical state restriction
and complete local-history bound to those windows. It does not assert a
unique infinite limit of the floating port routing.

## 1. An exact corner topology

This section concerns the exact oriented subdivision stencil. The literal
`_orient_faces_outward` operation uses numerical coordinates. Its orientation
and resulting addresses agree with this stencil only on the finite local
windows where the separate authentication below verifies that agreement.
The exact stencil theorem is not an all-level floating-orientation theorem.

The five original oriented faces incident to vertex0 are

```
(0,11,5), (0,5,1), (0,1,7), (0,7,10), (0,10,11).
```

They have vertex0 in their first position. The source replaces every
oriented triangle `(a,b,c)` by

```
(a,ab,ca), (b,bc,ab), (c,ca,bc), (ab,bc,ca),
```

where each unordered edge has one shared midpoint identifier. Inserting
child0 immediately after the original-face address therefore embeds each
coarse corner fan into the corresponding fine corner fan. A level-l cell
has literal address `f*4^l+d`, where `f` is its original face. The embedding
is `f*4^(l+1)+d`; the pair `(f,d)` names the identified cell.

Put `n=2^l`. Within each original face use the exact combinatorial grid
vertices `(i,j)`, with `i,j>=0` and `i+j<=n`. Its small triangles are
`[(i,j),(i+1,j),(i,j+1)]` and
`[(i+1,j),(i,j+1),(i+1,j+1)]` when their vertices fit. The five radial
edges are identified cyclically. A dual step increases the largest
radial coordinate `i+j` of a visited triangle by at most one. The root
triangle has largest coordinate one. Consequently a dual ball of radius
`r`, together with every incident dual edge, stays within this corner fan
when

\[
                         2^l>r+2.                         \tag{1}
\]

Insertion of child0 preserves all those incidences. Induction proves that
their exact addressed signatures agree at every pair of levels satisfying
(1). This is a full-sequence result for that particular corner root,
without an appeal to subsequence compactness. `infinite_ball(r)` evaluates
the signature at the smallest sufficient level. Restriction to a smaller
radius is consistent, so these signatures construct a countable rooted
dual graph. The graph contains the persistent five-face corner defect;
it is not asserted to be the root limit of an arbitrary address sequence.

The recursive neighbor evaluator follows the literal four-child stencil.
Child3 has the three corner children as neighbors. A corner child's middle
side meets child3; its other two sides inherit the adjacent parent and the
child at the shared endpoint. Exact midpoint names identify those endpoints.
This evaluates a finite rooted ball without enumerating the `20*4^l` cells.
No coordinate or physical dimension enters this proof.

## 2. Local evaluation of the numerical source

`LocalSource` additionally reproduces the source's finite coordinate
operations. It starts with the literal base vertices, keeps inherited
vertices, normalizes each fresh endpoint sum, applies `_clean_coordinates`
once in `_refine_level` and again in `_make_level`, and executes the literal
outward-face permutation. Re-cleaning inherited vertices at every level is
required; dropping it would change the source expression. Midpoints are
interned by their level and unordered parent identifiers. These internal
identifiers may differ from global array indices, but equality of vertices
and literal cell addresses are unchanged. The midpoint sum is independent
of its two endpoint orders.

Cell centers use the source's sum and row normalization. The local frame
uses the exact branch `abs(normal_z)>0.9`; its alternative reference axis
is retained. Neighbor directions use the source normalization floor.
The evaluator rejects a local direction that reaches that floor. Thus it
does not replace a small source vector by an ideal tangent direction.

The finite local-expression domain is `0<=l<=14`. Every source mesh edge
key `low*vertex_count+high` then fits a signed64-bit integer, since
`vertex_count=10*4^l+2` and `(10*4^14+2)^2<2^63`. Cell and port addresses
also fit. This conservative bound avoids using an unproved modular-overflow
interpretation of the source's global edge grouping. It is separate from
the native host's evaluated-level limit of two.

All scalar coordinate operations before the score contractions are the
same row operations on the same stored inputs, independent of remote rows.
The two three-component contractions may use different legal double
accumulation orders in a local array and a global array. Under the usual
IEEE nearest-rounding model, with unit roundoff `u=2^-53`, a length-three
dot product has absolute error bounded by `gamma_5` times the sum of
absolute products, where `gamma_5=5u/(1-5u)`. The evaluator explicitly checks
finite direction, frame and template components of magnitude at most two.
The first dot product has absolute-product sum at most12 and error at
most `12 gamma_5`. Its three component errors, multiplied by template
components at most two, contribute at most `72 gamma_5` to the second
contraction. That contraction has absolute-product sum at most
`72(1+gamma_5)` and rounding error at most `72 gamma_5(1+gamma_5)`.
Allowing both accumulation orders and an absolute `2^-900` allowance for
subnormal rounding therefore gives the exact bound
`288 gamma_5+144 gamma_5^2+2^-900 < 2^-40`, checked as a rational inequality.
The component bounds exclude overflow; the subnormal allowance exceeds
the sum of the finitely many minimum-subnormal rounding errors here.
The bound includes fused versus unfused
accumulation; it assumes finite standard floating operations and excludes
an arbitrary inaccurate numerical backend. The certificate rejects any
row whose decision is not fixed by this allowance.

Each score interval is rounded to binary32 by exact rational
nearest-even comparisons. A row's chosen port `p` is certified only if its
lower score is strictly larger than every earlier port's upper score and
at least as large as every later port's upper score. This proves the
literal first-argmax result, including stable binary32 ties. If the three
chosen ports collide, the evaluator enumerates all `12*11*10` injective
assignments and requires one assignment's lower total to exceed every
alternative's upper total. This proves the same unique Hungarian optimum
regardless of the global row ordering. An ambiguous collision is rejected,
rather than assigned a convenient port.

The original constructor is evaluated independently at levels0 through3.
Every coordinate and routed neighbor row in its radius-three cell ball
agrees with the local evaluation; the coordinates agree bit for bit.
These tests authenticate the finite attachment. The locality argument
explains why the evaluator computes the same neighborhood expressions
at other admitted finite levels if the full constructor returns with its
source postconditions. No remote graph allocation or high-level native
execution is inferred from a local query.

## 3. Finite routed windows and complete records

For a radius-r cell ball, keep all twelve ports of each included cell and
every source seam internal to their union. The source provides each cell's
thirty original seams and its three certified external attachment slots.
The evaluator compares the numerical cell incidences with the exact
corner topology before issuing a window. A port path of length at most r
crosses at most r cell boundaries, so this window contains the radius-r
port ball about all twelve root ports.

The canonical port name is `((f,d),p)`, with its explicit source address
`12*(f*4^l+d)+p`. The name order is the source address order. Hence the
identifications preserve seam endpoint order, native coin meaning, raw
slot meaning and signed current. Transport of the source tensor factors
therefore intertwines every fine native C/E Kraus map, including zero
slots and an arbitrary untouched reference. The record export convention
is unchanged. This is a finite instrument identification after address
transport; it is not an isomorphism of complete global archives.

The replay certifies matching routed windows for radii1,2,4,8 at their
sufficient corner levels and at levels6,10,14. These are finite tables
with explicit certificates. Routing equality between untested levels or
at every radius does not follow from those tables. A mutated address,
endpoint orientation, score decision, proof input or extra field causes
authentication to fail.

## 4. A consistent state and the finite record comparison

For any declared rational `p in [0,1]`, take

\[
 V_l=240\,4^l,\qquad q_l=\lfloor pV_l\rfloor,\qquad
 \tau_l=I_{q_l}/\binom{V_l}{q_l}.
\]

The total is read from the specification; repair does not select p.
On v named ports, the exact canonical restriction theorem gives the
all-number hypergeometric mixture, without conditioning the local total.
Its half-trace distance from the v-port Bernoulli product is at most

\[
 \epsilon_l\le \min\left(1,{v(v-1)\over2V_l}
                            +v|q_l/V_l-p|\right).             \tag{2}
\]

The product states restrict exactly under every consistent port inclusion.
Thus the exact corner topology has a consistent state family for any
declared p and any declared consistent routing extension. The finite
canonical restrictions converge whenever the specified windows exist
cofinally. A common independent external reference is allowed in this
specialization. A correlated exterior or reference requires its actual
joint restriction bound; marginal convergence does not supply it.

For each authenticated routed window, the
[complete local-history theorem](../native_local_history_limit/PROOF.md)
applies with degree6, root size12, `M_l=630*4^l`, v window ports and K
internal seams. Observe the records on the thirty seams with both endpoints
in the twelve-port root, together with the final state on those root ports.
All other window and remote record values are traced by this observation
map. Observing crossing seams or more final ports requires a larger root
and its corresponding cone bound. With independent uniform attempts, repair probability one,
and a common declared time partition aligned to `M_l`, it gives

\[
 D(\text{native local output},\text{common Poisson window})
 \le\min(1,\epsilon_l+a_l+tK^2/M_l),\qquad
 a_l\le12(6t)^{r+1}/(r+1)!.                                \tag{3}
\]

All visible raw records, their within-bin order and final quantum
disturbance are retained; only the stated remote records and exact tick
addresses are traced. This comparison uses the native per-seam attempt
clock. It makes no assertion about physical time. Two certified matching
windows have distance at most the sum of their two bounds in (3).

The displayed radius-four, level-fourteen example has348 ports and a
bound below `1/10^6` at `p=1/3`, `t=1/630`. Its description contains
`268435456` native attempts, within the single-program attempt cap.
Admission is verified. The expanded interpreter returns its unevaluated
carrier stop at this level; the comparison is a theorem about the declared
kernel and the certified finite geometry, not a reported execution of that
experiment. An infinite-process bound is recorded only conditionally on
an all-radius routed continuation.

## 5. Arithmetic continuation and its boundary

The unbounded nonnegative level grammar is distinct from geometry
evaluation. `execute` checks `level>2` and returns
`RESOURCE_STOP_UNEVALUATED_CARRIER` before density validation or graph
construction. The level200 control verifies zero attempts and this retained
stop. It is not a failed attempt to allocate an enormous graph.

The literal geometric implementation is not itself a total infinite
arithmetic specification. Its geometry hash writes a four-byte unsigned
level tag; a small-array call at `level=2^32` raises `OverflowError` at
that operation. Repeating the literal corner-child arithmetic64 times
also produces a zero spherical-area triangle in the bounded local control.
These tests execute small local operations only. They do not execute the
global source at those levels, prove a universal failure level for every
backend, or contradict the host's explicit resource stop.

A selected arithmetic continuation must specify vertex representation,
normalization, score comparison, tie handling and address encoding, and
prove matching routed balls at every radius. Exact real normalization,
arbitrary precision and a combinatorial routing rule are different choices;
this package selects none of them. The all-radius topological producer
and consistent state family reduce that remaining premise to the routed
arithmetic continuation. No unique full-sequence routed limit, constructive
routed subsequence, physical capacity or density, Lorentzian spacetime,
interacting chiral action or particle identification is asserted.
