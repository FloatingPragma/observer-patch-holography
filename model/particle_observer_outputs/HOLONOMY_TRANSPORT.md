# Creating and transporting a finite loop-holonomy defect

A declared link register on the actual icosahedral carrier gives an explicit
creation, transport and annihilation experiment. It retains a work register,
every detector label and the no-creation history. The protecting condition
is exact: the allowed operations preserve a specified region's boundary
holonomy. Arbitrary local link operations do not satisfy that condition,
and the closed sphere has no nontrivial flat global loop sector.

The construction is a finite extension of the native carrier. Its link
bits, work register and controls are specified program data. The original
C/E load-and-memory program is not asserted to create these spatial
curvature defects. No particle name, physical spin, mass or coupling is
assigned to them.

## Reused source and theorem domains

[`holonomy_transport.py`](holonomy_transport.py) reads the authenticated
source incidence and exact coordinates through
[`native_loop_readout.carrier`](../observer/native_loop_readout.py) and
[`slow_geometry`](../common_source/slow_geometry.py). It reuses
[`native_defect_stability_gate`](../../code/a5_closure/native_defect_stability_gate.py)
for the full sphere homology and the declared Z6 pair-transport paths.
The [finite flux certificate](../../code/a5_closure/flux_defect_criterion_certificate.py)
and [routed seam grammar](../../code/a5_closure/routed_seam_grammar_certificate.py)
establish the link-to-face map and its sum-zero image. Their existing local
path construction creates, moves and erases a Z6 flux pair. This calculation
uses its order-two subgroup `{0,3}`, represented by link signs.

The additions are a gauge-invariant quantum preparation, an explicit local
work exchange, a regional protection test and complete observer histories.
The sixty-state matter transport in
[`gauge_defect_transport`](../gauge_defect_transport/PROOF.md) is a different
carrier and is not identified with these curvature bits.

The native current record has a separate exact result: on all 78 admitted
load configurations, the sum of its two coin-conditioned currents on edge
`(u,v)` is `n_u-n_v`. This is a vertex gradient, so every face-loop sum is
zero. The executable reruns all 4,680 current branches to check it. The
[corrected native memory object](../common_source/native_gauge_object.py)
has its own protected charge/projector and complete C/E recovery. Protection
of that memory object is not a proof that a face holonomy is born or moved.

## The complete gauge-invariant space

There are 12 vertices, 30 edges and 20 triangular faces. Put a bit `a_e`
on each edge, with sign `z_e=(-1)^a_e`. For a face f define

\[
 B_f=\prod_{e\in\partial f}Z_e,\qquad
 b_f=\sum_{e\in\partial f}a_e\pmod2,\qquad
 B_f|a\rangle=(-1)^{b_f}|a\rangle.                    \tag{1}
\]

A defect is a negative face-loop product. It is an observable mismatch of
the link transport around that loop. It is not a marked vertex or an
arbitrary label attached to a flat configuration.

The vertex gauge operation is `A_v=product_(e incident to v) X_e`.
Each triangle meets a vertex star in zero or two edges, so `[A_v,B_f]=0`.
Every edge borders two faces, hence `product_f B_f=I` and the total number
of negative faces is even. The exact boundary ranks over F2 are 11 and 19,
giving `H1(S2;F2)=0`. Thus every flat link configuration is a vertex gauge,
and two link configurations with the same face curvature differ by a
vertex gauge. The code obtains that gauge through the existing Z6 solver.

The full link space has dimension `2^30`. Its simultaneous `A_v=+1`
sector has dimension `2^(30-11)=2^19`. A normalized basis is

\[
 |[a]\rangle=2^{-11/2}\sum_{g\in\operatorname{im}\delta_0}|a+g\rangle.
                                                               \tag{2}
\]

Each even 20-face curvature pattern labels exactly one basis vector. No
curvature sector is omitted. The flat preparation is (2) at a=0. It is
gauge invariant; a single computational link configuration would not be
the same preparation. A reversible preparation can generate the eleven
independent vertex stars coherently and uncompute their gauge coordinates
using the fixed-root spanning tree. This uses declared coherent controls,
rather than an assertion that bare load repair prepares (2).

## Creation, motion and the work balance

Flipping one edge with `X_e` toggles exactly its two incident face signs.
From two flat faces it creates a pair. If one incident face is negative,
the same flip moves that endpoint to the other face. From two negative
faces it annihilates the adjacent pair. The diagram certificate follows
an existing dual path from face 0 to its antipodal face 13 and then returns.

Choose a positive cost J and a two-level work register n with energy 2J:

\[
 H_{\rm cost}=J\sum_f(1-B_f)/2+2J n.                  \tag{3}
\]

J is a declared program coefficient, with no GeV conversion. Equation (3)
assigns energy to defects; it does not derive inertial mass, a dispersion
law or the cost of implementing a controller.

For the incident faces f,g let `P00` project onto both flat and `P11`
onto both negative. Define

\[
 T_e=X_eP_{00}\otimes|0\rangle\langle1|,\quad
 W_e=T_e+T_e^\dagger+I-T_e^\dagger T_e-T_eT_e^\dagger. \tag{4}
\]

The initial and final supports of `T_e` are orthogonal. Therefore
`W_e^2=I`, and it holds every inactive input unchanged. It trades a work
quantum for the two new defects, or returns it upon annihilation.
Consequently `[W_e,H_cost]=0` on the complete link/work space. The move
gate applies `X_e` only when exactly one of f,g is negative and holds the
zero/two-defect cases. It is another complete involution commuting with
(3). Both gates commute with every vertex gauge operation. Their truth
tables include all local curvature and work inputs, not only a prepared
creation branch.

Let `P=T_e†T_e+T_eT_e†`. A coherent creation pulse is

\[
 U(c,s)=I+(c-1)P-is(T_e+T_e^\dagger),\qquad c^2+s^2=1. \tag{5}
\]

It is unitary and preserves (3). The fixed rational choice `c=3/5,s=4/5`
sends flat/work1 to `3/5` times flat/work1 minus `4i/5` times
pair/work0. These are configured pulse coordinates, not experimental
particle inputs. Complete readout gives original probabilities `9/25`
and `16/25`. Record-controlled moves transport the created endpoint along
the full path; the no-creation branch holds. Return and annihilation
restore flat/work1 in both histories, while the distinct records remain.
Total field-plus-work energy is 2J at every stage. No history is postselected.
The apparatus selects transport only for the exact prepared-pair/work0
record; every other first record, including outside-Gauss input, selects
hold. This fixes the control policy on the complete outcome alphabet.

## What is protected, and what can remove it

For any set R of faces, cancellation of internal edges proves the exact
finite Stokes identity

\[
 Q_R:=\prod_{f\in R}B_f=\prod_{e\in\partial R}Z_e.     \tag{6}
\]

A link update `X(u)=product_(e:u_e=1)X_e` preserves `Q_R` precisely when
`u` intersects `boundary(R)` an even number of times. Vertex gauge updates
always satisfy this condition. Updates supported away from the boundary
also satisfy it. Arbitrary sequences and coherent combinations of such
operations preserve the sector, as do their complete commuting instrument
maps with arbitrary untouched references.

If `Q_R=-1`, no permitted operation satisfying that condition can make
every face in R flat. Defects may move within R or be created/annihilated
in pairs; the statement protects regional parity, not each individual
defect's position or lifetime. The example takes R to be its starting
triangle. Separation and return preserve its negative loop product.
Creation and the final annihilation cross its boundary and change it.
Those operations are visibly outside the stated protecting operation set.
Gauge invariance alone therefore does not give permanent protection.

The order-two charge is self-inverse: the two endpoints have equal
nontrivial Z2 labels. The larger reused Z6 construction has inverse labels
`k` and `6-k`. Neither fact determines physical particle/antiparticle
assignments, spatial spin or statistics.

## Complete readout and graphic data

The apparatus reads all twenty face bits and the work bit. Its pointer
alphabet has `2^21` labels; the `2^20` odd-total-curvature labels have zero
projectors. On the gauge-invariant link/work space the even labels give
`2^20` rank-one projectors whose sum is the identity. For the full ambient
space retain an additional outside-Gauss outcome, without assigning it a
physical curvature-basis state. These complete projectors define a unitary
pointer coupling by the [readout contract](READOUT_CONTRACT.md). Tensoring
by any untouched reference preserves all identities.

The JSON certificate exports the actual exact Q(sqrt5) vertex coordinates,
display coordinates, every edge and oriented face, both incident faces for
each edge, and every step of the creation/motion/annihilation episode. Each
step contains all twenty holonomy signs, the negative link edges, work and
defect energies, boundary sign and whether that step preserves the selected
region. The drawn negative-edge path is a gauge representative of (2),
not an additional measured physical string. Vertex regauging can change
that path while leaving every face record fixed. The clock is the declared
sequence of controls and readouts; no relativistic propagation time follows.
The separate complete experiment lists its two nonzero histories
and the zero-complement convention. It contains no claimed particle masses.

```sh
python -m model.particle_observer_outputs.holonomy_transport
python -m pytest -q model/particle_observer_outputs/test_holonomy_transport.py
```
