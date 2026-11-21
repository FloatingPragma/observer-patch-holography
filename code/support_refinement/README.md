# Equivariant rational support refinement

The declared spherical support has 12, 42, 162 and 642 vertices. Its
midpoint subdivisions admit an exact **rational chain deformation retract**
at each transition. The coarsening, subdivision and homotopy commute with
all sixty proper seed actions. Composing transitions gives the same
conclusion on every interval of the retained tower.

This makes the standard subdivision invariance of simplicial homology
explicit for the actual source. The added result is the specified local
equivariant contraction and its complete finite source checks.

This realizes a compatibility interface for the supplied support in CP-08.
It does not select a relative map from an observer's local ports to support
vertices, identify spatial and internal representations, or construct
responsive geometry, a physical clock or a continuum field theory.
The source's endpoint-valued process coarsening is unchanged. Rational
chains and a deterministic map of vertices are different objects.

## Actual source and exact domain

`construction.py` checks the bytes of the existing support producer
`code/source_selection_model/geometry.py`, the ordered Lean boundary,
the existing independent triangle-incidence search, and the source
realization receipt. It reconstructs the live support, compares all four
rows with that receipt, and infers midpoint identities from fine incidence.
Source modules are executed from freshly checked bytes on every load.
The finite executor accepts exact rational chains of dimensions zero,
one and two; floating coefficients, unknown simplices, invalid ancestry,
orientation damage and omitted levels are rejected.

The result concerns the **support** complexes. It supplies no map on the
extra fibre/slot carriers of the source nerve, no quantum channel, and no
new source preparation or operational controller.

## Maps and local proof

Let `K'` be the midpoint subdivision of an oriented triangular complex `K`.
Write `m_ab=m_ba` for the midpoint of the old edge `{a,b}`. Work over `Q`.
The following proof applies to any such finite complex with distinct edge
midpoints and the stated oriented subdivision; the executor checks the
four actual source levels.

The vertex map is a linear averaging map:

```
C_0(a) = a,                 C_0(m_ab) = (a+b)/2.
C_k([v_0,...,v_k]) = C_0(v_0) wedge ... wedge C_0(v_k).
```

Each fine simplex is supported in a single old simplex, so the exterior
formula lands in `C_k(K;Q)`. The augmented boundary is contraction with
the functional that sends every vertex to one. Since `C_0` preserves that
functional, `d C = C d`; degree zero has the usual zero boundary.

The subdivision map fixes old vertices and sends

```
S[a,b] = [a,m_ab] + [m_ab,b],
S[a,b,c] = T_a + T_b + T_c + T_m,
T_a = [a,m_ab,m_ac],       T_b = [m_ab,b,m_bc],
T_c = [m_ac,m_bc,c],       T_m = [m_ab,m_bc,m_ac].
```

Interior boundaries cancel, hence `d S=S d`. The two small edges each
coarsen to half the old edge, and each of the four triangles coarsens to
one quarter of the old triangle. Thus `C S=id` on every basis simplex.

There is also an explicit local homotopy `H:C_k(K';Q)->C_(k+1)(K';Q)`:

```
H_0(a) = 0,
H_0(m_ab) = ([a,m_ab]+[b,m_ab])/2,
H_1([a,m_ab]) = H_1([m_ab,b]) = 0,
H_1([m_ab,m_ac]) = T_a - S[a,b,c]/4,
H_2 = 0.
```

Extend by orientation and linearity. Every interior fine edge has the
displayed form in its unique parent triangle. Reversing that edge reverses
both terms, so the formula is well-defined and agrees across old faces.
For a midpoint, `dH_0=m_ab-(a+b)/2`. For a boundary fine edge the residual
`id-SC-H_0d` is zero; for an interior fine edge it is the boundary of the
displayed `H_1`. On a corner triangle, `H_1d=T_a-S[a,b,c]/4`.
On the central triangle, the three interior edges instead give
`-T_a-T_b-T_c+3S[a,b,c]/4=T_m-S[a,b,c]/4`. Therefore

```
dH + Hd = id - SC,       CS = id,
CH = 0,                 HS = 0,                 H^2 = 0.
```

The last identities follow locally: the two terms of `CH_0` cancel;
`C(T_a)=C(S[a,b,c])/4`; subdivided old edges use only boundary fine edges;
and `H_0` itself uses only those edges. This is a strong chain deformation
retract. In particular, rational homology is preserved, with explicit
inverse maps rather than only matching Euler characteristics.

For any simplicial symmetry `g`, its lift sends `m_ab` to `m_(ga,gb)`.
Every displayed formula commutes with that substitution, including its
orientation signs. Thus `C`, `S` and `H` are equivariant. The certificate
checks all basis simplices under all sixty actual proper actions on each
of the three transitions; this includes the homotopy, not only vertex
averaging or fundamental cycles.

## Tower composition

For consecutive retracts `(C_1,S_1,H_1)` and `(C_2,S_2,H_2)`, set

```
C = C_1 C_2,       S = S_2 S_1,       H = H_2 + S_2 H_1 C_2.
```

Substitution gives `dH+Hd=id-S_2 S_1 C_1 C_2` and `CS=id`.
The three strong side conditions also follow: cross terms contain
`H_2 S_2=0` or `C_2 H_2=0`. Iteration gives associative coarsening and
subdivision, compatible equivariant homotopies, and all interval retracts.
No infinite tower, geometric mesh limit or continuum response follows just
from these finite certificates.

## Endpoint obstruction and interpretation

An order-two proper seed action fixes an edge midpoint while interchanging
its two endpoints. Consequently an equivariant midpoint map whose value
must be one of those endpoints is impossible. The certificate checks such
an involution for every one of the thirty seed edges. Under the lifted
actions, only the identity commutes with the source's declared `max(edge)`
coarsening at each of the three transitions.

The average `(a+b)/2` is fixed by the interchange, which is why this chain
construction avoids that obstruction. Its two coefficients cannot be
silently replaced by a chosen endpoint or interpreted as a simplicial
vertex image. The algebra establishes a natural rational chain interface;
a pointwise geometric attachment and its physical realization are separate requirements.

## Replay

```
python -m pytest -q -W error code/support_refinement
PYTHONPATH=code python -m support_refinement --check code/support_refinement/results.json
```

All arithmetic in the certificate is exact. The retained result records
the source and implementation hashes, complete basis counts, the three
transition checks, all six interval retracts and an endpoint-obstruction
witness. It accepts no physical completion gate.
