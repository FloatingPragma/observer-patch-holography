# Continuous equivariant degree: proof proposed for acceptance

Acceptance status: **FROZEN_PAIR_AND_ALGEBRAIC_CROSSCHECK_ONLY**. The
argument below is a complete proof proposed for mathematical acceptance.
Its reconstruction and finite controls are not acceptance by a distinct
reviewer. The executable verifiers check finite data and representation
typing; they do not certify the differential-topological argument.
No `CONTINUOUS_TWIN_DEGREE_FLOOR_PROVED` gate is emitted.

This proof expands on the degree-floor statement and proof plan proposed by
B. Shatto (`dmobius3`) in the
[28 October follow-up to discussion #786](https://github.com/FloatingPragma/observer-patch-holography/discussions/786#discussioncomment-19246623).

## Oriented objects and theorem

Put `phi=(1+sqrt(5))/2`, with `sigma(phi)=1-phi`. The vertices `v_i` of the
convex realization, in port order 0 through 11, are

```text
(0,1,phi), (1,phi,0), (phi,0,1), (-1,phi,0),
(0,-1,phi), (phi,0,-1), (-phi,0,1), (0,1,-phi),
(1,-phi,0), (-phi,0,-1), (-1,-phi,0), (0,-1,-phi).
```

The twin vertices are `v'_i=sigma(v_i)`, with the same integer labels.
Let K be the abstract triangulated surface whose ordered faces are the
following list from `Lean/ObserverPatchHolography/CoreAxioms.lean`,
`CoreAxioms.orientedFaces`:

```text
(0,1,2), (0,3,1), (0,2,4), (0,6,3), (0,4,6),
(1,5,2), (1,3,7), (1,7,5), (2,8,4), (2,5,8),
(3,6,9), (3,9,7), (4,10,6), (4,8,10), (5,7,11),
(5,11,8), (6,10,9), (7,9,11), (8,11,10), (9,10,11).
```

The sum of these ordered simplices is the fundamental cycle `[K]`. Each
edge occurs twice with opposite boundary signs, every vertex link is a
five-cycle, and the connected complex has `(V,E,F)=(12,30,20)`. In the
plus coordinates these triangles are the boundary facets of a convex
icosahedron containing the origin. This also identifies K with S² without
using a classification of maps.

For barycentric coordinates on any listed face define

```text
r_+(sum_i t_i i) = (sum_i t_i v_i)/||sum_i t_i v_i||,
r_-(sum_i t_i i) = (sum_i t_i v'_i)/||sum_i t_i v'_i||.
```

Each face's three vectors are linearly independent, in both realizations,
so these formulas never normalize zero and agree on common edges. The map
`r_+` is a homeomorphism by convexity; `r_-` need not be injective.
Both target unit spheres have the boundary orientation induced by the
ambient orientation **minus** standard xyz. The determinants of the three
ordered vertex columns are `-2(phi+1)` on every plus face and
`-2(2-phi)` on every twin face. Thus `r_+` preserves the orientation `[K]`.

Let G be the sixty permutations preserving the oriented incidence, with
its fixed identification with A5. For each permutation g there is a unique
orthogonal matrix `rho(g)` taking every `v_i` to `v_(g i)`: the vectors span
R³ and their pairwise Gram values depend only on graph distance. Uniqueness
proves the multiplication law and faithfulness. The plus matrices are the
proper rotations of this oriented convex icosahedron. Their entries lie
in Q(phi). Applying sigma to their action equations gives
`rho'(g)=sigma(rho(g))`, acting on the twin vectors. Orthogonality and
determinant one survive these polynomial identities, so rho' is also a
faithful representation into SO(3). This argument transports algebraic
identities, not positivity through a field automorphism.

Consequently `r_±(g u)=rho_±(g) r_±(u)` for all `u in |K|`.
Write S for the plus unit sphere, S' for the twin unit sphere, and use
`r_+` to identify |K| with S, preserving orientation.

**Theorem.** Every continuous `f:S -> S'` satisfying

```text
f(rho(g)x) = rho'(g)f(x)    for every g in this same marked G and every x
```

has degree in `{7,13,17,23}` modulo 30. Hence `|deg(f)| >= 7`, with equality
attained by `r_- r_+^(-1)` of degree +7. In the untwisted problem, with rho
on both spheres, the necessary residues are `{1,11,19,29}` and the identity
attains +1. These are necessary conditions, with the stated witnesses;
no sufficient classification of all degrees is asserted. No channel,
linear-map or rational-map restriction is imposed on f.

## Lemma 1: the cyclic degree congruence

Let C_p, p prime, act faithfully by proper rotations on two oriented unit
spheres. Fix the **same generator** in both actions. Name the source poles
`x+,x-` and the target poles `y+,y-`. Suppose the continuous equivariant map
f sends `x+` to `y+` and `x-` to `y-`. In oriented tangent coordinates at
`x+` and `y+`, let the generator multiply by `exp(2πia/p)` and
`exp(2πib/p)`, where a and b are nonzero modulo p. Then

```text
deg(f) = b a^(-1) (mod p).
```

Only C_p-equivariance is required of the homotopies used to prove this
lemma. They need not preserve an ambient A5 action.

### Local replacement without an isolated initial fibre

Choose oriented orthonormal tangent bases and exponential charts at all
four poles, of radii less than π. A rotation fixing a pole satisfies
`R exp_x(v)=exp_x(Rv)`, so its action in each chart is exactly linear.
The chart disks are invariant, and target chart disks are convex in their
Euclidean coordinates. By continuity, disjoint invariant closed source
disks D+ and D- can be chosen with images in disjoint target chart disks
about the corresponding poles. Their images are compactly inside those
charts after shrinking the source disks.

Choose `1 <= k <= p-1` with `ak=b (mod p)`. On D+, let u(z) be the target
coordinate of f, and let chi(|z|) be a smooth radial cutoff, equal to one
on a smaller disk and zero near the boundary of D+. Choose epsilon>0 so
`m(z)=epsilon z^k` lies strictly inside the target coordinate disk for all
z in D+. The formula

```text
u_t(z) = (1-t chi(|z|))u(z) + t chi(|z|)m(z),   0 <= t <= 1,
```

stays in that convex disk. The relation `ak=b` makes both summands
equivariant. The formula agrees with f near the outer boundary and fixes
the value at x+, so it glues to an equivariant homotopy of the whole
sphere. Its endpoint is exactly `epsilon z^k` near x+ and has local degree
k there, including when k>1 and its derivative vanishes. The boundary loop
is `epsilon r^k exp(ik theta)`, of winding k. No property of the original
fibre, which could contain a neighborhood of x+, has been assumed.

On D-, perform the separate interpolation `u_t=(1-t chi)u` in the chart
centered at y-. It becomes constant at y- on a smaller disk. Throughout
this modification the image lies in the chart disjoint from y+, so it
introduces no preimage of y+ there. Call the resulting continuous map f0.
It is smooth on neighborhoods of two smaller closed protected disks P±.
Choose P+ strictly inside the polynomial region and P- strictly inside
the constant region. Choose boundary collars entirely within these
respective regions, so they avoid y+.

The cutoff interpolation can create additional preimages of y+ in the
transition annulus. This annulus lies outside P±, in the free source M
defined below, and is disjoint from their protected boundary collars.
No isolation or regularity of these additional preimages is assumed.
After relative smoothing, every preimage of y+ in this region belongs
to the compact set A used in the free-part transversality argument.

### Relative smoothing and the normalization bound

Here is a direct construction, without an equivariant approximation
theorem. For any delta>0 a smooth R³-valued map A uniformly delta-close to
f0 exists: take a finite cover with oscillation of f0 less than delta,
choose a smooth partition of unity subordinate to that cover, and at each
term use the constant vector f0 at one point of the corresponding set.
The triangle inequality gives the strict uniform bound. Smooth partitions
and cutoffs on S² are obtained from chart bump functions and division by
their positive sum.

Choose a smooth cutoff eta equal to one on neighborhoods of P± and their
boundary collars, supported inside the regions where f0 is smooth. Then
`F=eta f0+(1-eta)A` is smooth, equals f0 on those neighborhoods and satisfies
`||F-f0||_infinity < delta`. Average it in the ambient R³:

```text
Fbar(x) = (1/p) sum_g rho_target(g)^(-1) F(rho_domain(g)x).
```

The protected neighborhoods can be chosen invariant. On them every term
is f0(x); globally equivariance of f0 and orthogonality give
`||Fbar-f0||_infinity < delta`. Choose `0 < delta < 1/4`. Since `||f0||=1`,
each vector `(1-t)f0+t Fbar` has norm greater than `1-delta>0`. Normalizing
this straight segment gives a continuous equivariant homotopy, relative
to the protected neighborhoods. Its endpoint h is smooth and sphere-valued.
This proves the required relative smoothing with a strict nonzero bound.

### A transverse perturbation on the free source

Let `M=S² minus (int(P+) union int(P-))`. This is a compact smooth surface
with boundary. C_p acts freely on M, since every nonidentity element of
this prime rotation subgroup has exactly the two removed poles as its
fixed set. Choose a closed invariant collar C of the boundary on which
h avoids y+. It has a positive minimum distance from y+.

Take an invariant target exponential chart at y+, identified with a disk
`B_R(0) in R²`, and choose r>0 with `3r<R` so small that
`A=h^(-1)(exp_y+(closed B_r))` in M is disjoint from C. A is compact. If
A is empty no perturbation is required. Otherwise, every point of A has
a small open neighborhood U with pairwise disjoint C_p translates, with
closures away from C and with all images in `exp_y+(B_(2r))`. This follows
from freeness, finiteness of the group and continuity. Choose finitely many
such neighborhoods U_i and smooth nonnegative bumps beta_i with compact
support in U_i, whose positive sets and their translates cover A.

For independent parameters `t_i in R²`, define a coordinate displacement
on the translates by

```text
V_t(x) = sum_(i,g) beta_i(g^(-1)x) R_b(g) t_i,
```

where `R_b(g)` is the linear target action in the chart at y+. Each fixed
i contributes at most one g at any x, since its translated supports are
disjoint. The sum obeys `V_t(gx)=R_b(g)V_t(x)`. On the union of the supports
replace the coordinate w(x) of h by `w(x)+V_t(x)`; off that union leave h
unchanged. Compact supports inside the chart preimage make this a smooth
globally defined equivariant map h_t. Choose the parameter ball so small
that `||V_t||<r/2` everywhere. Coordinates on all supports then lie in
`B_(5r/2)`, strictly inside the target chart. Scaling t from zero to its
chosen value defines a homotopy fixed near the boundary of M.

Every zero of the perturbed coordinate has `||w||<r/2` before perturbation,
hence lies in A. At such a point some beta_i(g^(-1)x)>0. Differentiating
with respect to the corresponding pair t_i gives the invertible linear
map `beta_i(g^(-1)x) R_b(g)`. Thus the evaluation map
`E:int(M) x B -> S²`, `E(x,t)=h_t(x)`, is a submersion at every point of
`E^(-1)(y+)`. Its zero set Z is a smooth manifold of dimension `dim(B)`.

Apply Sard's theorem to the smooth projection `pi:Z -> B`. Its critical
values have measure zero, so every parameter ball contains a regular
value. At `(x,t) in Z`, the tangent equation is
`D_x E v + D_t E s=0`. Because `D_t E` is onto, surjectivity of `D pi` is
equivalent to surjectivity of `D_x E`. Therefore a regular parameter makes
h_t transverse to y+ on M. The inverse function and regular level-set
theorems used here apply to smooth maps between finite-dimensional,
second-countable manifolds without boundary; Z lies in `int(M) x B`.
Sard's Euclidean statement applies in a countable atlas. Boundary collars
are protected and have empty fibre.

The target action need not be free: at y+ it fixes the value and acts
invertibly on its tangent plane. Freedom of the **source** gives the
disjoint translates, so no isotropy constraint is imposed on the pair of
perturbation parameters. Overlaps between different U_i cause no problem,
since one nonzero pair already gives surjectivity.

The fibre on M is closed in compact M, avoids its boundary and is discrete
by transversality. It is therefore finite. There are no extra preimages
in the protected disks except x+. All modifications are homotopies.

### Adding the local degrees

The final fibre of y+ is x+ together with free orbits of p regular points.
Equivariance conjugates the local map at x by orientation-preserving
source and target rotations to the local map at gx. Their local degrees
are equal. Each free orbit therefore contributes a multiple of p.
The local degree at x+ is k. Additivity of degree over a finite fibre gives
`deg(f)=k+p m`, proving the congruence.

For specificity, the external topological results used are Sard's theorem
for smooth finite-dimensional maps, the inverse/regular level-set theorems,
homotopy invariance of degree, and the finite-fibre local-degree formula.
Sard and the level-set argument are given in Milnor,
[*Topology from the Differentiable Viewpoint*, §§2–3](https://people.dm.unipi.it/benedett/MILNOR-TDVPOINT.pdf).
The degree used here is defined by `f_*[S²]=deg(f)[S²]` in integral homology.
For a finite fibre, excision splits its local homology into one Z per
preimage; the fundamental class maps to the vector of local orientations,
giving the sum formula even at a critical isolated point. See Hatcher,
[*Algebraic Topology*, §2.2, Proposition 2.30](https://pi.math.cornell.edu/~hatcher/AT/ATch2.pdf).
The construction above verifies finiteness rather than assuming it for f.

### Independent local check and essential hypothesis

For any equivariant local map with an isolated zero, restrict to a small
invariant circle without zeros. Replace the generator by its a^(-1) power,
so its source angle is 2π/p and its target angle is `2πc/p`, where
`c=b a^(-1) mod p`. Lift the argument of the nonzero image loop to a real
function theta(s). Equivariance implies
`theta(s+2π/p)-theta(s)=2πc/p+2πn`, with n an integer constant by continuity.
Iterating p times gives winding `c+pn`. This checks the local congruence
without derivatives or transversality. It is independent of the smoothing
construction, but alone does not isolate the fibre of an arbitrary map
or establish the global congruence.

If x+ maps to the other target pole, name that pole y+; its oriented
tangent weight is the negative of the weight at the opposite pole. Pole
separation is essential: a constant map to a target pole is C_p-equivariant
of degree zero when both source poles have that same image. Also,
`z -> z^k` for k>1 has derivative zero at its isolated zero and local degree
k. An invertible derivative at a pole is not a permissible extra premise.

## Lemma 2: A5 forces the required pole separation

For any nonidentity g of order p, a point fixed by rho(g) maps under f to a
point fixed by rho'(g). Both fixed sets consist of exactly two poles.
It remains to prove that these two images are distinct, using the **same
abstract element** of A5 on the two sides.

For p=3, represent the subgroup in A5 by `<(123)>`; the even involution
`n=(23)(45)` conjugates (123) to its inverse. Its normalizer is S3 of order
6: a normalizing permutation preserves the three-letter support, and its
action on the complementary two letters is determined by parity.
For p=5 take `g=(12345)` and `n=(25)(34)`. Again `ngn^(-1)=g^(-1)`.
The normalizer in S5 consists of the twenty affine permutations of Z/5;
the even ones have multipliers ±1, giving the dihedral normalizer of order
10 in A5. Conjugating these descriptions covers every order-three or
order-five subgroup.

In either faithful SO(3) representation, n preserves the unordered axis
poles, by the conjugation relation. If it fixed a pole, it would be a
rotation about the same axis, hence commute with g. For odd p this
contradicts inversion of a nonidentity order-p rotation. It therefore
exchanges the poles in **both** representations.

For p=2 take `g=(12)(34)`. Its normalizer equals its centralizer and is
`{1,(12)(34),(13)(24),(14)(23)}`, of order 4. Choose either of the other
involutions as n. In a faithful SO(3) representation g and n are distinct
commuting half-turns. Write their unit axes as u and v. The half-turn
about u is `H_u=2u u^T-I`. Commutation implies that H_u preserves the
one-dimensional +1 eigenspace of H_v, so `H_u v=±v`. The plus case would
make their axes identical and their half-turn matrices equal, violating
faithfulness. Hence v is perpendicular to u, and `H_v u=-u`. Thus n swaps
the poles of g in both representations. The identity `g=g^(-1)` alone
supplies none of this geometry.

In all three cases, choose x-=rho(n)x+. If y+=f(x+), equivariance gives
`f(x-)=rho'(n)y+`, the other target pole. Lemma 1 applies.
The finite normalizer checks in `verify_port_support_source_audit.py`
reconstruct sizes 4,6,10 and pole swaps in the incidence action. These
finite calculations corroborate the group hypotheses; the SO(3) arguments
above establish the assertion simultaneously for the two marked actions.

## Weights, residues and sharpness

Use the order-five incidence permutation

```text
g=(0,2,4,1,6,8,3,5,10,7,9,11),
```

which fixes ports 0 and 11. Reconstruct its rotation by
`[v_(g0) v_(g1) v_(g2)] [v_0 v_1 v_2]^(-1)` and conjugate entries for the
twin. The traces are phi and `1-phi`. A proper order-five rotation has
trace `1+2cos(2πa/5)`, so the source weights at its two poles are ±1,
whereas the twin weights are ±2. Indeed
`2cos(2π/5)=phi-1` and `2cos(4π/5)=-phi`. These identities also follow
from the fifth-root relation `1+z+z²+z³+z⁴=0`, choosing the positive
root for `2cos(2π/5)`. The ratio b/a is therefore ±2 modulo 5,
regardless of which source pole is chosen. A trace fixes the unordered
pair of weights; it does not select the sign at f(x+).

At order three, all weights are units modulo 3, so b/a is ±1. At order
two the unique nonzero weight is 1 in both representations, so degree is
odd. In the untwisted order-five comparison the ratio is ±1. The exact
Chinese remainder calculation is

| Marked target | d mod 2 | d mod 3 | d mod 5 | d mod 30 | Least absolute degree allowed |
| --- | --- | --- | --- | --- | --- |
| twin | 1 | 1 or 2 | 2 or 3 | 7,13,17,23 | 7 |
| untwisted | 1 | 1 or 2 | 1 or 4 | 1,11,19,29 | 1 |

For example the twin residue representatives nearest zero are
`7,13,-13,-7`. All larger representatives differ by 30. These restrictions
concern the single global degree, with no sign chosen separately per face
or preimage.

The witnesses are continuous maps in this exact category. For each ordered
face and a direction y, solve `[v_a v_b v_c] c=y`, separately in the two
frames. A positive triple c gives exactly one interior preimage under the
radial map, with local sign `-sign(det[v_a v_b v_c])` for the prescribed
target orientation. All nonzero coefficients exclude edge preimages.
For y=(1,2,4) the plus face indices, counted from zero in the displayed
face list, are `{2}`; the twin indices are `{0,3,6,8,16,17,18}`.
All their local signs are +1. The directions (-3,5,2) and (7,-2,3)
give the same signed totals +1 and +7. The Cramer-rule reconstruction in
`verify_galois_port_frame_independent.py` verifies these data over Q(phi)
without importing the producer or using floating tolerances. Since
`deg(r_+)=1`, `deg(r_- r_+^(-1))=7`. The identity has untwisted degree 1.
No source-selected relative support map is used in either construction.

A continuous injection S² -> S² is open by invariance of domain in
two-dimensional charts and has closed image by compactness. Connectedness
of the target makes its image the whole sphere; compact-to-Hausdorff
bijectivity makes it a homeomorphism. Its induced isomorphism on H_2(S²; Z)
has degree ±1, excluded by the twin residues. The version of invariance
of domain used is the openness of a continuous injection between open
subsets of R², Hatcher, Theorem 2B.3 in the linked Chapter 2. Thus no
twin-equivariant continuous injection exists.

## Marking and orientation controls

An orientation-preserving recharting of either sphere conjugates its
representation and carries the map with it, preserving degree. Relabeling
the abstract group by the same automorphism in both actions merely
reindexes the equivariance equation. For any generator change `g -> g^u`
with u a unit modulo p, both weights multiply by u, preserving b/a.
Changing the target generator alone changes that ratio and the marked
problem. An outer automorphism inserted on one side exchanges the two
order-five classes; it is not a recharting of the fixed equation.

Reversing just the domain orientation or just the codomain orientation
negates degree and the corresponding tangent weight, consistently with
Lemma 1. Reversing both preserves degree. Reversing the common codomain
orientation for the two radial witnesses changes their degrees to -1 and
-7. Both residue sets are closed under negation, so their absolute floors
are unchanged. A convention varying by branch cannot be used to compare
the two signed witnesses.


## Independent algebraic cross-control and lift boundary

`equivariant_degree_crosscheck.py` constructs all nine irreducible
characters of 2I over Q(phi). Opposite lifts, including the central -1,
remain separate classes. Q is the defining SU(2) representation; Q' has
character sigma(chi_Q), with sigma(phi)=1-phi. Its symmetric square has the
twin SO(3) character. Character inner products give the McKay graph and the
two multiplicity sequences through degree 94. The independent verifier
starts with the affine E8 graph and uses `v_n=A v_(n-1)-v_(n-2)` instead.

The unique Q' in Sym^7(Q*) gives a nonzero covariant homogeneous polynomial
pair P of degree seven. Its common projective zero set, if nonempty, is
A5-invariant. A point stabilizer in the rotation action is cyclic of order
at most five, so an orbit has at least twelve points. Two degree-seven
polynomials cannot have twelve distinct common projective zeros. Hence
P has no base point and defines a degree-seven rational map. This existence
argument does not require matching HMS coordinates to the OPH labels.

Conversely a rational map with projective equivariance and a primitive
homogeneous pair of degree d satisfies `P(gv)=c(g) Q'(g)P(v)`. The multiplier
is constant: two primitive homogeneous pairs of the same degree defining
the same projective map differ by a nonzero constant, by coprimality in
C[x,y]. Composition makes c a character of 2I. For completeness its
abelianization vanishes: in the binary (2,3,5) presentation
`x²=y³=z⁵=xyz`, the abelian relation matrix

    [ 1 -1 -1 ]
    [-1  2 -1 ]
    [-1 -1  4 ]

has determinant -1. Thus c=1. The central -1 then forces d odd. This
justifies the homogeneous lift for this **rational** category; it says
nothing by itself about an arbitrary continuous map.

The HMS (3.18) pair has gcd one over Q and no common zero at infinity.
The producer uses rational Euclidean division; the independent verifier
uses coprimality modulo 101 with nonvanishing leading coefficients. This
checks base-point freeness of that explicit formula. Source-typed full
equivariance of that coordinate formula is not asserted by the gcd test.

## Literature and application boundary

[Houghton–Manton–Sutcliffe, hep-th/9705151v2](https://arxiv.org/html/hep-th/9705151)
gives degree-seven formulas in (3.18), the family in (3.19), and the
two-dimensional target representation in (4.6). These are rational-map
statements, not a proof of Lemma 1 for continuous maps.

[Kostant, math/0411142v1, Example 1.7](https://arxiv.org/pdf/math/0411142)
lists the near/far two-dimensional E8 numerators with denominator degrees
12 and 20. The checked numerators are {1,11,19,29} and {7,13,17,23}.

[Doyle–McMullen, §5, Proposition 5.3](https://math.dartmouth.edu/~doyle/docs/icos/icos/node5.html)
classifies the four commuting rational maps of degree below 31, with degrees
1,11,19,29. It does not classify twin-equivariant continuous maps. The
category difference prevents importing it as the universal degree theorem.

Even if the continuous theorem is accepted, applying it to the A1 support
requires a source-selected relative map, the same oriented realized target
sphere, compatible refinement maps, and a proved equivariance equation.
The seed/degree receipts do not supply these premises. Any degree-one twin
exclusion is conditional on them; physical identification is unfulfilled.
