# Galois port frames and the support attachment boundary

Baseline: `FloatingPragma/observer-patch-holography` at
`4511928b1d4bee6bca5601e309b44595a76f797f`.

Verdict: `GALOIS_FRAME_AMBIGUITY_PERSISTS_WITHOUT_SUPPORT_ATTACHMENT`.
The finite geometric equivalence is proved for the specified oriented pair.
The source-bound attachment to the A1 support bridge is not present in the
audited corpus. This is outcome B with missing interface C, not outcome A.

## Selection boundary audit

Roles classify mathematical use, not Python or Lean import dependencies.
The exact baseline above pins every surface in this inventory.

| Surface | Role | Evidence and boundary |
|---|---|---|
| `Lean/Screen/PortFrameGram.lean` | `DECLARES_POSITIVE_EMBEDDING` | `g5` assigns `(0,1)` to adjacency in the scaled `5G` table. `gram_sq` and `g5_trace` verify the declared table, not its selection. |
| `Lean/Geometry/ScreenCarrierMapCandidate.lean` | `CONSUMES_PRESELECTED_FRAME` | `PortCarrierCandidate.dot_table` requires `gramTargetZ`; `candidate_gram_bridge` ties it to `g5`. Its negated candidate has the same Gram, not the Galois-conjugate Gram. |
| `Lean/ObserverPatchHolography/CoreAxioms.lean` | `BRANCH_NEUTRAL` | `orientedFaces` and `BoundaryPacket` give incidence and orientation; `FederationSupportShadow` carries a separate nerve/support bridge. Its degree fields are propositions with supplied proofs, not a computed degree of port rays. |
| `docs/AXIOM_REFERENCE.md` | `BRANCH_NEUTRAL` | A1.5 declares no Euclidean placement; A1.7 supplies an oriented support homeomorphism and degree-one nerve map; the A1 non-implications explicitly include selection between the Gram pair. |
| `code/a5_closure/sl2f5_mckay_e8_certificate.py` | `USES_BOTH_GALOIS_EMBEDDINGS` | `build_certificate` checks faithfulness of the conjugate doublet, its distinct irreducible character, and independently recomputed fusion graph isomorphism to affine E8. |
| `code/a5_closure/sl2f5_mckay_e8_certificate.md` | `USES_BOTH_GALOIS_EMBEDDINGS` | The Galois-control section explicitly denies selection from graph type. |
| `code/source_selection_model/DERIVATION.md` | `USES_BOTH_GALOIS_EMBEDDINGS` | Section 2 places `v_p` and `sigma(v_p)` in the same response; section 1 distinguishes each local boundary from the global triangulation. |
| `Lean/Screen/A5FamilyBand.lean` | `USES_BOTH_GALOIS_EMBEDDINGS` | Exact `3`/`3'` projectors and costs; `family_band_selected` derives the positive-cost minimizer conditional on its admitted cost comparison. It does not derive an A1 geometric attachment. |
| `claims/selection_ledger.json` (baseline row 3) | `CONSUMES_PRESELECTED_FRAME` | The row cites the positive Gram identity while its two-element menu concerns orientation-preserving versus full automorphisms, not the two Gram branches. |
| `docs/SELECTION_LEDGER.md` (baseline row 3) | `CONSUMES_PRESELECTED_FRAME` | Generated projection of that row; the menus must not be conflated. |

The repository search included `Lean`, `code`, `docs` and `claims`, querying
Galois/conjugation, Gram/support links, degree and support-map declarations.
The strongest adjacent results are
`PortGramRepairBand.portGram_repair_band_packet`,
`PortGramRepairBand.selected_family_band_is_port_gram`, and
`PortGramRepairCovariance.normalizedKernel_tendsto_portGram`. They identify
the positive band under a declared Laplacian comparison or repair mean.
The covariance module explicitly says that its repair mean is not selected
by the axioms. These are conditional algebraic selectors, not the missing
source-bound identification with `b_r`. No competing repair-band construction
is introduced here. The independent executable reconstructs the two frames
from oriented incidence to expose precisely that distinction.

The direct Lean search for `FederationSupportShadow`, `supportMap` and
`orientedDegreeOne` resolves only to their definitions in `CoreAxioms`.
`CarrierDynamicsCompatibility` transports a supplied positive carrier map
and its barycentric refinements without a federation-support attachment.
The broader source receipts `code/rg_principle/DERIVATION.md` and
`code/source_operation_reads/DERIVATION.md` explicitly list a source-bound
degree-one support completion among their unproved obligations. The
degree-one phrase in `cone_whitney_bridge.py` concerns polynomial degree,
not topological degree of a local Gram-to-A1 bridge.

## Exact input and arithmetic

The primary input is the twenty ordered triples in
`CoreAxioms.orientedFaces`. The producer reads this declaration, constructs
adjacency from face edges and derives all graph distances. The reference
manifest pins the source file's SHA-256. The independent verifier rereads
the declaration, checks the oriented fundamental cycle, and reconstructs
distances by a separate breadth-first traversal. Neither executable imports
`PortFrameGram.g5` or a stored positive Gram target.

Write `s=sqrt(5)>0`, `phi=(1+s)/2`, and `sigma(s)=-s`. For distances
`0,1,2,3`, respectively,

```
G_plus  = 1,  s/5, -s/5, -1
G_minus = 1, -s/5,  s/5, -1.
```

The producer uses rational pairs in the basis `1,s`; the independent
verifier uses the basis `1,phi` with `phi^2=phi+1`. Both compute
`G^2=4G`, symmetry and `tr(G)=12`. Since `G/4` is a real symmetric
idempotent, its eigenvalues are zero or one; its trace and rank are three.
There is also a direct rank and positivity witness: `G=VV^T/r^2`,
`r^2>0`, and the first three row vectors have nonzero determinant.
Consequently `x^T(G/4)x=||V^T x||^2/(4r^2)>=0`. Positivity is verified
separately at both real embeddings; Galois conjugation is not assumed to
preserve arbitrary inequalities.

## Coordinates and point-set identification

The plus coordinates, in committed port order, are

```
0 (0,1,phi)     1 (1,phi,0)      2 (phi,0,1)
3 (-1,phi,0)    4 (0,-1,phi)     5 (phi,0,-1)
6 (-phi,0,1)    7 (0,1,-phi)     8 (1,-phi,0)
9 (-phi,0,-1)  10 (-1,-phi,0)   11 (0,-1,-phi).
```

The minus vectors apply `sigma(phi)=1-phi=-1/phi` to every coefficient.
The squared norms are `(5+s)/2` and `(5-s)/2`. Division of each inner
product by its branch's squared norm yields the incidence-derived Gram.

Let `S(x,y,z)=(y,x,z)`. The exact positive similarity
`phi S V_minus[i] = V_plus[p[i]]` has permutation

```
p = [5,3,4,10,9,11,0,2,1,7,8,6].
```

Since `r_plus^2=phi^2 r_minus^2` and `phi>0`, normalization gives
`S(V_minus[i]/r_minus)=V_plus[p[i]]/r_plus`. Thus the ordinary unlabeled
unit point sets are orthogonally congruent. `det(S)=-1`; this is not an
orientation-preserving identification of the committed face maps. The
permutation exchanges graph distances one and two and preserves antipodes.
Abstract A1 adjacency, Euclidean nearest neighbors, Euclidean second
neighbors and oriented face incidence remain separate data.

The Galois nearest/second-neighbor relation is established geometry, not a
novelty claim. See John Baez, [The Great Icosahedron](https://johncarlosbaez.wordpress.com/2026/05/27/the-great-icosahedron/)
and [From Pentagons to Pentagrams](https://johncarlosbaez.wordpress.com/2026/05/29/from-pentagons-to-pentagrams/),
which points to John H. Conway, Heidi Burgiel and Chaim Goodman-Strauss,
*The Symmetries of Things*, A K Peters, 2008. The finite proofs here do not
depend on those references.

## Four diagnostics and exact degree

The marked proper incidence permutation is
`g=[0,2,4,1,6,8,3,5,10,7,9,11]`. It fixes port zero and advances its
link by one step, sending port one to port two. Its order is exactly five.
The rotation is reconstructed from three independent vectors; orthogonality,
determinant one and the action on all twelve vectors are checked. The
minus rotation is both reconstructed directly and checked to equal the
coefficientwise conjugate of the plus rotation. Merely specifying order
five would not select a trace: `g^2` exchanges the two order-five classes.

| Diagnostic | Plus | Minus |
|---|---|---|
| Trace of this marked `g` | `phi` | `1-phi` |
| Adjacency coordinate eigenvalue | `s` | `-s` |
| Laplacian coordinate eigenvalue | `5-s` | `5+s` |
| Normalized squared chords at distances 1,2,3 | `2-2s/5, 2+2s/5, 4` | `2+2s/5, 2-2s/5, 4` |
| Rank of abstract-edge chord length | first | second |
| Degree in Cartesian ambient orientation | `-1` | `-7` |
| Degree in the fixed convex-positive ambient orientation | `1` | `7` |

The map on each ordered face is
`(a,b,c) -> (a v_i+b v_j+c v_k)/||a v_i+b v_j+c v_k||`,
for nonnegative barycentric coordinates with sum one. Each face's three
vectors are linearly independent. Thus its affine triangle avoids zero,
and radial maps on shared edges agree. The committed closed oriented
icosahedral boundary is the domain sphere, even when its image intersects
itself. Faces of the minus branch need not be convex-hull faces of its
point set.

For the rational regular value direction `y=(1,2,4)`, the producer solves
the twenty linear systems by Gaussian elimination. The independent
verifier uses Cramer's rule. Every coefficient is nonzero, so the direction
lies in none of the face-edge planes. A face has one preimage precisely
when its three coefficients are strictly positive; dividing them by their
sum produces its barycentric preimage. Its local orientation sign is
`sign det(v_i,v_j,v_k)`. This sign follows by differentiating radial
projection in the ordered simplex coordinates: positive normalization
factors leave the determinant sign unchanged. Hence the signed count is
the topological degree of the continuous map between the oriented spheres.
The degree theorem makes one regular value sufficient; several directions
are regression controls, not an estimate of a solid angle.

At `y=(1,2,4)`, the plus preimage is face index `2`; the minus preimages
are `0,3,6,8,16,17,18`, all with negative Cartesian sign. The additional
directions `(-3,5,2)` and `(7,-2,3)` give distinct preimage lists and the
same respective signed totals. The manifest records all sixty coefficient
triples per branch. The ambient orientation is fixed once, opposite the
Cartesian orientation, making the convex branch degree `+1`; it makes the
conjugate degree `+7` without a branch-dependent sign correction.

It follows by exhaustive evaluation of this two-element family that
degree one, trace `phi` for the marked `g`, Laplacian cost `5-s`, and
abstract edges being shortest chords are equivalent. Their conjugate
diagnostics are likewise equivalent. This is **equivalent within the two
Galois-conjugate rank-three realizations** specified here, with their fixed
orientation and generator conventions. It is not a classification of all
embeddings. In particular, applying a common orientation-reversing ambient
isometry flips both signed degrees while leaving both Grams unchanged.
A Gram alone does not choose that orientation. The trace diagnostic also
depends on the marked order-five class.

## A1 support attachment audit

| Object | Supplied data | Missing identification |
|---|---|---|
| Local twelve-port boundary `K_i` | `BoundaryPacket`, oriented triangles | Its response rays are not the federation nerve. |
| Carrier-map rays | `ScreenCarrierMapCandidate.candidateRayZ` and its Gram/rotation theorems | `PortCarrierCandidate` fixes the positive Gram; it does not consume `FederationSupportShadow`. |
| Barycentric local refinement | `DiscreteRefinement.Barycentric`, `baryCarrier_refine_meshRay` | Same-parent denominator scaling is not the global nerve/support bridge. |
| Federation overlap nerve `N_r` | Declared carrier overlaps and designated cycle `z_r` | No source-bound map carries the local fundamental cycle to `z_r` while identifying the radial rays. |
| Global spherical support `S_r` | A1.7 midpoint refinement and oriented `iota_r:|S_r| -> S^2` | An icosahedral seed does not identify `iota_r` with the local response frame. |
| A1 bridge `b_r:N_r -> S_r` | Source-bound degree-one map and refinement square | No theorem equates its oriented spherical realization with the local Gram face map. |

`source_selection_model/geometry.py` verifies the global subdivision chain
and degree-one coarsening; sections 1 and 4 of `DERIVATION.md` use `N_n=S_n`
and `b_n=id`. Each carrier has a separate local twelve-port boundary and
six-dimensional response. The positive geometric support chart in that
construction does not remove its conjugate local response component.
`PR-53` in `docs/PREMISE_REGISTER_V3.md` records no physical discharge
from carrier maps, finite cone constructions or supplied refinements. Its
finite geometric receipts are not a theorem identifying the A1 bridge with
these local Gram rays. `PortGramRepairBand` and the covariance limit supply
no such map either.

The absent interface is named `PORT-GRAM-SUPPORT-ATTACHMENT`. A sufficient
precise content, without a new axiom, is source-bound maps
`h_r:|K_i| -> |N_r|` with `(h_r)_*[K_i]=z_r`, together with a proof that
the local radial face map is homotopic to `iota_r |b_r| h_r` as maps to the
same oriented sphere. The comparison must pin the branch's source rays,
the fundamental cycles and orientation convention, and commute with the
relevant refinement data. Exact equality is stronger than necessary;
oriented homotopy plus the cycle identification suffices for degree.

Under that interface A1 would force degree one, which excludes the degree
seven member of the frozen pair. This conditional reasoning is not a
source selector in this packet: no source theorem inhabiting the interface
was found, and the reference has `source_theorem=null`. No interface is
added to the three axioms and PR-53 is neither consumed nor discharged.

## Integrated source negative control

`galois_source_response_control.py` pins `response.json`, its independent
checker, `geometry.py` and `DERIVATION.md` by SHA-256 at the baseline.
It applies field conjugation to every real and imaginary field coefficient
of the response tape while fixing the complex unit. Port indices, ordered
words, abstract adjacency and the global support construction are held fixed.
The independent formula is recomputed with the conjugate `v_p`, including
the symmetric quadrupole and its norm term.

Both replays check response rank twelve; observable rank twelve on `M_6`;
rank eleven on block-diagonal observables; six even symmetric directions;
six independent odd skew directions; 66 ordered brackets; 20 executed
Cayley factors; 60 proper actions; 720 naturality pairs; and all 3600 group
compositions. Skew-adjointness and injectivity give positive trace pairing
separately for each real source. The filled symmetric/skew tangent and the
declared primitive grammar give the same completeness argument. No scalar
direction is removed from the public tangent.

The A1 boundary and global support clauses do not depend on this exchange.
The carrier algebra, central record, normalized scalar seam restrictions,
trace-reset and measure/prepare instruments have the same definitions.
They are covariant under either proper unitary family by trace cyclicity
and unitary conjugation. Constant-family refinement intertwines either
identical descendant response; its instrument/refinement square is linear.
The A2 response-naturality and chart-composition equations are replayed
exactly above. The all-level support and instrument interpretations use
the same analytic arguments in sections 1, 3 and 4 of `DERIVATION.md`.

These are invariances of the stated clauses under replacement of the
constructed source law. They do not assert equality of arbitrary observed
probabilities between the two sources, or Galois positivity of arbitrary
states/effects. Literal exchange of the two matrix blocks does not equal
Galois conjugation: only one block carries the imaginary quadrupole.
The unchanged-quadrupole mutation fails the conjugate-source formula gate.
The original verifier's hard-coded positive chart is a fixture admission
check, not an A1/A2 selection theorem. The control keeps the abstract
adjacency fixed instead of testing edges against an un-conjugated `phi`.

The bounded conclusion is that response completeness plus same-response
proper holonomy does not by itself select `G_plus` over `G_minus` in this
audited construction. This is not a classification of every response source,
a re-verification of all record/causal programs, or a physical simulator
receipt. The independent global support continues to have degree one.

## Separate entropy root

`code/P_derivation/paper_math.py` computes its `phi` as `(1+sqrt(5))/2`;
`fine_structure_fixed_point_demo.py` does the same. Their README attributes
that root to an entropy balance and describes separate transport and
normalization inputs. This audit establishes the shared real-embedding
convention with `V_plus`, conditional on any geometric selection of that
branch. It does not merge those derivations, certify the entropy/transport
chain, run the numerical demo, or infer a physical numerical prediction.
Those files are not mathematical inputs of either frame executable.

## Replay and trust boundaries

```
python3 -B code/a5_closure/verify_galois_port_frame_independent.py
python3 -B code/a5_closure/galois_source_response_control.py
python3 -B -m pytest -q -p no:cacheprovider code/a5_closure/tests/test_galois_port_frame_certificate.py
(cd Lean && lake build PortGramRepairBand)
(cd Lean && lake build OPHScreen)
python3 tools/check_lean_docstring_style.py
python3 tools/check_axiom_consistency.py
```

Run the last two commands from the repository root. The producer emits JSON
to stdout; the checked reference is `manifests/galois_port_frame_reference.json`.
The independent verifier never imports the producer; its small AST gate
checks the producer's declared dependency boundary without execution.
This targeted gate is not a proof of noninterference for arbitrary Python.

Hostile tests cover a fixed square-root sign, missing branch, swapped
abstract adjacency, reversed face orientation, omitted/duplicated face,
relabeled chord, wrong-order generator, Laplacian sign swap, selected-Gram
input/import, unsupported support theorem and downstream numerical target.
They also cover forged degree, missing coordinate conjugation, false trace
and a boundary regular value. Lean proves the compact Gram and eigenvalue
identities, the marked permutation action and trace pair; degree, PSD and
rank are executable/mathematical results, not Lean claims in this module.

Additional controls reverse the live face input before production, reflect
the vectors while preserving the Gram, and replace the marked generator by
its square. These test the orientation and conjugacy-class scope directly.

No paper, core axiom, historical `PortFrameGram.g5` definition, invariant-mining
freeze or mandatory-runner registration is changed by this packet.

## Validation record

The Lean declarations are hosted in `Lean/Screen/PortGramRepairBand.lean` under
`OPH.GaloisPortFrames`. The independent frame replay and the complete source-response exchange
control pass. The packet's initial 21-test run passes; after adding the
three orientation/class scope controls, the 22 frame tests pass with the
two unchanged source tests deselected. This covers all 24 distinct tests.
The six McKay regression tests also pass. `lake build PortGramRepairBand`
and `lake build OPHScreen` pass, the latter with 8494 dependency/build jobs
including cache replays. Each of the twelve theorems in the module has an
explicit axiom receipt containing only standard Lean axioms.

`check_axiom_consistency.py`, `check_claim_registry.py`,
`build_selection_ledger.py --check`, `check_reader_style.py` and
`git diff --check` pass. Lean docstring style and explicit axiom coverage
were also checked for the tracked modules.

The full mandatory suite, a Lean topological-degree formalization, a positive
source attachment selector and numerical physics calculations are not
attempted. These local results do not establish remote CI validation.
