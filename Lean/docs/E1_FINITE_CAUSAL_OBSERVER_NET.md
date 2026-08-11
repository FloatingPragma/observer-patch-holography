# E1 finite causal observer-net handoff

## Attained finite theorem packet

`QFT/FiniteCausalObserverNet.lean` defines a proof-carrying enrichment of one
A3 `ConsensusTower`. At each finite regulator it carries:

- a finite region poset with declared overlaps and declared disjointness;
- region-indexed star subalgebras of the tower's private matrix algebra;
- isotony and exact commutation of observables in declared-disjoint regions;
- covariant region and local-algebra refinement;
- compatible regional expectations inherited from the A3 state-pairing law;
- idempotent complex-linear local repair, including explicit fixed-region and
  fixed-disjoint-observable receipts; and
- a generic B4 basis-reindexed partial-trace helper whose factorization and
  Kraus completeness remain explicit inputs.

The B2 relaxation of every local repair has the exact multiplicative law
`R_a R_b = R_(ab)`. A declared-disjoint observable is fixed for every
amplitude, so its expectation against every matrix is unchanged. This is an
exact Heisenberg-observable identity. It is not a CP/CPTP channel theorem.

`QFT/ObserverNetDescent.lean` defines compatible local families by equality
after restriction to each declared overlap. The API type `FiniteCover` is only
a nonempty finite family of declared subregions: the abstract region preorder
has no union or join operation, and the type contains no joint-coverage law.
A declared family has unique restriction gluing exactly when every compatible
family is the restriction of one unique element of the ambient local algebra.
The file constructs the glue, proves its restriction and uniqueness laws, and
proves that unique gluing forces the family restrictions to be jointly
injective.

All singleton identity families have unique restriction gluing. A parameterized
consistency model adapts any supplied A3 constant partition tower and B2
partition average. Its regions are subsets of two labels, it contains two
distinct nonempty declared-disjoint regions, and a two-member subregion family
has exact unique gluing. The code does not provide a closed concrete partition
and state in this file, and it does not formalize the two-member family's union
as a coverage theorem.

## Audit boundary

The contravariant star-homomorphic `restrict` maps, including their retraction
of isotony, are declared extra data. They are not supplied by an ordinary
isotonic AQFT net and are strong for noncommutative local algebras. Their role
is to make the requested overlap-descent statement well typed without hiding
a conditional expectation or split assumption.

The parameterized partition model is deliberately degenerate: every region is
assigned the same commutative public-record algebra and every restriction is
the identity. For each supplied partition and state it proves consistency of
the interface, but it is not a closed source witness, a noncommutative regional
factorization, or a realization of local quantum physics.

The `TensorSplitReceipt` is also weaker than a regional split: it stores a
basis equivalence and declared disjointness, but no field identifies the local
algebras of the two regions with the two matrix factors. Its marginal theorem
is the generic B4 identity after basis reindexing and does not use regional
disjointness. A regional application needs a separate factor-localization
theorem.

The negative controls prove that:

- assigning the full `M_2(C)` algebra to two labels declared disjoint violates
  locality, so finiteness and isotony do not imply the locality receipt;
- an idempotent repair can change the remote unit observable, so repair
  idempotence does not imply no-signalling; and
- a cover with two distinct global sections having identical restrictions
  cannot satisfy unique descent.

The imported B4 control additionally shows that Kraus or kernel normalization
cannot be dropped from the corresponding factorized no-signalling identity.

## Exact open boundary and issue status

This packet does not construct a closed noncommutative, region-separating inhabitant
from observer-patch source data. It does not derive a tensor split from
disjointness, identify regional algebras with the factors in its generic B4
helper, formalize genuine joint coverage, construct CP/CPTP repair, prove
scheduler locality for adaptive repairs, or supply spacetime causality, a
time-slice axiom, continuum QFT, physical clocks, fields, sectors, or
laboratory attachment.

Accordingly, the two requested Lean deliverable files exist and the
finite conditional interface is substantial; the remaining scope was the
noncommutative/source-attached model with coverage and factor receipts.

## The rich-fibre packet (2026-08-11)

The second preregistered bounded run substantially narrows that scope without
closing it. The pilot run `e1_prereg_64k_20260810` failed its gates
(one split fibre per truncation-12 window; overlapping windows) and is
closed negatively under its fail-closed clause. The second
preregistration froze two pilot-calibrated design changes before its
fresh-seed run: truncation 20 and target-blind greedy-disjoint observer
selection (support node identities only). On the fresh seed every gate
passed: split-fibre counts 3, 3, 4, 4 on pairwise-disjoint windows,
byte-replayed and custody-pinned
(`oph-physics-sim` `docs/E1_RICH_FIBRE_PAYLOAD_2.json`).

`QFT/RichFibreWitness.lean` mirrors the payload literally with
kernel-decided receipts (window lengths, pairwise node disjointness,
split census, per-fibre positions, ambient indexing on the 81-point
basis). `QFT/RichFibreRegionalNet.lean` builds the tower and the
`FiniteCausalObserverNet` instance over it and earns, with standard
axioms only:

- a genuinely noncommutative block algebra at **every** observer region
  (`richRegional_noncommutative_all`), with a constructed two-by-two
  matrix-unit corner associated with each designated split fibre
  (`richDesignatedFactor`);
- region separation: distinct windows are declared disjoint with the
  kernel node-disjointness receipt backing the declaration, and block
  members of distinct windows commute elementwise
  (`richBlock_commute_disjoint`, `richFibreNet_separation`);
- the character-and-block restriction system with the reflexivity,
  transitivity, and inclusion laws;
- **computed genuine coverage**: the four-window family generates the
  top algebra (`richWindowCover_coverageLaw`,
  `richWindowCover_familyJoin`), and dropping any one window breaks
  coverage because each window carries its own noncommutative block
  (`richDropCover_not_coverageLaw`).

The enrichment rule, region lattice, restriction system, base-point state,
and identity repair are declared postprocessors of the payload literals;
ambient window disjointness is constructional, and the witness receipt proves
the realized node windows are independently disjoint in the source. The
matrix-unit corner does not instantiate `TensorSplitReceipt` or identify the
regional algebra with a tensor factor. E1 and the B4 region-factor gate remain
open. CP/CPTP repair semantics, scheduler
locality, continuum causal and time-slice structure, and physical
clocks remain with #693, #700, and #703.

## Source-operator adapter preflight

`QFT/SourceOperatorGeneration.lean` mirrors the committed payload
`docs/E1_SOURCE_OPERATOR_PAYLOAD.json` of `oph-physics-sim` (schema v1,
sha256 `98d17ab4683eb967acd3b420ea96c14f6130b7802a61e7ef56561ef8f413ca6c`,
recounted post hoc from the retained gates-passing bundle
`runs/e1_prereg2_64k_20260810`): the 32-step support-modal label paths of
the four rich observers over the five packet fields, their realized
alphabets, counted transition operators, and field tables.  The
generation theorem is kernel-backed reachability: the field projectors
separate the realized states, compressing the counted transition operator
between diagonal units extracts each counted step as a matrix unit, and
the committed walk visits every realized state, so the generated star
algebra equals the full matrix algebra on each observer's alphabet.  The
diagonal-only negative control proves the transition operator
load-bearing.  This replaces the declared pair-groupoid enrichment with
algebras generated by transcribed source operators.

`QFT/JointSlotFactorisation.lean` assembles the designated pair
`(86, 88)` on the declared joint slot space `Fin 13 × Fin 14`.  The
lifted generated algebras are exactly the two tensor factors, the slots
commute, the checkpoint-projector Kraus family built from the committed
field table is complete inside the identified left slot, and its pinching
leaves the right marginal invariant through the committed B4 helper
`ptraceFst_local_kraus`.  The conditional expectation onto the right slot
carries its matrix-unit Kraus identity with positivity, trace,
fixed-slot, left-scalarisation, and idempotence receipts; it is the
restriction replacement for the excluded scalar-character interface.
Step-census decides mirror the payload (4 left-only, 4 right-only, 14
double, 9 unchanged of the 31 aligned steps).

The extraction is post-hoc and ineligible as validation; the step
alignment is the payload's checked conditional chain (shared engine
snapshot list, pinned history window of 32, exactly-32-step counts
forcing zero drops); the slot assembly and label conventions are declared
postprocessors; and the packet constructs no region lattice over the
tower.  The factor-identification corollaries recycle the generation
theorems through the slot embeddings and are not independent results.
Issue #692 stays open on the source-attached net packet.

## Verification

The targeted build order is:

```text
lake env lean -o .lake/build/lib/lean/QFT/FiniteCausalObserverNet.olean QFT/FiniteCausalObserverNet.lean
lake env lean QFT/ObserverNetDescent.lean
lake build QFT.RichFibreWitness QFT.RichFibreRegionalNet
```

The displayed `#print axioms` reports contain only `propext`,
`Classical.choice`, and `Quot.sound` (`richSplit_counts` is axiom-free).
No file declares a project axiom, an admission, an unsafe definition, or a
`native_decide` step. The rich-fibre packet additionally carries the
unique-restriction-gluing receipt on the four-window cover
(`richWindowCover_hasUniqueDescent`) and its reconstruction corollary
(`richWindowCover_reconstruction`) through the committed descent
interface.
