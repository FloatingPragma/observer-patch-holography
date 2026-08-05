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
finite conditional interface is substantial, but issue #692 should remain
open until at least one noncommutative/source-attached model discharges the
strong restriction/gluing and factor-localization receipts or the issue is
explicitly narrowed to the proof-carrying interface plus its parameterized
degenerate consistency model.

## Verification

The targeted build order is:

```text
lake env lean -o .lake/build/lib/lean/QFT/FiniteCausalObserverNet.olean QFT/FiniteCausalObserverNet.lean
lake env lean QFT/ObserverNetDescent.lean
```

The displayed `#print axioms` reports contain only `propext`,
`Classical.choice`, and `Quot.sound`. Neither file declares a project axiom,
an admission, or an unsafe definition.
