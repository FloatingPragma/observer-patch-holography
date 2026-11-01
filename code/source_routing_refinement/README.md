# M1: routing refinement and live local storage

This companion derives three requirements of the supplied source-read program:
a support-distance bound at every triangular refinement level, a local working
store with a fixed number of scalar slots per carrier, and a quantitative
obstruction to assigning the logical population's weight to every raw operation.
It does **not** derive M1's read/feedback law from the canonical source.

The [contract](CONTRACT.md) specifies the program, evidence and scientific
boundary. Routing components are formalized separately and composed
analytically; the kernel proves a 23-slot carrier layout, while the sharper
occupied-host census is enforced on the finite controls.

## 1. The routing bound is constructive

The simulator subdivides a triangle `(a,b,c)` into the three corner triangles
and the central midpoint triangle. All four children intersect pairwise in
a midpoint. Every old vertex survives in a corner child. Consequently:

1. Any two children of one face are separated by at most one fine edge.
2. A shared vertex of two coarse faces gives a fine edge between suitable
   children of those faces.
3. To lift a coarse path of `d` edges, move within each parent to the required
   child, cross the next parent boundary, and finally move to the requested
   target child. This takes at most `2d+1` edges.

The actual twenty base faces have diameter three. `base_paths` checks all
base pairs with kernel `decide`, without compiler-trusted reduction.
`stencil_refinement` proves the two lifting hypotheses for arbitrary vertex
labels and the literal four-child stencil. The analytic composition of
`base_diameter`, `stencil_refinement`,
`tower_diameter` and `envelope_closed` gives `D_l <= 4*2^l-1` for the
combinatorial subdivision. These four components are separately kernel
checked. Lean has no theorem identifying `baseEdge` with the base
`triangleEdge`, or `subdividedEdge` with the next level's `triangleEdge`,
and no theorem stating the composed all-level conclusion. No measured
asymptotic exponent is fitted.

The Python reconstruction uses integer vertex IDs and shared-edge midpoint IDs.
It exactly reproduces every face and undirected glued carrier pair in the
pinned L3, L4 and L5 captures. It checks all parent bridges through level five.
These finite comparisons connect the abstract construction to those captures;
the Lean files do not verify the simulator or its floating-point geometry.
Port assignments, carrier placement, read menus and the physical realization
of finer supports remain separate inputs. In particular, the finite degree
checks in the receipt are not an all-level geometric port-assignment proof.

There are `C_l=20*4^l` faces. Choosing the smallest level with `n <= C_l`
gives `n <= C_l < 4n` when `n>20`: the preceding level has fewer than `n`
faces. The deliberately oversized q=3 control support is not assigned this
minimal-capacity property. At q=13 and q=21 the bounds are respectively 63
and 127 hops, exceeding the inherited measured maxima 44 and 91.

## 2. Logical immutability does not require permanent physical allocation

The existing implementation gives every captured relay and every committed
layer version a fresh register ID. Its archive therefore grows with the entire
history. The program's actual reads are narrower:

- A layer consumes only its preceding payload layer.
- One source's multicast tree and all its recipient additions finish before
  the next source starts.
- Each pruned tree visits a carrier at most once, excluding its root.

Thus two payload banks and one relay slot per carrier suffice. All live logical
registers have distinct physical slots under this layout:

| Local slots | Use | Allocation |
|---|---|---|
| 0–11 | Mutable transport ports | Every carrier |
| 12 | Protected zero | Every carrier |
| 13–18 | Six protected address scalars | Occupied hosts |
| 19 | Mutable accumulation | Occupied hosts |
| 20–21 | Payload layer modulo two | Occupied hosts |
| 22 | Current source's captured relay | Every carrier as needed |

For `C` carriers and `n` occupied hosts, `storage.py` enforces the
`14C+9n` occupied-host census on the finite controls. Lean's
`occupied_slot_bound` proves the arithmetic inequality
`n <= C -> 14C+9n <= 23C`; it does not formalize occupancy. The kernel
lifetime predicate `Live` includes addresses, accumulators and both payload
banks on every carrier, yielding the `23C` layout. These are sufficient
bounds, with no necessity or optimality claim for two banks.
The native program finishes all old-layer reads before any commits;
its stronger barrier may allow one payload bank. We retain two banks because
the proof supports simultaneous live versions in two adjacent layers.

`SourceRoutingStorage.lean` defines those logical versions and their live set,
proves slot injectivity on the live set and preservation of the logical store
under writes and retirement. The abstract value type is arbitrary, so a value/writer
pair can be used without assuming that equal scalar values identify a writer.
The proof consumes the stated lifetime discipline. Its attachment to the
native loop uses code review and exhaustive finite tape checks, not a
kernel-verified interpreter. The retained binary format uses signed 64-bit
half-unit values; it is not an arbitrary-precision realization for every
possible input or refinement level.

`storage.py` supplies an allocation certificate for actual event tapes. A
separate last-use pass rejects reuse before a version's final admitted read.
At every event, the source and current physical writer must agree exactly.
Only register addresses change: every event, owner, consumed-writer ID and
scalar value remains in the mapped tape. A second replay uses only physical
addresses/current cell contents and checks local operand/destination kinds,
writer identities, the declared carrier bound and arithmetic. The retained
compressed tape is decoded and replayed, with an independent event-by-event
comparison of opcodes, values, owners, consumed writers and operand presence.
Transport ports must retain their original carrier and local port index, so
a store-consistent port permutation cannot change the captured wiring.
Preparation, permitted commit interventions,
seam validity, the complete read menu and induced logical order are checked by
the existing independent source oracle. The target replay does not claim to
replace that source oracle.

All four q=3 controls (baseline, source intervention, branch intervention and
scratch intervention) contain 22,818 events each. All **91,272** mapped events
are retained and replayed. Each source control checks all 6,561 logical pairs.
Each mapped control uses 16,919 physical slots against the bound 18,163, and
reuses 863 of its 17,782 original logical allocations. These are exact slot
counts, not measured process-memory savings.

The same formula gives the following conditional production bounds:

| Family | Carriers | Sites | Archived registers in inherited run | Derived working-slot upper bound |
|---|---:|---:|---:|---:|
| q=13 | 5,120 | 2,197 | 3,065,436 | 91,453 |
| q=21 | 20,480 | 9,261 | 46,674,948 | 370,069 |

These rows apply the proof to the pinned program/counters. They are not
q=13/q=21 mapped executions or measurements of their peak memory.

The following costs are outside the scalar-slot bound:

- Exact scalar precision, which can grow with layers and read fan-in.
- Writer IDs and tags; these require enough bits for the represented history.
- BFS queues, parent/recipient tables, schedule and read-menu control.
- The retained audit tape, semantic-ancestry analysis and verifier memory.

The offline allocator deliberately retains its input and last-use information;
it is not a 23-slot verifier. Old physical bindings may be discarded only
because the admitted program will not read them again. A later request for an
old retired writer is outside this interface and is rejected. The audit log
continues to retain history separately.

## 3. Raw operations do not inherit logical event mass

Let `K` be the number of rounds, `R` the total logical reads, and `H` the total
multicast tree hops. Counting the unchanged event program gives

```
E = 13C + 8n + 2Kn + R + 6H
N = (K+1)n
Kn <= R <= Kn²
H <= Kn(C-1)
```

Preparation contributes `13C+8n`; layer starts and commits contribute `2Kn`;
each logical read adds one accumulation, and each hop uses six transport
events. Self reads give the lower read bound. A pruned source tree captures
each nonroot carrier at most once, giving the hop bound without charging
shared prefixes repeatedly. These imply both the explicit event upper bound
in Lean and `E >= 3N + 13C + 5n`.

If a separately established logical population estimate gives
`|wN-V| <= epsilon` for nonnegative `w`, then the **same** weight on all raw
operations satisfies

```
|wE-V| >= 2V - 3 epsilon.
```

This bound can be vacuous for large error; it is an obstruction to transferring
a vanishing logical-count error at positive volume to the complete raw
population. It concerns the whole execution window. It does not compute the
population of every interval in the primitive semantic order, disprove all
possible weighting schemes, or select a physical measure. In particular, we
do not introduce inverse multiplicity weights and call them source-derived.

Similarly, if a serial layer has `cost>1` operations and a positive primitive
tick `tau` must fit in declared layer duration `Delta`, then
`cost*tau <= Delta` implies `tau < Delta`. The serial tick cannot simply be
identified with that logical duration. Parallel implementations require their
own scheduling/dependency analysis; neither model supplies a physical clock.

## 4. Supplied premises and scientific boundary

| Input | What this work establishes | Separate source-derivation obligation |
|---|---|---|
| Protected versions | Unbounded archival allocation is unnecessary for the admitted live reads | A physical local cell retaining an exact value/writer through unrelated interactions, and its precision/capacity cost |
| Feedback | Existing operations survive safe slot reuse | Selection/realization of local reset, export, capture and accumulation from the canonical source law |
| Support and placement | All-level combinatorial route bound; three exact captured-support matches | Source selection of those supports, hosts, port pairing and address/control transport |
| Read program | Complete finite menu and intervention behavior are preserved | Why that metric neighborhood and preceding-layer rule arise from source dynamics |
| Scheduling | Program lifetimes and every serial event are accounted for | Source-selected controller, fairness and physical time calibration |
| Population and measure | Added operations cannot silently inherit the logical event weight | Which primitive or distinguished events carry physical volume, with a same-world limit theorem |

The conditional six-operation hop and history theorem assumes the supplied
feedback and routing law. The closed pair-mean load invariant excludes resets
that change total load on a closed register set. An ancillary implementation
must account for changed ancillary state and reuse; supplied fresh or reset
ancillas do not derive the missing local law. Scalar precision, writer
metadata, controller state and retained audit history lie outside the
scalar-slot bound. The raw-mass result concerns the whole execution window;
it neither rules out every weighting nor supplies a physical clock or measure.

## Reproduction

From the repository root with the pinned Python requirements installed:

```sh
python code/source_routing_refinement/verify.py
python -m pytest -q code/source_routing_refinement/test_refinement.py
cd Lean
lake build Geometry.SourceRoutingRefinementAxiomAudit
```

`verify.py --write` regenerates the receipt and all four complete target tapes.
Ordinary verification only reads and compares them. The committed receipt
binds the implementation, Lean sources, captured supports, inherited run
receipts and complete decoded target tapes by SHA-256. It also pins the
inherited producer, source oracle, verifier and codecs used at this boundary.
Bare imports in those inherited scripts are scoped to their actual modules
and restored afterward; duplicate JSON keys and nonintegral numbers are
rejected. The module audit rejects
all axioms except `propext`, `Classical.choice`, and `Quot.sound`, with expected
failures for `sorryAx` and compiler-trusted reduction. Python controls enforce
coverage of every theorem; they do not substitute for running Lean.

The Python checks run in the dedicated Linux
`source-routing-refinement.yml` workflow. Lean CI imports the axiom audit via
`Geometry.lean`. The mandatory runner remains byte-identical to its registered
source freeze; this work does not refresh or weaken that separate freeze.
