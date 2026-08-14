# OPH-native W/Z source packet

This package defines the producer-side W/Z source frontier. It is separate
from the external Standard-Model calculation and records issue numbers only as
provenance pointers for the pinned scientific artifacts.

The committed inventory binds the finite source parents that can be consumed
without a physical W/Z claim:

| Parent | Bound result | Boundary retained |
| --- | --- | --- |
| #565 | Twelve-port carrier, incidence, inverse pairing, A5 frame, and finite refinement maps | No universal carrier selection or electroweak load identification |
| #566 | Finite source-model current algebra and response representation | No continuum or laboratory current, coupling, or kinetic normalization |
| #314 | Rank-fifteen chiral matter pair, finite Spin typing, anomaly checks, and conditional Yukawa channel lines | No scalar source, numerical Yukawa matrices, families, or pole masses |
| #567 | Finite diagonal Z6 form, common kernel, lattices, and screen sector transport | No continuum instanton, theta, monopole, or laboratory attachment |

The rank-three family-band result from #569 is carried as
`conditional_context_only`. Its finite screen-response resolvent is not a
physical matter pole. The finite domain attachment has complex rank 45, but
its declared signed matter operator is exactly the scalar signed operator
tensored with the identity on the 45-dimensional fiber. The extension is
conditional and not source-selected. It therefore supplies no
family-sensitive kinetic term, chiral gauge interaction, Yukawa action, or
physical family selection. The issue-314 spin packet lives on its separate
twelve-port support. No source, domain, or transport bridge attaches it to
the issue-634 local operator packet. Matter-pole identification, seam and
gauge action, physical Spin/locality, refinement, and laboratory-current
attachment are not supplied.

The #609 matter-menu receipt is bound as a negative boundary. It is complete
inside the declared exterior algebra and retains a source-invisible sterile
direct-sum countermodel. The #311 exact flux receipt supplies six twisted
adjacency spectra on the twelve-vertex support and the exact same-support
classical stiffness family `K_k = 5 I - A_k`. The regular Eisenstein blocks
use the positive hex-lattice metric `G = [[2,1],[1,2]]`; the verifier proves
phase isometry and the edge sum-of-squares Hessian exactly. A separately
pinned simulator receipt supplies a vector-spring realization for its own
issue-634 local-domain sector and scalar spectra. The latter uses an 8662-node main
spectral complex and a 1052-node complex from the 2048-carrier source
configuration. No identity bridge joins those spectra to the
twelve-vertex/forty-two-vertex exact flux packet. Neither finite result is a
complete-interface classical-versus-quantum indistinguishability theorem or
an extended-domain no-go. The
#616 scalar countermodels are transitively pinned by the #32 frontier. These
packets do not supply a complete physical census, quantum particle ontology,
or continuum pole.

The pinned #634 receipt supplies a finite causal, section, and local-operator
domain at finite bounded scope. It supplies no continuum Lorentzian spacetime
or quantum-EFT transfer. Such a transfer, or a stated-domain no-go, is not
supplied by this package.

The #32 RG frontier contributes exact per-copy representation indices and the
parametric one-loop gauge law after an explicit imported QFT functional. It
does not select family or scalar multiplicity and emits no interval,
threshold, finite map, Jacobian, term mask, or vector remainder.

The #630 frontier classifies the available renormalizable coefficient spaces
and retains exact two-completion witnesses for the scalar, Yukawa, and
`v_chart` to `v_F` choices. The corresponding required interfaces carry
provenance pointers #636 for the scalar action, #637 for complete Yukawa
matrices, and #638 for the source-to-FJ coordinate map. The integration
interface carries pointer #630.
The bound frontier emits no coefficient assignment, physical scalar action,
or coordinate identification.

The #631 carrier frontier proves a unique positive unital normalized-trace
isomorphism between the one-dimensional invariant screen order-unit line and
the four-copy weak multiplicity center. It does not identify that algebraic
line with one physical load or fix an electroweak normalization.

The pinned #633 clock contract proves interval arithmetic conditional on a
valid SI chart. The final #633 verdict is
`PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE`. That verdict
rests on a bounded schema and source scan plus one named channel-nonuse
experiment on the declared serialized interface. It is not a complete-domain
clock non-identifiability theorem. Dimensionless source verdicts remain admissible.
A GeV W/Z row is forbidden on this interface. Newton-G composition remains a
separate downstream use in #334.

The strongest machine status is:

```text
FINITE_BOUNDARIES_CLASSIFIED__DIMENSIONLESS_SOURCE_NOT_SUPPLIED__PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE
```

`promotion_allowed` is fixed to `false`. The inventory emits no pole
coordinate and no physical-unit value. Each resolved-boundary row hash-binds
an immutable scientific receipt and classifies its exact scope and
physical-promotion boundary. Issue numbers, where retained, are provenance
pointers only; live project state is not an input or authority. No bounded
receipt is promoted to a positive scientific source parent.

## Scientific boundary

The package is a deterministic, non-promoting inventory of the bounded parent
artifacts listed above. Its successful replay establishes path and content
closure for those inputs and enforces the separation between source inputs and
comparison data. It does not construct a physical W/Z source.

A promoting source packet would additionally require a complete coupled-field
action and operator census, physical matter and family attachment, complete
Yukawa matrices and coefficients, an explicit source-to-FJ coordinate map,
interval-by-interval matching and running data, and a source-selected joint law
with pole covariance. It would then have to pass an independent clean-room
substitution into the frozen pole calculator before any external comparison is
mounted. None of those objects is inferred from a successful inventory replay.

On the serialized interface represented here, the physical-unit verdict is
`PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE`; consequently
the package emits no GeV W/Z row.

## Files

- `data/source_parent_policy_v1.json`: trusted role, path, scope, and firewall
  policy.
- `schemas/source_parent_inventory_v1.schema.json`: strict output schema.
- `build_source_parent_inventory.py`: deterministic producer.
- `check_source_parent_inventory.py`: independent resolver. It imports no
  producer code and can replay all eleven allowlisted verifiers.
- `outputs/source_parent_inventory.json`: committed non-promoting frontier
  receipt.
- `code/a5_closure/manifests/{stage4_receipt,clock_unit_verdict,classical_realization_receipt}.json`:
  vendored, hash-bound receipts for the resolved #634, #633, and #311
  boundaries.
- `tests/test_source_parent_inventory.py`: deterministic, tamper, target-leak,
  coordinate, unit, family-boundary, partial-frontier, and DAG mutations.

## Run

Fast local replay:

```bash
python3 build_source_parent_inventory.py --check-byte-exact
python3 check_source_parent_inventory.py --skip-native-verifiers
python3 -m pytest -q tests/test_source_parent_inventory.py
```

Full parent replay:

```bash
python3 run_all.py
```

The full checker executes the production verifiers for #565, #566, #314,
#567, and the seven context packets from #569, #609, #311, #32, #630, #631,
and #633. It does not invoke or modify the #593 calculator.

## Coordinate and unit boundary

`v_chart` and `v_F` are different typed coordinates. The equality receipt is
absent and relabeling is rejected. The declared finite domain permits only
dimensionless coordinates normalized by `E_star`. Its physical-unit row is
`PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE` under the
bounded #633 interface audit.

The current external-validation schemas are pinned only as provisional
consumer shapes. They are not a frozen native interface, and their GeV fields
cannot be populated by this producer.
