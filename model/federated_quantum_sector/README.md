# Native quantum fields on joined observers

The existing [signed-difference theorem](../observer/QUANTUM_MIXING.md)
gives the quantum block formula on any finite connected weighted graph in
each fixed-number sector. This packet attaches it to the actual federations
and verifies individual current outcomes, full tensor support and continuing
joint records. Its added execution covers coherent interventions and exact
readout predictions while preserving arbitrary untouched references.

A complete four-attempt experiment crossing two elementary carriers gives:

| Declared experiment | Candidate C: plus | Candidate E: plus |
| --- | ---: | ---: |
| Coherent plus preparation | `5/8` | `1/2` |
| Opposite preparation phase | `3/8` | `1/2` |
| Phase flip controlled by the first zero-current record | `3/8` | `1/2` |
| Dephasing controlled by that record | `1/2` | `1/2` |

The contrast is produced by the seams inside the two carriers; the crossing
seam acts trivially on the probe coherence, so the cross-carrier content of
the experiment is the support of the superposition, not the joining dynamics.

Both candidates implement the same binary population repair. They differ in
the native quantum record: C resolves the current, while E resolves the
input endpoint imbalance. The graph, chosen schedule, coherent preparation,
interventions and phase analyzer are declared services. These are conditional
model predictions, without a laboratory comparison or a selected physical
clock.

## Full instrument reduction

For a finite graph with `V` ports and `M` seams, use
`H=(C²)⊗V`. The bit at each port is its binary occupation. The actual endpoint
rule supplies local Kraus maps `K_(coin,label)` with entries `0` or `1/√2`;
each map is extended by the identity on every other port. C uses the native
current as its label; E uses the input imbalance. Current remains a visible
outcome for both candidates. Coins and additional labels are summed as CP
alternatives. There is no coherent summation of forgotten outcomes.

The exact local audit checks all sixteen endpoint matrix units, each of the
three current outcomes, and both candidates. It proves Kraus completeness,
unitality, total-number conservation and self-adjointness of the unread
channel on operator space. Every Kraus map preserves the total number, so
each fixed-number sector is invariant. A tensor product of **local**
fixed-two sectors is not invariant across a joining seam.

For a matrix unit `|x><y|`, set `z=x−y` and `Z={v:z_v=0}`. Equality of the
committed record makes each outcome preserve `z`, not just the unread sum.
Consequently arbitrary sequences of native instruments, classical feedback
choosing their seams, and declared diagonal phase or dephasing interventions
can be propagated inside these exact operator blocks. This statement includes
every retained word and the unnormalized terminal state. Tensoring the
matrix-unit identity with an arbitrary reference gives the complete linear
comparison. The blocks are vector spaces of operators, not separate physical
state spaces or independent matrix algebras.

For the independently uniform seam/coin scheduler, the block is

```text
Phi_L[z] = (1−kappa_L(z)) I − L_exclusion(G[Z], common particles)/(2M),
kappa_C(z) = sum_edges |z_a−z_b|/(2M),
kappa_E(z) = sum_edges 1[z_a != z_b]/M.
```

The positive and negative difference sites are fixed. Common particles
exchange only along edges within `Z`. The formula follows by checking the
sixteen endpoint matrix units and summing edges; it does not enumerate an
exponentially large density matrix. The executable additionally compares
independent source enumeration with the formula on nonclassical columns of
all three pinned federation levels.

On a connected graph, nonzero `z` in one fixed-number operator space has
positive killing. There is an essential exception on the whole tensor
algebra: `|vacuum><full|` and its adjoint have constant `z` and are fixed.
Thus fixed-number uniqueness cannot be extended to unrestricted coherent
mixtures of number sectors.

## Exchange-field predictions

The two-port exchange observable is
`A_ab=|1><0|_a ⊗ |0><1|_b ⊗ I_other`. Summing the block's common-particle
configurations cancels the exclusion Laplacian. Self-adjointness then gives
the exact Heisenberg eigenrelation

```text
Phi_C(A_ab) = r_C A_ab,   r_C = 1−(deg(a)+deg(b))/(2M),
Phi_E(A_ab) = r_E A_ab,   r_E = 1−(deg(a)+deg(b)−1[a~b])/M.
```

Prepare the one-particle plus state on `a,b`, with all other ports empty,
and read the complete effects `(I±A_ab±A_ab*)/2` after `n` uniform attempts.
The exact plus probability is `(1+r_L^n)/2`. Other populations need not remain
inside the initial pair; the analyzer retains their fair plus/minus outcomes.
The report computes outward rational intervals by integer powers and exact
division, rather than trusting a floating exponential. On the first joining
seam of level zero, `r_C=104/105` and `r_E=619/630`.

These are quantum-coherence predictions in the native attempt parameter.
The uniform-scheduler forecasts are distinct from the fixed four-seam
experiment below. A nonzero coherence block has zero trace throughout native
record propagation, so native current histories alone cannot distinguish its
preparation phase. An explicitly coherent final measurement can.

## What executes

The shared [primitive observer runtime](../primitive_observer/runtime.py)
executes a source-derived four-seam tape on four ports across two elementary
carriers. Its one-particle support has dimension four. The other 236 ports
are exactly empty and every chosen seam stays within the active support;
the embedding therefore preserves all superpositions and references. Unused
global seams are not silently applied or discarded: this experiment declares
its four-operation schedule.

The initial state is the plus or minus superposition on the two internal
endpoints. Four native current reads follow, with the declared intervention
after the first read. The final analyzer first resolves its pair versus its
complement, measures X on the pair, and supplies a fair outcome on the
complement. This specifies the full measurement instrument, not just effects.
All `3⁴×2=162` visible words and all terminal matrices are retained, including
zero-probability words. A separate rational CP contraction predicts every
matrix entry. Serialized checkpoint continuation is compared with the
uninterrupted history. An entangled reference-qubit preparation verifies
complete outputs rather than probabilities alone.

The reference preparation is supplied through the runtime's joint-input
interface, whose preparation cost is left unpriced. Work-only runs retain
their declared preparation resource, and every branch retains all operation
costs.

Execution is finite channel semantics. It emits no compiled pulse words.
Floating residuals compare the runtime against the exact oracle; they are
not a rigorous hardware or source-pulse implementation error certificate.

## Existing theorem attachments and limits

The local binary factors and actual Kraus completeness instantiate the
finite split in
[`slot_locality_receipts`](../../Lean/EventAlgebra/SlotLocalitySurface.lean).
The theorem gives exact preservation of the remote marginal by a
**nonselective** local operation for every entangled input. Disjoint native
instruments commute outcome by outcome, and the tests check this on every
four-port matrix unit. Selected outcomes may steer a correlated remote
state; they do not satisfy the nonselective conclusion individually.

An edge operation enlarges Heisenberg support only by its endpoints. For a
fixed finite schedule, composing this fact attaches the existing
[dependency-cone theorem](../../Lean/ObserverPatchHolography/Locality/DependencyCone.lean)
to native quantum observables. Record-adaptive choices additionally require
their controller's consultation dependencies. The finite tensor split and
graph support do not supply physical spacelike separation or a relativistic
light cone.

The graph block formula is reused, not a general theorem; no Lean
proof is claimed. This packet provides interacting local quantum
instruments, noncommuting readouts and full finite record predictions. It
does not select a continuum state, compatible refinement dynamics, a
relativistic action, or the correspondence between native attempts and an
observer's physical time. The native CP channel is not replaced by a modular
automorphism or a reconstructed unitary evolution.

```sh
python -m model.federated_quantum_sector --check model/federated_quantum_sector/results.json
python -m pytest -q -W error model/federated_quantum_sector
```
