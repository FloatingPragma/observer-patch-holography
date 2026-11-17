# Charge-preserving field preparation

The common source uses `G=BE-n-r`. Its preparation service corrects electric
consistency while retaining native matter charge and the compensating
reference. It is a separate operation from free evolution under `H_total`.
The exact integer implementation is `audit_core.gauss_repair`.

## Local rule and termination

Supply a connected oriented graph, a root, protected integer charges
`q=n+r` and initial integer flux `E`. The sum of `q` must be zero on the
closed graph; an incompatible total is rejected. Choose a rooted spanning
tree. At a nonroot vertex `v`, with parent `p` and tree edge `e`, let
`g=BE-q`. The local batched correction is

```text
delta E_e = -B[v,e] g_v
g'_v = 0
g'_p = g_p + g_v.
```

Every other syndrome is unchanged. Native charge and every chord flux are
unchanged. The weighted score `sum_v depth(v)|g_v|` decreases by at least
`|g_v|` for a nonidle correction. Leaves-to-root postorder requires at most
`V-1` attempts. Repeated complete nonroot sweeps have at most
`sum_v depth(v)` genuine changes and a conservative attempt bound
`(V-1)sum_v depth(v)`. The root syndrome vanishes by total compatibility.
Integer arithmetic is batched; its bit cost depends on the represented
charges. The bound is not a claim that large binary charges require only
unit-cost arithmetic.

The routine returns the final flux, every nonidle history tuple
`(vertex, edge, removed_syndrome, delta_flux)`, and the attempt count.
The topology checks execute 768 schedules, all retained regional charge/flux
comparisons and a 257-bit charge example. These tests concern integer
constraint preparation, not the evolution or physical geometry of the
full native source.

## Quantum history and finite capacity

With protected matter/reference labels retained, store the removed syndrome
in a private orthogonal register. The map

```text
|matter, reference, E>|blank> -> |matter, reference, E+delta E>|g_v>
```

is invertible on its admitted input labels and therefore defines an
isometry on their span. Tracing out the history can dephase the reduced
source; arbitrary source-only coherence is not claimed to survive that
trace. A physical error-correction interpretation additionally specifies the
error model and syndrome reservoir.

The executable integer repair has no hard-wall cap. A finite implementation
must admit the initial flux, every intermediate flux reconstructed from the
history, and all syndrome labels within its declared resources. It rejects
or records capacity failure otherwise. Reusing the infinite-integer
termination bound at an arbitrary small cap is invalid. The complete
instrument retains unsuccessful preparation and capacity outcomes.

## Work and the physical interval

A repair may increase the physical Hamiltonian energy. Its preparation work
and apparatus history are accounted separately. The tree dressing fixes
chord flux and is not generally the compact interacting ground or the real
Coulomb minimum. A fractional minimum-energy solution is not relabelled as
an integer-flux eigenstate.

During free evolution every source and auxiliary term commutes with Gauss.
An initially physical state therefore needs no interleaved repair. The
charge and flux readouts use the same orientation and retained reference;
zero mismatch probability is a property of the physical preparation,
with all mismatch labels retained in the full instrument.
