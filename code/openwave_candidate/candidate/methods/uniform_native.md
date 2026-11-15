# Exact reduction of the complete local-coarse instrument

The reduction is exact for a separately specified uniform modular-record
preparation. It does not replace the definite-zero-record preparation of F1.
It includes each native outcome, its state, the no-commit branch, per-seam
GKSL evolution and arbitrary external references. It assumes one load, a
connected finite graph, Z_N registers with N>=2, and no record energy or other
operator beyond the declared port/transit generator and local-coarse C/E lift.

Let the original graph have P ports and M seams. Its subdivision has P+M
positions and 2M half-seams, with incidence A oriented from port to midpoint.
For the fixed Gauss residue g=delta_initial, define

```text
Fiber_z = {F in Z_N^(2M): A F = g-delta_z};
beta = 2M-(P+M)+1 = M-P+1.
```

Every fiber has N^beta elements, including composite N=6. To prove this,
choose a spanning tree of the subdivided graph. Assign all beta chord
registers arbitrarily. Remove tree leaves one at a time: each leaf equation
fixes its incident tree register uniquely, with coefficient +1 or -1.
The final equation is satisfied because both sides sum to zero modulo N.
No division by a nonunit in Z_N occurs. Thus the construction is a bijection
from Z_N^beta to each fiber. Closed walks of fundamental cycles generate its
kernel, so this is the complete reachable fixed-Gauss sector, not an arbitrary
subset of its record states.

Define an isometry from bare position space:

```text
V |z> = N^(-beta/2) sum_(F in Fiber_z) |z,F>.
```

The supports of its columns have different positions and are disjoint; their
norms are one. Each half-hop adds or subtracts one record coordinate and is
a bijection between the corresponding fibers. Its sign is a fixed incidence
sign, independent of F. The offset m is constant on each position type.
Consequently `H_full V = V H_bare`, and exponential evolution intertwines.

A complete native hop changes the two registers of its seam by (+j,-j),
modulo N. This is also a bijection between its source and destination fibers.
The C and E outcome classes depend only on the drawn seam, coin, law and
source position, never on F. Each outcome is an injective partial translation
on its position class. Let K_i be the full local-coarse Kraus operator and
k_i the same operator acting on bare positions. Then, separately for EVERY
outcome i,

```text
K_i V = V k_i.
```

The factor for fair coin and equal C/E mixture is 1/2 in each matrix; a
uniform seam choice adds 1/sqrt(M). The supplied transit rule treats every
transit position as elsewhere, without measuring its identity or records.
This is the existing local-coarse instrument, with no stronger readout added.

The equality holds before summing outcomes. It therefore preserves the full
classical-quantum ledger, conditional branches and repeated channels. The
no-commit branch follows from the Hamiltonian intertwiner. For any reference R,
replace V by V tensor I_R and each Kraus operator by K_i tensor I_R. The same
identity proves the extension for arbitrary entangled input/reference states.

For a rate r assigned to each seam,

```text
L_full(V rho V^dagger) = V L_bare(rho) V^dagger;
L = -i[H, .] + r sum_e (C_e-I).
```

Exponentiation gives the entire continuous semigroup on this invariant image.
`generator_action` executes the reduced action without forming a large
Liouville matrix. `evolve_per_seam` implements finite-time per-seam evolution through a matrix-free
LinearOperator and exponential action, with forward and adjoint operations and
an analytic supertrace. It assembles no dense Liouville matrix. The exponential
action is numerical; no certified numerical propagation error is claimed. The canonical run in this module instead declares its
Bernoulli split explicitly; its per-seam first-order rate is p/(M*t), and
is not silently identified with a different graph's rate at the same p.

The state assumption is substantive. A definite F=0 state at the prepared
port is one summand of V|port>, and has squared overlap 1/N^beta with it.
Individual cycle-record measurements or the fully monitored joint-source
instrument generally leave V's image. The rejecting tests demonstrate this.
Global POSITION-only measurement still intertwines, because position
projectors preserve the image; it is a stronger bare channel than the local
coarse instrument. It must not be confused with global position-and-record
monitoring. On trees beta=0, each fiber is a singleton, so this distinction
vanishes. Disconnected graphs require componentwise preparation and empty
fibers at unreachable positions; this API rejects them.

Record readouts are not silently discarded. For a half-edge on a cycle, its
coordinate is uniform on Z_N in a uniform fiber, since a fundamental cycle
can shift that coordinate by a unit. A bridge coordinate is fixed by the
regional charge on its cut. `half_record_marginals` implements both cases and
is checked against the full sector. A general regional boundary sum satisfies

```text
Q_region + boundary(F) = g_region modulo N.
```

`boundary_charge_readout` reports this joint net-flux/charge measurement without
enumerating record fibers. That coarse measurement is a function of position
in the fixed-Gauss image. It is compatible with V. Measuring an individual
cycle record generally disturbs the code and invalidates subsequent reduced
execution unless that disturbance is modeled. Seam histories and independent
winding counts are not recovered from a reduced position state.

Complete Z3 and Z6 triangle sectors verify H, every coarse Kraus operator,
eight-step density evolution, individual record marginals and external-reference
extension. A one-edge control verifies the bridge rule. The canonical
`flat_ring(1,1)` graph has 24 ports, 63 seams, beta=40 and bare dimension 87.
Its formal full sector has `87*6^40` states. Exact reduction permits repeated
native-coarse evolution on 87 positions, with a separately declared coherent
uniform superposition of `6^40` record states at each position. The preparation
is specified by the isometry; no executable finite-time native preparation
procedure is proved here.

The theorem licenses same-construction position and coarse local-record
observables for this named preparation, including compact transit-mode tests.
It does not license F1 definite-record substitution, many-load native balancing,
a physical particle/species or charge-unit identification, a continuum limit,
or an official OpenWave pass. No reserved Fable acceptance case is run.

Implementation: [uniform_native.py](../oph_openwave/uniform_native.py); tests:
[test_uniform_native.py](../tests/test_uniform_native.py). Run
`python -m oph_openwave.uniform_native --output PATH` for a compact source-bound
receipt. Add `--full-ledger` for every step's complete marginal local-outcome
ledger. This does not reconstruct a joint multistep quantum record history.
The canonical per-seam run uses the same r=0.29/(3*0.37) as the triangle,
for duration0.1, with no graph-size dilution. The eight-step Bernoulli run
separately declares p=0.29 and time0.37. The state preparation is a named
variant of the same modular local-coarse dynamics, not a default preparation
change. Record energy is unselected. The fully monitored hard-wall model is
a distinct comparison control.
