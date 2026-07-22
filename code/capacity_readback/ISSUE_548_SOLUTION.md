# Issue #548: Source-Derived Public-Checkpoint Capacity Packet

## Result

The first fixed-cutoff physical packet is the oriented edge-center register

```text
X_reach = P_12 x {write, check},
D = |X_reach| = 24.
```

The word *physical* is used in the repository's typed fixed-cutoff sense: the record atoms and checkpoint laws are generated from the declared finite screen carrier rather than arbitrary `record_i` labels. The artifact does not claim a laboratory realization, a capacity-indexed cosmic family, a unique finite-size slack zero, horizon saturation, or an electroweak load identification.

The exact certificate emits

```text
|Omega_tilde(r,D)| = 1,
|X_pub(q)| = |X_reach(q)| = 24,
|K(q)| = 40,
|<support(K)>| = 40,
E(G_q) = empty,
M_0(q) = alpha(G_q) = 24 = D.
```

## Acceptance mapping

| Acceptance item | Executable receipt |
|---|---|
| Frozen carrier; unclosed trials; output-blind membership | `build_terminal_fiber_manifest`: exact world plus all 30 single-edge deletions, 24 single-slot deletions, and 12 inverse fixed-point faults. `is_terminal_world` reads only structural carrier fields. |
| Complete terminal-fiber manifest | 67 fully materialized declared trials, SHA-256 constructor/candidate/manifest receipts, exactly one terminal ID, and `terminal_fiber_complete=true`. |
| Observer/interface atoms and total readouts | 12 observers, 24 local atoms each, 30 interfaces, 24 interface atoms each, total endpoint readout maps. |
| Endogenous histories and preregistered publicness | One semantic propagation/check/commit history per public section; universal twelve-port policy frozen before evaluation. |
| Complete joint kernels and support compositions | Full `D5 x C2_antipodal x C2_orientation` family of 40 global permutation kernels and a 40-by-40 composition table. |
| Local-marginal checks | Marginals are derived from every global row and compared against independently emitted observer packets: 11,520 row checks. |
| Injective reversible generators and exact model count | Every continuation is deterministic, injective, and surjective; the CSP backend counts exactly 24 public sections. |
| Compound graph, MIS, exact decoders | Empty 24-vertex graph, complete 24-vertex independent set, and inverse decoder for every continuation. |
| Approximate branch and TV robustness | Exact worst-input success 1; for rowwise TV distance `delta`, the same decoder gives success at least `1-delta`, hence `M_delta=24` by the carrier upper bound. |
| Carrier bound and rank-one saturation | Sparse projections `P_x=|i><i|`, 276 pairwise orthogonality checks, rank sum 24, identity resolution. |
| Empty/incomplete/ambiguous/singleton fibers | Explicit `classify_terminal_fiber` controls for all four cases. |
| Required controls | Isomorphic relabeling; cyclic permutation; same-marginal/different-joint coupling; tiny full-support noise; circular-definition rejection; target taint; identity family; erasure family; equal-finite-suffix nonpromotion. |
| Extension and refinement injections | Separate 24-to-48 capacity extension and fixed-24 refinement embeddings, both checked by `no_new_confusability`, plus negative controls that deliberately add a new edge and are rejected. |

## Why the CSP change is required

The previous global-section evaluator enumerated the Cartesian product of local atom sets. The source carrier has 24 local atoms at each of 12 observers, so that strategy starts from `24^12` candidates. `public_record_csp.py` performs the same exact enumeration with early interface propagation and a most-constrained-variable order. The existing finite tests verify extensional equality with Cartesian enumeration on a generic noninjective record diagram.

## Reproduction

```bash
cd code/capacity_readback
python3 source_derived_public_checkpoint_packet.py --output-dir runtime
python3 -m pytest -q
```

Expected test result for the supplied patch bundle:

```text
21 passed
```

The generated JSON files are:

```text
runtime/source_derived_terminal_fiber_manifest.json
runtime/source_derived_public_checkpoint_packet.json
runtime/source_derived_public_checkpoint_certificate.json
```

## Remaining theorem boundary

This issue closes the first source-derived reversible packet at frozen `D=24`. It does not supply the source family `D -> packet(D)` or the exact slack law selecting one regulator-stable cosmic dimension. Consequently it does not by itself establish the universe-level equation `N=log M_0(U_N)` at a selected `N`.
