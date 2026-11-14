# M10.2 proposed: native integer-current closure subgate

Proposed task ID, unallocated upstream. One concern: the carrier's literal native
current rule and integer record identity. No physical cell moves on this subgate.

## Task planning

Candidate: the native integer endpoint law retained in the proposed ring/H_m
construction. For d = n_a - n_b and coin in {0,1},
j = [d - (2 coin - 1)(abs(d) mod 2)] // 2. Update
n_a -= j, n_b += j, E_ab += j with div E positive at a.
Supplied inputs: integer endpoint registers, oriented integer seam register,
coin alphabet, arithmetic and local event service. No mass or charge anchor.

Gate: for every a,b in [-32,32] and both coins, exact integer current,
conserved total load and invariant n + div E; sum of the two currents equals d.
Control: adding one to the recorded current must violate the identity.
This is a complete local arithmetic gate; it does not identify an electron or
prove a topological physical unit. Row 1 remains open until those requirements
are met on the same candidate.

## Commands

```sh
python openwave/xperiments/m10_oph/research/scripts/m10_2_native_charge.py --prepare
python openwave/xperiments/m10_oph/research/scripts/m10_2_native_charge.py --run
```

No random seed is used; the finite domain is enumerated. Expected cost below one
second. Save environment, source commit, stdout hash, elapsed time and all failures
in an author-side receipt. The local preparatory execution is already known to
the author and cannot be relabelled as a fresh pre-registered physical result.
For a cell-claiming successor, file a new held-out gate before executing it.

## Definition of done

Independent derivation of the integer identities, a method note mapping equations
to code, complete run record and upstream review. The full physical row-1 gate is
in the 31-row packet, remains blocked and must not be inferred from this subgate.
