# Fixed physical noise, corrected source reads

Protection operates within an observer-like self-reading patch system:
local state and ports support readback, retained records and corrective
feedback, with public evidence bundles for the complete noisy instruments.

This package builds on the [all-sector fermionic source](../m1_fermionic_source/README.md).
It removes the need for a decreasing physical dephasing rate in that source's
constructed operational refinement. Four levels of error correction, finite
universal synthesis and a paid local control schedule give vanishing errors
in the named reads, clock and positive accounting at fixed finite rate.

The theorem-library premises A1--A3 are unchanged. This result uses the
proposed CCG source and its declared cofinal control capability; it does not promote CCG,
a microscopic source or a noise model into a consequence of A1--A3.
Read the [contract](CONTRACT.md) and the [proof](../../docs/research/FIXED_RATE_SOURCE_READS.md).

## Concrete results

- A complete 64-outcome seven-qubit decoder corrects every single-position
  error, including coherent errors and spectator entanglement.
- A bounded physical recovery uses verified cats, two candidates, four
  syndrome rounds and explicit idle slots. Independent forward and backward
  calculations check all **65,328** one-location Pauli faults on its **19,664**
  circuit locations. Every case is correctable. Removing verification or
  syndrome repetition produces explicit failures.
- Full M6 channels verify the continuous driven dephasing decomposition;
  noise during a noncommuting pulse is included.
- With the credited constructive universal fault-tolerance theorem, the
  complete prepared source experiments have decoded error
  O(q^-12 log^6(2q)) and additional accounting error O(q^-8 log^6(2q)).
  Recovery, software ancillas, physical flights and failures are charged.
- Raw microscopic fidelity can simultaneously vanish. Complete decoding
  preserves the intended accounting effect on correctable inputs, while
  all error diagnoses remain in the source history. An abort stays an output.

## Reproduce

With `PYTHONPATH=code` at the repository root:

```text
python -m m1_fixed_noise.build
python -m m1_fixed_noise.verify
python -m pytest -q code/m1_fixed_noise
```

The producer uses forward symbolic propagation. The independent verifier
propagates output sensitivities backward, constructs code projectors,
replays every decoder branch and evaluates driven channels by a different
matrix-exponential algorithm. A direct adaptive executor checks branch
selection. Mutation controls bypass source hashes to test physical rejection.
An independent contraction of the ideal quantum circuit checks all 256
Kraus maps of each of the six cat instruments, on every active data input.
Deliberately broken preparations, couplings and reads fail that check.
Adaptive replay rejects missing, out-of-range, duplicate and malformed
fault requests rather than silently treating them as a clean execution.

## Scope and costs

The accuracy theorem and Solovay--Kitaev synthesis are credited prior work.
The universal non-Clifford protection is an imported constructive theorem;
the finite Clifford recovery census is not a simulated universal machine.
No numerical threshold is inferred from that census or borrowed from a
different gadget schedule. Universal library constants are symbolic.

The local drive increases by a polylogarithmic factor relative to the
parent's particular pi c/a schedule. Each finite drive, time and inventory
is finite and counted. Parallel encoded flights preserve finite signal speed
and a uniform wall-time rescaling. This is not a fixed-strength device or a
practical hardware estimate. Central classical records retain their existing
reliable status in the specified quantum-noise model.

The theorem compares decoded logical outputs, named application records and
abort flags. All auxiliary records are retained, but their distributions
need not resemble noiseless records. Known preparation is included. A noisy
first encoding of arbitrary unprotected unknown inputs has a separate
unavoidable error and is explicitly not granted for free. Floquet accounting
is a declared positive reference, not physical drive energy.

A lower-level exhausted ancilla supply retains a complete local output and
a failure diagnostic, so an upper code can correct it. Only an outermost
failure aborts the application. Aborting the entire experiment on every
lower-level failure would not obey the proved concatenation bound.
