# Second independent joint spectral-continuation review

**Outcome: PASS within the declared finite domain.** This is a distinct Codex worker review of the producer and first reviewer, not an opposite-team Fable audit. The driver imports neither the producer constructor nor its grouped instrument. It uses the archived 384-state Hamiltonian, source labels and references, independent SciPy EVR number-block projectors, and direct endpoint-balancing/intervening-bit CAR parity for every native outcome.

The recomputed interacting ground is −2.7262233938630365 in N=3. Its projector agrees with the actual reference within 9.2e-17. Across 5,916 actual-input spectral rows, independently reconstructed degenerate projectors agree within 1.7e-13; full subnormalized branch matrices agree within 1.8e-15. Complete probabilities, including zeros, agree within 9.6e-15. The selected canonical native protocol is independently reproduced within 1.4e-16.

The actual populated degenerate branch is continued without renormalization. Direct same-H propagation and independently reconstructed complete MIX Kraus maps agree with its archived poststate within 3.6e-16 and every joint native-event probability within 1.4e-15. This checks continuation, rather than terminal probabilities alone.

A separate faithful complex density in an actual N=3 doublet is an explicit apparatus control, not a replacement canonical reference. A different complex unitary within that doublet leaves its reconstructed state invariant within 4.3e-17. Fine-basis measurement followed by forgetting the fine label instead changes the state by trace distance 0.36055513; after the native channel the difference remains 0.34260197. Full coefficient matrices inside degenerate groups are therefore necessary.

These are floating matrix checks without interval enclosures. The review does not independently rebuild the full interacting source, whose construction has its own separate review. It verifies the same archived H/reference, complete grouped number/energy apparatus and actual native continuation. It adds no physical species, spin, vacuum, mass, radiation, continuum identification or official OpenWave award. Every input source/archive snapshot is pinned before and after execution.

From the public bundle root:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python next/C10/review_joint_spectral_continuation.py
```
