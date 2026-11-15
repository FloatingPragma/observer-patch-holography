# Interacting observations on growing three-dimensional boxes

This runnable construction extends a supplied interacting quantum field over
27, 125 and 343 sites, executes its local operations through the existing
twelve-port source compiler, and retains its complete measurement/feedback
histories. Exact comparison bounds control evolution-step error and the
infinite-volume limit. The model is a repulsive **nonrelativistic lattice
Bose field at fixed spacing**. A relativistic continuum description and
selection of this action or geometry by native repair are open.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.spatial_interaction --check model/spatial_interaction/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/spatial_interaction/tests
```

The [derivation](DERIVATION.md) proves the comparison bounds and identifies
the existing mathematics used. The [receipt](results.json) retains all
1,056 native/quantum histories at each size, independent reference results,
matched controls, exact rational certificates and source-resource counts.
The [manifest](manifest.json) binds the implementation, tests and evidence.

## Action and observation dictionary

On an open cubic graph Lambda_R={-R,...,R}^3, choose

    H = m N + J sum_<xy> (a_x^dagger-a_y^dagger)(a_x-a_y)
        + (U/2) sum_x n_x(n_x-1),
    m=1, J=1/16, U=1.

This is a standard Bose--Hubbard family, not a new action; see
[Jaksch et al.](https://arxiv.org/abs/cond-mat/9805329).
Space, couplings, preparation, instrument locations, code routing and
Hamiltonian model time are specified apparatus. Their operational meaning
is fixed before the comparisons:

| Implementation object | Field/record meaning |
| --- | --- |
| Two occupation bits at each site | Code states 0, 1, 2; unused local state 3 has specified padding action |
| Basis pair i<=j | Normalized two-boson Fock configuration; i=j is a doublon |
| Local source-compiled bond word | Kinetic propagation on one nearest-neighbor bond |
| Local source-compiled contact phase | Energy U when both particles occupy the same site |
| Native first current j1 | Selects the relative preparation phase through the bit j1>0 |
| Pointer bit b | Binary outcome of the local projector onto (D+iS)/sqrt(2) |
| Committed first pointer bit | Selects a pi/2 phase on S before further evolution |
| Last pointer bit | Routes the final native repair to F for one, E for zero |
| Evolution parameter t | Supplied Hamiltonian model coordinate; no physical clock calibration |

Here D places two particles at the origin and S places one at the origin
and one at its +x neighbor. A source word prepares
(D+(-1)^(j1>0) S)/sqrt(2) from the supplied doublon. It evolves for t=1/16,
records a Ramsey outcome, applies the conditional phase, evolves for another
t=1/16 and records the second outcome. The disturbed field state continues
through both reads. Both pointer outcomes and both fair native repair coins
remain in the history; there is no branch normalization or postselection.
The native attachment is classical on the specified 66-state input prior.
It does not assert coherent realization on arbitrary native superpositions.

Total particle number two is an **exact invariant sector** of the action
and each completed local operation. Its dimensions are 378, 7,875 and
58,996. This is not a finite occupation approximation: no accessible
Hamiltonian branch can create a third particle. Intermediate compiler
pulses may visit local padding; all sixteen local code states are retained
while each word executes. Only completed words are compressed to the
invariant sector. Arbitrary creation operators or mid-word measurements
are outside this stated interface.

## Retained results

At 16 symmetric steps per evolution interval:

| Sites | Interaction-deletion history TV | Feedback-deletion history TV | Source/reference history TV | Certified source-to-infinite history TV bound |
| ---: | ---: | ---: | ---: | ---: |
| 27 | 0.0328241 | 0.2835851 | 2.70e-9 | 0.625010 |
| 125 | 0.0326510 | 0.2815759 | 2.71e-9 | 0.097667 |
| 343 | 0.0326510 | 0.2815759 | 2.69e-9 | 0.0101825 |

TV denotes total variation of the **complete two-bit law**, maximized over
the two native preparation phases. The native attachment is a common
record-controlled channel, so the same bounds apply to its complete joint
law. The controls delete only U or the bit-conditioned field phase. A
separate control disconnects the final field-bit/native-routing wire and
changes the native port-zero mean by approximately 0.0009953 at 343 sites.

The independent numerical reference exponentiates the complete sparse
Hamiltonian; it does not execute the source's split schedule. The analytic
temporal TV bound for the two-segment history is 125/12582912, about
9.94e-6. The finite/infinite-volume bound is
2[(24J+U)T]^R/R! with total model duration T=1/8. The last column adds
these two errors, capped at one. These are exact rational bounds for the
**exact-arithmetic ideal compilation**; angle-synthesis/compiler roundoff,
floating contractions and physical service errors are excluded. Float
comparisons and their 2e-8 guard are reproducibility
checks, not interval certificates or physical gate-precision guarantees.

A separate interaction witness makes the nontriviality claim rigorous.
Prepare (D+S)/sqrt(2), evolve for t=1/80 and read the same Ramsey projector,
without the intermediate read/feedback. Its interacting-minus-free
probability contrast is approximately 0.0062498. The continuous-Hamiltonian
contrast is at least 11/2400 in every box and on the infinite lattice.
After both 16-step splitting errors, the ideal executed contrast remains
at least 1802221/393216000, approximately 0.00458328. This positive bound
comes from the analytic response theorem, not a fitted numerical signal.

## What is characterized, and what is supplied

The infinite lattice has one bosonic species, conserved particle number,
repulsive on-site contact interaction and one-particle dispersion

    E(k) = m + 2J sum_a (1-cos(k_a)).

It has cubic rather than Lorentz symmetry and quadratic low-momentum
dispersion above m. The full-Fock extension satisfies H>=mN, with a unique
vacuum and vacuum spectral gap exactly m, uniformly in finite open boxes
and at the infinite-volume spectral threshold. This gap is supplied by m;
fixed-two-particle measurements cannot detect the common phase 2m and do
not measure that gap. The execution establishes the N=2 interface, not
arbitrary full-Fock instruments or a finite-density interacting vacuum.

The completed milestone is a spatially extended interacting construction
whose records admit a controlled lattice-field description. It supplies a
concrete step toward the [effective QFT target](../CONSTRUCTION.md), while
relativistic propagation, a controlled spatial continuum/EFT domain,
physical clocks and comparison with nature are additional work. It is
a separate spatial branch sharing the carrier controls, not yet the spatial
extension of the existing charged matter/flux/scalar common action.
Reproducing the deliberately compiled target is an implementation and
comparison result; it is not an independent empirical prediction.

## Source execution and resource accounting

Every bond is compiled in a four-qubit occupation code using the existing
[twelve-port source compiler](../observer/source_scalar/compiler.py). The
executor obtains the local matrices by contracting the emitted source
pulses, then checks them against independent target matrices on **every**
local input, including padding. The target matrices never replace source
execution. Preparation, feedback and both binary pointer instruments also
use emitted source words. Blank pointers, initial doublon preparation,
coherent code load/export, register routing and pointer readout are declared services. Ideal proper-code service-failure probability is zero
by the inherited compiler identities; no nonzero outcome is discarded.

Six disjoint edge matchings make each symmetric step local and give a
norm bound independent of volume. The simulator caches identical local
source contractions and applies them to exact Fock amplitudes. This reduces
simulation work, not the requested physical operation count. At 343 sites,
the unoptimized two-read history requests up to **1,320,713,826 source
pulses**, with their load/export service counts retained; it declares 687
lattice/pointer qubits. The source ledger excludes costs of initialization,
storage, physical readout and the two classical native repairs. No efficient
hardware or autonomous physical implementation is inferred from runtime.

Guarded polar corrections act only on already executed local matrices.
Roundoff-sized entries forbidden by the exact sector or diagonal identities
are compressed away; full-phase residuals and the onsite off-diagonal
residual are retained as numerical diagnostics. There is no normalization
of individual state or outcome branches. The bounds describe the ideal
compiled construction and keep this floating representation separate.
