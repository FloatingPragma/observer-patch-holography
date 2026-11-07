# The unread native channel and a coherent clock

The reduced native repair process cannot reproduce an arbitrary-state
unitary clock merely by choosing an encoder and decoder. This statement
concerns **unread carrier output**: the decoder receives no repair archive,
coherent environment, noiseless bypass memory or intermediate corrective
service. It does not exclude an observer using those resources, nor an
observer reading the simulator specification. A bound on every input is
also stronger than a comparison for one prepared state.

## The exact source mixing estimate

Use the binary 66-state invariant sector of the
[native instrument](../interacting_extension/source_selection/PROOF.md),
with independent uniform seam selection, fair placement coins, and either
the un-tilted C or E readout. Sum all repair outcomes to obtain `Phi`.
The 78-state carrier with transient double loads is not substituted into
this theorem, and an activity tilt is a different channel.

`source_contract_audit` reconstructs every signed-difference block from
the pinned native event columns. There is one population block of size 66,
132 coherence blocks of size 10, and 2970 of size 1: all 4356 matrix units are
present. The integer blocks, divided by 60, are nonnegative and symmetric.
Their exact Gershgorin lower bounds are 40/60 for C and 22/60 for E. The
coherence block row sums are at most 50/60 for C and 42/60 for E.

The population characteristic polynomial is checked exactly. Its eigenvalue
one is simple; its largest other eigenvalue is

```
r = (55+sqrt(5))/60.
```

Exact real-root isolation places the other polynomial factors below 57/60.
The smaller root of the exceptional quadratic is also below `r`. Positivity
of all quantum blocks and the coherence row-sum estimates therefore prove
`r` is the spectral radius on the full traceless operator space, not merely
on classical populations. The sole stationary density is `I_66/66`.

Let `R` replace its input by this density. Self-adjointness in the
Hilbert–Schmidt inner product gives `||Phi^n-R||_(2->2)<=r^n`. On a pure
input entangled with an arbitrary reference, its Schmidt support has
dimension at most 66. The output thus acts on dimension at most 66 squared;
`||X||_1<=66||X||_2`. The Hilbert–Schmidt bound is unchanged by tensoring
with an identity, and mixed inputs follow by convexity. Hence

```
half_diamond(Phi^n,R) <= min(1,33 r^n).
```

For `k` independently evolving carriers, channel telescoping gives the
upper bound `min(1,33 k r^n)` to their product replacement channel. Their
initial encoded input may be entangled. Independence here is a statement
about the product dynamics, not about that input. No claim is made for
coupled carriers or a common adaptive controller.

## The obstruction for every static encoding

Let `E` encode a `d`-dimensional logical system into those carriers and `D`
decode only their final reduced output. Both may be arbitrary CPTP maps.
They may not transmit any input-dependent side channel around the repair.
Their composition with replacement is another replacement `R_tau`, and
channel contractivity gives

```
half_diamond(D (Phi^n)^tensor(k) E, R_tau) <= min(1,33 k r^n).
```

Apply the target unitary channel `U` to half of a maximally entangled state
`Omega_d`. Its output is a pure maximally entangled state. Replacement
instead gives `tau tensor I_d/d`, whose probability in the projector onto
that pure state is `1/d²`, independently of `tau` and `U`. This two-outcome
measurement alone proves

```
half_diamond(R_tau,U) >= 1-1/d².
```

The triangle inequality gives the advertised exact lower bound

```
half_diamond(D (Phi^n)^tensor(k) E, U)
    >= max(0,1-1/d²-33 k r^n),       d>=2.
```

There is no need to simulate 66 to the power of the number of copies.
Since `sqrt(5)<9/4`, the rational replacement `r<229/240` supplies an
outward lower bound. For a logical qubit, one native carrier and 200
attempts, the error is greater than 0.74 for every such encoder, decoder
and unitary. The lower bound is zero when the estimate is inconclusive;
zero does not certify a successful encoding.

## What a logarithm or a reference can supply

The native stationary density has `-log(I_66/66)=log(66) I`, so its modular
automorphism is the identity. A unitary coordinate change cannot turn that
scalar into the nontrivial Higgs Hamiltonian. A normalized isometric
compression of a tracial density is again tracial. Nontracial ancillas or
selected conditional references add further preparation data.

The positive native **superoperator** instead has a nontrivial Euclidean
logarithm on the 4356-dimensional operator Hilbert space. Calling
`-log(Phi)` a Hamiltonian and applying `exp(-it H)` changes the evolution:
the original process is the dissipative channel `Phi^n`. The bound above
excludes identifying these two arbitrary-state processes by static coding
under its stated output restriction. It does not forbid implementing the
Euclidean logarithm as a separately controlled coherent program.

The [fixed Gibbs clock identity](../common_source/THERMAL_CLOCK.md) is a
positive conditional route. On the full matrix algebra, equality of
modular and Hamiltonian flows with `t=-beta s` holds exactly when the
faithful reference is `sigma=exp(-beta H)/Tr(exp(-beta H))`. That identity
does not derive the reference from the unconditioned tracial repair state.
The existing same-history scalar oscillator likewise has its declared
source pulse controls and is not already a common clock for every Higgs
sector.

Active coherent or record-assisted correction remains a separate route.
The tests explicitly exhibit both a noiseless-memory bypass and recovery
from a retained random-unitary record; each defeats the reduced-output
restriction and so escapes the theorem. A source-to-field implementation
must specify those resources and compare its complete process with the
field process. No claim that such an implementation is impossible, or that
physical OPH dynamics are dissipative on their full accessible algebra,
follows here.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/particle_identification/test_passive_clock.py
```
