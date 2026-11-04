# Full-carrier adaptive preparation

## The specified observer process

The source is the existing twelve-port native total-two carrier. Its binary
subspace has dimension 66. The full carrier, including double occupation at
one port, has dimension 78. The two source instruments are C (coin and current
records) and E (coin and contrast records). Every record has a Kraus operator
equal to an integer partial injection divided by the square root of two.
The common 17-slot raw alphabet includes zero slots and an outside slot
which is zero on the native carrier. The complete outcome sum preserves
trace, and on the binary carrier it also preserves the identity.

Working memory is the classical pair `(c,d)`. At a sweep boundary, `c=0`
selects C on all 30 seams in their source order; `c=1` selects E in the
reverse order. After the first event, its raw record `r` updates memory to
`(d XOR f(r),c)`, with `f` the existing
[balanced controller function](../adaptive_reference/controller.py).
The remaining 29 events use the selected law and order without another
memory update. All records are appended to the archive. The output's second
bit stores the selected law, so no hidden load or unrecorded selector is
needed. A 30-state classical phase counter selects the seam within a sweep.
The laws, ordering and memory rule are supplied simulator specifications.

The controller reads a source record: for C on the first seam, different
binary configurations give feedback probabilities `(1/2,1/2)` or `(0,1)`.
Its effect is not inferred from the field's hidden loads. Given the first
record, the update is invertible. The archive is retained and grows; the
stationary object below is a reduced field/controller state.

## Exact stationarity and support

For an output memory `(u,v)`, each incoming first-record term comes from
memory `(v,u XOR f(r))` and therefore uses the same instrument selected by
`v`. On a uniform field and memory input its output block is

```
sum_r L_r^(v) (I66/264) L_r^(v)* = I66/264.
```

The remaining events are complete unital instruments controlled by the
second memory bit, and preserve the same block. Thus the complete sweep
channel `Phi` on `direct_sum_4 M66` fixes `I264/264`. This uses the existing
[finite incoming-row identity](../../Lean/QFT/AdaptiveReference.lean),
attached to every concrete native branch by exact source enumeration.
The same reasoning applies at every elementary phase. Hence the homogeneous
phase-extended process fixes `I7920/7920` on `direct_sum_120 M66`.

A double occupation splits at its first incident repair. Every sweep visits
every seam, and a binary configuration never leaves the binary carrier.
Complete output support after one sweep is therefore binary for every
full-carrier input. Positivity also removes every row or column outside
this support, including inputs entangled with any untouched reference.
This is an absorbing support statement, with no outcome conditioning.

## Quantum minorization from exact path counts

For the E sweep, propagate all 4290 distinct basis pairs through common
fine-record branches. The computed set is empty at the end of the covering
word. All matrix-unit coefficients are nonnegative, so absence of a pair
means that no individual fine-record history transmits it. Every E-sweep
Kraus history is consequently a basis matrix unit times `2^(-15)`, or zero.
Any longer history containing an E sweep has the same rank-one property:
the other factors are partial injections and cannot create a second
nonzero entry.

Let `A` be the exact integer population numerator of one adaptive sweep,
including all four memory states. Its size is 264 and its denominator is
`2^30`. Both row and column sums are `2^30`. Let `A_C` retain only its
columns whose incoming selector is C. Then

```
B = A^4 - A_C^4
min_ij B_ij = m = 1586432651050269168489450503667712.
```

`B` counts exactly the complete four-sweep histories containing E. Every
counted path is a rank-one Kraus word. For every input/output basis pair,
subtract `m` such words; the remaining Kraus sum is completely positive.
Writing `D(X)=Tr(X) I264/264` on the classical-memory algebra gives

```
Phi^4 - alpha D is completely positive,
alpha = 264 m / 2^120,
3/10 < alpha < 1/3.
```

The exact fraction has decimal value approximately `0.31508380895`.
Classical positivity alone would not establish this quantum inequality;
the fine-history rank-one argument is essential. No eigenvalue calculation,
fit, floating tolerance or spectral-gap premise enters the proof.

Trace preservation and unitality give
`Phi^4 = alpha D + (1-alpha) R`, where `R` is complete positive, trace
preserving and unital. Both `D Phi = D` and `Phi D = D`. For any admitted
field/memory/reference density `rho`, including classical correlations with
memory and arbitrary quantum entanglement with the field,

```
distance((Phi^n tensor id)(rho), I264/264 tensor rho_reference)
 <= (1-alpha)^floor(n/4)
    distance(rho, I264/264 tensor rho_reference).
```

Distance is half trace norm. The initial distance is at most one. The same
bound times the input trace holds for subnormalized branches. For a full
78-state input, use one absorbing sweep and replace `n` by `n-1`.
This proves a unique faithful stationary boundary density and bounds every
other boundary eigenvalue in modulus by `(1-alpha)^(1/4)`. Source binding and
hashes of the full integer arrays are retained in the receipt; replay
reconstructs their entries rather than trusting those hashes as a proof.

The 30-phase process also has a unique stationary density: invariance gives
equal phase weights, each phase block is propagated from the boundary, and
the unique stationary boundary block is uniform. A fixed initial phase
cycles deterministically, however. Thus arbitrary initial phases do not
converge to the uniform-phase state. The contraction theorem compares
states at a common sweep phase, with the archive traced on both sides.

## Mathematical scope

The finite reference, complete instruments, continuing adaptive memory and
finite preparation budget are constructed on one source. The source repair
law is unchanged; the specified scheduling and memory law define the
observer used in this theorem. The admitted memory algebra is classical.
Coherent quantum controllers, additional field Hamiltonian steps and
arbitrary alternative feedback laws are outside the theorem. No records
are discarded from the physical execution, but the growing archive is
absent from the stationary marginal and is not asserted to mix or lose
information. This tracial reference has identity modular automorphism.
No physical time, particle mass, interacting continuum or geometric balance
law follows from the finite contraction.
