# A universal biased-core walk from the exact-code classification

This component strengthens the scope of the existing
[biased reset result](../biased_reset/README.md). Its source weights,
native instruments, reference, recovery service and reset policy stay as
declared there. The optimal asymptotic reset rate applies to **arbitrary
history-selected exact all-seam sharp-`q=1` source qubits**, not only the
60 displayed encoders. It does not prove the biased policy optimal for
every finite horizon.

The nearest results are the [native flagged-code theorem](../common_source/ENCODED_MEMORY.md)
and the [mobile exact-code classification](../particle_screening/mobile.py).
The latter's fair-coin walk is not substituted for the biased law. The new
step is that strictly positive endpoint probabilities multiply each retained
effect by a positive scalar, leaving the flagged Knill–Laflamme conditions
unchanged. The proof below therefore transfers the classification itself.

## Domain and proof

At every ready boundary the source is binary total-two H66, and `V` is an
isometry from a logical qubit into a sharp-`q=1` sector of a classical
two-port core `{a,b}`. Exact recovery must work for every complete retained
flag of all 30 native seam instruments, with arbitrary untouched references.
The unknown logical state is not stored in a persistent auxiliary. Coherent
recovery after every attempt is an explicit supplied service. The scheduler
and positive endpoint weights are fixed and source independent; C/E choices
and code choices may depend on the complete past. Native holds are allowed
when declared. Attempts cannot be avoided by halting or postselection.

For either C or E, positive coin factors can be divided out of the scalar
effect conditions. The existing theorem implies scalar compressions of
every occupation `n_j` and zero support on every occupied native edge. On
the sharp sector only source pairs `{a,j}` and `{b,j}` are possible. If `j`
is adjacent to either core port, at most one row survives; its rank-one
positive occupation effect cannot be a nonzero scalar on a qubit. Only
ports nonadjacent to both core ports survive.

For each such safe exterior port write the two source rows as a square
matrix `B_j`. Scalarity gives `B_j* B_j = w_j I2`, where `w_j >= 0`, and
isometry gives `sum_j w_j = 1`. A square matrix with this Gram also obeys
`B_j B_j* = w_j I2`, including `w_j=0`. Each core row therefore has squared
norm `w_j`. Since `V* n_a V` is scalar and its trace is `sum_j w_j = 1`,

```text
V* n_a V = V* n_b V = I2/2.
```

Conversely, safe-port blocks `B_j=sqrt(w_j) U_j` with
`sum_j w_j U_j* Z U_j=0` define exact codes: all native edge occupations
vanish, and the remaining complete effects are linear combinations of
`I,n_u,n_v,n_u n_v` with scalar compression. This includes complex frames
and unequal exterior weights. For a distance-two core there are two safe
ports; for an edge there are four. An antipodal core has none, contradicting
`sum_j w_j=1`. Thus all 66 core choices are covered: 60 admit codes and six
are impossible throughout this class.

For a seam from a core endpoint `a` to a port `b` outside the core, that
neighbor `b` is unoccupied on the code. The endpoint is occupied with
probability one half for every unknown logical state. With seam probability
`s_ab`, repair probability `r`, and positive weights `theta`, the exact move
probability is therefore

```text
P(core -> (core minus {a}) union {b})
    = r s_ab theta_b / [2 (theta_a + theta_b)].
```

All other probability stays at the core. This conditional kernel is
independent of the logical state, any reference, and the history-selected
code. C and E have the same core kernel although their full records differ.
The implementation compares the formula with the actual complete C/E
kernels, keeping all six forbidden destinations.

For the fixed biased law, the existing exact potential `h` solves
`h = 1 + Q h`, with zero on forbidden cores. Its unique maximum is at
`(0,1)`, approximately `203.0253781358` native attempts. Every attempt
lowers conditional expected `h` by one, for every admitted code choice.
Every counted reset raises it by at most `H=max(h)`. Thus
`E[N_n] >= max(0,(n-h(initial))/H)` for adaptive or randomized choices,
including voluntary early resets. Resetting only at the impossible boundary
to `(0,1)` attains asymptotic rate `1/H`, approximately `0.00492549261172`.
The compiled family supplies attainability; arbitrary exact codes
cannot lower this rate. This optimizes label interventions, not recovery
pulses, hardware errors, entropy production or energy.

## A clock can be a separate factor of the same process

Protection need not make the protected logical qubit itself a clock. If the
corrected memory branch obeys `M_h V_c = sqrt(p_h) V_c'`, a clock instrument
with retained Kraus label `j` on another factor obeys

```text
(M_h tensor L_j)(V_c tensor I_C)
    = sqrt(p_h) (V_c' tensor L_j).
```

This preserves the logical subsystem with arbitrary correlations among
logical qubit, clock and untouched reference; the clock's own conditional
measurement backaction is retained. Complete histories compose by induction,
and zero/abort outputs stay present. Memory and clock operations commute
because they act on separate factors, not because the memory is invariant
under its source reference's modular flow. The test attaches a nontrivial
qubit phase/readout instrument to the actual source-executed moving/reset
branches and checks every unnormalized joint output.

This is an operational composition theorem, with a supplied clock factor,
initialization, phase control, readout, shared-event order and storage. Its
test is an algebraic clock control, not a compiled native clock. A
record-inferred native clock could fill this interface only with its actual
preparation and complete archive contract retained. Neither the tensor
identity nor a native-attempt timestamp establishes a relation between
memory transport and clock phase, a Hamiltonian dispersion relation,
physical duration, a stationary product reference, or a particle mass.
The reference-like exterior and autonomous protection gaps are open.

```sh
python -m pytest -q -W error model/biased_code_universality
```
