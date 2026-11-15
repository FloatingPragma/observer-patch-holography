# Coherent readout at the observer's first native current

This construction joins the [native current counter](../native_clock/README.md)
and the [native coherent measurement](../native_quantum_history/README.md) in
**one executed observer history**. The phase readout occurs after the first
moving current, rather than after a supplied number of native attempts.
The same source process supplies both the event used as the clock tick and
the quantum evolution being measured.

The exact first-current readout predicts:

| Native candidate | Plus probability | Minus probability |
| --- | --- | --- |
| C | `883/1365 = 0.6468864469...` | `482/1365` |
| E | `5261/8855 = 0.5941276115...` | `3594/8855` |

The plus-probability difference is `3644/69069 = 0.0527588354...`.
Both candidates have **the same current waiting distribution**, with mean
`33/5` native attempts for this preparation. Their phase measurements
distinguish them. An independent clock, a supplied Hamiltonian and a conversion
to seconds are unnecessary for this conditional relational law.

The carrier, uniform seam scheduler, native candidate, coherent preparation,
phase detector and logical stopping controller are declared. No laboratory
measurement is identified with this protocol. These are exact laws
of the executable native models, not new empirical agreement claims.

## What executes

The preparation is the existing uniform mixture of ten coherent phase pairs
on the first seam `(u,v)`:

```text
rho_+ = (1/20) sum_(w outside {u,v})
            (|uw>+|vw>)(<uw|+<vw|).
```

Each attempt executes the actual native C or E instrument. Its visible outcome
is the absolute native current `b=|j|`, including zero. Seam, fair coin and the
remaining candidate repair label are probability-summed within each outcome.
After the first `b=1`, a retained flag chooses the identity operation for every
remaining controller slot. The latched flag stays one, and no further native
attempt is charged. A zero after a latched one is an impossible branch, which
is explicitly retained with zero density and probability.

The finite program allows **six native attempts** and then executes the complete
plus/minus phase instrument. It keeps all **128 flag/phase records**, including
the unresolved all-zero flag history and its two phase outcomes. Both phase
preparations are run under both native candidates. Every terminal density and
probability is compared with an exact source calculation; native attempt costs
are checked branch by branch. A checkpoint after two controller slots is
serialized, independently replayed and continued.

The first-current limit is calculated exactly by summing its convergent CP
series. It is **not** claimed to have been executed with a finite attempt
budget. The six-attempt unresolved tail is nonzero and remains in the receipt.
Conditioning on the successful finite trials would produce a different
measurement law. Tests also check a separate two-attempt program with an
entangled reference qubit against all complete 132-dimensional output states.
The six-attempt retained experiment has no external reference.

This is native channel-semantics execution with a declared logical halt
controller. It emits no compiled source pulse words and establishes no
physical implementation of the controller or apparatus.

## Exact first-current channel

Let `T0` and `T1` be the zero-current and moving-current CP maps. Every native
zero-current Kraus operator is diagonal: it leaves classical configurations
in place but can attenuate their coherences. The producer counts all source
branches and obtains the exact Schur multiplier

```text
T0(rho)_xy = W_xy rho_xy,        W_xy = Wcount_xy/60.
```

`W` is a Gram matrix of diagonal Kraus coefficients. Its **entries**, not
its matrix spectral norm, control these Schur powers. The native counts give
`0<=W_xy<=13/15`, with diagonal values `5/6` or `13/15`. Consequently, for every
density input, including an untouched entangled reference,

```text
Pr(no current through N attempts) <= (13/15)^N.
```

The exact first-current channel is the convergent CP sum

```text
S(rho) = sum_(n>=0) T1(T0^n(rho))
       = T1(R entrywise times rho),       R_xy = 1/(1-W_xy).
```

The resolvent is entrywise; it is not the matrix inverse `(I-W)^-1`.
Every recorded native Kraus map is checked to be a partial injection, so its
effect is diagonal. All 66 basis-state eventual-current trace identities are
checked exactly. Along with the convergent CP sum and vanishing no-event tail,
this establishes a complete trace-preserving stopped instrument on the
finite carrier. This reasoning does not require fractional native updates,
a continuous CP semigroup or a physical-time interpretation.

For a finite bound `N`, the producer separately calculates the first-current
sum through `N` and the remaining map `T0^N`. Their total trace is one.
Individual first-event branches use `T1 T0^(n-1)`, which is compared with
the state actually retained by the latched runtime.

## Why the two laws differ

The moving-current CP maps are exactly identical for C and E. On a moving
branch the current alone determines the endpoint imbalance. Their
zero-current instruments retain different committed information, so their
coherence attenuation differs before the first current arrives.

For the ten common-spectator phase pairs, the moving degrees have
multiplicities `{3:2, 4:4, 5:4}`. Direct native branch enumeration verifies
this census. Writing `kappa_C=1/6` and `kappa_E=3/10`, the remaining phase
coherence at the first current is

```text
r_L = (1/10) sum_d multiplicity(d) d/(60 kappa_L+d),
r_C = 401/1365,       r_E = 1667/8855,
p_L(plus) = (1+r_L)/2.
```

These formulas are checked against the complete exact stopped density,
rather than used as input to its construction. Reversing the initial phase
swaps the plus/minus probabilities. Native no-current outcomes cannot be
replaced with identity operations: they contribute measurable dephasing.
Simply conditioning on a move on the first attempt gives plus probability
`67/92` for both candidates, a different protocol that loses this distinction.

The diagonal preparation has class law `(2/5,1/2,1/10)`, exactly the native
current counter's stationary event-conditioned, or Palm, law. Here it is
obtained directly from the declared phase preparation; no unexecuted event
is postselected. The imported current-clock calculation therefore gives,
for `n>=0`,

```text
Pr(first current at n+1)
    = (4/75)(13/15)^n + (1/10)(5/6)^n,
Pr(first current at n+1, phase +/-)
    = [Pr(first current at n+1) +/- s_L(n)]/2,
s_L(n) = (1/600) sum_d multiplicity(d) d
                         (1-kappa_L-d/60)^n.
```

The receipt checks these joint laws through sixteen attempts and retains
their no-current tail. The first-attempt joint probabilities coincide:
`(67/600,1/24)`. Later joint phase outcomes distinguish C and E even though
their clock marginals continue to agree.

## What this removes, and what remains

If an independently inserted **identity** idle has probability `1-q` per
controller slot, then `T0'=(1-q)Id+q T0` and `T1'=q T1`. The factor `q`
cancels exactly in the entrywise stopped resolvent. Thus the entire
first-current quantum instrument, including phase outputs and reference
correlations, is invariant under this global idling. The producer checks the
cancellation on all 4,356 matrix units at `q=1/2`; the algebra proves it for
every `0<q<=1`. Attempt costs and zero-record histories change.
This control does not equate native no-current dephasing with identity idling,
and it does not remove arbitrary changes to seam weights.

The prediction is therefore relational: read one declared observable when
another observable of the same process first changes. It removes the chosen
attempt count from the ideal readout law. It does not derive physical proper
time, select C over E, supply the coherent apparatus, or establish spacetime
and interacting field recovery. The prior finite count and quantum-transfer
theorems remain at their stated domains; this packet adds their explicit
stopped-instrument composition and execution, with no new Lean theorem claimed.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.native_relational_readout --check model/native_relational_readout/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/native_relational_readout
```

The [receipt](results.json) binds the actual native code and both joined
packages. `--write PATH` creates a new receipt without overwriting an existing
one. Exact probabilities use rational strings; floating terminal-state
diagnostics use tolerance `2e-11` and are not outward-rounding certificates.
