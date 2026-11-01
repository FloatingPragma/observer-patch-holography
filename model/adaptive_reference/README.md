# Native feedback with a stationary finite controller

This closes a finite CP-02 scheduler premise that the earlier
`continuing_reference` construction left open. The new policy genuinely
changes the next repair according to the current raw native record. Its
two classical working bits and the quantum field have a stationary joint
marginal. The complete growing archive is retained and is **not** stationary.
The original repair, its zero slots, and its outside-sector outcome are unchanged.

## Balanced memory law

Let the working memory be `(c,d)` in `{0,1}^2`. Selector `c=0` requests the
complete C repair on native seam `(0,1)`; `c=1` requests the complete E repair
on `(0,5)`. Both use the unchanged charged-doublet full-Fock extension: all
native `L_r` on the admitted carrier, together with `Q=I-P` on its complement.
The common 17-slot alphabet pads C's unused raw slots with zero operators.
For a native record define `f(r)=coin XOR (raw mod 2)` and for `Q` define `f=0`.
After appending `r`, update

```
(c,d) -> (d XOR f(r), c).
```

The inverse given the retained record is `(u,v)->(v,u XOR f(r))`. This
is a supplied reversible **classical** controller law, executed exactly by
`controller.py`; no quantum pulse compilation of this memory law is claimed.
The admitted memory algebra is classical `C^4`, with arbitrary quantum field
states and arbitrary untouched reference systems. Controller coherences are
not silently admitted or silently measured.

For each output memory `(u,v)`, every incoming Kraus term uses selector `v`.
Thus, if each complete instrument is unital and commutes with a field weight
`sigma`, the corresponding output block is

```
sum_r L_r^(v) (sigma/4) L_r^(v)* = sigma/4.
```

Outgoing completeness follows separately from `sum L_r*L_r=I`. This proves
stationarity on the entire proper field, not just a sampled population.
More generally, on any finite cyclic controller alphabet, the update
`(c,d)->(d+f(c,r) mod q,c)` preserves `sigma tensor I_(q^2)/q^2` under the same
premises. Each tick may also be preceded or followed by a separate unitary
layer selected by the **current** working memory, using any original
number-preserving field Hamiltonian and real duration. The stationarity
theorem survives these layers; the finite mixing theorem below does not
include them.

The prior actual literal source gives `sigma_M=diag(a,b)^tensor M`, with
`a,b>0`, `b>a`, `gamma=log(b/a)` and trace `(a+b)^M`. Hence the executed proper
reference is `sigma_M tensor I4/4`, never a success-conditioned density.
Every initial local source-export failure is kept. A branch containing any
such failure freezes both field and memory and appends a distinct bypass
symbol. The complete field-plus-abort/controller marginal is stationary
after tracing source readouts, seed records, and subsequent archive records.

The modular generator on the successful CQ algebra is
`-gamma N-M log(a)+log(4)`. Its automorphism is implemented by
`exp(i t gamma N)` in every memory sector. Each complete adaptive instrument
intertwines that charge flow, including its output archive; the earlier
literal clock-word bounds apply. This remains a charge clock, scalar on each
fixed-number sector, and supplies no physical energy-temperature identity.

## Exact native initialization, finite memory, growing archive

Initialize a separate native simulator spectator in the explicitly supplied
binary total-two configuration occupied at ports 2 and 3. Apply two native
events on seam `(0,1)`. The spectator never changes, each raw current is zero,
and the two native coins give exactly four equiprobable working states.
Their records remain in the archive. No floating angle or postselection
enters this initialization. Preparation of that spectator as an interacting
field state by a hidden source word is **not** claimed: its deterministic
native basis initialization is an explicit supplied simulator input.

The two working bits are reused reversibly; the archive is not reused.
After `n` proper ticks there are `17^n` formal repair words, including zeros,
plus the repeated bypass word on initial abort sectors. Five bits per tick
suffice to store the 18-symbol union alphabet. Seed raw records require at
most six more bits. The field preparation still has all `6^M` local readout
records. A local-region bound traces complete outside outputs before local
proper-code projection, retaining the earlier volume-independent bounds.
The full literal preparation error does not accumulate under ideal balanced
adaptive evolution. Any implemented per-tick complete-isometry errors
`epsilon_j` add at most `2 sum epsilon_j` in output trace norm; no hardware
accuracy is assumed or inferred.

## Exact quantum mixing in a finite repair corner

Restrict only this additional theorem to the same two highest-spin charged
defects on the three ports `(0,1,5)`. There are three positions and four
classical memory states; the observable algebra is `direct_sum_4 M3`.
No original joined-field Hamiltonian step is inserted: that Hamiltonian
leaks this corner. Every E outcome has rank at most one. The exact integer
superoperator is `S/2` on 36 matrix entries; its population restriction is
`A/2` on 12 populations. Let `S_C,A_C` keep only columns with selector zero.
The packet checks

```
S^8-S_C^8 = diagonal-population embedding of (A^8-A_C^8),
min_ij (A^8-A_C^8)_ij = 1.
```

These are precisely histories that encounter an E event. Each is a retained
rank-one Kraus history, so the count proves complete positivity, including
off-diagonal field inputs, of

```
Phi^8 - (3/64) D,       D(X)=Tr(X) I12/12.
```

Consequently `Phi^8=(3/64)D+(61/64)R` with `R` CPTP and unital. The unique
stationary CQ density is `I12/12`; every other spectral value has modulus at
most `(61/64)^(1/8)`. For any untouched quantum reference and any admitted
joint density, the **archive-traced** error to `I12/12 tensor rho_reference`
contracts by `(61/64)^floor(n/8)` times its initial trace-norm error.
This is a bound, not a fitted decay rate. It neither asserts full-Fock mixing
nor erases information in the output archive. Exact reverse-current checks
recover the initial corner basis state and memory from each final state and
its full record, making that distinction operational.

## Negative controls and replay

A one-bit update `c->c XOR f(r)` is reversible for every fixed record but
fails joint balance. With the same native instruments, a uniform six-state
input maps to populations `(3,1,2,1,3,2)/12`, trace distance `1/6` from the
initial uniform state. Thus reversibility alone is insufficient. A fixed
zero memory is likewise not the stationary uniform seed. Omission of outside
outcomes, conditioning on successful readout, and identifying full archive
stationarity would all change the theorem.

```
.venv/bin/python -m pytest -q -W error model/adaptive_reference
.venv/bin/python -m model.adaptive_reference --check model/adaptive_reference/results.json
```

The receipt binds this exact implementation, source definitions, prior
literal-reference receipt, proof fragment, and tests. No canonical physical
premise, integrated CP gate, QFT completion, or Einstein limit is promoted.
