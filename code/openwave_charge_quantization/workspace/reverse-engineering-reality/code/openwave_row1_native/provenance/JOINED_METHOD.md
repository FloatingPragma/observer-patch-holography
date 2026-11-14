# Joined native source and winding readout

`joined_readout.py::execute_joined` is one complete execution path. Its
inputs are an actual known one-load preparation G, a source seam schedule,
the fair source coin word, a declared region, an auxiliary pointer in an
unknown one-load position, and its fair coin stream. It accepts no terminal
charge, expected winding or precomputed regional counter.

The actual `native()` function evolves the source. An observer initializes
twelve six-state counters to zero and updates only the two counters touched
by each returned signed current. The regional class is obtained from those
counters. The auxiliary pointer then performs the finite native discovery,
reset and cycle traversal defined in `NATIVE_CYCLE_METHOD.md`.

The source and pointer tapes are joined. Every letter retains the addressed
carrier stage, seam, coin and actual current. The primary source is idle
under the declared addressed scheduler while the terminal meter runs.
The reported quantity belongs to the retained source prefix, whose load
state is unchanged during the pointer phase. Clock counts include both
source and pointer attempts. The inspected source schedules use at most
three attempts and the bounded pointer at most 624, giving a maximum of
627 combined native attempts.

The native counters and the pointer winding use the spatial orientation

\[
w_A=\sum_{v\in A}(G_v-n_v).
\]

The corresponding relative source charge is `q_A=−w_A`. Both values are
reported explicitly on conclusive readout. Every reset or measurement
timeout reports no decoded charge. Zero charge, including the empty and
full region, goes through the same pointer-reset path; its possible reset
failure is retained.

The record replay uses all 456 complete source histories of the frozen
spatial construction. It retains 564 no-move source events and 324 moving
source events. For a source history h of length L, its original probability
is `p_h=2^(-L)`. Let p_x be the exact pointer-reset success probability for
its unknown initial port x, and let `s=1−2^(-R)`. On the bounded one-load
source, the conditional meter success probability is

\[
p_xs^{6|w_A(h)|}.
\]

For every configured source schedule, region and pointer initial column,
the complete joint law is

\[
\begin{aligned}
\Pr(\text{conclusive degree }w)
 &=\sum_{h:w_A(h)=w}p_h p_xs^{6|w|},\\
\Pr(\text{timeout})
 &=\sum_h p_h\bigl(1-p_xs^{6|w_A(h)|}\bigr).
\end{aligned}
\]

The sum is exactly one. The receipt checks 19,008 such probability budgets
over the 132 source schedules, twelve single-port regions and twelve
unknown pointer input columns. These are separate configured protocol
inputs, with no invented prior over them. Fair source coins, all pointer
coins and timeout complements are unconditional within each protocol.
The source-counter theorem separately covers every regional cut.

Every complete fine joined history has probability
`2^(-(source attempts + pointer attempts))`. The finite controller gives
the full stopping-tree semantics without storing an exponential list.
The receipt retains complete representative tapes for both charge signs,
native holds, zero charge, empty and full regions, and each of discovery,
routing and measurement timeout. A source with seven initial loads is
rejected at admission rather than interpreted through the bounded class-one
decoder. A deliberately omitted actual current gives a counter mismatch.

Simulator final loads are retained only for independent checking. The
controller does not read them. Source hashes are checked before and after
the joined run, including the actual native function, geometry, frozen
cycle meter, source history receipt and joined executor itself.

This supplies the complete operational path from the primary native
source's emitted currents to the spatial class readout. The geometric
interpretation is the separately checked dual-cell cohomology map. The
construction uses an additional one-load native pointer, finite observer
memory and an addressed scheduler. It uses no coherent transfer, analog
phase controller or supplied projective detector. Physical electromagnetic
units and platform status do not follow from this operational identity.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge/joined_readout.py
```

The complete output is `JOINED_READOUT.json`.
