# Compact row-1 review payload

`compact_review.py` is the complete standalone executable payload. Run:

```sh
python3 compact_review.py
```

It uses only Python's standard library and prints JSON. It contains no file,
network or process access and imports none of our packages. The thirty seams,
twelve-port domain, floor/ceiling endpoint repair, fair-coin branch alternatives
and swap/hold control are all visible in this one file. `OUTPUT.json` is the
saved output; `RUN.json` records an externally measured runtime and source hash.

This is a **same-group verification of already known results**, independently
implemented for full payload inspection. It is not a new prospectively
registered discovery or candidate. It checks the finite component of
`A-native-multiload-unit-selection-v1`, whose sealed full evidence remains in
[`../2026-11-15_codex_pr723_unit_selection/`](../2026-11-15_codex_pr723_unit_selection/).
That packet's manifest SHA-256 is
`35b67994df0c89d8478624815a95585293b99249da46b4dc90ba5fd27b72d5ae`.

The script enumerates all 1,820 nonnegative integer configurations with total
load 0 through 4, calls the disclosed native rule on all 109,200 seam/coin
branches, and checks binary-class closure. Reverse graph search establishes
that every configuration has a path to the binary class, of at most three
events. In this finite closed state space, IID scheduling that gives every
seam/coin branch positive probability therefore enters that class almost
surely. The script does not sample or claim an observed infinite-time run.

The control check runs all 720 branches from twelve actual double-load states
for each law: native repair enters the binary class on 120; swap/hold on zero.
Swap/hold preserves the multiset of loads, so subsequent control events cannot
remove the double. Both outcomes use the same full-state binary predicate.
A concrete paired event includes its oriented current and input/output state.
All 60 native branches of the uniform-two state remain fixed, explicitly
retaining the boundary to a universal 0/1 claim.

These are full-state mathematical audits. The script does **not** reimplement
the operational current-only terminal observer, replay the original long
histories, verify the unbounded signed-state proof, or inspect the complete
source package. Those claims require the original evidence and separate
checks. With the added double and uniform controls it executes 109,980 native
calls and 720 swap/hold calls; this count includes only literal calls here.

Falsifying outcomes for the finite statement would be a binary exit, a state
with no path to the binary class, or failure of the stated counts. Larger
integer loads and double-preserving control outcomes are representable. The
integer alphabet is supplied, not derived. The uniform-two counterexample and
the need for a density promise remain; this payload establishes neither
physical electric/topological charge nor a row award. It needs no coherent
memory. It makes no energy/composition claim and does not resolve row 2's
separated-pair energy excess. G6 and G9 remain unawarded.
