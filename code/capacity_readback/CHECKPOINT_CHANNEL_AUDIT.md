# Exact public checkpoint capacity

This repair under issue #1033 makes the finite evaluator implement the
existing Pro5 readback definition. It selects no physical source, capacity,
temperature or clock. The mathematical support criterion in
`paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex` already distinguishes exact
zero-error recovery from approximate recovery; the old executable did not
preserve that distinction for sufficiently small positive errors.

## Reproduced defects

The regression-only commit precedes the implementation repair. All thirteen
cases in `test_checkpoint_channel_precision.py` fail on main `7f99c8cc`.

| Input or operation | Old result | Required result |
| --- | --- | --- |
| Binary crossover with error `2^-40` or `2^-60`, tolerance zero | Approximate evaluator returns two codewords | One: both inputs can produce both outputs |
| Binary crossover with error `10^-400` | Conversion to float deletes positive support, returning two codewords | Preserve the overlap; zero-error capacity is one |
| Binary crossover with error `1/8`, tolerance `1/8 - 10^-14` | Acceptance tolerance admits two codewords | Reject the two-word code |
| NaN checkpoint rows in an otherwise valid packet | Empty supports manufacture saturation and `PASS` | Reject the kernel |
| Boolean/string probabilities, duplicate inputs or coerced output labels | Accepted or silently reinterpreted | Reject malformed channels and alphabets |
| Positive transition mass outside a continuation alphabet | Silently removed from support closure | Reject the same-alphabet continuation claim |
| Floating addition for an error-transfer bound | Can round below the exact sum | Return an upper bound |

These are errors in finite evidence, not objections to a valid stochastic
channel model. A further algorithmic defect made the generic decoder search
enumerate exponentially many assignments even for an invertible channel.

## Numerical input and output contract

`checkpoint_channels.py` is the common parser for support graphs, approximate
capacity and continuation relations. It accepts finite `int`, `float`,
`Fraction` and `Decimal` probabilities in `[0,1]`, rejecting Booleans,
strings, invalid numbers and unsupported types. Inputs and outputs are
nonempty string labels; input labels must be distinct. Every channel has
exactly one row per input, and a declared family must be nonempty.

A float denotes its exact represented binary value. A decimal or fraction
retains its exact value. The parser preserves the historical row-mass
acceptance boundary `abs(sum(row)-1) <= 10^-12`, but makes its interpretation
explicit: accepted entries are **relative weights**, normalized by their
exact positive row sum. Every reported decoder discloses the original mass
for every codeword. There is no additional tolerance on the decoding error
or on whether a support entry is positive. Normalization defines the
evaluated kernel; it is not an uncertainty bound on an unknown physical
channel. Supply exact unit-mass fractions when that distinction matters.

Exact errors and epsilon are serialized as canonical numerator/denominator
strings, including after a JSON round trip. Floating success summaries are
lower bounds and floating error summaries are upper bounds. For example,
error `10^-400` stays positive in the exact certificate, has a positive
outward display bound, and is never displayed as perfect success.

The historical private `_channel_rows` adapter remains for the reversible
packet producers. It accepts only losslessly representable normalized
binary64 probabilities; it refuses lossy conversion. General capacity and
support computations bypass that adapter and use exact rows.

## Decoder simplification and correctness

For a fixed code `C` and channel `K`, define

```text
e_x(d) = sum(K(y|x) for y with d(y) != x)
e_C(K) = min_d max_(x in C) e_x(d).
```

The common code is admissible precisely when `e_C(K) <= epsilon` for every
declared channel. Each channel may have its own decoder, as in Pro5's `d_K`.
The channel is known to the decoder. This does not assert the stronger
unknown-channel decoding contract: identity and cyclic shift need different
inverse decoders.

The exact search removes only dominated assignments:

1. If an output is possible for just one codeword, assigning it to that
   codeword decreases its error without increasing another codeword's error.
2. If several codewords can produce an output, assigning it to a codeword
   with zero probability, or leaving it unassigned, is dominated by assigning
   it to any positive-probability owner. All such owners are searched.
3. Partial errors are sums of nonnegative probabilities. They are lower
   bounds on the errors of every completed decoder. A branch whose maximum
   partial error is at least the best completed error cannot improve it.

Thus these reductions preserve the minimax optimum. An explicit stack avoids
a recursion-depth limit. Codes are searched in decreasing size and every
rejected code has a channel with optimal error strictly above epsilon. The
first accepted code therefore has maximum cardinality. A shared allocation
budget raises an error if exhausted; it never returns a smaller code as a
certified maximum. Defaults remain receipt-scale: at most twelve inputs and
100,000 decoder nodes. The input limit can be raised explicitly for simple
families; arbitrary noisy decoding remains combinatorial.

At epsilon zero, two codewords share an output precisely when a deterministic
decoder cannot recover both without error. Consequently the exact search
agrees with the maximum independent set of the union support graph.
Positive rescaling of a row preserves this graph. Arbitrarily small
full-support noise can still destroy every nontrivial zero-error code.

For a reversible family, every positive output has just one owner. Each
channel needs one search node for the full code. The complete retained
24-record, 40-channel source packet now has independently replayable generic
decoder witnesses with a budget of forty nodes; no special capacity formula
is substituted for the generic evaluator.

Support closure also preserves the difference between one checkpoint and
indefinite continuation. For `0 -> 1 -> 2 -> 2`, the one-step capacity is two,
but the squared relation maps every input to `2`. The continuation closure
includes that erasure. Any positive output outside the declared continuation
alphabet is rejected, rather than removed by implicit postselection.

Finally, a fixed decoder's error changes by at most the row total-variation
distance. Exact addition with upward rounding implements the existing
transfer bound `min(1, epsilon+delta)` without narrowing it accidentally.

## Independent evidence and limitations

`verify_checkpoint_decoder.py` imports neither the optimizer nor the parser.
It independently parses the supplied channels, reconstructs normalized
probabilities and errors, checks the full input alphabet, and replays each
decoder. It also checks the exact ratios and outward summaries. It runs with
the producer absent under optimized Python. Thirteen forged-witness controls
include constant `PASS`, false errors, false row masses, altered decoders,
missing channels and a false scope claim.

The replay certificate proves **attainability only**. It does not prove that
a larger code is impossible. Maximality follows from the search argument
above and is independently tested with a flat integer-count implementation
that enumerates all decoders, including unassigned outputs:

- All 325 unordered pairs with replacement of the 25 binary quarter-count
  channels, at five error tolerances: 1,625 capacity comparisons.
- All 216 three-input, three-output half-count channels, at three error
  tolerances: 648 comparisons, also checking monotonicity and the zero-error
  support-graph identity.

Together these are 2,273 exact comparisons, not a proof by testing on arbitrary
channels. The search proof supplies that general finite justification.
Ten isolated implementation mutations restore false threshold slack,
float conversion, support truncation, downward TV rounding, acceptance of
escaping mass, zeroed decoder errors, an always-accepting verifier, omitted
normalization, greedy assignment or ignored invalid rows outside the code. Each
is caught by an ordinary assertion or a missing required exception, not by
an import or collection failure. CI executes these controls and the entire
capacity directory on Linux.

The subsequent maintainer-style audit also compares the actual minimax error
against flat decoder enumeration for 125 pairs of three-input, four-output
channels with unequal row denominators. It tests each pair at the exact full
code threshold and immediately on either side when that remains in `[0,1]`.
A three-channel counterexample has capacity two for every individual channel
but common capacity one: the calculation must select one common code, not
take the minimum of separately optimized code sizes. Other controls exercise
both exact row-mass boundaries, low Decimal precision, invalid rows outside
the returned code, zero-probability external labels, the legacy float
adapter's explicit refusal, and the TV transfer through the capacity caller.
These additional controls found no further production defect in the repair.

## Downstream impact

The corrected generic noisy evaluator can change answers previously accepted
within `10^-12` of a boundary, and exact support now retains probabilities
below float range. Malformed rows that previously manufactured a graph or
packet `PASS` are rejected. The extra exact fields extend the approximate
evaluator's return object; existing field names remain available, with
floating success values now conservatively rounded.

The retained source-derived packet, terminal-fiber manifest and certificate
were regenerated separately using main's evaluator and the repaired
evaluator. All three outputs are byte-identical to each other and to the
tracked receipts. Their respective SHA-256 values are:

```text
packet      3604147f715327abdc96ad56160fa5ff5384af767d0a57c49beae5b92f3ed96a
manifest    34e22b6c554c657ca0573516a8c543c2bb995e581cd072c5592f3a4b36c94f69
certificate bf48d67735360ca3d54d4c2b58bb14f0884bb917cc0f4d8d66e7229a75a3ff93
```

The existing capacity suite also replays the bounded continuation,
scalarization, source-family, carrier and horizon controls. The Pro5 support,
TV-stability and zero-error discontinuity statements already have the required
mathematical scope. No paper, public claim, claim-registry payload, frozen
registration, pinned receipt or mandatory-runner byte changes are needed.
`F_READBACK_SPEC.md` is itself a frozen source input and remains byte-identical;
the numerical interpretation is documented here and linked from the directory
README instead of editing that frozen artifact.
In particular, correcting finite capacity evidence supplies neither the
complete physical source antecedent nor a universe-level value of `N`.

Reproduce the complete affected suite, including isolated mutations, with:

```sh
python -m pytest -q code/capacity_readback -W error
```
