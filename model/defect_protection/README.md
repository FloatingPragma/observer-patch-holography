# A protected exchange charge in an explicit observer memory

This construction joins the [native exchange topology](../exchange_sector/STATISTICS.md)
to the complete C/E repair instrument and an explicit one-bit memory law. It
obtains a nonconstant protected **quantum** charge on the enlarged system.
The charge has not yet been embedded in a spatially localized source defect.

The source is the existing binary total-two carrier, with 66 unordered pairs
of distinct native ports. Its ordered double cover has 132 states. A bit `r`
records which ordering of the occupied pair is used. For every complete raw
branch `h`, let `f_h(x)` be the actual endpoint repair and `p_h(x)` its
transported-order parity. The joint Kraus map is

```
K_h = 1/sqrt(2) sum_x |f_h(x)><x| tensor X^p_h(x).
```

The parity is derived by following the two named tokens through the actual
move. For a seam `(a,b)` with `a<b`, it is one precisely when the moving
token crosses its spectator in the sorted-pair chart. Holding columns have
parity zero. The binary domain is closed under C/E. The construction makes
no claim about collision states or refined carriers. It retains every raw
coin/record slot, including zero and abort slots.

Every branch commutes with `Q=I66 tensor X`. Thus both spectral projectors
`(I +/- Q)/2` are preserved by the complete instrument for any finite
classical adaptive schedule; no independence, mixing or all-seam premise is
needed for this preservation. In the X eigenbasis the two blocks are exactly
the existing positive and negative exchange-character instruments. This
does not promise preservation of coherence between the blocks.

The memory is coupled to transport. Starting from pair `(0,11)`, the actual
three-move history on seams `(0,5),(0,11),(5,11)` returns the source and flips
the computational orientation bit. That complete history has probability
`1/8` on its initial pair. The full experiment includes its other records;
the history is not assigned unit probability. Its X eigenstates acquire
characters `+1` and `-1`. Changing the ordering chart at any source vertex
conjugates each map by `G=sum_x |x><x| tensor X^g(x)`; the closed-loop action
and Q are unchanged.

There is also an exact classification for the **all-seam** branch algebra.
The native effects separate all 66 source addresses, hence generate
`P_x tensor I`. A connected tree of native movements gives transport matrix
units. The odd loop gives `P_x tensor X` at one vertex, and the tree transports
it to every vertex. Products therefore give every `E_xy tensor I` and
`E_xy tensor X`. All generators commute with Q, so the generated algebra is
exactly `M66 tensor span{I,X}`, or `M66 direct_sum M66`, of complex dimension
8712. Its commutant is exactly `span{I,Q}`. For uniform seam/coin selection,
the channel is bistochastic on all 132 states. The usual finite
Hilbert-Schmidt sum-of-squared-commutators identity then identifies its
Heisenberg fixed space with this commutant. These are discrete constructive
certificates, not numerical rank estimates or new Lean proofs.

In contrast, the classical 132-state ordered-cover transition graph is
connected. It has no nonconstant classical configuration invariant. The
protected charge uses coherence of the ordering bit; it is not a second
classical connected component or an inert extra label.

## What the memory law must supply

The bit must remain coherent and every additional service must commute with
Q at the Kraus level. For a trace-preserving service, this condition is also
necessary for exact conservation of this sharp two-valued charge: preservation
of each eigenspace forces all off-diagonal Kraus blocks to vanish, since their
leakage probabilities are nonnegative and sum to zero. A computational-basis
measurement/dephasing destroys the distinction between `|+>` and `|->`.
Resetting the bit to `|0>` also destroys it. Independent Z error probability
`p` per service multiplies charge expectation after `n` services by
`(1-2p)^n`. This explicitly binds protection to a memory implementation.

With source reference `sigma=I66/66`, the sharp states
`rho_+/-=sigma tensor |+/-><+/−|` have exactly the same source marginal.
Every complete classical native record history also has the same law, and
the same conditional source output, in both sectors. Their joint trace
distance is one; after discarding memory it is zero. Thus the source exterior
is reference-like, but the distinguishing charge resides in the observer
bit. The missing particle obligation is a local core realization of that
charge algebra, together with a repair/memory law preserving it and a
compatible moving energy response. No physical species, Lorentz spin,
autonomous spatial defect or mass is identified here.
The equal-tape statement uses committed classical records. Keeping different
histories coherent and erasing their records is a different readout and can
reveal the exchange character; it is not excluded by this statement.

## Execution and reproduction

The [service](service.py) runs the existing proper-code source compiler.
Nine padded source bits and one orientation bit are retained. Record-selected
controlled-X feedback is compiled by conjugating the controlled phase word
with two Y rotations. The integer Walsh audit checks all 1024 addresses for
all 90 seam/current words, including scalar phases, restored parity wires,
padding and collision identity extensions. A separate fresh execution tests
all 132 input columns of C and E on seam `(0,5)`, retaining all 1024 terminal
addresses and every raw output. The combined Frobenius error bounds the
operator error for arbitrary untouched references; its fixed admission
ceiling is `2e-7`, not an outward interval or hardware guarantee.

Native routing, preparation of the coin and record, all feedback pulses,
proper-code load/export calls, and their zero-probability failure complements
are counted per branch in the receipt. The full-column check accepts an
arbitrary supplied source/memory input. A separate executed test prepares
either memory character with one Y pulse and reads it after complete C/E
continuation with one Y pulse and both computational outcomes. Each attempt
adds three C tape bits or four E tape bits; the protected orientation bit is
one persistent bit, not a hidden growing history. Physical initialization,
storage-error and energy costs remain uncalibrated.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_protection --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/defect_protection
```
