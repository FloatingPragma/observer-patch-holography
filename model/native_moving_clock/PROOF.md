# Native clocks in a transported source code

Let `W=S(2AA*-I)` be the charged coherent event of the declared two-load
program, and let `V` be the six-word matching code. The original clock has
a complete corrected source channel `Phi` on the 512-address source. Its
reachable matching-support space has dimension 36. The event register and
active memory are retained reference systems for this clock.

The transported encoding is `C_n=W^n A V`. Since both factors of `W` are
unitary, `C_n*C_n=I` and `W C_n=C_(n+1)`. The controlled clock instrument
at frame `n` is the actual original clock instrument conjugated by `W^n`.
Its implementation applies `W^-n`, the native source clock and quadrature,
and `W^n`. Every inverse transport is coherent; it resets neither source,
memory, event, nor a reference register. Each complementary clock outcome
is conjugated with the two admitted outcomes, so their effect sum is the
identity on the full source and event carrier.

For two frames `n1,n2`, interleave the source clock blocks with
`W^(n2-n1)`. Cancelling adjacent exact inverse words proves that the complete
two-time quadrature probabilities equal those of the original clock. This
identity holds on arbitrary source/reference matrix units. Unitary
conjugation preserves the inherited trace-norm error bound for initial
states in the six matching-code words, with arbitrary untouched references.
The bound is not asserted for arbitrary initial states in the 36-row
matching-support space. Complete-channel
contraction bounds additional implementation errors by their sum. With
frames 1 and 4 the word uses 14 forward or inverse coherent events, so the
additional walk contribution is at most `14 epsilon_W` if each implemented
walk or inverse has complete-channel error at most `epsilon_W`. The native
readout, preparation and compensated charge-feedback errors are additional.

The calculation executes all 36 source-support columns on all 3960
source/event coordinates for six character sectors and four frames. The
full padded program also executes its stored reflection on
`512 x 8 x 64` addresses and compares complete native clock records after
forward and inverse transport. Numerical diagnostics and exact structural
unitary identities have separate roles. The reflection compiler supplies
its full-word implementation bound. Controlled swaps, source-clock steps,
compensated feedback, preparation and readout require their additional
complete-channel implementation bounds. No end-to-end hardware precision
certificate is asserted.

The primary program preparation has character 1, spectator 0 and native
law C. Two of the 36 repeated-record cases use that preparation and law.
The other 34 cases are controls over the inherited law and spectator family;
they are not relabeled as the fixed manifest's primary preparation.

The fixed six-word code has population `67/90` after one coherent event
from its first word; `23/90` occupies its complement. A fixed-code clock
cannot be attached to the moving state without an interface proof. The
controlled transported encoding supplies such an interface with an explicit
tracking apparatus. It supplies no spatially local, autonomous clock. Its
preparation and coherent phase are specified services. The work calculation
uses the declared two-hop source observable, with zero reference energy
assigned to memory and event registers; it supplies no physical work scale.
No stationary joint generator, physical particle pole or de Broglie-clock
criterion follows from this controlled experiment.
