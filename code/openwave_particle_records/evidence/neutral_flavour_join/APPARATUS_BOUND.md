# Complete apparatus and stored-source error envelopes

The body laws in RECORDS.json.gz are conditional on ready auxiliary pointer
registers. They are not presented as probabilities for every original
preparation attempt. The full instrument has the following additional
factorization.

Each pointer service starts with a fresh independent uniform mixture on the
full78 total-two configurations. Its current-only reset has exact average
success probability p, stored as a rational in the receipt. A successful
reset prepares the pure ready auxiliary independently of the work/reference.
Every failed reset is a terminal flag retaining its native tape and the work
state. Three successful bit reads, each followed by the conditional native X,
prepare the three-bit blank. The actual source X words prepare the input
family label; the complete source frame word then prepares its tag state.
The receiver uses the inverse frame word and three pointer bits.

The six reset stages therefore have exact first-failure masses
`p^m(1-p)`, m=0,...,5, and all-success mass p^6. For terminal-read failures,
the preceding measured bits and remaining quantum state are retained in the
factorized instrument. The stages after source release additionally carry
that source's stopped success/timeout branch. This describes complete
instruments without enumerating exponentially many fine strings. The archive
retains every fine native tape used in the independent78-input reset check;
its empirical78 successes are not substituted for p.

The ideal native bit instrument is exactly Lüders on arbitrary work/reference
inputs. The saved raw word is compared with its complete eight-column
permutation target. Work-blank preparation is independently checked on all
four one-bit matrix units before tensoring. The full frame words retain all
eight columns; every padding outcome remains. Proper-code load/export,
coherent routing and source control precision remain declared services.

Let t=35/2147483648 be the uniform reset timeout bound. The frozen word
operator envelopes are eta_p=10^-11 for the pointer, eta_X=10^-10 for the
basis/conditional X, and eta_F=2×10^-10 for the frame and inverse frame.
For an approximate word A and exact unitary U with ||A-U||<=eta, the
unhalved completely bounded trace-norm difference is at most
`eta(2+eta)`. Thus each half-diamond term is delta(eta)=eta(2+eta)/2,
and the approximate word's channel norm is at most(1+eta)^2.

The main instrument has six bit services, at most six X words and two
frame words. A conservative complete comparison bound is

    G [6(t+delta_p)+6delta_X+2delta_F],
    G=(1+eta_p)^12 (1+eta_X)^12 (1+eta_F)^4.

An arm with k intermediate full-family reads has6+3k bit services and2+2k
frame words, with the corresponding larger G. The post-run arithmetic audit
checks k=0,...,6. The source's internal cap failures are already part of the
compared instrument; their probability is not hidden in this apparatus error.

The audit treats every stored floating-point entry as an exact binary
rational. It rechecks the charged Y polynomial's three root intervals using
an outward rational sqrt(5) interval. For each proposed eigenvector v and
root approximation lambda, the full residual bounds its component outside
the correct simple eigenline by residual/separation. Its norm defect and
that projection bound give a phase-aligned eigenvector distance. Summing
those distances and exact entrywise bounds on the saved source words proves
the admitted frame envelopes against the exact source eigenframe. Separate
exact entrywise checks verify the X and pointer envelopes. No numerical
polar projection is used to change the saved output matrices.

These bounds concern the stored finite-source realization and its stated
word envelopes. They are not authenticated physical hardware error rates
or a physical pulse-duration/energy calibration. The repeated-read arm in
the frozen receipt is the exact ideal factorized channel; its implementation
comparison has the larger budget computed by the audit. This separates the
reported ideal control from the actually executed single terminal readout
matrices and does not silently rename a diagnostic as a second source run.

The exact successful source-clock mass at k=3 is(255/256)^12. Its conditional
coherent return is one ideally. The original-attempt correct-return weight
also includes p^6 and the apparatus comparison bound. A conditional return
near one must not be reported as near-unit success for every attempt.
