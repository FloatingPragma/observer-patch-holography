# Complete finite outcome encoding

No event is dropped from the observer alphabet. The body table records every
native signed-unit input, seam and coin, plus vacuum on every body/reference
seam. A current reconstructs the next full signed load vector. Exact suffix
histories retain every eight-coin word, current and source position under each
of the three tested body laws, for both endpoint starts and every signed,
ordered discovery context. Repeated table evaluation is explicitly not counted
as another native call.

SEARCH has a finite directed acyclic history graph. Its nodes are
`(initial signed port, source-attempt count, selector retry count, partial
five-bit ticket)`. Each reference bit is an independent fair native coin on
the blank seam, with zero current and its full bit recorded. Tickets 0–29
address their named body seam; 30 and 31 are distinct retained rejection
records. At most 16 complete tickets are drawn. A selector timeout is terminal.
After a body event with current zero, the source remains at its original port
and the next source attempt begins; a nonzero current terminates SEARCH with
the ordered departure/arrival record. After 256 unmoved attempts, search timeout
is terminal. The table and these finite transition rules retain a lossless
description of every finite tape; they do not combine distinct record symbols.

For one selector, `q=16^-16` is its failure probability and `a=1-q` its success
probability. Each accepted seam has probability `a/30`. A signed source has five
moving and fifty-five nonmoving body seam/coin alternatives. If
`r=11a/12`, the mass of each of its five directed first moves at attempt `t` is
`r^(t-1) a/60`; selector failure there has mass `r^(t-1)q`. Final search timeout
has mass `r^256`. These explicit retained masses sum to one. All five neighbour
outcomes have equal probability, and both charges have the same discovery law.
No conditioning on a preferred axis is performed.

`RECEIPT.json` holds all attempt-indexed terminal masses as exact rationals.
`RECORDS.json.gz` holds the complete source table and suffix histories, plus
every event of the 96 finite controller trajectories. Those trajectories are
diagnostics; the finite graph above is the proof for all possible search tapes.

The source-field law is compared conditional on the *recorded discovered
frame*, with unconditional discovery, selector failure and timeout masses
reported beside it. This conditioning does not remove a failure of source
apolarization: both source starts and every suffix coin outcome remain in the
field comparison. Matched controls use the identical native discovery prefix,
then change only the holding law. Separate whole-run HOLD trajectories retain
their discovery failures and do not silently condition them away.
