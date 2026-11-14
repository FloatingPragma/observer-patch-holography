# Opposing review of the balanced native occupancy meter

The frozen `balanced_occupancy.py` source at SHA256
`a23385a0bdea821ec47b41c353ddbcbc45088ffb99ae6b78d259e782105b9732`
passes an independent endpoint and complete-probability review. The review
uses the earlier frozen floor/ceiling implementation, its direct absorbing
probability enumeration and an independent shortest-path calculation.
It does not reuse the producer's native transition helper for replay.

There are756 actual meter replays,24 joined source-to-occupancy replays and
108 exact class/input/retry laws. Both signed cycles and the six-hop
out-and-back zero path have their measured degrees. All discovery, routing
and measurement timeouts remain separate. The six successful hops give
the same factor (1-2^-R)^6 for every class, conditional on fixed initial
pointer position. Combined fine histories retain their powers of1/2 and
the primary source remains idle during its terminal apparatus.

The independence qualification is necessary. At R=1, reset plus six hops
has success probability31/4096 from pointer port0 and91/8192 from port1.
If a source survives with prior probability1/2, and the pointer is prepared
at0 on survival but at1 on movement, the conclusive survival fraction is
62/153. It is not1/2. An independent mixture of those pointer preparations
does preserve1/2. Fresh pointer coins must also be conditionally fair
given the source history; marginal fairness alone is insufficient.

The owning method states the independent-preparation condition correctly.
The review adds this explicit counterexample and reports no implementation
defect. `BALANCED_OCCUPANCY_REVIEW.json` binds the exact source and receipt
bytes. This review leaves all frozen row-one files unchanged.
