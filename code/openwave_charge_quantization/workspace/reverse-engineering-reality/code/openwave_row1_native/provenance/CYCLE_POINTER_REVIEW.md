# Independent review of the literal-native winding pointer

The bounded pointer passes independent review. The target is
`../auditor_bridge/native_cycle_meter.py`, with its complete receipt
`NATIVE_CYCLE_METER.json`. The reviewer
[audit_cycle_pointer.py](audit_cycle_pointer.py) reconstructs the one-load
floor/ceiling rule without importing the producer's endpoint, route-law or
success-probability routines. It calls the target executor only to compare
its actual output records.

All 36 combinations of unknown initial pointer port and requested winding
`−1,0,+1` have independently reconstructed exact unconditional output laws.
The separate discovery-timeout, routing-timeout, measurement-timeout and
success probabilities sum to one. Their reset probabilities and worst
combined timeout agree exactly with the producer's rational numbers.
At 16 attempts per desired hop, the worst timeout is
`0.00013732048715552606`; the complete attempt cap is 624. No success weight
is divided by its probability and no timeout is decoded as charge zero.

The implementation was also executed on 576 independently selected coin
tapes, covering all twelve unknown initial ports and all three charges.
Every actual seam, coin and current was replayed against the independent
rule. Each successful output counted the actual signed wraps of the native
pointer on `0→1→8→6→2→10→0`; its final winding agreed with the requested
spatial class. The full fine-history probability was `2^(−attempts)` in
every case. Holds remain in those histories. The pointer's initial hidden
position is not a controller input.

The reset has the required reference scope. Before the first nonzero current,
every fixed coin/current-zero branch is diagonal in the one-load position
basis. Each nonzero-current branch is rank one. After that record, the pointer
location is identified by the record itself; every successful route ends at
port zero. Thus each successful fine reset operator has the form
`c |0><x|`, for its identified initial x. This proves the reset statement on
all source/reference matrix units, while retaining the classical discovery
record and every failure branch. It does not assume an unknown coherent
pointer has a pre-existing classical position.

The earlier simulator-only endpoint contrast was removed from the declared
records before this review. Current records now contain only the native
coin/current, controller stage, selected seam and requested direction.

This review closes the actual pointer/reset component. The spatial class is
identified separately by the integral dual-cell theorem. A single joined
source-history→counter→pointer execution must retain both source and pointer
clocks and their combined unconditional law; it is reviewed separately when
its frozen artifact is available. No platform status is awarded here.
