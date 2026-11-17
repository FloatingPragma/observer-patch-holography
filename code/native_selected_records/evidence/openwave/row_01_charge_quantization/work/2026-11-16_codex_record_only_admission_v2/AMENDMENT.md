# Strict saved-receipt verification

The first registered scientific run passed all gates. Before independent
checking and sealing, source inspection identified that the default producer
compared newly calculated scientific JSON with saved evidence but did not also
require the saved receipt's Boolean verdict and evidence digest. A false receipt
or stale digest therefore needed an explicit rejection gate of its own.

The amendment adds those unconditional checks to the default producer and
independent checker. Explicit --evidence on the checker remains a labelled
scientific mutation-test interface and does not claim the original custody.
Default paths are read-only. Only --record-to NEW_DIRECTORY creates outputs,
exclusively; scientific failure remains nonzero. No law, controller, decoder,
input, record, outcome, finite gate or theorem changed. The original source and
registry are preserved in amendments/receipt_custody. The first record remains
unchanged. This verification-only amendment is frozen before its replay and
mutation checks.
