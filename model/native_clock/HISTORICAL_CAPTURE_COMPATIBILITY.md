# Replaying the retained host capture

The original `host_clock_records.json` remains byte-for-byte unchanged, with
SHA-256 `8a7590e9cd5aa13c7a703651c1851a46c51dfcffd49b125d072914434890900f`.
Its first baseline block took 32,356,500 ns. The timestamps describe that
historical execution. Replaying its native events does not measure today's
host or establish portable execution times.

The captured source manifest differs from the reviewed analysis environment
in the dependency manifest and the added compatibility check in
`host_clock.py`. The original declared requirements and complete capture
source are retained under `historical/` with their capture-time hashes.
The capture did not record Python's version, an installed-package inventory,
or hardware identity; these missing facts are not reconstructed.

The compatibility rule is specific to this exact historical source manifest
and the exact reviewed old/current requirements hashes. All other captured
source-file hashes must match. An AST comparison retains the complete original
module except `analyse` and the one explicitly named compatibility import.
Thus the imports, constants, policy, capture method, native timer loop, full
native event replay, scheduler replay and CLI are unchanged. The updated
analysis checks provenance and reports it separately; the event replay and
calibration arithmetic run on the original records.

The timed loop uses the standard-library seeded generator, bytearray/list
operations, integer arithmetic, the unchanged native repair function and the
recorded monotonic timer. Source/geometry initialization occurs before the
timer starts. Third-party version changes do not justify a claim about timer
performance or future execution. They also do not prevent exact structural
replay of the retained native choices, currents, workload sink and counts.
Every such event is checked; timestamps are read, never regenerated or fitted.

Countercontrols may modify a historical input before analysis. They retain
the exact source declaration and undergo the complete event/metadata checks,
but the receipt explicitly marks them as different from the retained measured
capture. A different historical manifest, a different current requirements
transition, changed native dependency, or changed capture/replay syntax fails.
The current receipt binds this rule and both snapshots in addition to its
current analysis sources. No capture hash or timestamp is rewritten to make
it look current.
