# Prospective seed and replay correction

The first frozen execution and its common consumer are preserved byte-for-byte
under `attempts/first_seeded_run_v1`. All its numerical gates passed. Rotation's
independent review identified two real presentation/interface defects before
the packet was sealed:

1. The harness passed `72361000+K+12` as the controller's RNG seed. Actual
   decisions never decoded that seed and the all-input proof assumes fair
   native choices, but the seed nevertheless encoded K. This left an avoidable
   information channel inconsistent with the strongest unknown-total claim.
2. The independent checker rewrote saved handoffs and results on every replay.
   That could silently repair changed evidence instead of rejecting it.

V2 supplies exactly the same public seed, 72361000, to every controller and
control arm. The seed is independent of K. The checker now verifies saved
outputs by default and records them only through explicit first-write mode.
The source law, inputs, cap, mathematical proof, controls and numeric gates
are unchanged. The new seed may change individual native histories, so the
full run is frozen and executed anew; no old outcomes are discarded.

The known-vacuum datums and bounded neutral preparation remain explicit
assumptions. Choosing to retain one of each selected sign is a controller
choice, not credited as source-generated particle number. The latest master
CAUTIONS and SUBMISSION_CHECKLIST are assessed separately from the original
frozen copies; none of those older bytes are rewritten.
