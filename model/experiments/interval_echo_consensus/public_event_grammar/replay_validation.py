"""Verify both versions of the source and compare the complete grammar record."""

from ..signed_effects.sequential import compare as compare_records

REPLAY = "model/candidate_laws/interval_echo/replay.py"
OBSERVATIONS = "model/candidate_laws/interval_echo/observations.py"
# Deposited in a2f9e9147. The successor adds optional attendance, service
# histories and scheduling. Default preparation and draw order are preserved.
ORIGINAL = ("ff310f4d0a13d4ef9b060e2f2144b7e4da821741edcb5092baaa76a853daba7e",
            "e8c883004831a466c7ddf0030e1e3e50c18c2b008d2eb3a62b9953a47981e2cb")
# The producer change delegates only the post-acquisition comparison.
REVIEWED = ("fa2a0dfc26506e28cb8e9eeef7d8de6213a631192cbc4813c9138954ecf20293",
            "d4d0d825bef3b643cf13b72f000e79206a86d30ecd5e615c4a3eedbb92292255")

# The publication repair preserves the entire measurement generator and changes
# the verification and file handling around it.
DEPOSITED_REPLAY = "3eaeeb1601a825dd281c545886cb0c1781147c3232e7d1abb5e632dc0c6198e5"
TRANSACTIONAL_REPLAY = "5a38c59ef3d5c8c2293d15e0127f6c0ea60da3b50bdc7ba9085ff5cb7bb6b1c4"
REVIEWED_BINDINGS = {
    (*ORIGINAL, DEPOSITED_REPLAY),
    (*REVIEWED, DEPOSITED_REPLAY),
    (*REVIEWED, TRANSACTIONAL_REPLAY),
}


def compare(saved, actual):
    for record in (saved, actual):
        sources = record["source_audit"]["source_sha256"]
        binding = (record["producer_sha256"], sources[OBSERVATIONS], sources[REPLAY])
        if binding not in REVIEWED_BINDINGS:
            raise ValueError("unreviewed grammar producer, observation source or replay source")
    expected = saved["source_audit"]["source_sha256"]
    current = actual["source_audit"]["source_sha256"]
    replay_audit = {**actual["source_audit"],
                    "source_sha256": {**current, OBSERVATIONS: expected[OBSERVATIONS], REPLAY: expected[REPLAY]}}
    compare_records(saved, {**actual, "source_audit": replay_audit,
                            "producer_sha256": saved["producer_sha256"]})
