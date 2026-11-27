"""Verify both versions of the source and compare the complete grammar record."""

from ..signed_effects.sequential import compare as compare_records

OBSERVATIONS = "model/candidate_laws/interval_echo/observations.py"
# Deposited in a2f9e9147. The successor adds optional attendance, service
# histories and scheduling. Default preparation and draw order are preserved.
ORIGINAL = ("ff310f4d0a13d4ef9b060e2f2144b7e4da821741edcb5092baaa76a853daba7e",
            "e8c883004831a466c7ddf0030e1e3e50c18c2b008d2eb3a62b9953a47981e2cb")
# The producer change delegates only the post-acquisition comparison.
REVIEWED = ("fa2a0dfc26506e28cb8e9eeef7d8de6213a631192cbc4813c9138954ecf20293",
            "d4d0d825bef3b643cf13b72f000e79206a86d30ecd5e615c4a3eedbb92292255")


def compare(saved, actual):
    for record in (saved, actual):
        source = record["source_audit"]["source_sha256"][OBSERVATIONS]
        if (record["producer_sha256"], source) not in (ORIGINAL, REVIEWED):
            raise ValueError("unreviewed grammar producer or observation source")
    expected = saved["source_audit"]["source_sha256"]
    current = actual["source_audit"]["source_sha256"]
    replay_audit = {**actual["source_audit"],
                    "source_sha256": {**current, OBSERVATIONS: expected[OBSERVATIONS]}}
    compare_records(saved, {**actual, "source_audit": replay_audit,
                            "producer_sha256": saved["producer_sha256"]})
