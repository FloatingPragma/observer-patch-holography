"""Check the deposited material streams independently of the validator code."""

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey

from .measure import run
from .transport import layout

HERE = Path(__file__).parent


def saved():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def wire_bytes(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":")).encode("ascii")


def check_certificate(c, kind, public, q):
    if not c or c[0] != kind or len(c) != 5:
        return False
    tag, view, val, parent, votes = c
    levels = {"DC": ("LC", "commit"), "LC": ("PC", "accept"), "PC": (None, "prepare")}
    parent_kind, phase = levels[tag]
    if parent_kind:
        if parent[1:3] != [view, val] or not check_certificate(
            parent, parent_kind, public, q
        ):
            return False
        ref = hashlib.sha256(wire_bytes(parent)).hexdigest()
    else:
        if parent is not None:
            return False
        ref = ""
    if len(votes) != q or len({v[1] for v in votes}) != q:
        return False
    for v in votes:
        stage, signer, vview, value, proof, signature = v
        if [stage, vview, value, proof] != [phase, view, val, ref]:
            return False
        try:
            public[signer].verify(bytes.fromhex(signature), wire_bytes(v[:-1]))
        except InvalidSignature:
            return False
    return True


def test_public_signatures_and_independent_reader_outcomes():
    data = saved()
    for name, trial in data["studies"].items():
        q = {"low": 2, "high": 4}.get(name, 3)
        public = [
            Ed25519PublicKey.from_public_bytes(bytes.fromhex(k))
            for k in data["authentication"][name]["keys"]
        ]
        decisions = iter(e for e in trial["observer_events"] if e[0] == "decision")
        for _, receiver, path, packet_id, length in trial["wire"]:
            p = trial["packets"][packet_id]
            raw = wire_bytes(p)
            assert (
                len(raw) * 8 == length and hashlib.sha256(raw).hexdigest() == packet_id
            )
            if p[0] != "DC":
                continue
            table = trial["links"][path]["raw"]
            # Material word 10 is source zero, 01 source one. Each byte
            # below is rebuilt from its eight complete balanced pointer pairs.
            assert all(sorted(pair) == [0, 1] for pair in table)
            octets = bytes(
                sum(table[(v >> j) & 1][1] << j for j in range(8)) for v in range(256)
            )
            try:
                observed = json.loads(raw.translate(octets))
            except (ValueError, UnicodeError):
                observed = None
            accepted = check_certificate(observed, "DC", public, q)
            assert next(decisions)[:3] == ["decision", receiver, accepted]
        assert next(decisions, None) is None
    # Each exercise is a separate enrollment. An adverse history supplies no
    # reusable signature for any of the five other public-key rings.
    keys = [k for item in data["authentication"].values() for k in item["keys"]]
    assert len(keys) == len(set(keys)) == 24


def test_outcomes_lock_tapes_and_held_controls():
    r = saved()
    ss = r["studies"]
    assert ss["low"]["final"][4:] == [0, 1] and ss["equivocation"]["final"][4:] == [
        0,
        0,
    ]
    assert ss["silent"]["final"][1:] == [1] * 5 and ss["high"]["final"] == [None] * 6
    for name in ("orphan", "mixed"):
        t = ss[name]
        assert t["final"][4:] == [1, 1]
        for i, tape in enumerate(t["lock_reads"]):
            far = [row for row in tape if row[:2] == ["read", 5]]
            assert len(far) == 1
            assert far[0][3] == (t["locks"][i][2] if t["locks"][i] else None)
            if t["locks"][i]:
                assert far[0][2][1] == t["locks"][i][1]
    events = ss["orphan"]["observer_events"]
    node = [e for e in events if e[0] in ("accept", "prepare", "commit") and e[1] == 1]
    assert any(
        e[0:3] == ["prepare", 1, 1] and e[-1] and e[3] == e[4] and e[3][:2] == [0, 0]
        for e in node
    )
    assert any(
        e[0:3] == ["commit", 1, 1]
        and e[-1]
        and e[3][:2] == [0, 0]
        and e[4][:2] == [1, 1]
        for e in node
    )
    assert Counter(e[0] for e in node if e[2] == 1 and not e[-1]) == {
        "prepare": 1,
        "accept": 1,
        "commit": 1,
    }
    # Every true transfer excludes the dedicated persistent lock rails.
    scratch = {port for row in layout().rows for port in row[1]}
    for t in ss.values():
        for link in t["links"].values():
            assert all(a not in scratch and b not in scratch for a, b in link["path"])
    assert r["intersections"]["3"] == {
        "cases": 64,
        "minimum_honest_overlap": 1,
        "without_honest_witness": 0,
    }
    assert r["intersections"]["2"]["without_honest_witness"] == 48


def test_public_replay_and_invalid_signature():
    r = saved()
    assert wire_bytes(run(r["authentication"])) == wire_bytes(r)
    # A holder of the receipt can repeat the whole execution without a signing key.
    broken = r["authentication"]["mixed"]["signatures"]
    first = next(iter(broken))
    broken[first] = "00" * 64
    import pytest

    with pytest.raises(ValueError, match="Invalid signature"):
        run(r["authentication"])


def test_collectors_with_extra_duplicate_and_misaddressed_arrivals():
    from .certificates import Signers
    from .measure import Trial

    trial = Trial(Signers())
    all_four = [0, 1, 2, 3]
    changes = trial.reports(0, all_four)
    inputs = [changes[2], None, changes[2], *changes]
    proposal = trial.propose(0, 1, inputs)
    assert trial.rules.proposal(proposal)
    assert [r[0][1] for r in proposal[1]] == [2, 0, 1]
    pc = trial.collect("prepare", proposal, [2, 2, 0, 1, 3], 0)
    assert trial.rules.cert(pc, "PC") and [v[1] for v in pc[4]] == [2, 0, 1]
