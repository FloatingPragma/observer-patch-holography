"""Independent checks on the routed terminals, signed bytes and local epochs."""

from collections import Counter
import gzip
from itertools import permutations
import json

import networkx as nx
import pytest

from ..certificates import Rules, Signers
from ..transport import encode
from .measure import HERE
from .routes import owner_graph


def archive(name):
    return json.loads(gzip.decompress((HERE / name).read_bytes()))


def test_four_singletons_and_five_cell_witness():
    records = archive("records.json.gz")
    for level, b in records["levels"].items():
        g = owner_graph(int(level))
        # Simplicity, cubic degree and no 3-/4-cycle make every connected
        # four-set a tree, either one star or a unique length-three path.
        assert all(d == 3 for _, d in g.degree)
        assert not any(nx.triangles(g).values())
        assert all(
            len(set(g[a]) & set(g[c])) == 1
            for a in g
            for b0 in g[a]
            for c in g[b0]
            if c != a
        )
        top = b["topology"]
        assert (
            top["connected_quadruples"] == 7 * len(g)
            and top["all_triples_connected_singletons"] == []
        )
        assert len(top["pentagons"]) == 12
        for case in b["cases"].values():
            triple = case["members"]
            cells = set().union(*(case["regions"][i] for i in triple))
            homes = [s[0] for s in case["regions"]]
            for s, r in permutations(triple, 2):
                col = case["columns"]["native"][f"{s}/{r}"]
                assert {j // 12 for pair in col["path"] for j in pair} <= cells
                crossing = {
                    frozenset([a // 12, c // 12])
                    for a, c in col["path"]
                    if a // 12 != c // 12
                }
                assert nx.has_path(
                    nx.Graph([tuple(e) for e in crossing]), homes[s], homes[r]
                )
                assert col["raw"] == [[1, 0], [0, 1]]
                assert (
                    case["columns"]["parked"][f"{s}/{r}"]["raw"][0]
                    == case["columns"]["parked"][f"{s}/{r}"]["raw"][1]
                )


def test_material_clock_and_signed_terminals():
    for board in archive("records.json.gz")["levels"].values():
        for c in board["cases"].values():
            q = Rules(Signers(c["signature_bank"]), 3)
            assert c["authenticated"] == [True] * 3
            assert c["decisions"] == [c["silent"] % 2] * 3 and not c["parked_accepted"]
            assert all(raw == [0, 0] for raw in c["silent_view_mailboxes"].values())
            assert (
                c["timeout_clock"]["tick"]["terminal"]
                == c["timeout_clock"]["parked"]["terminal"]
            )
            assert (
                c["timeout_clock"]["tick"]["period"]
                != c["timeout_clock"]["parked"]["period"]
            )
            for e in c["traffic"]:
                packet = c["packets"][e["packet"]]
                word = encode(packet)
                ones = sum(x.bit_count() for x in word)
                if not e["parked"]:
                    assert e["raw_hist"] == [0, ones, 8 * len(word) - ones, 0]
                if e["phase"] == "decision":
                    assert q.cert(packet, "DC")
                clk = e["clock"]
                assert (
                    len(clk["period"]) * clk["repeats"] + len(clk["suffix"])
                    == e["stop"] - e["start"]
                )
                assert clk["period"] == [[0, 1], [1, 0]]
            assert max(c["phase_slots"].values()) <= c["finite_phase_wire_bound"]


def test_duplicate_fixed_schedule_and_fail_silent_cuts():
    bank = archive("duplication.json.gz")
    counts = Counter()
    for level, cases in bank["levels"].items():
        for f, c in cases.items():
            for key, trial in c["deliveries"].items():
                s, r = map(int, key.split("/"))
                copies = trial["copies"]
                assert len(copies) == 2
                assert all(s in p["via"] and r in p["via"] for p in copies)
                assert sorted(set(copies[0]["via"]) & set(copies[1]["via"])) == sorted(
                    (s, r)
                )
                assert sum(p["slots"] for p in copies) <= c["phase_bound"]
                n = sum(p["admitted_value"] is not None for p in copies)
                counts[r != int(f), n] += 1
                assert (trial["first_value"] is not None) == (r != int(f))
                if r != int(f):
                    assert trial["first_value"] == int(f) % 2
            assert all(
                p["admitted_value"] is None
                for copies in c["both_interiors_disabled"].values()
                for p in copies
            )
    assert counts == {(True, 2): 48, (True, 1): 24, (False, 0): 36}


def test_local_epoch_reconstructed_by_each_receiver():
    from .local_view import views

    bank = archive("local_view.json.gz")
    for cases in bank.values():
        patterns = set()
        for v, c in cases.items():
            if "insufficient_stock" in c:
                assert any(
                    len(r["digits"]) < 10 for r in c["insufficient_stock"].values()
                )
                continue
            regs = c["registers"]
            for i, r in regs.items():
                values = [row[1] for row in r["raw"]]
                assert sum(x << j for j, x in enumerate(values)) == int(v)
                assert len(set(p for d in r["digits"] for p in d)) == 20
                assert all(
                    set(r["cells"]).isdisjoint(other["cells"])
                    for j, other in regs.items()
                    if j != i
                )
            for event in c["events"]:
                assert event["source_raw"] == regs[str(event["source"])]["raw"]
                assert event["receiver_raw"] == regs[str(event["destination"])]["raw"]
                assert event["live_data_raw"] == [[1, 0], [0, 1]]
            rules = Rules(Signers(c["signature_bank"]), 3)
            assert c["decisions"] == [1] * 3 and not any(
                c["wrong_view_decisions"].values()
            )
            assert c["prior_lock_refused_by_single_view_codec"]
            cert = c["packets"][c["decision"]]
            compact = views(cert, int(v), restore=False)
            assert rules.cert(views(compact, int(v), restore=True), "DC")
            assert not rules.cert(views(compact, int(v) + 1, restore=True), "DC")
            with pytest.raises(ValueError):
                views(c["prior_lock_report"], int(v) + 1, restore=False)
            patterns.add(tuple(e["compact_bytes"] for e in c["events"]))
        assert len(patterns) <= 1
    assert all("insufficient_stock" not in c for c in bank["2"].values())
