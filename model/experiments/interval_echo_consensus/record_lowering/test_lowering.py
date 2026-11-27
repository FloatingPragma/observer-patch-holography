"""Independent address, intervention, cost and arbitrary finite truth-table checks."""

from collections import Counter
import gzip
import hashlib
from itertools import product
import json

from model.executable_model.carrier import federation
from .compiler import compile_table
from .hardware import Lane, domains
from .measure import HERE, play


def archive():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def test_all_binary_functions_on_material_detector_with_both_prior_targets():
    board = federation(0)
    owners = domains(board)[0]
    for truth in range(16):
        lane = Lane(board, owners)
        lane.gate(((2, 1),), 3)
        lane.table((0, 1), 3, tuple((truth >> k) & 1 for k in range(4)))
        lane.survey("done", (0, 1, 2, 3, 8))
        for w in range(8):
            got = play(board, [lane], [(w, 0)])["branches"][0]["readbacks"]["done"]
            a, b, t = ((w >> k) & 1 for k in range(3))
            assert got == [a, b, t, t ^ ((truth >> (w % 4)) & 1), 0]


def test_compiler_all_three_bit_functions_and_vector_outputs():
    for truth in range(256):
        rows = [((truth >> i) & 1) + 2 * (1 - ((truth >> i) & 1)) for i in range(8)]
        gates = compile_table((0, 1, 2), (3, 4), rows, (5,))
        for x, old in product(range(8), range(4)):
            record = [(x >> j) & 1 for j in range(3)] + [old & 1, old >> 1, 0]
            initial = record.copy()
            for controls, addr in gates:
                record[addr] ^= all(record[i] == bit for i, bit in controls)
            assert record[:3] == initial[:3] and record[5] == 0
            assert record[3] + 2 * record[4] == old ^ rows[x]
            for controls, addr in reversed(gates):
                record[addr] ^= all(record[i] == bit for i, bit in controls)
            assert record == initial


def test_complete_word_factorization_and_resource_accounting():
    r = archive()
    for result in r["records"]:
        board = federation(result["level"])
        physical = {tuple(sorted(edge)) for edge in board.seams}
        owners = [set(lane["owners"]) for lane in result["lanes"]]
        assert owners[0].isdisjoint(owners[1])
        marginals = [{}, {}]
        words = [len(lane["word"]) for lane in result["lanes"]]
        for launch in result["bank"]:
            assert launch["calendar"] == max(words)
            for i, (case, branch, lane) in enumerate(
                zip(launch["inputs"], launch["branches"], result["lanes"], strict=True)
            ):
                key = tuple(case)
                assert marginals[i].setdefault(key, branch) == branch
                raw = iter(map(int, branch["flags"]))
                ports = lane["word"]
                fired = 0
                slots = {}
                for t, (op, a, b) in enumerate(ports):
                    assert a // 12 in owners[i] and b // 12 in owners[i]
                    assert (a == b and op == "wait") or tuple(
                        sorted((a, b))
                    ) in physical
                    if op in ("project", "intervene", "shutter"):
                        slots[t] = next(raw)
                        if op != "project":
                            fired += slots[t]
                assert next(raw, None) is None
                tally = Counter(op for op, _, _ in ports)
                assert (
                    branch["costs"]["native_exchanges"]
                    == tally["swap"] + tally["flip"] + fired
                )
                assert branch["costs"]["echo_or_inverse"] == 284
                assert branch["costs"]["projections"] == 128
                taps = {}
                for label, at, reg, a, b in lane["reads"]:
                    assert slots[a] + slots[b] == 1
                    taps.setdefault((label, at), []).append(slots[a])
                assert all(
                    taps[k, 0] == taps[k, 1] == v
                    for k, v in branch["readbacks"].items()
                )
        assert [len(x) for x in marginals] == [16, 16]
        serial = result["serial"]
        assert serial["calendar"] == sum(words)
        assert all(
            marginals[i][tuple(case)] == b
            for i, (case, b) in enumerate(
                zip(serial["inputs"], serial["branches"], strict=True)
            )
        )
    a, b = r["records"]
    assert [[v["readbacks"] for v in x["branches"]] for x in a["bank"]] == [
        [v["readbacks"] for v in x["branches"]] for x in b["bank"]
    ]
    assert a["bank"][0]["calendar"] != b["bank"][0]["calendar"]


def test_prepublication_intervention_stale_consumer_and_ancilla_controls():
    r = archive()
    for name, sha in r["sources"].items():
        assert hashlib.sha256((HERE.parents[3] / name).read_bytes()).hexdigest() == sha
    for board in r["records"]:
        native = [board["bank"][17 * j]["branches"][0]["readbacks"] for j in range(16)]
        for j, item in enumerate(native):
            word, flip = divmod(j, 2)
            initial = int(word.bit_count() >= 2)
            assert item["cached"] == [initial, initial, 0]
            assert item["edited"] == [initial ^ flip, initial, 0]
            assert item["sealed"] == [initial ^ flip, initial, initial ^ flip, 1]
        changed = {}
        for mode, runs in board["controls"].items():
            changed[mode] = [
                j
                for j, launch in enumerate(runs)
                if launch["branches"][0]["readbacks"]["end"][7] != native[j]["end"][7]
            ]
        assert changed["stale"] == list(range(1, 16, 2))
        assert len(changed["parked"]) == 8
        assert changed["dirty"] == list(range(16))
