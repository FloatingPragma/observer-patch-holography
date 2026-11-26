import gzip
import hashlib
from itertools import permutations, product
import json

import numpy as np

from ..apparatus import COMMIT, FAR, MASKS, hardware, layout
from .measure import HERE, extract_state


def load(name):
    return json.loads(gzip.decompress((HERE / name).read_bytes()))


def test_snapshot_tapes_and_fault_latch_have_the_declared_physical_addresses():
    data = load("records.json.gz")
    root = HERE.parents[5]
    for source, digest in data["source_sha256"].items():
        assert hashlib.sha256((root / source).read_bytes()).hexdigest() == digest

    for bank in data["banks"]:
        full = (1 << bank["leaves"]) - 1
        for (start, end), snapshot in zip(
            bank["snapshot_ranges"], bank["state"], strict=True
        ):
            stream = iter(bank["reads"][start:end])

            def click(rail, mask):
                station, u, v, active, c = next(stream)
                assert (station, u, v, int(active, 16)) == (
                    FAR,
                    rail,
                    hardware().rows[FAR][2][0],
                    mask,
                )
                read = int(c, 16)
                assert not (read & ~mask)
                return read

            def digit(pair, mask):
                u, v = [click(rail, mask) for rail in pair]
                assert u ^ v == mask
                return v

            for side in (0, 1):
                p = layout(side)
                active = [full ^ click(port, full) for port in p["selector"]]
                assert sum(x.bit_count() for x in active) == bank["leaves"]
                b = digit(p["boundary"], full)
                word = [0] * 3
                for k, select in enumerate(active):
                    if select:
                        for i, pair in enumerate(p["banks"][k]):
                            word[i] |= digit(pair, select)
                assert [b, *word, *active] == [int(x, 16) for x in snapshot[side]]
            assert next(stream, None) is None
        for step, (start, end, *masks) in enumerate(bank["selector_records"]):
            raw = bank["reads"][end - 4 : end]
            for line, rail in zip(
                raw,
                [v for pair in data["selector"]["latch"] for v in pair],
                strict=True,
            ):
                assert line[:3] == [COMMIT, rail, hardware().rows[COMMIT][2][0]]
                assert int(line[3], 16) == full
            a, b, c, d = [int(line[4], 16) for line in raw]
            assert a ^ b == c ^ d == full
            expected = [full ^ b, b ^ d, d]
            if bank["reverse_faults"]:
                expected[1:] = reversed(expected[1:])
            assert expected == [int(x, 16) for x in masks]
            if bank["held"] and step:
                assert end - start == 4


def test_all_permutations_prefixes_and_exact_normal_distance():
    cert = load("records.json.gz")["certificate"]
    masks = [{3 * s + i for i in MASKS[m]} for m in range(4) for s in (0, 1)]
    subsets = {tuple(sorted(mask)) for mask in masks}
    assert len(subsets) == 8
    prefixes = set()
    for order in permutations(range(8)):
        used = 0
        for item in order:
            used |= 1 << item
            prefixes.add(used)
    assert prefixes | {0} == set(range(256))
    corrections = [(4096, 0), (768, 256), (144, 304), (27, 313)]
    largest = 0
    for row in cert["subsets"]:
        enrolled = row["offers"]
        n = [
            sum(i in mask for j, mask in enumerate(masks) if enrolled & (1 << j))
            for i in range(6)
        ]
        assert n == row["hits"]
        for q, y in enumerate(row["expected_4096"]):
            assert y == sum(
                corrections[k][1] + ((q >> i) & 1) * corrections[k][0]
                for i, k in enumerate(n)
            )
            if enrolled == 255:
                assert y == 1878 + 27 * q.bit_count()
            else:
                largest = max(largest, y - 4096 * q.bit_count())
    assert largest == 1860
    for b, q in product(range(4), range(64)):
        # Hamming on the chart of protected bits and six error flags.
        assert (
            min((b ^ normal_b).bit_count() + q.bit_count() for normal_b in range(4))
            == q.bit_count()
        )


def test_independent_scalar_receivers_and_all_complete_selected_branches():
    scalar = load("scalar.json.gz")["runs"]
    banks = load("records.json.gz")["banks"]
    for run in scalar:
        error = run["initial_flags"].copy()
        for step, (dual, parity) in enumerate(run["checkpoints"]):
            for side, (b, word, rank) in enumerate(dual):
                assert b == run["boundary"][side]
                assert [v ^ b for v in word] == error[side * 3 : side * 3 + 3]
                assert rank == sum(v ^ b for v in word)
            assert parity == (run["boundary"][0] ^ run["boundary"][1])
            if step < 8:
                side, m = banks[0]["order"][step]
                chosen = run["branch"][step]
                assert run["selectors"][step][2] == chosen
                for i in MASKS[m]:
                    if chosen != 1:
                        error[side * 3 + i] = int(chosen == 2)
    for source, bank in zip(scalar[-2:], banks[:2], strict=True):
        address = sum(x * 3 ** (7 - i) for i, x in enumerate(source["branch"]))
        for (record, _), snap in zip(source["checkpoints"], bank["state"], strict=True):
            assert [
                [(int(word, 16) >> address) & 1 for word in side[1:4]] for side in snap
            ] == [s[1] for s in record]
    for bank in banks:
        n = bank["leaves"]
        branch = np.array(
            list(product(range(3), repeat=(1 if bank["held"] else 8))), dtype=np.uint8
        )
        s = extract_state(bank["state"][0], n)
        for k, (side, m) in enumerate(bank["order"]):
            mask = sum(1 << (3 * side + j) for j in MASKS[m])
            mode = branch[:, 0 if bank["held"] else k]
            if bank["reverse_faults"]:
                mode = np.where(mode == 0, 0, 3 - mode)
            s = np.where(mode == 0, s & ~mask, np.where(mode == 2, s | mask, s))
            np.testing.assert_array_equal(s, extract_state(bank["state"][k + 1], n))


def test_same_one_offer_marginals_do_not_certify_the_latched_block():
    banks = load("records.json.gz")["banks"]
    fresh = banks[1]["receiver_census"]
    held = banks[4]["receiver_census"]
    assert all(
        a * held["denominator"] == b * fresh["denominator"]
        for a, b in zip(fresh["joint_rows"][1], held["joint_rows"][1], strict=True)
    )
    assert (
        fresh["expected_distances"][-1] == "255/512"
        and held["expected_distances"][-1] == "3/2"
    )
    assert banks[2]["receiver_census"]["expected_distances"][7] == "465/1024"
    assert [banks[i]["receiver_census"]["joint_rows"][-1][0] for i in (0, 2)] == [
        3565164681,
        3130738209,
    ]
    assert banks[5]["receiver_census"]["joint_rows"][-1][0] == 0
