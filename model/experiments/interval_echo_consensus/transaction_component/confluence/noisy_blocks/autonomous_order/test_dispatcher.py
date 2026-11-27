"""Independent pointer decoding and finite-kernel checks of the held runs."""

from functools import cache
import gzip
import json

from ...apparatus import layout, hardware, COMMIT, FAR
from ..instrument import NOMINAL
from ..measure import flags, one
from . import measure
from .dispatcher import SHELVES


@cache
def archive():
    return json.loads(gzip.decompress((measure.HERE / "records.json.gz").read_bytes()))


class Tapes:
    def __init__(self, rows, cursor, full):
        self.rows, self.cursor, self.full = rows, cursor, full

    def rail(self, site, station, active):
        at, a, b, mask, output = self.rows[self.cursor]
        self.cursor += 1
        assert (at, a, b) == (station, site, hardware().rows[station][2][0])
        assert int(mask, 16) == active
        return int(output, 16)

    def bit(self, edge, station, active):
        first, second = (self.rail(v, station, active) for v in edge)
        assert first ^ second == active
        return second

    def shelf(self):
        both = [
            [self.bit(pair, at, self.full) for pair in SHELVES] for at in (COMMIT, FAR)
        ]
        assert both[0] == both[1]
        return both[0]

    def flag_words(self):
        dual = []
        for at in (COMMIT, FAR):
            cols = []
            for s in range(2):
                p = layout(s)
                ranks = [self.full ^ self.rail(x, at, self.full) for x in p["selector"]]
                assert sum(v.bit_count() for v in ranks) == self.full.bit_count()
                b = self.bit(p["boundary"], at, self.full)
                values = [0, 0, 0]
                for i, enable in enumerate(ranks):
                    if enable:
                        for j, edge in enumerate(p["banks"][i]):
                            values[j] |= self.bit(edge, at, enable)
                cols.extend(v ^ b for v in values)
            dual.append(cols)
        assert dual[0] == dual[1]
        return dual[0]


def test_readouts_are_independently_decodable():
    source = archive()
    for bank in (
        source["tree"],
        source["two_blocks"],
        source["parked"],
        *source["mask_probes"],
    ):
        steps = bank.get("steps") or [bank["proposal"]]
        full = (1 << bank.get("leaves", 256)) - 1
        for s in steps:
            begin, end = s["read_span"]
            tapes = Tapes(bank["reads"], begin, full)
            for columns, record in (
                (tapes.shelf(), "before"),
                (tapes.flag_words(), "flags"),
                (tapes.shelf(), "retired"),
            ):
                assert columns == [int(x, 16) for x in s[record]]
            assert tapes.cursor == end


def test_block_mass_for_every_initial_word_and_held_alternative():
    @cache
    def numerator(q, mask, held):
        if not mask:
            return q.bit_count()
        j = measure.expected_label(q, mask)
        next_mask = mask ^ (1 << j)
        cover = flags(*NOMINAL[j])
        if held is not None:
            return numerator(one(q, cover, held), next_mask, held)
        return sum(
            w * numerator(one(q, cover, m), next_mask, None)
            for m, w in enumerate(measure.WEIGHTS)
        )

    for q in range(64):
        assert numerator(q, 255, None) * 4096 == (27 * q.bit_count() + 1878) * 16**8
    assert sum(w * numerator(63, 255, m) for m, w in enumerate(measure.WEIGHTS)) == 24


def test_complete_joint_law_and_actual_adaptive_fork():
    source = archive()
    result = measure.validate(source)
    assert result == source["summary"] and result["distinct_orders"] == 20
    tree = source["tree"]
    first, second = tree["steps"][:2]
    assert int(first["groups"][7], 16) == (1 << 3**8) - 1
    assert second["groups"][0] == first["fault"][0]
    assert int(second["groups"][6], 16) == int(first["fault"][1], 16) | int(
        first["fault"][2], 16
    )
    scalar = source["scalar"]["stages"]
    for k, s in enumerate(scalar):
        assert s["final"][-1] == ((k + 1) // 8) % 2
    assert scalar[7]["retired"][:8] == scalar[15]["retired"][:8] == [0] * 8


def test_both_native_scans_on_all_enrolled_branches():
    from .scanner import SCRATCH

    records = json.loads(
        gzip.decompress((measure.HERE / "scan_records.json.gz").read_bytes())
    )
    for key, banks in [
        ("tree", [records["tree"]]),
        ("masks", records["masks"]),
        ("parked", [records["parked"]]),
        ("parity_off", [records["parity_off"]]),
    ]:
        for bank in banks:
            full = (1 << bank.get("leaves", 256)) - 1
            for s in bank.get("steps") or [bank["proposal"]]:
                data = s["scanner"]
                begin, end = data["span"]
                tapes = Tapes(bank["reads"], begin, full)

                def both(edge, active):
                    a, b = (tapes.bit(edge, at, active) for at in (COMMIT, FAR))
                    assert a == b
                    return a

                assert [both(p, full) for p in SCRATCH] == [0] * 3
                seen = []
                for at in (COMMIT, FAR):
                    qs = []
                    for side in (0, 1):
                        stock = layout(side)
                        ranks = [
                            full ^ tapes.rail(port, at, full)
                            for port in stock["selector"]
                        ]
                        for i in range(3):
                            out = 0
                            for r, live in enumerate(ranks):
                                if not live:
                                    continue
                                station, a, b, enabled, click = tapes.rows[tapes.cursor]
                                tapes.cursor += 1
                                assert (station, a, b) == (
                                    at,
                                    stock["boundary"][0],
                                    stock["banks"][r][i][0],
                                )
                                assert (
                                    int(enabled, 16) == live
                                    and int(click, 16) & ~live == 0
                                )
                                out |= int(click, 16)
                            qs.append(out)
                    seen.append(qs)
                assert seen[0] == seen[1] == [int(h, 16) for h in s["flags"]]
                odd = 0
                if key != "parity_off":
                    for v in seen[0]:
                        odd ^= v
                exits = [tapes.rail(p, COMMIT, full) for p in SCRATCH[0]]
                assert exits == [full ^ odd, odd] and both(SCRATCH[0], full) == odd
                current = [int(h, 16) for h in s["before"]][:8]
                groups = [0] * 8
                for edge, accept, walk in [
                    (SCRATCH[1], full ^ odd, range(8)),
                    (SCRATCH[2], odd, reversed(range(8))),
                ]:
                    for label in walk:
                        assert both(edge, full) == accept
                        fired = both(SHELVES[label], accept)
                        assert fired == accept & current[label]
                        if key != "parked":
                            current[label] ^= fired
                        accept ^= fired
                        groups[label] |= fired
                    assert both(edge, full) == accept
                assert [both(p, full) for p in SCRATCH] == [0] * 3
                assert current == [int(h, 16) for h in s["retired"][:8]]
                assert groups == [int(h, 16) for h in s["groups"]]
                assert tapes.cursor == end
    assert records["summary"] == archive()["summary"]


def test_native_rail_inversion():
    import numpy as np
    from ..scalar import seam

    pi_gate = abs(seam.pulse(np.pi)) ** 2
    for word in ((2, 0), (0, 2)):
        predicted = np.zeros(9)
        predicted[seam.LAW.SRC.index(word[::-1])] = 1
        np.testing.assert_allclose(
            pi_gate[:, seam.LAW.SRC.index(word)], predicted, atol=seam.TOL
        )
