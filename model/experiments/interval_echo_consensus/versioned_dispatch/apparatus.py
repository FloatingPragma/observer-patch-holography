"""Sharp bit pairs, independently routed queries and persistent writer stamps."""

from ..boolean_nodes.repair_ballots.measure import Apparatus

NAMES = ("permit", "bank", "left", "right", "onward", "spare")
AUDITORS = (16, 19)
DISPATCH = range(24, 34)
FINAL = range(34, 40)


def decode(raw):
    # Each digit encodes two pointers with the first pointer most significant.
    if any(c not in "12" for c in raw):
        raise ValueError("A physical bit requires exactly one bright rail")
    return tuple(int(c == "1") for c in raw)


def integer(bits):
    return sum(b << n for n, b in enumerate(bits))


def source_record(bits):
    return bits[0], integer(bits[1:])


def dispatch_record(bits):
    a, x, b, y, first, second = bits[:6]
    tag = integer(bits[6:])
    if not any(bits):
        return None, 0
    if a + b != 1 or first + second != 1 or (x and not a) or (y and not b):
        raise ValueError("Incomplete receiver or continuation record")
    return (b, y if b else x, second), tag


def final_record(bits):
    if not any(bits):
        return None, 0
    if not bits[0]:
        raise ValueError("Unpublished final response")
    return bits[1], integer(bits[2:])


def registers(raw):
    word = decode(raw)
    return [source_record(word[4 * i : 4 * i + 4]) for i in range(6)] + [
        dispatch_record(word[24:34]),
        final_record(word[34:40]),
    ]


class Board:
    def __init__(self):
        self.reader = Apparatus(20)
        self.pairs = [self.reader.device.rows[i][j] for i in range(20) for j in (0, 1)]

    def reset(self):
        self.loads = self.reader.device.fresh([0] * 20)
        self.reader.probes = self.reader.exchanges = 0

    def raw(self, slots, at, mode="native"):
        standard = self.reader.device.rows[at][2][0]
        values = []
        for slot in slots:
            a, b = self.pairs[slot]
            left = self.reader.probe(self.loads, a, standard, at, mode)
            right = self.reader.probe(self.loads, b, standard, at, mode)
            values.append(str(2 * left + right))
        return "".join(values)

    def put(self, slot, value, at):
        (old,) = decode(self.raw([slot], at))
        if old != value:
            self.reader.flip(self.loads, self.pairs[slot])

    def stamp(self, first, width, identifier, at):
        if identifier not in range(1, 1 << width):
            raise ValueError("Identifier outside the material stamp field")
        for i in range(width):
            self.put(first + i, (identifier >> i) & 1, at)

    def snapshot(self, at):
        return self.raw(range(40), at)

    def run(self, word, overwrite=None, mode="native"):
        self.reset()
        audit = [self.snapshot(AUDITORS[0])]
        operations = []

        certified = {}

        def finish(identifier, register, versions=()):
            certified[str(identifier)] = list(versions)
            operations.append((identifier, register))
            audit.append(self.snapshot(AUDITORS[0]))

        preps = list(enumerate(word))
        if overwrite is not None:
            preps.append((overwrite, word[overwrite]))
        for i, (reg, bit) in enumerate(preps, 1):
            at = 2 * reg
            self.put(4 * reg, bit, at)
            self.stamp(4 * reg + 1, 3, i, at)
            finish(i, reg)

        sensor = "parked" if mode == "parked" else "native"
        source_raw = self.raw(range(24), 12, sensor)
        accepted = False
        declared = None
        extra = []
        if all(c in "12" for c in source_raw):
            values = decode(source_raw)[::4]
            gate, bank, _, _, turn, spare = values
            declared = [0, 1, 2 + bank, 4] if gate else [0]
            if gate:
                accepted = True
                dest = 12 + bank
                # The addressed receiving bank rereads the selected physical data rail.
                payload = self.raw([4 * (2 + bank)], dest)
                extra.append(payload)
                self.put(24 + 2 * bank + 1, decode(payload)[0], dest)
                self.put(24 + 2 * bank, 1, dest)
                nxt = self.raw([16, 20] if mode == "hidden" else [16], 14)
                extra.append(nxt)
                onward = decode(nxt)[0] ^ (decode(nxt)[1] if mode == "hidden" else 0)
                if onward != (turn ^ (spare if mode == "hidden" else 0)):
                    raise ValueError(
                        "The continuation receiver disagrees with the dispatcher"
                    )
                self.put(28 + onward, 1, 14 + onward)
                self.stamp(30, 4, 8, 15)
                bits = decode(source_raw)
                finish(
                    8,
                    6,
                    [(r, *source_record(bits[4 * r : 4 * r + 4])) for r in declared],
                )

        # Both observers read the actual composite register, including empty outputs.
        public = [self.raw(DISPATCH, at) for at in AUDITORS]
        record, tag = dispatch_record(decode(public[1]))
        if (record is not None) != accepted:
            raise ValueError("Decision and delivered record differ")
        continuation_read = None
        if record is not None:
            station = 17 + record[2]
            continuation_read = self.raw(DISPATCH, station)
            consumed, writer = dispatch_record(decode(continuation_read))
            if consumed != record or tag != writer:
                raise ValueError(
                    "Continuation station and auditor read different versions"
                )
            _, payload, onward = consumed
            self.put(35, payload ^ onward, 17 + onward)
            self.put(34, 1, 17 + onward)
            self.stamp(36, 4, 9, 18)
            finish(9, 7, [(6, consumed, writer)])
        terminal = [self.raw(FINAL, at) for at in AUDITORS]
        peers = self.snapshot(AUDITORS[1])
        if peers != audit[-1] or public[0] != public[1] or terminal[0] != terminal[1]:
            raise ValueError("The two native auditors disagree")
        return {
            "source": list(word),
            "overwrite": overwrite,
            "mode": mode,
            "operations": operations,
            "certified_versions": certified,
            "continuation_read": continuation_read,
            "audit_snapshots": audit,
            "declared_offer_support": declared,
            "dispatcher_raw": source_raw,
            "receiver_data": extra,
            "dispatch_raw": public,
            "continuation_raw": terminal,
            "independent_final_snapshot": peers,
            "probes": self.reader.probes,
            "exchanges": self.reader.exchanges,
        }
