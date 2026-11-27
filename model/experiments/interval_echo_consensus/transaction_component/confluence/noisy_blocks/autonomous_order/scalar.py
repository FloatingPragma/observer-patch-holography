"""One independently routed scalar trajectory, using native projected noise pulses."""

from ..scalar import scalar_quarter
from ...apparatus import COMMIT, FAR, MASKS, Experiment, layout
from ..instrument import BOUNDARY, NOMINAL
from .dispatcher import SHELVES


def acquire():
    dev = Experiment(BOUNDARY, ((0, 0, 0), (0, 1, 1)))
    dev.initialize(dev.m.rows[14][0], 1)
    for i, pair in enumerate(SHELVES):
        dev.initialize(pair, int(i < 8))

    def shelf():
        rows = [[dev.read_bit(pair, at) for pair in SHELVES] for at in (COMMIT, FAR)]
        assert rows[0] == rows[1]
        return rows[0]

    def mismatch(at):
        out = []
        for s in (0, 1):
            b, word, _ = dev.state(s, at, "dispatch")
            out.extend(x ^ b for x in word)
        return out

    stages = []
    for tick in range(16):
        before = shelf()
        seen = mismatch(COMMIT)
        assert seen == mismatch(FAR)
        labels = [i for i, v in enumerate(before[:8]) if v]
        selected = (max if sum(seen) % 2 else min)(labels)
        dev.write(SHELVES[selected], 0)
        retired = shelf()
        assert [i for i in range(9) if before[i] != retired[i]] == [selected]
        held = tick % 3
        first = scalar_quarter(dev, held != 0)
        mode = 1 + scalar_quarter(dev, held == 2) if first else 0
        side, m = NOMINAL[selected]
        b, word, rank = dev.state(side, COMMIT, "prepare")
        assert dev.guard(side, rank, MASKS[m], FAR) == [word[i] ^ b for i in MASKS[m]]
        target = list(word)
        for i in MASKS[m]:
            target[i] = (b, word[i], 1 - b)[mode]
        new = sum(x ^ b for x in target)
        packet = layout(side)
        for pair, v in zip(packet["banks"][new], target):
            dev.write(pair, v)
        if new != rank:
            for a, c in dev.m.move(packet["selector"][rank], packet["selector"][new]):
                dev.loads[a], dev.loads[c] = dev.loads[c], dev.loads[a]
        assert dev.state(side, FAR, "post") == (b, tuple(target), new)
        if not any(retired[:8]):
            for edge in SHELVES[:8]:
                dev.write(edge, 1)
            dev.write(SHELVES[-1], 1 - retired[-1])
        stages.append(
            dict(
                mask=before,
                flags=seen,
                label=selected,
                retired=retired,
                mode=mode,
                final=shelf(),
                flags_after=mismatch(FAR),
            )
        )
    return dict(stages=stages, raw=dev.raw, observations=dev.observations)
