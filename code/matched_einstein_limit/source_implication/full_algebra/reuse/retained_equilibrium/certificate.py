"""Split the certified capped head channel at the registered pre-deadline return bit."""

from itertools import product
from pathlib import Path

import numpy as np

from .. import measure as reuse

HERE = Path(__file__).parent
BANK = HERE / "channels.npz"
GRID = 2**38


def sector(h, ports):
    routes = reuse.ROUTES
    d = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert d == 3840 and row <= 43 and h.dtype == np.int64 and not (h - h.T).nnz
    state = np.zeros((d, len(routes), 4, 2), dtype=np.int64)
    for k, (i, _) in enumerate(routes):
        state[4 * ports[i] : 4 * ports[i] + 4, k, :, 0] = reuse.UNIT * np.eye(
            4, dtype=np.int64
        )
    hits = {route: [] for route in routes}
    identity = np.stack(
        (np.eye(16, dtype=object) * reuse.UNIT**2, np.zeros((16, 16), dtype=object))
    )
    alive = {route: [identity] for route in routes}
    for _ in range(1, reuse.TICKS + 1):
        flat = state.reshape(d, -1)
        for _ in range(reuse.SUBSTEPS):
            flat = reuse.substep(h, flat)
        state = flat.reshape(d, len(routes), 4, 2)
        for k, route in enumerate(routes):
            end = 4 * ports[route[1]]
            at = state[end : end + 4, k]
            hits[route].append(reuse.moment(at[None]))
            at.fill(0)
            alive[route].append(reuse.moment(state[:, k].reshape(-1, 4, 4, 2)))
    result = []
    for route in routes:
        slow = alive[route][12] * reuse.UNIT**2
        fast = np.zeros_like(slow)
        for a, first in enumerate(hits[route], 1):
            reverse = route[::-1]
            slow += reuse.compose(alive[reverse][12 - a], first)
            for b, returned in enumerate(hits[reverse][: 12 - a], 1):
                piece = reuse.compose(returned, first)
                if a + b < 12:
                    fast += piece
                else:
                    slow += piece
        # Round nearest on a dyadic grid, by integer arithmetic.
        panels = np.array([slow, fast], dtype=object)
        panels = (panels * GRID + reuse.UNIT**4 // 2) // reuse.UNIT**4
        result.append(panels.astype(np.int64))
    return np.array(result), float(reuse.channel_error(d, row))


def build():
    form = reuse.station.iso.pairing.formation
    board = form.federation(1)
    ports = form.region(board)
    channel, errors = [], []
    for charge, odd in product((-1, 1), (False, True)):
        h = reuse.native_precision.generator(board, ports, charge, odd)
        block, err = sector(h, ports)
        channel.append(block)
        errors.append(err)
    return np.array(channel).reshape(2, 2, 12, 2, 2, 16, 16), max(errors)


def channels():
    with np.load(BANK) as bank:
        a = bank["numerator"] / GRID
    return a[:, :, :, :, 0] + 1j * a[:, :, :, :, 1]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    bank, error = build()
    if args.write:
        np.savez_compressed(BANK, numerator=bank)
    else:
        with np.load(BANK) as data:
            np.testing.assert_array_equal(bank, data["numerator"])
    print("Complete grouped channel error before grid rounding:", error)
