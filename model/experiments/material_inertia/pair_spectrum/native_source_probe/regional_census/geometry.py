"""Regional meter sites on the actual 240-port board and the public record alphabet."""

from hashlib import sha256
import json

import numpy as np

from .. import physical
from ..measure import check, regions, routes

READ_WORDS = 17
PREFIXES = 1 + 8 + 8**2
EXHAUSTION = PREFIXES * READ_WORDS
FILTER_MISS = EXHAUSTION + 8**3
BINS = FILTER_MISS + 1


def survey(board):
    profiles = regions(board)
    sites = []
    for owner in sorted(
        {
            k
            for p in profiles
            for k in sorted(set(range(board.cells)) - set(p["owners"]))[:3]
        }
    ):
        field = routes(board, set(range(board.cells)) - {owner})
        paths = dict(
            zip(
                field["target_ports"],
                field["endpoint_transposition_paths"],
                strict=True,
            )
        )
        for path in paths.values():
            physical.validate_route(
                board, path, field["source_park"], field["source_cross"]
            )
        sites.append(
            dict(
                owner=owner,
                park=field["source_park"],
                cross=field["source_cross"],
                contact=field["interaction_port"],
                paths=paths,
            )
        )
    by_owner = {s["owner"]: s for s in sites}
    result = []
    for p in profiles:
        owners = p["owners"]
        candidates = sorted(set(range(board.cells)) - set(owners))[:3]
        heads = [by_owner[d] for d in candidates]
        check(
            len({site[k] for site in heads for k in ("park", "cross")}) == 6,
            "overlapping meter seeds",
        )
        target = p["apparatus"]["target_ports"]
        schedule = [[s["paths"][t] for t in target] for s in heads]
        result.append(
            dict(
                center=p["center"],
                radius=p["radius"],
                owners=owners,
                sites=[
                    {k: s[k] for k in ("park", "cross", "contact", "owner")}
                    for s in heads
                ],
                routed_paths_sha256=sha256(
                    json.dumps(schedule, separators=(",", ":")).encode()
                ).hexdigest(),
            )
        )
    return result


def labels(x, y, profile, decoder):
    region = np.isin(np.arange(240) // 12, profile["owners"]).astype(int)
    q = region[x] - region[y]
    prefix = np.zeros(x.size, dtype=int)
    slot = np.full(x.size, -1)
    for depth, s in enumerate(profile["sites"]):

        def load(port):
            return 1 + (x == port).astype(int) - (y == port).astype(int)

        word = 3 * load(s["park"]) + load(s["cross"])
        pending = slot < 0
        accepted = pending & (word == 4)
        slot[accepted] = (0, 1, 9)[depth] + prefix[accepted]
        rejected = pending & ~accepted
        prefix[rejected] = 8 * prefix[rejected] + word[rejected] - (word[rejected] > 4)
    check(np.all(slot >= 0), "three disjoint unit-host attempts exhausted")
    terminal = np.array([decoder.index(i) for i in (-1, 0, 1)])[q + 1]
    return slot * READ_WORDS + terminal, q


def event_masks(decoder):
    word = np.arange(BINS)
    valid = word < EXHAUSTION
    signed = np.array([decoder[i % READ_WORDS] in (-1, 1) for i in word])
    return valid & signed, valid & (word >= READ_WORDS)
