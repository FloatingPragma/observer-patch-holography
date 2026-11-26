"""Twin material readback of versioned, one-publication repairs."""

import argparse
from collections import Counter
from dataclasses import replace
import gzip
import json
from itertools import product, permutations
from pathlib import Path

from .board import Board, allocation
from .contract import (
    EVENTS,
    acceptance_supports,
    assemble,
    batch,
    check_genesis,
    launch,
    mismatch,
    normal,
    parents,
    proposals,
    require,
    validate,
)

READER = 13
SOURCE = launch(1, 0, (0, 0, 0), (1, 0, 1))
ORDERS = [(0, 1), (1, 0)]


class Service:
    def __init__(self, states, hidden=0):
        self.board = Board(states, hidden)
        self.events = []
        self.attempts = []
        check_genesis(self.board.survey())

    def commit(self, tx, challenge=None):
        side = tx.side
        a = allocation(side)
        u = self.board.snapshot(side, READER, "commit-start")
        outcome = validate(tx, u)
        if outcome != "ok":
            self.attempts.append([tx.identifier, outcome])
            return outcome
        target = tx.result()
        if not tx.is_seal:
            for j, bit in enumerate(target[2:5]):
                self.board.write(a["bank1"][j], bit, READER)
                self.board.audit(side)
                if j == 0 and challenge == "interrupted":
                    self.attempts.append([tx.identifier, "interrupted"])
                    return "interrupted"
                if j == 0 and challenge == "early_publication":
                    self.board.write(a["p"], 1, READER)
                    visible = self.board.audit(side)
                    require(
                        visible != u and visible != target,
                        "partial-word control missing",
                    )
                    self.attempts.append([tx.identifier, "mixed word exposed"])
                    return "mixed word exposed"
            if challenge == "stale":
                self.board.write(a["bank0"][0], 1 - u[2], READER)
            if challenge == "boundary_version":
                self.board.write(a["b"], 1 - u[0], READER)
                self.board.write(a["b"], u[0], READER)
                self.board.write(a["g"], 1, READER)
            if challenge == "bad_shadow":
                self.board.write(a["bank1"][2], 1 - target[4], READER)
        # The following check and single publication are one isolated commit
        # appointment. Interrupted staging has never selected an inactive bank.
        final = self.board.snapshot(side, READER, "commit-check")
        outcome = validate(tx, final)
        if outcome == "ok" and not tx.is_seal:
            if tuple(self.board.read(p, READER) for p in a["bank1"]) != target[2:5]:
                outcome = "payload"
        if outcome != "ok":
            self.attempts.append([tx.identifier, outcome])
            return outcome
        self.board.write(a["s" if tx.is_seal else "p"], 1, READER)
        output = self.board.audit(side)
        require(output == target, "published word or version differs from the payload")
        edge = tuple(
            after - before for before, after in zip(final[5:], output[5:], strict=True)
        )
        observed_id = {(1, 0): 1 + side, (0, 1): 3 + side}.get(edge)
        require(
            observed_id == tx.identifier,
            "material publication has the wrong event identity",
        )
        self.events.append(
            [
                observed_id,
                side,
                list(tx.members),
                sorted(tx.writes),
                parents(side, final),
                list(final),
                list(output),
            ]
        )
        self.attempts.append([tx.identifier, "published"])
        return "published"

    def record(self):
        return dict(
            reads=self.board.records,
            raw=self.board.raw,
            events=self.events,
            attempts=self.attempts,
        )


def authenticate(events):
    commits = {}
    for event, side, _, writes, cited, before, after in events:
        require(event not in commits, "reused event")
        if event > 2:
            require(cited == [side + 1] and side + 1 in commits, "missing parent")
            parent_side, parent_writes, parent_after = commits[side + 1]
            # Middle register y belongs to the seal's full local read set.
            # Check the actual parent's write, the child's visible writer and
            # the child's payload against that parent's post-commit middle word.
            require(
                parent_side == side
                and "y" in parent_writes
                and before[5] == 1
                and parents(side, before) == cited
                and before[3] == parent_after[3],
                "unauthenticated value-version",
            )
        else:
            require(not cited, "unexpected genesis writer")
        commits[event] = (side, writes, after)


def orders_trace(data_order, seal_order, hidden=0):
    service = Service(SOURCE, hidden)
    for order in (data_order, seal_order):
        frozen = {t.side: t for t in batch(service.board.survey())}
        require(len(frozen) == 2, "expected a pair of separately prepared components")
        for side in order:
            require(
                service.commit(frozen[side]) == "published",
                "lost source-batch component",
            )
            other = 1 - side
            if side == order[0]:
                require(
                    service.board.audit(other) == frozen[other].snapshot,
                    "disjoint snapshot lost",
                )
    final = service.board.survey()
    require(
        all(normal(u) for u in final) and not batch(final), "terminal reader mismatch"
    )
    require(
        sorted(e[0] for e in service.events) == list(EVENTS), "event carrier incomplete"
    )
    authenticate(service.events)
    return service.record()


def controls():
    result = {}
    for name in (
        "interrupted",
        "early_publication",
        "stale",
        "boundary_version",
        "bad_shadow",
        "split",
        "disagreement",
        "replay",
    ):
        service = Service(SOURCE)
        tx = batch(service.board.survey())[0]
        if name == "split":
            outcome = service.commit(replace(tx, members=tx.members[:-1]))
            require(
                outcome == "component" and not service.events, "split commit admitted"
            )
        elif name == "disagreement":
            members = proposals(0, tx.snapshot)
            p = members[-1]
            key, original = p.payload[0]
            changed = replace(p, payload=((key, 1 - original), p.payload[1]))
            try:
                assemble([*members[:-1], changed], tx.snapshot)
            except ValueError as ex:
                require(str(ex) == "incompatible overlap", "different rejection")
                service.attempts.append([tx.identifier, str(ex)])
            else:
                raise ValueError("incompatible gluing admitted")
        elif name == "replay":
            require(
                service.commit(tx) == "published" and service.commit(tx) == "stale",
                "duplicate commit",
            )
        else:
            outcome = service.commit(tx, name)
            expected = {
                "early_publication": "mixed word exposed",
                "interrupted": "interrupted",
                "stale": "stale",
                "boundary_version": "stale",
                "bad_shadow": "payload",
            }[name]
            require(outcome == expected and not service.events, "failed control")
            if name == "interrupted":
                require(
                    service.board.audit(0) == tx.snapshot
                    and service.commit(tx) == "published",
                    "interruption changed the public state or prevented retry",
                )
        service.board.audit(0)
        result[name] = service.record()
    for name in ("duplicate_ids", "genesis_collision"):
        service = Service(SOURCE)
        if name == "genesis_collision":
            service.board.write(allocation(0)["p"], 1, READER)
        try:
            check_genesis(
                service.board.survey(),
                (1, 1, 3, 4) if name == "duplicate_ids" else EVENTS,
            )
        except ValueError as ex:
            service.attempts.append([0, str(ex)])
        else:
            raise ValueError("Invalid namespace admitted")
        result[name] = service.record()
    return result


def census():
    counts = Counter()
    terminal = {}
    words = list(product((0, 1), repeat=3))
    supports = acceptance_supports()
    for a, b, x, y in product((0, 1), (0, 1), words, words):
        initial = launch(a, b, x, y)
        if not all(mismatch(u) for u in initial):
            continue
        check_genesis(initial)
        counts["sources"] += 1
        for data, seal in product(ORDERS, repeat=2):
            states = initial.copy()
            ids = []
            for order in (data, seal):
                pending = {t.side: t for t in batch(states)}
                for side in order:
                    t = pending[side]
                    require(
                        validate(t, states[side]) == "ok", "snapshot did not survive"
                    )
                    affected = set().union(
                        *(s for s in map(set, supports.values()) if s & t.writes)
                    )
                    require(affected <= t.reads, "missing semantic dependency")
                    members = proposals(side, states[side])
                    for perm in permutations(members):
                        require(
                            assemble(perm, states[side]) == t, "glue-order dependence"
                        )
                    require(
                        parents(side, states[side])
                        == ([1 + side] if t.is_seal else []),
                        "incorrect causal support",
                    )
                    states[side] = t.result()
                    ids.append(t.identifier)
            key = f"{a}{b}"
            require(
                terminal.setdefault(key, states) == states
                and all(normal(u) for u in states),
                "distinct same-boundary normal forms",
            )
            require(sorted(ids) == list(EVENTS), "duplicate or absent label")
            counts["orders"] += 1
    return dict(counts=dict(counts), terminal=terminal, supports=supports)


def calibrate():
    board = Board(SOURCE)
    for side in range(2):
        layout = allocation(side)
        pairs = (
            [layout[key] for key in ("b", "g", "p", "s")]
            + layout["bank0"]
            + layout["bank1"]
        )
        for pair in pairs:
            for bit in (0, 1):
                board.initialize(pair, bit)
                before = board.loads.copy()
                for at in (board.local(side), READER, 15):
                    require(board.read(pair, at) == bit, "incorrect calibration")
                require(before == board.loads, "routed read disturbed a held register")
                board.write(pair, 1 - bit, READER)
                require(board.read(pair, 15) == 1 - bit, "unseen write")
                expected = before.copy()
                i, j = pair
                expected[i], expected[j] = expected[j], expected[i]
                require(
                    expected == board.loads,
                    "one-register write changed another address",
                )
    return board.raw


def run():
    traces = {f"{a[0]}{b[0]}": orders_trace(a, b) for a, b in product(ORDERS, repeat=2)}
    return dict(
        schema=1,
        encoding="Each raw row is [reader, first rail address, 2*left_click+right_click]. The pointer pair 10 denotes zero; 01 denotes one.",
        calibration=calibrate(),
        allocation=[allocation(i) for i in range(2)],
        source=SOURCE,
        orders=traces,
        hidden_control=orders_trace(ORDERS[0], ORDERS[0], hidden=1),
        controls=controls(),
        census=census(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args()
    result = run()
    raw = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    path = Path(__file__).with_name("records.json.gz")
    if args.acquire:
        with path.open("xb") as f:
            f.write(gzip.compress(raw, mtime=0))
    else:
        require(
            gzip.decompress(path.read_bytes()) == raw, "transaction receipt differs"
        )
    print(
        "Verified all material publications and refusals, 196 sources and 784 source-batch schedules."
    )


if __name__ == "__main__":
    main()
