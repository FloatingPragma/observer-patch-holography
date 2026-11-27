"""An addressed finite dyadic workspace with a two-observer native digit mean."""

from functools import cache
from itertools import chain
from math import pi

import networkx as nx
import numpy as np

from .. import apparatus
from ..chart import exchange_word

LAW = apparatus.LAW
WIDTH = 7
SIZE = 12
CENTER = 64
WRITER = (180, 181)
HEADS = (216, 228)
REFS = ((217, 221), (229, 233))
READ_CODES = (1, 0, 2)  # First bright flag (zero standard) is the high bit.


@cache
def digit_instruments():
    """Integer supports, probabilities and reversal for sharp native writes."""
    pattern = LAW._write_pattern()
    if not np.array_equal(pattern.T @ pattern, 2 * np.eye(9, dtype=int)):
        raise ArithmeticError("The written pair cannot be restored")
    bright = np.array([pair == (1, 1) for pair in LAW.SRC] * LAW.N_LAB)
    read = {}
    for v in range(3):
        for ref in (0, 2):
            col = pattern[:, LAW.SRC.index((v, ref))]
            weight = int(col[bright] @ col[bright])
            if weight not in (0, 2):
                raise ArithmeticError("Reading a sharp load has unresolved outcomes")
            read[v, ref] = weight // 2
    write = {}
    for a in (0, 2):
        for b in (0, 2):
            col = pattern[:, LAW.SRC.index((a, b))]
            support = [divmod(int(i), LAW.N_SRC) for i in np.flatnonzero(col)]
            ends = {LAW.SRC[j] for _, j in support}
            if len(ends) != 1:
                raise ArithmeticError(
                    "Writing separated the coin and loads incorrectly"
                )
            write[a, b] = ends.pop(), [LAW.LABELS[i] for i, _ in support]
    swap = apparatus.pulse(pi)
    np.testing.assert_allclose(
        swap,
        np.eye(9)[:, [LAW.SRC.index((b, a)) for a, b in LAW.SRC]],
        rtol=0,
        atol=apparatus.TOL,
    )
    return read, write


def template_matching(remaining):
    if not remaining:
        return ()
    first = min(remaining)
    for edge in sorted(apparatus.SEAMS):
        if edge[0] == first and edge[1] in remaining:
            suffix = template_matching(remaining - set(edge))
            if suffix is not None:
                return (edge,) + suffix
    return None


def decode_word(raw):
    value = 0
    for k, pair in enumerate(raw):
        if tuple(pair) == (READ_CODES[2], READ_CODES[0]):
            value += 1 << k
        elif tuple(pair) != (READ_CODES[0], READ_CODES[2]):
            return None
    return value


def midpoint(trits):
    """Normalize the independently observed sum digits, then drop the zero low place."""
    carry = 0
    digits = []
    for trit in trits:
        carry, bit = divmod(trit + carry, 2)
        digits.append(bit)
    if digits[0]:
        return None
    return sum(bit << k for k, bit in enumerate(digits[1:] + [carry]))


class Bench:
    def __init__(
        self,
        *,
        board=None,
        words=None,
        flags=None,
        writer=WRITER,
        heads=HEADS,
        references=REFS,
    ):
        self.graph = nx.Graph(
            (apparatus.federation(0) if board is None else board).seams
        )
        matching = template_matching(set(range(12)))
        if words is None:
            dimers = [(12 * i + a, 12 * i + b) for i in range(14) for a, b in matching]
            words = [tuple(dimers[i : i + WIDTH]) for i in range(0, len(dimers), WIDTH)]
        self.words = tuple(tuple(word) for word in words)
        self.flags = (
            tuple((168 + a, 168 + b) for a, b in matching[:2])
            if flags is None
            else tuple(flags)
        )
        self.writer, self.heads, self.references = writer, heads, references
        self.width = len(self.words[0]) if self.words else 0
        pairs = list(chain.from_iterable(self.words)) + list(self.flags)
        reserved = [
            *chain.from_iterable(pairs),
            *writer,
            *heads,
            *chain.from_iterable(references),
        ]
        seams = [
            *pairs,
            writer,
            *((h, r) for h, rr in zip(heads, references, strict=True) for r in rr),
        ]
        if (
            self.width < 1
            or any(len(w) != self.width for w in self.words)
            or len(self.flags) != 2
            or len(heads) != 2
            or any(len(rr) != 2 for rr in references)
            or len(set(reserved)) != len(reserved)
            or not all(self.graph.has_edge(*e) for e in seams)
        ):
            raise ValueError(
                "The supplied material shelves and service ports need disjoint supported seams"
            )
        self.read_table, self.write_table = digit_instruments()

    @cache
    def route(self, source, end):
        excluded = (
            set(chain.from_iterable(self.references))
            | set(self.writer)
            | set(self.heads)
        ) - {source, end}
        graph = self.graph.subgraph(set(self.graph) - excluded)
        route = exchange_word(nx.shortest_path(graph, source, end))
        if any(set(e) == set(self.writer) for e in route):
            raise ArithmeticError("Transit crosses the active write memory")
        return route

    @staticmethod
    def swaps(loads, steps):
        for a, b in steps:
            loads[a], loads[b] = loads[b], loads[a]

    def start(self, values):
        if len(values) != len(self.words) or any(
            not 0 <= x < 1 << self.width for x in values
        ):
            raise ValueError("Bad seven-place enrollment")
        self.loads = [1] * self.graph.number_of_nodes()
        for r0, r2 in self.references:
            self.loads[r0], self.loads[r2] = 0, 2
        for pairs, value in [*zip(self.words, values), (self.flags, 0)]:
            for k, (a, b) in enumerate(pairs):
                self.loads[a] = 2 * ((value >> k) & 1)
                self.loads[b] = 2 - self.loads[a]
        self.tape = []
        self.sensor = "native"

    def observe(self, source, station):
        head = self.heads[station]
        outward = self.route(source, head)
        self.swaps(self.loads, outward)
        code = 0
        for reference in self.references[station]:
            a, b = self.loads[head], self.loads[reference]
            click = (
                self.read_table[a, b]
                if self.sensor == "native"
                else int((a, b) == (1, 1))
            )
            code = 2 * code + click
        self.swaps(self.loads, reversed(outward))
        return code

    def survey(self, pairs, station):
        return [[self.observe(a, station), self.observe(b, station)] for a, b in pairs]

    def read(self, reg, station):
        raw = self.survey(self.words[reg], station)
        value = decode_word(raw)
        self.tape.append(dict(op="read", register=reg, station=station, raw=raw))
        return value

    def store(self, pairs, value, station, name):
        raw = self.survey(pairs, station)
        old = decode_word(raw)
        if old is None or not 0 <= value < 1 << len(pairs):
            raise ArithmeticError("Rejected a store")
        for bit, (a, b) in enumerate(pairs):
            if (old ^ value) & (1 << bit):
                self.swaps(self.loads, [(a, b)])
        self.tape.append(
            dict(
                op="store",
                target=name,
                station=station,
                command=value,
                before=raw,
                after=self.survey(pairs, station),
            )
        )

    def put(self, reg, value, station):
        self.store(self.words[reg], value, station, reg)

    def prepare_mean(self, left, right, parked=False):
        if left == right:
            raise ValueError("The two endpoints must have separate supported shelves")
        raw = [[], []]
        for (a, _), (b, _) in zip(self.words[left], self.words[right], strict=True):
            schedule = self.route(a, self.writer[0]) + self.route(b, self.writer[1])
            self.swaps(self.loads, schedule)
            u, v = self.writer
            incoming = self.loads[u], self.loads[v]
            if not parked:
                outgoing, _ = self.write_table[incoming]
                self.loads[u], self.loads[v] = outgoing
            for station, port in enumerate(self.writer):
                raw[station].append(self.observe(port, station))
            # Both observed loads are sharp on the entire written superposition.
            # The certified inverse echo restores the two incoming loads and blank.
            self.loads[u], self.loads[v] = incoming
            self.swaps(self.loads, reversed(schedule))
        self.tape.append(
            dict(op="mean", registers=[left, right], parked=parked, raw=raw)
        )
        if raw[0] != raw[1] or any(flag not in READ_CODES for flag in raw[0]):
            return "read_disagreement", None
        if self.read(left, 0) is None or self.read(right, 1) is None:
            return "invalid_standard", None
        results = [midpoint([READ_CODES.index(flag) for flag in bank]) for bank in raw]
        if results[0] is None:
            return "off_grid", None
        return "accepted", results

    def mean_axes(self, axes, parked=False):
        axes = tuple(axes)
        registers = tuple(chain.from_iterable(axes))
        if len(set(registers)) != len(registers):
            raise ValueError("Coefficient axes must use disjoint physical shelf pairs")
        pending = []
        for left, right in axes:
            status, results = self.prepare_mean(left, right, parked)
            if status != "accepted":
                return status
            pending.append((left, right, results))
        for left, right, results in pending:
            self.put(left, results[0], 0)
            self.put(right, results[1], 1)
        return "accepted"

    def mean(self, left, right, parked=False):
        return self.mean_axes([(left, right)], parked)

    def public(self, choice):
        self.store(self.flags, {None: 0, "low": 1, "high": 2}[choice], 0, "flags")
        return [self.survey(self.flags, station) for station in (0, 1)]
