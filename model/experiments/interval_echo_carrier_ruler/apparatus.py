"""Conventional material splitters, private phase dwells and independent remote sign heads."""

from functools import cache
from math import pi

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_29_maxwell.seam_symbol import (
    measure as symbol,
)
from model.experiments.interval_echo_consensus.seam_word_meter.chart import (
    exchange_word,
)

EDGES = tuple(symbol.carrier.CELL_SEAMS)
ROD = (12, 13)


def stamp(x):
    return format(float(x), ".12g")


class Ruler:
    def __init__(self):
        self.board = symbol.carrier.federation(0)
        self.net = nx.Graph(self.board.seams)
        self.points, self.axes = symbol.chart()
        self.h = symbol.local_block(1)
        self.exchange = expm(-1j * pi * self.h)
        np.testing.assert_allclose(self.exchange, [[0, 1], [1, 0]], atol=3e-13, rtol=0)
        exchange, contrast = symbol.native_law.local_term_parts()
        vacuum = symbol.native_law.SRC.index((0, 0))
        np.testing.assert_allclose(
            (exchange + symbol.KAPPA * contrast)[:, vacuum], 0, atol=1e-14, rtol=0
        )
        self.phase_port = next(
            p for p in range(12 * 17, 12 * 18) if self.board.degree[p] == 5
        )
        self.heads = []
        for owner in (18, 19):
            pair = next(
                (a, b)
                for a, b in self.board.seams
                if a // 12 == b // 12 == owner
                and self.board.degree[a] == self.board.degree[b] == 5
            )
            refs = [
                next(v for v in self.board.neighbours[u] if v not in pair) for u in pair
            ]
            self.heads.append(dict(owner=owner, pair=pair, references=refs))
        self.source_path = list(symbol.path())
        symbol.pointer()
        self.calendars = {}

    @cache
    def transposition(self, source, dest):
        return exchange_word(nx.shortest_path(self.net, source, dest))

    @staticmethod
    def swap(psi, moves):
        for a, b in moves:
            psi[a], psi[b] = psi[b], psi[a]

    def prepare(self, label, support, phases):
        """Synthesize a fixed equal-amplitude exposure word, before any receiver read."""
        if len(support) == 2:
            seq = list(support)
        else:
            seq = self.source_path
        amp = 1 / np.sqrt(len(support))
        psi = np.zeros(self.board.n, complex)
        psi[seq[0]] = 1
        splits = []
        for j, (a, b) in enumerate(zip(seq, seq[1:])):
            if not self.net.has_edge(a, b):
                raise ValueError("Nonlocal splitter")
            duration = 2 * np.arccos(1 / np.sqrt(len(seq) - j))
            psi[[a, b]] = expm(-1j * duration * self.h) @ psi[[a, b]]
            splits.append([a, b, stamp(duration)])
        dwells = []
        for u, phase in zip(support, phases, strict=True):
            delta = (np.angle(psi[u]) - phase) % (2 * pi)
            if min(delta, 2 * pi - delta) < 1e-12:
                delta = 0
            dwell = delta / (4 * symbol.KAPPA)
            route = self.transposition(u, self.phase_port)
            self.swap(psi, route)
            psi[self.phase_port] *= np.exp(-1j * 4 * symbol.KAPPA * dwell)
            self.swap(psi, reversed(route))
            dwells.append([u, stamp(dwell)])
        np.testing.assert_allclose(
            psi[list(support)], amp * np.exp(1j * np.array(phases)), rtol=0, atol=4e-13
        )
        self.calendars[label] = dict(split=splits, private_dwells=dwells)
        return psi

    def read(self, psi, edge, at, receiver):
        local = psi.copy()
        head = self.heads[at]
        a, b = head["pair"]
        for source, dest in zip(edge, (a, b), strict=True):
            self.swap(local, self.transposition(source, dest))
        analyser = symbol.gates(receiver)[1][:2, :2]
        local[[a, b]] = analyser @ local[[a, b]]
        if any(abs(local[r]) > 3e-13 for r in head["references"]):
            raise ArithmeticError("Remote pointer standard contains a load")
        prob = abs(local) ** 2
        other = np.ones(len(prob), bool)
        other[[a, b]] = False
        return [prob[a], prob[b], prob[other].sum()]

    def layout(self):
        return dict(
            phase=self.phase_port,
            heads=self.heads,
            rod=ROD,
            positions=[[stamp(v) for v in x] for x in self.points],
            private_routes={
                str(u): self.transposition(u, self.phase_port)
                for u in (*range(12), *ROD)
            },
            observer_routes=[
                [
                    dict(
                        source=e,
                        first=self.transposition(e[0], h["pair"][0]),
                        second=self.transposition(e[1], h["pair"][1]),
                    )
                    for e in (ROD, *EDGES)
                ]
                for h in self.heads
            ],
        )
