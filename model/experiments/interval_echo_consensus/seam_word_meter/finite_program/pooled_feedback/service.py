"""Fresh physical replicas, shared adaptive words, interval decisions and material flags."""

from dataclasses import dataclass
from functools import cache
from math import log, pi, sqrt

import networkx as nx
import numpy as np

from ... import apparatus as app, chart
from .publication import WIDTH, write_scores


@dataclass(frozen=True)
class Node:
    edge: tuple
    port: int
    station: int
    children: tuple = ()  # Low, high; a leaf reads two independent terminal cohorts.


def program():
    terminal_low = Node((0, 1), 0, 0)
    terminal_high = Node((0, 5), 0, 0)
    return Node((0, 1), 0, 0, (
        Node((0, 5), 0, 1, (terminal_low, terminal_high)),
        Node((0, 7), 0, 1, (terminal_low, terminal_high))))


def requirements(node):
    if (node.edge not in app.SEAMS or node.port not in range(app.N) or node.station not in (0, 1)
            or len(node.children) not in (0, 2)):
        raise ValueError('unsupported_layout_or_read')
    if not node.children:
        return 2
    return 1 + max(requirements(child) for child in node.children)


@cache
def exchange_probability(angle):
    indices = [app.LAW.SRC.index((0, 2)), app.LAW.SRC.index((2, 0))]
    u = app.pulse(angle)
    block = abs(u[np.ix_(indices, indices)])**2
    p = float(block[1, 0])
    np.testing.assert_allclose(block, [[1-p, p], [p, 1-p]], atol=app.TOL, rtol=0)
    for n in (0, 2):
        idx = app.LAW.SRC.index((n, n))
        np.testing.assert_allclose(abs(u[:, idx])**2, np.eye(9)[idx], atol=app.TOL, rtol=0)
    return p


class Cohort:
    """Sharp trajectory sampler of the exact diagonal path effects, with private audit history."""
    def __init__(self, index, n, probabilities, seed):
        self.index = index
        self.rng = np.random.default_rng(seed)
        self.loads = 2 * (self.rng.random((n, app.N)) < probabilities).astype(np.int8)
        self.private_prefix = np.zeros(n, np.int64)
        self.word = []

    def offer(self, edge, p):
        a, b = edge
        self.private_prefix = 3 * self.private_prefix + self.loads[:, a]
        swap = self.rng.random(len(self.loads)) < p
        self.loads[swap, a], self.loads[swap, b] = self.loads[swap, b], self.loads[swap, a]
        self.word.append(edge)

    def consume(self, port, station):
        terminal = self.loads[:, port]
        codes = 3 * self.private_prefix + terminal
        joined = np.bincount(codes, minlength=3 ** (len(self.word) + 1))
        return dict(cohort=self.index, station_owner=(1, 3)[station], port=port, word=list(self.word),
                    private_prefix_then_terminal_counts=joined.tolist(),
                    terminal_counts=np.bincount(terminal, minlength=3).tolist())


def estimate(counts, k, delta, numerical):
    n = sum(counts)
    mean = (counts[2] - counts[0]) / n
    radius = sqrt(2 * log(2*k / delta) / n) + numerical
    interval = [max(-1., mean-radius), min(1., mean+radius)]
    sign = 'low' if interval[1] < 0 else 'high' if interval[0] > 0 else None
    return dict(estimate=mean, radius=radius, interval=interval, sign=sign)


@cache
def flag_layout():
    pairs = ((72, 73), (74, 75))
    board = app.federation(0)
    graph = nx.Graph(board.seams)
    stations = app.layout()['stations']
    paths = []
    if not all(graph.has_edge(*pair) for pair in pairs):
        raise ValueError('missing publication dimers')
    for s in stations:
        forbidden = {p for other in stations for p in other['standards']}
        for port in (p for pair in pairs for p in pair):
            path = nx.shortest_path(graph.subgraph(set(graph) - forbidden), port, s['head'])
            labels = np.arange(board.n)
            for a, b in chart.exchange_word(path):
                labels[a], labels[b] = labels[b], labels[a]
            expected = np.arange(board.n)
            expected[port], expected[s['head']] = expected[s['head']], expected[port]
            np.testing.assert_equal(labels, expected)
            paths.append(path)
    return dict(pairs=pairs, observer_owners=[1, 3], paths=paths)


def publish(sign):
    layout = flag_layout()
    swap = app.pulse(pi)
    projectors = app.projectors()
    reads = []
    for active in (sign == 'low', sign == 'high'):
        v = np.eye(9)[:, app.LAW.SRC.index((0, 2))]
        if active:
            v = swap @ v
        values = []
        for first in (True, False):
            probs = [np.linalg.norm((np.kron(p, np.eye(3)) if first else np.kron(np.eye(3), p)) @ v)**2
                     for p in projectors]
            value = int(np.argmax(probs))
            np.testing.assert_allclose(probs, np.eye(3)[value], atol=app.TOL, rtol=0)
            values.append(value)
        reads.append(values)
    # Repeated visits from both stations use a restored, sharp, nondemolition pointer.
    pointer = ['01', '00', '10']
    return dict(sign=sign, material_word=reads, flags_at_both_owners=[
        [[pointer[n] for n in pair] for pair in reads] for _ in layout['observer_owners']])


def run(node, values, seed, *, n=32768, replicas=4*32768, k=72, delta=.01, numerical=1e-9, angle=pi/2):
    try:
        groups = requirements(node)
    except ValueError as error:
        return dict(refusal=str(error), measured_copies=0, public=publish(None))
    if (n < 1 or 2*n >= 1 << WIDTH or replicas < n*groups or k < groups or not 0 < delta < 1 or numerical < 0
            or groups > 9):
        return dict(refusal='unsupported_budget_or_audit_depth', measured_copies=0, public=publish(None))
    try:
        preparation = chart.prepare_chart(values, -1, 1)
    except ValueError:
        return dict(refusal='outside_chart', measured_copies=0, public=publish(None))
    probabilities = np.array([row['ternary_read'][2] / sum(row['ternary_read']) for row in preparation['coordinates']])
    p = exchange_probability(angle)
    cohorts = [Cohort(i, n, probabilities, child) for i, child in enumerate(seed.spawn(groups))]
    live = list(cohorts)
    reads, actions = [], []
    def observe(port, station):
        row = live.pop(0).consume(port, station)
        # Only the three terminal bin totals cross the controller interface.
        row['decision_read'] = estimate(row['terminal_counts'], k, delta, numerical)
        reads.append(row)
        return row['decision_read']
    final_interval = None
    while True:
        actions.append(dict(edge=node.edge, receiving_cohorts=[g.index for g in live]))
        for group in live:
            group.offer(node.edge, p)
        if not node.children:
            a, b = [observe(node.port, s) for s in (0, 1)]
            lo = max(a['interval'][0], b['interval'][0])
            hi = min(a['interval'][1], b['interval'][1])
            sign = a['sign'] if a['sign'] == b['sign'] and lo <= hi else None
            if sign:
                final_interval = [lo, hi]
            reason = None if sign else 'terminal_unresolved'
            break
        observation = observe(node.port, node.station)
        sign = observation['sign']
        if sign is None:
            reason = 'decision_unresolved'
            break
        node = node.children[int(sign == 'high')]
    public = publish(sign)
    public['score_words_at_both_owners'] = write_scores(sign, reads)
    return dict(refusal=reason, public=public, final_interval=final_interval, source_preparation=preparation,
                native_exchange_probability=p, actions=actions, records=reads, reserved_copies=n*groups,
                measured_copies=n*len(reads), survivor_copies=n*len(live),
                survivor_common_words=[dict(cohort=g.index, word=g.word) for g in live],
                scalar_repair_executions=n*sum(len(row['receiving_cohorts']) for row in actions))
