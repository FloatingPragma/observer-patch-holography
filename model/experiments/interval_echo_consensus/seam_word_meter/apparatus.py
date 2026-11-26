"""Derive the working instruments and their readback effects from the native seam."""

from functools import lru_cache
from math import asin, sqrt, pi

import networkx as nx
import numpy as np
import sympy as sp
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import CELL_SEAMS, federation

LAW = candidate.native_law
SEAMS = tuple(CELL_SEAMS)
N = 12
TOL = 4e-12


@lru_cache(None)
def pulse(t):
    ex, con = LAW.local_term_parts()
    full = expm(-1j * t * (ex + candidate.KAPPA * con))
    np.testing.assert_allclose(full[9:, :9], 0, atol=TOL)
    return full[:9, :9]


@lru_cache(None)
def projectors():
    """Two extreme references implement three Lüders projectors, including their phases."""
    u = LAW.echo()
    bright = np.tile([s == (1, 1) for s in LAW.SRC], 11)
    flags = {}
    for ref in (0, 2):
        cols = [LAW.SRC.index((v, ref)) for v in range(3)]
        for click in (False, True):
            back = u.T.conj() @ (u[:, cols] * (bright == click)[:, None])
            mask = np.ones(LAW.DIM, dtype=bool)
            mask[cols] = False
            np.testing.assert_allclose(back[mask], 0, atol=TOL)
            flags[ref, click] = back[cols]
    outcomes = {}
    for first in (False, True):
        for second in (False, True):
            outcomes[int(first), int(second)] = flags[2, second] @ flags[0, first]
    p = [outcomes[key] for key in [(0, 1), (0, 0), (1, 0)]]
    for i, m in enumerate(p):
        np.testing.assert_allclose(m, np.diag(np.eye(3)[i]), atol=TOL)
    np.testing.assert_allclose(outcomes[1, 1], 0, atol=TOL)
    return p


def layout():
    """Both stations transpose the actual queried port into the sensor and restore it."""
    board = federation(0)
    graph = nx.Graph(board.seams)
    actual = {tuple(sorted((a, b))) for a, b in graph.edges if a < 12 and b < 12}
    if actual != set(SEAMS) or len(actual) != len(SEAMS):
        raise ValueError("The source alphabet differs from the cell seams")
    swap = pulse(pi)
    np.testing.assert_allclose(
        swap, np.eye(9)[:, [LAW.SRC.index((b, a)) for a, b in LAW.SRC]], atol=TOL
    )
    stations = []
    calendar = (24, 25)
    if not graph.has_edge(*calendar):
        raise ValueError("Calendar ports do not share a native seam")
    for owner in (1, 3):
        base = 12 * owner
        ref0, ref2 = (base + j for j in (1, 5))
        if not graph.has_edge(base, ref0) or not graph.has_edge(base, ref2):
            raise ValueError("Missing reference seam")
        routes = []
        for port in [*range(N), calendar[0]]:
            net = graph.subgraph(set(graph) - {ref0, ref2, calendar[1]})
            path = nx.shortest_path(net, port, base)
            pairs = list(zip(path, path[1:]))
            swaps = pairs + pairs[-2::-1]
            labels = list(range(board.n))
            for a, b in swaps:
                labels[a], labels[b] = labels[b], labels[a]
            ideal = list(range(board.n))
            ideal[port], ideal[base] = ideal[base], ideal[port]
            if labels != ideal:
                raise ValueError("Sensing route changes a spectator")
            routes.append(path)
        stations.append(dict(head=base, standards=[ref0, ref2], paths=routes))
    return {"seams": SEAMS, "calendar": calendar, "stations": stations}


def local_effects():
    """Full nine-column Heisenberg test, retaining coherences and every local bin."""
    ps = projectors()
    gate = pulse(pi / 2)
    pinches = [np.kron(p, np.eye(3)) for p in ps]
    rows = []
    for k, p in enumerate(ps):
        left, right = np.kron(p, np.eye(3)), np.kron(np.eye(3), p)
        target = (left + right) / 2
        for side, op in (("first", left), ("second", right)):
            free = gate.T.conj() @ op @ gate
            actual = sum(d @ free @ d for d in pinches)
            np.testing.assert_allclose(actual, target, atol=TOL)
            # Querying the other endpoint preserves the entire scalar alphabet.
            other = [np.kron(np.eye(3), q) for q in ps]
            np.testing.assert_allclose(
                sum(d @ free @ d for d in other), target, atol=TOL
            )
            rows.append([k, side, float(np.linalg.norm(free - target, ord=2))])
    core = np.eye(99)[:, LAW.SRC.index((0, 2))]
    echo = LAW.echo()
    num = np.tile([a for a, b in LAW.SRC], 11)
    after = [
        int(round(abs(w) ** 2 @ num)) for w in (core, echo @ core, echo @ (echo @ core))
    ]
    witness = {}
    for sign in (1, -1):
        state = np.zeros(9, dtype=complex)
        state[LAW.SRC.index((1, 0))] = 1 / np.sqrt(2)
        state[LAW.SRC.index((0, 1))] = sign * 1j / np.sqrt(2)
        n = np.diag([a for a, b in LAW.SRC])
        held = gate.T.conj() @ n @ gate
        witness[str(sign)] = [
            float(np.vdot(state, obs @ state).real)
            for obs in (n, held, sum(d @ held @ d for d in pinches))
        ]
    return {
        "unread_effect_errors": rows,
        "phase_inputs_initial_unread_read": witness,
        "initial_first_port_then_two_echoes": after,
    }


@lru_cache(None)
def flight(edge, t):
    """One-load block in the physical convention, including the untouched empty pair."""
    a, b = SEAMS[edge]
    gate = pulse(t)
    pairs = [LAW.SRC.index(s) for s in ((1, 0), (0, 1))]
    out = np.eye(N, dtype=complex) * gate[0, 0]
    out[np.ix_([a, b], [a, b])] = gate[np.ix_(pairs, pairs)]
    return out


def observed_word(symbols, starts, t=pi / 2, read=True):
    """Unnormalized amplitudes per complete ordered read transcript; preparations are columns."""
    branches = {(): starts}
    for edge in symbols:
        a, _ = SEAMS[edge]
        gate = flight(edge, t)
        if read:
            local_values = np.zeros(N, dtype=int)
            local_values[a] = 1
            split = {}
            for value, p in enumerate(projectors()[:2]):
                diag = p.diagonal()[local_values, None]
                for transcript, amps in branches.items():
                    split[transcript + (value,)] = gate @ (diag * amps)
            branches = split
        else:
            branches = {key: gate @ amps for key, amps in branches.items()}
    return {key: abs(amplitude) ** 2 for key, amplitude in branches.items()}


@lru_cache(None)
def calendar_bit(p):
    """Tuned blank dimer, two-reference load read, and branch-conditioned reset."""
    angle = 2 * asin(sqrt(p))
    start = np.eye(9)[:, LAW.SRC.index((1, 0))]
    out = pulse(angle) @ start
    probs = []
    for outcome in (0, 1):
        # The left endpoint reads one on a stay and zero on a crossing.
        read = np.kron(projectors()[1 - outcome], np.eye(3)) @ out
        weight = float(np.vdot(read, read).real)
        back = pulse(pi) @ read if outcome else read
        np.testing.assert_allclose(np.delete(back, LAW.SRC.index((1, 0))), 0, atol=TOL)
        probs.append(weight)
    np.testing.assert_allclose(probs, [1 - p, p], atol=TOL)
    return probs


def schedule(weights):
    """Return exact emitted masses and independently evaluated numerical dimer diagnostics.

    A deterministic continuation consumes its full word without growing the Python
    call stack. The exact controller carries a positive branch even when a binary64
    propagation cannot numerically distinguish its probability from an endpoint.
    """
    law = {
        w: sp.Rational(v) if isinstance(v, float) and np.isfinite(v) else sp.sympify(v)
        for w, v in weights.items()
    }
    if not law or any(
        p.is_nonnegative is not True or p.is_finite is not True for p in law.values()
    ):
        raise ValueError("Finite, known nonnegative weights are required")
    if sp.simplify(sum(law.values()) - 1) != 0:
        raise ValueError("Word weights must total one")
    if len({len(w) for w in law}) != 1 or any(
        e not in range(len(SEAMS)) for w in law for e in w
    ):
        raise ValueError("Inconsistent word length or symbol")
    choices = {w: p for w, p in law.items() if p != 0}
    exact, numeric, nodes = {}, {}, []
    pending = [(0, choices, sp.Integer(1), 1.0)]
    while pending:
        depth, branch, mass, sampled = pending.pop()
        words = list(branch)
        if len(words) == 1:
            word = words[0]
            exact[word] = mass
            numeric[word] = sampled
            continue
        while all(w[depth] == words[0][depth] for w in words):
            depth += 1
        groups = {}
        for w, p in branch.items():
            groups.setdefault(w[depth], {})[w] = p
        splits = [(sorted(groups), mass, sampled)]
        while splits:
            options, total, approximate = splits.pop()
            if len(options) == 1:
                pending.append((depth + 1, groups[options[0]], total, approximate))
                continue
            mid = len(options) // 2
            halves = options[:mid], options[mid:]
            sums = [sum(sum(groups[e].values()) for e in half) for half in halves]
            p = sp.cancel(sums[1] / sum(sums))
            nominal = float(p)
            observed = None
            if min(nominal, 1 - nominal) > TOL:
                observed = calendar_bit(nominal)
            nodes.append(
                {
                    "emitted_prefix": list(words[0][:depth]),
                    "options": options,
                    "exact_crossing": str(p),
                    "native_binary64_probs": observed,
                }
            )
            for i, half in reversed(list(enumerate(halves))):
                result = (
                    None
                    if approximate is None or observed is None
                    else approximate * observed[i]
                )
                if result == 0:
                    result = None
                splits.append((half, sp.cancel(total * (p if i else 1 - p)), result))
    return exact, numeric, nodes
