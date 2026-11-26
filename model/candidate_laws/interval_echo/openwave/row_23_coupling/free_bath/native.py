"""Guarded source contacts, free vacancy propagation and calibrated terminal heads."""

from collections import deque
from functools import cache
import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from model.candidate_laws.interval_echo import candidate

LAW = candidate.native_law


@cache
def contacts(k):
    h0, h1 = LAW.local_term_parts()
    h = h0 + k * h1
    reduced = LAW.seam_term(k)
    np.testing.assert_allclose(h[:9, :9], reduced, atol=2e-13)
    np.testing.assert_allclose(h[9:, :9], 0, atol=2e-13)
    for load in range(3):
        col = np.eye(99)[:, LAW.SRC.index((load, 0))]
        np.testing.assert_allclose(h1 @ col, load**2 * col, atol=2e-13)
    return reduced


@cache
def instruments():
    from ...electron_field.lepton_family.face_frame.shared_report import certify

    reports = certify()
    for load, expected in enumerate(((0, 1), (0, 0), (1, 0))):
        for flip in (True, False):
            if tuple(int(reports[load, r, flip]) for r in (0, 2)) != expected:
                raise ArithmeticError(
                    "Ternary head cannot distinguish the active material"
                )
    timer = []
    for t in (0, 1 / 8, 1 / 2, 1, 3):
        gate = expm(-1j * t * contacts(8))
        col = gate[:, LAW.SRC.index((1, 0))]
        p = abs(col[LAW.SRC.index((0, 1))]) ** 2
        np.testing.assert_allclose(p, np.sin(t / 2) ** 2, atol=2e-13)
        timer.append([t, float(p)])
    return timer


def distances(board, start):
    d = np.full(board.n, -1, dtype=int)
    d[start] = 0
    pending = deque([start])
    while pending:
        v = pending.popleft()
        for u in board.neighbours[v]:
            if d[u] < 0:
                d[u] = d[v] + 1
                pending.append(u)
    if d.min() != 0:
        raise ValueError("Disconnected screen")
    return d


def vacancy_operator(board, k, pins):
    """Project the actual blank generator onto one mobile vacancy and fixed source signs."""
    e = contacts(k).diagonal().reshape(3, 3)
    free = np.array([v for v in range(board.n) if v not in pins])
    tag = np.ones(board.n, dtype=int)
    for v, load in pins.items():
        tag[v] = load
    diag = np.array(
        [
            sum(e[0, tag[u]] - e[1, tag[u]] for u in board.neighbours[v])
            - k * (6 - board.degree[v])
            for v in free
        ]
    )
    # The filled-host energy in this fixed guard setting is a common phase.
    mat = (
        board.adjacency().tocsr()[free][:, free]
        * contacts(k)[LAW.SRC.index((0, 1)), LAW.SRC.index((1, 0))]
    )
    return (mat + sparse.diags(diag)).tocsr(), free


def propagation(h, initial, t):
    if t == 0:
        return initial.copy()
    u = expm_multiply((-1j * t) * h, initial, traceA=(-1j * t) * h.diagonal().sum())
    np.testing.assert_allclose(np.vdot(u, u).real, 1, atol=3e-11)
    return u


def uniform_control(board, blocked):
    free = set(range(board.n)) - {blocked}
    root = min(free)
    order = [root]
    seen = {root: None}
    children = {v: [] for v in free}
    for v in order:
        for u in board.neighbours[v]:
            if u in free and u not in seen:
                order.append(u)
                seen[u] = v
                children[v].append(u)
    if len(order) != len(free):
        raise ValueError("The vacancy preparation tree is disconnected")
    sizes = {}
    for v in reversed(order):
        sizes[v] = 1 + sum(sizes[u] for u in children[v])
    amps = np.zeros(board.n, dtype=complex)
    amps[root] = 1
    for v in order:
        remaining = sizes[v]
        for u in children[v]:
            angle = 2 * np.arcsin(np.sqrt(sizes[u] / remaining))
            phase = np.exp(-1j * angle * (8 + 0.5))
            old = amps[v]
            amps[v] = old * phase * np.cos(angle / 2)
            amps[u] = 1j * old * phase * np.sin(angle / 2)
            remaining -= sizes[u]
    # Private phase pulses, routed by full swaps where needed, align the
    # accumulated phases. Only one vacancy is present, so each selects one component.
    if not any(board.degree[v] == 5 for v in free):
        raise ValueError("No available private phase port")
    corrected = amps * np.exp(-1j * np.angle(amps))
    np.testing.assert_allclose(
        corrected[sorted(free)], 1 / np.sqrt(len(free)), atol=3e-13
    )
    return np.full(len(free), 1 / np.sqrt(len(free)))


@cache
def clock_bins(k):
    """A coherent static four-sign comparison resolves 4K on an energy clock."""
    e = contacts(k).diagonal().reshape(3, 3)
    order = ((2, 2, 1), (2, 0, -1), (0, 2, -1), (0, 0, 1))
    size = 256
    times = 2 * np.pi * np.arange(size) / size
    values = []
    for adjacent in (False, True):
        amp = np.ones(size, dtype=complex)
        for a, b, sign in order:
            # Three arbitrary occupied/vacant/double spectators. Their phases
            # depend on at most one guard and disappear in the four-arm comparison.
            energy = e[a, 0] + e[a, 2] + e[b, 1] + (e[a, b] if adjacent else 0)
            amp *= np.exp(-1j * sign * times * energy)
        read = abs(np.fft.ifft(amp)) ** 2
        peak = int(np.argmax(read))
        value = peak if peak < size // 2 else peak - size
        np.testing.assert_allclose(read, np.eye(1, size, peak).ravel(), atol=1e-11)
        values.append(value)
    # Independently read the single-boundary difference using the two sign
    # settings on an ordinary or vacant, sharp active neighbour.
    boundary = []
    for n in (0, 1, 2):
        pair = np.exp(-2j * times * e[2, n]) * np.exp(2j * times * e[0, n])
        spec = abs(np.fft.ifft(pair)) ** 2
        peak = int(np.argmax(spec))
        np.testing.assert_allclose(spec, np.eye(1, size, peak).ravel(), atol=1e-11)
        boundary.append(peak if peak < size // 2 else peak - size)
    return values, boundary


def classification(board, start, shell):
    """For each possible captured source, score its guard relative to a uniformly addressed probe."""
    aperture = np.flatnonzero(distances(board, start) == shell)
    weight = np.zeros(board.n)
    refused = np.zeros(board.n)
    for q in aperture:
        weight[board.neighbours[q]] += 1 / len(aperture)
        refused[q] = 1 / len(aperture)
    return aperture, weight, refused


def location_instrument(board):
    """Every owner reconstructs its sparse defect word from eleven restored signed writes."""
    std = LAW._write_pattern()
    for i, (a, b) in enumerate(LAW.SRC):
        outcomes = {LAW.LABELS[j // 9][1] for j in np.flatnonzero(std[:, i])}
        if outcomes != {a - b}:
            raise ArithmeticError("Native write does not measure the signed difference")
    trees = set()
    for owner in range(board.cells):
        base = 12 * owner
        queue = [base]
        known = {base}
        tree = []
        for a in queue:
            for b in board.neighbours[a]:
                if b // 12 == owner and b not in known:
                    known.add(b)
                    queue.append(b)
                    tree.append((a - base, b - base))
        if len(tree) != 11:
            raise ArithmeticError("Disconnected local word reader")
        trees.add(tuple(tree))
    for tree in trees:
        for mobile in range(-1, 12):
            for pin in range(-1, 12):
                if pin >= 0 and pin == mobile:
                    continue
                for load in (0, 2):
                    word = np.ones(12, dtype=int)
                    if mobile >= 0:
                        word[mobile] = 0
                    if pin >= 0:
                        word[pin] = load
                    relative = np.zeros(12, dtype=int)
                    for a, b in tree:
                        relative[b] = relative[a] - (word[a] - word[b])
                    decoded = relative - int(np.median(relative)) + 1
                    np.testing.assert_array_equal(decoded, word)
    return dict(
        owners=board.cells,
        distinct_local_trees=[[list(e) for e in t] for t in sorted(trees)],
        rule="Restore each signed write. Set the median relative load to the ordinary load one.",
    )
