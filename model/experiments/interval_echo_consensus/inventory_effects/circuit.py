"""Compile a two-input, six-address instrument on the native ternary triangle."""

from functools import lru_cache
from itertools import combinations, product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation


# Consecutive sharp configurations differ by one swap.
WORDS = ((0, 1, 2), (0, 2, 1), (2, 0, 1), (2, 1, 0), (1, 2, 0), (1, 0, 2))
INDEX = [np.ravel_multi_index(word + (1, 1, 1), (3,) * 6) for word in WORDS]
K = candidate.KAPPA
TAU = 2 * np.pi


@lru_cache(None)
def ports():
    board = federation(0)
    for data in combinations(range(12), 3):
        if any(board.degree[i] != 5 for i in data):
            continue
        if any(j not in board.neighbours[i] for i, j in combinations(data, 2)):
            continue
        options = [sorted(set(map(int, board.neighbours[i])) - set(data)) for i in data]
        for guards in product(*options):
            if len(set(guards)) == 3:
                return data + guards
    raise ValueError("No private triangle with three distinct guard neighbors")


def wrap(angle):
    """Reduce to [0, 2pi) and send roundoff below a multiple of 2pi to zero.

    Compiled phases that are exact multiples of 2pi arrive with a 1e-16
    roundoff whose sign differs between linear algebra libraries. Under the
    registered law both representatives give the same operation. The
    exchange-only control would otherwise run either a zero or a full-period
    exposure depending on the platform.
    """
    r = angle % TAU
    return 0.0 if TAU - r < 1e-12 else float(r)


def argument(z):
    """np.angle on (-pi, pi], with values within roundoff of -pi sent to +pi.

    A negative real entry carries an imaginary roundoff whose sign differs
    between platforms. Keeping one branch keeps the compiled phases canonical.
    """
    a = float(np.angle(z))
    return np.pi if a < 1e-12 - np.pi else a


def private(i, phase):
    return ("private", i, i, wrap(-phase) / K)


def contact(i, j, time):
    return ("contact", i, j, time)


def sign_echo(i, j, time, flip):
    return [contact(i, j, time / 2), *flip, contact(i, j, time / 2), *flip]


def label_phase(i, phases):
    """Diagonal qutrit phases up to a common phase when guard i+3 starts in 1."""
    alpha = (phases[2] - phases[0]) / 4
    beta = phases[1] - phases[0] - alpha
    guard_time = wrap(beta) / (K + 0.5)
    flip = [private(i, np.pi)]
    return [private(i, alpha), *sign_echo(i, i + 3, guard_time, flip)]


def word_rotation(row, theta, lower_phase):
    """Givens matrix with upper off-diagonal -conj(z)sin(theta), lower z sin(theta)."""
    a, b = WORDS[row : row + 2]
    i, j = [u for u in range(3) if a[u] != b[u]]
    remaining = 3 - i - j
    time = 2 * theta
    # The third label determines the scalar part of the active contact in this inventory.
    energies = []
    for label in range(3):
        u, v = sorted(set(range(3)) - {label})
        energies.append(K * (u - v) ** 2 + 0.5)
    toggle = np.zeros(3)
    toggle[a[remaining]] = np.pi
    phi = np.zeros(3)
    phi[b[j]] = lower_phase - np.pi / 2
    return [
        *label_phase(j, -phi),
        *sign_echo(i, j, time, label_phase(j, toggle)),
        *label_phase(remaining, time * np.array(energies)),
        *label_phase(j, phi),
    ]


def target(menu):
    y = np.array([[0, -1j], [1j, 0]])
    h = np.eye(2) / 4 + y / 8
    unit = np.eye(2)
    effects = {
        "fine": (h, h, unit - 2 * h),
        "union": (2 * h, unit - 2 * h, 0 * h),
        "single": (h, unit - h, 0 * h),
    }[menu]
    blocks = []
    for j, effect in enumerate(effects):
        vals, basis = np.linalg.eigh(effect)
        block = (basis * np.sqrt(vals)) @ basis.conj().T
        if menu != "fine":
            angle = (j + 1) * np.pi / (5 if menu == "single" else 7)
            block = (
                np.array(
                    [
                        [np.cos(angle), 1j * np.sin(angle)],
                        [1j * np.sin(angle), np.cos(angle)],
                    ]
                )
                @ block
            )
        blocks.append(block)
    return np.vstack(blocks)


def compile_reader(iso):
    reduced = iso.copy()
    pulses = []
    for col in (0, 1):
        for row in reversed(range(col, 5)):
            a, b = reduced[row : row + 2, col]
            if abs(b) < 1e-14:
                continue
            theta = np.arctan2(abs(b), abs(a))
            phase = argument(b) - argument(a)
            z = np.exp(1j * phase)
            c, s = np.cos(theta), np.sin(theta)
            reduced[row : row + 2] = (
                np.array([[c, z.conjugate() * s], [-z * s, c]]) @ reduced[row : row + 2]
            )
            pulses.append(word_rotation(row, theta, phase))
    phases = [argument(reduced[j, j]) for j in (0, 1)]
    i = next(i for i in range(3) if WORDS[0][i] != WORDS[1][i])
    fix = np.zeros(3)
    fix[WORDS[0][i]], fix[WORDS[1][i]] = phases
    return label_phase(i, fix) + [p for group in reversed(pulses) for p in group]


@lru_cache(None)
def local_h(k):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + k * contrast
    ids = [law.SRC.index(pair) for pair in product(range(3), repeat=2)]
    blank = np.eye(law.DIM)[:, ids]
    reduced = blank.T @ h @ blank
    np.testing.assert_allclose(h @ blank, blank @ reduced, atol=1e-13)
    return reduced


@lru_cache(None)
def gate(kind, t, k):
    if kind == "private":
        return np.diag(np.exp(-1j * t * k * np.array([0, 1, 4])))
    return expm(-1j * t * local_h(k))


def play(pulses, initial, k):
    state = initial.reshape((3,) * 6 + (-1,))
    board = federation(0)
    nodes = ports()
    for kind, i, j, time in pulses:
        if kind == "private":
            if board.degree[nodes[i]] != 5:
                raise ValueError("Unavailable private reference")
            axes = (i,)
        else:
            if nodes[j] not in board.neighbours[nodes[i]]:
                raise ValueError("Unavailable addressed contact")
            axes = (i, j)
        leading = tuple(range(len(axes)))
        moved = np.moveaxis(state, axes, leading)
        matrix = gate(kind, time, k)
        out = matrix @ moved.reshape(len(matrix), -1)
        state = np.moveaxis(out.reshape(moved.shape), leading, axes)
    return state.reshape(729, -1)


def prepare(name, k):
    source = np.zeros((729, 1), complex)
    source[INDEX[int(name == "one")], 0] = 1
    if name in ("zero", "one"):
        return source
    phi = {"plus_x": 0, "plus_y": np.pi / 2, "minus_y": -np.pi / 2}[name]
    return play(word_rotation(0, np.pi / 4, phi), source, k)


def read(state):
    mass = abs(state.flatten()) ** 2
    mask = np.ones(len(mass), bool)
    mask[INDEX] = False
    prob = np.r_[mass[INDEX], mass[mask].sum()]
    np.testing.assert_allclose(prob.sum(), 1, rtol=0, atol=2e-12)
    return prob
