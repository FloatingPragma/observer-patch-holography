"""Native 012 junction, independently signed six-position reduction, and spectator."""

from itertools import permutations, product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

LAW = candidate.native_law
WORDS = tuple(sorted(set(permutations((0, 1, 1, 2)))))
START = (0, 1, 1, 2)
LEFT, RIGHT = (1, 3, 2), (2, 3, 1)
BONDS = (1, 2, 3)
# Cyclic order of unordered core positions, starting from the two occupied leaves.
RING = ((1, 2), (0, 2), (2, 3), (0, 3), (1, 3), (0, 1))
TIMES = (1, 8, 32, 64, 128)
SHOTS = 2**24


def swap(word, i, j):
    out = list(word)
    out[i], out[j] = out[j], out[i]
    return tuple(out)


def ports():
    board = federation(0)
    # All other contacts, including any bonds among these leaves, are shuttered.
    return (0, *board.neighbours[0][:4])


def positions(word):
    return tuple(i for i, load in enumerate(word) if load == 1)


def word_matrix(k, bonds=BONDS, kinetic=True):
    local = LAW.seam_term(k)
    if not kinetic:
        local = local - LAW.seam_term(0)
    h = np.zeros((len(WORDS), len(WORDS)))
    for col, w in enumerate(WORDS):
        for b in bonds:
            p = (w[0], w[b])
            a = LAW.SRC.index(p)
            h[col, col] += local[a, a]
            if p != p[::-1]:
                h[WORDS.index(swap(w, 0, b)), col] += local[LAW.SRC.index(p[::-1]), a]
    return h


def sharp():
    return np.eye(len(WORDS))[:, WORDS.index(START)]


def pulse(h):
    return expm(-1j * np.pi * h)


def prepare(sign, k=8):
    # The reference is an idle exposure of equal duration. The contrast phase
    # on the 02 word is exp(-4 pi i k), exactly one at both enrolled contrasts.
    src = sharp()
    projected = (src + sign * (pulse(word_matrix(k, (3,))) @ src)) / 2
    prob = float(np.vdot(projected, projected).real)
    return prob, projected / np.sqrt(prob)


def evolved(src, sequence, k=8):
    for b in sequence:
        src = pulse(word_matrix(k, (b,))) @ src
    return src


def dark(left, right):
    return float(np.linalg.norm(left - right) ** 2 / 4)


def permutation(seq, n=4):
    labels = tuple(range(n))
    for b in seq:
        labels = swap(labels, 0, b)
    return labels


def addressed():
    banks = [np.random.default_rng(seed) for seed in (702, 1702)]
    spin_flip = np.array([np.exp(-1j * np.pi * w[0] ** 2 / 4) for w in WORDS])
    result = {}
    for sign in (1, -1):
        p, src = prepare(sign)
        closed = LEFT + RIGHT[::-1]
        retraced = LEFT + LEFT[::-1]
        left, right = evolved(src, LEFT), evolved(src, RIGHT)
        result[str(sign)] = dict(
            source_probability=p,
            at_common_endpoint=float(
                sum(
                    abs(left[i]) ** 2
                    for i, w in enumerate(WORDS)
                    if positions(w) == (0, 3)
                )
            ),
            half_paths_dark=dark(left, right),
            retraced_dark=dark(src, evolved(src, retraced)),
            once_dark=dark(evolved(src, closed), evolved(src, retraced)),
            twice_dark=dark(evolved(src, closed * 2), evolved(src, retraced * 2)),
            fixed_cores_host_swap_dark=dark(src, evolved(src, (3,))),
            matched_host_swap_dark=dark(evolved(src, closed), evolved(src, (3,))),
            reads_by_step=[
                [
                    list(positions(WORDS[int(np.argmax(abs(evolved(src, seq[:i]))))]))
                    for i in range(len(seq) + 1)
                ]
                for seq in (LEFT, RIGHT)
            ],
            # A diagnostic phase plate on a host, never used on the native flight.
            changed_host_preparation_dark=dark(
                evolved(spin_flip * src, LEFT), evolved(spin_flip * src, RIGHT)
            ),
        )
    result["sharp_unheralded_dark"] = dark(
        evolved(sharp(), LEFT), evolved(sharp(), RIGHT)
    )
    result["half_path_tags"] = [permutation(LEFT), permutation(RIGHT)]
    result["closed_tags"] = permutation(LEFT + RIGHT[::-1])
    result["joint_launch_counts"] = {}
    for name in (
        "half_paths_dark",
        "retraced_dark",
        "once_dark",
        "twice_dark",
        "fixed_cores_host_swap_dark",
        "matched_host_swap_dark",
        "changed_host_preparation_dark",
    ):
        law = []
        for branch in ("1", "-1"):
            p = result[branch]["source_probability"]
            dark_exit = float(np.clip(result[branch][name], 0, 1))
            law.extend((p * (1 - dark_exit), p * dark_exit))
        result["joint_launch_counts"][name] = []
        for bank in banks:
            result["joint_launch_counts"][name].append(
                bank.multinomial(SHOTS, probabilities(law)).tolist()
            )
    return result


def effective(sign, k):
    # Compare against only the core coordinates, with the connection placed at
    # an arbitrary cut. All dynamical phases from the host contact are retained.
    v = [2 * k + 1 if 0 in pair else 6 * k + 1 + (1 - sign) / 2 for pair in RING]
    h = np.diag(v)
    for i in range(6):
        j = (i + 1) % 6
        h[i, j] = h[j, i] = -0.5 * (sign if j == 0 else 1)
    return h


def law_from_words(sign, k, t, kinetic=True):
    _, src = prepare(sign, k)
    wave = expm(-1j * t * word_matrix(k, kinetic=kinetic)) @ src
    return np.array(
        [
            sum(abs(wave[j]) ** 2 for j, w in enumerate(WORDS) if positions(w) == pair)
            for pair in RING
        ]
    )


def probabilities(p):
    p = np.asarray(p, dtype=float)
    if min(p) < -1e-11 or abs(sum(p) - 1) > 1e-9:
        raise ArithmeticError("Invalid normalized detector law")
    return p.clip(0) / p.clip(0).sum()


def release():
    records = []
    banks = [np.random.default_rng(seed) for seed in (701, 1701)]
    for k, sign, t in product((8, 4), (1, -1), TIMES):
        p = probabilities(law_from_words(sign, k, t))
        q = abs(expm(-1j * t * effective(sign, k))[:, 0]) ** 2
        if np.max(abs(p - q)) > 1e-10:
            raise ArithmeticError(
                "Free position law disagrees with the independently signed ring"
            )
        sampled = [rng.multinomial(SHOTS, p).tolist() for rng in banks]
        records.append(
            dict(k=k, host_sign=sign, time=t, law=p.tolist(), counts=sampled)
        )
    return records


def spectator():
    # Three two-level hosts in port order O,C,D; the two unit loads start on A,B.
    host_words = tuple(product((0, 2), repeat=3))
    host_swaps = []
    routes = []
    for spare in (3, 4):
        route = (1, spare, 2, 1, spare, 2)
        tags = permutation(route, n=5)
        routes.append(tags)
        m = np.zeros((8, 8))
        host_ports = (0, 3, 4)
        for j, w in enumerate(host_words):
            final = tuple(w[host_ports.index(tags[p])] for p in host_ports)
            m[host_words.index(final), j] = 1
        host_swaps.append(m)
    identity = np.eye(8)
    exits = [(identity - u) / 2 for u in host_swaps]
    # Kernel intersection by full matrix diagonalization. Also test an extra
    # sharp 0 at D next to the four-site odd source, with both paths available.
    odd = np.zeros(8)
    odd[host_words.index((0, 2, 0))] = 1 / np.sqrt(2)
    odd[host_words.index((2, 0, 0))] = -1 / np.sqrt(2)
    return dict(
        route_tags=routes,
        odd_source_with_spectator=[float(odd @ m @ odd) for m in exits],
        mean_exit_spectrum=np.linalg.eigvalsh(sum(exits) / 2).tolist(),
    )


def run():
    return dict(
        ports=list(ports()),
        core_word_order=RING,
        addressed=addressed(),
        free=release(),
        free_parked=law_from_words(-1, 8, 128, kinetic=False).tolist(),
        simultaneous_spectator=spectator(),
        trials_per_free_bank=SHOTS,
    )
