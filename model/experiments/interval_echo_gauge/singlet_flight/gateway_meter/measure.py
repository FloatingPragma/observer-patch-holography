"""Read gateway overlaps by sequential time comparisons and a root-gated loop."""

from collections import Counter
from itertools import combinations, permutations
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation

HERE = Path(__file__).resolve().parent
STARTS = 2**20
SEED = 150177
SCALES = {"native": 1.0, "half_exchange": 0.5}
ACQUISITION = HERE / "acquisition.json.gz"
# The normalization tolerance of every detector law.
PROBABILITY_TOLERANCE = 1e-12


def incidence(car):
    a = car.adjacency().toarray().astype(int)
    identity = np.eye(12, dtype=int)
    # The unique unreachable address in two steps is the antipode.
    j = ((identity + a + a @ a) == 0).astype(int)
    b = identity - j
    c = b @ a
    ones = np.ones_like(a)
    facts = [
        (j @ j, identity),
        (a @ j, j @ a),
        (b @ b, 2 * b),
        (b @ c, 2 * c),
        (c @ c, 10 * b),
        (a @ c, 5 * b),
        ((identity + j) @ (a + identity), ones),
    ]
    for left, right in facts:
        np.testing.assert_array_equal(left, right)
    assert np.trace(b) == 12 and np.trace(c) == 0
    assert np.all((c == 0) == (b != 0)) and set(c.flat) == {-1, 0, 1}
    return b, c


def mounted_h(car, scale):
    swap, contrast = law._load_ops()
    seam = scale * (np.eye(9) - swap) / 2 + candidate.KAPPA * contrast
    h = np.zeros((12, 12))
    for x, y in car.seams:
        h[np.ix_([x, y], [x, y])] += seam[np.ix_([3, 1], [3, 1])]
    return h


def bright_filters(car, scale):
    h = mounted_h(car, scale)
    ident = np.eye(12)
    offset = 5 * candidate.KAPPA
    energies = offset + (5 - np.array([5, np.sqrt(5), -1, -np.sqrt(5)])) / 2
    target = energies[1]
    operators = []
    for other in energies[[0, 2, 3]]:
        delay = np.pi / abs(target - other)
        operators.append((ident + expm(-1j * (h - target * ident) * delay)) / 2)
    return operators


def clock_probs(filters, state):
    alive = state
    probabilities = []
    for bright in filters:
        before = np.vdot(alive, alive).real
        alive = bright @ alive
        probabilities.append(before - np.vdot(alive, alive).real)
    probabilities.append(np.vdot(alive, alive).real)
    return normalized(probabilities)


def normalized(probabilities):
    p = np.asarray(probabilities)
    if min(p) < -1e-12 or abs(p.sum() - 1) > 1e-12:
        raise ArithmeticError("detector probabilities must cover all launches")
    p = np.maximum(p, 0)
    return p / p.sum()


def loop_probs(scale):
    # Local basis contains every arrangement of the vacancy, excess and host.
    words = list(permutations((0, 1, 2)))
    ident = np.eye(len(words))
    loop = ident.copy()
    for a, b in [(0, 1), (1, 2), (2, 0)]:
        s = np.zeros_like(ident)
        costs = []
        for i, word in enumerate(words):
            updated = list(word)
            updated[a], updated[b] = word[b], word[a]
            s[words.index(tuple(updated)), i] = 1
            costs.append((word[a] - word[b]) ** 2)
        contact = scale * (ident - s) / 2 + candidate.KAPPA * np.diag(costs)
        loop = expm(-1j * np.pi * contact) @ loop
    gate = np.diag([int(w[0] == 0) for w in words])
    vectors = np.eye(len(words))
    u, v, off = [vectors[words.index(w)] for w in ((0, 2, 1), (0, 1, 2), (2, 0, 1))]
    settings = [(u + v) / np.sqrt(2), (u - v) / np.sqrt(2), u, off]
    values = []
    for state in settings:
        entered = gate @ state
        dark = np.linalg.norm((ident - loop) @ entered / 2) ** 2
        prob = np.vdot(entered, entered).real
        values.append(normalized([dark, prob - dark, 1 - prob]))
    return values


def board_sites(board, car):
    sites = [[] for _ in range(board.cells)]
    for edge in board.external:
        for port in edge:
            sites[port // 12].append(port % 12)
    for owner, site in enumerate(sites):
        site.sort()
        assert len(set(site)) == 3
        for port in range(12):
            neighbours = board.neighbours[12 * owner + port]
            assert (
                sorted(p % 12 for p in neighbours if p // 12 == owner)
                == car.neighbours[port]
            )
    return sites


def gateway_states(ports):
    eye = np.eye(12)
    states = list(eye[ports])
    for x, y in combinations(ports, 2):
        states.extend([(eye[x] + eye[y]) / np.sqrt(2), (eye[x] - eye[y]) / np.sqrt(2)])
    return states


def acquisition_law(fresh, stored):
    stored = np.asarray(stored, dtype=float)
    gap = np.max(np.abs(fresh - stored)) if stored.shape == fresh.shape else np.inf
    if gap > PROBABILITY_TOLERANCE:
        raise ArithmeticError("detector probabilities differ from acquisition law")
    return stored


def acquisition():
    return json.loads(gzip.decompress(ACQUISITION.read_bytes()))


def replay(acquired=None):
    """Return the seeded counts and the detector laws they were drawn from.

    With an acquisition record the freshly computed laws must agree with the
    stored binary64 probabilities, and the counts are drawn from the stored ones.
    """
    car = cell()
    incidence(car)
    np.testing.assert_allclose(
        mounted_h(car, 1) + candidate.KAPPA * np.eye(12),
        candidate.response(car),
        atol=1e-13,
        rtol=0,
    )
    native_g0, native_g1 = law.local_term_parts()
    native = native_g0 + candidate.KAPPA * native_g1
    np.testing.assert_allclose(
        native[:, :9],
        np.vstack([law.seam_term(candidate.KAPPA), np.zeros((90, 9))]),
        atol=1e-13,
        rtol=0,
    )
    rng = np.random.default_rng(SEED)
    records = {"seed": SEED, "starts_per_setting": STARTS, "boards": []}
    laws = {"seed": SEED, "starts_per_setting": STARTS, "boards": []}
    for level in range(3):
        board = federation(level)
        sites = board_sites(board, car)
        faces = sorted(board.faces())
        assert len(faces) == 20 * board.cells
        bank = {"owners": board.cells, "gateways": sites, "arms": {}}
        laws["boards"].append({"arms": {}})
        for name, scale in SCALES.items():
            filters = bright_filters(car, scale)
            loop = np.array(loop_probs(scale))
            gateway = np.array(
                [[clock_probs(filters, s) for s in gateway_states(p)] for p in sites]
            )
            laws["boards"][-1]["arms"][name] = {
                "gateway_probabilities": gateway.tolist(),
                "loop_probabilities": loop.tolist(),
            }
            if acquired is not None:
                stored = acquired["boards"][level]["arms"][name]
                gateway = acquisition_law(gateway, stored["gateway_probabilities"])
                loop = acquisition_law(loop, stored["loop_probabilities"])
            # Every consecutive group of three root settings follows the sorted face.
            rooted = [[i, *sorted(set(face) - {i})] for face in faces for i in face]
            face_counts = np.stack(
                [rng.multinomial(STARTS, p, size=len(rooted)) for p in loop], axis=1
            )
            owners = [[rng.multinomial(STARTS, p).tolist() for p in o] for o in gateway]
            bank["arms"][name] = {
                "gateway_counts": owners,
                "loop_counts": face_counts.tolist(),
            }
        bank["rooted_faces"] = rooted
        records["boards"].append(bank)
    return records, laws


def decode(data):
    b, c = incidence(cell())
    target = (b + c / np.sqrt(5)) / 4
    report = {}
    for name in SCALES:
        errors, shifts = [], []
        classes = Counter()
        root_sums = np.zeros((4, 3), dtype=np.int64)
        for board in data["boards"]:
            arm = board["arms"][name]
            for ports, counts in zip(
                board["gateways"], arm["gateway_counts"], strict=True
            ):
                obs = np.array(counts)[:, -1] / STARTS
                for p, reading in zip(ports, obs[:3], strict=True):
                    shifts.append(abs(reading - target[p, p]))
                for i, (p, q) in enumerate(combinations(ports, 2)):
                    signed = (obs[3 + 2 * i] - obs[4 + 2 * i]) / 2
                    errors.append(abs(signed - target[p, q]))
                    classes[str((int(b[p, q]), int(c[p, q])))] += 1
            root_sums += np.asarray(arm["loop_counts"]).sum(axis=0)
        report[name] = {
            "max_pair_error": max(errors),
            "max_diagonal_error": max(shifts),
            "pair_classes": dict(classes),
            "summed_loop_outcomes": root_sums.tolist(),
        }
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "counts.json.gz"
    if args.write:
        records, laws = replay()
        for target, payload in ((path, records), (ACQUISITION, laws)):
            encoded = json.dumps(payload, separators=(",", ":"), allow_nan=False)
            target.write_bytes(gzip.compress(encoded.encode(), mtime=0))
    else:
        records, _ = replay(acquisition())
        assert records == json.loads(gzip.decompress(path.read_bytes()))
    print(json.dumps(decode(records), indent=2))
