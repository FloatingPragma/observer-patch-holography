"""Independent assembly of all sharp calibration and joint probe probabilities."""

import json
from functools import cache
from itertools import combinations, product

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.special import betainccinv, betaincinv

from .certificate import BOARDS, HERE, apertures, blank_terms, rail_matrix
from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate

READS = ("phase", "rail")
STRENGTHS = {
    "native": (8, 8),
    "distant": (8, 8),
    "uncoupled": (8, 0),
    "flight_off": (0, 8),
    "half": (4, 4),
}


@cache
def propagated(n, amplitude, strength):
    """From bit words and raw nearest-neighbour swaps, then independent Krylov flight."""
    words = np.array(
        sorted(sum(1 << j for j in c) for c in combinations(range(n), n // 2 + 1))
    )
    idx = np.arange(len(words))
    diag = np.zeros(len(words))
    rows, cols, values = [], [], []
    for j in range(n):
        k = (j + 1) % n
        unequal = ((words >> j) ^ (words >> k)) & 1
        active = idx[unequal.astype(bool)]
        rows.extend(active)
        cols.extend(np.searchsorted(words, words[active] ^ ((1 << j) | (1 << k))))
        values.extend(np.full(len(active), -0.5))
        diag += (0.5 + strength * amplitude**2) * unequal
    h = sparse.csr_matrix(
        (values, (rows, cols)), shape=(len(idx), len(idx))
    ) + sparse.diags(diag)
    psi = np.zeros(len(idx))
    seed = 1 | sum(1 << j for j in range(1, n, 2))
    psi[np.searchsorted(words, seed)] = 1
    amplitudes = expm_multiply(-3j * h, psi, traceA=-3j * diag.sum())
    rates = abs(amplitudes) ** 2
    np.testing.assert_allclose(rates.sum(), 1, rtol=0, atol=2e-11)
    return words, rates


@cache
def source_probabilities(level, n, amount, sign, flight):
    words, probs = propagated(n, amount, flight)
    if sign < 0:
        words = ((1 << n) - 1) ^ words
    _, sensors = apertures(level, n)
    bits = (words[:, None] >> np.arange(n)) & 1
    walls = bits + np.roll(bits, -1, axis=1) - 1
    left, right = (s[0] for s in sensors)
    herald = (
        ((walls != 0).sum(axis=1) == 2)
        & (walls[:, left] == sign)
        & (walls[:, right] == sign)
    )
    sums = bits + np.roll(bits, -1, axis=1)
    return np.bincount(
        9 * herald + 3 * sums[:, left] + sums[:, right], weights=probs, minlength=18
    ).reshape(2, 3, 3)


@cache
def probability(level, n, side, x, y, probe_sign, kappa, time, setting, distant=None):
    ring, local = apertures(level, n)
    board = federation(level)
    j, a, b = local[side]
    stock = {ring[j]: x, ring[(j + 1) % n]: y}
    if distant is not None:
        a, b = distant
    matrix = rail_matrix(board, (a, b), stock, probe_sign, kappa)
    psi = expm(-1j * matrix * ((1, 2)[time] / 32)) @ (np.ones(2) / np.sqrt(2))
    if setting == "phase":
        idx = [
            candidate.native_law.SRC.index(w)
            for w in ((1, 1 - probe_sign), (1 - probe_sign, 1))
        ]
        exchange, _ = blank_terms()
        analyzer = expm(-1j * np.pi * exchange[np.ix_(idx, idx)] / 2)
        psi = analyzer @ psi
    return abs(psi[0]) ** 2


def check(raw):
    size = 2**24
    if raw["shots"] != size or raw["seed"] != 26011262:
        raise ValueError("Unexpected enrollment")
    if len(raw["mounts"]) != 4:
        raise ValueError("Missing instrument")
    remote = {}
    for (lvl, n), mount in zip(BOARDS, raw["mounts"], strict=True):
        ring, stations = apertures(lvl, n)
        assert (lvl, n) == (mount["level"], mount["length"])
        assert (
            tuple(mount["ring"]) == ring
            and tuple(map(tuple, mount["contacts"])) == stations
        )
        board = federation(lvl)
        far = mount["distant"]
        remote[lvl, n] = tuple(map(tuple, far))
        for a, b in far:
            assert a // 12 == b // 12 and a in board.neighbours[b]
            assert not set(ring).intersection(
                set(board.neighbours[a]) | set(board.neighbours[b]) | {a, b}
            )
        assert not any(v in board.neighbours[u] for u in far[0] for v in far[1])

    def index(rows, keys, width):
        data = {tuple(k): np.array(v) for k, v in rows}
        if len(rows) != len(keys) or set(data) != keys:
            raise ValueError("Missing or duplicate bank")
        for v in data.values():
            if (
                len(v) != width
                or not np.issubdtype(v.dtype, np.integer)
                or min(v) < 0
                or v.sum() != size
            ):
                raise ValueError("Invalid bank")
        return data

    prep = index(
        raw["standard"],
        {
            (*g, *k)
            for g in BOARDS
            for k in product(
                (4, 8), range(2), range(3), range(3), (-1, 1), range(2), READS
            )
        },
        2,
    )
    held = index(
        raw["measurements"],
        {
            (*g, *k)
            for g in BOARDS
            for k in product(STRENGTHS, (1, 2), (-1, 1), (-1, 1), range(2), READS)
        },
        8,
    )
    dataset = []
    for (lvl, n, k, s, x, y, sg, t, read), counts in prep.items():
        p = probability(lvl, n, s, x, y, sg, k, t, read)
        dataset.append((counts, np.array([1 - p, p])))
    for (lvl, n, law, amount, src, sg, t, read), counts in held.items():
        flight, det = STRENGTHS[law]
        mass = source_probabilities(lvl, n, amount, src, flight)
        byside = []
        for side in (0, 1):
            qs = []
            for a, b in ((0, 0), (0, amount), (amount, amount)):
                p = probability(
                    lvl,
                    n,
                    side,
                    a,
                    b,
                    sg,
                    det,
                    t,
                    read,
                    distant=remote[lvl, n][side] if law == "distant" else None,
                )
                qs.append([1 - p, p])
            byside.append(qs)
        prob = np.einsum("hab,ai,bj->hij", mass, *byside).ravel()
        dataset.append((counts, prob))
    counts = np.concatenate([r[0] for r in dataset])
    p = np.concatenate([r[1] for r in dataset])
    tail = 0.01 / (2 * 32768)
    lower = np.zeros(len(p))
    upper = np.ones(len(p))
    i = counts > 0
    j = counts < size
    lower[i] = betaincinv(counts[i], size - counts[i] + 1, tail)
    upper[j] = betainccinv(counts[j] + 1, size - counts[j], tail)
    covered = (lower <= p + 1e-12) & (p <= upper + 1e-12)
    return dict(
        calibration_banks=len(prep),
        joint_banks=len(held),
        predictions=len(counts),
        covered=int(sum(covered)),
        worst_absolute_difference=float(max(abs(counts / size - p))),
    )


if __name__ == "__main__":
    report = check(json.loads((HERE.parent / "material_probe/counts.json").read_text()))
    print(json.dumps(report, indent=2))
    if report["covered"] != report["predictions"]:
        raise SystemExit("Observed bank outside the simultaneous counting bounds")
