"""Build the moving ring and both active probes from physical contacts and trits."""

import json
from functools import cache
from itertools import combinations, product

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .certificate import HERE, apertures, blank_terms
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation


@cache
def initial_ring(k):
    n = 14
    words = np.array(sorted(sum(1 << i for i in c) for c in combinations(range(n), 8)))
    bits = (words[:, None] >> np.arange(n)) & 1
    unequal = bits != np.roll(bits, -1, axis=1)
    diag = (0.5 + k) * unequal.sum(axis=1)
    h = sparse.diags(diag, format="csr")
    for j in range(n):
        active = np.flatnonzero(unequal[:, j])
        other = np.searchsorted(
            words, words[active] ^ ((1 << j) | (1 << ((j + 1) % n)))
        )
        h += sparse.csr_matrix(
            (np.full(len(active), -0.5), (other, active)), shape=h.shape
        )
    psi = np.zeros(len(words))
    psi[np.searchsorted(words, 1 | sum(1 << j for j in range(1, n, 2)))] = 1
    psi = expm_multiply(-3j * h, psi, traceA=-3j * sum(diag))
    return bits, h, psi


def apparatus(k, src, probe, mode, mount):
    board = federation(0)
    ring, stations = apertures(0, 14)
    rails = (
        tuple(tuple(s[1:]) for s in stations)
        if mode != "detached"
        else tuple(map(tuple, mount["distant"]))
    )
    bits, h, psi = initial_ring(k)
    if src < 0:
        bits = 1 - bits
    bonds = tuple(s[0] for s in stations)
    signs = bits + np.roll(bits, -1, axis=1) - 1
    gate = ((signs != 0).sum(axis=1) == 2) & (signs[:, bonds] == src).all(axis=1)
    sites = {port: values[:, None] for port, values in zip(ring, bits.T, strict=True)}
    orientations = np.array(list(product(range(2), repeat=2)), int).T
    for (a, b), orient in zip(rails, orientations, strict=True):
        sites[a] = np.where(orient, 1 - probe, 1)[None, :]
        sites[b] = np.where(orient, 1, 1 - probe)[None, :]
    # Trace each physical edge and private zero reference. Retain full (n_p-n_q)^2.
    c = blank_terms()[1].diagonal().reshape(3, 3)
    diag = np.zeros((len(bits), 4))
    # The literal cut removes every term of each squared contrast at its fixed tip.
    cut = (
        {
            frozenset((station[1], ring[i % len(ring)]))
            for station in stations
            for i in (station[0], station[0] + 1)
        }
        if mode == "cut"
        else set()
    )
    for a, b in board.seams:
        if frozenset((a, b)) in cut:
            continue
        if a in sites or b in sites:
            diag += k * c[sites.get(a, 0), sites.get(b, 0)]
    for a, load in sites.items():
        diag += k * (6 - board.degree[a]) * c[load, 0]
    exchange = blank_terms()[0]
    aux = np.array(
        [candidate.native_law.SRC.index(w) for w in ((1, 1 - probe), (1 - probe, 1))]
    )
    plate = exchange[np.ix_(aux, aux)]
    T = sparse.kron(sparse.csr_matrix(plate), np.eye(2)) + sparse.kron(
        sparse.eye(2), plate
    )
    unequal = (bits != np.roll(bits, -1, axis=1)).sum(axis=1)
    H = sparse.kron(sparse.eye(len(bits)), T, format="csr")
    off = h - sparse.diags(h.diagonal())
    # In stopped mode the planned instrument parks source hopping at the given appointment.
    if mode != "stopped":
        H += sparse.kron(off, sparse.eye(4), format="csr")
    diag += unequal[:, None] / 2
    # A common scalar changes no local reading.
    diag -= diag.mean()
    return H + sparse.diags(diag.ravel()), np.kron(psi, np.ones(4) / 2), gate


def check(filename, modes):
    p = HERE.parent / "material_probe"
    raw = json.loads((p / filename).read_text())
    mount = next(
        s
        for s in json.loads((p / "counts.json").read_text())["mounts"]
        if s["level"] == 0 and s["length"] == 14
    )
    banks = {tuple(key): v for key, v in raw["banks"]}
    size = raw["shots"]
    labels = list(product((8, 4), (-1, 1), (-1, 1), modes))
    assert size == 2**28 and len(raw["banks"]) == len(banks) == 96
    eps = np.sqrt(np.log(2 * 2048 / 0.01) / (2 * size))
    worst = 0
    covered = 0
    total = 0
    q = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
    for k, src, s, mode in labels:
        h, psi, gate = apparatus(k, src, s, mode, mount)
        for i, time in enumerate((1 / 16, 1 / 4)):
            matrix = -1j * time * h
            out = expm_multiply(matrix, psi, traceA=matrix.diagonal().sum()).reshape(
                len(gate), 4
            )
            for name, field in (("rail", out), ("phase", out @ np.kron(q, q).T)):
                weights = abs(field[:, ::-1]) ** 2
                pred = np.concatenate(
                    [weights[~gate].sum(axis=0), weights[gate].sum(axis=0)]
                )
                counts = np.array(banks[k, src, s, mode, i, name])
                assert counts.sum() == size and len(counts) == 8 and min(counts) >= 0
                errors = abs(counts / size - pred)
                covered += int((errors <= eps).sum())
                total += len(counts)
                worst = max(worst, float(max(errors)))
    return dict(
        predictions=total,
        covered=covered,
        worst_absolute_difference=worst,
        simultaneous_half_width=eps,
    )


if __name__ == "__main__":
    for filename, modes in (
        ("live_counts.json", ("live", "detached", "stopped")),
        ("cut_counts.json", ("live", "detached", "cut")),
    ):
        result = check(filename, modes)
        print(json.dumps({filename: result}, indent=2))
        if result["covered"] != result["predictions"]:
            raise SystemExit("Moving-source probabilities failed held counting bounds")
