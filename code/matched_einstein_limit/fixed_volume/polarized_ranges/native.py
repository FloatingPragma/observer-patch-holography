"""Cross-check measured angular word positions against reversible exchange appointments."""
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.equilibrium_refinement import measure as source
from model.experiments.interval_echo_einstein_station.proper_ball_ruler.measure import pulse_calibration


DIRECTORY = Path(source.__file__).parent / "tracer_flux/bearings"


def routes(adjacent, words, start):
    paths = {start: ()}
    queue = deque([start])
    while queue:
        i = queue.popleft()
        for j in adjacent[i]:
            if j not in paths:
                edge = tuple(sorted(set(words[i]) ^ set(words[j])))
                paths[j] = (*paths[i], edge)
                queue.append(j)
    return paths


def timing_read(reg, path, start, end, gateway, exterior):
    owner = {0}
    reg.loads.fill(0)
    reg.loads[list(start)] = 2
    for a, b in path:
        reg.swap(a, b, owner)
    # Two resident flags at the nominated endpoint give its full sharp word
    # since the register always contains exactly two residents.
    for p in end:
        steps = reg.route(p, gateway, owner)
        for a, b in steps:
            reg.swap(a, b, owner)
        if not reg.flag(gateway, exterior, {0, exterior // 12}):
            raise ValueError("The timed transfer misses its registered destination")
        for a, b in reversed(steps):
            reg.swap(a, b, owner)
    for a, b in reversed(path):
        reg.swap(a, b, owner)
    if tuple(np.flatnonzero(reg.loads)) != start:
        raise ValueError("Inverse transfer did not restore the source word")


def gram_coefficients():
    """Integer parts A+B/sqrt(5) of the counted angular Gram matrix."""
    records = json.loads((DIRECTORY / "results.json").read_text())
    labels = records["anonymous_labels"]
    exposure = records["exposures"][0]
    bank = np.array(exposure["counts"]["uniform_ball"]["counts"])[:, :12]
    edges = np.zeros((12, 12), dtype=int)
    edges[np.ix_(labels, labels)] = bank / bank.sum(axis=1)[:, None] > exposure["decision_probability"]
    np.fill_diagonal(edges, 0)
    antipodes = np.linalg.matrix_power(edges, 2) == 0
    eye = np.eye(12, dtype=int)
    return eye - antipodes, 2 * edges + eye + antipodes - 1


def perpendicular(pair, other, gram):
    a, b = pair[0], other[0]
    return all(m[np.ix_(a, b)].sum() == 0 for m in gram)


def run():
    board = federation(0)
    ports, contacts, border, words = source.apparatus(board, [0])
    h = source.generator(board, words, contacts)
    assert np.all(np.diag(h) >= 0) and set(h[h < 0]) == {-.5}
    idx = {word: i for i, word in enumerate(words)}
    adjacent = [np.flatnonzero(row < 0) for row in h]
    paths = [routes(adjacent, words, i) for i in range(len(words))]

    def distance(a, b):
        return len(paths[idx[a]][idx[b]])

    placement = json.loads((DIRECTORY / "placements.json").read_text())
    origins = [tuple(v) for v in placement["shared_origin_records"]]
    axes = [(tuple(a), tuple(b)) for row in placement["shells"].values()
            for a, b in row["antipodal_word_pairs"]]
    gram = gram_coefficients()
    reader = Registers(board, seed=source.ROOT_SEED)
    gateway = min(border)
    surveyed = set()
    for a, b in combinations(sorted({*origins, *(p for axis in axes for p in axis)}), 2):
        timing_read(reader, paths[idx[a]][idx[b]], a, b, gateway, border[gateway])
        surveyed.add((a, b))

    centers = []
    for origin in origins:
        mid = [(a, b) for a, b in axes if distance(origin, a) == distance(origin, b)
               and distance(a, b) == 2 * distance(origin, a)]
        triples = []
        for triple in combinations(mid, 3):
            if not all(perpendicular(a, b, gram) for a, b in combinations(triple, 2)):
                continue
            rad = [distance(origin, a) for a, _ in triple]
            contractions, pairs, forks = [], [], []
            for i, j in combinations(range(3), 2):
                seg = [distance(p, q) for p, q in product(triple[i], triple[j])]
                r, s = rad[i], rad[j]
                # Squared-range average cancels the independent leading bearing.
                x = F(r**2 + s**2) - F(sum(n*n for n in seg), 4)
                k = 3 * x / (r*r*s*s)
                contractions.append(k)
                pairs.append(dict(axes=[i, j], ranges=seg, contraction=str(k)))
                a, b, c, d = seg
                if (a == b == r+s or c == d == r+s or a == c == r+s or b == d == r+s):
                    forks.append([i, j])
            triples.append(dict(axes=triple, radii=rad, cross_pairs=pairs,
                                fork_pairs=forks, scalar=str(2 * sum(contractions))))
        centers.append(dict(origin=origin, antipodal_pairs_at_true_midpoint=[
            dict(radius=rad, pairs=n) for rad,n in sorted(Counter(distance(origin, a) for a,b in mid).items())],
            orthogonal_triads=triples))
    return dict(source_words=len(words), native_exchange_hand=pulse_calibration()["decoded_single_step"],
                distinct_roundtrips=len(surveyed), read_gates=dict(reader.gates), centers=centers)
