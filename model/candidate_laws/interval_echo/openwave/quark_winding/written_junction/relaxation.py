"""Exact short-time screening, finite recurrent averages and stationary pair weights."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import currents
from . import ring


def panel(length, reference, extension, stationary):
    j = np.arange(length)
    seed = (j % 3 == 0).astype(int)
    gauss = np.roll(seed, 1)
    w, e, n = ring.basis(gauss)
    h = ring.generator(w, n, reference, extension)
    mask = sum(1 << int(i) for i in j[seed == 1])
    idx = np.searchsorted(w, [mask, 0, (1 << length) - 1])
    source = np.zeros((len(w), 3))
    source[idx, np.arange(3)] = 1
    head = e.mean(axis=1)

    def difference(values):
        return [values[1] - values[0], values[2] - values[0]]

    coefficients = [
        sum(
            abs(h[k, i]) ** 2 * (head[k] - head[i])
            for k in h.getrow(i).indices
            if k != i
        )
        for i in idx
    ]
    psi = source.astype(complex)
    previous = 0.0
    reads = []
    for t in (0.0, 1 / 64, 1 / 16, 1 / 8, 0.25, 1.0, 4.0):
        if t > previous:
            psi = expm_multiply(-1j * (t - previous) * h, psi)
        previous = t
        reads.append(dict(time=t, charges=difference(head @ (abs(psi) ** 2))))
    result = dict(
        length=length,
        dimension=len(w),
        law="positive_laplacian" if extension else "native",
        supplied_uniform_reference=not stationary,
        quadratic=difference(coefficients),
        reads=reads,
    )
    if stationary:
        levels, v = np.linalg.eigh(h.toarray())
        groups = np.split(
            np.arange(len(levels)), np.where(np.diff(levels) > 1e-8)[0] + 1
        )
        # Keep coherence between all states of equal energy.
        mean = np.zeros(3)
        for ix in groups:
            if len(ix) == 1:
                col = v[:, ix[0]]
                mean += (head @ (col**2)) * col[idx] ** 2
            else:
                weight = v[:, ix] @ v[idx][:, ix].T
                mean += head @ (abs(weight) ** 2)
        result["infinite_time_mean_difference"] = difference(mean)
        gaps = []
        for config in n:
            sites = np.flatnonzero(config)
            gaps.append(tuple(sorted(np.diff(np.r_[sites, sites[0] + length]))))
        maxima = {}
        for pair in ((2, 4), (1, 5)):
            read = np.array(
                [s == tuple(sorted([*pair] + [3] * (length // 3 - 2))) for s in gaps]
            )
            projected = v[read]
            peak = np.max(np.sum(abs(projected) ** 2, axis=0))
            for ix in groups:
                if len(ix) > 1:
                    peak = max(peak, np.linalg.norm(projected[:, ix], 2) ** 2)
            maxima[str(pair)] = peak
        result["stationary_pair_maximum"] = maxima
    return result


def run():
    board, cycles = ring.physical_cycles()
    return currents._portable(
        [
            panel(
                l,
                6 - board.degree[cycles[l]] if l in cycles else np.zeros(l, int),
                extension,
                l in cycles,
            )
            for l in (12, 15, 18, 21, 24)
            for extension in (False, True)
        ]
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
