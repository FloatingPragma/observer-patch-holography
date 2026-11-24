"""Response of candidate baryon lines to a local reference, and their radial material records."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import currents
from . import backgrounds, reservoir


def shells(weights, loads):
    length = loads.shape[1]
    mean = weights @ loads
    excess = (mean + np.roll(mean, 1) + np.roll(mean, -1) - weights.sum()) / 3
    result = []
    for centre in range(length):
        d = abs(np.arange(length) - centre)
        result.append(np.bincount(np.minimum(d, length - d), weights=excess))
    return np.array(result)


def panel(r, site, ext):
    length = r["length"]
    _, spokes, edges, faces, geo = reservoir.patch(length, site)
    words, material, entry, tight = reservoir.configurations(
        backgrounds.bits(r["gauss"], length), spokes
    )
    h = reservoir.generator(
        words, material, edges, faces, np.array(geo["reference"]), ext
    )
    loads = (material[:, None] >> np.arange(length)) & 1
    masks = [tight & (entry == key) for key in (-1, site)]
    result = dict(length=length, entry_site=site)
    if not all(m.any() for m in masks):
        return result | {"lines_available": False}
    e, v = eigh(h.toarray(), driver="evd", check_finite=False)
    winners = [backgrounds.peak(e, v, m) for m in masks]
    step = 1e-3
    response = []
    for scale in (-step, step):
        a = h.toarray()
        a[np.diag_indices_from(a)] += scale * (entry == -1)
        ev, vec = eigh(a, driver="evd", check_finite=False)
        match = []
        for _, state in winners:
            overlap = abs(state @ vec) ** 2
            best = np.argmax(overlap)
            match.append(dict(energy=ev[best], overlap=overlap[best]))
        response.append(match)
    energies = []
    slopes = []
    for i, ((energy, state), k, m) in enumerate(
        zip(winners, ("neutral", "positive"), masks)
    ):
        p = state**2
        avg = shells(p, loads)
        conditional = shells(p * m / sum(p[m]), loads)
        np.testing.assert_allclose(avg.sum(axis=1), sum(p[entry != -1]), atol=1e-10)
        np.testing.assert_allclose(conditional.sum(axis=1), i, atol=1e-10)
        slope = p[entry == -1].sum()
        diff = (response[1][i]["energy"] - response[0][i]["energy"]) / (2 * step)
        result[k] = dict(
            energy=energy,
            bulk_reservoir_probability=slope,
            measured_slope=diff,
            tracked=[a[i] for a in response],
            unconditional_all_centres=avg,
            tight_all_centres=conditional,
        )
        energies.append(energy)
        slopes.append(slope)
    return result | dict(
        lines_available=True,
        neutral_minus_positive=energies[0] - energies[1],
        difference_slope=slopes[0] - slopes[1],
        step=step,
    )


def run(scans):
    return currents._portable(
        {
            k: [
                panel(r, v["additional_site"], k == "extension")
                for r in v["selected"]
                if r["n"] == r["length"] // 3
            ]
            for k, v in scans["paired"].items()
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    scans = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    data = run(scans)
    with args.output.open("x") as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
