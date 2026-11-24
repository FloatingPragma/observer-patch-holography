"""Stationary composite and calibrated signed constituent reads on one reservoir graph."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import currents
from . import backgrounds, field_probe, joint_field, loop_colour, reservoir


def read(weights, cores, entry, site, edges, heads, ref, rng):
    result = {}
    for channel, select in (("neutral", entry == -1), ("positive", entry == site)):
        for tag in ("negative", "positive"):
            rows = [(i, c) for i, c in cores.items() if select[i]]
            p = sum(weights[i] for i, _ in rows)
            if p < 1e-14:
                result[channel + "/" + tag] = {"detection_probability": p}
                continue
            totals = np.zeros((2, 2, 4))
            blind_joint = np.zeros_like(totals)
            for i, items in rows:
                selected = [(f, q) for k, f, q in items if k == tag]
                for f, q in selected:
                    aperture = np.r_[f, 0]
                    grad = np.array([aperture[b] - aperture[a] for a, b in edges])
                    delta = grad @ (ref - heads[i])
                    np.testing.assert_allclose(delta, q, atol=1e-12)
                    weight = weights[i] / len(selected)
                    for j, sign in enumerate((-1, 1)):
                        for target, after in ((totals, heads[i]), (blind_joint, ref)):
                            trial = field_probe.edge_response(
                                ref, after, grad, np.pi / 2, sign
                            )
                            target[j] += weight * np.array(
                                [trial["quadratures"][b] for b in ("x", "y")]
                            )
            with_quiet = np.concatenate((totals, np.full((2, 2, 1), 1 - p)), axis=2)
            shots = 2**19
            result[channel + "/" + tag] = dict(
                detection_probability=p,
                signed_joint=with_quiet,
                phase_units=[
                    np.angle(complex(*(v[:, 0] - v[:, 1]))) / (sign * np.pi / 2)
                    for v, sign in zip(totals, (-1, 1))
                ],
                head_independent_imaginary=blind_joint[:, 1, 0] - blind_joint[:, 1, 1],
                counts=[
                    [
                        rng.multinomial(
                            shots, np.maximum(a, 0) / np.maximum(a, 0).sum()
                        ).tolist()
                        for a in y
                    ]
                    for y in with_quiet
                ],
            )
    return result


def panel(r, site, extension, rng):
    length = r["length"]
    z, spokes, edges, faces, geo = reservoir.patch(length, site)
    words, material, entry, tight = reservoir.configurations(
        backgrounds.bits(r["gauss"], length), spokes
    )
    h = reservoir.generator(
        words, material, edges, faces, np.array(geo["reference"]), extension
    )
    bits = (material[:, None] >> np.arange(length)) & 1
    heads = (words[:, None] >> np.arange(len(edges))) & 1
    patch = SimpleNamespace(bit={edge: 1 << i for i, edge in enumerate(edges)})
    sizes, _, _, components = loop_colour.refine(
        patch, words, material, include_witness=False
    )
    result = dict(
        length=length,
        reservoir_port=z,
        indistinguishability=sizes,
        hop_components=components,
    )
    ref = None
    for phase in range(3):
        ids = np.flatnonzero(
            (entry == -1) & (bits == (np.arange(length) % 3 == phase)).all(axis=1)
        )
        if len(ids):
            ref = int(ids[0])
            break
    if (
        ref is None
        or not any(tight & (entry == -1))
        or not any(tight & (entry == site))
    ):
        return result | {"common_field_read_available": False}
    result.update(
        common_field_read_available=True,
        reference_phase=phase,
        reference_heads=int(words[ref]),
    )
    far = np.zeros(len(words), bool)
    nonwall = np.zeros_like(far)
    for n in np.unique(bits.sum(axis=1)):
        mask = bits.sum(axis=1) == n
        chart = backgrounds.arcs(length, int(n))
        a = chart[(material[mask] & ((1 << length) - 1))]
        far[mask] = a >= (np.nanmin(chart) + np.nanmax(chart)) / 2
        nonwall[mask] = ~np.isfinite(a)
    span = np.linalg.eigvalsh(h[far][:, far].toarray())
    trials = h.diagonal()[nonwall]
    result.update(
        stretched_interval=[span[0], span[-1]],
        nonwall_diagonal_interval=[min(trials), max(trials)],
    )
    e, v = eigh(h.toarray(), driver="evd", check_finite=False)
    masks = [(entry == -1) & tight, (entry == site) & tight]
    winners = [backgrounds.peak(e, v, m) for m in masks]
    scores = np.minimum(*(np.sum(v[m] ** 2, axis=0) for m in masks))
    i = int(np.flatnonzero(scores >= scores.max() - 1e-12)[0])
    winners.append((e[i], v[:, i]))
    cores = {int(i): list(joint_field.cores(bits[i])) for i in np.flatnonzero(tight)}
    for key, (energy, state) in zip(("neutral", "positive", "compromise"), winners):
        weights = state**2
        result[key] = dict(
            energy=energy,
            residual=np.linalg.norm(h @ state - energy * state),
            full_read=reservoir.outcomes(weights, entry, tight, site),
            field=read(weights, cores, entry, site, edges, heads, heads[ref], rng),
        )
    return result


def run(scans):
    rng = np.random.default_rng(12013004)
    return currents._portable(
        {
            k: [
                panel(r, v["additional_site"], k == "extension", rng)
                for r in v["selected"]
                if r["n"] == r["length"] // 3
            ]
            for k, v in scans["paired"].items()
        }
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    d = run(s)
    with args.output.open("x") as f:
        f.write(json.dumps(d, indent=2, allow_nan=False) + "\n")
