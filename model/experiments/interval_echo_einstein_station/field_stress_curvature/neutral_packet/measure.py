"""Read the interference of two neutral field packets from independent ring clicks."""

import argparse
from hashlib import sha256
import json
from math import ceil, log, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import measure as field

HERE = Path(__file__).resolve().parent
DELAYS = (0, 3 / 16, 1 / 2, 1)
SIGNS = (-1, 1)
ARMS = (
    "minus_minus",
    "minus_plus",
    "plus_minus",
    "plus_plus",
    "vacuum_minus",
    "vacuum_plus",
    "source_alone",
    "probe_alone",
)
SHOTS = 2**27
SEED = 612283


def fly(h, states, dwell):
    return expm_multiply(-1j * dwell * h, states) if dwell else states


def laws():
    previous = json.loads((HERE.parent / "results.json").read_text())
    select = previous["selected"]
    assert select["head"] == [4, 11] and select["ray"] == [10, 11]
    assert (
        previous["field_dwell_exposure"] == 0.5 and previous["ray_exposure"] == 3 / 16
    )
    patch = field.search.Patch()
    words = np.array(patch.field_component())
    assert len(words) == previous["field_words"] == 458
    vacuum = np.ones(len(words)) / sqrt(len(words))
    h = field.ring_laplacian(patch, words)
    rings = field.components(patch, words)

    def head(edge):
        return np.where(words & patch.bit[tuple(edge)], -1, 1)

    z, w = head(select["head"]), head(select["ray"])
    v = fly(h, z * vacuum, 0.5)
    c, s = np.cos(field.THETA), np.sin(field.THETA)
    assert abs(2 * np.sin(2 * field.THETA) ** 2 - 1) < 1e-14
    states = np.column_stack(
        [
            (c + 1j * t * s * w) * (c * vacuum + 1j * sig * s * v)
            for sig in SIGNS
            for t in SIGNS
        ]
        + [(c + 1j * t * s * w) * vacuum for t in SIGNS]
        + [v, w * vacuum]
    )
    assert np.allclose(np.linalg.norm(states, axis=0), 1)
    probabilities, predicted = [], []
    for dwell in DELAYS:
        current = fly(h, states, dwell)
        projected = np.asarray([op @ current for op in rings])
        ps = np.sum(abs(projected) ** 2, axis=1)
        cross = np.einsum(
            "ij,ij->i", projected[:, :, 6].conj(), projected[:, :, 7]
        ).real
        np.testing.assert_allclose(
            ps[:, 0] - ps[:, 1] - ps[:, 2] + ps[:, 3], cross, atol=1e-13, rtol=0
        )
        np.testing.assert_allclose(ps[:, 4], ps[:, 5], atol=1e-13, rtol=0)
        probabilities.append(ps)
        predicted.append(cross)
    return patch.faces, np.array(probabilities), np.array(predicted)


def infer(faces, counts, predicted):
    # Simultaneous raw-bin coverage, rounded upwards to one part in a billion.
    eps = ceil(1e9 * sqrt(log(2 * counts.size / 0.01) / (2 * SHOTS))) / 1e9
    rates = counts / SHOTS
    contrast = rates[:, :, 0] - rates[:, :, 1] - rates[:, :, 2] + rates[:, :, 3]
    vacuum = rates[:, :, 5] - rates[:, :, 4]
    assert np.max(abs(vacuum)) <= 2 * eps
    assert np.max(abs(contrast - predicted)) <= 4 * eps
    upper = np.sqrt((rates[:, :, 6] + eps) * (rates[:, :, 7] + eps))
    assert np.all(abs(contrast) - 4 * eps <= upper)
    probe = {10, 11}
    remote = np.array([not probe.issubset(set(face)) for face in faces])
    separated = np.array([probe.isdisjoint(face) for face in faces])
    stats = []
    for j, d in enumerate(DELAYS):
        clear = abs(contrast[j]) > 4 * eps
        k = int(np.argmax(np.where(separated, abs(contrast[j]), -1)))
        value = float(contrast[j, k])
        stats.append(
            dict(
                delay=d,
                resolved_unaddressed_faces=int(sum(clear & remote)),
                resolved_vertex_disjoint_faces=int(sum(clear & separated)),
                strongest_vertex_disjoint_face=list(faces[k]),
                strongest_contrast=value,
                simultaneous_interval=[value - 4 * eps, value + 4 * eps],
                worst_field_forecast_error=float(max(abs(contrast[j] - predicted[j]))),
                largest_vacuum_contrast=float(max(abs(vacuum[j]))),
            )
        )
    return dict(
        total_banks=int(counts.size),
        shots_per_bank=SHOTS,
        raw_probability_radius=eps,
        contrast_radius=4 * eps,
        faces_without_probe_head=int(sum(remote)),
        faces_without_probe_vertex=int(sum(separated)),
        by_delay=stats,
    )


def acquire(write):
    faces, probabilities, predicted = laws()
    archive = HERE / "counts.json"
    if write:
        counts = np.random.default_rng(SEED).binomial(
            SHOTS, np.clip(probabilities, 0, 1)
        )
        data = dict(
            faces=[list(f) for f in faces],
            delays=DELAYS,
            arms=ARMS,
            trials_per_setting=SHOTS,
            seed=SEED,
            clicks=counts.tolist(),
            protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
            charged_comparison_sha256=sha256(
                (HERE.parent / "results.json").read_bytes()
            ).hexdigest(),
        )
        archive.write_text(json.dumps(data, indent=1) + "\n")
    saved = json.loads(archive.read_text())
    assert (
        saved["protocol_sha256"]
        == sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert (
        saved["charged_comparison_sha256"]
        == sha256((HERE.parent / "results.json").read_bytes()).hexdigest()
    )
    assert saved["faces"] == [list(f) for f in faces] and saved["arms"] == list(ARMS)
    assert saved["delays"] == list(DELAYS) and saved["trials_per_setting"] == SHOTS
    counts = np.asarray(saved["clicks"], dtype=np.int64)
    assert counts.shape == probabilities.shape
    summary = infer(faces, counts, predicted)
    assert (
        np.max(abs(counts / SHOTS - probabilities)) <= summary["raw_probability_radius"]
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--replay", action="store_true", help="verify the deposited counts"
    )
    args = parser.parse_args()
    print(json.dumps(acquire(args.write), indent=2))
