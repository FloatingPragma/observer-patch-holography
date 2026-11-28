"""Read the archived two-charge outcomes without importing their producer."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = Path(__file__).resolve().parent


def decode(z, times):
    angle = np.angle(z)
    # Innovations remove the integer windings between consecutive dyadic times.
    diff = np.angle(np.exp(1j * (angle[1:] - 2 * angle[:-1])))
    return angle[0] / times[0] + np.sum(diff / times[1:])


def test_all_pair_frequencies_from_raw_joint_outcomes():
    data = json.loads((HERE / "copies.json").read_text())
    path = HERE / "copies_counts.npz"
    assert sha256(path.read_bytes()).hexdigest() == data["counts_sha256"]
    assert not set(data["stock"]["routes"][0]) & set(data["stock"]["routes"][1])
    dt = np.array(data["dyadic_times"])
    S = data["shots_per_joint_bank"]
    errs = []
    with np.load(path) as raw:
        for row in data["rows"]:
            i = row["index"]
            count = raw[f"pair_{i}"]
            assert (count.sum(-1) == S).all()
            same = (count[..., 0] + count[..., 3]) / S
            # Equal outcomes have parity +; free phases add across the two axes.
            z = 2 * (same[0, 0] - same[1, 1]) + 2j * (same[0, 1] + same[1, 0] - 1)
            value = decode(z, dt) / 2
            np.testing.assert_allclose(value, row["energy"], atol=1e-13)
            errs.append(value - np.sqrt(3.5**2 + 8 * np.sin(row["q"] / 2) ** 2))
            if not i:
                continue
            plus, minus = raw[f"space_{i}"]
            edges = (plus - minus)[2:-4:2] / float(sum(plus))
            u = edges[:-2] + edges[2:]
            v = 2 * edges[1:-1]
            # With s_n a sinusoid this recurrence gives cos(q) locally.
            projected = np.clip(v @ u / (v @ v), -1, 1)
            assert round(np.arccos(projected) * len(data["rows"]) / np.pi) == i
    assert max(abs(np.array(errs))) < data["paired_frequency_radius"]


def test_prior_selector_over_every_interior_source_and_end_controls():
    data = json.loads((HERE / "copy_source.json").read_text())
    n = 84
    diag = np.full(n, 4.5)
    diag[[0, -1]] = 8.5
    hops = np.tile([-0.5, -4.0], n // 2)[:-1]
    vals, vec = eigh_tridiagonal(diag, hops)
    np.testing.assert_allclose(vals[np.argmin(abs(vals - 4.5))], 8, atol=1e-13, rtol=0)
    frequencies = vals[:, None] + vals
    mask = abs(frequencies - 9) < 1e-10
    assert np.sum(mask) == 82
    weights = (vec**2) @ mask.astype(float) @ (vec**2).T
    inside = weights[1:-1, 1:-1]
    rec = data["paired_sharp"][0]
    np.testing.assert_allclose(
        [inside.min(), inside.mean(), inside.max()],
        [rec[k] for k in ("min_input_weight", "mean_input_weight", "max_input_weight")],
    )
    assert inside.min() > 0 and inside.size == 6724
    assert [r["selected"] for r in data["paired_sharp"]] == [9.0, 16.0, 16.0, 17.0]
    assert [r["pass_both"] for r in data["filters"][0]["records"]] == [
        False,
        False,
        True,
        True,
    ]
    assert [r["pass_both"] for r in data["filters"][1]["records"]] == [
        False,
        False,
        True,
        True,
    ]
    for src in data["filters"]:
        for r in src["records"]:
            for t in r["phase_reads"]:
                assert np.all(np.sum(t["counts"], axis=1) == 2**22)
    scale = data["fixed_law"]
    assert scale["bulk_compton_length_cells"] < 0.41
    assert scale["continuum_energy"] - scale["pair_energy"] > 0.54
    assert abs(scale["measured_lattice_energy"] - scale["pair_energy"]) < 1e-5
