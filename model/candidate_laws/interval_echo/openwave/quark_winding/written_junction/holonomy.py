"""Peierls flux on the matter hops of a written ring, with two flux-sensitive controls."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import currents
from . import backgrounds, ring


def twist(h, loads, potential):
    """Assign phases from transported matter alone. Field-only flips keep their amplitude."""
    h = h.tocoo()
    values = h.data.astype(complex)
    for k, (r, c) in enumerate(zip(h.row, h.col)):
        diff = loads[r] - loads[c]
        if np.any(diff):
            a = int(np.flatnonzero(diff == -1)[0])
            b = int(np.flatnonzero(diff == 1)[0])
            if (a + 1) % len(potential) == b:
                values[k] *= np.exp(1j * potential[a])
            else:
                assert (b + 1) % len(potential) == a
                values[k] *= np.exp(-1j * potential[b])
    return coo_matrix((values, (h.row, h.col)), shape=h.shape).tocsr()


def ring_panel(row, geometry, extension, angle):
    length = row["length"]
    data = backgrounds.block(
        backgrounds.bits(row["gauss"], length), np.array(geometry["reference"])
    )
    h = backgrounds.generator(data, extension)
    a = np.zeros(length)
    path = geometry["ports"]
    seam = next(
        i for i in range(length) if path[i] // 12 != path[(i + 1) % length] // 12
    )
    a[seam] = angle
    new = twist(h, data[2], a)
    phase = np.exp(-1j * (data[1] @ a))
    expected = h.multiply(phase[:, None] * phase.conj()[None, :]).tocsr()
    residual = abs(new - expected).max()
    basis = np.eye(1, h.shape[0], 0).ravel()
    p = abs(expm_multiply(-0.5j * h, basis)) ** 2
    q = abs(expm_multiply(-0.5j * new, basis)) ** 2
    return dict(
        length=length,
        gauss=row["gauss"],
        seam=[path[seam], path[(seam + 1) % length]],
        twist=angle,
        conjugacy_residual=float(residual),
        port_probability_residual=max(abs(p - q)),
    )


def controls(angle):
    # Head-independent single-particle ring on the same twelve cyclic coordinates.
    length = 12
    loads = np.eye(length, dtype=int)
    h = coo_matrix(-4 * (np.roll(loads, 1, axis=0) + np.roll(loads, -1, axis=0)))
    a = np.zeros(length)
    a[0] = angle
    b = twist(h, loads, a)
    ea, eb = [np.linalg.eigvalsh(z.toarray()) for z in (h, b)]
    ps = [abs(expm_multiply(-1.5j * z.tocsr(), loads[0])) ** 2 for z in (h, b)]
    # A real active face. The extra can make a full tour, then a field-only
    # flip restores all three heads. Matter and field acquire different phases.
    w, _, n = ring.basis(np.array([1, 0, 0]))
    face = ring.generator(w, n, np.zeros(3), True).tolil()
    i, j = 0, len(w) - 1
    assert w[i] == 0 and w[j] == 7
    face[i, j] = face[j, i] = -0.5
    face[i, i] += 0.5
    face[j, j] += 0.5
    bt = twist(face, n, np.array([angle, 0, 0]))

    def moment(z):
        v = z.tocsr().getcol(j)
        v = z @ v
        return float(abs(v).power(2).sum())

    fourth = moment(bt) - moment(face)
    return dict(
        twist=angle,
        ordinary_spectrum_shift=max(abs(ea - eb)),
        ordinary_port_probability_shift=max(abs(ps[0] - ps[1])),
        active_face_fourth_moment_difference=fourth,
        single_triangle_loop_prediction=64 * (np.cos(angle) - 1),
    )


def run(scans):
    angles = [np.pi / 2, 2 * np.pi / 3]
    return currents._portable(
        dict(
            written={
                law: [
                    ring_panel(
                        r, scans["geometry"][str(r["length"])], law == "extension", a
                    )
                    for r in pair["selected"]
                    for a in angles
                ]
                for law, pair in scans["paired"].items()
            },
            controls=[controls(a) for a in angles],
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    with args.output.open("x") as f:
        f.write(json.dumps(run(s), indent=2, allow_nan=False) + "\n")
