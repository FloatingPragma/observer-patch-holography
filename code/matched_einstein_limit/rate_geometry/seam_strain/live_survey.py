"""Five-station intrinsic ball and boundary surveys on the retained seam source."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import live_radar, live_source, live_sphere
from .source import measure as source

HERE = Path(__file__).parent
FACES = tuple(combinations(range(4), 3))


def stations():
    board = source.formation.federation(1)
    outer = source.formation.region(board)
    hop = shortest_path(
        board.adjacency(), indices=outer, directed=False, unweighted=True
    )
    candidates = [
        (int(hop[:, p].sum()), int(hop[:, p].max()), p)
        for p in range(board.n)
        if p not in outer
    ]
    best = min(candidates)[:2]
    return tuple(outer) + tuple(
        p for total, worst, p in candidates if (total, worst) == best
    )


def square(distances, n):
    matrix = np.zeros((n, n), dtype=distances.dtype)
    a, b = np.triu_indices(n, 1)
    matrix[a, b] = matrix[b, a] = distances
    return matrix


def star(d):
    """Exterior area, sum of tetrahedron volumes, and the apex solid angle."""
    radii = d[4, :4]
    gram = (radii[:, None] ** 2 + radii[None, :] ** 2 - d[:4, :4] ** 2) / 2
    areas = []
    volumes = []
    angles = []
    for i, j, k in FACES:
        a, b, c = d[i, j], d[i, k], d[j, k]
        area = np.sqrt(4 * a * a * b * b - (a * a + b * b - c * c) ** 2) / 4
        three = gram[np.ix_([i, j, k], [i, j, k])]
        determinant = np.sqrt(np.linalg.det(three))
        denominator = (
            radii[i] * radii[j] * radii[k]
            + radii[i] * gram[j, k]
            + radii[j] * gram[k, i]
            + radii[k] * gram[i, j]
        )
        areas.append(area)
        volumes.append(determinant / 6)
        angles.append(
            4
            * np.arctan(
                determinant / (np.sqrt(denominator**2 + determinant**2) + denominator)
            )
        )
    return np.array([sum(areas), sum(volumes), sum(angles)])


def admission(d):
    r = d[4, :4]
    gram = (r[:, None] ** 2 + r[None, :] ** 2 - d[:4, :4] ** 2) / 2
    eigen = np.linalg.eigvalsh(gram)
    if eigen[0] <= 0:
        raise ValueError(
            "the five input distances do not form a nondegenerate 4-simplex"
        )
    heights = []
    for ids in FACES:
        face = gram[np.ix_(ids, ids)]
        # Distance to affine plane through three points of the face.
        heights.append(1 / np.sqrt(np.linalg.solve(face, np.ones(3)).sum()))
    return dict(gram_eigenvalues=eigen.tolist(), safe_cone_radius=min(heights))


def slopes(d, change):
    return star(d + 1e-15j * change).imag / 1e-15


def rotations(rho):
    result = {}
    for side, marginal in enumerate(source.formation.marginals(rho)):
        traceless = marginal - np.eye(2) / 2
        generator = traceless / np.ptp(np.linalg.eigvalsh(traceless))
        op = (
            np.kron(generator, np.eye(2))
            if side == 0
            else np.kron(np.eye(2), generator)
        )
        result[("left", "right")[side]] = op
    return result


def modular_response(ref, tangents):
    # The uncoupled strain remains in a pure state. Its null space is a
    # physically absent sector shared by every incoming preparation.
    values, axes = np.linalg.eigh(ref)
    keep = values > 1e-12
    basis = axes[:, keep]
    log = -(basis * np.log(values[keep])) @ basis.conj().T
    out = {}
    for key, t in tangents.items():
        np.testing.assert_allclose(
            basis @ (basis.conj().T @ t @ basis) @ basis.conj().T, t, atol=1e-11
        )
        out[key] = float(np.trace(log @ t).real)
    return out


def retained(coupled):
    arm = live_source.postpump(coupled)
    ref = arm["reference"]
    changes = {key: arm[key] - ref for key in ("ordinary", "decorrelated")}
    changes["hidden"] = arm["hidden"]
    board = source.formation.federation(1)
    heads, _ = source.form(board, source.formation.region(board), -1)
    for name, generator in rotations(heads).items():
        variation = -1j * (generator @ heads - heads @ generator)
        changes[name] = live_source.propagate_heads(variation, coupled)
    k = live_source.source_calibration()
    charges = {
        "incoming": {
            tag: k[tag + "_modular"] for tag in ("ordinary", "hidden", "decorrelated")
        }
        | {"left": 0.0, "right": 0.0},
        "retained": modular_response(ref, changes),
    }
    return ref, changes, charges


def read(effect, state):
    return np.einsum("eij,ji->e", effect, state).real


def panel(ref, changes, charges, effects, ports, deadline, tag, coupled):
    n = len(ports)
    distances = square(read(effects, ref), n)
    dr = {k: square(read(effects, v), n) for k, v in changes.items()}
    centers = []
    for idx in range(4, n):
        sl = np.ix_([0, 1, 2, 3, idx], [0, 1, 2, 3, idx])
        d = distances[sl]
        admission_record = admission(d)
        area, volume, angle = star(d)
        values = {k: slopes(d, v[sl]) for k, v in dr.items()}
        matched = {}
        for adjuster in ("left", "right"):
            a, v, _ = values[adjuster]
            if abs(v) < 1e-12:
                raise ValueError("rotation cannot set the measured volume")
            group = {}
            for key in ("ordinary", "decorrelated", "hidden"):
                da, dv, _ = values[key]
                alpha = -dv / v
                group[key] = dict(
                    rotation=alpha,
                    area_slope=da + alpha * a,
                    modular={
                        s: c[key] + alpha * c[adjuster] for s, c in charges.items()
                    },
                )
            matched[adjuster] = group
        # The implicit radius is chosen from volume. At reference radius r,
        # dA|V = r^2 dW/3 holds identically throughout this star's cone.
        centers.append(
            dict(
                apex=ports[idx],
                **admission_record,
                outer_area=area,
                interior_volume=volume,
                solid_angle=angle,
                matched_boundary=matched,
                slopes={
                    k: dict(
                        area=v[0],
                        volume=v[1],
                        angle=v[2],
                        matched_cone_area_per_radius_squared=v[2] / 3,
                    )
                    for k, v in values.items()
                },
            )
        )
    smooth = live_sphere.sphere_fit(distances)
    smooth["curvature_slopes"] = {
        k: live_sphere.sphere_curvature_slope(distances, m, smooth["curvature"])
        for k, m in dr.items()
    }
    return dict(
        read=tag,
        coupled=coupled,
        deadline=deadline,
        source_slopes=charges,
        round_comparison=smooth,
        lengths=distances.tolist(),
        length_slopes={k: v.tolist() for k, v in dr.items()},
        centers=centers,
    )


def run():
    ports = stations()
    panels = []
    for coupled in (True, False):
        ref, changes, charges = retained(coupled)
        for deadline in (8, 10, 12):
            effects = [
                live_radar.range_effects(1, q, coupled, ports, deadline)
                for q in (-1, 1)
            ]
            for tag, e in zip(
                ("vacancy", "excess", "pooled"),
                [*effects, (effects[0] + effects[1]) / 2],
                strict=True,
            ):
                panels.append(
                    panel(ref, changes, charges, e, ports, deadline, tag, coupled)
                )
    return dict(stations=list(ports), panels=panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "live_survey.json").write_text(result)
    else:
        print(result, end="")
