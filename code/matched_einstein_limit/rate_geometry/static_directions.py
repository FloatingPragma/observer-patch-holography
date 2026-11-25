"""Angular apertures of a native cell and a separate static curvature control."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, null_space

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
    compile_wave,
)
from model.executable_model import cell
from model.executable_model.carrier import icosahedron_from_coordinates, isomorphism


def design():
    board = cell()
    directions, edges = icosahedron_from_coordinates()
    mapping = isomorphism(board.seams, edges)
    assert mapping is not None
    directions = directions[mapping] / np.linalg.norm(directions[0])
    dots = directions @ directions.T
    opposite = dots.argmin(axis=1)
    poles = np.flatnonzero(np.arange(12) < opposite)
    axis = np.maximum(np.arange(12), opposite)
    stars = np.eye(12, dtype=int) + board.adjacency().toarray().astype(int)
    for row in stars:
        np.testing.assert_array_equal(np.sort(axis[row > 0]), np.sort(axis[poles]))
    raw = stars[poles]
    weighted = raw + np.eye(12, dtype=int)[poles]
    # Collapse antipodal pairs for an even response.
    expand = np.eye(12, dtype=int)[:, poles] + np.eye(12, dtype=int)[:, opposite[poles]]
    np.testing.assert_array_equal(raw @ expand, np.ones((6, 6), int))
    np.testing.assert_array_equal(
        weighted @ expand, np.ones((6, 6), int) + np.eye(6, dtype=int)
    )
    return directions, poles, raw, weighted


def native_reads(directions, raw, weighted):
    board = cell()
    pole = directions[0]
    wavevectors = np.vstack([pole, null_space(pole[None, :]).T])
    waves = np.column_stack(
        [compile_wave(board, directions @ v / 2) for v in wavevectors]
    )
    h = candidate.generator(board).toarray()
    free = expm(-1.75j * h) @ waves
    probabilities = abs(free) ** 2
    np.testing.assert_allclose(probabilities, abs(waves) ** 2, atol=2e-13)
    laws = [probabilities[:, 0], probabilities[:, 1:].mean(axis=1)]
    rng = np.random.default_rng(76153)
    record = []
    axial = np.ones(12, dtype=int)
    axial[[0, (directions @ pole).argmin()]] = 5
    transverse = (axial == 1).astype(int)
    for name, law, weights in zip(("axial", "transverse_mixture"), laws, (axial, transverse), strict=True):
        # The two poles have squared axial direction cosine 1; the others 1/5.
        # The mean over two orthogonal transverse polarizations is basis independent.
        canonical = weights / weights.sum()
        np.testing.assert_allclose(law, canonical, atol=2e-13, rtol=0)
        counts = rng.multinomial(2**20, canonical)
        np.testing.assert_allclose(raw @ law, 0.5, atol=3e-13)
        # The score is one at the pole, one half at each rim port, zero outside.
        record.append(
            dict(
                preparation=name,
                counts_by_port=counts.tolist(),
                shots=int(counts.sum()),
                equal_apertures=(raw @ counts / counts.sum()).tolist(),
                pole_weighted_apertures=(
                    weighted @ counts / (2 * counts.sum())
                ).tolist(),
                reference_pole_exact=float((weighted @ law)[0] / 2),
            )
        )
    return record


def comparison(directions, poles, weighted, sectional, lapse_hessian):
    # On Sigma^2_K x R the two radial normal Jacobi sectional curvatures
    # are 0 and K*(1-n_z^2). The even central lapse is 1+c*z^2/2.
    z2 = directions[:, 2] ** 2
    sectional_normal = sectional * (1 - z2)
    target = (sectional + lapse_hessian) * (1 - z2[poles])
    rows = []
    for radius in (0.5, 0.25, 0.125, 0.0625):
        x = radius * np.sqrt(abs(sectional_normal))
        factor = np.ones(12)
        use = x > 0
        fn = np.sin if sectional > 0 else np.sinh
        factor[use] = fn(x[use]) / x[use]
        lapse = 1 + lapse_hessian * radius**2 * z2 / 2
        # factor gives area per calibrated solid angle divided by radius squared.
        area = -6 * (weighted @ factor / 7 - 1) / radius**2
        clock = 2 * (weighted @ lapse / 7 - 1) / radius**2
        ric = area + 3 * clock.mean() - clock
        ray = 7 * ric - 6 * ric.mean()
        rows.append(
            dict(
                radius=radius,
                reconstructed=ray.tolist(),
                max_absolute_error=float(max(abs(ray - target))),
            )
        )
    return dict(
        sectional_curvature=sectional,
        lapse_hessian_zz=lapse_hessian,
        independent_null_ricci=target.tolist(),
        radii=rows,
    )


def report():
    directions, poles, raw, weighted = design()
    return dict(
        axis_ports=poles.tolist(),
        ordinary_star_even_rank=1,
        weighted_star_even_rank=6,
        native_single_cell=native_reads(directions, raw, weighted),
        static_product_comparisons=[
            comparison(directions, poles, weighted, k, c)
            for k, c in ((1.0, -2.0), (-1.0, 2.0), (0.0, 0.0))
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result, end="")
