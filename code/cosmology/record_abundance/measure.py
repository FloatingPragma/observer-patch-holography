"""Test how much homogeneous signed-record charge individual source measurements fix."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from matched_einstein_limit.collar_identifiability.refinement import (
    admission,
    transport,
)
from matched_einstein_limit.collar_identifiability.refinement.composition import (
    measure as array,
)
from model.executable_model import federation
from model.executable_model.refinement import load_refinement

HERE = Path(__file__).parent
RECEIPT = HERE / "results.json"
WIDTH = 25
SURVIVAL = 7 / 8
UNIFORM_PLACEMENT = 1 / 128
PAIRING = 1 / 2
RADII = (1, 2, 4)


def placements():
    c = np.cos(2 * np.pi * np.fft.fftfreq(WIDTH))
    neighbour = c[:, None, None] + c[None, :, None] + c[None, None, :]
    mark = np.zeros((WIDTH,) * 3)
    mark[(WIDTH // 2,) * 3] = 1 / 4
    visits = np.fft.ifftn(np.fft.fftn(mark) / (1 - SURVIVAL * (1 + neighbour) / 4)).real
    np.testing.assert_allclose(
        visits - SURVIVAL * array.periodic_step(visits), mark, atol=1e-15, rtol=0
    )
    assert 0 < visits.min() < visits.max() < 1
    return visits, 2 * (3 - neighbour)


def policies(level):
    board = (
        load_refinement(HERE / "FEDERATION_LEVEL3.npz", level=level)
        if level == 3
        else federation(level)
    )
    owners = np.asarray(board.external) // 12
    a, b = owners.T
    graph = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(board.cells,) * 2
    )
    face_size = 4**level
    faces = np.arange(board.cells).reshape(-1, face_size)
    first = faces[0]
    second = next(face for face in faces[1:] if not graph[first, :][:, face].nnz)
    region = np.union1d(first, second)
    # The galaxy and the switched domains must have no common comparator seam.
    leave_out = np.zeros(board.cells, dtype=bool)
    leave_out[region] = True
    leave_out[graph[region, :].indices] = True
    galaxy = int(np.flatnonzero(~leave_out)[0])

    def greater_color(face):
        local = graph[face, :][:, face]
        distances = shortest_path(local, directed=False, unweighted=True, indices=0)
        assert np.isfinite(distances).all()
        group = distances.astype(int) % 2 == 0
        row, col = local.nonzero()
        assert np.all(group[row] != group[col])
        return group if group.sum() >= len(group) / 2 else ~group

    staggered = np.r_[first[greater_color(first)], second[~greater_color(second)]]
    complements = {}
    for label, selected in (("whole_face", first), ("staggered", staggered)):
        choices = np.stack(
            [
                np.isin(np.arange(board.cells), group)
                for group in (selected, np.setdiff1d(region, selected))
            ]
        )
        np.testing.assert_array_equal(choices.sum(axis=1), [face_size, face_size])
        np.testing.assert_array_equal(
            choices.sum(axis=0), np.isin(np.arange(board.cells), region)
        )
        complements[label] = choices
    return owners, region, galaxy, complements


def charge(rate, contrast):
    signal = rate * contrast / 3
    return signal * np.arctanh(PAIRING * signal)


def gravity(density, eigen):
    # The scale and conversion are fixed to one in this chart. Set the potential zero mode to zero.
    inverse = np.zeros_like(eigen)
    np.divide(1, eigen, out=inverse, where=eigen != 0)
    phi = np.fft.ifftn(
        -4 * np.pi * np.fft.fftn(density - density.mean()) * inverse
    ).real
    acceleration = np.stack(
        [(np.roll(phi, 1, axis) - np.roll(phi, -1, axis)) / 2 for axis in range(3)]
    )
    sheet = density.sum(axis=2)
    coords = np.arange(WIDTH) - WIDTH // 2
    squared = coords[:, None] ** 2 + coords[None, :] ** 2
    esd = np.array(
        [sheet[squared <= r * r].mean() - sheet[squared == r * r].mean() for r in RADII]
    )
    c = WIDTH // 2
    force = [-acceleration[0, c + r, c, c] for r in RADII]
    return acceleration, esd, force


def explicit_probability_slope(owners, mask, galaxy, visits, fresh, formed):
    placements = UNIFORM_PLACEMENT * mask.astype(float)
    placements[galaxy] = visits
    left, right = (placements[owners] / 3).T
    p = array.edge_probabilities(left, right, PAIRING, fresh, formed)
    dp = array.edge_probabilities(
        left, right, 1, fresh, formed
    ) - array.edge_probabilities(left, right, 0, fresh, formed)
    return float(dp @ (np.log(p) - np.log1p(-p)))


def panel(level, visits, eigen, fresh, formed):
    contrast = formed - fresh
    owners, region, galaxy, banks = policies(level)
    membership = np.isin(owners, region)
    outer = int((membership.sum(axis=1) == 1).sum())
    inner = int((membership.sum(axis=1) == 2).sum())
    background = charge(UNIFORM_PLACEMENT, contrast)
    galaxy_charge = 3 * charge(visits, contrast)
    rows, signals = {}, []
    c = WIDTH // 2
    for label, masks in banks.items():
        cuts = [int(np.sum(mask[owners[:, 0]] != mask[owners[:, 1]])) for mask in masks]
        density = galaxy_charge + np.mean(cuts) * background
        # Four complete log variations, including binary outcome reversals at negative seams.
        for r in (0, *RADII):
            direct = np.mean(
                [
                    explicit_probability_slope(
                        owners, mask, galaxy, visits[c + r, c, c], fresh, formed
                    )
                    for mask in masks
                ]
            )
            np.testing.assert_allclose(direct, density[c + r, c, c], atol=5e-17, rtol=0)
        field, esd, inward = gravity(density, eigen)
        assert np.all(esd > 0) and np.all(np.array(inward) > 0)
        signals.append((field, esd))
        rows[label] = {
            "cut_counts_for_two_words": cuts,
            "active_cells_each_word": 4**level,
            "mean_background_charge": float(np.mean(cuts) * background),
            "mean_total_charge": float(density.mean()),
            "inward_peculiar_acceleration": list(map(float, inward)),
            "surface_density_contrast": esd.tolist(),
            "mean_number_of_disagreeing_internal_neighbours": float(
                np.mean(cuts) - outer / 2
            ),
        }
    np.testing.assert_allclose(
        np.mean(rows["whole_face"]["cut_counts_for_two_words"]), outer / 2
    )
    np.testing.assert_allclose(
        np.mean(rows["staggered"]["cut_counts_for_two_words"]), outer / 2 + inner
    )
    ratio = (
        rows["staggered"]["mean_background_charge"]
        / rows["whole_face"]["mean_background_charge"]
    )
    np.testing.assert_allclose(ratio, 2**level)
    residuals = [float(np.max(np.abs(a - b))) for a, b in zip(*signals, strict=True)]
    assert max(residuals) < 1e-17
    return {
        "board_cells": int(20 * 4**level),
        "subdivision_length": 2**level,
        "switched_region_cells": region.tolist(),
        "galaxy_cell": galaxy,
        "outside_boundary_edges": outer,
        "interior_edges": inner,
        "policies": rows,
        "background_ratio": float(ratio),
        "full_mean_ratio": rows["staggered"]["mean_total_charge"]
        / rows["whole_face"]["mean_total_charge"],
        "maximum_peculiar_acceleration_difference": residuals[0],
        "maximum_enrolled_shear_difference": residuals[1],
    }


def run():
    transport.native_fair_bit()
    reference = admission.material_panel(1, 1, 8)["records"][-1]
    fresh, contrast = reference["q"][-1], reference["d"][-1]
    visits, eigen = placements()
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "chart_side": WIDTH,
        "survival": SURVIVAL,
        "uniform_injection_per_site_per_active_cell": (1 - SURVIVAL)
        * UNIFORM_PLACEMENT,
        "galaxy_injection_at_one_site": 0.25,
        "reference_pairing": PAIRING,
        "native_return_contrast": contrast,
        "chart_radii": list(RADII),
        "panels": [panel(i, visits, eigen, fresh, fresh + contrast) for i in (1, 2, 3)],
    }


if __name__ == "__main__":
    from cosmology.rar_deep_regime.lensing_profile import check

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RECEIPT.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        expected = json.loads(RECEIPT.read_text())
        for row in (*expected["panels"], *result["panels"]):
            for field in (
                "maximum_peculiar_acceleration_difference",
                "maximum_enrolled_shear_difference",
            ):
                if not 0 <= row.pop(field) < 1e-17:
                    raise ValueError(
                        "Model residual exceeds the shift-invariance tolerance"
                    )
        check(expected, result)
    for row in result["panels"]:
        print(row["board_cells"], row["background_ratio"], row["full_mean_ratio"])
