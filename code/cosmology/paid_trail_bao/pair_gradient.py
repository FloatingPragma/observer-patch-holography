"""Price the observed squared field differences on the original clean material chart."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    pair_field as field,
)
from . import closed_inventory as inventory

HERE = inventory.HERE


def run():
    pair = field.pair
    _, w, cells, _, codes, h, ids, groups = pair.setup()
    board, home, marks, *_ = field.lab.energy.returned(w)
    clean = set(range(board.n)) - set(home) - set(np.ravel(marks))
    ports = np.array([e for e in board.seams if set(e) <= clean])
    # Each row is reconstructed from the four singly placed ordinary donors.
    # Divide the signed field difference by -16 before forming its square.
    adjacent = board.adjacency().toarray()[:, home].astype(int)
    gradients = adjacent[ports[:, 0]] - adjacent[ports[:, 1]]
    gram = gradients.T @ gradients
    np.testing.assert_array_equal(gram, np.diag(np.diag(gram)))
    assert all(gram.diagonal() > 0)
    homogeneous = ((w[:, :4] @ gradients.T) ** 2).sum(1)
    weighted = (w[:, :4] ** 2).sum(1)
    np.testing.assert_array_equal(homogeneous, (w[:, :4] ** 2) @ gram.diagonal())
    # Fixed reference and sharp stationary reweightings matching the head
    # marginals, measured source displacement, compact energies and mean loads.
    rows = np.vstack(
        (pair.constraints(w, cells, codes, h, "energy"), w[:, :4].T, (w[:, 1] - 1) ** 2)
    )
    A = pair.group_matrix(rows, ids, groups)
    b = A.sum(1)
    b[-1] = 80
    z = field.lab.energy.least_change(cells)
    result, ranges = {}, {}
    for key, score, calibration in (
        ("same_link_units", homogeneous, 4 * gram.diagonal()),
        ("each_home_calibrated", weighted, [4]),
    ):
        c = pair.group_matrix(score[None, :], ids, groups)[0]
        canonical = (z - 1) @ score / len(w)
        extreme = []
        witnesses = []
        for sign in (1, -1):
            opt = linprog(sign * c, A_eq=A, b_eq=b, bounds=(0, None))
            assert opt.success
            primal = np.array(pair.rational(opt.x), object)
            dual = np.array(pair.rational(opt.eqlin.marginals), object)
            assert all(primal >= 0) and np.all(A @ primal == b)
            assert np.all(dual @ A <= sign * c) and dual @ b == sign * (c @ primal)
            extreme.append((c @ primal - int(score.sum())) / len(w))
            witnesses.append(
                dict(primal=list(map(str, primal)), dual=list(map(str, dual)))
            )
        result[key] = dict(
            canonical=str(canonical),
            extremes=list(map(str, extreme)),
            witnesses=witnesses,
        )
        for tag, ends in (("repair", [canonical]), ("all_matched", extreme)):
            vals = [v / F(int(u)) for v in ends for u in calibration]
            ranges[key + "_" + tag] = min(vals), max(vals)
    src = HERE / "closed_inventory_maps.json"
    maps = json.loads(src.read_bytes())
    inputs, *_ = inventory.rp.galaxy.transfer.budget.load_inputs()
    forecasts = {}
    for tag, (lo, hi) in ranges.items():
        fits = inventory.mass.forecast(maps, [float(lo), float(hi)], inputs)
        settings = {}
        for tracer, cap in fits.items():
            for k, tab in cap.items():
                row = tab["closed"]
                settings[tracer + "/" + k] = dict(
                    within_two=row["within_two"],
                    complete=sum(
                        v == len(inputs[n]["radius_kpc"])
                        for n, v in row["galaxies"].items()
                    ),
                )
        forecasts[tag] = dict(
            per_accepted_patch_in_double_units=[str(lo), str(hi)],
            settings=settings,
            maximum_positive_population_mass_ratio=float(
                max(0, hi) * F(2, 9) * F(448, 2110)
            ),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                field.HERE / "pair_field.json",
                src,
                HERE / "record_profile_maps.json",
            )
        },
        clean_edges=len(ports),
        edges_with_donor_contrast=int(np.any(gradients, axis=1).sum()),
        donor_gram=gram.tolist(),
        native_scores=result,
        galaxy_forecasts=forecasts,
        galaxies=len(inputs),
        rings=sum(len(v["radius_kpc"]) for v in inputs.values()),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pair_gradient.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(dst.read_bytes()), result)
    for tag, v in result["galaxy_forecasts"].items():
        print(
            tag,
            v["per_accepted_patch_in_double_units"],
            v["maximum_positive_population_mass_ratio"],
            max(r["within_two"] for r in v["settings"].values()),
        )
