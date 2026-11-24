"""Fixed-marginal perturbations preserving both material rulers' inferred volumes."""

import argparse
import json

import numpy as np

from .measure import HERE, NAMES, geometry


def cycle_basis():
    contrasts = np.array([[1, 0, -1], [0, 1, -1]], float).T
    axes, _ = np.linalg.qr(contrasts)
    return np.array([np.outer(a, b) for a in axes.T for b in axes.T])


def gradients(laws, ref, delta):
    read = [geometry(laws[:, :, :, j, :], ref, delta) for j in range(2)]
    return np.array([item[key] for key in ("dlog_volume", "dlog_area") for item in read])


def null_transfer(operator, area, basis, ref):
    rows = operator[[0, 1, 2 + area]]
    rows = rows / np.linalg.norm(rows, axis=1)[:, None]
    _, singular, right = np.linalg.svd(rows, full_matrices=True)
    if singular[-1] < 1e-8:
        raise ValueError("the three reconstructed constraints are numerically dependent")
    delta = np.einsum("a,aij->ij", right[-1], basis)
    if -np.sum(delta * np.log(ref)) < 0:
        delta = -delta
    return delta, singular


def entropy_maximum_check(laws, basis):
    uniform = np.ones((3, 3)) / 9
    operator = np.stack([gradients(laws, uniform, d) for d in basis], axis=1)
    volume = operator[:2]
    volume = volume / np.linalg.norm(volume, axis=1)[:, None]
    _, singular, right = np.linalg.svd(volume, full_matrices=True)
    if singular[-1] < 1e-8:
        raise ValueError("the two volume constraints are numerically dependent")
    kernel = right[2:]
    rows = []
    for area, name in enumerate(NAMES):
        projected = (operator[2 + area] @ kernel.T) @ kernel
        norm = np.linalg.norm(projected)
        if norm < 1e-14:
            rows.append(dict(area_ruler=name, projected_gradient_norm=norm))
            continue
        delta = np.einsum("a,aij->ij", projected / norm, basis)
        rows.append(dict(area_ruler=name, projected_gradient_norm=norm,
                         transfer=delta.tolist(),
                         derivatives=gradients(laws, uniform, delta).tolist()))
    return rows


def run():
    bank = json.loads((HERE / "results.json").read_text())
    ref = np.array(bank["reference"])
    basis = cycle_basis()
    transfers = []
    boards = []
    for board in bank["screens"]:
        laws = np.array(board["word_moments"])
        operator = np.stack([gradients(laws, ref, d) for d in basis], axis=1)
        rows = []
        for area, name in enumerate(NAMES):
            delta, singular = null_transfer(operator, area, basis, ref)
            if len(transfers) <= area:
                transfers.append(delta)
            held = transfers[area]
            slope = -np.sum(held * np.log(ref))
            for sign in (-1, 1):
                assert np.min(ref + sign * held / 32) > 0
            rows.append(dict(area_ruler=name, gradient_singular_values=singular.tolist(),
                             local_entropy_slope=-np.sum(delta * np.log(ref)),
                             local_derivatives=gradients(laws, ref, delta).tolist(),
                             transferred_entropy_slope=slope,
                             transferred_derivatives=gradients(laws, ref, held).tolist()))
        boards.append(dict(cells=board["cells"], records=rows,
                           uniform_reference=entropy_maximum_check(laws, basis)))
    return dict(derivative_order=["first_return_log_volume", "free_pulse_log_volume",
                                  "first_return_log_area", "free_pulse_log_area"],
                twenty_cell_transfers=np.array(transfers).tolist(), boards=boards)


if __name__ == "__main__":
    from ..collar_pairing.measure import compare

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    path = HERE / "both_volumes.json"
    if args.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        compare(output, json.loads(path.read_text()))
        print("All fixed-volume derivatives agree with the receipt.")
