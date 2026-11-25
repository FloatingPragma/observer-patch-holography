"""A stationary binary record at one selected reference, with an independent stability test."""

import json
from itertools import combinations

import mpmath as mp
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq

from .intrinsic_area import (
    geometry_slopes,
    head_source,
    pairing,
    radar_banks,
    ranging,
    survey_corners,
)


def bridge_corners():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    p, d = head_source()
    corners = []
    for odd in (False, True):
        h = ranging.generator(board, route, 1, 8, odd)
        (arm,) = pairing.material_read(h, route[1:3], {1: (p + d, p)})
        corners.append(
            [
                [r["complete"] for r in arm[prep]["appointments"]]
                for prep in ("independent", "formed")
            ]
        )
    return np.asarray(corners)


def gains(corners, click, u, s):
    area, volume, _, _, _, a, adjust = geometry_slopes(corners, u, s)
    c_read = (click[0, 1] - click[0, 0]) - (click[1, 1] - click[1, 0])
    direct = np.array(
        [
            (click[1] - click[0]) @ np.array([1 - s, s]),
            (click[:, 1] - click[:, 0]) @ np.array([1 - u, u]),
        ]
    )
    # A comparison with an independent, even fresh return and a fair tie bit.
    mean = np.array([1 - u, u]) @ click @ np.array([1 - s, s]) - click[0, 0]
    entropy = -np.arctanh(mean) * (direct + c_read * adjust)
    eta = -entropy[0] / a[0]
    return dict(
        ordinary_eta=float(eta),
        fraction=float(-eta * a[1] / entropy[1]),
        area=a.tolist(),
        entropy=entropy.tolist(),
        volume=volume,
        survey_area=area,
    )


def stationary_hessian(corners, click, u, s, eta):
    """Differentiate the area, volume and binary entropy in three affine mixture coordinates."""
    with mp.workdps(45):
        matrices = [
            [mp.matrix(corners[i, j].tolist()) for j in range(2)] for i in range(2)
        ]
        read = mp.matrix(click.tolist())

        def observables(u, s, v):
            # v is the probability of odd parity AND the formed pairing arm.
            weights = (1 - u - s + v, s - v, u - v, v)
            ranges = mp.zeros(4)
            mean = -read[0, 0]
            for weight, (i, j) in zip(
                weights, ((0, 0), (0, 1), (1, 0), (1, 1)), strict=True
            ):
                ranges += weight * matrices[i][j]
                mean += weight * read[i, j]
            logs = mp.matrix(4, 4)
            for i, j in combinations(range(4), 2):
                logs[i, j], logs[j, i] = mp.log(ranges[i, j]), mp.log(ranges[j, i])
            clock = [
                sum(logs[i, j] - logs[j, i] for j in range(4)) / 4 for i in range(4)
            ]
            sq = mp.matrix(4, 4)
            for i, j in combinations(range(4), 2):
                sq[i, j] = sq[j, i] = mp.exp(
                    logs[i, j] + logs[j, i] - clock[i] - clock[j]
                )

            def gram(indices):
                i, *rest = indices
                g = mp.matrix(
                    [[(sq[i, j] + sq[i, k] - sq[j, k]) / 2 for k in rest] for j in rest]
                )
                return mp.sqrt(mp.det(g)) / mp.factorial(len(rest))

            area = sum(gram(face) for face in combinations(range(4), 3))
            volume = gram(range(4))
            entropy = -sum(y * mp.log(y) for y in ((1 + mean) / 2, (1 - mean) / 2))
            return area, volume, entropy

        x = tuple(map(mp.mpf, (u, s, u * s)))

        def total(*args):
            a, v, k = observables(*args)
            return k + eta * a

        def vol(*args):
            return observables(*args)[1]

        columns = np.eye(3, dtype=int)

        def gradient(f):
            return np.array(
                [float(mp.diff(f, x, tuple(int(v) for v in axis))) for axis in columns]
            )

        grad_v, grad_f = gradient(vol), gradient(total)
        multiplier = np.dot(grad_v, grad_f) / np.dot(grad_v, grad_v)

        def lagrangian(*args):
            a, v, k = observables(*args)
            return k + eta * a - multiplier * v

        hess = np.array(
            [
                [
                    float(mp.diff(lagrangian, x, tuple(int(v) for v in (i + j))))
                    for j in columns
                ]
                for i in columns
            ]
        )
    tangent = null_space(grad_v.reshape(1, 3))
    return dict(
        fixed_volume_gradient_norm=float(np.linalg.norm(tangent.T @ grad_f)),
        eigenvalues=np.linalg.eigvalsh(tangent.T @ hess @ tangent).tolist(),
        volume_multiplier=multiplier,
    )


def run():
    banks = radar_banks()
    corners = survey_corners(banks, 0, 6)
    clicks = bridge_corners()
    click = clicks[:, :, 0]
    u = 0.25

    def mismatch(s):
        return gains(corners, click, u, s)["fraction"] - 1

    lo, hi = 0.25, 0.5
    s = brentq(mismatch, lo, hi, xtol=1e-13)
    station = gains(corners, click, u, s)
    eta = station["ordinary_eta"]
    held = []
    for tick, horizon, leg, odd, paired in (
        (0, 6, 0, u, 0.25),
        (0, 6, 0, u, 0.5),
        (0, 6, 0, 0.375, s),
        (0, 4, 0, u, s),
        (1, 6, 0, u, s),
        (0, 6, 2, u, s),
    ):
        m = gains(survey_corners(banks, tick, horizon), clicks[:, :, leg], odd, paired)
        held.append(
            dict(
                tick=banks[1, False][tick]["step"],
                horizon=horizon,
                leg=pairing.TIMES[leg],
                odd=odd,
                pairing=paired,
                transfer_at_frozen_eta=-eta * m["area"][1] / m["entropy"][1],
                ordinary_at_frozen_eta=-eta * m["area"][0] / m["entropy"][0],
            )
        )
    return dict(
        odd=u,
        selected_pairing=s,
        bracket=[lo, hi],
        endpoint_mismatches=[mismatch(lo), mismatch(hi)],
        calibrated=station,
        constrained=stationary_hessian(corners, click, u, s, eta),
        held_at_frozen_calibration=held,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
