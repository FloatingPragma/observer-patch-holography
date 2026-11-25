"""Increase entropy by adding ordinary energy or by erasing pairing on one selected head source."""

import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry import (
    measure as phase,
)
from .full_algebra import station


def right_coin(x):
    return np.trace(x.sum(axis=0).reshape(2, 2, 2, 2), axis1=0, axis2=2)


def rotation(axis, theta):
    """Rotate the right head about its independently measured marginal Bloch vector."""
    gate = np.kron(np.eye(2), expm(-0.5j * theta * axis))
    return gate


def check_rail(axis, theta):
    identity, x, y, z = phase.AXES
    v = np.array([np.trace(a @ axis).real / 2 for a in (x, y, z)])
    phi, alpha = np.arctan2(v[1], v[0]), np.arccos(v[2])

    def rz(t):
        return expm(-0.5j * t * z)

    ry = rz(np.pi / 2) @ expm(-0.5j * alpha * x) @ rz(-np.pi / 2)
    orient = rz(phi) @ ry
    wanted = orient @ rz(theta) @ orient.conj().T
    # One calibrated native contact transfers the coin to a rail with Z
    # conjugation. Apply the conjugated Euler word, then reverse the transfer.
    exchange = phase.marked_gate(1, 8).reshape(4, 4)
    circuit = exchange.conj().T @ np.kron(z @ wanted @ z, identity) @ exchange
    return float(
        np.linalg.norm(
            circuit[:, :2] - np.vstack([expm(-0.5j * theta * axis), np.zeros((2, 2))])
        )
    )


def run():
    path = Path(station.__file__).with_name("results.json")
    record = json.loads(path.read_text())
    chi = np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    eta = record["ordinary_eta"]
    effects = station.radar_menus()[6]

    def read(x):
        r = station.pair(effects, x)
        return np.r_[station.geometry(r), station.entropy(x), np.trace(x[1]).real + 16]

    baseline = read(chi)
    r = right_coin(chi)
    bloch = np.array([np.trace(r @ a).real for a in phase.AXES[1:]])
    axis = sum(v * a for v, a in zip(bloch, phase.AXES[1:])) / np.linalg.norm(bloch)
    h = np.kron(np.eye(2), axis / 2)
    adjust = -1j * (h @ chi - chi @ h)
    left = np.trace(chi.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    product = np.stack([np.kron(block, r) for block in left])
    # A left sharp sign read with result forgotten equalizes parity blocks.
    dephased = np.stack([chi.sum(axis=0) / 2] * 2)
    jac = station.range_gradient(station.pair(effects, chi))

    def differential(d):
        return np.r_[
            jac @ station.pair(effects, d),
            -station.pair(d[None], station.log_blocks(chi))[0],
            np.trace(d[1]).real,
        ]

    dv = differential(adjust)
    tests = []
    for label, target in [("ordinary", dephased), ("collar", product)]:
        d = target - chi
        slope = differential(d)
        angle = -slope[1] / dv[1]
        matched = slope + angle * dv
        finite = []
        for t in (2e-7, 1e-7, 5e-8):
            source = chi + t * d

            def rotated(theta):
                u = rotation(axis, theta)
                return u @ source @ u.conj().T

            interval = sorted((0.0, 2 * angle * t))
            theta = brentq(
                lambda x: read(rotated(x))[1] - baseline[1], *interval, xtol=8e-15
            )
            active = rotated(theta)
            meas = read(active)
            diff = (meas - baseline) / t
            finite.append(
                dict(
                    probability=t,
                    angle=theta,
                    volume_difference=meas[1] - baseline[1],
                    signed_per_unit=diff.tolist(),
                    residual=diff[2] + eta * diff[0],
                    right_marginal_error=float(np.linalg.norm(right_coin(active) - r)),
                    transfer_gate_error=check_rail(axis, theta),
                )
            )
        tests.append(
            dict(
                source=label,
                raw=slope.tolist(),
                angle_per_unit=angle,
                matched=matched.tolist(),
                effective_eta=-matched[2] / matched[0],
                finite=finite,
            )
        )
    return dict(
        eta=eta,
        order=["survey_area", "survey_volume", "full_head_entropy", "parked_energy"],
        baseline=baseline.tolist(),
        unitary_control=dv.tolist(),
        tests=tests,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
