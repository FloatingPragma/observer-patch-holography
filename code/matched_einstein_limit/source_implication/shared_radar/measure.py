"""Independent recovery of the registered survey from each probe's return histogram."""

import json
from itertools import permutations
from pathlib import Path

import numpy as np

from .. import signed_transfers
from ..full_algebra import station
from ..full_algebra.enclosure import geometry_jet, iv
from .certificate import source as enclosure


def effects():
    """Complete two-leg arrival probabilities on fresh copies. Quiet is the complement."""
    iso = station.iso
    board = iso.pairing.formation.federation(1)
    ports = iso.pairing.formation.region(board)
    menus = {}
    for charge in (-1, 1):
        by_parity = []
        for odd in (False, True):
            h = iso.ranging.generator(board, ports, charge, 8, odd)
            u = iso.ranging.kernels(h, ports, 0.5, 24)
            by_edge = []
            for i, j in permutations(range(4), 2):
                out = iso.ranging.first_hits(u, i, j)
                back = iso.ranging.first_hits(u, j, i)
                by_time = []
                for n in range(2, 25):
                    amplitudes = np.array([back[n - k] @ out[k] for k in range(1, n)])
                    by_time.append(
                        np.einsum("kai,kaj->ij", amplitudes.conj(), amplitudes)
                    )
                by_edge.append(by_time)
            by_parity.append(by_edge)
        menus[charge] = np.stack(by_parity, axis=2)
    times = np.arange(2, 25)
    target = 3 * np.eye(4) + np.einsum(
        "t,ptbij->pbij", np.minimum(times / 4 - 3, 0), (menus[-1] + menus[1]) / 2
    )
    return menus, target


def calibrate(menu, target):
    """A fixed affine dual minimizing a bound on the per-bin variance."""
    effects = menu.reshape(-1, 2, 4, 4)
    directions = station.directions()
    base = np.stack([np.eye(4) / 8] * 2)
    design = np.einsum("ibjk,dbkj->id", effects, directions).real
    target_design = np.einsum("ibjk,dbkj->id", target, directions).real
    max_prob = np.linalg.eigvalsh(effects).max(axis=(1, 2))
    scale = np.sqrt(np.maximum(max_prob, 1e-15))
    dual = (target_design @ np.linalg.pinv(design / scale[:, None])) / scale[None, :]
    offset = station.pair(target, base) - dual @ station.pair(effects, base)
    residual = np.linalg.norm(dual @ design - target_design, axis=1)
    # The traceless part of a normalized density has HS norm less than one.
    return offset, dual, max_prob, residual


def bounds(weights, max_prob, shots, log_level):
    """Bernstein for independent edge launches; time bins on one edge are exclusive."""
    weights = weights.reshape(-1, 12, 23)
    width = (
        np.maximum(weights.max(axis=2), 0) - np.minimum(weights.min(axis=2), 0)
    ).max(axis=1)
    variance = np.sum(weights**2 * max_prob.reshape(1, 12, 23), axis=(1, 2))
    return np.sqrt(2 * variance * log_level / shots) + 2 * width * log_level / (
        3 * shots
    )


def preparations():
    data = json.loads(Path(station.__file__).with_name("results.json").read_text())
    x = np.array(data["x_real"]) + 1j * np.array(data["x_imag"])
    right = signed_transfers.right_coin(x)
    pauli = station.phase.AXES[1:]
    vec = np.array([np.trace(right @ a).real for a in pauli])
    axis = sum(v * a for v, a in zip(vec, pauli)) / np.linalg.norm(vec)
    left = np.trace(x.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    states = dict(
        equilibrium=x,
        ordinary=np.stack([x.sum(axis=0) / 2] * 2),
        pairing=np.stack([np.kron(b, right) for b in left]),
    )
    for label, theta in (("plus", np.pi / 2), ("minus", -np.pi / 2)):
        u = signed_transfers.rotation(axis, theta)
        states[label] = u @ x @ u.conj().T
    return states


def experiment():
    menus, target = effects()
    states = preparations()
    target_r = station.pair(target, states["equilibrium"])
    # One independently certified survey volume and its local matching rotations.
    archived = json.loads(
        Path(signed_transfers.__file__).with_name("signed_transfers.json").read_text()
    )
    angles = {arm["source"]: arm["angle_per_unit"] for arm in archived["tests"]}
    seed = 894031
    rng = np.random.default_rng(seed)
    # For each dose cover two species, the calibration bank's twelve means,
    # and twelve means plus two projected values in each of five fresh banks.
    log_level = np.log(2 * 2 * 2 * (12 + 5 * 14) / 0.01)
    out = {}
    for q, menu in menus.items():
        offset, dual, pmax, error = calibrate(menu, target)
        runs = {}
        for shots in (2**42, 2**50):
            counts = {}

            def observe(label, x):
                p = station.pair(menu, x)
                complete = np.column_stack([p, 1 - p.sum(axis=1)])
                assert complete.min() >= 0
                hits = np.array([rng.multinomial(shots, edge) for edge in complete])
                counts[label] = hits.tolist()
                return offset + dual @ (hits[:, :-1].ravel() / shots)

            range_bound = bounds(dual, pmax, shots, log_level) + error
            calibration = observe("reference", states["equilibrium"])
            gradient = station.range_gradient(calibration)
            # These intervals cover the unknown range Jacobian. The separate
            # reference bank also fixes the two projection rows before the held draws.
            intervals = [
                enclosure.ri_band(float(z), float(e))
                for z, e in zip(calibration, range_bound)
            ]
            iv.dps = 30
            geometric = geometry_jet(intervals)
            grad_error = np.array(
                [
                    [float(abs(u - enclosure.ri(float(v))).b) for u, v in zip(jet.g, g)]
                    for jet, g in zip(geometric, gradient)
                ]
            )
            shape_bound = (
                bounds(gradient @ dual, pmax, shots, log_level)
                + np.abs(gradient) @ error
            )
            observations = {label: observe(label, x) for label, x in states.items()}
            baseline = observations["equilibrium"]
            rotation = (observations["plus"] - observations["minus"]) / 2
            matched = {}
            for label, key in (("ordinary", "ordinary"), ("pairing", "collar")):
                alpha = angles[key]
                dr = observations[label] - baseline + alpha * rotation
                projected = gradient @ dr
                width = (2 + abs(alpha)) * shape_bound + grad_error @ (
                    abs(dr) + (2 + abs(alpha)) * range_bound
                )
                matched[label] = dict(
                    survey_area_slope=projected[0],
                    survey_volume_slope=projected[1],
                    simultaneous_sampling_halfwidth=width.tolist(),
                )
            runs[str(shots)] = dict(
                counts=counts,
                inferred_equilibrium_ranges=calibration.tolist(),
                range_halfwidths=range_bound.tolist(),
                inferred_equilibrium_area_volume=station.geometry(calibration).tolist(),
                matched=matched,
            )
        out[str(q)] = dict(numerical_design_error=float(max(error)), exposures=runs)
    return dict(
        seed=seed,
        clock_step=0.5,
        timeout=12,
        target_timeout=6,
        common_population_area_volume=station.geometry(target_r).tolist(),
        charges=out,
    )


if __name__ == "__main__":
    print(json.dumps(experiment(), indent=2, allow_nan=False))
