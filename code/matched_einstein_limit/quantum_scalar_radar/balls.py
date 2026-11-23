"""Measure equal-volume area deficits on a local spacetime-geodesic rest slice."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import measure

TICK = 2.0**-48
RADII = (0.4, 0.2, 0.1)
ANGULAR_STEP = 0.00025
APPOINTMENT = 0.60


def grid(order):
    z, weight = leggauss(order)
    longitude = (np.arange(2 * order) + 0.5) * math.pi / order
    z, angle = np.broadcast_arrays(z[:, None], longitude)
    radial = np.sqrt(1 - z * z)
    n = np.stack((radial * np.cos(angle), radial * np.sin(angle), z), axis=-1).reshape(
        -1, 3
    )
    east = np.stack((-np.sin(angle), np.cos(angle), np.zeros_like(z)), axis=-1).reshape(
        -1, 3
    )
    north = np.cross(n, east)
    return n, east, north, np.repeat(weight, 2 * order) * math.pi / order


def congruence(source, order):
    n, east, north, weights = grid(order)
    c, s = math.cos(ANGULAR_STEP), math.sin(ANGULAR_STEP)
    rays = np.concatenate(
        [n, c * n + s * east, c * n - s * east, c * n + s * north, c * n - s * north]
    )
    start = source(APPOINTMENT)
    momentum = rays * np.exp(start[:3])
    values = np.zeros((len(rays), 5))
    values[:, 0] = APPOINTMENT

    def rates(_, flat):
        state = flat.reshape(-1, 5)
        metric = source(state[:, 0])
        length = np.exp(metric[:3].T)
        # The forward source supplies the connection to the measuring apparatus.
        conjugate = metric[3:6].T
        hubble = (
            measure.classical.COUPLING
            / (measure.VOL0 * length.prod(axis=1)[:, None])
            * (conjugate - 0.5 * conjugate.sum(axis=1)[:, None])
        )
        v = momentum / length**2
        return np.column_stack(
            (state[:, 1], -np.sum(hubble * momentum * v, axis=1), v)
        ).ravel()

    path = solve_ivp(
        rates,
        (0, 0.45),
        values.ravel(),
        rtol=2e-12,
        atol=2e-13,
        method="DOP853",
        dense_output=True,
    )
    if not path.success:
        raise RuntimeError(path.message)

    def meters(radius, interior):
        coordinates = path.sol(radius).reshape(5, len(n), 5)
        centre, ep, em, np_, nm = coordinates
        lengths = np.exp(source(centre[:, 0])[:3, :].T)
        derivatives = []
        for plus, minus in ((ep, em), (np_, nm)):
            delta = (plus - minus) / (2 * s)
            derivatives.append(np.column_stack((delta[:, 0], delta[:, 2:] * lengths)))
        if interior:
            derivatives.append(
                np.column_stack((centre[:, 1], momentum[: len(n)] / lengths))
            )
        jacobi = np.stack(derivatives, axis=1)
        gram = np.einsum("nac,nbc,c->nab", jacobi, jacobi, [-1.0, 1.0, 1.0, 1.0])
        determinant = np.linalg.det(gram)
        if np.any(determinant < 0):
            raise ValueError(
                "Geodesic rest slice left the instrument's spacelike domain"
            )
        return np.rint(np.sqrt(determinant) / TICK).astype(np.int64).tolist()

    return meters, weights


def acquire(preparation):
    source = measure.propagate(preparation)
    interior, weights = congruence(source, 4)
    rim, boundary_weights = congruence(source, 5)
    nodes, radial_weights = leggauss(6)

    def volume(radius):
        samples = [interior(s, True) for s in radius * (1 + nodes) / 2]
        reading = radius / 2 * TICK * radial_weights @ (np.array(samples) @ weights)
        return float(reading), samples

    spheres = []
    for reference in RADII:
        target = 4 * math.pi * reference**3 / 3
        matched = brentq(
            lambda radius: volume(radius)[0] - target,
            0.97 * reference,
            1.03 * reference,
            xtol=1e-13,
        )
        interior_value, samples = volume(matched)
        faces = rim(matched, False)
        boundary = TICK * np.asarray(faces) @ boundary_weights
        estimated = 15 * (1 - boundary / (4 * math.pi * reference**2)) / reference**2
        spheres.append(
            dict(
                reference=reference,
                matched_radius=matched,
                interior_ticks=samples,
                boundary_ticks=faces,
                volume=interior_value,
                area=float(boundary),
                estimated_Guu=float(estimated),
            )
        )
    return spheres


def run(primary):
    result = {}
    fine = primary["surveys"][-1]
    for name in ("ground", "y", "z", "even_x"):
        reading = measure.read(fine["records"][name][1], primary["recovered_mass2"])
        prediction = (
            fine["calibration"]["vacuum"]
            + fine["calibration"]["coupling"] * reading["energy"]
        )
        result[name] = dict(
            spheres=acquire(measure.SOURCES[name]),
            stress_prediction=float(prediction),
            radar_Guu=float(reading["temporal"]),
        )
    # Constant-t slices of this metric are Euclidean after rescaling the axes.
    # Run the identical three-rod/two-rod cubatures on their exact inertial rulers.
    flat = []
    for r in RADII:
        nodes, wr = leggauss(6)
        density = np.rint(((1 + nodes) * r / 2) ** 2 / TICK) * TICK
        value = float(r / 2 * wr @ density * sum(grid(4)[-1]))
        area = float(round(r * r / TICK) * TICK * sum(grid(5)[-1]))
        flat.append(
            dict(
                reference=r,
                volume=value,
                area=area,
                Guu=15 * (1 - area / (4 * math.pi * r * r)) / r**2,
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "BALL_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        local_meter_tick=TICK,
        angular_step=ANGULAR_STEP,
        flat_control=flat,
        records=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = run(json.loads((measure.HERE / "results.json").read_text()))
    with args.output.open("x") as stream:
        json.dump(data, stream, separators=(",", ":"), allow_nan=False)
        stream.write("\n")
    for n, row in data["records"].items():
        print(
            n,
            row["stress_prediction"],
            row["radar_Guu"],
            [s["estimated_Guu"] for s in row["spheres"]],
        )
