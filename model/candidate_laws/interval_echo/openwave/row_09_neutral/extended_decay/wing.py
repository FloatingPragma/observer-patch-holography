"""Read the downstream material energy and the power crossing into it."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from ..joined_decay.measure import core_pulse
from ..passive_decay.material import columns, edge
from ..passive_decay.measure import SHOTS, compare_records, hamiltonians, physical
from ..passive_decay.readers import Z, charge, circuit, clocks, rotate
from .material import route
from .measure import HERE, rebuild

WAITS = (0, 4, 8, 16, 32, 64)
PRIMARY = 32
BINS = 1024
SEED = 26112709


def installation():
    board, dimers, heads, _, _ = route((3,), crossing=True)
    native = rebuild(board, dimers, heads, False, True)
    crossing = edge(dimers[2][1], dimers[3][0])
    cut = {**heads, crossing: 32}
    changed = dict.fromkeys(heads, 32)
    settings = {"native": native}
    for name, costs in (("boundary_cut", cut), ("no_reversal", changed)):
        potential = physical(board, dimers, costs)
        settings[name] = (
            native - np.diag(np.diag(native)) + np.diag(potential - potential[0])
        )
    return board, dimers, settings


def boundary_power(h, stock, stop=PRIMARY, bins=BINS):
    """Measure the C,D rails in two product frames, averaged over the midpoint clock."""
    _, wing, _, _ = hamiltonians(4)
    local = np.kron(np.eye(8), wing)
    occupied = np.diag([0.0, 1.0])
    vacant = np.eye(2) - occupied
    boundary = -np.kron(np.eye(4), np.kron(vacant, np.kron(occupied, np.eye(2))))
    y = np.array([[0.0, -1j], [1j, 0.0]])
    y_frame = np.diag([1.0, 1j]) @ rotate(np.pi / 4) @ Z
    np.testing.assert_allclose(y_frame.conj().T @ y @ y_frame, Z, atol=1e-13)
    times = (np.arange(bins) + 0.5) * stop / bins
    values, vectors = eigh(h)
    waves = vectors @ (
        np.exp(-1j * values[:, None] * times) * (vectors.T @ stock)[:, None]
    )
    waves = waves.reshape(4, 2, 2, 2, bins)
    output = {}
    for name, energy, measured in (
        ("bare", local, np.kron(vacant, y / 2)),
        ("with_boundary", local + boundary, np.kron(y / 2, occupied)),
    ):
        current = np.kron(np.eye(4), np.kron(measured, np.eye(2)))
        np.testing.assert_allclose(1j * (h @ energy - energy @ h), current, atol=1e-13)
        if name == "bare":
            amplitudes = np.einsum("dj,acdet->acjet", y_frame.conj(), waves)
            joint = (abs(amplitudes) ** 2).sum(axis=(0, 3)).mean(axis=-1)
            outcomes = [joint[0, 1], joint[1].sum(), joint[0, 0]]
        else:
            amplitudes = np.einsum("cj,acdet->ajdet", y_frame.conj(), waves)
            joint = (abs(amplitudes) ** 2).sum(axis=(0, 3)).mean(axis=-1)
            outcomes = [joint[1, 1], joint[:, 0].sum(), joint[0, 1]]
        slope = h @ current - current @ h
        curvature = h @ slope - slope @ h
        error = stop * (stop / bins) ** 2 * np.linalg.norm(curvature, 2) / 24
        output[name] = dict(
            negative_zero_positive=np.round(outcomes, 13).tolist(),
            midpoint_error_bound=float(error),
        )
    return output


def acquire():
    board, dimers, settings = installation()
    core, wing, _, _ = hamiltonians(4)
    core_values, core_basis = eigh(core)
    wing_values, wing_basis = eigh(wing)
    core_reader = circuit(core_basis, 32)
    wing_reader = circuit(wing_basis, 33)
    reader = np.kron(core_reader, np.kron(np.eye(2), wing_reader))
    _, native_basis = eigh(settings["native"])
    cold = native_basis[:, 0]
    pulse = np.kron(circuit(core_pulse(core_basis), 32), np.eye(8))
    stocks = {"pulse": pulse @ cold, "idle": cold}
    charge_meters = charge(board, (*dimers[:2], *dimers[-2:]))
    assert len(charge_meters) == 2
    for meter in charge_meters:
        labels = meter["word_00_to_22"]
        assert (labels[0], labels[1], labels[3], labels[4]) == (7, 0, 0, 1)
        assert meter["coherent_amplitude_error"] < 1e-12
    panels = {}
    for name, h in settings.items():
        values, vectors = eigh(h)
        panels[name] = {}
        for arm, stock in stocks.items():
            rows = []
            for wait in WAITS:
                evolved = vectors @ (np.exp(-1j * values * wait) * (vectors.T @ stock))
                pair = abs(reader.T @ evolved).reshape(4, 2, 4) ** 2
                words = (abs(evolved.reshape(4, 2, 2, 2)) ** 2).sum(axis=(0, 3))
                rows.append(
                    dict(
                        core_by_wing=pair.sum(axis=1).round(13).tolist(),
                        crossing_c_by_d=words.round(13).tolist(),
                    )
                )
            panels[name][arm] = rows
    return dict(
        protocol_sha256=sha256((HERE / "WING.md").read_bytes()).hexdigest(),
        dimers=dimers,
        waits=list(WAITS),
        native_and_changed_column_errors=[columns(), columns(True)],
        parent_gap=float(np.ptp(core_values)),
        isolated_wing_gaps=(wing_values - wing_values[0]).tolist(),
        clock_simulation=dict(
            seed=SEED,
            shots_per_quadrature=SHOTS,
            fringes=clocks(wing, 33, np.random.default_rng(SEED)),
        ),
        charge_readers=charge_meters,
        panels=panels,
        power=dict(
            stop=PRIMARY,
            bins=BINS,
            arms={a: boundary_power(settings["native"], s) for a, s in stocks.items()},
        ),
    )


def contrast(record, setting, wait):
    index = record["waits"].index(wait)
    pulse, idle = (record["panels"][setting][arm][index] for arm in ("pulse", "idle"))
    levels = (np.asarray(pulse["core_by_wing"]) - idle["core_by_wing"]).sum(axis=0)
    occupancies = np.asarray(pulse["crossing_c_by_d"]) - idle["crossing_c_by_d"]
    measured_gaps = [
        0.0,
        *(line["gap"] for line in record["clock_simulation"]["fringes"]),
    ]
    bare = levels @ measured_gaps
    boundary = -occupancies[0, 1] if setting == "native" else 0.0
    return levels, np.array([bare, bare + boundary])


def decode(record):
    for name, stocks in record["panels"].items():
        for rows in stocks.values():
            assert len(rows) == len(record["waits"])
            for row in rows:
                for key, shape in (
                    ("core_by_wing", (4, 4)),
                    ("crossing_c_by_d", (2, 2)),
                ):
                    values = np.asarray(row[key])
                    assert values.shape == shape and values.min() >= 0
                    np.testing.assert_allclose(values.sum(), 1, atol=1e-11, rtol=0)
        levels, energies = contrast(record, name, PRIMARY)
        if name != "native":
            np.testing.assert_allclose(levels, 0, atol=1e-11, rtol=0)
            np.testing.assert_allclose(energies, 0, atol=1e-11, rtol=0)
        print(
            name,
            "level excess",
            np.round(levels, 7),
            "energy excess",
            np.round(energies, 7),
        )
    arms = record["power"]["arms"]
    final_levels, final = contrast(record, "native", PRIMARY)
    initial_levels, initial = contrast(record, "native", 0)
    endpoints = final - initial
    radii = [line["radius"] for line in record["clock_simulation"]["fringes"]]
    clock_error = np.abs((final_levels - initial_levels)[1:]) @ radii
    assert np.min(endpoints) > clock_error, "The downstream energy sign is unresolved."
    print(
        "excited levels above a one percentage point excess",
        int(sum(final_levels[1:] > 0.01)),
    )
    for i, allocation in enumerate(("bare", "with_boundary")):
        integral = 0.0
        bound = 0.0
        for arm, sign in (("pulse", 1), ("idle", -1)):
            panel = arms[arm][allocation]
            probs = np.asarray(panel["negative_zero_positive"])
            assert probs.shape == (3,) and min(probs) >= 0
            np.testing.assert_allclose(probs.sum(), 1, atol=1e-11, rtol=0)
            integral += sign * (probs[2] - probs[0]) * record["power"]["stop"] / 2
            bound += panel["midpoint_error_bound"]
        assert abs(integral - endpoints[i]) <= bound + clock_error + 1e-10
        print(
            allocation,
            "integrated pulse-minus-idle power",
            integral,
            "quadrature bound",
            bound,
            "endpoint clock half-width",
            clock_error,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    measured = acquire()
    if args.acquire:
        with args.acquire.open("x") as stream:
            json.dump(measured, stream, allow_nan=False, separators=(",", ":"))
            stream.write("\n")
    else:
        compare_records(json.loads((HERE / "wing.json").read_text()), measured)
    decode(measured)
