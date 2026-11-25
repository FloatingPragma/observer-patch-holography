"""Two equal-survey-volume source pairs and a held material chart for each chart pair."""

import argparse
import json

import numpy as np
from mpmath import iv, mp
from scipy.optimize import brentq

from .... import station
from . import mounts, transfer

cuts = mounts.cuts


def preparations():
    raw = json.loads((cuts.survey.HERE.parent / "results.json").read_text())
    x = np.array(raw["x_real"]) + 1j * np.array(raw["x_imag"])
    radar = station.radar_menus()[6]
    ranges = station.pair(radar, x)
    volume = station.geometry(ranges)[1]
    basis = station.directions()
    grad = station.range_gradient(ranges)[1]
    ds = np.stack([station.pair(radar, d) for d in basis], axis=1)
    g = grad @ ds
    adjust = np.einsum("a,abij->bij", g / np.linalg.norm(g), basis)
    da = station.pair(radar, adjust)
    ds_a = g @ (g / np.linalg.norm(g))
    variations = []
    obs = []
    for indices, coeff in (
        ([0, 8, 22, 29, 30], [1.0, 0.8, -0.5, 0.4, 0.3]),
        ([1, 6, 19, 23, 30], [-0.35, 0.7, -1.0, -0.6, 0.8]),
    ):
        raw_dir = sum(c * basis[i] for i, c in zip(indices, coeff))
        v = raw_dir - adjust * (grad @ station.pair(radar, raw_dir)) / ds_a
        v /= np.linalg.norm(v)
        arms = []
        vols = []
        areas = []
        for sign in (-1, 1):
            step = 4e-5
            change = sign * step * v
            arc = ranges + station.pair(radar, change)
            shift = brentq(
                lambda t: station.geometry(arc + t * da)[1] - volume,
                -1e-6,
                1e-6,
                xtol=1e-16,
            )
            delta = change + shift * adjust
            assert (
                np.linalg.norm(delta) < 1e-4
                and np.linalg.eigvalsh(x + delta).min() > 0.02
            )
            arms.append(delta)
            a, vol = station.geometry(station.pair(radar, x + delta))
            vols.append(float(vol))
            areas.append(float(a))
        variation = (arms[1] - arms[0]) / (2 * step)
        variations.append(variation)
        direct = (station.entropy(x + arms[1]) - station.entropy(x + arms[0])) / (
            2 * step
        )
        areal = -raw["ordinary_eta"] * (areas[1] - areas[0]) / (2 * step)
        assert max(abs(u - volume) for u in vols) < 2e-13 and abs(areal - direct) < 1e-7
        obs.append(
            dict(
                arm_step=step,
                arm_volumes=vols,
                entropy=float(direct),
                ordinary_calibrated_area=float(areal),
            )
        )
    return x, np.asarray(variations), obs, radar


def blocks(coefficients):
    # Inverse of source_components in its trace-zero Hilbert-Schmidt basis.
    joint = station.iso.basis(True)
    local = np.asarray(mounts.full_source.local_basis(), dtype=complex)
    parity = np.stack((np.eye(4), -np.eye(4))) / np.sqrt(8)
    return (
        np.einsum("ba,aij->bij", coefficients[:18].reshape(2, 9), joint)
        + np.einsum("ba,aij->bij", coefficients[18:30].reshape(2, 6), local)
        + coefficients[-1] * parity
    )


def material_reference(inverse, means, keys, fair):
    # The two inputs are absolute material expectation records from separate stocks.
    difference = mp.matrix(
        [mp.mpf(float(means[a, j, 0])) - mp.mpf(float(means[a, j, 1])) for a, j in keys]
    )
    x = fair + blocks(np.array((inverse * difference).tolist(), dtype=float).ravel())
    if np.linalg.eigvalsh(x).min() <= 0:
        raise ValueError(
            "the absolute material records do not resolve a faithful reference"
        )
    return x


def run(report):
    iv.dps = 75
    mp.dps = 80
    data = json.loads((mounts.HERE / "effects.json").read_text())
    labels = [
        [(r["charge"], r["reading"]) for r in arm["rows"]] for arm in data["charts"]
    ]
    reference, dirs, pairs, radar = preparations()
    fair = np.tile(np.eye(4), (2, 1, 1)) / 8
    readings = mounts.float_effects(np.concatenate(([reference, fair], dirs)))
    changes = np.array(
        [[record[k] for k in keys] for keys, record in zip(labels, readings)]
    )
    absolute, changes = changes[:, :, :2], changes[:, :, 2:]
    intervals = transfer.read_native()
    direct = [mounts.source_components(d) for d in dirs]
    reports = []
    for instance in report["instruments"]:
        pool = instance["chart_indices"]
        ids = instance["physical_rows"]
        h = transfer.center(np.array([intervals[a][j] for a, j in ids]))
        inverse = h**-1
        weights = mp.matrix([mp.mpf(w) for w in instance["prices"]])
        measured = material_reference(inverse, absolute, ids, fair)
        log_meter = -station.log_blocks(measured)
        rebuilt_ranges = station.pair(radar, measured)
        reference_read = dict(
            min_eigenvalue=float(np.linalg.eigvalsh(measured).min()),
            source_HS_error=float(np.linalg.norm(measured - reference)),
            entropy=float(station.entropy(measured)),
            twelve_directed_range_max_error=float(
                np.max(abs(rebuilt_ranges - station.pair(radar, reference)))
            ),
        )
        arms = []
        for t in range(len(pairs)):
            received = mp.matrix([float(changes[a, j, t]) for a, j in ids])
            reader = float(mp.fdot(weights, received))
            d = np.array((inverse * received).tolist(), dtype=float).ravel()
            source_error = np.linalg.norm(d - direct[t])
            out = []
            if len(pool) == 2:
                for other in set(range(3)) - set(pool):
                    prediction = np.array(
                        (
                            transfer.center(intervals[other]) * inverse * received
                        ).tolist(),
                        dtype=float,
                    ).ravel()
                    actual = changes[other, :, t]
                    out.append(
                        dict(
                            chart=other,
                            absolute_l2_error=float(
                                np.linalg.norm(prediction - actual)
                            ),
                            actual_l2=float(np.linalg.norm(actual)),
                            relative_l2_error=float(
                                np.linalg.norm(prediction - actual)
                                / np.linalg.norm(actual)
                            ),
                        )
                    )
            arms.append(
                dict(
                    modular=reader,
                    material_reference_modular=float(
                        station.pair(log_meter, blocks(d))
                    ),
                    reconstructed_source_HS_error=float(source_error),
                    held_readings=out,
                )
            )
        reports.append(
            dict(chart_indices=list(pool), reference=reference_read, cases=arms)
        )
    return dict(
        preparations=pairs,
        readings_by_chart_and_preparation=changes.tolist(),
        reference_and_fair_stock_means=absolute.tolist(),
        instruments=reports,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    r = run(json.loads((mounts.HERE / "transfer.json").read_text()))
    if a.write:
        (mounts.HERE / "held.json").write_text(json.dumps(r, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in r.items() if k != "readings_by_chart_and_preparation"}
        )
    )
