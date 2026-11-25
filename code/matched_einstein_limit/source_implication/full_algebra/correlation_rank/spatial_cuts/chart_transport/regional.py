"""Whole-region transfer at two nested domain radii and the same survey volume."""

import argparse
import json
from math import floor

import numpy as np
from mpmath import iv, mp
from scipy.linalg import qr

from .... import station
from . import apertures, held, mounts, transfer

cuts = mounts.cuts
POOLS = transfer.POOLS[3:]


def read():
    data = json.loads(apertures.EFFECTS.read_text())
    matrices = []
    labels = []
    for arm, patch in zip(data["charts"], apertures.apparatus()[-1], strict=True):
        assert arm["cells"] == patch["cells"] and arm["images"] == patch["images"]
        keys = [(r["charge"], r["reading"]) for r in arm["rows"]]
        assert keys == [(q, n) for q in (-1, 1) for n in mounts.REGIONAL]
        matrices.append(transfer.interval_rows(arm["rows"], data["scale"]))
        labels.append(keys)
    return matrices, labels


def study(pool, matrices, held_data, intervals, target):
    source, actual, base, means, changes, compared, radar = held_data
    all_keys = [(a, i) for c in pool for a in (2 * c, 2 * c + 1) for i in range(8)]
    m = np.array([matrices[a][i] for a, i in all_keys])
    # Select physical detectors from the material response only; neither the
    # stationary density nor any area or volume reading enters the selection.
    pivots = qr(cuts.mid(m).T, mode="economic", pivoting=True)[2][:31]
    keys = [all_keys[i] for i in sorted(pivots)]
    b = np.array([matrices[a][i] for a, i in keys])
    inv = transfer.center(b) ** -1
    R = transfer.intervals(transfer.publish(inv))
    defect = cuts.norm(np.eye(31, dtype=int) - R @ b)
    assert defect.b < cuts.ri("1e-35")
    norm = cuts.norm(R)
    sigma = (1 - defect) / norm
    # A fixed independent modular functional supplements the entirely
    # material reconstruction of the reference below.
    w = transfer.intervals(transfer.publish(target * inv)[0])
    gain = sum(abs(x) for x in w)
    source_hat = held.material_reference(inv, means, keys, base)
    log = -station.log_blocks(source_hat)
    reference_error = float(np.linalg.norm(source_hat - source))
    radar_error = float(max(abs(station.pair(radar, source_hat - source))))
    omitted = set(range(3)) - set(pool)
    cases = []
    for t, d in enumerate(actual):
        incoming = mp.matrix([float(changes[a, i, t]) for a, i in keys])
        inferred = np.array((inv * incoming).tolist(), dtype=float).ravel()
        physical = float(station.pair(log, held.blocks(inferred)))
        controls = []
        for third in omitted:
            # Compare the original full outer mount, including its local cuts,
            # with a prediction based only on the other mounts' nested regions.
            prediction = np.array(
                (transfer.center(intervals[third]) * inv * incoming).tolist(),
                dtype=float,
            ).ravel()
            output = compared[third, :, t]
            center = transfer.intervals(
                transfer.publish(transfer.center(intervals[third]))
            )
            errors = cuts.norm(center) * defect + cuts.norm(intervals[third] - center)
            controls.append(
                dict(
                    held_chart=third,
                    relative_l2_error=float(
                        np.linalg.norm(prediction - output) / np.linalg.norm(output)
                    ),
                    deterministic_vector_error_per_unit_HS_upper=transfer.upper(errors),
                )
            )
        cases.append(
            dict(
                material_only_modular=physical,
                source_HS_error=float(np.linalg.norm(inferred - d)),
                full_held_charts=controls,
            )
        )
    return dict(
        mounts=list(pool),
        rank=31,
        physical_rows=keys,
        inverse_defect_upper=transfer.upper(defect),
        min_singular_lower_in_1e20=floor(transfer.endpoint(sigma, 0) * 10**20),
        source_noise_gain_upper=transfer.upper(norm / (1 - defect)),
        modular_price_l1_upper=transfer.upper(gain),
        reference_reconstruction_HS_error=reference_error,
        maximum_radar_range_error=radar_error,
        cases=cases,
    )


def run():
    iv.dps = 75
    mp.dps = 80
    matrices, labels = read()
    intervals = transfer.read_native()
    source, dirs, pairs, radar = held.preparations()
    fair = np.tile(np.eye(4), (2, 1, 1)) / 8
    flow = mounts.float_effects(
        np.concatenate(([source, fair], dirs)), apertures.apparatus, regional_only=True
    )
    observations = np.array([[arm[k] for k in keys] for keys, arm in zip(labels, flow)])
    old = json.loads((mounts.HERE / "held.json").read_text())
    assert (
        max(
            abs(a["entropy"] - b["entropy"]) for a, b in zip(pairs, old["preparations"])
        )
        < 1e-8
    )
    direct = [mounts.source_components(d) for d in dirs]
    previous = np.asarray(old["readings_by_chart_and_preparation"])
    target = transfer.center(transfer.archive_operator()[None, :])
    measurements = (
        source,
        direct,
        fair,
        observations[:, :, :2],
        observations[:, :, 2:],
        previous,
        radar,
    )
    return dict(
        preparations=pairs,
        absolute_means_and_held_slopes=observations.tolist(),
        instruments=[
            study(pool, matrices, measurements, intervals, target) for pool in POOLS
        ],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    r = run()
    if a.write:
        (mounts.HERE / "regional.json").write_text(json.dumps(r, indent=2) + "\n")
    print(json.dumps(r["instruments"]))
