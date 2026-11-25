"""Blind regional reconstruction, positive analyzer synthesis and reciprocal optical flights."""

import argparse
import json

import numpy as np
from mpmath import iv, mp

from .. import mounts, regional
from . import counts, crosscheck, optics, receiver as rx
import random


def unknown_preparations():
    # The emitter supplies states only to the two native flight simulators.
    # All calls to the receiver below use material detector records instead.
    p = json.loads((mounts.cuts.survey.HERE.parent / "results.json").read_text())
    x = np.array(p["x_real"]) + 1j * np.array(p["x_imag"])
    u = np.kron(np.diag([1.0, 1j]), np.eye(2))
    y = u @ x @ u.conj().T
    return {"stationary": x, "nearby": x + (y - x) / 8192, "quarter_turn": y}


def calibrations():
    matrices, _ = regional.read()
    sources = unknown_preparations()
    choices = counts.choose(matrices)
    costs = counts.score_ranges()
    streams = [random.Random(280003), random.Random(930047)]
    banks = {}
    programs = {}
    eff = {}
    source_error = {}
    for name, truth in sources.items():
        arms = []
        status = []
        errors = []
        banks[name] = []
        for position, keys in enumerate(choices):
            observations, bank = counts.acquire(
                matrices, keys, costs, truth, streams[position]
            )
            banks[name].append(bank)
            count, analyzer, report = rx.calibrate(
                keys, matrices, observations, (-1 if position == 0 else 1) * rx.SHIFT
            )
            arms.append(analyzer)
            status.append(report)
            errors.append(float(np.linalg.norm(analyzer - truth.transpose(0, 2, 1))))
        programs[name] = status
        eff[name] = arms
        source_error[name] = errors
    return (
        sources,
        eff,
        dict(
            material_record_order=["tracial", *sources],
            material_banks=banks,
            receivers=programs,
            reconstructed_analyzer_HS_errors=source_error,
            nearby_distance_HS=float(
                np.linalg.norm(sources["nearby"] - sources["stationary"])
            ),
            endpoint_material_rows=choices,
        ),
    )


def run():
    iv.dps = 60
    mp.dps = 80
    sources, effects, record = calibrations()
    rng = np.random.default_rng(60043)
    samples = {}
    for charge, bounded, cap, numerical in optics.kernels():
        groups = {}
        for name, source in sources.items():
            rows = {}
            trials = ("adaptive",) if name == "stationary" else ("adaptive", "frozen")
            for mode in trials:
                stock = effects[name if mode == "adaptive" else "stationary"]
                forward, reverse = optics.directional(bounded, source, stock)
                forward_f = optics.camera.transmission(numerical, source, stock[1])[
                    5:, :5
                ]
                reverse_f = optics.camera.transmission(numerical, source, stock[0])[
                    :5, 5:
                ]
                err = max(
                    np.max(abs(forward - forward_f)), np.max(abs(reverse - reverse_f))
                )
                assert err < 1e-13
                measured = optics.counts(forward_f, reverse_f, rng)
                asym = float(np.max(abs(forward - reverse.T)))
                # Bound includes native columns, final binary conversion,
                # and arithmetic in the stored nominal emission preparation.
                gain = []
                for f in (forward, reverse):
                    gain.append(
                        optics.camera.gain_bounds(np.maximum(f - 1e-14, 0), f + 1e-14)
                    )
                rows[mode] = dict(
                    max_pair_mismatch_upper=asym + 2e-14,
                    max_pair_mismatch_lower=max(0.0, asym - 2e-14),
                    native_directional_area_gain=gain,
                    counts=measured,
                )
            groups[name] = rows
        samples[str(charge)] = dict(
            unconditional_pixel_operator_norm_squared_upper=rx.transfer.upper(cap),
            sources=groups,
            archived_camera_audit=crosscheck.previous_banks(charge, bounded),
        )
    # Each endpoint has a separate finite bank. On their simultaneous Hoeffding
    # event the exact trace-one source error is bounded by the certified operator
    # defect, the implemented lottery residual and the native noise gain.
    worst = []
    for position in range(2):
        errors = []
        for name in sources:
            model = record["receivers"][name][position]
            data = record["material_banks"][name][position]
            errors.append(
                iv.mpf(model["implementation_and_inverse_error_upper"])
                + iv.mpf(model["source_gain"]) * iv.mpf(data["contrast_l2_error_upper"])
            )
        worst.append(max(x.b for x in errors))
    floor = min(min(p["min_count"] for p in rr) for rr in record["receivers"].values())
    # The reference lottery determines a strictly positive ball: each kept or
    # dephased target changes by at most twice the state HS norm. The product
    # Pauli frame has smallest singular value one.
    open_radius = iv.mpf(".001")
    margin = (
        iv.mpf(floor) / rx.TOTAL
        - 2 * (open_radius + iv.mpf("2.01e-7") + 2 * max(worst))
        - iv.mpf(rx.SHIFT + 1) / rx.TOTAL
    )
    assert margin.a > iv.mpf(".0005")
    for q, v in samples.items():
        bound = iv.mpf(v["unconditional_pixel_operator_norm_squared_upper"]) * sum(
            worst
        )
        v["simultaneous_mismatch_upper_from_material_counts"] = rx.transfer.upper(bound)
    return dict(
        **record,
        shots_per_emitter=optics.camera.SHOTS,
        cameras=samples,
        state_error_upper=[rx.transfer.upper(e) for e in worst],
        guaranteed_unknown_neighborhood_HS=".001",
        positive_weight_lower_millionths=int(rx.transfer.endpoint(margin, 0) * 10**6),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    if args.write:
        (rx.HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "endpoint_material_rows",
                    "reconstructed_analyzer_HS_errors",
                    "nearby_distance_HS",
                )
            }
        )
    )
