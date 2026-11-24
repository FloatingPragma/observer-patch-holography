"""Test the released, anchored vacancy on the connected one-vacancy generators."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import anchored_vacancy as anchor

SHOTS = 2**24
SEED = 178036
STEP = 0.25


def band_edges(board, h):
    # For E between these bounds H55-E is negative and the Schur complement
    # H66-E-H65(H55-E)^(-1)H56 is positive. Off-block hopping cannot fill the gap.
    blocks = []
    for degree in (5, 6):
        indices = np.flatnonzero(board.degree == degree)
        block = h[indices][:, indices]
        center = block.diagonal()
        radius = np.asarray(abs(block).sum(axis=1)).ravel() - abs(center)
        blocks.append([float(min(center - radius)), float(max(center + radius))])
    if blocks[0][1] >= blocks[1][0]:
        raise ValueError("Degree blocks do not certify a gap")
    return blocks


def count(wave, success, rng):
    # The signed vacancy echo identifies its occupied port on the unit-filled host.
    inside, outside = (float(np.vdot(v, v).real) for v in (wave[:12], wave[12:]))
    bins = np.array(
        [
            success * inside,
            success * outside,
            success * max(0, 1 - inside - outside),
            1 - success,
        ]
    )
    bins = np.round(bins, 13)
    return rng.multinomial(SHOTS, bins / bins.sum()).tolist()


def panel(board, h, source, success, launch, anchor_read, rng):
    e, v = np.linalg.eigh(h.toarray())
    rest = float(np.vdot(source, h @ source).real)
    np.testing.assert_allclose(h @ source, rest * source, atol=2e-12)
    target = float(np.expm1(np.log(2) / 2) * anchor_read)
    delta = e - rest
    kicked = source.copy()
    kicked[launch] *= -1
    coeff = v.T @ kicked
    distance = abs(delta - target)
    idx = int(np.argmin(distance))
    tol = 0.01 * target
    # Projection onto the entire degenerate line is basis independent.
    flat = abs(delta) < 1e-10
    records = []
    for time in (0, 8, 32, 128):
        phase = np.exp(-1j * time * delta)
        for impulse, state in ((False, source), (True, kicked)):
            wave = v @ (phase * (v.T @ state))
            records.append(
                {
                    "exposure": time,
                    "pulse": impulse,
                    "counts": count(wave, success, rng),
                }
            )
    filters = []
    for bins in (512, 4096):
        # All histories have equal total duration. The source accumulates rest
        # phase until its one local impulse. The pointer supplies a Fourier port.
        x = (delta - target) * STEP / (2 * np.pi)
        kernel = np.exp(-1j * np.pi * (bins - 1) * x) * np.sinc(bins * x) / np.sinc(x)
        if bins == 512:
            phases = np.exp(-2j * np.pi * np.outer(x, np.arange(bins)))
            np.testing.assert_allclose(kernel, phases.mean(axis=1), atol=2e-13)
        filters.append(
            {
                "appointments": bins,
                "counts": count(v @ (kernel * coeff), success, rng),
                "accepted_probability_given_source": float(
                    np.linalg.norm(kernel * coeff) ** 2
                ),
            }
        )
    return {
        "held_energy": rest,
        "rest_anchor_clock_read": anchor_read,
        "target_kinetic": target,
        "target_total_in_host_reference": rest + target,
        "closest_spectral_distance": float(distance[idx]),
        "closest_kinetic_line": float(delta[idx]),
        "number_of_eigenvalues_in_one_percent_window": int(sum(distance <= tol)),
        "kicked_probability_in_one_percent_window": float(
            sum(abs(coeff[distance <= tol]) ** 2)
        ),
        "kicked_weight_remaining_at_rest_energy": float(sum(abs(coeff[flat]) ** 2)),
        "rest_energy_component_outside_home": float(
            np.linalg.norm((v[:, flat] @ coeff[flat])[12:]) ** 2
        ),
        "homogeneous_energy_shift": "Any common sector shift cancels from delta",
        "flight": records,
        "coherent_kick_time_filter": filters,
    }


def run():
    original = json.loads(
        Path(anchor.__file__).with_name("material_reference.json").read_text()
    )
    clock = original["native"][0]["panels"]["6"]["vacancy"]
    board = anchor.federation(0)
    (source, _), success, prep = anchor.preparation(board)
    rng = np.random.default_rng(SEED)
    launch = prep["launch_local_port"]
    studies = {}
    for name, strength, read in (
        ("registered", anchor.K, clock["frequency_per_clock_gap"]),
        ("half_contrast", anchor.K / 2, None),
    ):
        h = anchor.defect_generator(board, -1, strength, 0)
        ref = float(np.vdot(source, h @ source).real) if read is None else read
        studies[name] = panel(board, h, source, success, launch, ref, rng)
    bounds = []
    for level in range(3):
        board = anchor.federation(level)
        h = anchor.defect_generator(board, -1, anchor.K, 0)
        bounds.append(
            {
                "cells": board.cells,
                "degree_block_enclosures": band_edges(board, h),
                "vacancy_floor_lower_bound": float(
                    min(
                        h.diagonal()
                        - (np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal()))
                    )
                ),
                "vacancy_floor_trial_upper_bound": float(min(h.diagonal())),
            }
        )
    return {
        "shots": SHOTS,
        "seed": SEED,
        "source_success": float(success),
        "launch_port": launch,
        "clock_read_standard_error": clock["frequency_standard_error"],
        "time_step": STEP,
        "bins": [
            "accepted Fourier output, home",
            "accepted Fourier output, outside",
            "other timing output",
            "source failed",
        ],
        "flight_convention": "Unfiltered flights accept every timing output",
        "boards": bounds,
        "studies": studies,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        from ..electron_field.replay import compare

        compare(
            result,
            json.loads(Path(__file__).with_suffix(".json").read_text()),
            "localized_release",
        )
        print(
            "Reproduced the native vacancy release, spectral enclosures and half-contrast control."
        )
