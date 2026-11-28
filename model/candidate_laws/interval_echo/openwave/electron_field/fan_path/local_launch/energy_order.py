"""Read the same signed face arrival on either side of the protected energy test."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from runpy import run_path

import numpy as np

from ...field_charge_records.asymmetric_clock import (
    STEP,
    TICKS,
    UNDERFLOW_PAD,
    checked_spectrum,
    interval,
    step_error,
    taylor_step,
)
from ..fixture import FanPath, archived_fixture
from ..measure import conjugacy
from .measure import HERE, band, pauli_record


def weight(columns):
    return float(np.vdot(columns, columns).real / 2)


def case(path, negative):
    conjugacy(path, negative)
    frame, delta, _ = band(path)
    block = path.blocks[3, 0]
    words = block.words[3]
    assert path.cells[3][0] == 1045
    indices = np.arange(block.rows.start, block.rows.stop).reshape(-1, len(words))
    mask = indices[:, [words.index(w) for w in (334, 365)]]
    assert [sum(path.field.bits.values()) ^ w for w in (334, 365)] == [177, 146]
    eta = float(Fraction(2**41, 38**21))
    proj_err = (2 + delta) * delta + eta + 1e-11
    vf = frame[mask.ravel()]
    bound = max(checked_spectrum(vf.T @ vf)) + 2 * proj_err
    records = []
    for rate in (4, 2):
        h = (path.parked + (rate / 4) * path.bridge_term).tocsr()
        _, _, truncation, rounding = step_error(h)
        wave = path.source().astype(complex)
        for tick in range(TICKS[-1] + 1):
            if tick in (0, *TICKS):
                wave_err = (
                    4 * tick * float(truncation + rounding + UNDERFLOW_PAD) + 1e-14
                )
                # All norms below refer to the reference-purified two-input launch.
                obs_err = 2 * (wave_err + proj_err) + (wave_err + proj_err) ** 2 + 1e-11
                yes = frame @ (frame.T @ wave)
                no = wave - yes
                local = wave[mask]
                positive = yes[mask]
                negative_band = no[mask]
                after = vf.T @ local.reshape(-1, 2)
                p_a, p_p, p_ap, p_pa = map(weight, (local, yes, after, positive))
                p_an = p_a - p_ap
                assert -obs_err < p_an < 1 + obs_err
                assert p_ap <= bound * p_a + obs_err
                # Distinguishing the two low lines leaves a mixture of two
                # rank-one maps from the reference into the receiver.
                overlaps = frame.T @ wave
                lines = np.concatenate(
                    [(frame[mask, i, None] * overlaps[i]) for i in range(2)], axis=0
                )
                reads = {
                    name: pauli_record(kraus, obs_err)
                    for name, kraus in (
                        ("local_only", local),
                        ("energy_yes_then_local", positive),
                        ("energy_no_then_local", negative_band),
                        ("fine_low_line_then_local", lines),
                    )
                }
                assert (
                    reads["fine_low_line_then_local"][
                        "partial_transpose_minimum_interval"
                    ][1]
                    >= 0
                )
                ap_order = [1 - p_a, p_an, p_ap]
                pa_order = [
                    p_p - p_pa,
                    p_pa,
                    1 - p_p - weight(negative_band),
                    weight(negative_band),
                ]
                for outcomes in (ap_order, pa_order):
                    assert min(outcomes) >= -obs_err
                    np.testing.assert_allclose(sum(outcomes), 1, atol=1e-12)
                records.append(
                    {
                        "flight_bridge_rate": rate,
                        "wait": float(tick * STEP),
                        "local_first_outcomes": np.maximum(ap_order, 0).tolist(),
                        "energy_first_outcomes": np.maximum(pa_order, 0).tolist(),
                        "later_band_pass_given_local": None
                        if p_a == 0
                        else interval(p_ap / p_a, 2 * obs_err / (p_a - obs_err)),
                        "absolute_error_single_pass_or_reference_moment": obs_err,
                        "reads": reads,
                    }
                )
            if tick < TICKS[-1]:
                wave = taylor_step(h, wave)
    return {
        "every_locally_selected_face_maximum_band_acceptance_upper": float(bound),
        "operator_error_finite_clock": eta,
        "finite_clock_reset_failure_probability_upper": 2 * eta,
        "flights": records,
    }


def run():
    with archived_fixture() as fixture:
        data = case(FanPath(fixture, 4, 1, True), FanPath(fixture, 4, -1, True))
    return {
        "protocol_sha256": sha256((HERE / "ENERGY_ORDER.md").read_bytes()).hexdigest(),
        "local_first_outcome_order": [
            "local_refused",
            "local_then_energy_no",
            "local_then_energy_yes",
        ],
        "energy_first_outcome_order": [
            "energy_yes_local_refused",
            "energy_yes_local_pass",
            "energy_no_local_refused",
            "energy_no_local_pass",
        ],
        "ideal_clock_records_with_finite_error": data,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parents[1] / "replay.py"))["compare"]
        compare(
            json.loads(json.dumps(data, allow_nan=False)),
            json.loads((HERE / "energy_order.json").read_text()),
            "energy_order",
        )
        print(
            "Reproduced both orders of the local charge–face and whole-guide energy reads."
        )


if __name__ == "__main__":
    main()
