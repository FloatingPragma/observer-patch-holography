"""Read a localized vacancy's internal period and energy independently."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    anchored_vacancy as anchor,
)
from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    experiment as compact,
)
from model.candidate_laws.interval_echo.openwave.row_04_compact.gateway_acquisition import (
    prepare_layout,
)

from .sharp import unmonitored

SHOTS = 2**18
SEED = 178137
WAITS = np.linspace(0.125, 16.125, 129)
SEGMENTS = 8
HERE = Path(__file__).resolve().parent


def counts(rng, probabilities):
    values = np.round(probabilities, 12)
    if min(values) < -1e-10 or abs(sum(values) - 1) > 1e-9:
        raise ArithmeticError("Detector record has invalid normalization")
    values = np.clip(values, 0, 1)
    return rng.multinomial(SHOTS, values / sum(values)).tolist()


def fit_trace(rows):
    """Estimate the clock only from timer swaps and the two enrolled port counts."""
    timer, sample = np.array(rows)[:, 0], np.array(rows)[:, 1:3].sum(axis=1)
    times = 2 * SEGMENTS * np.arcsin(np.sqrt(timer / SHOTS))
    observed = sample / SHOTS

    def design(freq):
        return np.column_stack(
            (np.ones(len(times)), np.cos(freq * times), np.sin(freq * times))
        )

    def fit(freq):
        matrix = design(freq)
        coeff = np.linalg.lstsq(matrix[:65], observed[:65], rcond=None)[0]
        return matrix @ coeff, coeff

    def loss(freq):
        predicted, _ = fit(freq)
        return np.sum((predicted[:65] - observed[:65]) ** 2)

    grid = np.arange(0.15, 24.51, 0.05)  # below the nominal Nyquist frequency
    best = grid[np.argmin([loss(f) for f in grid])]
    freq = minimize_scalar(loss, bounds=(best - 0.05, best + 0.05), method="bounded").x
    forecast, coeff = fit(freq)
    amp = 2 * np.linalg.norm(coeff[1:])
    rms = np.sqrt(np.mean((forecast[65:] - observed[65:]) ** 2))
    return {
        "frequency_per_clock_gap_if_resolved": float(freq) if amp > 0.03 else None,
        "fitted_peak_to_trough": float(amp),
        "held_rms": float(rms),
        "held_rms_constant": float(
            np.sqrt(np.mean((observed[65:] - observed[:65].mean()) ** 2))
        ),
    }


def energy_meter(local, state, success, strength, rate, rng):
    # The controlled full exchange on an occupied seam reads the two swap parities.
    law = anchor.candidate.native_law
    indices = [law.SRC.index(p) for p in ((0, 1), (1, 0))]
    seam = anchor.free_seam(strength)[np.ix_(indices, indices)]
    swap = expm(-1j * np.pi * seam)
    minus = (np.eye(2) - swap) / 2
    np.testing.assert_allclose(minus, [[0.5, -0.5], [-0.5, 0.5]], atol=2e-13)
    anti = present = 0.0
    for edge in local.seams:
        pair = state[list(edge)]
        anti += np.linalg.norm(minus @ pair) ** 2
        present += np.linalg.norm(pair) ** 2
    n = len(local.seams)
    clicks = counts(
        rng, [anti / n, (present - anti) / n, success - present / n, 1 - success]
    )
    admitted = sum(clicks[:3])
    prob = clicks[0] / admitted
    kinetic = rate * n * prob
    return {
        "counts_antisymmetric_symmetric_elsewhere_source_failed": clicks,
        "exchange_estimate": float(kinetic),
        "exchange_standard_error": float(
            rate * n * np.sqrt(prob * (1 - prob) / admitted)
        ),
        "estimated_host_referenced_mean_energy": float(4 * strength + kinetic),
    }


def panel(level, axis, prep_pulse, strength, rate, rng):
    local, axes = compact.geometry()
    board = anchor.federation(level)
    charts = compact.gateway_table(board)
    home = next(
        (
            i
            for i, tips in enumerate(charts)
            if any(set(pair) <= set(tips) for pair in axes)
        ),
        0,
    )
    tips = charts[home]
    pair = [p for p in axes if not set(p).intersection(tips)][axis]
    prep = prepare_layout(tips, pulse=prep_pulse)
    launch = pair[0]
    source = prep.accepted_vectors[:, launch]
    success = float(np.vdot(source, source).real)
    if success < 0.01:
        raise ArithmeticError("Insufficient source acceptance")
    local_slice = slice(12 * home, 12 * (home + 1))
    vacancy = np.zeros(board.n, complex)
    vacancy[local_slice] = source
    h = anchor.defect_generator(board, -1, strength, 0)
    kinetic = anchor.defect_generator(board, -1, 0, 0)
    h = h + (rate - 1) * kinetic
    # [32,54] etc. can be centered arbitrarily for a single-history observation.
    centered = h - 4 * strength * sparse.eye(board.n, format="csr")
    waves = expm_multiply(
        -1j * centered, vacancy, start=WAITS[0], stop=WAITS[-1], num=len(WAITS)
    )
    outside = (abs(waves[:, : local_slice.start]) ** 2).sum(axis=1)
    outside += (abs(waves[:, local_slice.stop :]) ** 2).sum(axis=1)
    occupancies = abs(waves[:, local_slice]) ** 2
    probs = occupancies[:, list(pair)]
    rows = []
    dimer_ids = [anchor.candidate.native_law.SRC.index(p) for p in ((0, 1), (1, 0))]
    timer_h = anchor.free_seam(anchor.K)[np.ix_(dimer_ids, dimer_ids)]
    for t, (a, b) in zip(WAITS, probs, strict=True):
        swaps = abs(expm(-1j * timer_h * t / SEGMENTS)[1, 0]) ** 2
        dimer_swaps = rng.binomial(SHOTS, round(swaps, 12))
        rows.append(
            [int(dimer_swaps), *counts(rng, [a, b, success - a - b, 1 - success])]
        )
    bands = compact.dark_bands(compact.candidate.generator(local).toarray(), tips)
    pieces = [projector @ source for _, projector in bands]
    # Each piece is checked on the full connected vacancy Hamiltonian, including external seams.
    energies = []
    residuals = []
    for piece in pieces:
        v = np.zeros_like(vacancy)
        v[local_slice] = piece
        hv = h @ v
        e = (np.vdot(v, hv) / np.vdot(v, v)).real
        energies.append(float(e))
        residuals.append(float(np.linalg.norm(hv - e * v)))
    return {
        "cells": board.cells,
        "cell": home,
        "gateways": tips,
        "detector_ports": list(pair),
        "preparation_pulse": prep_pulse,
        "strength": strength,
        "exchange_rate": rate,
        "source_acceptance": success,
        "max_outside_probability": float(outside.max()),
        "compact_energies": energies,
        "maximum_eigen_residual": max(residuals),
        "source_norm_outside_compact_space": float(
            np.linalg.norm(source - sum(pieces))
        ),
        "measurements": rows,
        "frequency_fit": fit_trace(rows),
        "energy_read": energy_meter(local, source, success, strength, rate, rng),
    }


def spectral_bounds(board, sign):
    h = anchor.defect_generator(board, sign, anchor.K, 0)
    radius = np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal())
    lo, hi = float(min(h.diagonal() - radius)), float(max(h.diagonal() + radius))
    return {"energy_interval": [lo, hi], "largest_possible_internal_frequency": hi - lo}


def run():
    rng = np.random.default_rng(SEED)
    reference = json.loads(
        Path(anchor.__file__).with_name("material_reference.json").read_text()
    )
    transfers = []
    for cells, entry in zip((20, 320), reference["native"], strict=True):
        read = entry["panels"]["6"]["vacancy"]
        transfers.append(
            {
                "cells": cells,
                "frequency": read["frequency_per_clock_gap"],
                "standard_error": read["frequency_standard_error"],
            }
        )
    cases = {
        "no_antipodal_gateways": (0, 0, 1.0, 8, 1),
        "antipodal_80": (1, 0, 1.0, 8, 1),
        "antipodal_320": (2, 0, 1.0, 8, 1),
        "alternate_axis": (1, 1, 1.0, 8, 1),
        "shorter_preparation": (2, 0, 0.5, 8, 1),
        "half_contrast": (1, 0, 1.0, 4, 1),
        "exchange_off_at_release": (1, 0, 1.0, 8, 0),
    }
    return {
        "seed": SEED,
        "shots": SHOTS,
        "equal_timer_segments": SEGMENTS,
        "record_columns": [
            "timer swaps of enrolled timer copies",
            "launch port",
            "opposite port",
            "another port",
            "source failure",
        ],
        "first_wait": float(WAITS[0]),
        "last_wait": float(WAITS[-1]),
        "number_of_waits": len(WAITS),
        "independent_stationary_vacancy_six_seam_transfers": transfers,
        "single_defect_bounds": [
            {
                "cells": (board := anchor.federation(level)).cells,
                "vacancy": spectral_bounds(board, -1),
                "excess": spectral_bounds(board, 1),
            }
            for level in range(3)
        ],
        "panels": {name: panel(*case, rng) for name, case in cases.items()},
        "unmonitored_sharp_returns": [
            unmonitored(level, SHOTS, rng) for level in (1, 2)
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        compare(result, json.loads((HERE / "results.json").read_text()), "row_03_clock")
        print(
            "Reproduced seven free-vacancy panels and the full one-defect clock bounds."
        )


if __name__ == "__main__":
    main()
