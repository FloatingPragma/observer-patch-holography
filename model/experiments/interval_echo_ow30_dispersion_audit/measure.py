"""Compare two dispersion forecasts with counted sign and clock records."""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_collective import (
    experiment as strand,
    formation,
)
from model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation import (
    acquire,
    spatial,
)


HERE = Path(__file__).resolve().parent
LIMIT = 0.0048
VALIDATION = strand.HERE.parent / "submission/row_30_confirmation/validation.json"


class Counts:
    def __init__(self, directory):
        self.directory = directory

    def keep(self, key, clicks):
        np.save(self.directory / key, np.asarray(clicks, dtype="<u8"))


def forecast(wave_numbers, meter):
    z = meter["c_squared"] * (2 - 2 * np.cos(wave_numbers)) / meter["m"] ** 2
    root = meter["m"] * np.sqrt(1 + z)
    hopping = meter["m"] * (1 + z / 2 - z**2 / 8)
    return root, hopping


def transfer(n, coupling, times, profiles):
    u, energies, _ = strand.dual(n, coupling)
    phases = np.cos(times[:, None] * energies)
    return [(phases * (u.T @ profile)) @ u.T for profile in profiles]


def board_run(plan, meter, directory):
    n = len(plan["route"])
    m = n + 1
    channels = min(8, max(4, (m + 1) // 4))
    labels = np.arange(1, channels + 1)
    launch = np.pi * (m + 1 - labels) / (m + 1)
    profiles = np.sqrt(2 / (m + 1)) * np.sin(np.outer(launch, np.arange(1, m + 1)))
    times = strand.DT * np.arange(strand.SAMPLES)
    rng = np.random.default_rng(plan["seed"])
    counts = Counts(directory)
    clock = strand.read_clock(times, rng, counts)
    step = clock["dt_observed"]
    observed_times = step * np.arange(strand.SAMPLES)
    window = np.hanning(strand.SAMPLES)
    native = transfer(n, strand.K, times, profiles)
    off = transfer(n, 0.0, times, profiles)

    def amplitude(signal, frequency):
        return float(
            2
            * abs(np.dot(window * signal, np.exp(-1j * frequency * observed_times)))
            / window.sum()
        )

    lines, controls = {}, []
    for first in (1, 2):
        for second in range(first + 1, channels + 1):
            fa, fb = profiles[first - 1], profiles[second - 1]
            pa, pb = native[first - 1], native[second - 1]
            ca, cb = off[first - 1], off[second - 1]
            norm = strand.wedge_norm(fa, fb)
            aa, ab, bb = pa @ fa, pa @ fb, pb @ fb
            key = f"p{first}_{second}_"
            signal = strand.counts(aa * bb - ab**2, norm**2, rng, counts, key + "sum")
            changed = strand.counts(
                (ca @ fa) * (cb @ fb) - (ca @ fb) ** 2,
                norm**2,
                rng,
                counts,
                key + "changed",
            )
            frequency = strand.projection(signal, None, step)
            lines[f"{first},{second}"] = frequency
            controls.append(
                [
                    amplitude(signal, frequency),
                    amplitude(changed, frequency),
                    amplitude(changed, 2),
                ]
            )
            for i, (ref, response) in enumerate(
                (
                    (fa, aa[:, None] * pb - ab[:, None] * pa),
                    (fb, ab[:, None] * pb - bb[:, None] * pa),
                )
            ):
                strand.counts(
                    response,
                    norm * (abs(ref).sum() - abs(ref)),
                    rng,
                    counts,
                    key + f"scan{i}",
                )

    # Reflected sine profiles have the same absolute weights as the low launches
    # used by the frozen count decoder. Its demodulation and sign normalization
    # therefore apply directly. The source-only foil below uses the high labels.
    observer = {
        "pair_sum_lines": lines,
        "modes": [{"spatial_launch": int(i), "k": 0} for i in labels],
    }
    reading = spatial.score({"dimers": n, "observer": observer}, directory)
    if "even_trained_k" not in reading:
        return {
            "seed": plan["seed"],
            "dimers": n,
            "status": reading["status"],
            "reading": reading,
        }
    foil = spatial.score(
        {"dimers": n, "observer": observer}, directory, candidate_k=launch
    )
    k = np.array(reading["even_trained_k"])
    omega = np.array(reading["omega_from_pair_counts"])
    root, hopping = forecast(k, meter)
    closure = max(
        abs(lines[f"2,{j}"] - omega[1] - omega[j - 1]) for j in range(4, channels + 1)
    )
    cutoff = reading["quadratic_cutoff"]
    contrast, quiet, calibrated = np.array(controls).T
    resolved = bool(
        reading["even_trained_q"] <= cutoff < foil["candidate_q"]
        and closure <= LIMIT
        and min(contrast) > 0.4
        and max(quiet) < 0.02
        and min(calibrated) > 0.4
    )
    result = {
        "level": plan["level"],
        "seed": plan["seed"],
        "dimers": n,
        "status": "resolved" if resolved else "unresolved",
        "clock_step": step,
        "pair_sum_closure": float(closure),
        "min_contrast": float(min(contrast)),
        "max_changed_contrast": float(max(quiet)),
        "min_changed_at_two": float(min(calibrated)),
        "odd_profile_q": reading["even_trained_q"],
        "launch_label_q": foil["candidate_q"],
        "profile_limit": cutoff,
        "pair_sums": lines,
        "modes": [
            {
                "launch": int(m + 1 - i),
                "k": float(p),
                "omega": float(w),
                "root_residual": float(w - r),
                "hopping_residual": float(w - h),
            }
            for i, p, w, r, h in zip(labels, k, omega, root, hopping)
        ],
    }
    return result


def run(folder):
    published = json.loads(VALIDATION.read_text())
    meter = acquire.measured_meter()
    sources = {
        "seam_columns": strand.native_checks(),
        "commissioning_column_error": formation.native(),
        "dense_strands": [strand.exact_check(n) for n in (5, 6)],
        "boards": [],
    }
    if sources["commissioning_column_error"] > 1e-11:
        raise ValueError("commissioning columns differ from the native binary sector")
    boards = []
    for saved in published["boards"]:
        plan = saved["enrollment"]
        if acquire.enrollment(plan["level"], plan["seed"]) != plan:
            raise ValueError(
                "enrolled geometry differs from archived validation geometry"
            )
        setup = formation.prepare(plan["level"], plan["seed"])
        if setup["dimers_enrolled"] != len(plan["route"]):
            raise ValueError("prepared and measured strands have different lengths")
        sources["boards"].append(setup)
        target = folder / str(plan["seed"])
        target.mkdir()
        boards.append(board_run(plan, meter, target))
    comparison = []
    for saved in published["boards"]:
        modes = saved["dispersion"]["modes"]
        k = np.array([line["k"] for line in modes])
        w = np.array([line["omega"] for line in modes])
        root, hopping = forecast(k, meter)
        comparison.append(
            {
                "seed": saved["enrollment"]["seed"],
                "n": len(k),
                "max_root_error": float(max(abs(w - root))),
                "max_hopping_error": float(max(abs(w - hopping))),
            }
        )
    return {
        "validation_sha256": hashlib.sha256(VALIDATION.read_bytes()).hexdigest(),
        "source_checks": sources,
        "meter": meter,
        "limit": LIMIT,
        "published_frequencies": comparison,
        "high_band_boards": boards,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--retain", type=Path, help="new directory for every count array"
    )
    parser.add_argument(
        "--record",
        type=Path,
        help="write a new result (default compares with results.json)",
    )
    args = parser.parse_args()
    with TemporaryDirectory() as tmp:
        directory = Path(tmp)
        if args.retain:
            args.retain.mkdir()
            directory = args.retain
        result = run(directory)
    if args.record:
        with args.record.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        strand.verify(result, json.loads((HERE / "results.json").read_text()))
    for board in result["high_band_boards"]:
        errors = (
            {
                key: max(abs(line[key]) for line in board.get("modes", []))
                for key in ("root_residual", "hopping_residual")
            }
            if board.get("modes")
            else {}
        )
        print(board["seed"], board["status"], errors)
    if any(b["status"] != "resolved" for b in result["high_band_boards"]):
        raise SystemExit(3)
    if any(
        abs(line["root_residual"]) > LIMIT or abs(line["hopping_residual"]) <= LIMIT
        for b in result["high_band_boards"]
        for line in b["modes"]
    ):
        raise SystemExit(2)
