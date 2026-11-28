"""Keep the two archived preparations and remove one tip reference in flight.

Run with the OPH checkout containing the three frozen inputs as the positional
argument. The frozen author headline may be next to this file, in its ordinary
checkout path, or supplied by --headline. This control uses no OPH modules.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import isclose
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.linalg import eigh


HERE = Path(__file__).resolve().parent
HEADLINE = "openwave/interval_echo/row_04/spectroscopy_release/headline.py"
HEADLINE_SHA256 = "b39cc85d5198989dd45c5b05e516a5144b832caf5ee76d551200bf36b404415d"
TIP = 9
APERTURE = (9, 11)
WAITS = (0, 1, 7, 37, 150)
GATE = 0.98
ALLOWANCE = 2e-9
RECORDS = ("signed_weight", "raw_signed_click", "two_tip_occupancy")


def digests(root, author):
    for relative, expected in author["DIGESTS"].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input digest mismatch for {relative}")


def native_and_controlled(raw, author):
    paths, positions, laplacian, native = author["physical_histories"](raw)
    physical_degree = raw["degree"]
    assert tuple(int(physical_degree[p]) for p in APERTURE) == (5, 5)
    tip = (positions == TIP).astype(int)
    assert np.count_nonzero(tip) == 6
    # A blank zero-load reference contributes κ n_p². There is exactly one
    # such reference at a degree-five port. M is expressed in units κ/2.
    altered = native - 2 * np.diag(tip)
    np.testing.assert_array_equal(np.diag(native - altered), 2 * tip)
    selected = author["select"](paths, laplacian, native)
    geometry = []
    for column, ray in selected:
        indices = np.flatnonzero(ray)
        assert tuple(sorted(positions[indices])) == APERTURE
        direction = np.sign(ray).astype(int)
        twist = altered @ direction
        assert int(direction @ twist) == 0 and int(twist @ twist) == 2
        # On the ideal signed ray, the exact variance of the altered M/2
        # is 1/4. This obstruction concerns the coherent law itself.
        geometry.append(
            {"column": column, "two_tip_histories": [paths[i] for i in indices]}
        )
    return paths, positions, selected, native, altered, geometry


def evolve(mode, scale, vector):
    # At zero duration both arms use the same supplied vector exactly,
    # without a distinct round trip through approximate bases.
    if scale == 0:
        return vector
    levels, basis = mode
    return basis @ (np.exp(-1j * scale * levels) * (basis.T @ vector))


def recorded(vector, ray, aperture_mask, commission):
    assert abs(float(np.vdot(vector, vector).real) - 1) < ALLOWANCE
    weight = commission * abs(ray @ vector) ** 2
    pocket = commission * float(np.sum(abs(vector[aperture_mask]) ** 2))
    assert 0 <= weight <= pocket + ALLOWANCE and pocket <= commission + ALLOWANCE
    return [round(float(weight), 12), round(float(weight / 8), 12), round(pocket, 12)]


def record(root, headline):
    if not __debug__:
        raise RuntimeError("Replay requires Python without -O")
    headline_sha256 = hashlib.sha256(headline.read_bytes()).hexdigest()
    if headline_sha256 != HEADLINE_SHA256:
        raise ValueError("The author headline does not match the archived SHA-256")
    author = run_path(str(headline))
    digests(root, author)
    boards = root / author["BOARD"]
    with np.load(boards, allow_pickle=False) as board:
        paths, positions, selected, native, altered, geometry = native_and_controlled(
            board, author
        )
        stations = len(board["edges"]) + 1
    assert stations == 631 and len(paths) == 32
    spectra = (eigh(native), eigh(altered))
    mask = np.isin(positions, APERTURE)
    clock = json.loads((root / author["CLOCK"]).read_bytes())
    words = json.loads((root / author["WORDS"]).read_bytes())
    commission_exact = (1 - Fraction(1, 2**96)) ** stations
    commission = float(commission_exact)
    source_refusal = float(1 - commission_exact)
    arms = []
    for bank in clock["runs"]:
        strength = bank["strength_in_simulator"]
        calibrated = author["read_clock"](
            bank["counts"],
            clock["appointments"],
            clock["shots_per_appointment"],
            len(clock["runs"]),
        )
        row = {"strength_in_simulator": strength, "original_clock_record": calibrated}
        arms.append(row)
        if calibrated is None:
            assert strength == 0
            row["disposition"] = "inconclusive_clock"
            continue
        estimate, lower, upper = calibrated
        assert lower <= abs(strength) <= upper
        ratio = strength / estimate
        pulse = evolve_matrix(spectra[0], 4 * ratio * words["step"])
        row["disposition"] = "calibrated"
        row["launches"] = []
        for column, ray in selected:
            source = np.eye(len(paths), dtype=complex)[0]
            for angles in words["pulses"][column]["angles"]:
                source = pulse @ source
                for port, angle in zip(words["ports"], angles, strict=True):
                    source[positions == port] *= np.exp(-1j * ratio * angle)
            observed = []
            for time in WAITS:
                phase = 0.5 * strength * time
                base, cut = (
                    recorded(evolve(mode, phase, source), ray, mask, commission)
                    for mode in spectra
                )
                if time == 0:
                    assert base == cut
                if base[0] - ALLOWANCE <= GATE or base[2] - ALLOWANCE <= GATE:
                    raise ArithmeticError(
                        "The native reference did not retain the excitation"
                    )
                observed.append(
                    {"exposure": time, "native": base, "reference_9_removed": cut}
                )
            row["launches"].append({"column": column, "outcomes": observed})

    calibrated_rows = [row for row in arms if row["disposition"] == "calibrated"]
    unresolved = [
        row["strength_in_simulator"] for row in arms if row not in calibrated_rows
    ]
    runs = [run for row in calibrated_rows for run in row["launches"]]
    assert len(runs) == 32 and unresolved == [0]
    contrasts = {}
    for time in WAITS:
        outcomes = [run["outcomes"][WAITS.index(time)] for run in runs]
        contrasts[str(time)] = {}
        for j, reading in ((0, "signed_weight"), (2, "two_tip_occupancy")):
            altered_weights = [entry["reference_9_removed"][j] for entry in outcomes]
            contrasts[str(time)][reading] = {
                "observed_minimum": min(altered_weights),
                "observed_maximum": max(altered_weights),
                "launches_rejecting_0p98": sum(
                    weight + ALLOWANCE < GATE for weight in altered_weights
                ),
            }
    if (
        any(
            contrasts[str(time)]["signed_weight"]["launches_rejecting_0p98"]
            != len(runs)
            for time in WAITS[1:]
        )
        or contrasts["0"]["signed_weight"]["launches_rejecting_0p98"]
    ):
        raise ArithmeticError(
            "The coherent intervention failed its distinction control"
        )
    return {
        "scope": "target-aware coherent flight-only control on the archived preparation",
        "headline_sha256": headline_sha256,
        "physical_tip": TIP,
        "physical_degrees_of_two_tips": [5, 5],
        "physical_reference_count_removed": 1,
        "squared_contrast_term": "kappa * n_9**2 with a blank reference at load zero",
        "eligible_history_count": 6,
        "eligible_histories": [path for path in paths if path[-1] == TIP],
        "original_detector_and_prefix": geometry,
        "ideal_prepared_ray_energy_variance_over_kappa_squared": [1, 4],
        "resolved_releases": len(runs),
        "unresolved_strengths": unresolved,
        "records_per_arm": RECORDS,
        "preparation_exhaustion_per_original": source_refusal,
        "failure_gate": GATE,
        "probability_allowance": ALLOWANCE,
        "summary_by_exposure": contrasts,
        "full_release_records": arms,
    }


def evolve_matrix(mode, scale):
    levels, basis = mode
    return (basis * np.exp(-1j * scale * levels)) @ basis.T


def check_receipt(actual, saved, field="receipt"):
    if type(actual) is not type(saved):
        raise ArithmeticError(f"Different types at {field}")
    if isinstance(actual, dict) and actual.keys() == saved.keys():
        for key in actual:
            check_receipt(actual[key], saved[key], f"{field}.{key}")
        return
    if isinstance(actual, list) and len(actual) == len(saved):
        for i, (read, expected) in enumerate(zip(actual, saved, strict=True)):
            check_receipt(read, expected, f"{field}[{i}]")
        return
    if isinstance(actual, float) and isclose(
        actual, saved, abs_tol=2e-10, rel_tol=2e-13
    ):
        return
    if actual != saved:
        raise ArithmeticError(
            f"The coherent reference-control replay differs at {field}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("oph_root", type=Path)
    parser.add_argument(
        "--headline", type=Path, help="Path to the frozen author headline"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--record", type=Path, help="Create a new evidence file")
    group.add_argument("--check", type=Path, help="Verify a saved evidence file")
    args = parser.parse_args()
    headline = args.headline
    if headline is None:
        adjacent = HERE / "headline.py"
        headline = adjacent if adjacent.is_file() else args.oph_root / HEADLINE
    data = json.loads(json.dumps(record(args.oph_root, headline), allow_nan=False))
    if args.record:
        with args.record.open("x") as output:
            output.write(json.dumps(data, separators=(",", ":")) + "\n")
    if args.check:
        check_receipt(data, json.loads(args.check.read_text()))
    print(json.dumps(data["summary_by_exposure"], indent=2))


if __name__ == "__main__":
    main()
