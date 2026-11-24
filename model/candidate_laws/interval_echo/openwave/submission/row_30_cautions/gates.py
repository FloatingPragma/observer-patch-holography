"""Exercise the proposed OW-30 gates on exposed material, without a public draw."""

import argparse
import json
from math import sqrt
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from ...row_30_collective import experiment, holdout
from ...row_30_collective.registration import verify as registration


def local_bank():
    # Supply only the published pilot selector. This never calls settings(),
    # official_anchor() or the guarded command for a prospective acquisition.
    with TemporaryDirectory() as temporary, patch.object(
        registration, "settings", return_value=[holdout.SEED] * 3
    ):
        release = Path(temporary) / "audit"
        passed = registration.run_bank("exposed-fixture", release)
        result = json.loads((release / "result.json").read_text())
        boards = []
        for board in result["boards"]:
            level = board["level"]
            summary = {
                "level": level,
                "failed_gates": sorted(k for k, ok in board["gates"].items() if not ok),
                "count_arrays": len(list((release / str(level)).glob("*.npy"))),
            }
            if "score" in board:
                summary["lattice_max"] = board["score"]["lattice_max"]
                summary["far_count_fraction"] = board["score"]["far_count_difference"]
            if "error" in board:
                summary["acquisition_error"] = board["error"]
            boards.append(summary)
    return {"bank_return": passed, "archived_passed": result["passed"], "boards": boards}


def changed_flight():
    native_spectrum = experiment.dual
    native_spectroscopy = experiment.spectroscopy
    flight_started = False

    def start_flight(*args, **kwargs):
        nonlocal flight_started
        flight_started = True
        return native_spectroscopy(*args, **kwargs)

    def flight_spectrum(length, coupling):
        # Keep the meter, commissioning, exact preflight checks and zero arm.
        # The first spectroscopy call starts the altered flight for both
        # spectroscopy and transport, including raw terminal-count replay.
        if flight_started and coupling == experiment.K:
            coupling = 7.5
        return native_spectrum(length, coupling)

    with patch.object(experiment, "spectroscopy", start_flight), patch.object(
        experiment, "dual", flight_spectrum
    ):
        return local_bank()


def ordinary_wave_number():
    archive = json.loads((experiment.HERE / "holdout.json").read_text())
    m, c2 = archive["m"], archive["c_squared"]
    rows = []
    for board in archive["boards"]:
        pairs = [
            (r["k"], r["omega"])
            for r in board["observer"]["modes"]
            if 0 < r["k"] <= 0.50
        ]
        square_root = max(abs(sqrt(m * m + c2 * k * k) - w) for k, w in pairs)
        tangent = max(abs(m + c2 * k * k / (2 * m) - w) for k, w in pairs)
        span = max(k for k, _ in pairs) - min(k for k, _ in pairs)
        coverage = len(pairs) >= 2 and span >= 0.02
        rows.append(
            {
                "level": board["level"],
                "lines": len(pairs),
                "k_span": span,
                "square_root_max": square_root,
                "tangent_max": tangent,
                "square_root_passes": coverage and square_root < 0.0048,
                "tangent_passes": coverage and tangent < 0.0048,
            }
        )
    return rows


def audit():
    with patch.object(
        experiment, "spectroscopy", side_effect=ArithmeticError("injected acquisition outage")
    ):
        missing = local_bank()
    return {
        "scope": "Local gate audit on previously exposed material",
        "material_seed_on_each_level": holdout.SEED,
        "frozen_manifest_sha256": registration.file_hash(registration.HERE / "files.json"),
        "verifier_sha256": registration.file_hash(Path(registration.__file__)),
        "calibrated_native_coupling": experiment.K,
        "perturbed_flight_coupling": 7.5,
        "perturbed_flight": changed_flight(),
        "missing_acquisition": missing,
        "ordinary_k_tolerance": 0.0048,
        "ordinary_k_models": ordinary_wave_number(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write a separate audit receipt")
    target = parser.parse_args().output
    receipt = audit()
    if target is None:
        saved = Path(__file__).with_suffix(".json")
        experiment.verify(receipt, json.loads(saved.read_text()))
        print("Reproduced the gate audit on exposed material.")
    else:
        target.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
