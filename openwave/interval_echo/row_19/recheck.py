"""Recompute the finite collar example, the five-station Gram rank and the disk holdout."""

import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = (
    ROOT / "code/matched_einstein_limit/source_implication/full_geometry/results.json"
)
DISK = ROOT / "code/cosmology/rar_deep_regime/cut_spacing/disk_domain.json"
TRANSFERS = (
    ROOT / "code/matched_einstein_limit/source_implication/signed_transfers.json"
)


def run():
    # FiniteCapGeneratorSplit.exampleState: two bulk sectors and one collar.
    weights = np.array([1.0, 1.0, math.exp(-1)])
    modular = -np.log(weights)
    anomalous = modular - np.average(modular[:2], weights=weights[:2])
    charge = np.average(anomalous, weights=weights)
    np.testing.assert_allclose(anomalous, [0, 0, 1], atol=1e-15, rtol=0)
    np.testing.assert_allclose(charge, 1 / (2 * math.e + 1), atol=1e-15, rtol=0)

    transfers = json.loads(TRANSFERS.read_bytes())
    if transfers["order"] != [
        "survey_area",
        "survey_volume",
        "full_head_entropy",
        "parked_energy",
    ]:
        raise ValueError("Signed transfer observables changed order")
    matched = {t["source"]: np.array(t["matched"]) for t in transfers["tests"]}
    ordinary, collar = matched["ordinary"], matched["collar"]
    np.testing.assert_allclose([ordinary[1], collar[1], collar[3]], 0, atol=1e-12)
    eta = -ordinary[2] / ordinary[0]
    prediction = -collar[2] / eta
    np.testing.assert_allclose(prediction, collar[0], atol=1e-12)
    geometry = json.loads(SOURCE.read_bytes())
    dimensions = {}
    for name in ("original_four", "extended_five"):
        block = geometry[name]
        d = np.array(block["distances"])
        center = np.eye(len(d)) - np.ones_like(d) / len(d)
        gram = -center @ d**2 @ center / 2
        positive = np.linalg.eigvalsh(gram)[::-1][:-1]
        np.testing.assert_allclose(
            positive, block["centered_eigenvalues_descending"], atol=1e-10, rtol=0
        )
        dimensions[name] = dict(
            positive_eigenvalues=positive.tolist(),
            flat_embedding_rank=int(sum(positive > 1e-8)),
        )
    outer = json.loads(DISK.read_bytes())
    ring = outer["outer"]["figure8"]
    return dict(
        finite_example_collar_expectation=float(charge),
        four_station_ordinary_ratio=float(eta),
        predicted_collar_area=float(prediction),
        observed_collar_area=float(collar[0]),
        change_in_parked_energy=float(collar[3]),
        original_and_extended_ranging=dimensions,
        agc114905_radius_kpc=ring["radius_kpc"],
        measured_kms=ring["circular_kms"],
        minimum_disk_model_kms=ring["min_phase_model_kms"],
        minimum_excess_speed_kms=ring["min_phase_model_kms"] - ring["circular_kms"],
        enclosed_gas_fraction=outer["gas_mass_fraction_in_outer_cylinder"],
        input_sha256={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (SOURCE, DISK, TRANSFERS)
        },
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
