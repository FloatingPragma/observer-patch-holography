"""Transfer the radial boundary law from NGC 3198 to eight other galaxies."""

import hashlib
import json

import numpy as np

from .measure import run as primary
from .source import HERE, enroll, profile

CALIBRATOR = "NGC3198"
OTHER_PAIR_MEMBER = "NGC3521"


def run():
    kept, _ = enroll()
    targets = primary()["targets"]
    held = sorted(set(kept) - {CALIBRATOR, OTHER_PAIR_MEMBER})
    output = {}
    for name, options in (
        ("stars_atomic_molecular", {}),
        ("atomic_molecular", {"star": 0}),
        ("atomic", {"star": 0, "molecular": 0}),
    ):
        src = {key: profile(*item[:4], **options) for key, item in kept.items()}
        for law, use_finite in (("compact", False), ("finite", True)):
            x = {
                key: s["S_msun_kpc"] * (s["F"] if use_finite else 1)
                for key, s in src.items()
            }
            gain = targets[CALIBRATOR]["excess_v2"] / x[CALIBRATOR]
            rows = {}
            for key in held:
                target = targets[key]
                excess = gain * x[key]
                speed = np.sqrt(excess + target["baryon_v2"])
                rows[key] = {
                    "number_measured_CO_rings": src[key]["n_measured_star_HI_H2"][2],
                    "source": x[key],
                    "measured_excess_v2": target["excess_v2"],
                    "predicted_excess_v2": excess,
                    "predicted_over_measured_excess": excess / target["excess_v2"]
                    if target["excess_v2"]
                    else None,
                    "excess_residual_v2": excess - target["excess_v2"],
                    "observed_kms": target["observed_kms"],
                    "predicted_total_kms": speed,
                    "total_residual_kms": speed - target["observed_kms"],
                }
            result = {
                "gain": gain,
                "calibration_source": x[CALIBRATOR],
                "calibration_excess_v2": targets[CALIBRATOR]["excess_v2"],
                "held": rows,
            }
            for metric in ("total_residual_kms", "excess_residual_v2"):
                errors = np.array([row[metric] for row in rows.values()])
                result[metric] = {
                    "median_abs": np.median(abs(errors)),
                    "rms": np.sqrt(np.mean(errors**2)),
                    "worst_to_best": sorted(
                        rows, key=lambda key: abs(rows[key][metric]), reverse=True
                    ),
                }
            output[f"{name}_{law}"] = result
    return {
        "protocol_sha256": hashlib.sha256(
            (HERE / "GAS_TRANSFER_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "calibration": CALIBRATOR,
        "outside_motivating_pair": held,
        "predictions": output,
    }


if __name__ == "__main__":
    (HERE / "gas_transfer.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
