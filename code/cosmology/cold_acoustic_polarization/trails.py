"""Test all paid, retained native trail stocks against the same Planck bands."""

import hashlib
import json
from fractions import Fraction

import numpy as np

from .boltzmann import BARYONS, RATIOS, compute
from .chains import CHANNELS, weighted_samples
from .measure import HERE, Likelihood, norm

ROOT = HERE.parents[2]
RECEIPT = ROOT / (
    "model/experiments/interval_echo_shell_flux/product_address/green_response/"
    "material_recoil/overlap_source/native_boundary/stock/trail.json"
)


def prepare(board):
    initial = Fraction(board["initial_energy"])
    final = Fraction(board["simultaneous_terminal_energy"])
    echo = Fraction(board["work_from_addressed_echoes"])
    read = Fraction(board["work_from_sharp_reads"])
    work = Fraction(board["net_work"])
    if work != final - initial or work != echo + read:
        raise ValueError("Mismatch between the paid work and conserved board excess")
    if Fraction(board["terminal_matched_blank"]) != initial:
        raise ValueError("Initial and terminal blank controls have different energies")
    return {
        "owners": board["owners"],
        "live_seams": board["live_seams"],
        "mean_read_echo_appointments": board["supplied_mean_echoes"],
        "initial_one_load": str(initial),
        "terminal_total": str(final),
        "work_from_addressed_echoes": str(echo),
        "work_from_sharp_reads": str(read),
        "net_work": str(work),
        "ratio_exact": str(work / initial),
        "ratio": float(work / initial),
        "local_passive_excess_floor": board["passive_excess_vs_any_one_load_blank"],
    }


def rank(ratio, samples):
    q = np.quantile(samples, [0.025, 0.16, 0.84, 0.975])
    return {
        "fraction_of_weight_at_or_below": float(np.mean(samples <= ratio)),
        "within_central_68": bool(q[1] <= ratio <= q[2]),
        "within_central_95": bool(q[0] <= ratio <= q[3]),
    }


def run():
    import camb

    raw = RECEIPT.read_bytes()
    native = json.loads(raw)
    boards = [prepare(board) for board in native["boards"]]
    with np.load(HERE / "chain_density_columns.npz") as posteriors:
        for channel in CHANNELS:
            samples = weighted_samples(
                np.vstack([posteriors[f"{channel}_{i}"] for i in range(1, 5)])
            )
            for board in boards:
                board.setdefault("published_chain", {})[channel] = rank(
                    board["ratio"], samples
                )

    lk = Likelihood()
    control = json.loads((HERE / "results.json").read_bytes())["models"]
    fac = np.tile(lk.eff * (lk.eff + 1) / (2 * np.pi), 3)
    for baryon in BARYONS:
        standard = control[str(baryon)]["standard_control"]["amplitude_scale_tilt"]
        standard_theory, _ = compute(baryon, RATIOS["standard_control"])
        standard_mean = lk.transformed(
            standard_theory, standard["temperature_calibration"]
        )
        standard_score = lk.check_polarization(standard_mean)["conditional_P"]
        np.testing.assert_allclose(
            standard_score, standard["scores"]["conditional_P"], rtol=3e-5
        )
        for board in boards:
            theory, zeq = compute(baryon, board["ratio"])
            main = lk.evaluate(theory, True)
            fitted = lk.transformed(theory, main["temperature_calibration"])
            diff = fitted - standard_mean
            separation = norm(lk.cond, diff[lk.n :] - lk.B @ diff[: lk.n])
            board.setdefault("baryon_density_runs", {})[str(baryon)] = {
                "z_equality": zeq,
                "amplitude_scale_tilt": main,
                "amplitude_scale": lk.evaluate(theory, False),
                "broad_tilt_sensitivity": lk.evaluate(theory, True, 1.0),
                "conditional_P_excess_score_over_standard": main["scores"][
                    "conditional_P"
                ]
                - standard_score,
                "conditional_prediction_squared_separation_from_standard": separation,
            }
            if baryon == BARYONS[0]:
                board["predicted_bands_uK2"] = (fitted * fac).reshape(3, -1).tolist()

    return {
        "CAMB": camb.__version__,
        "native_receipt": str(RECEIPT.relative_to(ROOT)),
        "native_receipt_sha256": hashlib.sha256(raw).hexdigest(),
        "native_protocol_sha256": native["protocol_sha256"],
        "TT_training_bins": lk.n,
        "held_polarization_bins": 2 * lk.n,
        "boards": boards,
    }


if __name__ == "__main__":
    (HERE / "trail_results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
