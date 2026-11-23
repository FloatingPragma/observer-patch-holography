"""Test curvature inferred from native route volumes against unused physical range records."""

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "proper_ball_ruler" / "coarse.json"


def profile(x):
    """Twice the integral of sin(t)^2 on [0,x]; Taylor form avoids cancellation at zero."""
    if abs(x) < 0.01:
        return (2 / 3) * x**3 * (1 - x * x / 5 + 2 * x**4 / 105 - x**6 / 945)
    return x - np.sin(2 * x) / 2


def doubling(x):
    return 8.0 if x == 0 else profile(2 * x) / profile(x)


def infer_angle(ratio):
    if not 2 <= ratio <= 8:
        raise ValueError("Volume ratio outside the positive-curvature pre-cut domain")
    return brentq(lambda x: doubling(x) - ratio, 0, np.pi / 2, xtol=1e-14)


def gram_check(distance, frequency, landmarks):
    angle = frequency * distance
    gram = np.cos(angle)
    best = dict(relative_error_lower_bound=0.0)
    indices = range(len(distance))
    for size in range(5, len(distance) + 1):
        for subset in combinations(indices, size):
            sl = np.ix_(subset, subset)
            spectrum = np.linalg.eigvalsh(gram[sl])
            norm_bound = np.linalg.eigvalsh(angle[sl])[-1]
            necessary = max(-spectrum[0], spectrum[-5], 0) / norm_bound
            if necessary > best["relative_error_lower_bound"]:
                best = dict(
                    relative_error_lower_bound=float(necessary),
                    indices=list(subset),
                    addresses=[landmarks[i] for i in subset],
                    least_eigenvalue=float(spectrum[0]),
                    fifth_largest_eigenvalue=float(spectrum[-5]),
                    distance_angle_norm=float(norm_bound),
                )
    # Supplementary interval certificate for the same witness, chosen by the registered bound.
    if best.get("least_eigenvalue", 0) < -1e-12:
        selected = np.ix_(best["indices"], best["indices"])
        angles = angle[selected]
        _, vectors = np.linalg.eigh(gram[selected])
        vector = vectors[:, 0]
        weights = np.outer(vector, vector)

        def quadratic_form_ceiling(error):
            lower = np.maximum((1 - error) * angles, 0)
            upper = np.minimum((1 + error) * angles, np.pi)
            entries = np.where(weights >= 0, np.cos(lower), np.cos(upper))
            return np.sum(entries * weights)

        best["interval_witness_vector"] = vector.tolist()
        best["quadratic_form_ceiling_at_five_percent"] = quadratic_form_ceiling(0.05)
        best["interval_relative_error_lower_bound"] = brentq(
            quadratic_form_ceiling, 0, 1
        )
    return dict(
        full_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
        strongest_submatrix=best,
        excludes_five_percent_ranges=bool(best["relative_error_lower_bound"] > 0.05),
    )


def run():
    archive = json.loads(ARCHIVE.read_bytes())
    results = []
    for panel in archive["panels"]:
        radius = panel["microscopic_ball_radius"]
        relative = panel["relative_ball_radius"]
        n = panel["ball_count"]
        ratio = panel["doubled_ball_count"] / n
        x = infer_angle(ratio)
        releases = panel["launch_classes"]
        acquired_x = infer_angle((releases[0] + releases[1]) / releases[0])
        distance = np.asarray(panel["ball_landmark_ranges"], dtype=float)
        shell_fraction = 1 - profile(x * (radius - 1) / radius) / profile(x)
        flat_radius_angle = (1.5 * profile(x)) ** (1 / 3)
        results.append(
            dict(
                cells=panel["cells"],
                radius_swaps=radius,
                relative_radius=relative,
                home_range_swaps=panel["home_max_range"],
                volume_doubling_ratio=ratio,
                fitted_angle=x,
                independently_sampled_angle=acquired_x,
                curvature_in_home_range_units=(x / relative) ** 2,
                predicted_fixed_volume_area_over_flat=(np.sin(x) / flat_radius_angle)
                ** 2,
                predicted_fraction_on_last_shell=shell_fraction,
                observed_fraction_on_last_shell=panel["range_front"] / n,
                distance_test=gram_check(distance, x / radius, panel["ball_landmarks"]),
            )
        )
    k = results[1]["curvature_in_home_range_units"]
    forecast = []
    for panel in results[2:]:
        estimate = doubling(np.sqrt(k) * panel["relative_radius"])
        observed = panel["volume_doubling_ratio"]
        forecast.append(
            dict(
                cells=panel["cells"],
                predicted=estimate,
                observed=observed,
                relative_error=abs(estimate - observed) / observed,
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        native_archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
        panels=results,
        training_cells=results[1]["cells"],
        fixed_curvature_transfer=forecast,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    with parser.parse_args().output.open("x") as out:
        out.write(json.dumps(run(), indent=2, allow_nan=False) + "\n")
