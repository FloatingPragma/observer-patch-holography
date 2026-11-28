"""Find and certify rare band directions with a remote ring at receiver capture."""

from fractions import Fraction
import json
from math import inf, nextafter

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

from ..fan_path.channel_bound import capture_floor, dyadic_vectors, sqrt_bounds
from ..fan_path.fixture import FanPath, archived_fixture
from ..fan_path.spectral import separated_reference, verified_inertia


def probability_bounds(squared, error):
    amplitude_low, amplitude_high = sqrt_bounds(squared)
    return max(Fraction(0), amplitude_low - error) ** 2, (amplitude_high + error) ** 2


def detector_grams(path, trial):
    """Exact terminal capture and terminal capture with a source-ring excitation."""
    receiver = path.blocks[path.size - 1, 0]
    assert receiver.shape[-1] == 4
    code = (2 * path.code).astype(int).astype(object)
    np.testing.assert_array_equal(code, 2 * path.code)
    np.testing.assert_array_equal(code.T @ code, 4 * np.eye(2, dtype=int))
    # The physical terminal hub is the last four-word factor. Dot products
    # use Python integers. The terminal code contributes the common factor 1/4.
    words = trial[receiver.rows].reshape(-1, 4, 2).astype(object)
    codes = (code.T @ words).reshape(*receiver.shape[:-1], 2, 2)
    flat = codes.reshape(-1, 2)
    capture = flat.T @ flat
    # At the greatest distance from the last hub only the first fan counts.
    # Rejecting its normalized uniform mode is a reader on a separate factor.
    source_vacuum = codes.sum(axis=0).reshape(-1, 2)
    remote = capture - Fraction(1, receiver.shape[0]) * (
        source_vacuum.T @ source_vacuum
    )
    return capture, remote


def run_case(path):
    _, bracket, _ = separated_reference(path)
    assert [verified_inertia(path.h, cut)["count_below"] for cut in bracket] == [2, 2]
    start = np.random.default_rng(923).standard_normal(path.h.shape[0])
    values, vectors = eigsh(path.h, k=4, which="SA", v0=start, tol=1e-13)
    order = values.argsort()
    trial, norms, _, _, error = dyadic_vectors(
        path.h, values[order[:2]], vectors[:, order[:2]], bracket[0]
    )
    source_admission, _ = capture_floor(path, 0, trial, norms, error)
    assert source_admission > 0
    capture, remote = detector_grams(path, trial)
    _, directions = eigh(remote.astype(float), capture.astype(float))
    direction = directions[:, -1]
    # The generalized eigensolve only proposes a direction. Verify the rounded
    # coefficient pair by exact comparisons with the certified frame error.
    direction /= direction[np.argmax(abs(direction))]
    q = np.rint(direction * 2**26).astype(int).astype(object)
    norm = sum(q[i] ** 2 * norms[i] for i in range(2))
    assert norm > 0
    denominator = 4 * norm
    capture_lower, capture_upper = probability_bounds(
        Fraction(q @ capture @ q, denominator), error
    )
    remote_lower, _ = probability_bounds(Fraction(q @ remote @ q, denominator), error)
    assert capture_lower > 0
    conditional = remote_lower / capture_upper
    threshold = {3: Fraction(994, 1000), 4: Fraction(984, 1000)}
    assert conditional > threshold[path.size]
    return {
        "fans": path.size,
        "source_to_terminal_bridges": [list(edge) for edge in path.bridges],
        "charge_sign": path.sign,
        "source_ring_distance": path.size - 1,
        "integer_trial_direction": q.tolist(),
        "joint_terminal_capture_given_band_probability_bounds": [
            nextafter(float(capture_lower), -inf),
            nextafter(float(capture_upper), inf),
        ],
        "source_ring_positive_given_band_and_terminal_capture_lower": nextafter(
            float(conditional), -inf
        ),
        "source_to_band_smallest_singular_value_lower": nextafter(
            float(source_admission), -inf
        ),
    }


def run():
    with archived_fixture() as archive:
        return [run_case(FanPath(archive, n, 1, asymmetric=True)) for n in (3, 4)]


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
