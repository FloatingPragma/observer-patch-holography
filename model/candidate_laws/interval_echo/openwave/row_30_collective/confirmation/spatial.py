"""Cross-check inferred wave numbers against disjoint spatial detector trials.

Even holds supply training profiles. Odd holds at the same addresses test the
two-wave spatial span. Both splits use pair frequencies measured on separate
temporal counters. All hypotheses use the very same odd-hold signs and weights.
"""

from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from ..experiment import projection

SIGN_SHOTS = 2**32
CLOCK_SHOTS = 2**20
SAMPLES = 16384
TAIL_PROBABILITY = 1e-5


class PhysicalFailure(ValueError):
    """The physical records are present but do not resolve the enrolled waves."""


def _counts(directory, name, shape, shots):
    data = np.load(directory / (name + ".npy"), mmap_mode="r", allow_pickle=False)
    if data.shape != shape or data.dtype.kind not in "ui":
        raise ValueError(f"invalid shape or type for {name}")
    if np.any(data < 0) or np.any(data > shots):
        raise ValueError(f"count outside shot budget for {name}")
    return data


def _coordinates(board, directory):
    modes = board["observer"]["modes"]
    if "route" in board:
        n = len(board["route"])
    else:
        n = int(board["dimers"])
    m = n + 1
    ids = [entry["spatial_launch"] for entry in modes]
    channels = min(8, max(4, (m + 1) // 4))
    if ids != list(range(1, channels + 1)):
        raise ValueError("enrolled spatial launch labels are missing or reordered")
    pairs = [(0, b) for b in range(1, channels)] + [(1, b) for b in range(2, channels)]
    clock = [
        2.0 * _counts(directory, "clock_" + key, (SAMPLES,), CLOCK_SHOTS) / CLOCK_SHOTS
        - 1
        for key in ("cos", "sin")
    ]
    # Validate the entire spatial record inventory before returning a physical
    # failure. An early weak pair cannot conceal a later absent detector.
    scans = {
        (a, b, j): _counts(
            directory, f"p{a + 1}_{b + 1}_scan{j}", (SAMPLES, m), SIGN_SHOTS
        )
        for a, b in pairs
        for j in (0, 1)
    }
    sums = {
        (a, b): _counts(directory, f"p{a + 1}_{b + 1}_sum", (SAMPLES,), SIGN_SHOTS)
        for a, b in pairs
    }
    if m <= channels:
        raise PhysicalFailure("too few distinct receiver addresses")
    samples = np.arange(SAMPLES)
    phase = np.unwrap(np.arctan2(clock[1], clock[0]))
    step = np.polyfit(samples, phase, 1)[0]
    if not 0 < step < np.pi:
        raise PhysicalFailure("clock does not resolve increasing phase")
    launch_k = np.pi * np.arange(1, channels + 1) / (m + 1)
    templates = np.sqrt(2 / (m + 1)) * np.sin(np.outer(np.arange(1, m + 1), launch_k))
    projectors = [[] for _ in range(channels)]
    training_data = []
    tests = []
    pair_lines = {}
    for a, b in pairs:
        label = f"{a + 1},{b + 1}"
        temporal_sign = 2.0 * sums[a, b] / SIGN_SHOTS - 1
        if not np.ptp(temporal_sign):
            raise PhysicalFailure(f"no temporal signal for pair {label}")
        freq = projection(temporal_sign, None, step)
        if not np.isfinite(freq) or not 0 < freq < np.pi / step:
            raise PhysicalFailure(f"unresolved temporal pair {label}")
        if not np.isclose(
            freq, float(board["observer"]["pair_sum_lines"][label]), atol=1e-8, rtol=0
        ):
            raise ValueError(
                f"reported temporal line differs from counts for pair {label}"
            )
        pair_lines[a + 1, b + 1] = freq
        # A real demodulator fixes the quadrature independently of all scans.
        # Alternating the time indices retains the complete observation span.
        weight = np.hanning(SAMPLES) * np.cos(freq * step * samples)
        training = []
        for j, ref in enumerate((templates[:, a], templates[:, b])):
            # The absolute wedge weight at the source is a common scalar and
            # cancels in both QR and the noise-standardized spatial subspace fit.
            scale = np.sum(abs(ref)) - abs(ref)
            data = scans[a, b, j]
            for parity in (0, 1):
                selected = weight[parity::2]
                sign = 2.0 * data[parity::2] / SIGN_SHOTS - 1
                response = scale * (selected @ sign)
                sd = scale * np.linalg.norm(selected) / np.sqrt(SIGN_SHOTS)
                if parity == 0:
                    training.append(response)
                    training_data.append((a, b, response / sd, sd))
                else:
                    tests.append((a, b, response / sd, sd))
        basis = np.stack(training, axis=1)
        if np.linalg.matrix_rank(basis) < 2:
            raise PhysicalFailure(f"rank-deficient training scans on pair {label}")
        q = np.linalg.qr(basis)[0]
        p = q @ q.T
        projectors[a].append(p)
        projectors[b].append(p)
    trained = []
    for planes in projectors:
        _, vec = np.linalg.eigh(sum(planes))
        w = vec[:, -1]
        middle = w[1:-1]
        cosine = (middle @ (w[:-2] + w[2:])) / (2 * (middle @ middle))
        if not -1 < cosine < 1:
            raise PhysicalFailure(
                "training profiles do not yield an interior wave number"
            )
        trained.append(np.arccos(cosine))
    first = (pair_lines[1, 2] + pair_lines[1, 3] - pair_lines[2, 3]) / 2
    frequencies = [first] + [pair_lines[1, i] - first for i in range(2, channels + 1)]
    return m, np.asarray(trained), launch_k, training_data, tests, frequencies


def _residuals(m, wave_numbers, tests):
    # The open-chain receiver at address m is the fixed zero endpoint of the
    # standing-wave recurrence. With receiver x indexed from zero, each mode
    # has the shape sin((m-x)k); the two amplitudes are unrestricted nuisances.
    waves = np.sin(np.outer(np.arange(m, 0, -1), wave_numbers))
    residuals = []
    for a, b, data, sd in tests:
        design = waves[:, [a, b]] / sd[:, None]
        coeff, _, rank, _ = np.linalg.lstsq(design, data, rcond=None)
        if rank != 2:
            raise PhysicalFailure("candidate waves do not span a receiver pair")
        residuals.append(data - design @ coeff)
    return np.concatenate(residuals)


def _quadratic(m, wave_numbers, tests):
    residual = _residuals(m, wave_numbers, tests)
    return float(residual @ residual)


def score(board, directory, candidate_k=None):
    """Fit k on even holds, validate it on odd holds, and test any supplied foil."""
    result = dict(
        records_complete=False, profile_pass=False, launch_label_rejected=False
    )
    try:
        m, initial, launch, training, tests, frequencies = _coordinates(
            board, Path(directory)
        )
        # The interior recurrence is a consistent initializer; the final fit
        # uses every even-hold address to estimate k before consuming odd holds.
        fit = least_squares(
            lambda k: _residuals(m, k, training),
            initial,
            bounds=(0, np.pi),
            max_nfev=100,
        )
        if not fit.success:
            raise PhysicalFailure("even-hold spatial fit did not converge")
        trained = fit.x
        trained_q = _quadratic(m, trained, tests)
        if candidate_k is None:
            fitted_q = trained_q
        else:
            candidate = np.asarray(candidate_k, dtype=float)
            if candidate.shape != trained.shape:
                raise ValueError("wrong number of candidate wave numbers")
            if not np.all(np.isfinite(candidate)) or np.any(
                (candidate <= 0) | (candidate >= np.pi)
            ):
                raise PhysicalFailure(
                    "substitute includes an unresolved or endpoint wave number"
                )
            fitted_q = _quadratic(m, candidate, tests)
        # The full-record recurrence used in v2 is diagnostic. A degenerate
        # old reading cannot veto the prospectively specified v3 estimator.
        legacy_k = np.asarray(
            [row["k"] for row in board["observer"]["modes"]], dtype=float
        )
        legacy_q = None
        if np.all(np.isfinite(legacy_k)) and np.all(
            (legacy_k > 0) & (legacy_k < np.pi)
        ):
            try:
                legacy_q = _quadratic(m, legacy_k, tests)
            except PhysicalFailure:
                pass
        recurrence_q = _quadratic(m, initial, tests)
        launch_q = _quadratic(m, launch, tests)
        df = (m - 2) * len(tests)
    except FileNotFoundError as exc:
        return dict(result, status="missing_records", detail=Path(exc.filename).name)
    except PhysicalFailure as exc:
        return dict(
            result, records_complete=True, status="physical_fail", detail=str(exc)
        )
    except (ValueError, OSError, KeyError, TypeError, EOFError) as exc:
        return dict(result, status="invalid_records", detail=str(exc))
    # Each normalized sign sum has sub-Gaussian variance proxy <=1 by
    # Hoeffding. After an orthogonal least-squares projection of rank df the
    # usual quadratic-form upper tail is df + 2 sqrt(df t) + 2t. Calibration
    # of this bound applies to a fixed, correct spatial mean. Training noise,
    # finite-window leakage or inadequate signal can cause physical failure.
    t = -np.log(TAIL_PROBABILITY)
    cutoff = float(df + 2 * np.sqrt(df * t) + 2 * t)
    pass_profile = bool(fitted_q <= cutoff and trained_q <= cutoff)
    reject_label = bool(launch_q > cutoff)
    return dict(
        records_complete=True,
        status="pass" if pass_profile and reject_label else "physical_fail",
        profile_pass=pass_profile,
        launch_label_rejected=reject_label,
        projected_degrees_of_freedom=df,
        quadratic_cutoff=cutoff,
        candidate_q=fitted_q,
        even_trained_q=trained_q,
        even_recurrence_only_q=recurrence_q,
        full_recurrence_diagnostic_q=legacy_q,
        omega_from_pair_counts=frequencies,
        launch_label_q=launch_q,
        even_trained_k=trained.tolist(),
    )
