"""A first-read mass tag on a shared, possibly interacting, finite source."""

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


TOLERANCE = 2e-10
SIGNS = (-1, 1)
ANALYZERS = (1, -1j)


def _propagate(source, initial, wait):
    step = -1j * wait * source
    return expm_multiply(step, initial, traceA=step.diagonal().sum())


def _record(minus, plus, groups):
    """Retain both reads and every terminal position group for each analyzer."""
    labels, bins = np.unique(groups, return_inverse=True)
    probabilities = np.empty((len(ANALYZERS), len(SIGNS), len(SIGNS), len(labels)))
    for axis, analyzer in enumerate(ANALYZERS):
        for first_index, first in enumerate(SIGNS):
            for final_index, final in enumerate(SIGNS):
                # The first X read acts on a definite tag. The second
                # probe contributes the remaining factor 1/sqrt(2).
                occupation = abs(minus + first * final * analyzer * plus) ** 2 / 8
                probabilities[axis, first_index, final_index] = np.bincount(
                    bins, weights=occupation, minlength=len(labels),
                )
    np.testing.assert_allclose(probabilities.sum(axis=(2, 3)), 0.5,
                               atol=TOLERANCE, rtol=0)
    parity = np.outer(SIGNS, SIGNS)
    real, imag = (probabilities.sum(axis=-1) * parity).sum(axis=(1, 2))
    return labels, probabilities, complex(real, imag)


def _grade_integral(source, grade, initial, wait):
    # The lower block solves i dχ/dt = Hχ + βψ from χ(0) = 0.
    tangent = sparse.bmat(((source, None), (grade, source)), format="csr")
    wave, derivative = np.split(_propagate(tangent, np.r_[initial, np.zeros_like(initial)], wait), 2)
    return float(-np.vdot(wave, derivative).imag)


def audit(H, beta, initial, wait, *, gap=1 / 16, groups=None):
    """Read two mass arms and certify their inclusive clock response.

    The same H acts in both arms. Only the signed mass term differs.
    The terminal groups report the full joint outcomes, while clock
    estimation sums over them before forming the phase.
    """
    source = sparse.csr_matrix(H)
    if source.shape[0] == 0 or source.shape[0] != source.shape[1]:
        raise ValueError("the clock needs a nonempty square generator")
    skew = source - source.getH()
    if not np.all(np.isfinite(source.data)) or np.max(abs(skew.data), initial=0) > TOLERANCE:
        raise ValueError("the clock generator must be finite and Hermitian")
    size = source.shape[0]
    diagonal = np.asarray(beta)
    if diagonal.shape != (size,) or not np.all(np.isin(diagonal, (-1, 1))):
        raise ValueError("the mass grading must give one sign per basis state")
    initial = np.asarray(initial, dtype=complex)
    if initial.shape != (size,) or not np.isclose(np.linalg.norm(initial), 1.0, atol=TOLERANCE, rtol=0):
        raise ValueError("the clock needs a normalized source state")
    if not math.isfinite(wait) or not math.isfinite(gap) or wait <= 0 or gap <= 0 or gap * wait > 0.5:
        raise ValueError("the clock certificate requires 0 < gap * wait <= 1/2")
    groups = np.zeros(size, dtype=int) if groups is None else np.asarray(groups)
    if groups.shape != (size,) or not np.issubdtype(groups.dtype, np.integer):
        raise ValueError("the terminal read must give one integer group per basis state")

    grade = sparse.diags(diagonal.real.astype(float), format="csr")
    minus = _propagate(source - (gap / 2) * grade, initial, wait)
    plus = _propagate(source + (gap / 2) * grade, initial, wait)
    labels, probabilities, overlap = _record(minus, plus, groups)
    np.testing.assert_allclose(overlap, np.vdot(minus, plus), atol=TOLERANCE, rtol=0)

    integral = _grade_integral(source, grade, initial, wait)
    x = gap * wait
    phase = -math.atan2(overlap.imag, overlap.real)
    # The Dyson expansion gives a real-parameter derivative bound wait**n.
    # Conjugacy under gap reversal supplies the even/odd Taylor bounds.
    if not 1 - x ** 2 / 2 - TOLERANCE <= overlap.real <= 1 + TOLERANCE:
        raise ValueError("the clock overlap lies outside the phase certificate")
    if abs(overlap.imag + gap * integral) > x ** 3 / 6 + TOLERANCE:
        raise ValueError("the clock sine response exceeds the centered remainder")
    if abs(phase - gap * integral) > 2 * x ** 3 + TOLERANCE:
        raise ValueError("the clock phase exceeds the centered remainder")
    return {
        "wait": float(wait),
        "mass_gap": float(gap),
        "analyzers": ["X", "Y"],
        "first_and_last_read_signs": list(SIGNS),
        "terminal_groups": labels.tolist(),
        "joint_records": probabilities.tolist(),
        "ramsey_overlap": [overlap.real, overlap.imag],
        "visibility": abs(overlap),
        "phase": phase,
        "phase_clock_rate": phase / x,
        "sine_clock_time": -overlap.imag / gap,
        "exact_integrated_grade": integral,
        "exact_grade_rate": integral / wait,
        "normalized_phase_error_bound": 2 * x ** 2,
        "sine_clock_time_error_bound": gap ** 2 * wait ** 3 / 6,
    }
