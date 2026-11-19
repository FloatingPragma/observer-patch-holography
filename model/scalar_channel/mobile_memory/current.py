"""Signed seam-current records and record-controlled pair conversion."""

from fractions import Fraction
from math import isfinite

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .law import EXTERNAL_WEIGHT, MASS, pair_generator


CURRENT_RECORDS = (-2, -1, 0, 1, 2)
TERMINAL_RECORDS = (
    "charged_pair_at_target",
    "other_charged_pair",
    "neutral_pair_at_target",
    "other_neutral_pair",
)
PREPARATIONS = (-1, 1)

# The columns are the +1 and -1 oriented number-current modes for a
# difference-square hop. A positive seam weight multiplies both eigenvalues.
_CURRENT_MODES = np.array([[1, 1], [1j, -1j]]) / np.sqrt(2)

# These are the actual local propagators, including the phase relative to
# untouched ports. Forward and reverse use positive pulses of relative
# duration pi/4 and 3*pi/4 under the same difference square.
_FORWARD = (1 - 1j) / 2 * np.array([[1, 1j], [1j, 1]])
_REVERSE = _FORWARD.conj().T


def _ports(nports, seam):
    first, second = seam
    if not 0 <= first < nports or not 0 <= second < nports or first == second:
        raise ValueError("current read needs two distinct ports on the carrier")
    return [first, second]


def one_charge_read(state, seam, charge):
    """Return the unnormalized projections labeled by signed current -1, 0, +1."""
    state = np.asarray(state, dtype=complex)
    if state.ndim != 1 or charge not in (-1, 1):
        raise ValueError("current read needs one charged port vector and charge +1 or -1")
    ports = _ports(len(state), seam)
    local = state[ports]
    elsewhere = state.copy()
    elsewhere[ports] = 0
    records = {0: elsewhere}
    for direction, mode in zip((1, -1), _CURRENT_MODES.T, strict=True):
        branch = np.zeros_like(state)
        branch[ports] = mode * np.vdot(mode, local)
        records[charge * direction] = branch
    return {record: records[record] for record in (-1, 0, 1)}


def _rotate_axes(amplitudes, ports, positive, negative):
    """Apply two two-port operators to the positive and negative charge coordinates."""
    result = amplitudes.copy()
    result[ports, :] = positive @ result[ports, :]
    result[:, ports] = result[:, ports] @ negative.T
    return result


def pair_current_read(wave, nports, seam):
    """Read the signed current while preserving coherence within each outcome."""
    ports = _ports(nports, seam)
    charged, neutral = np.asarray(wave, dtype=complex).reshape(2, nports, nports)
    basis_change = _CURRENT_MODES.conj().T
    read_basis = _rotate_axes(charged, ports, basis_change, basis_change)
    directions = np.zeros(nports, dtype=int)
    directions[ports] = (1, -1)
    signed_current = directions[:, None] - directions[None, :]
    branches = {}
    for record in CURRENT_RECORDS:
        state = np.zeros((2, nports, nports), dtype=complex)
        state[0] = _rotate_axes(
            np.where(signed_current == record, read_basis, 0),
            ports,
            _CURRENT_MODES,
            _CURRENT_MODES,
        )
        if record == 0:
            state[1] = neutral
        branches[record] = state.ravel()
    return branches


def _steer_pair(wave, nports, ports, sign):
    result = wave.copy().reshape(2, nports, nports)
    if sign:
        left, right = (_REVERSE, _FORWARD) if sign > 0 else (_FORWARD, _REVERSE)
        result[0] = _rotate_axes(result[0], ports, left, right)
    return result.ravel()


def _terminal_record(wave, nports, target):
    charged, neutral = np.abs(wave.reshape(2, nports, nports)) ** 2
    charged_here = charged[target, target]
    neutral_here = neutral[target, target]
    charged[target, target] = 0
    neutral[target, target] = 0
    return np.array([charged_here, charged.sum(), neutral_here, neutral.sum()])


def current_bound(wait=Fraction(1, 128), spectral_bound=28):
    """Exact lower bound on feedback minus fixed-arm neutral detection."""
    wait = Fraction(wait)
    spectral_bound = Fraction(spectral_bound)
    if wait <= 0 or spectral_bound <= 0:
        raise ValueError("the wait and spectral bound must be positive")
    remainder = (spectral_bound * wait) ** 3 / 6
    # The target amplitude has imaginary part t + O(t^3). The other
    # endpoint has zero constant, first and second coefficients.
    at_target_lower = max(Fraction(0), wait - remainder) ** 2
    from_other_upper = remainder**2
    return (at_target_lower - from_other_upper) / 2


def current_experiment(car, seam, wait=1 / 128, external_weight=EXTERNAL_WEIGHT):
    """Retain full joint records for two balanced phase preparations and three controls."""
    nports = car.n
    ports = _ports(nports, seam)
    if tuple(sorted(ports)) not in car.seams:
        raise ValueError("the current detector must use a carrier seam")
    if not (isfinite(wait) and wait > 0 and isfinite(external_weight) and 0 < external_weight <= 1):
        raise ValueError("current experiment needs a finite positive wait and weight in (0, 1]")

    prepared = {}
    for direction in PREPARATIONS:
        mode = (1 - direction) // 2
        initial = np.zeros((2, nports, nports), dtype=complex)
        initial[0][np.ix_(ports, ports)] = np.outer(
            _CURRENT_MODES[:, mode], _CURRENT_MODES[:, 1 - mode]
        )
        prepared[direction] = pair_current_read(initial.ravel(), nports, seam)

    # Both phases and every nonzero branch use the same source and wait. Reuse
    # the two pulse arms in the independent fair-arm mixture and fixed-arm control.
    indexed_branches = {}
    for direction, records in prepared.items():
        for record, branch in records.items():
            if np.any(branch):
                actions = (-1, 1) if record else (0,)
                for action in actions:
                    indexed_branches[direction, record, action] = _steer_pair(
                        branch, nports, ports, action
                    )
    generator = pair_generator(car, external_weight)
    step = -1j * wait * generator
    evolved = expm_multiply(
        step,
        np.column_stack(list(indexed_branches.values())),
        traceA=step.diagonal().sum(),
    )
    terminals = {
        key: _terminal_record(state, nports, ports[0])
        for key, state in zip(indexed_branches, evolved.T, strict=True)
    }
    zero = np.zeros(len(TERMINAL_RECORDS))
    results = {name: [] for name in ("feedback", "fixed", "shuffled")}
    for direction in PREPARATIONS:
        feedback, fixed, shuffled = [], [], []
        for record in CURRENT_RECORDS:
            if record == 0:
                row = terminals.get((direction, record, 0), zero)
                feedback.append(row)
                fixed.append(row)
                shuffled.append(row)
                continue
            forward = terminals.get((direction, record, 1), zero)
            reverse = terminals.get((direction, record, -1), zero)
            feedback.append(forward if record > 0 else reverse)
            fixed.append(forward)
            shuffled.append((forward + reverse) / 2)
        results["feedback"].append(feedback)
        results["fixed"].append(fixed)
        results["shuffled"].append(shuffled)

    joint = {name: np.array(rows) for name, rows in results.items()}
    marginals = {
        name: float(rows[:, :, 2].sum() / len(PREPARATIONS))
        for name, rows in joint.items()
    }
    spectral_bound = 4 * int(max(car.degree)) + 2 * MASS + 2
    return {
        "oriented_seam": list(ports),
        "wait_in_inverse_J": float(wait),
        "external_weight": float(external_weight),
        "phase_preparations": list(PREPARATIONS),
        "current_records": list(CURRENT_RECORDS),
        "terminal_records": list(TERMINAL_RECORDS),
        "joint_probabilities": {name: rows.tolist() for name, rows in joint.items()},
        "neutral_pair_at_target": marginals,
        "feedback_advantage": marginals["feedback"] - marginals["fixed"],
        "source_norm_upper_bound": spectral_bound,
        "feedback_advantage_lower_bound": str(current_bound(wait, spectral_bound)),
    }
