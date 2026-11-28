"""Compare outside-head frequencies and material-arrival counters on fresh copies."""

from itertools import product
from math import log, sqrt

import numpy as np

from ..calibration import reported_interval
from ..local_read import PAULI, partial_transpose, reconstructed
from ......fan_recoil.self_start.measure import DRAW_PAD, PROBABILITY_TICKS, bound


CURRENT_SHOTS = 2**42
WITNESS_SHOTS = 2**45
SIGNS = np.array(list(product((1, -1), repeat=2)))


def sample(probabilities, n, rng):
    p = np.asarray(probabilities)
    assert abs(p.sum() - 1) < 2e-13 and p.min() > -2e-14
    ticks = np.rint(p * PROBABILITY_TICKS).astype(np.int64)
    ticks[ticks.argmax()] += PROBABILITY_TICKS - ticks.sum()
    assert ticks.min() >= 0 and abs(ticks / PROBABILITY_TICKS - p).max() < DRAW_PAD
    return rng.multinomial(n, ticks / PROBABILITY_TICKS)


def mean(effect, rho):
    return float(np.trace(effect @ rho).real)


def refusal_bins(incoming, passed, collector):
    if collector is None:
        return ["energy_refused", "capture_refused"], [incoming - passed, 1 - incoming]
    pole, returned = collector
    return (
        [
            "wheel_energy_refused",
            "collector_refused",
            "rest_filter_refused",
            "capture_refused",
        ],
        [
            returned * (incoming - passed),
            incoming * (pole - returned),
            incoming * (1 - pole),
            1 - incoming,
        ],
    )


def current_bank(
    bank, embed, plus, minus, local, pad, rng, *, collector=None, shots=CURRENT_SHOTS
):
    sub = {s: embed @ r @ embed.conj().T for s, r in bank.items()}
    returned = 1 if collector is None else collector[1]
    labels, refusal = refusal_bins(
        bank["unread"].trace().real, sub["unread"].trace().real, collector
    )
    outside, inside = [], []
    for key in ("plus", "minus"):
        rho = sub[key]
        p = rho.trace().real
        a, b = (mean(e, rho) for e in (plus, minus))
        outside.extend([a / 2, (p - a) / 2, b / 2, (p - b) / 2])
        src, dest = (mean(e, rho) for e in local)
        inside.extend([src, dest, p - src - dest])
    ec = sample([*(returned * np.array(outside)), *refusal], shots, rng)
    mc = sample([*(returned * np.array(inside)), *refusal], shots, rng)
    rows = {}
    for i, key in enumerate(("plus", "minus", "unread")):
        a, b = ec[:8].reshape(2, 4), mc[:6].reshape(2, 3)
        ext = a[i] if i < 2 else a.sum(axis=0)
        loc = b[i] if i < 2 else b.sum(axis=0)
        rho = sub[key]
        passed = float(rho.trace().real)
        # The archived two-word uncertainty bounds each Pauli moment. Both
        # preparations and each unresolved projection are contractions.
        delta = 2 * pad + 2e-10
        prob_pad = (2 * delta + 12 * DRAW_PAD / returned) / (passed - delta)
        u, v = bound(ext[0] + ext[3], ext.sum(), prob_pad)
        l0, h0 = bound(loc[0], loc.sum(), prob_pad)
        l1, h1 = bound(loc[1], loc.sum(), prob_pad)
        target = 4 * mean(plus - minus, rho) / passed
        np.testing.assert_allclose(
            target, 4 * mean(local[1] - local[0], rho) / passed, atol=4e-13, rtol=0
        )
        outside_interval = reported_interval(4 * (2 * u - 1), 4 * (2 * v - 1))
        inside_interval = reported_interval(4 * (l1 - h0), 4 * (h1 - l0))
        assert max(outside_interval[0], inside_interval[0]) <= min(
            outside_interval[1], inside_interval[1]
        )
        rows[key] = {
            "source_gate_per_original": returned * passed,
            "short_long_outside_flip_given_gate": [
                mean(e, rho) / passed for e in (plus, minus)
            ],
            "agreement_probability_given_gate": (1 + target / 4) / 2,
            "outward_current_native": target,
            "outside_current_simultaneous_99pct": outside_interval,
            "independent_material_current_simultaneous_99pct": inside_interval,
        }
    return {
        "originals_per_apparatus": shots,
        "outside_bins_per_report": [
            "short_changed",
            "short_preserved",
            "long_changed",
            "long_preserved",
        ],
        "material_bins_per_report": ["source", "destination", "neither"],
        "bin_order": ["plus_report", "minus_report", *labels],
        "outside_counts": ec.tolist(),
        "independent_material_counts": mc.tolist(),
        "reports": rows,
    }


def frame(effects, v, err):
    eye = np.eye(v.shape[0])
    projected = {name: v.conj().T @ (2 * a - eye) @ v for name, a in effects.items()}
    # Bounds two columns, both occurrences of the angular frame and small matrix products.
    guard = 4 * sqrt(2) * err + 1e-11
    ops = [np.eye(2)]
    names = ["identity"]
    for name, op in projected.items():
        design = np.array([[mean(o, a) / 2 for a in PAULI] for o in [*ops, op]])
        # For at most four rows the perturbation has norm <= 2*guard.
        if len(ops) < 4 and np.linalg.svd(design, compute_uv=False)[-1] > 2 * guard:
            ops.append(op)
            names.append(name)
    assert len(ops) == 4
    design = np.array([[mean(o, a) / 2 for a in PAULI] for o in ops])
    singular = np.linalg.svd(design, compute_uv=False)
    return (
        ops,
        design,
        {
            "outside_settings": names,
            "calibration_operator_norm_guard": guard,
            "smallest_singular_lower": float(singular[-1] - 2 * guard),
            "condition_number_upper": float(
                (singular[0] + 2 * guard) / (singular[-1] - 2 * guard)
            ),
            "closed_aperture_span_dimension_including_identity": 1,
        },
    )


def fixed_witness(design, transfer, ops, guard):
    # This direction was published for the 1045 incoming word/reference pair.
    # Its transferred direction uses the analyzer calibration, before taking new counts.
    direction = np.linalg.solve(transfer.conj().T, [2, 1, 1, -1 - 1j])
    direction /= np.linalg.norm(direction)
    w = partial_transpose(np.outer(direction, direction.conj()))
    pauli = np.array([[mean(w, np.kron(a, b)) / 4 for b in PAULI] for a in PAULI])
    coeff = np.linalg.solve(design.T, pauli)
    assembled = sum(
        coeff[i, j] * np.kron(op, b)
        for i, op in enumerate(ops)
        for j, b in enumerate(PAULI)
    )
    assert np.linalg.norm(assembled - w, ord=2) < 2e-12
    return w, coeff, guard * abs(coeff).sum() + 2e-11


def witness_bank(
    t_input,
    transfer,
    ops,
    coeff,
    w,
    pad,
    op_error,
    randomize,
    rng,
    *,
    collector=None,
    shots=WITNESS_SHOTS,
    comparisons=3,
    receiver="outside",
):
    raw = reconstructed(t_input)
    joint = transfer @ raw @ transfer.conj().T
    if randomize:
        reduced = np.trace(joint.reshape(2, 2, 2, 2), axis1=1, axis2=3)
        joint = np.kron(reduced, np.eye(2) / 2)
    p = float(joint.trace().real)
    incoming = float(t_input[0, 0])
    returned = 1 if collector is None else collector[1]
    labels, refusal = refusal_bins(incoming, p, collector)
    if receiver == "face":
        assert collector is not None and abs(p - incoming) < 1e-14
        labels, refusal = labels[1:], refusal[1:]
    sample_score = coeff[0, 0]
    span = 0.0
    banks = []
    x, y = SIGNS.T
    for i, j in product((1, 2, 3), repeat=2):
        o, b = ops[i], PAULI[j]
        mx, my = mean(np.kron(o, np.eye(2)), joint), mean(np.kron(np.eye(2), b), joint)
        mxy = mean(np.kron(o, b), joint)
        sharp = (p + x * mx + y * my + x * y * mxy) / 4
        n = sample([*(returned * sharp), *refusal], shots, rng)
        scores = coeff[i, j] * x * y
        if i == j:
            scores += coeff[i, 0] * x + coeff[0, j] * y
        accepted = n[:4].sum()
        sample_score += n[:4] @ scores / accepted
        span += np.ptp(scores) ** 2 / accepted
        setting = (
            {"outside_setting": i}
            if receiver == "outside"
            else {"receiving_face_axis": "IXYZ"[i]}
        )
        banks.append({**setting, "reference_axis": "IXYZ"[j], "counts": n.tolist()})
    # Conditional on each setting's accepted count, its four joint outcomes
    # are independent bounded trials. Union over the enrolled comparisons.
    stat = sqrt(log(2 * comparisons / 0.01) * span / 2)
    # Sixteen input Pauli errors sum in trace norm. inputs()'s pad is four
    # times a single-moment bound plus 4e-10. Include frame and dyadic errors.
    delta = 4 * pad + 4e-10
    draw_pad = DRAW_PAD / returned
    model = abs(coeff).sum() * (2 * delta + 16 * draw_pad) / (p - delta - 4 * draw_pad)
    radius = stat + model + op_error
    return {
        "originals_per_setting": shots,
        "bins": ["++", "+-", "-+", "--", *labels],
        "source_gate_per_original": returned * p,
        "expected_witness": mean(w, joint) / p,
        "count_score": float(sample_score),
        "count_radius": stat,
        "model_and_dyadic_radius": float(model),
        "apparatus_witness_norm_guard": float(op_error),
        "simultaneous_99pct": reported_interval(
            sample_score - radius, sample_score + radius
        ),
        "settings": banks,
    }
