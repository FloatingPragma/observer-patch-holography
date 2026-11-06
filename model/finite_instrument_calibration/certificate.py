"""Exact finite probes and reference-uniform instrument error bounds.

The quantum channel model and trusted preparation/readback are hypotheses.
Examples below are algebraic controls, not observed laboratory records.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from math import isqrt
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = (
    "Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean",
    "Lean/EventAlgebra/LuedersPhaseInstrument.lean",
    "Lean/EventAlgebra/FiniteProbeInstrumentSelection.lean",
    "model/finite_instrument_calibration/certificate.py",
)


def exact_fraction(value):
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("exact rational inputs are required")
    if not isinstance(value, (int, str, F, s.Rational)):
        raise TypeError("exact rational inputs are required")
    return F(str(value))


def sqrt_upper(value, bits=64):
    """Outward rational square root; exact squares stay exact."""
    x = exact_fraction(value)
    if x < 0:
        raise ValueError("a nonnegative radicand is required")
    if isinstance(bits, bool) or not isinstance(bits, int) or bits < 1:
        raise ValueError("positive integer precision is required")
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    if a*a == x.numerator and b*b == x.denominator:
        return F(a, b)
    scale = 1 << bits
    k = isqrt((x.numerator*scale*scale)//x.denominator)
    if F(k*k, scale*scale) < x:
        k += 1
    return F(k, scale)


def blocks(ranks):
    ranks = tuple(ranks)
    if not ranks or any(isinstance(d, bool) or not isinstance(d, int) or d < 1
                        for d in ranks):
        raise ValueError("positive integer block ranks are required")
    start, out = 0, []
    for d in ranks:
        out.append(tuple(range(start, start+d)))
        start += d
    return tuple(out)


def calibration_bound(ranks, basis_errors, coherent_errors):
    """Half-diamond bound from joint correct-label/survival probabilities.

    coherent_errors has one entry for each block of rank greater than one.
    Upper bounds are unconditional; failures and wrong labels count as errors.
    """
    groups = blocks(ranks)
    errors = tuple(exact_fraction(x) for x in basis_errors)
    coherent = tuple(exact_fraction(x) for x in coherent_errors)
    if len(errors) != sum(map(len, groups)):
        raise ValueError("one basis error per input dimension is required")
    if len(coherent) != sum(len(g) > 1 for g in groups):
        raise ValueError("one coherent error per nontrivial block is required")
    if any(not 0 <= x <= 1 for x in errors+coherent):
        raise ValueError("probability errors must lie in [0,1]")
    spread_sq = F(0)
    j = 0
    for g in groups:
        if len(g) == 1:
            continue
        block_error = sum((errors[i] for i in g), F(0))
        spread_sq += len(g)*(sqrt_upper(coherent[j])+sqrt_upper(block_error))**2
        j += 1
    eta = sqrt_upper(sum(errors, F(0)))+sqrt_upper(spread_sq)
    return min(F(1), 2*eta)


def history_bound(instrument_error, uses, preparation_error=0):
    delta, epsilon = map(exact_fraction, (instrument_error, preparation_error))
    if not 0 <= delta <= 1 or not 0 <= epsilon <= 1:
        raise ValueError("half trace/distance bounds must lie in [0,1]")
    if isinstance(uses, bool) or not isinstance(uses, int) or uses < 0:
        raise ValueError("a nonnegative integer use count is required")
    return min(F(1), epsilon+uses*delta)


def exact_matrix(matrix):
    a = s.Matrix(matrix)
    if any(x.atoms(s.Float) or x.free_symbols for x in a):
        raise TypeError("fully specified exact matrix entries are required")
    return a.applyfunc(s.simplify)


def zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)


def projectors(ranks, basis=None):
    groups = blocks(ranks)
    n = sum(map(len, groups))
    u = s.eye(n) if basis is None else exact_matrix(basis)
    if u.shape != (n, n) or not zero(u.H*u-s.eye(n)):
        raise ValueError("an exact orthonormal basis is required")
    return tuple(sum((u[:, i]*u[:, i].H for i in g), s.zeros(n)) for g in groups)


def check_instrument(kraus, ranks):
    groups = blocks(ranks)
    n = sum(map(len, groups))
    if len(kraus) != len(groups):
        raise ValueError("all labeled outcome branches are required")
    out = tuple(tuple(exact_matrix(k) for k in branch) for branch in kraus)
    if any(k.shape != (n, n) for branch in out for k in branch):
        raise ValueError("all Kraus maps must have the declared input/output size")
    effect = sum((k.H*k for branch in out for k in branch), s.zeros(n))
    if not zero(effect-s.eye(n)):
        raise ValueError("the complete instrument must be trace preserving")
    return out


def probe_scores(kraus, ranks, basis=None):
    """Exact joint probability of the correct label and pure-state readback."""
    groups = blocks(ranks)
    n = sum(map(len, groups))
    projectors(ranks, basis)  # Validate the preparation/readback basis.
    u = s.eye(n) if basis is None else exact_matrix(basis)
    operators = check_instrument(kraus, ranks)

    def score(label, vector):
        return s.simplify(sum((s.conjugate((vector.H*k*vector)[0]) *
                              (vector.H*k*vector)[0]
                              for k in operators[label]), s.Integer(0)))

    basis_scores, coherent_scores = [], []
    for label, group in enumerate(groups):
        basis_scores.extend(score(label, u[:, i]) for i in group)
        if len(group) > 1:
            v = sum((u[:, i] for i in group), s.zeros(n, 1))/s.sqrt(len(group))
            coherent_scores.append(score(label, v))
    return tuple(basis_scores), tuple(coherent_scores)


def outcome_maps_equal(left, right, ranks):
    """Independent exact equality check on all matrix units and labels."""
    a, b = check_instrument(left, ranks), check_instrument(right, ranks)
    n = sum(ranks)
    for i in range(n):
        for j in range(n):
            unit = s.zeros(n)
            unit[i, j] = 1
            for ax, bx in zip(a, b):
                diff = sum((k*unit*k.H for k in ax), s.zeros(n))
                diff -= sum((k*unit*k.H for k in bx), s.zeros(n))
                if not zero(diff):
                    return False
    return True


def controls():
    ranks = (2, 1)
    p = projectors(ranks)
    ideal = tuple((x,) for x in p)
    e = tuple(s.diag(*[int(i == j) for i in range(3)]) for j in range(3))
    dephased = ((e[0], e[1]), (e[2],))
    swap = s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    twisted = tuple((swap*x,) for x in p)
    noisy = tuple((s.Rational(3, 5)*x, s.Rational(4, 5)*swap*x) for x in p)
    # The phase context is the committed two-dimensional Pauli-Y context.
    y_basis = s.Matrix([[1, 1], [s.I, -s.I]])/s.sqrt(2)
    yp = projectors((1, 1), y_basis)
    y_ideal = tuple((x,) for x in yp)
    x_swap = s.Matrix([[0, 1], [1, 0]])
    y_twisted = tuple((x_swap*x,) for x in yp)
    t = s.Rational(1, 1024)
    a, b = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    cross_rotation = s.Matrix([[a, 0, -b], [0, 1, 0], [b, 0, a]])
    nonprojective = tuple((s.Rational(3, 5)*x,
                          s.Rational(4, 5)*x*cross_rotation) for x in p)
    return (
        ("lueders_rank_two_plus_one", ideal, ranks, None),
        ("dephasing_passes_basis_but_fails_coherent_probe", dephased, ranks, None),
        ("same_effect_repeatable_twist", twisted, ranks, None),
        ("trace_preserving_rational_kraus_mixture", noisy, ranks, None),
        ("committed_pauli_y_lueders", y_ideal, (1, 1), y_basis),
        ("committed_pauli_y_swap_twist", y_twisted, (1, 1), y_basis),
        ("nonprojective_effect_pre_rotation_mixture", nonprojective, ranks, None),
    )


def run():
    rows = []
    for name, operators, ranks, basis in controls():
        bs, cs = probe_scores(operators, ranks, basis)
        ideal = tuple((p,) for p in projectors(ranks, basis))
        passes = all(x == 1 for x in bs+cs)
        equal = outcome_maps_equal(operators, ideal, ranks)
        effect_match = all(zero(sum((k.H*k for k in branch), s.zeros(sum(ranks)))-p[0])
                           for branch, p in zip(operators, ideal))
        if passes != equal:
            raise ArithmeticError("finite selection and independent map comparison disagree")
        errors = tuple(1-x for x in bs)
        cerrors = tuple(1-x for x in cs)
        rows.append({"case": name, "ranks": list(ranks),
                     "basis_joint_survival": list(map(str, bs)),
                     "coherent_joint_survival": list(map(str, cs)),
                     "passes_exact_protocol": passes,
                     "all_matrix_units_equal_lueders": equal,
                     "induced_effects_equal_target_projectors": effect_match,
                     "half_diamond_upper": str(calibration_bound(ranks, errors, cerrors))})
    delta = calibration_bound((2, 1), [F(1, 2**40)]*3, [F(1, 2**40)])
    return {"schema": "oph.finite_instrument_calibration.v1",
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in PARENTS},
            "controls": rows,
            "rational_error_example": {
                "ranks": [2, 1], "each_probability_error_upper": str(F(1, 2**40)),
                "half_diamond_upper": str(delta),
                "hundred_use_complete_history_upper": str(history_bound(delta, 100))},
            "theorem": "complete CPTP instrument selected by basis and block coherent survival",
            "trusted_preparation_and_readback_assumed": True,
            "quantum_probability_and_CPTP_semantics_assumed": True,
            "model_examples_are_not_observed_records": True,
            "lean_scope": "single full-space channel from basis and uniform return probabilities",
            "multi_label_selection_and_stability": "analytic proof and exact finite controls",
            "physical_readout_selected": False,
            "all_pass": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = run()
    path = HERE/"results.json"
    if args.check:
        if json.loads(path.read_text()) != data:
            raise SystemExit("finite instrument receipt drift")
        print("finite instrument calibration certificate OK")
    else:
        path.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
