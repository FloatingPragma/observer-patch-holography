"""Exact-factor source propagation of the calibrated free Weyl--Higgs action.

For one helicity pair H=e*(n_a+n_b)+m*(a†b†+ba), E=sqrt(m²+e²).
With chi=atan2(e,m) and C=exp(-i*chi*sigma_y/2) on the even sector,

    exp(-itH)=exp(-iet) C exp(-iEt*sigma_x) C†.

The two mixing rotations stay bounded while the central phase is reduced
with exact intervals. This is a source word, not a Trotter approximation or
an eigenvalue substituted for a detector record. Odd sectors and the whole
operator's scalar phase are retained. All calibrations remain supplied.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import calibrated_matter as matter, calibrated_source_frame as family_source
from . import calibrated_spectrum_observer as observer
from .chiral_candidate import gaussian_source as g, kinetic, protected_momentum as protected
from .chiral_candidate import dispersion_identifiability as bounds
from ..calibrated_measurement.bounds import root_interval


ROOT = Path(__file__).resolve().parents[2]
COMPILER = protected.compiler
ADMISSION = F(1, 10**10)
MOMENTUM_MAGNITUDES = (F(1), F(2))
PROBE_MAGNITUDES = (F(), *MOMENTUM_MAGNITUDES)
DEFAULT_LEVELS = 31
SEED = 20261027
QUADRATURES = ("X", "Y", "Z")


def _freeze(value):
    array = np.asarray(value, complex)
    return np.frombuffer(array.tobytes(), complex).reshape(array.shape)


def _exact(value, *, nonnegative=False):
    if isinstance(value, (bool, float, np.floating)):
        raise ValueError("an original exact finite real parameter required")
    result = s.sympify(value)
    if result.free_symbols or result.is_real is not True or result.is_finite is not True or result.atoms(s.Float):
        raise ValueError("an original exact finite real parameter required")
    if nonnegative and result.is_nonnegative is not True:
        raise ValueError("a nonnegative exact source parameter required")
    return result


def _time(value):
    if type(value) is not F or value < 0 or value > 2**48:
        raise ValueError("a nonnegative exact duration at most2**48 in inverse GeV required")
    return value


def energy(mass, kinetic_energy):
    m, e = _exact(mass, nonnegative=True), _exact(kinetic_energy)
    return s.sqrt(m*m+e*e)


@lru_cache(maxsize=4096)
def parameter_interval(value):
    """Exact algebraic enclosure with phase-safe root precision.

    The older rest receipt's1e-18 absolute square-root grid is deliberately
    unchanged. Multiplying it by a long programmed duration need not meet
    this source admission. Here roots use1e-40 and arithmetic rounds
    outward at1e-28 before phase multiplication.
    """
    value = s.sympify(value)
    if value.is_Rational:
        return (F(int(value.p), int(value.q)),)*2
    if value.is_Add:
        out = (F(), F())
        for term in value.args:
            out = bounds._add(out, parameter_interval(term))
        return bounds._compact(out)
    if value.is_Mul:
        out = (F(1), F(1))
        for term in value.args:
            out = bounds._mul(out, parameter_interval(term))
        return bounds._compact(out)
    if value.is_Pow and value.exp == s.Rational(1, 2):
        low, high = parameter_interval(value.base)
        if low < 0:
            raise ValueError("nonnegative exact algebraic radicand required")
        grid = 10**40
        def sqrt_floor(x):
            return F(math.isqrt(x.numerator*grid*grid//x.denominator), grid)
        lo, hi = sqrt_floor(low), sqrt_floor(high)
        if hi*hi < high:
            hi += F(1, grid)
        return lo, hi
    if value.is_Pow and value.exp.is_Integer:
        exponent = int(value.exp)
        box = parameter_interval(value.base)
        if exponent < 0:
            exponent, box = -exponent, bounds._inverse(box)
        out = (F(1), F(1))
        for _ in range(exponent):
            out = bounds._mul(out, box)
        return bounds._compact(out)
    raise ValueError("unsupported exact transport parameter: "+str(value))


def _complex_multiply(a, b):
    return bounds._cmul(a, b)


@lru_cache(maxsize=16384)
def _trig(box, cosine=False):
    # Pi reduction returns exact Machin-series rational endpoints. Compact
    # them outward before Taylor powers, preserving the enclosure while
    # preventing huge denominator powers at long clock windings.
    return bounds._trig(bounds._compact(box), cosine=cosine)


def _target_boxes(mass, kinetic_energy, duration):
    """Directed exact4-column propagator, including the odd-sector phase."""
    m, e, t = _exact(mass, nonnegative=True), _exact(kinetic_energy), _time(duration)
    E = energy(m, e)
    zero, one = (F(), F()), (F(1), F(1))
    result = [[(zero, zero) for _ in range(4)] for _ in range(4)]
    phase_box = tuple(x*t for x in parameter_interval(e))
    _, phase_box, _ = observer._reduced(phase_box)
    phase = _trig(phase_box, True), bounds._neg(_trig(phase_box))
    if E == 0:
        for i in range(4):
            result[i][i] = one, zero
        return result
    argument = tuple(x*t for x in parameter_interval(E))
    _, reduced, _ = observer._reduced(argument)
    co, si = _trig(reduced, True), _trig(reduced)
    inv = bounds._inverse(parameter_interval(E))
    me = bounds._mul(parameter_interval(m), inv)
    ee = bounds._mul(parameter_interval(e), inv)
    off = zero, bounds._neg(bounds._mul(me, si))
    result[0][0] = _complex_multiply(phase, (co, bounds._mul(ee, si)))
    result[3][3] = _complex_multiply(phase, (co, bounds._neg(bounds._mul(ee, si))))
    result[0][3] = result[3][0] = _complex_multiply(phase, off)
    result[1][1] = result[2][2] = phase
    return result


def moving_pair(mass, kinetic_energy, duration):
    """Execute all4 source columns of one exact helicity-pair evolution."""
    m, e, t = _exact(mass, nonnegative=True), _exact(kinetic_energy), _time(duration)
    E = energy(m, e)
    if max(abs(float(m)), abs(float(e))) > 256:
        raise ValueError("the fixed256GeV source window was exceeded")
    chi = math.atan2(float(e), float(m)) if E != 0 else 0.
    prepare = g.local_gate("pairing", 1j, chi/2)
    inverse = COMPILER.inverse(prepare.word)
    interval = tuple(x*t for x in parameter_interval(E))
    central, _, _, _ = observer._source_pair(interval)
    scalar_angle = tuple(-x*t for x in parameter_interval(e))
    reduced_phase, _, _ = observer._reduced(scalar_angle)
    scalar = [COMPILER.pulse("I", reduced_phase, [0])]
    word = list(inverse)+list(central.word)+list(prepare.word)+scalar
    unitary = _freeze(COMPILER.execute(word, 2))
    prep = _freeze(COMPILER.execute(prepare.word, 2))
    inv = _freeze(COMPILER.execute(inverse, 2))
    target = _target_boxes(m, e, t)
    error = observer._operator_error(unitary, target)
    if error > ADMISSION:
        raise ArithmeticError("exact moving-pair all-column source error exceeded admission")
    # Actual preparation/evolution/read-frame composition, not a target
    # evolution replaced after diagonalization. Global phase remains.
    effective = _freeze(inv@unitary@prep)
    expected_effective = [[entry for entry in row] for row in observer._pair_boxes(observer._reduced(interval)[1])]
    phase_interval = tuple(x*t for x in parameter_interval(e))
    _, phase_interval, _ = observer._reduced(phase_interval)
    phase = _trig(phase_interval, True), bounds._neg(_trig(phase_interval))
    expected_effective = [[_complex_multiply(phase, entry) for entry in row] for row in expected_effective]
    effective_error = observer._operator_error(effective, expected_effective)
    if effective_error > ADMISSION:
        raise ArithmeticError("source energy-frame preparation/readout composition failed")
    return {"unitary": unitary, "preparation": prep, "inverse_preparation": inv,
            "effective_read_unitary": effective, "energy_GeV": E,
            "word": tuple(word), "prepare_word": tuple(prepare.word), "inverse_prepare_word": tuple(inverse),
            "resources": COMPILER.resources(word),
            "preparation_resources": COMPILER.resources(prepare.word),
            "inverse_preparation_resources": COMPILER.resources(inverse),
            "exact_binary_operator_error_upper": error,
            "exact_binary_effective_read_operator_error_upper": effective_error,
            "bounded_mixing_angle": chi, "global_phase_removed_or_fitted": False,
            "time_parameter": "programmed inverse-GeV action parameter; same global clock dictionary as calibrated rest observer"}


def helicity_frame(axis):
    """Actual source SU2 spin frame S, S†sigma_axis S=sigma_z."""
    if type(axis) is not str or axis not in ("x", "y", "z"):
        raise ValueError("one of the three fixed spatial axes required")
    if axis == "x":
        gate = g.local_gate("hopping", -1j, math.pi/4)
        target = np.array([[1, -1], [1, 1]], complex)/math.sqrt(2)
    elif axis == "y":
        gate = g.local_gate("hopping", -1., math.pi/4)
        target = np.array([[1, 1j], [1j, 1]], complex)/math.sqrt(2)
    else:
        gate = g.local_gate("hopping", 1j, 0.)
        target = np.eye(2, dtype=complex)
    actual = _freeze(COMPILER.execute(gate.word, 2))
    spin = actual[np.ix_((1, 2), (1, 2))]
    expected = np.zeros((4, 4), complex)
    expected[0, 0], expected[3, 3] = 1, np.linalg.det(target)
    expected[np.ix_((1, 2), (1, 2))] = target
    residual = float(np.linalg.norm(actual-expected, 2))
    sigma = kinetic.pauli_matrices()[("x", "y", "z").index(axis)]
    epsilon = kinetic.spin_epsilon()
    if (residual > float(ADMISSION) or np.linalg.norm(spin.conj().T@sigma@spin-np.diag([1., -1.]), 2) > float(ADMISSION)
            or np.linalg.norm(spin.T@epsilon@spin-epsilon, 2) > float(ADMISSION)):
        raise ArithmeticError("source helicity frame violates spin, chirality or epsilon convention")
    return {"spinor_unitary": _freeze(spin), "Fock_unitary": actual,
            "word": tuple(gate.word), "resources": dict(gate.resources),
            "source_Fock_operator_difference_diagnostic": residual}


def momentum_menu():
    """One finite spectral family in a supplied2pi coordinate cube, in GeV."""
    values = [(F(), F(), F())]
    for magnitude in MOMENTUM_MAGNITUDES:
        for axis in range(3):
            for sign in (1, -1):
                p = [F(), F(), F()]
                p[axis] = sign*magnitude
                values.append(tuple(p))
    return tuple(values)


def source_layout():
    """All1170 original CAR slots and their disjoint canonical Fock factors."""
    momenta = momentum_menu()
    indices = {p: 90*i for i, p in enumerate(momenta)}
    factors, order = [], []
    labels = observer.catalogue()
    def add(original, row, e, sign, block):
        start = len(order)
        order.extend(original)
        factors.append({"label": row["label"], "channel": row["channel"], "family": row["family"],
                        "mass_GeV": row["mass_GeV"], "e_GeV": e, "pair_sign": sign,
                        "momentum_block": block, "original_CAR_indices": list(original),
                        "canonical_CAR_indices": [start, start+1]})
    for row in labels:
        for (left, right) in matter._oriented_components(row["channel"])[1].todok():
            i, j = 3*left+row["family"], 3*right+row["family"]
            add((i, 45+j), row, F(), 1, "zero")
            add((j, 45+i), row, F(), 1, "zero")
    neutral = []
    used = set(order)
    for i in range(90):
        if i not in used:
            neutral.append({"original_CAR_index": i, "e_GeV": F(), "canonical_CAR_index": len(order)})
            order.append(i)
    for magnitude in MOMENTUM_MAGNITUDES:
        for axis, axis_name in enumerate(("x", "y", "z")):
            p = [F(), F(), F()]; p[axis] = magnitude
            plus, minus = indices[tuple(p)], indices[tuple(-x for x in p)]
            before = set(order)
            for row in labels:
                for left, right in matter._oriented_components(row["channel"])[1].todok():
                    i, j = 3*left+row["family"], 3*right+row["family"]
                    for a, b, e, sign in ((i, 45+j, -magnitude, 1), (j, 45+i, -magnitude, 1),
                                           (45+i, j, magnitude, -1), (45+j, i, magnitude, -1)):
                        add((plus+a, minus+b), row, e, sign, f"{axis_name}:{magnitude}")
            used = set(order)-before
            for offset, orientation in ((plus, 1), (minus, -1)):
                for i in range(90):
                    original = offset+i
                    if original not in used:
                        e = orientation*magnitude*(-1 if i < 45 else 1)
                        neutral.append({"original_CAR_index": original, "e_GeV": e,
                                        "canonical_CAR_index": len(order)})
                        order.append(original)
    if sorted(order) != list(range(1170)) or len(factors) != 546 or len(neutral) != 78:
        raise ArithmeticError("finite momentum family discarded or duplicated a physical CAR mode")
    inverse = {original: canonical for canonical, original in enumerate(order)}
    return {"momenta": momenta, "pairs": tuple(factors), "neutral": tuple(neutral),
            "canonical_to_original": tuple(order), "original_to_canonical": inverse}


def _route_plan(sites):
    """Compact reversible adjacent-FSWAP route in declared factor order."""
    start = min(sites)
    ordering = list(range(start, max(sites)+1))
    legs = []
    for i, site in enumerate(sites):
        position = start+ordering.index(site)
        target = start+i
        if position < target:
            raise ArithmeticError("fermion-factor extraction changed an earlier selected mode")
        legs.append([target, position])
        value = ordering.pop(position-start)
        ordering.insert(target-start, value)
    if ordering[:len(sites)] != list(sites):
        raise ArithmeticError("source route failed its declared mode permutation")
    return {"contiguous_start": start, "left_bubble_ranges": legs,
            "forward_and_inverse_FSWAP_calls": 2*sum(b-a for a, b in legs)}


def frame_routing(layout=None):
    """One compact source-word grammar for all weak/helicity field frames.

    Every local Fock factor is actually contracted. Reversible FSWAP routes
    place that ordered mode subset contiguously, apply the local word and
    reverse the route. The mode map and all routing multiplicities are
    retained; no qubit embedding silently drops fermionic spectator signs.
    """
    layout = source_layout() if layout is None else layout
    inverse = layout["original_to_canonical"]
    family = family_source.source_frame()
    local_admission = F(1, 5*10**11)
    inverse_family = COMPILER.execute(COMPILER.inverse(family["sourceword"]), 3)
    if (family["diagnostics"]["CKM_source_operator_difference"] > float(local_admission)
            or family["diagnostics"]["singlet_phase_source_difference"] > float(local_admission)
            or np.linalg.norm(inverse_family-family["CKM_fock_target"].conj().T, 2) > float(local_admission)):
        raise ArithmeticError("local calibrated family frame exceeds its tighter field admission")
    swap = g.local_gate("fswap")
    # Exact FSWAP entries in the same little-endian four-state basis.
    z, p, n = (F(), F()), (F(1), F(1)), (F(-1), F(-1))
    swap_boxes = [[(z, z) for _ in range(4)] for _ in range(4)]
    for i, j, value in ((0, 0, p), (1, 2, p), (2, 1, p), (3, 3, n)):
        swap_boxes[i][j] = value, z
    swap_error = observer._operator_error(swap.fock_unitary, swap_boxes)
    swap_admission = F(1, 10**13)
    if swap_error > swap_admission:
        raise ArithmeticError("binary FSWAP exceeds its fixed operator admission")
    factors, swap_calls, pulses = [], 0, 0
    local_calls = 0
    for block, momentum in enumerate(layout["momenta"]):
        offset = 90*block
        axis = next((name for i, name in enumerate(("x", "y", "z")) if momentum[i]), "z")
        spin = helicity_frame(axis)
        inverse_spin = COMPILER.execute(COMPILER.inverse(spin["word"]), 2)
        if (spin["source_Fock_operator_difference_diagnostic"]
                +np.linalg.norm(inverse_spin-spin["Fock_unitary"].conj().T, 2) > float(local_admission)):
            raise ArithmeticError("local spin frame exceeds its tighter field admission")
        for internal in range(45):
            sites = tuple(inverse[offset+i] for i in (internal, 45+internal))
            route = _route_plan(sites)
            factors.append({"kind": "spin_"+axis, "sites": list(sites), "route": route})
            calls = route["forward_and_inverse_FSWAP_calls"]
            swap_calls += calls
            pulses += spin["resources"]["native_pulses"]+calls*swap.resources["native_pulses"]
            local_calls += 1
        for row in family["mapping"]["CKM_blocks"]:
            sites = tuple(inverse[offset+i] for i in row["CAR_indices"])
            route = _route_plan(sites)
            factors.append({"kind": "CKM", "sites": list(sites), "route": route})
            calls = route["forward_and_inverse_FSWAP_calls"]
            swap_calls += calls
            pulses += family["CKM_block_resources"]["native_pulses"]+calls*swap.resources["native_pulses"]
            local_calls += 1
        for row in family["mapping"]["negative_singlet_phases"]:
            site = inverse[offset+row["CAR_index"]]
            factors.append({"kind": "singlet_sign", "sites": [site]})
            pulses += family["singlet_phase_resources"]["native_pulses"]
            local_calls += 1
    summed = local_calls*F(1, 5*10**11)+swap_calls*swap_admission
    if summed >= F(1, 2):
        raise ArithmeticError("routed source frame is outside its finite product-error domain")
    # prod(1+epsilon_j)-1 <= S/(1-S), S=sum(epsilon_j)<1;
    # this rational bound avoids constructing giant intermediate powers.
    error = summed/(1-summed)
    # Four frame operations account for preparation, weak evolution and
    # returning to the read coordinates. Local admissions here are numerical
    # all-input conditions; the routed FSWAP binary-array bound is exact.
    four = (1+error)**4-1
    score = 2*four+four*four
    if score > observer.WEAK_FRAME_SCORE_ALLOWANCE:
        raise ArithmeticError("complete routed weak/spin-frame score budget exceeded")
    return {"factor_plan": factors, "factor_count": len(factors), "FSWAP_calls_per_frame": swap_calls,
            "native_pulses_per_frame": pulses, "source_native_occupation_wire_rule": "canonicalCARj -> decodedwire9*j",
            "local_spin_CKM_sign_operator_admission_assumed": "1/500000000000",
            "exact_binary_FSWAP_operator_error_upper": str(swap_admission),
            "backend_FSWAP_error_diagnostic": str(swap_error),
            "four_frame_score_error_upper_conditional": str(observer._outward((F(), score))[1]),
            "all_field_frame_arrays_materialized": False,
            "local_admissions_are_physical_or_outward_certificates": False}


def bind_source_word(word, sites):
    """Validate literal occupation-wire relocation before factor composition."""
    literal = protected.decoded_source_word(word, sites)
    lowered = protected._lower_decoded_word(literal, sites)
    if lowered != list(word):
        raise ArithmeticError("decoded occupation wiring changed the ordered CAR source word")
    return literal


def _signed_pair(mass, e, sign, duration):
    base = moving_pair(mass, e, duration)
    if sign == 1:
        return base
    phase_word = g.local_gate("onsite", 1., math.pi).word
    phase = COMPILER.execute(phase_word, 2)
    inverse_phase = COMPILER.execute(COMPILER.inverse(phase_word), 2)
    unitary = _freeze(phase@base["unitary"]@inverse_phase)
    prep = _freeze(phase@base["preparation"])
    inverse = _freeze(base["inverse_preparation"]@inverse_phase)
    effective = _freeze(inverse@unitary@prep)
    result = dict(base)
    result.update(unitary=unitary, preparation=prep, inverse_preparation=inverse,
                  effective_read_unitary=effective,
                  word=tuple(COMPILER.inverse(phase_word))+base["word"]+tuple(phase_word),
                  prepare_word=base["prepare_word"]+tuple(phase_word),
                  inverse_prepare_word=tuple(COMPILER.inverse(phase_word))+base["inverse_prepare_word"])
    result["resources"] = COMPILER.resources(result["word"])
    result["preparation_resources"] = COMPILER.resources(result["prepare_word"])
    result["inverse_preparation_resources"] = COMPILER.resources(result["inverse_prepare_word"])
    target = _target_boxes(mass, e, duration)
    for i, j in ((0, 3), (3, 0)):
        real, imaginary = target[i][j]
        target[i][j] = bounds._neg(real), bounds._neg(imaginary)
    result["exact_binary_operator_error_upper"] = observer._operator_error(unitary, target)
    E = energy(mass, e)
    box = tuple(x*duration for x in parameter_interval(E))
    even = observer._pair_boxes(observer._reduced(box)[1])
    _, scalar, _ = observer._reduced(tuple(x*duration for x in parameter_interval(e)))
    phase_box = _trig(scalar, True), bounds._neg(_trig(scalar))
    expected = [[_complex_multiply(phase_box, entry) for entry in row] for row in even]
    result["exact_binary_effective_read_operator_error_upper"] = observer._operator_error(effective, expected)
    if max(result["exact_binary_operator_error_upper"], result["exact_binary_effective_read_operator_error_upper"]) > ADMISSION:
        raise ArithmeticError("signed source factor exceeded its exact binary admission")
    return result


def _groups(layout):
    groups, lookup = [], {}
    for row in layout["pairs"]:
        key = (row["label"], row["e_GeV"], row["pair_sign"])
        if key not in lookup:
            lookup[key] = len(groups)
            groups.append({key: row[key] for key in ("label", "channel", "family", "mass_GeV", "e_GeV", "pair_sign")}
                          | {"multiplicity": 0})
        groups[lookup[key]]["multiplicity"] += 1
    return tuple(groups), lookup


def _field_at(duration, groups, neutral_groups):
    pairs, outputs, inputs, effective = [], [], [], []
    norms, initial_norms = [], []
    resources = {"pair_evolution": 0, "pair_preparation": 0, "neutral_evolution": 0}
    max_error = F()
    for row in groups:
        pair = _signed_pair(row["mass_GeV"], row["e_GeV"], row["pair_sign"], duration)
        pairs.append(pair)
        fixed = fixed_probe()
        initial = _freeze(fixed["unitary"][:, 0])
        final = _freeze(pair["unitary"]@initial)
        inputs.append(initial); outputs.append(final)
        norms.append(observer._norm(final)); initial_norms.append(observer._norm(initial))
        effective.append(_freeze(pair["unitary"]@fixed["unitary"]))
        resources["pair_evolution"] += row["multiplicity"]*pair["resources"]["native_pulses"]
        resources["pair_preparation"] += row["multiplicity"]*fixed["resources"]["native_pulses"]
        max_error = max(max_error, pair["exact_binary_operator_error_upper"], pair["exact_binary_effective_read_operator_error_upper"])
    neutrals, neutral_norm = [], F(1)
    for e, multiplicity in neutral_groups:
        angle = e*duration
        reduced, _, _ = observer._reduced((angle, angle))
        gate = g.local_gate("onsite", 1., reduced)
        matrix = _freeze(COMPILER.execute(gate.word, 1))
        neutrals.append(matrix)
        neutral_norm *= observer._norm(matrix[:, 0])**multiplicity
        resources["neutral_evolution"] += multiplicity*gate.resources["native_pulses"]
    return {"pairs": pairs, "effective": effective, "inputs": inputs, "outputs": outputs,
            "norms": norms, "initial_norms": initial_norms, "neutrals": neutrals,
            "neutral_norm": neutral_norm, "resources": resources, "max_pair_error": max_error}


def fixed_probe():
    """One common+Y even-pair state, independent of every action parameter."""
    gate = g.local_gate("pairing", 1., -math.pi/4)
    actual = _freeze(COMPILER.execute(gate.word, 2))
    return {"unitary": actual, "word": tuple(gate.word), "resources": dict(gate.resources)}


def fixed_read(field, quadrature):
    """Three fixed even-pair Pauli analyzers with both number-pointer outputs."""
    if quadrature not in QUADRATURES or type(quadrature) is not str:
        raise ValueError("one of the fixedX/Y/Z pair quadratures required")
    coefficient, angle = {"X": (1j, -math.pi/4), "Y": (1., math.pi/4), "Z": (1., 0.)}[quadrature]
    gate = g.local_gate("pairing", coefficient, angle)
    changed = gate.fock_unitary@field
    word, kraus, _ = g._occupation_pointer()
    outputs = tuple(protected.model.apply_subsystem(changed, k, [0], 2) for k in kraus)
    return outputs, gate.resources["native_pulses"]+COMPILER.resources(word)["native_pulses"]


def signed_sine_interval(minus_x, z):
    """From(−X,Z)=A*(p/E,m/E)*sin(theta), reconstruct signedA*sin.

    The sign is certified by either component. Conflicting certified signs
    reject the declared model. If neither sign resolves, the symmetric
    interval remains valid through the zero crossing; no division by a
    small mass, visibility or sine is performed.
    """
    def valid(box):
        if type(box) is not tuple or len(box) != 2 or any(type(x) is not F for x in box) or box[0] > box[1]:
            raise ValueError("ordered original Fraction quadrature boxes required")
        return box
    minus_x, z = valid(minus_x), valid(z)
    positive = minus_x[0] > 0 or z[0] > 0
    negative = minus_x[1] < 0 or z[1] < 0
    if positive and negative:
        raise ValueError("opposite certified signs contradict the oriented source model")
    def square(box):
        low = F() if box[0] <= 0 <= box[1] else min(x*x for x in box)
        return low, max(x*x for x in box)
    a, b = square(minus_x), square(z)
    lo = root_interval(a[0]+b[0])[0]
    hi = root_interval(a[1]+b[1])[1]
    return (lo, hi) if positive else (-hi, -lo) if negative else (-hi, hi)


def rotate_intervals(cosine, sine, offset):
    """Classical interval rotation; the quantum probe/analyzers stay fixed."""
    _, reduced, _ = observer._reduced((offset, offset))
    c, q = _trig(reduced, True), _trig(reduced)
    return (bounds._add(bounds._mul(c, cosine), bounds._mul(q, sine)),
            bounds._add(bounds._neg(bounds._mul(q, cosine)), bounds._mul(c, sine)))


def _record_means(counts):
    patterns = np.arange(512)
    means = {label: {} for label in observer.RECORD_LABELS}
    for setting, axis in enumerate(QUADRATURES):
        accepted = counts[setting, :30].sum(axis=0)
        for index, label in enumerate(observer.RECORD_LABELS[:-1]):
            signs = 1-2*((patterns >> index) & 1)
            means[label][axis] = F(int(np.sum(accepted*signs[:, None])), observer.TRIALS_PER_SETTING)
        if axis == "Y":
            for bit, key in ((1, "C"), (0, "S")):
                signs = np.asarray([-1 if ((record >> bit) & 1) == 0 else 1 for record in range(4)])
                means["clock"][key] = F(int(np.sum(accepted*signs)), observer.TRIALS_PER_SETTING)
    return means


def _update(priors, means, level):
    from .calibrated_phase_inference import phase_update
    result = {}
    radius = observer._score_radius()
    for label in observer.RECORD_LABELS:
        values = means[label]
        if label == "clock":
            cosine = values["C"]-radius, values["C"]+radius
            sine = values["S"]-radius, values["S"]+radius
        else:
            cosine = values["Y"]-radius, values["Y"]+radius
            minus_x = -values["X"]-radius, -values["X"]+radius
            z = values["Z"]-radius, values["Z"]+radius
            sine = signed_sine_interval(minus_x, z)
        offset = sum(priors[label])*2**level/2
        residual = tuple(observer._outward(box) for box in rotate_intervals(cosine, sine, offset))
        update = phase_update(priors[label], level, *residual)
        if not update["certified"]:
            raise ValueError("fixed-probe phase is unresolved: "+update["reason"])
        result[label] = update["phase_base_interval"]
    return result


def _counts(value):
    if (type(value) is not np.ndarray or value.dtype != np.dtype("int64") or value.ndim != 6
            or value.shape[0] != 3 or value.shape[2:] != (3, 32, 512, 4)
            or not 1 <= value.shape[1] <= 40 or np.any(value < 0)
            or np.any(value > observer.TRIALS_PER_SETTING)
            or np.any(value.sum(axis=(3, 4, 5)) != observer.TRIALS_PER_SETTING)):
        raise ValueError("complete original3protocol/level/3setting/address/field/clock counts required")
    return value


def decode_records(counts):
    """Packet-free fixed-probe inverse: counts, global units and momenta only."""
    counts = _counts(counts)
    recovered = []
    for data in counts:
        priors = {label: (F(), F(1)) for label in observer.RECORD_LABELS}
        for level, settings in enumerate(data):
            priors = _update(priors, _record_means(settings), level)
        clock = priors["clock"]
        if clock[0] <= 0:
            raise ValueError("the common source clock remains unresolved")
        recovered.append({"phase_intervals": priors,
            "energy_GeV_intervals": {label: (priors[label][0]/(2*clock[1]), priors[label][1]/(2*clock[0]))
                                     for label in observer.RECORD_LABELS[:-1]}})
    return tuple(recovered)


def execute(*, levels=DEFAULT_LEVELS, seed=SEED, retained_counts=None):
    """One fixed-probe p=0,1,2 record experiment on the common finite field."""
    if type(levels) is not int or not 1 <= levels <= 40 or type(seed) is not int or seed < 0:
        raise ValueError("one through40 dyadic levels and a nonnegative integer seed required")
    if retained_counts is not None:
        retained_counts = _counts(retained_counts)
        if retained_counts.shape[1] != levels:
            raise ValueError("retained count horizon changed")
    layout = source_layout()
    frames = frame_routing(layout)
    groups, lookup = _groups(layout)
    neutral_table = {}
    for row in layout["neutral"]:
        neutral_table[row["e_GeV"]] = neutral_table.get(row["e_GeV"], 0)+1
    neutral_groups = tuple(sorted(neutral_table.items()))
    rows = observer.catalogue()
    energies = [[parameter_interval(energy(row["mass_GeV"], p)) for row in rows] for p in PROBE_MAGNITUDES]
    if any(high >= 256 for protocol in energies for _, high in protocol):
        raise ValueError("the common256GeV phase window was exceeded")
    selected = [[lookup[(row["label"], -p, 1)] for row in rows] for p in PROBE_MAGNITUDES]
    clocks_initial, clock_calibration, clock_cost = observer.history._clock_preparation()
    contrast = -(F(clock_calibration["populations"][0])-F(clock_calibration["populations"][1]))
    coin = COMPILER.execute([COMPILER.pulse("H", math.pi/2, [0])], 1)[:, 0]
    selector_read = observer.instruments.pauli_read(coin, 0, axis="Z")
    coin_weights = [observer._norm(value) for value in selector_read.outputs]
    selector = tuple(math.prod(coin_weights[(a >> bit) & 1] for bit in range(5)) for a in range(32))
    accepted, failed = sum(selector[:30], F()), sum(selector[30:], F())
    priors = [{label: (F(), F(1)) for label in observer.RECORD_LABELS} for _ in PROBE_MAGNITUDES]
    rngs = [np.random.Generator(np.random.PCG64(seed+i)) for i in range(3)]
    arrays, pair_unitaries, pair_outputs, neutral_unitaries = {}, [], [], []
    success_arrays, failure_arrays, clock_arrays, clock_failed_arrays = [], [], [], []
    count_rows = [[] for _ in PROBE_MAGNITUDES]
    max_bias, max_operator = F(), F()
    resources = {"charged_evolution_native_pulses": 0, "neutral_evolution_native_pulses": 0,
                 "fixed_pair_preparation_native_pulses": 0, "fixed_pair_read_native_pulses": 0,
                 "full_routed_spin_family_frame_native_pulses": 4*3*3*levels*frames["native_pulses_per_frame"],
                 "clock_preparation_native_pulses": 9*levels*clock_cost["native_pulses"],
                 "clock_advance_native_pulses": 36*levels,
                 "clock_pointer_native_pulses": 0,
                 "fair_selector_native_pulses": 45*levels*(1+selector_read.resources["native_pulses"])}
    patterns = np.arange(512)
    for level in range(levels):
        time = F(2**level, observer.BASE_DENOMINATOR)
        field = _field_at(time, groups, neutral_groups)
        pair_unitaries.append(np.stack([value["unitary"] for value in field["pairs"]]))
        pair_outputs.append(np.stack(field["outputs"]))
        neutral_unitaries.append(np.stack(field["neutrals"]))
        if level == 0:
            arrays["pair_initial_preparations"] = np.stack([fixed_probe()["unitary"] for _ in groups])
            for factor in layout["pairs"]:
                index = lookup[(factor["label"], factor["e_GeV"], factor["pair_sign"])]
                bind_source_word(field["pairs"][index]["word"], tuple(factor["canonical_CAR_indices"]))
        max_operator = max(max_operator, field["max_pair_error"])
        for key, raw in (("charged_evolution_native_pulses", "pair_evolution"),
                         ("neutral_evolution_native_pulses", "neutral_evolution"),
                         ("fixed_pair_preparation_native_pulses", "pair_preparation")):
            resources[key] += 9*field["resources"][raw]
        clocks, cp, clock_means, clock_pointer = observer._clock(clocks_initial, time, F())
        clocks_failed, cpf, _, _ = observer._clock(clocks_initial, F(), F())
        clock_arrays.append(np.stack(clocks)); clock_failed_arrays.append(np.stack(clocks_failed))
        resources["clock_pointer_native_pulses"] += 9*sum(row["native_pulses"] for row in clock_pointer[:2])
        clock_mass, failed_clock_mass = sum(cp.values(), F()), sum(cpf.values(), F())
        level_success, level_failed = [], []
        for protocol, magnitude in enumerate(PROBE_MAGNITUDES):
            chosen = selected[protocol]
            spectators = field["neutral_norm"]
            initial_spectators = F(1)
            for index, row in enumerate(groups):
                power = row["multiplicity"]-int(index in chosen)
                spectators *= field["norms"][index]**power
                initial_spectators *= field["initial_norms"][index]**power
            setting_counts, success, unsuccessful, exact_means = [], [], [], {label: {} for label in observer.RECORD_LABELS}
            for setting, axis in enumerate(QUADRATURES):
                outputs, aborts, probabilities, failed_probabilities = [], [], [], []
                for index in chosen:
                    measured, cost = fixed_read(field["effective"][index], axis)
                    failure, _ = fixed_read(fixed_probe()["unitary"], axis)
                    outputs.append(np.stack(measured)); aborts.append(np.stack(failure))
                    probabilities.append(tuple(observer._norm(v[:, 0]) for v in measured))
                    failed_probabilities.append(tuple(observer._norm(v[:, 0]) for v in failure))
                    resources["fixed_pair_read_native_pulses"] += cost
                totals = [sum(v, F()) for v in probabilities]
                failed_totals = [sum(v, F()) for v in failed_probabilities]
                field_mass = spectators*math.prod(totals)
                failed_mass = initial_spectators*math.prod(failed_totals)
                if abs(accepted*field_mass*clock_mass+failed*failed_mass*failed_clock_mass-1) > F(1, 10**7):
                    raise ArithmeticError("complete moving field/clock/selector instrument lost raw mass")
                for lane, row in enumerate(rows):
                    mean = accepted*spectators*clock_mass*(probabilities[lane][0]-probabilities[lane][1])
                    mean *= math.prod(totals[j] for j in range(9) if j != lane)
                    exact_means[row["label"]][axis] = mean
                    phase = tuple(2*x*time for x in energies[protocol][lane])
                    _, phase, _ = observer._reduced(phase)
                    if axis == "Y":
                        ideal = _trig(phase, True)
                    else:
                        ratio = (-magnitude if axis == "X" else row["mass_GeV"])/energy(row["mass_GeV"], magnitude)
                        ideal = bounds._mul(parameter_interval(ratio), _trig(phase))
                    ideal = tuple(accepted_ideal*x for accepted_ideal in (observer.IDEAL_ACCEPTED,) for x in ideal)
                    bias = observer._distance_to_interval(mean, ideal)
                    max_bias = max(max_bias, bias)
                    if bias > observer.NUMERIC_SCORE_ADMISSION:
                        raise ArithmeticError("fixed-probe source Bloch-score error exceeds its numerical allowance")
                field_law, abort_law = np.ones(512), np.ones(512)
                for lane in range(9):
                    bits = (patterns >> lane) & 1
                    field_law *= np.asarray([float(x/totals[lane]) for x in probabilities[lane]])[bits]
                    abort_law *= np.asarray([float(x/failed_totals[lane]) for x in failed_probabilities[lane]])[bits]
                clock_law = np.asarray([float(x/clock_mass) for x in cp.values()])
                failed_clock_law = np.asarray([float(x/failed_clock_mass) for x in cpf.values()])
                total_selector = accepted+failed
                law = np.stack([float(weight/total_selector)*(field_law[:, None]*clock_law if address < 30
                              else abort_law[:, None]*failed_clock_law) for address, weight in enumerate(selector)])
                law /= law.sum()
                if retained_counts is None:
                    counts = rngs[protocol].multinomial(observer.TRIALS_PER_SETTING, law.ravel()).reshape(32, 512, 4)
                else:
                    counts = retained_counts[protocol, level, setting]
                setting_counts.append(counts)
                for lane, row in enumerate(rows):
                    normalized = accepted/total_selector*(probabilities[lane][0]-probabilities[lane][1])/totals[lane]
                    if abs(normalized-exact_means[row["label"]][axis]) > F(1, 10**9):
                        raise ArithmeticError("moving sampler normalization exceeds the separate numerical budget")
                if axis == "Y":
                    for bit, key in enumerate(("C", "S")):
                        mean = accepted*field_mass*clock_means[bit]
                        exact_means["clock"][key] = mean
                        ideal = tuple(observer.IDEAL_ACCEPTED*contrast*x for x in _trig((time, time), key == "C")) if time <= 4 else None
                        if ideal is None:
                            _, phase, _ = observer._reduced((time, time))
                            ideal = tuple(observer.IDEAL_ACCEPTED*contrast*x for x in _trig(phase, key == "C"))
                        bias = observer._distance_to_interval(mean, ideal)
                        max_bias = max(max_bias, bias)
                        if bias > observer.NUMERIC_SCORE_ADMISSION:
                            raise ArithmeticError("common fixed clock source error exceeded numerical allowance")
                        normalized = accepted/total_selector*clock_means[bit]/clock_mass
                        if abs(normalized-mean) > F(1, 10**9):
                            raise ArithmeticError("clock sampler normalization exceeds its separate numerical budget")
                success.append(np.stack(outputs)); unsuccessful.append(np.stack(aborts))
            count_rows[protocol].append(np.stack(setting_counts))
            priors[protocol] = _update(priors[protocol], _record_means(np.stack(setting_counts)), level)
            level_success.append(np.stack(success)); level_failed.append(np.stack(unsuccessful))
        success_arrays.append(np.stack(level_success)); failure_arrays.append(np.stack(level_failed))
    counts = np.stack([np.stack(values) for values in count_rows])
    decoded = decode_records(counts)
    if any(row["phase_intervals"] != prior for row, prior in zip(decoded, priors, strict=True)):
        raise ArithmeticError("independent fixed-probe count-only inverse disagrees")
    recovered = []
    for protocol, magnitude in enumerate(PROBE_MAGNITUDES):
        for lane, row in enumerate(rows):
            box = decoded[protocol]["energy_GeV_intervals"][row["label"]]
            target = energies[protocol][lane]
            if not box[0] <= target[0] <= target[1] <= box[1]:
                raise ArithmeticError("fixed-probe record energy interval missed the supplied action")
            square = max(F(), box[0]**2-magnitude*magnitude), box[1]**2-magnitude*magnitude
            if square[1] < 0:
                raise ArithmeticError("record dispersion contradicts nonnegative invariant mass squared")
            mass_box = root_interval(square[0])[0], root_interval(square[1])[1]
            recovered.append({"label": row["label"], "momentum_GeV": str(magnitude),
                              "energy_interval_GeV": list(map(str, box)),
                              "energy_interval_numeric_GeV": list(map(float, box)),
                              "invariant_mass_interval_GeV": list(map(str, mass_box)),
                              "invariant_mass_interval_numeric_GeV": list(map(float, mass_box)),
                              "positive_invariant_mass_resolved": mass_box[0] > 0})
    for lane, row in enumerate(rows):
        rest = decoded[0]["energy_GeV_intervals"][row["label"]]
        for protocol in (1, 2):
            found = recovered[9*protocol+lane]["invariant_mass_interval_GeV"]
            if max(rest[0], F(found[0])) > min(rest[1], F(found[1])):
                raise ArithmeticError("moving record-derived invariant mass is inconsistent with rest records")
    sampling = observer.concentration_certificate(3*29*levels, observer.TRIALS_PER_SETTING,
                                                 math.nextafter(float(observer.RECORD_RADIUS), 0.))
    if F(sampling["familywise_error_upper"]) >= F(1, 100):
        raise ArithmeticError("simultaneous fixed-probe sampling certificate failed")
    arrays.update(pair_unitaries=np.stack(pair_unitaries), pair_outputs=np.stack(pair_outputs),
                  neutral_unitaries=np.stack(neutral_unitaries), complete_record_counts=counts,
                  success_read_matrices=np.stack(success_arrays), failed_read_matrices=np.stack(failure_arrays),
                  clock_reference_outputs=np.stack(clock_arrays), failed_clock_reference_outputs=np.stack(clock_failed_arrays),
                  selector_coin=coin, selector_read_outputs=np.stack(selector_read.outputs))
    arrays["helicity_Fock_frames"] = np.stack([helicity_frame(axis)["Fock_unitary"] for axis in ("x", "y", "z")])
    arrays["weak_CKM_Fock_frame"] = family_source.source_frame()["CKM_fock_unitary"]
    arrays["opposite_momentum_source_controls"] = np.stack([
        charged_momentum_pair(F(7, 10), F(1), axis, F(1, 8))["unitary"] for axis in ("x", "y", "z")])
    cost = sum(resources.values())
    resources.update(native_pulses_per_one_trial_setting_schedule=cost,
                     charged_load_export_requirements=2*cost,
                     hypothetical_full_campaign_native_pulses=observer.TRIALS_PER_SETTING*cost,
                     hypothetical_full_campaign_service_events=2*observer.TRIALS_PER_SETTING*cost)
    report = {"schema": "oph.calibrated_transport.v1", "dyadic_levels": levels, "sampler_seed": seed,
        "field": {"coordinate_box_side_GeV_inverse": "2*pi", "momentum_menu_GeV": [[str(x) for x in p] for p in layout["momenta"]],
                  "CAR_modes": 1170, "charged_pair_factors": 546, "neutral_modes": 78,
                  "zero_momentum_duplicated": False,
                  "pair_groups": [{key: str(value) if isinstance(value, (s.Basic, F)) else value for key, value in row.items()} for row in groups],
                  "canonical_to_original_CAR_indices": list(layout["canonical_to_original"]),
                  "neutral_groups": [[str(e), count] for e, count in neutral_groups]},
        "complete_symbol_and_current_controls": symbol_intertwining(),
        "source_wiring": {"decoded_occupation_wire_rule": "canonicalCARj -> wire9*j; syndrome wires excluded",
                          "disjoint_pair_factor_bindings_checked": len(layout["pairs"]),
                          "native_code_load_recovery_export_executed": False},
        "fixed_probe": {"preparation": "+Y even-pair state, same for every mass and momentum",
                        "analyzers": "fixed even-pair PauliX/Y/Z; no parameter-dependent preparation/read rotation",
                        "Hamiltonian_compilation": "bounded mass/momentum-dependent C rotations occur only inside the actual propagator",
                        "classical_phase_compensation": "rotate interval quadratures using prior-count midpoint; no adaptive quantum analyzer",
                        "signed_sine_identity": "for selected e=-p,(-X,Z)=A*(p/E,m/E)*sin(2Et); oriented Euclidean norm gives A*sin",
                        "inverse_receives_action_parameters": False},
        "common_clock": {"native_gap": "log(5)", "global_energy_dictionary": "1GeV/log(5)",
                         "source_modular_parameter": "2**level/(512*log(5))", "per_particle_fit": False},
        "records": {"kind": "complete fixed-seed simulated counts, not laboratory data",
                    "history_rule": "separate full-field preparation/evolution/read histories per setting and trial; one common action, momentum family and clock",
                    "trials_per_setting": observer.TRIALS_PER_SETTING, "settings_per_level": 9,
                    "all_attempt_trials": 9*levels*observer.TRIALS_PER_SETTING,
                    "outcomes_per_setting": 65536, "selector_timeouts_retained": True,
                    "complete_unnormalized_quantum_factors_retained": True,
                    "quantum_factor_rule": "charged+Y factors and neutral vacuum with declared multiplicities; tensor product in canonical CAR order, then retained routed family/spin frame grammar",
                    "failure_rule": "selector30/31 have no field/clock advance; fixed initial factors and both pointer outputs retained",
                    "sampler_normalizes_born_law_separately_from_source_arrays": True,
                    "simultaneous_sampling_certificate": sampling},
        "dispersion_from_records": recovered, "source_resources": resources,
        "resource_contract": "charged full successful schedule upper bound for every attempted trial, including aborts; source C/E services required, not executed",
        "frame_source_word_grammar": frames,
        "error_contract": {"binary_all_column_pair_error_admission": str(ADMISSION),
                           "maximum_binary_pair_error_upper": str(ADMISSION),
                           "maximum_binary_complete_score_bias_upper": str(observer.NUMERIC_SCORE_ADMISSION),
                           "backend_pair_error_diagnostic": str(observer._outward((F(), max_operator))[1]),
                           "backend_complete_score_bias_diagnostic": str(observer._outward((F(), max_bias))[1]),
                           "backend_diagnostics_authenticated_by_replay": False,
                           "complete_source_array_replay_tolerance_is_physical_error_bound": False,
                           "source_sampling_radius": str(observer.RECORD_RADIUS),
                           "physical_score_allowance_required": str(observer.PHYSICAL_SCORE_ALLOWANCE),
                           "separate_frame_score_allowance_conditional": str(observer.WEAK_FRAME_SCORE_ALLOWANCE),
                           "separate_sampler_rounding_allowance": str(observer.SAMPLER_SCORE_ALLOWANCE),
                           "physical_service_authentication": False},
        "scope": {"calibrated_finite_free3D_momentum_field_source_executed": True,
                  "record_energy_inversion_uses_mass_conditioned_probe_controls": False,
                  "empirical_mass_Yukawa_inputs_derived": False,
                  "native_CE_protection_history_or_physical_service_error_certified": False,
                  "local_interacting_continuum_or_physical_momentum_clock_dictionary_derived": False,
                  "quantum_gauge_Higgs_dynamics_or_confinement_executed": False}}
    return report, arrays


def symbol_intertwining():
    """All180 Nambu columns of the calibrated symbol in each source frame."""
    weak = np.asarray(matter.mass_kernel().evalf(20), complex)
    target_family = np.asarray(matter.mass_basis_change().evalf(20), complex)
    mass = target_family.T@weak@target_family
    family = family_source.source_frame()["one_body_basis_change90"]
    rows = []
    for magnitude in MOMENTUM_MAGNITUDES:
        for axis_number, axis in enumerate(("x", "y", "z")):
            spin = helicity_frame(axis)["spinor_unitary"]
            A = np.kron(spin, np.eye(45))@family
            transform = np.block([[A, np.zeros_like(A)], [np.zeros_like(A), A.conj()]])
            p = [0., 0., 0.]; p[axis_number] = float(magnitude)
            target = kinetic.bdg_symbol((0., 0., float(magnitude)), mass, speed=1.)
            full = matter.weyl_symbol(tuple(p))
            residual = float(np.linalg.norm(transform.conj().T@full@transform-target, 2))
            q = kinetic.nambu_charge(matter.matrices()["Q45"])
            charge = float(np.linalg.norm(transform.conj().T@q@transform-q, 2))
            if residual > 2e-10 or charge > 2e-10:
                raise ArithmeticError("actual source spin/family frame does not intertwine the calibrated moving action/current")
            rows.append({"axis": axis, "momentum_GeV": str(magnitude),
                         "complete_Nambu_columns": 180,
                         "source_frame_symbol_difference_diagnostic": residual,
                         "source_frame_electric_charge_difference_diagnostic": charge})
    return rows


def _route_word(local_word, sites):
    """Literal even two-mode operator routed with actual fermionic swaps."""
    a, b = sites
    if not 0 <= a < b:
        raise ValueError("ordered distinct source factor sites required")
    route = [(j-1, j) for j in range(b, a+1, -1)]
    swap = g.local_gate("fswap")
    word = []
    for pair in route:
        word += COMPILER.remap(swap.word, pair)
    word += COMPILER.remap(local_word, (a, a+1))
    for pair in reversed(route):
        word += COMPILER.remap(swap.word, pair)
    return word


def charged_momentum_pair(mass, magnitude, axis, duration):
    """Actual256-column8-CAR Dirac pair at two distinct opposite momenta.

    Order is (p spin-up L,Rc; p spin-down L,Rc; -p likewise). The four
    nested pair terms and spin-frame routing are independent source words.
    This explicit small Fock check anchors the factor theorem used for the
    larger calibrated family; Nambu doubling is not extra physical modes.
    """
    m, p, t = _exact(mass, nonnegative=True), _exact(magnitude, nonnegative=True), _time(duration)
    if p <= 0:
        raise ValueError("distinct opposite momenta require strictly positive magnitude")
    spin_word = helicity_frame(axis)["word"]
    spin_frame_word = []
    for pair in ((0, 2), (1, 3), (4, 6), (5, 7)):
        spin_frame_word += _route_word(spin_word, pair)
    mass_word = []
    for (a, b), e, sign in (((0, 7), -p, 1), ((1, 6), -p, 1),
                            ((2, 5), p, -1), ((3, 4), p, -1)):
        pair = moving_pair(m, e, t)
        local = list(pair["word"])
        if sign == -1:
            phase = g.local_gate("onsite", 1., math.pi).word
            local = list(phase)+local+COMPILER.inverse(phase)
        mass_word += _route_word(local, (a, b))
    word = COMPILER.inverse(spin_frame_word)+mass_word+spin_frame_word
    actual = _freeze(COMPILER.execute(word, 8))
    vector = [0., 0., 0.]; vector[("x", "y", "z").index(axis)] = float(p)
    kernel = np.array([[0., float(m)], [float(m), 0.]])
    h, delta = kinetic.momentum_pair_source(tuple(vector), kernel, speed=1.)
    target = expm(-1j*float(t)*protected._fock_hamiltonian(h, delta))
    residual = float(np.linalg.norm(actual-target, 2))
    if not math.isfinite(residual) or residual > 2e-10:
        raise ArithmeticError("complete moving8-CAR source differs from independent action target")
    return {"unitary": actual, "target": _freeze(target), "word": tuple(word),
            "resources": COMPILER.resources(word), "all_Fock_columns": 256,
            "physical_CAR_modes": 8, "opposite_momentum_is_distinct": True,
            "source_action_difference_diagnostic": residual,
            "global_phase_removed_or_fitted": False}


def source_bindings():
    """Closed code/data custody, independent of interpreter import order."""
    names = set(observer.source_bindings()) | set(bounds.source_bindings()) | {
        "model/common_source/calibrated_transport.py",
        "model/common_source/test_calibrated_transport.py"}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def canonical_report(report):
    """Check fixed admissions; do not demand backend residual-byte identity."""
    result = deepcopy(report)
    fields = ((result["error_contract"], "backend_pair_error_diagnostic", ADMISSION),
              (result["error_contract"], "backend_complete_score_bias_diagnostic", observer.NUMERIC_SCORE_ADMISSION),
              (result["frame_source_word_grammar"], "backend_FSWAP_error_diagnostic", F(1, 10**13)))
    for owner, key, admission in fields:
        value = owner[key]
        if type(value) is not str or not 0 <= F(value) <= admission:
            raise ValueError("retained numerical source diagnostic exceeds fixed admission: "+key)
        owner[key] = "backend diagnostic within fixed admission"
    return result


def compare_execution(report, arrays, fresh_report, fresh_arrays):
    """Exact schemas/counts and portable full complex-output comparison."""
    observer._same_report(canonical_report(fresh_report), canonical_report(report))
    if set(arrays) != set(fresh_arrays):
        raise ValueError("complete transport artifact array keys differ")
    for name, expected in fresh_arrays.items():
        actual = arrays[name]
        if type(actual) is not np.ndarray or actual.dtype != expected.dtype or actual.shape != expected.shape:
            raise ValueError("transport artifact array schema differs: "+name)
        if actual.dtype.kind in "Ui":
            if not np.array_equal(actual, expected):
                raise ValueError("complete transport records differ: "+name)
        elif (not np.isfinite(actual).all()
              or not np.allclose(actual, expected, rtol=2e-11, atol=2e-12)):
            raise ValueError("fresh complete complex transport source output differs: "+name)


def _finite_json(value):
    if type(value) is dict:
        for member in value.values():
            _finite_json(member)
    elif type(value) is list:
        for member in value:
            _finite_json(member)
    elif type(value) is float and not math.isfinite(value):
        raise ValueError("nonfinite retained transport JSON is inadmissible")


def check(path):
    """Own a fresh source replay driven by the exact retained full counts."""
    path = Path(path)
    report = strict_json.load(path)
    _finite_json(report)
    if type(report) is not dict or report.get("source_bindings") != source_bindings():
        raise ValueError("retained transport closed source bindings differ")
    artifact = path.with_suffix(".npz")
    if report.get("artifact") != {"filename": artifact.name,
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError("retained complete transport artifact custody differs")
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError("duplicate transport NPZ array names are inadmissible")
        arrays = {name: archive[name] for name in archive.files}
    counts = _counts(arrays.get("complete_record_counts"))
    fresh, fresh_arrays = execute(levels=counts.shape[1], seed=report["sampler_seed"], retained_counts=counts)
    scientific = {key: value for key, value in report.items() if key not in ("source_bindings", "artifact")}
    compare_execution(scientific, arrays, fresh, fresh_arrays)
    return scientific


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    parser.add_argument("--levels", type=int, default=DEFAULT_LEVELS)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print("calibrated finite moving-field complete-record/source replay: PASS")
        return
    artifact = args.output.with_suffix(".npz")
    if args.output.exists() or artifact.exists():
        parser.error("choose fresh JSON and complete artifact paths")
    initial_bindings = source_bindings()
    report, arrays = execute(levels=args.levels)
    if source_bindings() != initial_bindings:
        raise RuntimeError("transport sources changed during execution; receipt was not written")
    report["source_bindings"] = initial_bindings
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps([{key: row[key] for key in ("label", "momentum_GeV", "invariant_mass_interval_numeric_GeV")}
                      for row in report["dispersion_from_records"]]))


if __name__ == "__main__":
    main()
