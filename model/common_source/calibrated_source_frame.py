"""Phase-exact source implementation of the calibrated weak/mass frame.

The CKM Euler factors execute as actual number-conserving fermion source
words. Nonadjacent modes0,2 are routed by FSWAP, including the doubly
occupied sign. All eight columns and the vacuum phase are compared with
the independent exterior-power lift of the supplied calibrated CKM matrix.
Six disjoint down-left color/spin registers and the actual singlet signs
then specify the complete90-mode weak-to-mass frame as exact tensor factors.
This compiles a calibrated basis; it does not derive its empirical values.
"""
from __future__ import annotations

import hashlib
import math
from numbers import Real
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import calibrated_matter as matter
from .chiral_candidate import gaussian_source as g
from .chiral_candidate import protected_momentum as protected


ROOT = Path(__file__).resolve().parents[2]
ADMISSION = 2e-10


def _freeze(value):
    array = np.asarray(value, complex)
    return np.frombuffer(array.tobytes(), complex).reshape(array.shape)


def exterior(unitary):
    """Independent full Fock lift by minors, in little-endian CAR order."""
    value = np.asarray(unitary, complex)
    if value.shape != (3, 3) or not np.isfinite(value).all():
        raise ValueError("one finite3x3 one-particle matrix required")
    subsets = [tuple(j for j in range(3) if mask >> j & 1) for mask in range(8)]
    output = np.zeros((8, 8), complex)
    for i, left in enumerate(subsets):
        for j, right in enumerate(subsets):
            if len(left) == len(right):
                output[i, j] = np.linalg.det(value[np.ix_(left, right)]) if left else 1.
    return _freeze(output)


def _logical_word():
    packet = matter.parameters()["CKM"]
    angles = {ij: math.asin(float(packet["sin_theta"+ij])) for ij in ("12", "13", "23")}
    delta = float(packet["delta_radians"])
    calls = (("hopping", 1j, angles["12"], (0, 1)),
             ("fswap", 0., 0., (1, 2)),
             ("hopping", 1j*np.exp(-1j*delta), angles["13"], (0, 1)),
             ("fswap", 0., 0., (1, 2)),
             ("hopping", 1j, angles["23"], (1, 2)))
    word = []
    for kind, coefficient, duration, sites in calls:
        gate = g.local_gate(kind, coefficient, duration)
        word.extend(protected.compiler.remap(gate.word, sites))
    return protected.compiler._validate(word, 3)


def _global_word(word, modes):
    """Emit and validate the decoded occupations of the disjoint native codes."""
    emitted = protected.decoded_source_word(word, tuple(modes))
    lowered = protected._lower_decoded_word(emitted, tuple(modes))
    if lowered != protected.compiler._validate(word, len(modes)):
        raise ArithmeticError("calibrated frame changes decoded occupation/factor order")
    return emitted


def _sum_resources(*rows):
    return {name: sum(row[name] for row in rows) for name in g._ADDITIVE_RESOURCES}


def source_frame():
    """Complete source CKM factor and disjoint full90-mode frame mapping.

    Arrays have immutable byte backing. All words and mutable metadata are
    fresh objects on every call, including after source-gate caches warm.
    `resources` counts one complete90-mode basis change; `sourceword` is
    the single three-mode CKM factor, and `decoded_sourceword` is the full
    field word on decoded native occupation wires9*j.
    """
    word = _logical_word()
    actual = _freeze(protected.compiler.execute(word, 3))
    ckm = np.asarray(matter.ckm_matrix().evalf(20), complex)
    target = exterior(ckm)
    residual = float(np.linalg.norm(actual-target, 2))
    unitarity = float(np.linalg.norm(actual.conj().T@actual-np.eye(8), 2))
    if (not math.isfinite(residual) or residual > ADMISSION
            or not math.isfinite(unitarity) or unitarity > ADMISSION):
        raise ArithmeticError("source CKM Fock lift differs from phase-exact exterior target")
    # The actual same-source sign pulse is not a freely discarded global
    # phase: its vacuum column remains1 and occupied column is-1.
    phase_word = list(g.local_gate("onsite", 1., math.pi).word)
    phase = _freeze(protected.compiler.execute(phase_word, 1))
    phase_error = float(np.linalg.norm(phase-np.diag([1., -1.]), 2))
    if not math.isfinite(phase_error) or phase_error > ADMISSION:
        raise ArithmeticError("source singlet sign has a wrong vacuum/occupied phase")
    source_one = np.asarray(actual[np.ix_((1, 2, 4), (1, 2, 4))])
    target45 = np.asarray(matter.mass_basis_change().evalf(20), complex)
    target90 = np.kron(np.eye(2), target45)
    onebody = np.eye(90, dtype=complex)
    blocks, signs, emitted = [], [], []
    covered = set()
    for spin in (0, 1):
        for internal in (3, 6, 8):
            modes = tuple(45*spin+3*internal+family for family in range(3))
            onebody[np.ix_(modes, modes)] = source_one
            block_word = _global_word(word, modes)
            # Check this literal's action even after the local source word
            # and numerical native contractions have been cached.
            lowered = protected._lower_decoded_word(block_word, modes)
            if np.max(abs(protected.compiler.execute(lowered, 3)-actual)) > ADMISSION:
                raise ArithmeticError("emitted full-field CKM block changes source action")
            emitted.extend(block_word)
            blocks.append({"spin": spin, "down_left_internal_component": internal,
                           "CAR_indices": list(modes), "decoded_occupation_wires": [9*j for j in modes]})
            if covered.intersection(modes):
                raise ArithmeticError("calibrated CKM factors overlap")
            covered.update(modes)
        for internal in range(15):
            indices = tuple(45*spin+3*internal+family for family in range(3))
            if not np.array_equal(target90[np.ix_(indices, indices)], -np.eye(3)):
                continue
            for mode in indices:
                if mode in covered:
                    raise ArithmeticError("a singlet phase overlaps the CKM factor")
                phase_global = _global_word(phase_word, (mode,))
                emitted.extend(phase_global)
                onebody[mode, mode] = phase[1, 1]
                signs.append({"CAR_index": mode, "decoded_occupation_wire": 9*mode,
                              "phase": "-1", "internal_component": internal, "spin": spin})
                covered.add(mode)
    error90 = float(np.linalg.norm(onebody-target90, 2))
    if len(blocks) != 6 or len(signs) != 18 or not math.isfinite(error90) or error90 > ADMISSION:
        raise ArithmeticError("actual factor mapping differs from the complete calibrated M15xSpin frame")
    block_resources = protected.compiler.resources(word)
    phase_resources = protected.compiler.resources(phase_word)
    additive_resources = _sum_resources(*([block_resources]*6+[phase_resources]*18))
    emitted_resources = protected.compiler.resources(emitted)
    if not protected._same_resources(additive_resources, emitted_resources):
        raise ArithmeticError("global source basis-change word has inconsistent resource counts")
    return {"CKM_fock_unitary": actual, "CKM_fock_target": target,
            "singlet_phase_unitary": phase, "one_body_basis_change90": _freeze(onebody),
            "sourceword": tuple(word), "decoded_sourceword": tuple(emitted),
            "resources": emitted_resources, "CKM_block_resources": block_resources,
            "singlet_phase_resources": phase_resources,
            "mapping": {"basis_convention": "chi_weak=W chi_mass; creation-Fock lift Gamma(W)",
                        "CKM_blocks": blocks, "negative_singlet_phases": signs,
                        "identity_CAR_indices": sorted(set(range(90))-covered),
                        "all90_CAR_modes_retained": True,
                        "full_Fock_representation": "six disjoint complete8-column CKM factors,18 complete2-column phases,and identity on all other modes",
                        "native_occupation_wire_rule": "CAR modej -> decoded occupation wire9*j"},
            "diagnostics": {"CKM_source_operator_difference": residual,
                            "CKM_source_unitarity_residual": unitarity,
                            "singlet_phase_source_difference": phase_error,
                            "full90_mode_onebody_difference": error90,
                            "vacuum_phase_difference": float(abs(actual[0, 0]-1)),
                            "full_CKM_Fock_columns_executed": 8,
                            "numerical_operator_admission": ADMISSION,
                            "global_phase_removed_or_fitted": False,
                            "polar_norm_or_probability_repair_applied": False,
                            "physical_native_service_error_certified": False}}


def source_bindings():
    paths = set(matter.source_bindings()) | {
        "model/common_source/calibrated_source_frame.py",
        "model/common_source/test_calibrated_source_frame.py",
        "model/common_source/chiral_candidate/gaussian_source.py",
        "model/common_source/chiral_candidate/protected_momentum.py",
        "model/observer/source_scalar/compiler.py"}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def source_defined_down_sector(angles, *, singlet_sign=1):
    """Actual64-column weak-frame propagation of three charged rest pairs.

    Local CAR order is left families0..2, conjugated-right families3..5.
    Pair(f,f+3) rotations commute but each literal uses its full FSWAP
    route. The chronological word is inverse-frame, mass evolution, frame.
    Three opposite-charge all-left mass pairs therefore implement the weak
    creation coefficient W6 Delta_mass W6.T, independently checked here.

    The numerical all-input admissions are not outward arithmetic bounds
    or physical service certificates. A caller propagating them must state
    that conditional numerical contract and retain the complete frame word.
    """
    if (type(angles) is not tuple or len(angles) != 3
            or any(isinstance(x, (bool, np.bool_)) or not isinstance(x, Real)
                   or not math.isfinite(float(x)) or abs(float(x)) > 4 for x in angles)):
        raise ValueError("three finite real reduced pair angles with absolute value<=4 required")
    if type(singlet_sign) is not int or singlet_sign not in (-1, 1):
        raise ValueError("source singlet sign+1 or-1 required")
    angles = tuple(map(float, angles))
    compiler = protected.compiler
    frame_word = _logical_word()
    if singlet_sign == -1:
        phase_word = g.local_gate("onsite", 1., math.pi).word
        for site in (3, 4, 5):
            frame_word += compiler.remap(phase_word, (site,))
    inverse_word = compiler.inverse(frame_word)
    mass_word = []
    for family, angle in enumerate(angles):
        a, b = family, family+3
        route = [(j-1, j) for j in range(b, a+1, -1)]
        swap = g.local_gate("fswap")
        pair = g.local_gate("pairing", 1., angle)
        for sites in route:
            mass_word += compiler.remap(swap.word, sites)
        mass_word += compiler.remap(pair.word, (a, a+1))
        for sites in reversed(route):
            mass_word += compiler.remap(swap.word, sites)
    full_word = inverse_word+mass_word+frame_word
    # Execute every input, including odd sectors and absolute vacuum phase.
    actual_frame = compiler.execute(frame_word, 6)
    actual_inverse = compiler.execute(inverse_word, 6)
    actual_mass = compiler.execute(mass_word, 6)
    actual_weak = compiler.execute(full_word, 6)
    ckm = np.asarray(matter.ckm_matrix().evalf(20), complex)
    w6 = np.zeros((6, 6), complex)
    w6[:3, :3], w6[3:, 3:] = ckm, singlet_sign*np.eye(3)
    delta = np.zeros((6, 6), complex)
    for family, angle in enumerate(angles):
        delta[family, family+3], delta[family+3, family] = angle, -angle
    target_mass = expm(-1j*protected._fock_hamiltonian(np.zeros((6, 6)), delta))
    target_weak = expm(-1j*protected._fock_hamiltonian(np.zeros((6, 6)), w6@delta@w6.T))
    expected_frame = np.kron(np.eye(8), exterior(ckm))
    if singlet_sign == -1:
        right_sign = np.diag([(-1.)**mask.bit_count() for mask in range(8)])
        expected_frame = np.kron(right_sign, exterior(ckm))
    diagnostics = {
        "all_input_mass_source_target_difference": float(np.linalg.norm(actual_mass-target_mass, 2)),
        "all_input_weak_source_target_difference": float(np.linalg.norm(actual_weak-target_weak, 2)),
        "frame_source_target_difference": float(np.linalg.norm(actual_frame-expected_frame, 2)),
        "inverse_frame_source_target_difference": float(np.linalg.norm(actual_inverse-expected_frame.conj().T, 2)),
        "literal_composition_difference": float(np.linalg.norm(actual_weak-actual_frame@actual_mass@actual_inverse, 2)),
        "weak_source_unitarity_residual": float(np.linalg.norm(actual_weak.conj().T@actual_weak-np.eye(64), 2)),
    }
    if any(not math.isfinite(x) or x > ADMISSION for x in diagnostics.values()):
        raise ArithmeticError("complete source-defined weak/mass conjugation failed numerical all-input admission")
    literal = _global_word(full_word, tuple(range(6)))
    return {"weak_unitary": _freeze(actual_weak), "mass_pair_product": _freeze(actual_mass),
            "frame_unitary": _freeze(actual_frame), "inverse_frame_unitary": _freeze(actual_inverse),
            "weak_target": _freeze(target_weak), "mass_target": _freeze(target_mass),
            "resources": compiler.resources(full_word), "sourceword": tuple(full_word),
            "decoded_sourceword": tuple(literal), "errors": diagnostics,
            "numerical_operator_admission": ADMISSION,
            "numerical_admission_is_outward_or_physical_certificate": False,
            "source_all_input_columns": 64, "source_global_phase_removed_or_fitted": False,
            "factor_local_CAR_order": "left families0,1,2; conjugated-right families3,4,5",
            "singlet_sign": singlet_sign,
            "full_field_embedding": "one declared down-quark color/spin block of the complete90-mode source frame"}
