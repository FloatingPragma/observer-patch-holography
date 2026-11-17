"""Charge-conserving record tomography of a finite charged CAR response.

This is a new apparatus protocol for the existing supplied quadratic action.
It reconstructs a two-mode charge block from detector probabilities and a
programmed clock, not from an eigenvalue oracle. The source action, reference
fermion, preparation, spatial/Higgs coordinates and alias window are premises.
It establishes neither a selected mass law nor an interacting physical field.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm

from tools import strict_json
from . import action, history, protected_momentum as protected

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INTERVAL = F(1, 4)
SECOND_INTERVAL = F(3, 8)
PHASES = {"X+": 1., "X-": -1., "Y+": 1j, "Y-": -1j}
TRAINING = ((0., 1.), (-1/3, 1.), (1/3, 1.), (0., .75), (0., 1.25))
HELD_OUT = ((2/3, .9), (-.5, 1.1), (1/6, 1.4))
ADMISSION = 2e-8


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, F)):
        raise ValueError(name+" must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(name+" must be finite")
    return value


def _pack(matrix):
    a = np.asarray(matrix, complex)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def _unpack(value):
    a = np.asarray(value["real"], float)+1j*np.asarray(value["imag"], float)
    if a.shape != (2, 2) or not np.isfinite(a).all():
        raise ValueError("a finite two-by-two response is required")
    return a


def _freeze(value):
    value = np.asarray(value)
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


@lru_cache(maxsize=1)
def _apparatus():
    """Actual neutral source links in a fixed-total-charge reference code.

    d=(c_-,c_+^dagger,c_ref), with physical charges(-1,+1,-1).
    The d vacuum is physical |010>; its one-d-particle code has Q=0.
    The basis signs below follow the bitwise CAR ordering, not a bosonic swap.
    """
    basis = np.eye(8, dtype=complex)[:, [3, 0, 6]]
    basis[:, 2] *= -1
    charge = np.diag([-((i >> 0) & 1)+((i >> 1) & 1)-((i >> 2) & 1)
                      for i in range(8)])
    assert np.max(abs(charge@basis)) == 0
    rotations, resources = {}, {}
    for site, phase_name in itertools.product(range(2), PHASES):
        z = 1j*PHASES[phase_name]
        h, delta = np.zeros((3, 3), complex), np.zeros((3, 3), complex)
        if site == 0:
            h[0, 2], h[2, 0] = z, z.conjugate()
        else:
            # z*d_1^dagger*d_ref+h.c. is neutral physical pair creation.
            delta[1, 2], delta[2, 1] = -z.conjugate(), z.conjugate()
        gate = protected._small_evolution(h, delta, math.pi/4, 1)
        reverse = protected._small_evolution(h, delta, -math.pi/4, 1)
        expected = (basis[:, 2]+PHASES[phase_name]*basis[:, site])/math.sqrt(2)
        if (np.max(abs(gate.unitary@basis[:, 2]-expected)) > 2e-10
                or np.max(abs(reverse.unitary@gate.unitary-np.eye(8))) > 2e-10
                or np.max(abs(charge@gate.unitary-gate.unitary@charge)) > 2e-10):
            raise ArithmeticError("source reference rotation lost its CAR sign, inverse or charge")
        rotations[site, phase_name] = (_freeze(gate.unitary), _freeze(reverse.unitary))
        resources[site, phase_name] = dict(gate.resources)
    return _freeze(basis), _freeze(charge), rotations, resources


def apparatus():
    basis, charge, rotations, resources = _apparatus()
    return basis, charge, dict(rotations), {key: dict(value) for key, value in resources.items()}


def source_pair(momentum, higgs_scale, *, law="continuum"):
    """An exactly decoupled two-mode factor of the existing z-axis field.

    The opposite-momentum field remains the owner of species multiplicities.
    This component does not interpret two Nambu entries as two species.
    """
    p, scale = _real(momentum, "coordinate momentum"), _real(higgs_scale, "Higgs multiplier")
    if abs(p) > 1 or not .5 <= scale <= 1.5:
        raise ValueError("the declared coordinate/Higgs intervention domain was exceeded")
    if law == "continuum":
        effective = p
    elif law == "two_hop_control":
        total = 2/(3*math.sqrt(3))
        effective = (2*total-1/math.pi)*math.sin(math.pi*p)+(1/math.pi-total)*math.sin(2*math.pi*p)
    else:
        raise ValueError("unknown source kinetic control")
    full_h, full_delta, layout = protected.field_layout((0., 0., effective))
    partners = np.flatnonzero(full_delta[0] != 0)
    if len(partners) != 1:
        raise ArithmeticError("selected charged source mode does not have a unique paired factor")
    indices = (0, int(partners[0]))
    original = [layout["original_fermionic_wire_order"][i] for i in indices]
    labels = action.mode_labels()
    charges = [F(labels[i % 90]["electric_Q"]) for i in original]
    if charges != [-1, 1]:
        raise ArithmeticError("selected field factor does not have the declared opposite physical charges")
    outside = [j for j in range(len(full_h)) if j not in indices]
    if (np.any(full_h[np.ix_(indices, outside)]) or np.any(full_delta[np.ix_(indices, outside)])):
        raise ArithmeticError("the proposed component is not a decoupled source factor")
    h = full_h[np.ix_(indices, indices)]
    delta = scale*full_delta[np.ix_(indices, indices)]
    return h, delta, {"source_CAR_modes": layout["CAR_modes"], "source_factor_indices": list(indices),
        "original_source_Fock_indices": original, "source_factor_electric_charges": list(map(str, charges)),
        "selected_internal_types": [labels[i % 90]["internal_type"] for i in original],
        "native_family_axes": [labels[i % 90]["family"] for i in original],
        "source_factor_decouples_exactly": True, "coordinate_momentum": p,
        "neutral_Higgs_multiplier": scale, "kinetic_law": law}


def source_unitary(momentum, higgs_scale, *, steps=64, law="continuum", duration=INTERVAL):
    h, delta, domain = source_pair(momentum, higgs_scale, law=law)
    work = protected._small_evolution(h, delta, float(duration), steps)
    joint = np.zeros((8, 8), complex)
    joint[:4, :4] = joint[4:, 4:] = work.unitary
    _, charge, _, _ = apparatus()
    if np.max(abs(charge@joint-joint@charge)) > 2e-10:
        raise ArithmeticError("source propagation does not conserve the declared charge")
    return joint, {**domain, "steps": steps, "source_resources_per_interval": work.resources,
        "source_to_continuous_operator_error": work.source_error}


def _read(vector, site):
    return protected._number_vectors(vector, site)


def measure(unitary, *, duration=INTERVAL):
    """Prepare once, read a neutral quadrature, continue, and read again.

    Each of 16 settings retains 16 outcomes, including detector complements.
    The clock purifiers and unnormalized work states are retained separately.
    No probability or post-measurement state is renormalized or repaired.
    """
    basis, charge, rotations, resources = apparatus()
    unitary = np.asarray(unitary, complex)
    if unitary.shape != (8, 8) or not np.isfinite(unitary).all():
        raise ValueError("a three-mode work/reference source is required")
    if np.linalg.norm(unitary.conj().T@unitary-np.eye(8), 2) > ADMISSION:
        raise ValueError("the apparatus requires a unitary source on the complete work/reference space")
    clock, calibration, clock_cost = history._clock_preparation()
    clocks = protected._clock_outputs(clock, 2, float(duration))
    vacuum = np.eye(8, dtype=complex)[:, 0]
    # Fill the positive work mode and negative reference together, using two
    # actual neutral inverse quarter rotations. Blank initialization is the
    # supplied service; no charged reference state is injected.
    reference = rotations[1, "X+"][1]@rotations[1, "X+"][1]@vacuum
    if np.linalg.norm(reference-basis[:, 2]) > 2e-10:
        raise ArithmeticError("neutral source reference preparation failed")
    contexts, arrays = [], {}
    maximum_charge = 0.
    for i, j, phase in itertools.product(range(2), range(2), PHASES):
        key = f"i{i}_j{j}_{phase}"
        state = rotations[j, "X+"][0]@reference
        state = unitary@state
        forward, reverse = rotations[i, phase]
        first = tuple(forward@v for v in _read(reverse@state, 2))
        terminal, outputs = [], []
        for bit, branch in enumerate(first):
            # An observer's first read changes the work state used next.
            for continued, vector in enumerate(_read(unitary@branch, 0)):
                maximum_charge = max(maximum_charge, float(np.linalg.norm(charge@vector)))
                for (cx, cy), clock_columns in clocks.items():
                    probability = float(np.vdot(vector, vector).real*np.vdot(clock_columns, clock_columns).real)
                    terminal.append({"first": bit, "continued": continued, "clock_X": cx, "clock_Y": cy,
                                     "probability": probability})
                    outputs.append(vector)
        contexts.append({"output": i, "input": j, "phase": phase, "terminal": terminal})
        arrays[key+"_work"] = np.asarray(outputs)
    arrays["clock_reference_columns"] = np.asarray(list(clocks.values()))
    if maximum_charge > 2e-10:
        raise ArithmeticError("a retained detector branch leaves total charge zero")
    packet = {"schema": "oph.charged_record_two_point.v1", "contexts": contexts,
        "clock_calibration": calibration,
        "instrument_contract": {"physical_mode_charges": [-1, 1, -1], "fixed_total_charge": 0,
            "charge_embedding": "supplied model electric charge on physical CAR modes; no empirical identification",
            "clock_reference_charge_representation": "trivial neutral apparatus premise",
            "charged_external_phase_reference_required": False, "reference_fermion_and_neutral_phase_controls_supplied": True,
            "work_preparations_per_context": 1, "work_resets_per_context": 0,
            "complete_detector_complements_retained": True, "clock_purifiers_retained": True,
            "source_preparation_native_pulses_by_input": {
                str(j): 2*resources[1, "X+"]["native_pulses"]+resources[j, "X+"]["native_pulses"] for j in range(2)},
            "reference_prepared_from_supplied_blank_by_neutral_source_pair_pulses": True,
            "clock_preparation_resources": clock_cost,
            "source_detector_rotation_pulses": {f"{i}_{phase}": r["native_pulses"] for (i, phase), r in resources.items()},
            "detector_basis_rotations_per_context": 2,
            "source_propagation_intervals_per_context": 2,
            "resource_entries_are_component_counts_not_total_maintenance_cost": True,
            "parity_and_charge_checked_on_actual_source_operators": True,
            "domain": "completed compiled coherent source service and complete detector instrument",
            "native_attempt_abort_and_export_failure_tree_included_here": False,
            "native_service_success_probabilities_not_set_to_one_by_this_protocol": True},
        "branch_contract": {"generator_eigenphases_strictly_inside_principal_interval": True,
            "maximum_absolute_eigenphase": math.pi/2, "clock_two_intervals_before_final_read": True,
            "this_prior_is_not_inferred_from_single_time_data": True}}
    return packet, arrays


def decode(packet, *, zero_time=False):
    """Blind complex transfer inverse using complete detector records only."""
    if packet.get("schema") != "oph.charged_record_two_point.v1":
        raise ValueError("wrong record schema")
    branch = packet.get("branch_contract", {})
    if (branch.get("generator_eigenphases_strictly_inside_principal_interval") is not True
            or branch.get("clock_two_intervals_before_final_read") is not True):
        raise ValueError("an independently supplied generator/clock branch contract is required")
    phase_bound = _real(branch["maximum_absolute_eigenphase"], "generator phase bound")
    if not 0 < phase_bound < math.pi:
        raise ValueError("a strict independently supplied principal logarithm window is required")
    contrast = _real(packet["clock_calibration"]["contrast"], "clock contrast")
    if not 0 < abs(contrast) <= 1:
        raise ValueError("invalid calibrated clock contrast")
    expected = set(itertools.product(range(2), range(2), PHASES))
    probabilities, angles = {}, []
    for context in packet["contexts"]:
        key = context["output"], context["input"], context["phase"]
        if key not in expected or key in probabilities:
            raise ValueError("missing, duplicate or unknown tomographic context")
        rows = context["terminal"]
        labels = [(r["first"], r["continued"], r["clock_X"], r["clock_Y"]) for r in rows]
        if len(labels) != 16 or set(labels) != set(itertools.product((0, 1), repeat=4)):
            raise ValueError("complete uniquely labelled detector/continuation records required")
        weights = [_real(r["probability"], "record probability") for r in rows]
        if min(weights) < 0 or max(weights) > 1+ADMISSION or abs(math.fsum(weights)-1) > ADMISSION:
            raise ValueError("complete unnormalized instrument probabilities required")
        probabilities[key] = math.fsum(w for r, w in zip(rows, weights) if r["first"] == 1)
        means = [math.fsum((1 if r[name] == 0 else -1)*w for r, w in zip(rows, weights))/contrast
                 for name in ("clock_X", "clock_Y")]
        if abs(math.hypot(*means)-1) > ADMISSION:
            raise ValueError("clock records leave the calibrated circle")
        angles.append(math.atan2(means[1], means[0]))
    if set(probabilities) != expected or max(angles)-min(angles) > ADMISSION:
        raise ValueError("incomplete tomography or incompatible same-process clock")
    transfer = np.array([[probabilities[i, j, "X+"]-probabilities[i, j, "X-"]
                         +1j*(probabilities[i, j, "Y+"]-probabilities[i, j, "Y-"])
                         for j in range(2)] for i in range(2)])
    defect = float(np.linalg.norm(transfer.conj().T@transfer-np.eye(2), 2))
    if defect > ADMISSION:
        raise ValueError("records do not identify a closed canonical two-mode CAR transfer")
    result = {"transfer": _pack(transfer), "unitarity_defect": defect,
        "decoder_reads_source_action_or_simulator_duration": False,
        "finite_shot_or_outward_uncertainty_certified": False}
    if zero_time:
        if abs(angles[0]) > ADMISSION or np.linalg.norm(transfer-np.eye(2), 2) > ADMISSION:
            raise ValueError("equal-time interpolator calibration is not identity")
        return result
    if not 0 < angles[0] < math.pi:
        raise ValueError("clock phase outside the declared positive branch")
    time = angles[0]/2
    if max(abs(np.angle(np.linalg.eigvals(transfer)))) >= phase_bound:
        raise ValueError("generator logarithm exceeds the declared phase window")
    generator = 1j*logm(transfer)/time
    if np.linalg.norm(generator-generator.conj().T, 2) > ADMISSION:
        raise ValueError("reconstructed CAR generator is not Hermitian")
    values, vectors = np.linalg.eigh(generator)
    if values[1]-values[0] <= ADMISSION:
        raise ValueError("unresolved poles require whole-subspace residues, outside this simple-pole inverse")
    residues = [np.outer(vectors[:, i], vectors[:, i].conj()) for i in range(2)]
    result.update({"clock_interval": time, "generator_over_clock": _pack(generator),
        "poles_over_clock": values.tolist(), "canonical_residues": list(map(_pack, residues)),
        "residue_sum_error": float(np.linalg.norm(sum(residues)-np.eye(2), 2)),
        "inverse_response": "omega*I - B_record; retarded resolvent at omega+i*eta, eta>0",
        "canonical_interpolator_and_quadratic_closure_premises": True,
        "simple_poles_separated_above_numerical_admission": True,
        "generator_identified_only_within_supplied_alias_window": True,
        "interacting_wavefunction_renormalization_recovered": False})
    return result


def fit_joint_response(settings, decoded):
    design = np.array([[1., _real(p, "momentum"), _real(h, "Higgs multiplier")] for p, h in settings])
    if len(settings) != len(decoded) or np.linalg.matrix_rank(design) != 3:
        raise ValueError("independent momentum and Higgs interventions are required")
    matrices = np.asarray([_unpack(row["generator_over_clock"]) for row in decoded])
    coefficients = np.linalg.lstsq(design, matrices.reshape(len(settings), 4), rcond=None)[0].reshape(3, 2, 2)
    return {"constant": _pack(coefficients[0]), "momentum_derivative": _pack(coefficients[1]),
        "Higgs_multiplier_derivative": _pack(coefficients[2]),
        "maximum_training_operator_residual": float(max(np.linalg.norm(a-b, 2) for a, b in
            zip(matrices, np.einsum('sa,aij->sij', design, coefficients)))),
        "scope": "z-axis kinetic and neutral-Higgs-multiplier derivatives of one canonical finite charge block"}


def predict_generator(fit, setting):
    p, h = setting
    return _unpack(fit["constant"])+p*_unpack(fit["momentum_derivative"])+h*_unpack(fit["Higgs_multiplier_derivative"])


def embedded_transfer(transfer):
    basis = apparatus()[0][:, :2]
    return np.eye(8)+basis@(transfer-np.eye(2))@basis.conj().T


def validate_outputs(packet, arrays):
    """Tie each probability to its full unnormalized output, without conditioning."""
    clocks = np.asarray(arrays["clock_reference_columns"])
    if clocks.shape != (4, 4, 4) or not np.isfinite(clocks).all():
        raise ValueError("complete finite clock purifier columns are required")
    charge = apparatus()[1]
    for context in packet["contexts"]:
        key = f'i{context["output"]}_j{context["input"]}_{context["phase"]}_work'
        work = np.asarray(arrays[key])
        if work.shape != (16, 8) or not np.isfinite(work).all():
            raise ValueError("complete finite unnormalized continuation vectors are required")
        for row, vector in zip(context["terminal"], work, strict=True):
            columns = clocks[2*row["clock_X"]+row["clock_Y"]]
            probability = float(np.vdot(vector, vector).real*np.vdot(columns, columns).real)
            if (abs(probability-row["probability"]) > ADMISSION
                    or np.linalg.norm(charge@vector) > ADMISSION):
                raise ValueError("record probability/state or retained charge mismatch")


def complete_output_distance(left, a, right, b):
    """Largest context CQ trace distance, including measured continuation."""
    validate_outputs(left, a)
    validate_outputs(right, b)
    maximum = 0.
    for x, y in zip(left["contexts"], right["contexts"], strict=True):
        if (x["output"], x["input"], x["phase"]) != (y["output"], y["input"], y["phase"]):
            raise ValueError("incompatible record contexts")
        key = f'i{x["output"]}_j{x["input"]}_{x["phase"]}_work'
        distance = 0.
        for index, (v, w) in enumerate(zip(a[key], b[key], strict=True)):
            ca, cb = a["clock_reference_columns"][index % 4], b["clock_reference_columns"][index % 4]
            mass_a, mass_b = np.vdot(v, v).real*np.vdot(ca, ca).real, np.vdot(w, w).real*np.vdot(cb, cb).real
            # Orthogonal projection avoids subtracting two nearly equal fourth
            # powers when the retained and predicted branches almost coincide.
            av, bv = np.kron(v, ca.ravel()), np.kron(w, cb.ravel())
            if mass_a > 0:
                perpendicular = bv-av*(np.vdot(av, bv)/mass_a)
                radicand = float((mass_a-mass_b)**2+4*mass_a*np.vdot(perpendicular, perpendicular).real)
                distance += math.sqrt(radicand)/2
            else:
                distance += mass_b/2
        maximum = max(maximum, distance)
    return maximum


def ambiguity_controls(packets, decoded, *, steps):
    """Competing supplied actions expose exactly what the records identify."""
    control_packets, control_arrays, training_errors = {}, {}, []
    for index, setting in enumerate((*TRAINING, HELD_OUT[0])):
        unitary, _ = source_unitary(*setting, steps=steps, law="two_hop_control")
        packet, arrays = measure(unitary)
        key = "two_hop_"+str(index)
        control_packets[key] = packet
        control_arrays.update({key+"_"+name: value for name, value in arrays.items()})
        recovered = decode(packet)
        error = float(np.linalg.norm(_unpack(recovered["generator_over_clock"])
                                    -_unpack(decoded[str(index)]["generator_over_clock"]), 2))
        if index < len(TRAINING):
            training_errors.append(error)
        else:
            separation = error
    h, delta, _ = source_pair(*TRAINING[0])
    negative = protected._small_evolution(h, -delta, float(INTERVAL), steps).unitary
    positive = protected._small_evolution(h, delta, float(INTERVAL), steps).unitary
    packet, arrays = measure(np.kron(np.eye(2), negative))
    control_packets["negative_pair"] = packet
    control_arrays.update({"negative_pair_"+name: value for name, value in arrays.items()})
    negative_decoded = decode(packet)
    # These are native occupation reads of the empty pair, not eigenvalue data.
    marginal = [[float(np.vdot(v, v).real) for v in _read(u[:, 0], 0)] for u in (positive, negative)]
    reference = decoded["0"]
    generator = _unpack(reference["generator_over_clock"])
    time = reference["clock_interval"]
    projector = _unpack(reference["canonical_residues"][0])
    aliased = generator+4*math.pi/time*projector
    alias_error = max(float(np.linalg.norm(expm(-1j*t*generator)-expm(-1j*t*aliased), 2))
                      for t in (time, 1.5*time))
    return {"two_hop_same_training_maximum_generator_difference": max(training_errors),
        "two_hop_held_out_generator_separation": separation,
        "opposite_pair_sign_empty_pair_number_probabilities": marginal,
        "opposite_pair_sign_number_marginal_difference": float(np.max(abs(np.diff(marginal, axis=0)))),
        "opposite_pair_sign_record_generator_separation": float(np.linalg.norm(
            generator-_unpack(negative_decoded["generator_over_clock"]), 2)),
        "unbounded_alias_shift_over_clock": 4*math.pi/time,
        "unbounded_alias_two_duration_operator_difference": alias_error,
        "unbounded_alias_is_matrix_countermodel_not_native_source_execution": True,
        "unbounded_alias_excluded_only_by_supplied_generator_prior": True}, control_packets, control_arrays


def execute(*, steps=64):
    if isinstance(steps, bool) or not isinstance(steps, int) or steps <= 0 or steps % 2:
        raise ValueError("an even positive source step count is required for the fixed-step second duration")
    zero, zero_arrays = measure(np.eye(8), duration=F(0))
    decode(zero, zero_time=True)
    packets, decoded, domains = {}, {}, {}
    arrays = {"calibration_"+name: value for name, value in zero_arrays.items()}
    for index, setting in enumerate((*TRAINING, *HELD_OUT)):
        unitary, domain = source_unitary(*setting, steps=steps)
        packet, work = measure(unitary)
        key = str(index)
        packets[key], decoded[key], domains[key] = packet, decode(packet), domain
        arrays.update({key+"_"+name: value for name, value in work.items()})
    fit = fit_joint_response(TRAINING, [decoded[str(i)] for i in range(len(TRAINING))])
    validation, second_records, second_decoded = [], {}, {}
    for index, setting in enumerate(HELD_OUT, len(TRAINING)):
        key = str(index)
        predicted = predict_generator(fit, setting)
        time = decoded[key]["clock_interval"]
        predicted_packet, predicted_arrays = measure(embedded_transfer(expm(-1j*time*predicted)))
        original = {name[len(key)+1:]: value for name, value in arrays.items() if name.startswith(key+"_")}
        validation.append({"setting": list(setting),
            "operator_residual": float(np.linalg.norm(predicted-_unpack(decoded[key]["generator_over_clock"]), 2)),
            "complete_continued_record_quantum_output_distance": complete_output_distance(packets[key], original, predicted_packet, predicted_arrays)})
        second_unitary, second_domain = source_unitary(*setting, steps=3*steps//2, duration=SECOND_INTERVAL)
        first_unitary, _ = source_unitary(*setting, steps=steps)
        composition_error = float(np.linalg.norm(np.linalg.matrix_power(first_unitary, 3)
                                                -np.linalg.matrix_power(second_unitary, 2), 2))
        second_packet, second_arrays = measure(second_unitary, duration=SECOND_INTERVAL)
        recovered = decode(second_packet)
        second_records[key], second_decoded[key] = second_packet, recovered
        arrays.update({"second_"+key+"_"+name: value for name, value in second_arrays.items()})
        predicted_second, predicted_second_arrays = measure(
            embedded_transfer(expm(-1j*recovered["clock_interval"]*predicted)), duration=SECOND_INTERVAL)
        validation[-1].update({"second_duration_clock_interval": recovered["clock_interval"],
            "second_duration_operator_residual": float(np.linalg.norm(predicted-_unpack(recovered["generator_over_clock"]), 2)),
            "second_duration_complete_continued_output_distance": complete_output_distance(
                second_packet, second_arrays, predicted_second, predicted_second_arrays),
            "second_duration_source_domain": second_domain})
        validation[-1]["same_elementary_word_three_short_vs_two_long_operator_error"] = composition_error
    controls, control_packets, control_arrays = ambiguity_controls(packets, decoded, steps=steps)
    arrays.update({"control_"+name: value for name, value in control_arrays.items()})
    return {"schema": "oph.record_two_point_joint_response.v1", "source_steps": steps,
        "training_settings": list(map(list, TRAINING)), "held_out_settings": list(map(list, HELD_OUT)),
        "equal_time_calibration": zero, "domains": domains, "records": packets, "decoded": decoded,
        "held_out_second_duration_records": second_records, "held_out_second_duration_decoded": second_decoded,
        "joint_response": fit, "held_out_validation": validation,
        "ambiguity_controls": controls, "control_records": control_packets,
        "source_interval_relation": {"short_duration": str(INTERVAL), "long_duration": str(SECOND_INTERVAL),
            "short_steps": steps, "long_steps": 3*steps//2,
            "same_elementary_source_word_per_setting": True,
            "two_commensurate_times_do_not_remove_unbounded_generator_aliases": True},
        "scope": {"actual_neutral_source_probes_and_continuation_executed": True,
            "charge_homogeneous_Nambu_block": "(c_minus,c_plus_dagger)",
            "same_supplied_candidate_component_finite_generator_and_joint_derivatives_recovered": True,
            "source_mass_coefficient_supplied_not_selected_by_repair_maintenance_cost": True,
            "dynamical_Higgs_or_current_vertex_reconstructed": False,
            "new_particle_candidate_introduced": False,
            "source_action_clock_and_apparatus_supplied": True,
            "canonical_response_not_interacting_physical_propagator": True,
            "unmeasured_generator_branch_prior_required": True,
            "empirical_mass_fit": False, "physical_chiral_interacting_continuum_established": False}}, arrays


def source_bindings():
    paths = {ROOT/path for path in protected.source_bindings()}
    paths.add(Path(__file__))
    paths.update(ROOT/name for name in ("model/common_source/chiral_candidate/test_record_two_point.py",
        "model/common_source/hierarchy_diagnosis/audit.py", "model/common_source/replay.py", "tools/strict_json.py"))
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--check", type=Path)
    parser.add_argument("--steps", type=int, default=64)
    args = parser.parse_args()
    if args.check:
        retained = strict_json.load(args.check)
        for name, digest in retained["source_bindings"].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError("changed record-response source: "+name)
        steps = retained["source_steps"]
    else:
        if args.output.exists() or args.output.with_suffix(".npz").exists():
            parser.error("choose fresh JSON and NPZ paths")
        steps = args.steps
    initial_bindings = source_bindings()
    report, arrays = execute(steps=steps)
    report["source_bindings"] = source_bindings()
    if initial_bindings != report["source_bindings"]:
        raise ValueError("record-response sources changed during execution")
    if args.check:
        from ..hierarchy_diagnosis.audit import compare_replay
        artifact = args.check.with_name(retained["artifact"]["filename"])
        if hashlib.sha256(artifact.read_bytes()).hexdigest() != retained["artifact"]["sha256"]:
            raise ValueError("changed record-response quantum artifact")
        report["artifact"] = retained["artifact"]
        compare_replay(retained, report)
        with np.load(artifact, allow_pickle=False) as stored:
            from ..replay import compare_arrays
            compare_arrays(dict(stored), arrays)
        print("Recomputed complete neutral record tomography and every retained continuation output")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        artifact = args.output.with_suffix(".npz")
        np.savez_compressed(artifact, **arrays)
        report["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        print(args.output)


if __name__ == "__main__":
    main()
