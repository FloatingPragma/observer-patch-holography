"""Protected source propagation of the existing finite Weyl--Higgs candidate.

The seven retained momenta are one zero mode and three opposite pairs of
the same periodic spectral field.  At each momentum all M15, family and
Spin components remain.  Actual source Fock matrices on independent blocks
of at most eight modes retain the complete complex field output.  Their
tensor product, not a 2**630 allocation, defines the complete finite field.

Every fermion occupation is carried by record_code's qubit factor in its
full six-state code. Native correction is an exact tensor identity, with a
separate finite source error/service budget. Geometry, the cutoff, coherent
controls, reference clock and the unsuccessful trial Yukawa law remain
declared inputs. None is identified with observed Standard Model physics.
"""
from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

from . import action, gaussian_source as g, history, kinetic
from .. import record_code, record_field
from ..mass_response import history as instruments, model
from ...observer.source_scalar import compiler


MOMENTUM = F(1, 3)
SHORT_INTERVAL = F(1, 4)
ACCEPTED = F(15, 16)


def momentum_menu():
    """One periodic box of coordinate side6*pi; the cutoff retains7 modes."""
    p = float(MOMENTUM)
    return ((0., 0., 0.), (p, 0., 0.), (0., p, 0.), (0., 0., p))


def _steps(value):
    if type(value) is not int or not 1 <= value <= 4096:
        raise ValueError("a positive bounded source split count required")
    return value


def _duration(value):
    if type(value) is not F or value not in (SHORT_INTERVAL, 2*SHORT_INTERVAL):
        raise ValueError("one of the two declared exact probe intervals required")
    return float(value)


def _same_resources(left, right):
    return all((abs(left[key]-right[key]) <= 5e-12*max(1., abs(left[key]))
                if key == "integrated_processor_norm_upper" else left[key] == right[key])
               for key in left)


def decoded_source_word(word, sites):
    """Literal emitted source requests after the native code decoders."""
    return compiler.remap(word, tuple(9*site for site in sites))


def _lower_decoded_word(word, sites):
    literal = compiler._validate(word)
    inverse = {9*site: index for index, site in enumerate(sites)}
    if any(target not in inverse for event in literal for target in event["targets"]):
        raise ArithmeticError("protected source word touches a wrong occupation or syndrome wire")
    return [compiler.pulse(event["kind"], event["angle"],
                           [inverse[target] for target in event["targets"]]) for event in literal]


@lru_cache(maxsize=4096)
def _decoded_matrix(frozen_word, sites):
    word = [compiler.pulse(kind, angle, targets) for kind, angle, targets in frozen_word]
    return g._freeze(compiler.execute(_lower_decoded_word(word, sites), len(sites)))


def _executed_decoded_gate(gate, sites):
    # The actual emitted literal, not an expected target set, is the cache key.
    literal = decoded_source_word(gate.word, sites)
    frozen = tuple((event["kind"], event["angle"], tuple(event["targets"])) for event in literal)
    actual = _decoded_matrix(frozen, sites)
    expected = gate.fock_unitary[:2, :2] if len(sites) == 1 else gate.fock_unitary
    if np.max(abs(actual-expected)) > 2e-11:
        raise ArithmeticError("literal protected source word changes the declared field factor order")
    return actual, literal


@dataclass(frozen=True)
class SmallEvolution:
    unitary: np.ndarray
    target: np.ndarray
    majorana: np.ndarray
    resources: dict
    source_error: float
    local_pulse_menu: tuple
    local_calls: tuple


def _fock_hamiltonian(h, delta):
    """Independent bitwise CAR construction, including every Fock sector."""
    modes = len(h)
    result = np.zeros((1 << modes, 1 << modes), complex)

    def apply(mask, site, create):
        occupied = bool(mask & (1 << site))
        if occupied == create:
            return None
        return mask ^ (1 << site), (-1)**((mask & ((1 << site)-1)).bit_count())

    for initial in range(1 << modes):
        for i in range(modes):
            for j in range(modes):
                if h[i, j] != 0:
                    first = apply(initial, j, False)
                    second = None if first is None else apply(first[0], i, True)
                    if second is not None:
                        result[second[0], initial] += h[i, j]*first[1]*second[1]
            for j in range(i+1, modes):
                if delta[i, j] == 0:
                    continue
                first = apply(initial, j, True)
                second = None if first is None else apply(first[0], i, True)
                if second is not None:
                    result[second[0], initial] += delta[i, j]*first[1]*second[1]
                first = apply(initial, i, False)
                second = None if first is None else apply(first[0], j, False)
                if second is not None:
                    result[second[0], initial] += delta[i, j].conjugate()*first[1]*second[1]
    if np.max(abs(result-result.conj().T)) > 2e-13:
        raise ArithmeticError("independent Fock Hamiltonian lost Hermiticity")
    return result


@lru_cache(maxsize=192)
def _small(hbytes, dbytes, modes, duration, steps):
    h = np.frombuffer(hbytes, complex).reshape(modes, modes)
    delta = np.frombuffer(dbytes, complex).reshape(modes, modes)
    if not 1 <= modes <= 8:
        raise ValueError("the literal comparison retains at most eight Fock modes")
    _, terms, _ = g._quadratic(h, delta)
    state = np.eye(1 << modes, dtype=complex)
    resources = {key: 0 for key in g._ADDITIVE_RESOURCES}
    menu, calls = set(), {}
    swap = g.local_gate("fswap")

    def apply(gate, sites):
        nonlocal state
        matrix, literal = _executed_decoded_gate(gate, tuple(sites))
        state = model.apply_subsystem(state, matrix, sites, modes)
        actual_resources = compiler.resources(literal)
        for key in resources:
            resources[key] += actual_resources[key]
        menu.update((event["kind"], event["angle"]) for event in literal)
        key = gate.kind, gate.coefficient.real, gate.coefficient.imag, gate.duration, tuple(sites)
        calls[key] = calls.get(key, 0)+1

    half = duration/(2*steps)
    prepared = [(kind, sites, g.local_gate(kind, coefficient, half))
                for kind, sites, coefficient in terms]
    for _, sites, gate in (*prepared, *reversed(prepared)):
        if len(sites) == 1:
            apply(gate, sites)
            continue
        a, b = sites
        route = [(j-1, j) for j in range(b, a+1, -1)]
        for pair in route:
            apply(swap, pair)
        apply(gate, (a, a+1))
        for pair in reversed(route):
            apply(swap, pair)
    # The literal source step repeats unchanged. Contracting its computed
    # matrix by integer power saves host work, never physical source requests.
    # No target exponential replaces a source gate or its intermediate word.
    state = np.linalg.matrix_power(state, steps)
    resources = {key: value*steps for key, value in resources.items()}
    calls = {key: count*steps for key, count in calls.items()}
    # Source and independent comparison are separate paths. No target
    # exponential, polar factor or normalized state supplies an executed gate.
    target = expm(-1j*duration*_fock_hamiltonian(h, delta))
    gaussian = g.evolve(h, delta, duration, steps=steps)
    if not _same_resources(resources, gaussian.resources):
        raise ArithmeticError("literal and Gaussian source schedules disagree")
    error = float(np.linalg.norm(state-target, 2))
    if not math.isfinite(error):
        raise ArithmeticError("nonfinite complete small-block source comparison")
    return SmallEvolution(g._freeze(state), g._freeze(target), gaussian.orthogonal,
                          dict(resources), error, tuple(sorted(menu)),
                          tuple((*key, count) for key, count in calls.items()))


def _small_evolution(h, delta, duration, steps):
    return _small(np.asarray(h, complex).tobytes(), np.asarray(delta, complex).tobytes(),
                  len(h), float(duration), _steps(steps))


def field_layout(momentum):
    """All CAR modes; charged connected components precede neutral ones.

    Nonzero p retains both p and -p. Zero is a single self-conjugate
    momentum and has90 CAR modes, never two copies of the zero mode.
    """
    p = kinetic._momentum(momentum)
    kernel = np.asarray(action.mass_kernel().evalf(17), complex)
    moving = any(p)
    if moving:
        h, delta = kinetic.momentum_pair_source(p)
    else:
        h, delta = kinetic.weyl_kinetic(p), kinetic.hamiltonian_pairing()
    chosen = next(row["index"] for row in action.mode_labels()
                  if row["spin"] == 0 and row["internal_type"] == "L"
                  and row["family"] == 0 and F(row["electric_Q"]) == -1)
    partner = np.flatnonzero(kernel[chosen] != 0)
    if len(partner) != 1:
        raise ArithmeticError("the selected Higgs line is not a charged pair")
    pairs = [(chosen, int(partner[0]))]
    used = set(pairs[0])
    for i in range(45):
        if i in used:
            continue
        partners = np.flatnonzero(kernel[i] != 0)
        if len(partners) == 1:
            j = int(partners[0])
            if j in used:
                raise ArithmeticError("overlapping internal Higgs components")
            pairs.append((i, j)); used.update((i, j))
        elif len(partners) != 0:
            raise ArithmeticError("nonfactorizing trial Higgs kernel")
    internal_groups = [*pairs, *((i,) for i in range(45) if i not in used)]
    groups = [tuple(90*q+45*spin+i for q in range(2 if moving else 1)
                    for spin in (0, 1) for i in group) for group in internal_groups]
    order = tuple(i for group in groups for i in group)
    if len(order) != (180 if moving else 90) or len(set(order)) != len(order):
        raise ArithmeticError("a momentum or internal field component was discarded")
    h, delta = h[np.ix_(order, order)], delta[np.ix_(order, order)]
    sizes = tuple(map(len, groups))
    cut = 0
    for size in sizes:
        other = list(range(cut))+list(range(cut+size, len(h)))
        if np.any(h[cut:cut+size, other]) or np.any(delta[cut:cut+size, other]):
            raise ArithmeticError("the retained Fock factorization mixes components")
        cut += size
    return h, delta, {"momentum": list(p), "opposite_momentum_retained": moving,
        "CAR_modes": len(h), "original_fermionic_wire_order": list(order),
        "factor_sizes": list(sizes), "selected_probe_Fock_modes": sizes[0],
        "all_M15_family_and_Spin_components_retained": True,
        "zero_momentum_duplicated": False, "modes_discarded": 0}


@dataclass(frozen=True)
class FieldBlock:
    layout: dict
    factors: tuple[SmallEvolution, ...]
    resources: dict
    tensor_operator_error: float
    majorana_factorization_residual: float


@lru_cache(maxsize=16)
def _field_block(momentum, duration, steps):
    h, delta, layout = field_layout(momentum)
    factors, cut = [], 0
    composed = np.zeros((2*len(h), 2*len(h)))
    resources = {key: 0 for key in g._ADDITIVE_RESOURCES}
    logerror = 0.
    for size in layout["factor_sizes"]:
        selected = slice(cut, cut+size)
        factor = _small_evolution(h[selected, selected], delta[selected, selected], duration, steps)
        factors.append(factor)
        composed[2*cut:2*(cut+size), 2*cut:2*(cut+size)] = factor.majorana
        for key in resources:
            resources[key] += factor.resources[key]
        logerror += math.log1p(factor.source_error)
        cut += size
    full = g.evolve(h, delta, duration, steps=steps)
    residual = float(np.linalg.norm(full.orthogonal-composed, 'fro'))
    if (not math.isfinite(residual) or residual > 2e-9
            or not _same_resources(resources, full.resources)):
        raise ArithmeticError(f"full momentum source and literal Fock factors disagree: {residual}")
    return FieldBlock(layout, tuple(factors), resources, math.expm1(logerror), residual)


def field_block(momentum, duration=SHORT_INTERVAL, *, steps=16):
    """Actual complete finite field factorization, after original-input checks."""
    p = kinetic._momentum(momentum)
    time = _duration(duration)
    result = _field_block(p, time, _steps(steps))
    # Validate the current literal producer even after warming numerical
    # caches, and detach mutable metadata from those caches.
    for factor in result.factors:
        for kind, real, imag, local_time, sites, _ in factor.local_calls:
            _executed_decoded_gate(g.local_gate(kind, complex(real, imag), local_time), sites)
    layout = deepcopy(result.layout)
    layout["source_split_count"] = steps
    factors = tuple(SmallEvolution(f.unitary, f.target, f.majorana, dict(f.resources),
                                  f.source_error, f.local_pulse_menu, f.local_calls)
                    for f in result.factors)
    return FieldBlock(layout, factors, dict(result.resources), result.tensor_operator_error,
                      result.majorana_factorization_residual)


def _number_vectors(vector, site, *, continuous=False):
    if not continuous:
        _, kraus, _ = g._occupation_pointer()
        modes = len(vector).bit_length()-1
        return tuple(model.apply_subsystem(vector[:, None], k, [site], modes)[:, 0]
                     for k in kraus)
    indices = np.arange(len(vector))
    empty, occupied = vector.copy(), vector.copy()
    empty[(indices >> site) & 1 == 1] = 0
    occupied[(indices >> site) & 1 == 0] = 0
    return empty, occupied


def _clock_outputs(clock, elapsed, duration, *, continuous=False):
    state = history._clock_advance(clock, elapsed*duration, qubits=2,
                                   backend="continuous" if continuous else "source")
    branches = {(): state}
    for site, axis in ((0, "X"), (1, "Y")):
        outputs = {}
        for record, vector in branches.items():
            if continuous:
                sigma = np.array([[0, 1], [1, 0]]) if axis == "X" else np.array([[0, -1j], [1j, 0]])
                changed = model.apply_subsystem(vector, sigma, [site], 2)
                states = ((vector+changed)/2, (vector-changed)/2)
            else:
                states = instruments.pauli_read(vector, site, axis=axis).outputs
            outputs.update({(*record, value): output for value, output in enumerate(states)})
        branches = outputs
    return branches


def protocol(blocks, selected, duration, *, continuous=False):
    """One preparation, two native/field cycles and every stopped output.

    Complete field outputs are pure complex tensor factors. Other momentum
    and internal factors are retained at each elapsed-cycle value. Only the
    selected factor is disturbed by reads/feedback. Clock purifiers remain
    explicit columns; no probability or norm repair is used.
    """
    time = _duration(duration)
    if type(selected) is not int or not 0 <= selected < len(blocks):
        raise ValueError("a retained momentum block must be selected")
    probe = blocks[selected].factors[0]
    unitary = probe.target if continuous else probe.unitary
    clock, calibration, clock_cost = history._clock_preparation()
    clock_branches = {n: _clock_outputs(clock, n, time, continuous=continuous) for n in (0, 1, 2)}
    spectator_states = {}
    spectator_mass = {n: 1. for n in (0, 1, 2)}
    for block_index, block in enumerate(blocks):
        for factor_index, factor in enumerate(block.factors):
            if (block_index, factor_index) == (selected, 0):
                continue
            matrix = factor.target if continuous else factor.unitary
            vector = np.eye(len(matrix), dtype=complex)[:, 0]
            for elapsed in (0, 1, 2):
                if elapsed:
                    vector = matrix@vector
                spectator_states[block_index, factor_index, elapsed] = g._freeze(vector)
                spectator_mass[elapsed] *= float(np.vdot(vector, vector).real)
    blank = np.eye(len(unitary), dtype=complex)[:, 0]
    after_first = _number_vectors(unitary@blank, 0, continuous=continuous)
    phase = g.local_gate("onsite", -math.pi/2, 1.)
    phase_matrix = np.diag([1., 1j]) if continuous else phase.fock_unitary[:2, :2]
    continued = []
    for first, vector in enumerate(after_first):
        if first:
            vector = model.apply_subsystem(vector[:, None], phase_matrix, [1],
                                          probe.unitary.shape[0].bit_length()-1)[:, 0]
        continued.append(vector)
    terminal = []

    def add(flag, elapsed, first, second, selector_raw, weight, vector):
        for record, clock_state in clock_branches[elapsed].items():
            probability = (float(weight)*spectator_mass[elapsed]
                           *float(np.vdot(vector, vector).real)
                           *float(np.vdot(clock_state, clock_state).real))
            terminal.append({"flag": flag, "elapsed_cycles": elapsed,
                "first": first, "second": second, "selector_raw": selector_raw,
                "clock_X": record[0], "clock_Y": record[1],
                "selector_weight_exact": str(weight), "probability": probability,
                "probe_vector": g._freeze(vector), "clock_reference_columns": g._freeze(clock_state)})

    for raw in (30, 31):
        add("first_selector_timeout", 0, None, None, raw, F(1, 32), blank)
    for first, vector in enumerate(continued):
        for raw in (30, 31):
            add("second_selector_timeout", 1, first, None, raw, ACCEPTED/32, vector)
        for second, output in enumerate(_number_vectors(unitary@vector, 0, continuous=continuous)):
            add("completed", 2, first, second, None, ACCEPTED**2, output)
    mass = math.fsum(row["probability"] for row in terminal)
    if not math.isfinite(mass) or abs(mass-1) > 2e-7:
        raise ArithmeticError("complete protected momentum history lost probability")
    return {"terminal": terminal, "spectator_states": spectator_states,
        "clock_calibration": calibration, "clock_preparation_resources": clock_cost,
        "total_probability": mass, "work_preparations": 1, "work_resets": 0,
        "native_and_selector_records_are_separate_exact_factors": True,
        "zero_and_abort_outcomes_retained": True,
        "field_pointer_contract": "two source-resolved number projectors; a first1 applies a neutral source phase on mode1",
        "proper_code_export_failures": "retained first-failure ordinal blocks, bounded separately from the complete successful field factors"}


def infer_frequency(short, long, calibration, phase_contract):
    """Cancel unknown Rabi visibility using actual first-read records.

    P(t)=A sin²(E t), P(2t)=A sin²(2E t). The clock supplies Omega*t.
    Neither Hamiltonian coefficients, momenta, eigensolves nor simulator
    durations enter this decoder. Later read backaction remains in each
    complete history and cannot alter the earlier recorded marginal.
    """
    if (type(phase_contract) is not dict
            or phase_contract.get("same_law_and_preparation_at_doubled_interval") is not True
            or phase_contract.get("short_energy_phase_strictly_between_zero_and_pi_over_two") is not True):
        raise ValueError("a declared same-law doubled-interval alias window is required")
    values, angles = [], []
    for packet in (short, long):
        terminal = packet["terminal"]
        if any(not math.isfinite(row["probability"]) or row["probability"] < 0 for row in terminal):
            raise ValueError("nonnegative finite complete recorded probabilities required")
        if abs(math.fsum(row["probability"] for row in terminal)-1) > 2e-7:
            raise ValueError("complete unnormalized selector/detector mass is required")
        values.append(math.fsum(row["probability"] for row in terminal if row["first"] == 1)/float(ACCEPTED))
        means = [math.fsum((1 if row[key] == 0 else -1)*row["probability"]
                          for row in terminal if row["flag"] == "completed")
                 /float(ACCEPTED**2)/calibration["contrast"] for key in ("clock_X", "clock_Y")]
        if abs(math.hypot(*means)-1) > 2e-7:
            raise ValueError("source clock leaves its declared calibrated circle")
        angle = math.atan2(means[1], means[0])
        if not 0 < angle < math.pi:
            raise ValueError("clock phase is outside the admitted window")
        angles.append(angle)
    if abs(angles[1]-2*angles[0]) > 2e-7 or not 0 < values[0] < 1:
        raise ValueError("unresolved visibility or incompatible doubled-interval clock")
    ratio = values[1]/(4*values[0])
    if not 0 < ratio < 1:
        raise ValueError("no resolved principal-branch Rabi frequency")
    phase = math.acos(math.sqrt(ratio))
    visibility = values[0]/math.sin(phase)**2
    if not 0 < visibility <= 1+2e-5:
        raise ValueError("recorded Rabi visibility is outside the declared two-level law")
    inverse_derivative = 1/(angles[0]*math.sqrt(ratio*(1-ratio)))
    return {"energy_over_clock": 2*phase/angles[0], "short_energy_phase": phase,
        "measured_visibility": visibility, "short_first_occupation_probability": values[0],
        "long_first_occupation_probability": values[1], "short_final_clock_phase": angles[0],
        "local_inverse_sensitivity_diagnostics": {
            "absolute_derivative_wrt_short_probability": inverse_derivative*ratio/values[0],
            "absolute_derivative_wrt_long_probability": inverse_derivative/(4*values[0]),
            "absolute_derivative_wrt_short_final_clock_phase": 2*phase/angles[0]**2},
        "decoder_uses_Hamiltonian_momentum_or_simulator_duration": False,
        "recovered_frequency_interval_or_physical_pole_certified": False,
        "probability_clock_and_floating_error_propagation_outward_certified": False,
        "finite_shot_or_physical_service_error_certified": False}


def protection_budget(blocks, duration):
    """Reference-uniform complete-source budget, scaled to all630 carriers.

    Small Fock comparisons include splitting and successful source-word
    errors on every input column. Native code error covers every seam and
    recorded recovery. All first export failures are orthogonal ordinal
    blocks; their computed amplitudes are bounded before channel conversion.
    These floating calculations are not outward or hardware certificates.
    """
    time = _duration(duration)
    if tuple(tuple(block.layout["momentum"]) for block in blocks) != momentum_menu():
        raise ValueError("the fixed momentum menu and its complete source factors are required")
    owned = tuple(field_block(momentum, duration, steps=block.layout["source_split_count"])
                  for momentum, block in zip(momentum_menu(), blocks, strict=True))
    for supplied, current in zip(blocks, owned, strict=True):
        if (supplied.layout != current.layout or supplied.resources != current.resources
                or supplied.tensor_operator_error != current.tensor_operator_error
                or len(supplied.factors) != len(current.factors)):
            raise ValueError("supplied protected field budget differs from owned source execution")
        for a, b in zip(supplied.factors, current.factors, strict=True):
            if (a.resources != b.resources or a.source_error != b.source_error
                    or not np.array_equal(a.unitary, b.unitary)
                    or not np.array_equal(a.target, b.target)):
                raise ValueError("supplied protected factor differs from owned source execution")
    blocks = owned
    modes = sum(block.layout["CAR_modes"] for block in blocks)
    if modes != 630:
        raise ValueError("the complete fixed seven-momentum field is required")
    native = record_code.source_composition_report()
    code_exact = record_code.exact_certificate()
    encoder, decoder = record_code.encoding_resources(), record_code.decoding_resources()
    clock, calibration, clock_cost = history._clock_preparation()
    pop = np.array([float(F(value)) for value in calibration["populations"]])
    hadamard = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    clock_target = -np.kron(hadamard, hadamard)@np.diag(np.sqrt(np.kron(pop, pop)))
    preparation_error = float(np.linalg.norm(clock-clock_target, 'fro'))
    field_log = math.fsum(math.log1p(block.tensor_operator_error) for block in blocks)
    pointer_word, pointer_kraus, _ = g._occupation_pointer()
    pointer_error = float(np.linalg.norm(np.concatenate(pointer_kraus)-np.vstack((np.diag([1., 0.]), np.diag([0., 1.]))), 2))
    reads = [instruments.pauli_read(np.eye(2, dtype=complex), 0, axis=axis) for axis in ("X", "Y")]
    read_errors = []
    for axis, read in zip(("X", "Y"), reads, strict=True):
        pauli = np.array([[0., 1.], [1., 0.]]) if axis == "X" else np.array([[0., -1j], [1j, 0.]])
        target = np.vstack(((np.eye(2)+pauli)/2, (np.eye(2)-pauli)/2))
        read_errors.append(float(np.linalg.norm(np.concatenate(read.outputs)-target, 2)))
    feedback = g.local_gate("onsite", -math.pi/2, 1.)
    feedback_error = float(np.linalg.norm(feedback.fock_unitary-np.diag([1., 1j, 1., 1j]), 2))
    clock_angle = -model.parameters().clock*time/2
    clock_error = record_code._primitive_error("Z", clock_angle)
    auxiliary_log = (math.log1p(preparation_error)+2*math.log1p(pointer_error)
                     +sum(map(math.log1p, read_errors))+math.log1p(feedback_error)
                     +4*math.log1p(clock_error))
    pulse_menu = {item for block in blocks for factor in block.factors for item in factor.local_pulse_menu}
    for word in (pointer_word, feedback.word, *(read.word for read in reads)):
        pulse_menu.update((event["kind"], event["angle"]) for event in word)
    amplitudes = np.sqrt(pop)
    preparation_word = []
    for bit, purifier in ((0, 2), (1, 3)):
        preparation_word += compiler.remap(instruments.preparation_word(amplitudes), [bit])
        preparation_word += compiler.cnot(bit, purifier)
        preparation_word.append(compiler.pulse("H", math.pi/2, [bit]))
    pulse_menu.update((event["kind"], event["angle"]) for event in preparation_word)
    pulse_menu.add(("Z", clock_angle))
    export = dict(record_field._proper_code_export_diagnostic())
    export["native_base_pulse_menu_count"] = export.pop("actualpulse_menu_count")
    export["native_base_pulse_family_counts"] = export.pop("pulse_family_counts")
    export["additional_field_clock_detector_pulse_menu_count"] = len(pulse_menu)
    for kind, angle in pulse_menu:
        full = compiler.native_pulse(kind, angle)
        size = 4 if kind == "J" else 2
        if full.shape != (6, 6) or not np.isfinite(full).all():
            raise ArithmeticError("nonfinite full source pulse in protected momentum program")
        export["maximum_complement_operator_norm"] = max(
            export["maximum_complement_operator_norm"], float(np.linalg.norm(full[size:, :size], 2)))
        export["maximum_success_norm_excess"] = max(
            export["maximum_success_norm_excess"], float(np.linalg.norm(full[:size, :size], 2))-1)
    if any(not math.isfinite(export[key]) or not 0 <= export[key] <= 2e-11
           for key in ("maximum_complement_operator_norm", "maximum_success_norm_excess")):
        raise ArithmeticError("full momentum pulse export exceeds its fixed numerical guard")
    resource_keys = ("native_pulses", "proper_code_load_events", "proper_code_export_events", "proper_code_service_events")
    field_resources = {key: sum(block.resources[key] for block in blocks) for key in resource_keys}
    auxiliary_pulses = (clock_cost["native_pulses"]+2*len(pointer_word)
                        +sum(read.resources["native_pulses"] for read in reads)
                        +feedback.resources["native_pulses"]+4)
    laws = {}
    for law in ("C", "E"):
        local = native["laws"][law]
        native_log = (modes*native["encoder_log_error"]
            +2*modes*(local["maximum_corrected_attempt_log_error"]
                      +native["decoder_log_error"]+native["encoder_log_error"]))
        resources = {key: modes*encoder[key]+2*modes*(
            local["maximum_corrected_attempt_resources"][key]+decoder[key]+encoder[key])
            +2*field_resources[key] for key in resource_keys}
        resources["native_pulses"] += auxiliary_pulses
        resources["proper_code_load_events"] += auxiliary_pulses
        resources["proper_code_export_events"] += auxiliary_pulses
        resources["proper_code_service_events"] += 2*auxiliary_pulses
        pulses = resources["native_pulses"]
        successful = math.expm1(native_log+2*field_log+auxiliary_log)
        failure = record_field._first_export_failure_bound(pulses, export)
        epsilon = successful+failure
        # A fixed numerical admission criterion for this larger program;
        # the old six-mode experiment's 1e-6 gate is not transferred here.
        if not math.isfinite(epsilon) or not 0 <= epsilon < .01:
            raise ArithmeticError("complete630-mode source comparison exceeds its fixed1e-2 numerical admission")
        raw_bits = 2 if law == "C" else 3
        archive = 10+2*modes*(1+raw_bits)+4
        export_bits = pulses.bit_length()
        laws[law] = {"complete_successful_operator_error_diagnostic": successful,
            "first_export_failure_amplitude_diagnostic": failure,
            "complete_flagged_operator_error_diagnostic": epsilon,
            "complete_cq_half_diamond_error_diagnostic": epsilon*(1+epsilon/2),
            "fixed_complete_operator_admission": .01, "source_resources": resources,
            "native_attempts_per_two_cycle_history": 2*modes,
            "maximum_raw_archive_bits": archive,
            "first_export_failure_ordinal_bits": export_bits,
            "maximum_archive_and_export_ordinal_bits": archive+export_bits,
            "physical_service_invocations": {
                "native_blank_code_register": modes, "blank_qubits_per_register": 9,
                "native_coin_and_raw_blank_pair": 2*modes,
                "native_coin_or_raw_record_commit": 4*modes,
                "native_seam_route_and_record_recovery_control": 2*modes,
                "native_logical_storage_for_complete_window": modes,
                "reference_clock_preparation": 1, "clock_purification_reference_storage": 1,
                "fair_selector_bit": 10, "selector_five_bit_commit": 2,
                "source_advance_or_halt_control": 2, "number_read_and_commit": 2,
                "clock_quadrature_read_and_commit": 2, "conditional_neutral_phase_control": 1,
                "first_export_failure_ordinal_commit": 1,
                "classical_archive_bit_storage_for_complete_window": archive+export_bits,
                "source_analogue_pulse": pulses,
                "proper_code_coherent_load": resources["proper_code_load_events"],
                "proper_code_complete_export": resources["proper_code_export_events"]},
            "physical_service_error_values_supplied": False}
    native_factors = {}
    for law in ("C", "E"):
        row = record_field.record_summary(law)
        row["native_carrier_factors_per_accepted_seam"] = modes
        row["complete_native_record_addresses_per_accepted_cycle"] = 30*(8 if law == "C" else 16)**modes
        row["raw_abort_policy"] = "complete all630 local record attempts; retain their work/records; any raw abort halts before field decode/advance; ideal code weight is exactly zero"
        row["corrected_branch_action"] = "sqrt(joint native branch weight)*I on all630 logical factors and arbitrary references"
        row["representation"] = "complete exact product of630 local factors; no512**630 state array"
        native_factors[law] = row
    return {"native_carriers": modes, "quantum_work_code_bits": 9*modes,
        "occupation_factor": "code6=C3_spectator tensor C2_fermion, label=2*s+n",
        "code_spectator_multiplicity": str(3**modes),
        "all_native_seams_and_raw_records_certified": True,
        "native_exact_capacity_certificate": code_exact,
        "native_record_factors": native_factors,
        "native_tensor_branch_action": "sqrt(product native weights)*I on the full630-mode CAR field and arbitrary references",
        "decoded_source_word_wiring": "factor-local source words are relocated into the checked disjoint global slot map; CAR index j is decoded native wire9*j; other local code wires are identity",
        "zero_native_raw_abort_and_padding_outputs_retained": True,
        "computed_export_failures_are_orthogonal_first_failure_ordinals": True,
        "native_record_tree_or_full_native_tensor_materialized": False,
        "source_field_tensor_operator_error_per_cycle": math.expm1(field_log),
        "source_clock_preparation_vector_error": preparation_error,
        "source_detector_complete_operator_errors": [pointer_error, *read_errors],
        "proper_code_export_diagnostic": export, "laws": laws,
        "physical_duration_energy_error_or_Standard_Model_spectrum_certified": False,
        "floating_diagnostics_are_outward_certificates": False}


def _conditional_probe_comparison(actual, target):
    if len(actual["terminal"]) != 40 or len(target["terminal"]) != 40:
        raise ArithmeticError("a stopped selector/detector branch was omitted")
    maximum = probability = 0.
    for a, b in zip(actual["terminal"], target["terminal"], strict=True):
        for key in ("flag", "elapsed_cycles", "first", "second", "selector_raw", "clock_X", "clock_Y"):
            if a[key] != b[key]:
                raise ArithmeticError("source and target terminal records differ")
        # Full complex probe and clock/purifier columns, including zeros.
        left = np.kron(a["clock_reference_columns"], a["probe_vector"][:, None])
        right = np.kron(b["clock_reference_columns"], b["probe_vector"][:, None])
        maximum = max(maximum, math.sqrt(float(F(a["selector_weight_exact"])))
                      *float(np.linalg.norm(left-right, 'fro')))
        probability = max(probability, abs(a["probability"]-b["probability"]))
    if not math.isfinite(maximum) or maximum > .01:
        raise ArithmeticError("complete complex probe/reference output exceeds fixed numerical admission")
    return {"maximum_complex_probe_and_reference_output_error": maximum,
            "maximum_complete_record_probability_difference": probability,
            "probabilities_alone_used_as_quantum_output_oracle": False}


def _current_checks(blocks):
    """The same decoded M15 current/charge acts on every retained momentum."""
    currents = [np.kron(np.eye(2), np.asarray(t.evalf(17), complex))
                for t in action.generators()]
    charge = np.diag(np.asarray(action.spin_electric_charge().diagonal(), float).ravel())
    maximum_kinetic = maximum_charge = maximum_source = 0.
    for block in blocks:
        p = block.layout["momentum"]
        h, delta, layout = field_layout(p)
        moving = any(p)
        original_h = kinetic.weyl_kinetic(p)
        maximum_kinetic = max(maximum_kinetic, *(float(np.max(abs(t@original_h-original_h@t))) for t in currents))
        original_q = np.kron(np.eye(2 if moving else 1), charge)
        order = layout["original_fermionic_wire_order"]
        q = np.diag(original_q)[order]
        maximum_charge = max(maximum_charge, float(np.max(abs(q[:, None]*h-h*q[None, :]))),
                             float(np.max(abs(q[:, None]*delta+delta*q[None, :]))))
        cut = 0
        for size, factor in zip(layout["factor_sizes"], block.factors, strict=True):
            number_charge = np.array([sum(q[cut+j] for j in range(size) if mask >> j & 1)
                                      for mask in range(1 << size)])
            maximum_source = max(maximum_source, float(np.max(abs(
                number_charge[:, None]*factor.unitary-factor.unitary*number_charge[None, :]))))
            cut += size
    if max(maximum_kinetic, maximum_charge, maximum_source) > 2e-10:
        raise ArithmeticError("the same M15 kinetic current or completed source electric charge changed")
    return {"same_retained_internal_generators_per_momentum": 12,
        "all_kinetic_current_commutator_residual": maximum_kinetic,
        "neutral_Higgs_model_charge_covariance_residual": maximum_charge,
        "all_Fock_column_source_model_electric_charge_residual": maximum_source,
        "broken_electroweak_generators_are_asserted_unbroken": False,
        "physical_electromagnetic_current_selected": False}


def execute(*, steps=16, with_protection_budget=True):
    """Execute the fixed field at both intervals; return compact exact factors.

    A protocol selects which existing field factor is read, while every
    factor evolves. It does not select a different fitted action. Host caches
    reuse contractions only; every represented native request is charged.
    """
    _steps(steps)
    if type(with_protection_budget) is not bool:
        raise ValueError("an explicit source-budget switch required")
    menu = momentum_menu()
    intervals = (SHORT_INTERVAL, 2*SHORT_INTERVAL)
    blocks = {time: tuple(field_block(p, time, steps=steps*(1 if time == SHORT_INTERVAL else 2))
                         for p in menu) for time in intervals}
    doubled_error = max(float(np.linalg.norm(long.unitary-short.unitary@short.unitary, 2))
                        for a, b in zip(blocks[SHORT_INTERVAL], blocks[2*SHORT_INTERVAL], strict=True)
                        for short, long in zip(a.factors, b.factors, strict=True))
    if doubled_error > 2e-8:
        raise ArithmeticError("doubled-duration protocol changed the elementary source law")
    phase_contract = {"same_law_and_preparation_at_doubled_interval": True,
        "short_energy_phase_strictly_between_zero_and_pi_over_two": True,
        "proof": "selected lepton m<1/2, c<1, |p|<=1/3 imply E²<13/36<1; t=1/4; pi>3",
        "other_or_unbounded_frequency_laws_admitted": False}
    protocols, decoded, comparisons = {}, {}, {}
    arrays = {}
    for index, momentum in enumerate(menu):
        actual = [protocol(blocks[time], index, time) for time in intervals]
        target = [protocol(blocks[time], index, time, continuous=True) for time in intervals]
        decoded[str(index)] = infer_frequency(*actual, actual[0]["clock_calibration"], phase_contract)
        exact = infer_frequency(*target, target[0]["clock_calibration"], phase_contract)
        comparisons[str(index)] = {"independent_continuous_record_decoder": exact,
            "source_frequency_over_clock_difference": abs(decoded[str(index)]["energy_over_clock"]-exact["energy_over_clock"]),
            "interval_outputs": [_conditional_probe_comparison(a, b) for a, b in zip(actual, target, strict=True)]}
        for interval_index, packet in enumerate(actual):
            name = f"p{index}_t{interval_index}"
            protocols[name] = {key: value for key, value in packet.items() if key not in ("terminal", "spectator_states")}
            protocols[name]["terminal"] = [{key: value for key, value in row.items()
                if key not in ("probe_vector", "clock_reference_columns")} for row in packet["terminal"]]
            arrays[name+"_probe_vectors"] = np.stack([row["probe_vector"] for row in packet["terminal"]])
            arrays[name+"_clock_reference_columns"] = np.stack([row["clock_reference_columns"] for row in packet["terminal"]])
    # Store every shared field factor once, not once per record/protocol.
    factor_records = {}
    global_bindings_checked = 0
    for interval_index, time in enumerate(intervals):
        block_offset = 0
        for block_index, block in enumerate(blocks[time]):
            factor_offset = block_offset
            for factor_index, factor in enumerate(block.factors):
                name = f"t{interval_index}_p{block_index}_f{factor_index}"
                blank = np.eye(len(factor.unitary), dtype=complex)[:, 0]
                arrays[name] = np.stack((blank, factor.unitary@blank, factor.unitary@factor.unitary@blank))
                size = len(factor.unitary).bit_length()-1
                for kind, real, imag, local_time, sites, count in factor.local_calls:
                    gate = g.local_gate(kind, complex(real, imag), local_time)
                    global_sites = tuple(factor_offset+site for site in sites)
                    lowered = _lower_decoded_word(decoded_source_word(gate.word, global_sites), global_sites)
                    if lowered != list(gate.word):
                        raise ArithmeticError("global protected field word changes its declared local source action")
                    global_bindings_checked += 1
                factor_records[name] = {"Fock_modes": size,
                    "CAR_modes": block.layout["CAR_modes"], "elapsed_cycles": [0, 1, 2],
                    "global_native_CAR_indices": list(range(factor_offset, factor_offset+size)),
                    "literal_decoded_occupation_wires": [9*j for j in range(factor_offset, factor_offset+size)],
                    "all_Fock_input_columns_executed": len(factor.unitary),
                    "complete_source_operator_error": factor.source_error,
                    "source_resources_per_cycle": factor.resources,
                    "complex_factor_sha256": hashlib.sha256(arrays[name].tobytes()).hexdigest()}
                factor_offset += size
            if factor_offset != block_offset+block.layout["CAR_modes"]:
                raise ArithmeticError("global native carrier allocation omitted a factor")
            block_offset = factor_offset
        if block_offset != 630:
            raise ArithmeticError("global native carrier allocation is incomplete")
    rest = decoded["0"]["energy_over_clock"]
    shell = {str(i): {"measured_energy_squared_minus_rest_squared_over_clock_squared":
                       decoded[str(i)]["energy_over_clock"]**2-rest**2,
                     "declared_common_symbol_value": float(MOMENTUM)**2/4}
             for i in range(1, 4)}
    return {"schema": "oph.chiral_candidate.protected_momentum.v1",
        "field": {"retained_momentum_representatives": [list(p) for p in menu],
            "opposite_pairs_in_addition_to_single_zero": 3, "total_retained_momenta": 7,
            "CAR_modes_per_momentum": 90, "total_CAR_modes": 630,
            "periodic_coordinate_side": "6*pi", "nonzero_coordinate_momentum": "1/3",
            "clock_parameter": "the same programmed native-reference modular parameter",
            "coordinate_speed": "Omega/2", "same_supplied_trial_Higgs_and_Yukawa_law_in_every_block": True,
            "finite_spectral_cutoff_is_a_local_gauge_regulator": False,
            "spatial_coordinate_volume_or_clock_unit_physically_calibrated": False},
        "layouts": [block.layout for block in blocks[SHORT_INTERVAL]],
        "same_field_current_attachment": _current_checks(blocks[SHORT_INTERVAL]),
        "source_interval_relation": {"short_steps": steps, "long_steps": 2*steps,
            "same_elementary_source_word_repeated_twice": True,
            "all_Fock_column_doubled_word_operator_error": doubled_error},
        "phase_contract": phase_contract, "protocols": protocols, "decoded": decoded,
        "dispersion_readout": shell, "independent_comparisons": comparisons,
        "complete_complex_tensor_factors": factor_records,
        "literal_global_native_factor_word_bindings_checked": global_bindings_checked,
        "factor_output_rule": "each successful-code terminal is sqrt(selector weight) times its retained probe vector, clock/purifier matrix, and every other field factor at its elapsed cycle; multiply exact native record factors. Ideal zero raw/padding blocks remain zero; nonzero computed export-failure quantum outputs are bounded flagged blocks, not explicitly retained vectors",
        "Gaussian_and_full_Fock_factorization_residuals": {
            str(time): [block.majorana_factorization_residual for block in blocks[time]] for time in intervals},
        "protection_and_service_budgets": {str(time): protection_budget(blocks[time], time)
                                           for time in intervals} if with_protection_budget else {},
        "scope": {"same_declared_field_has_source_executed_three_axis_nonzero_momentum_response": True,
            "local_record_code_correction_preserves_complete_field_and_references_exactly": True,
            "observed_particle_masses_or_spacetime_derived": False,
            "interacting_quantum_gauge_Higgs_or_local_chiral_regulator_constructed": False,
            "native_records_select_the_trial_transport_or_Yukawa_law": False,
            "physical_service_errors_or_finite_shot_uncertainty_supplied": False,
            "operator_admission_certifies_frequency_uncertainty_or_a_physical_pole": False}}, arrays


def source_bindings():
    root = Path(__file__).resolve().parents[3]
    from ..record_observer import SOURCE_PATHS
    names = set(SOURCE_PATHS)|set(action.report()["source_bindings"])
    names.update(compiler.source_summary()["source_sha256"])
    names.update({"model/common_source/chiral_candidate/"+name for name in (
        "protected_momentum.py", "test_protected_momentum.py", "README.md", "action.py",
        "gaussian_source.py", "history.py", "kinetic.py", "__init__.py")})
    names.update("model/common_source/mass_response/"+name for name in (
        "model.py", "history.py", "experiment.py", "__init__.py"))
    paths = {root/name for name in names}
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def canonical_report(report):
    """Retain fixed admissions; raw backend residuals remain fresh diagnostics.

    This does not authenticate stored microscopic residual values. Actual
    complex amplitudes, probabilities, decoded readouts and discrete costs
    remain evidence and are compared with a fresh complete replay.
    """
    result = deepcopy(report)
    result["backend_contract"] = {"raw_backend_residuals_authenticated": False,
        "source_output_comparison_absolute_tolerance": 2e-9,
        "decoded_observable_comparison_absolute_tolerance": 2e-7,
        "retained_admissions_checked_against_owned_fresh_execution": True,
        "comparison_tolerances_are_physical_error_bounds": False}
    relation = result["source_interval_relation"]
    value = relation.pop("all_Fock_column_doubled_word_operator_error")
    if not math.isfinite(value) or not 0 <= value <= 2e-8:
        raise ArithmeticError("doubled source-word diagnostic exceeds admission")
    relation["all_Fock_column_doubled_word_operator_error_admission"] = 2e-8
    current = result["same_field_current_attachment"]
    for key in tuple(current):
        if key.endswith("residual"):
            value = current.pop(key)
            if not math.isfinite(value) or not 0 <= value <= 2e-10:
                raise ArithmeticError("current diagnostic exceeds admission")
            current[key+"_admission"] = 2e-10
    residuals = result.pop("Gaussian_and_full_Fock_factorization_residuals")
    if any(not math.isfinite(value) or not 0 <= value <= 2e-9 for rows in residuals.values() for value in rows):
        raise ArithmeticError("complete Gaussian/Fock factorization exceeds admission")
    result["Gaussian_and_full_Fock_factorization_admission"] = 2e-9
    for factor in result["complete_complex_tensor_factors"].values():
        value = factor.pop("complete_source_operator_error")
        if not math.isfinite(value) or not 0 <= value < .01:
            raise ArithmeticError("complete local source factor exceeds admission")
        factor["complete_source_operator_error_admission"] = .01
    for budget in result["protection_and_service_budgets"].values():
        for key, admission in (("source_field_tensor_operator_error_per_cycle", .01),
                               ("source_clock_preparation_vector_error", 1e-10)):
            value = budget.pop(key)
            if not math.isfinite(value) or not 0 <= value <= admission:
                raise ArithmeticError("complete source budget diagnostic exceeds admission")
            budget[key+"_admission"] = admission
        values = budget.pop("source_detector_complete_operator_errors")
        if any(not math.isfinite(v) or not 0 <= v <= 1e-10 for v in values):
            raise ArithmeticError("complete source detector exceeds admission")
        budget["source_detector_complete_operator_error_admission"] = 1e-10
        export = budget["proper_code_export_diagnostic"]
        for key in ("maximum_complement_operator_norm", "maximum_success_norm_excess"):
            value = export.pop(key)
            if not math.isfinite(value) or not 0 <= value <= 2e-11:
                raise ArithmeticError("source export exceeds admission")
            export[key+"_admission"] = 2e-11
        for row in budget["laws"].values():
            epsilon = row["complete_flagged_operator_error_diagnostic"]
            if not math.isfinite(epsilon) or not 0 <= epsilon < .01 or row["fixed_complete_operator_admission"] != .01:
                raise ArithmeticError("complete630-mode source output exceeds admission")
            for key in tuple(row):
                if key.endswith("_diagnostic"):
                    del row[key]
            row["complete_cq_half_diamond_error_admission"] = .01*(1+.01/2)
    return json.loads(json.dumps(result, allow_nan=False))


def _finite_json(value):
    # strict_json rejects duplicate keys; this also rejects JSON NaN/Infinity.
    json.dumps(value, allow_nan=False)


def verify_retained(path):
    """Check original custody and every compact output/probability relation."""
    from tools import strict_json
    path = Path(path)
    document = strict_json.load(path)
    _finite_json(document)
    if document.get("schema") != "oph.chiral_candidate.protected_momentum.v1":
        raise ValueError("wrong protected momentum schema")
    if document.get("source_bindings") != source_bindings():
        raise ValueError("protected momentum source custody changed")
    artifact = document["artifact"]
    if Path(artifact["filename"]).name != artifact["filename"]:
        raise ValueError("protected quantum artifact must be adjacent")
    array_path = path.with_name(artifact["filename"])
    if hashlib.sha256(array_path.read_bytes()).hexdigest() != artifact["sha256"]:
        raise ValueError("protected quantum artifact custody changed")
    with np.load(array_path, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError("duplicate protected quantum array names")
        arrays = dict(archive)
    expected = set(document["complete_complex_tensor_factors"])
    expected.update(name+suffix for name in document["protocols"]
                    for suffix in ("_probe_vectors", "_clock_reference_columns"))
    if set(arrays) != expected or any(a.dtype.kind != "c" or not np.isfinite(a).all() for a in arrays.values()):
        raise ValueError("incomplete or nonfinite protected quantum factors")
    for name, factor in document["complete_complex_tensor_factors"].items():
        array = arrays[name]
        if (array.shape != (3, 1 << factor["Fock_modes"])
                or hashlib.sha256(array.tobytes()).hexdigest() != factor["complex_factor_sha256"]):
            raise ValueError("retained complete Fock factor changed")
    for name, protocol in document["protocols"].items():
        rows = protocol["terminal"]
        probes, clocks = arrays[name+"_probe_vectors"], arrays[name+"_clock_reference_columns"]
        selected, interval = (int(part[1:]) for part in name.split("_"))
        if (len(rows) != 40 or probes.shape != (40, 16 if selected == 0 else 256)
                or clocks.shape != (40, 4, 4)):
            raise ValueError("complete stopped clock/reference output changed")
        for row, probe, clock in zip(rows, probes, clocks, strict=True):
            mass = float(F(row["selector_weight_exact"]))*float(np.vdot(probe, probe).real)*float(np.vdot(clock, clock).real)
            for key, factor in document["complete_complex_tensor_factors"].items():
                if key.startswith(f"t{interval}_") and key != f"t{interval}_p{selected}_f0":
                    vector = arrays[key][row["elapsed_cycles"]]
                    mass *= float(np.vdot(vector, vector).real)
            if not 0 <= row["probability"] <= 1 or not math.isclose(mass, row["probability"], rel_tol=2e-8, abs_tol=2e-10):
                raise ValueError("protected complete quantum output and recorded probability disagree")
        if abs(math.fsum(row["probability"] for row in rows)-1) > 2e-7:
            raise ValueError("complete stopped protected instrument lost mass")
    return document, arrays


def compare_execution(retained, stored_arrays, current, current_arrays):
    """Compare all complex outputs; validate hashes against stored bytes only."""
    left, right = deepcopy(retained), canonical_report(current)
    right["source_bindings"] = source_bindings()
    right["artifact"] = left["artifact"]
    for report in (left, right):
        for factor in report["complete_complex_tensor_factors"].values():
            factor.pop("complex_factor_sha256")
    def compare(a, b, path="root"):
        if type(a) is not type(b):
            raise ValueError("protected receipt type mismatch: "+path)
        if type(b) is dict:
            if a.keys() != b.keys():
                raise ValueError("protected receipt fields mismatch: "+path)
            for key in b:
                compare(a[key], b[key], path+"."+key)
        elif type(b) is list:
            if len(a) != len(b):
                raise ValueError("protected receipt length mismatch: "+path)
            for index, (x, y) in enumerate(zip(a, b, strict=True)):
                compare(x, y, path+f"[{index}]")
        elif type(b) is float:
            if not math.isfinite(a) or not math.isfinite(b):
                raise ValueError("nonfinite protected evidence: "+path)
            if "admission" in path or path.endswith("comparison_absolute_tolerance"):
                if a != b:
                    raise ValueError("protected numerical admission changed: "+path)
            elif not math.isclose(a, b, rel_tol=2e-7, abs_tol=2e-7):
                raise ValueError("protected numerical evidence differs: "+path)
        elif a != b:
            raise ValueError("protected exact evidence differs: "+path)
    compare(left, right)
    if stored_arrays.keys() != current_arrays.keys():
        raise ValueError("protected quantum output keys changed")
    for name, actual in current_arrays.items():
        before = stored_arrays[name]
        if (before.shape != actual.shape or before.dtype.kind != actual.dtype.kind
                or not np.isfinite(before).all() or not np.isfinite(actual).all()
                or not np.allclose(before, actual, rtol=2e-9, atol=2e-9)):
            raise ValueError("complete complex protected output differs: "+name)


def check(path):
    retained, arrays = verify_retained(path)
    current, fresh_arrays = execute(steps=retained["source_interval_relation"]["short_steps"])
    compare_execution(retained, arrays, current, fresh_arrays)
    print("Recomputed protected630-mode field and checked every complete complex factor, record and admission", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--steps", type=int, default=16)
    args = parser.parse_args()
    if args.check:
        if args.output is not None:
            parser.error("checking does not write outputs")
        check(args.check)
        return
    if args.output is None:
        parser.error("--output or --check is required")
    artifact = args.output.with_suffix(".npz")
    if args.output.exists() or artifact.exists():
        parser.error("choose fresh JSON and companion NPZ output paths")
    report, arrays = execute(steps=args.steps)
    report = canonical_report(report)
    report["source_bindings"] = source_bindings()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    report["artifact"] = {"filename": artifact.name,
        "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
        "contract": "all complete complex Fock tensor factors and all40 probe/clock/reference outputs per protocol; computed export failures have norm bounds, not retained vectors"}
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"output": str(args.output), "artifact_bytes": artifact.stat().st_size,
        "decoded_energy_over_clock": {key: value["energy_over_clock"] for key, value in report["decoded"].items()},
        "dispersion": report["dispersion_readout"]}), flush=True)


if __name__ == "__main__":
    main()
