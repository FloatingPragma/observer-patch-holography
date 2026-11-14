"""Executed source-port compilation of two exactly invariant field sectors.

This is an added coherent-controller architecture, not the passive simulator
or a universality claim for Observer.repair. No target exponential is installed
as the implemented gate: full six-level source generators are exponentiated.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
from itertools import cycle
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RER = HERE.parents[1]

import numpy as np
import sympy as s
from scipy.linalg import expm
from model.coherent_defect_wave.chiral_higgs.model import quantum_lepton_block
from model.native_clock import configuration
from model.native_interface import pointer
from model.native_interface import construction as native_interface
from model.observer import native
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler as source


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def complex_array(array):
    return [[[float(z.real), float(z.imag)] for z in row]
            for row in np.atleast_2d(array)]


def load_link():
    """Load the original producer; its sibling current import stays original."""
    directory = RER / 'code/sm_fermion_current'
    previous = sys.modules.get('current')
    if previous is not None and Path(previous.__file__).resolve() != (directory / 'current.py').resolve():
        raise ValueError('a foreign quantum-link current module is already loaded')
    sys.path.insert(0, str(directory))
    try:
        spec = importlib.util.spec_from_file_location('joined_quantum_link', directory / 'quantum_link.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(directory))
    return module


@lru_cache(maxsize=1)
def yukawa_code():
    """Prove HV=Vh and Q_a V=0 using the complete native 320-state source."""
    full, charges, support, _, ward = quantum_lepton_block()
    vacuum = s.zeros(320, 1)
    vacuum[0] = 1
    bright = -full * vacuum / 2
    code = s.ImmutableMatrix(s.Matrix.hstack(vacuum, bright))
    small = s.ImmutableMatrix([[0, -2], [-2, 1]])
    require(code.H * code == s.eye(2), 'Yukawa code is not an isometry')
    require(full * code == code * small, 'Yukawa code is not invariant')
    require(all(q * code == s.zeros(320, 2) for q in charges), 'Yukawa code is not gauge singlet')
    require(ward == 12 and len(support) == 5, 'native source changed')
    return code, small, support


def source_pulse(weights):
    """One simultaneous addressed port pulse; coherent code service is supplied.

    The existing exact source checker verifies all ten coefficient identities
    in Q(sqrt(5)). We combine those *coefficients*, then exponentiate the full
    six-level response. The requested target never enters this execution.
    """
    if not weights or any(k not in ('X', 'Y', 'Z', 'I', 'J') for k in weights):
        raise ValueError('nonempty supported source-control weights required')
    if any(isinstance(v, bool) or not np.isfinite(float(v)) for v in weights.values()):
        raise ValueError('finite real integrated control weights required')
    coefficients = sum((source.integrated_controls(k, float(v)) for k, v in weights.items()), np.zeros(12))
    generators = source._source()[0]  # Pinned actual full source, not logical matrices.
    generator = np.tensordot(coefficients, generators, axes=(0, 0))
    full = expm(generator)
    return full, {
        'kind': 'simultaneous_source_port_pulse',
        'base_weights': {k: str(v) for k, v in weights.items()},
        'integrated_port_coefficients': coefficients.tolist(),
        'proper_code_load_events': 1, 'proper_code_export_events': 1,
        'source_pulses': 1,
        'generator_operator_norm_diagnostic': float(np.linalg.norm(generator, 2)),
        'raw_unitarity_defect': float(np.linalg.norm(full.conj().T @ full - np.eye(6), 2)),
        'valid_code_leakage_diagnostic': float(np.linalg.norm(full[2:, :2], 2)),
        'complete_export_kraus_count': len(source.transfer_kraus(2)),
        'failure_outcomes': ['export_complement_2', 'export_complement_3',
                             'export_complement_4', 'export_complement_5'],
        'exact_valid_code_failure_probability': 0,
        'physical_service_error_measured': False,
    }


def yukawa_pulse(time=F(1, 3)):
    """exp(-it h), h=I/2-2X-Z/2; no splitting or dropped scalar phase."""
    if type(time) not in (int, F):
        raise ValueError('exact rational model time required')
    return source_pulse({'I': -time / 2, 'X': 2 * time, 'Z': time / 2})


def coefficient_precision(bits=32, time=F(1, 3), port_nudges=None):
    """Exact rational per-port quantization enclosure plus complete raw replay.

    All bounds below concern the exact supplied response generators and
    quantized controls. Binary64 exponentiation is reported separately.
    No experimentally authenticated hardware error is supplied by this bound.
    """
    if type(bits) is not int or not 8 <= bits <= 128 or type(time) is not F:
        raise ValueError('8..128 precision bits and a Fraction model time required')
    port_nudges = {} if port_nudges is None else port_nudges
    if type(port_nudges) is not dict or any(type(k) is not int or not 0 <= k < 12 or type(v) is not int
                                          for k, v in port_nudges.items()):
        raise ValueError('integer port-address to integer grid-nudge mapping required')
    scale = 1 << (bits + 32)
    sqrt_lo = F(math.isqrt(5 * scale * scale), scale)
    sqrt_hi = sqrt_lo + F(1, scale)

    def interval(pair):
        a, b = map(F, pair)
        ends = a + b * sqrt_lo, a + b * sqrt_hi
        return min(ends), max(ends)

    weights = {'I': -time / 2, 'X': 2 * time, 'Z': time / 2}
    selected = source._source()[2]
    exact = [tuple(sum((weights[k] * selected[k][p][part] for k in weights), F())
                   for part in range(2)) for p in range(12)]
    quantized, radii = [], []
    for port, pair in enumerate(exact):
        lo, hi = interval(pair)
        grid_value = round((lo + hi) * (1 << bits) / 2)
        value = F(grid_value + port_nudges.get(port, 0), 1 << bits)
        quantized.append(value)
        radii.append(max(abs(value - lo), abs(value - hi)))
    response = json.loads((RER / source.RESPONSE).read_text())
    norm_bounds = []
    for generator in response['generators']:
        bound = F()
        for block in generator:
            for row in block:
                for cell in row:
                    # Entrywise absolute sum dominates operator norm. The
                    # complex magnitude is bounded by |real|+|imaginary|.
                    for pair in (cell[:2], cell[2:]):
                        lo, hi = interval(pair)
                        bound += max(abs(lo), abs(hi))
        norm_bounds.append(bound)
    eta = sum((radius * norm for radius, norm in zip(radii, norm_bounds, strict=True)), F())
    generator = np.tensordot(np.array(list(map(float, quantized))), source._source()[0], axes=(0, 0))
    actual = expm(generator)
    ideal, _ = yukawa_pulse(time)
    residual = float(np.linalg.norm(actual - ideal, 2))
    require(residual <= float(eta) + 1e-14, 'control precision enclosure failed')
    state = actual[:, 0]
    export = []
    success = None
    for index, kraus in enumerate(source.transfer_kraus(2)):
        value = kraus @ state
        rho = np.outer(value, value.conj())
        if index == 0:
            success = rho[:2, :2]
        export.append({'record': 'export_success' if index == 0 else f'export_complement_{index + 1}',
                       'probability': float(np.trace(rho).real), 'output_code_plus_failure_density': complex_array(rho)})
    pointer_branches = fine_pointer(success, 'C')
    total = sum(x['probability'] for x in export[1:]) + sum(x['probability'] for x in pointer_branches)
    require(abs(total - 1) < 1e-12, 'quantized complete output lost probability')
    return {
        'precision_bits_per_integrated_port_coefficient': bits,
        'explicit_extra_grid_nudges': {str(k): v for k, v in port_nudges.items()},
        'sqrt5_exact_rational_enclosure': [str(sqrt_lo), str(sqrt_hi)],
        'exact_coefficients_Qsqrt5': [[str(a), str(b)] for a, b in exact],
        'quantized_coefficients': list(map(str, quantized)),
        'coefficient_error_upper': list(map(str, radii)),
        'generator_operator_norm_upper': list(map(str, norm_bounds)),
        'single_pulse_operator_error_upper': str(eta),
        'complete_export_and_pointer_half_diamond_error_upper_from_quantization': str(eta),
        'error_upper_decimal_display': float(eta),
        'source_roundoff_residual_separate_diagnostic': residual,
        'complete_export_outputs': export,
        'C_pointer_outputs_on_subnormalized_success': [{**b, 'output': complex_array(b['output'])} for b in pointer_branches],
        'complete_raw_output_mass': total,
        'export_failure_work_retained': True,
        'no_success_renormalization': True,
        'exact_proof': 'unitary Duhamel: ||exp(sum cprime D)-exp(sum c D)|| <= sum_p ||D_p|| |cprime_p-c_p|; complete export/readout are contractions',
        'physical_generator_drift_or_code_service_error_included': False,
    }


def fine_pointer(rho, law='C'):
    """Actual six-pulse source coupling followed by every native C/E record.

    rho may include a reference as its second tensor factor. Auxiliary output
    configurations are Kraus-environment indices, not extra observed records.
    Both records and their unnormalized work/reference outputs are retained.
    Called after a conclusive auxiliary reset; reset failures are separate.
    """
    rho = np.asarray(rho, dtype=complex)
    if law not in ('C', 'E') or rho.ndim != 2 or rho.shape[0] != rho.shape[1] or rho.shape[0] % 2:
        raise ValueError('C/E and a work-qubit/reference density required')
    rd = rho.shape[0] // 2
    matrix = pointer.local_matrix()
    groups = defaultdict(lambda: defaultdict(lambda: np.zeros((2, 2), complex)))
    currents = {}
    for z in range(4):
        state = circuit.decode(z)
        coupling = matrix[np.ix_([z, z + 4], [pointer.READY_CODE, pointer.READY_CODE + 4])]
        for coin in (0, 1):
            out, (_, current) = native.native(state, pointer.SEAM, coin)
            label = current if law == 'C' else state[0] - state[1]
            record = (coin, label)
            require(record not in currents or currents[record] == current, 'record/current ambiguity')
            currents[record] = current
            groups[record][out] += coupling / np.sqrt(2)
    branches = []
    completeness = np.zeros((2, 2), complex)
    for (coin, label), environment in sorted(groups.items()):
        output = np.zeros_like(rho)
        for kraus in environment.values():
            full = np.kron(kraus, np.eye(rd))
            output += full @ rho @ full.conj().T
            completeness += kraus.conj().T @ kraus
        current = currents[coin, label]
        branches.append({'record': {'law': law, 'coin': coin, 'label': label, 'current': current,
                                    'bit': int(2 * current + 2 * coin - 1 < 0)},
                         'probability': float(np.trace(output).real), 'output': output,
                         'environment_configurations': [list(x) for x in environment]})
    require(np.linalg.norm(completeness - np.eye(2), 2) < 1e-12, 'incomplete native pointer')
    return branches


def reset_receipt(target_code, seed_scans=4, retries=32):
    """Actual finite adaptive controller, with full success/failure recipe."""
    target = circuit.decode(target_code)
    cert = configuration.certificate(seed_scans, retries, full=True, mode='reset', target=target)
    executions = []
    # Deliberate inputs/tapes are replay examples, never probability samples.
    for law, initial, coins in (
        ('C', native.states()[0], cycle((0, 1, 1, 0))),
        ('E', native.states()[-1], cycle((1, 0))),
        ('C', native.states()[0], cycle((0,))),
    ):
        execution = configuration.execute(initial, coins, seed_scans, retries,
                                          full=True, mode='reset', target=target, law=law)
        support = configuration.word_support(execution, full=True)
        require(initial in support and support[initial] == execution.output, 'reset replay mismatch')
        if execution.outcome.success:
            require(execution.output == target and len(support) == 1, 'reset success not conclusive')
        executions.append({
            'law': law, 'input': list(initial), 'output': list(execution.output),
            'success': execution.outcome.success, 'reason': execution.outcome.reason,
            'attempts': execution.outcome.attempts,
            'fine_record_tape': [[list(e), c, j] for e, c, j in execution.records],
            'native_fine_word_squared_amplitude': str(F(1, 2 ** execution.outcome.attempts)),
            'word_support_columns': len(support),
        })
    return {'certificate': cert, 'replay_examples': executions,
            'all_paths_recipe': {
                'controller': 'model.native_clock.configuration.policy',
                'executor': 'model.native_clock.configuration.execute',
                'kraus_support': 'model.native_clock.configuration.word_support',
                'kraus_formula': '2^(-attempts/2) sum_{x in word_support} |final(x)><x|',
                'termination': 'finite certified attempt cap; every terminal tape is success or explicit timeout',
                'records_retained': ['seam', 'coin', 'C-current or E-contrast', 'timeout reason', 'final quantum work'],
                'recipe_is_complete_enumeration_semantics_not_a_list_of_all_exponentially_many_tapes': True,
                'private_reset_records_retained_in_audit_and_traced_only_in_stated_channel_comparison': True,
            }}


def car_step(mask, mode, create):
    if bool(mask & (1 << mode)) == create:
        return None
    return mask ^ (1 << mode), (-1) ** ((mask & ((1 << mode) - 1)).bit_count())


def link_audit():
    """All 491,520 CAR/Gauss columns plus actual 15-pulse state execution."""
    link = load_link()
    rows = link.channels()
    n = len(rows)
    require(n == 15 and sum(r['integer_charge'] for r in rows) == 0, 'generation changed')
    size = 1 << n
    tape = []
    for row in rows:
        r = F(row['step']) * F(row['kappa']) / 2
        tape.append(source.pulse('X', -2 * np.arctan(float(r)), [row['id']]))
    initial = np.zeros(size, complex)
    initial[0] = 1
    actual = source.execute(tape, n, initial)
    expected = np.ones(size, complex)
    histogram = defaultdict(F)
    car_checks = 0
    for mask in range(size):
        occupied = sum(1 << (2 * i + ((mask >> i) & 1)) for i in range(n))
        moved = sum(row['integer_charge'] for row in rows if mask >> row['id'] & 1)
        ell = -moved
        charge_left, charge_right = -moved, moved
        require(ell - charge_left == 0 and -ell - charge_right == 0, 'Gauss failed')
        exact = link.parent.C(1)
        for row in rows:
            i = row['id']
            src = 2 * i + ((mask >> i) & 1)
            dst = 4 * i + 1 - src
            m, sign1 = car_step(occupied, src, False)
            m, sign2 = car_step(m, dst, True)
            require(sign1 * sign2 == 1, 'CAR hop acquired wrong sign')
            wanted = occupied ^ (1 << src) ^ (1 << dst)
            require(m == wanted, 'CAR occupation mismatch')
            next_mask = mask ^ (1 << i)
            next_flux = -sum(r['integer_charge'] for r in rows if next_mask >> r['id'] & 1)
            require(next_flux - ell == (row['integer_charge'] if mask >> i & 1 else -row['integer_charge']),
                    'electric compensation mismatch')
            exact *= link.amplitudes(row)[(mask >> i) & 1]
            car_checks += 1
        expected[mask] = complex(float(exact.real), float(exact.imag))
        histogram[ell] += exact.norm()
    error = float(np.linalg.norm(actual - expected))
    require(error < 2e-13, 'source word differs from exact quantum-link state')
    census, parent_histogram = link.branch_census(rows)
    require(dict(histogram) == parent_histogram, 'quantum-link histogram mismatch')
    return {
        'encoding': '|mask> -> ordered CAR occupations (2i+bit_i) tensor |ell=-sum q_i bit_i>',
        'encoding_is_declared_memory_semantics_not_a_native_repair_operation': True,
        'sector_dimension': size, 'ambient_link': 'ell2(Z); no cyclic electric cutoff',
        'operator_intertwining_columns': car_checks,
        'every_H_i_restricts_to': 'kappa_i X_i, including CAR sign +1 and compensating electric shift',
        'all_encoded_states_exactly_Gauss_invariant': True,
        'source_tape': tape, 'resources': source.resources(tape),
        'prepared_state_error_norm_diagnostic': error,
        'raw_state_norm_squared': float(np.vdot(actual, actual).real),
        'exact_branch_census': census,
        'exact_electric_readout': link.histogram_readout(histogram),
        'pointer_scope': 'measuring all work bits gives the complete configuration instrument; coarse flux labels do not make it the flux-Lueders instrument',
        'continuous_hamiltonian_compiler': 'replace each angle by -t*kappa_i; all 15 encoded X_i commute',
        'Cayley_global_operator_error_upper_vs_time_step': 'sum_i abs(2 atan(step*kappa_i/2)-step*kappa_i)',
        'Cayley_error_at_shipped_step': float(sum(abs(2*np.arctan(float(F(r['step'])*F(r['kappa'])/2))
                                                  -float(F(r['step'])*F(r['kappa']))) for r in rows)),
        'physical_time_calibrated': False,
        'physical_electric_memory_or_15_qubit_hardware_built': False,
        'logical_fifteen_bit_blank_is_supplied_for_this_execution': True,
        'fifteen_bit_blank_not_derived_from_one_78_state_reset': True,
    }


def repair_boundary():
    """Read the actual Observer repair block, without inferring gate control."""
    from model.coherent_seam_repair.program import Observer
    observer = Observer()
    a, b = observer.addresses[0]
    coordinates = [(a, 0), (b, 0), (a, 1), (b, 1)]
    actual = s.zeros(4)
    for j, (port, parity) in enumerate(coordinates):
        output = observer.repair(observer.basis_state(port, parity, 0))
        for i, (p, r) in enumerate(coordinates):
            actual[i, j] = output.amplitudes[observer.index(p, r, 0)]
    signs = s.Matrix([1, -1, -1, -1])
    require(actual == s.eye(4) - signs * signs.T / 2, 'Observer repair block changed')
    require(actual * actual == s.eye(4) and actual == actual.conjugate(), 'repair reflection lost its real involution')
    audits = {law: native_interface.exact_native_audit(law) for law in ('C', 'E')}
    return {
        'actual_Observer_repair_block': [[str(x) for x in row] for row in actual.tolist()],
        'real_involution': True,
        'actual_Observer_address_reflection': '2 AA^dagger-I, uniform 60-address A; a separately declared coherent operation',
        'optional_Observer_quarter_turn': 'address-selective multiplication by i; optional apparatus, not idle repair',
        'native_committed_C_E_audits': audits,
        'exact_real_control_counterexample': {
            'target_Cayley_state': ['15/17', '-8i/17'],
            'density_off_diagonal': '120i/289',
            'real_only_words_on_same_real_code_cannot_produce_it': True,
            'alternative_rebit_encoding_excluded': False,
            'quarter_turn_plus_reflections_universality_claimed': False,
        },
        'minimum_claim': 'two independently stated resources are used: complex port actuation and reference-preserving code transfer; no strict primitive-count lower bound',
    }


def run():
    code, h, support = yukawa_code()
    full, receipt = yukawa_pulse()
    target = expm(-1j * np.asarray(h.tolist(), complex) / 3)
    error = float(np.linalg.norm(full[:2, :2] - target, 2))
    require(error < 1e-13, 'source Yukawa pulse mismatch')
    state = full[:2, 0]
    rho = np.outer(state, state.conj())
    all_branches = {}
    for law in ('C', 'E'):
        branches = fine_pointer(rho, law)
        all_branches[law] = [{**b, 'output': complex_array(b['output'])} for b in branches]
    source_files = [
        'model/coherent_seam_repair/program.py', 'model/observer/native.py',
        'model/native_interface/construction.py', 'model/native_interface/pointer.py',
        'model/native_clock/configuration.py', 'model/observer/source_discrimination/circuit.py',
        'model/observer/source_scalar/compiler.py', 'code/source_selection_model/response.json',
        'code/m1_source_realization/receipt.json', 'code/sm_fermion_current/quantum_link.py',
        'code/sm_fermion_current/current.py', 'model/coherent_defect_wave/chiral_higgs/model.py',
        'model/coherent_defect_wave/physical_attachment/chiral.py',
    ]
    return {
        'schema': 'oph.join.compiled_field_services.v1',
        'scope': {
            'added_coherent_port_controller': True,
            'derived_from_Observer_real_reflections': False,
            'passive_classical_run_is_this_program': False,
            'whole_primary_nonabelian_spatial_program_compiled_and_run': False,
            'physical_service_noise_authenticated': False,
        },
        'smallest_sufficient_interface_used': {
            'operation': 'coherently load one/pair memory code; apply exp(sum c_p D_p); complete export with explicit complement failures',
            'continuous_port_actuation': 'declared architecture; supplied real c_p and precision',
            'coherent_proper_code_transfer': 'additional CCG requirement; source pulse alone does not prove it',
            'address_controller': 'finite retained classical schedule with one- and two-code routing',
            'strict_resource_minimality_proved': False,
            'additional_Householder_or_target_matrix_oracle_used': False,
        },
        'source_controls': source.source_summary(),
        'actual_repair_boundary': repair_boundary(),
        'finite_precision_execution': coefficient_precision(),
        'finite_precision_leakage_control': coefficient_precision(port_nudges={0: 1}),
        'native_yukawa': {
            'full_source_dimension': 320, 'full_CAR_modes': 6, 'scalar_levels': 5,
            'vacuum_cyclic_support': support,
            'code_basis': ['|0>', '(|53>+|64>+|88>+|169>)/2'],
            'restricted_H': [[str(x) for x in row] for row in h.tolist()],
            'exact_isometry': True, 'exact_H_intertwining': True,
            'exact_gauge_charge_annihilations': 24,
            'individual_four_leaf_basis_states_are_not_all_gauge_singlets': True,
            'model_time': '1/3', 'source_pulse': receipt,
            'source_vs_target_operator_error_diagnostic': error,
            'raw_output_state': complex_array(state.reshape(-1, 1)),
            'native_C_and_E_complete_pointer_outputs': all_branches,
            'pointer_certification': {k: str(v) if isinstance(v, F) else v for k, v in pointer.certificate().items()},
            'reset_work': reset_receipt(0), 'reset_auxiliary': reset_receipt(pointer.READY_CODE),
            'successful_service_source_pulses': 7,
            'successful_service_proper_code_events': 14,
            'readout_native_attempts_after_reset': 1,
            'scalar_occupation_readout': 'N_scalar restricts to |1><1|; complete bit instrument preserves arbitrary references',
            'no_source_state_renormalization': True,
            'no_Trotter_error_on_this_sector': True,
        },
        'compensated_quantum_link': link_audit(),
        'error_contract': {
            'metric': 'half diamond distance for complete output instrument, with stated marginalization of private reset records',
            'bound': 'min(1, tau_work + tau_aux + sum_pulse eta_op + sum_service delta_halfdiamond + delta_native_read)',
            'pulse_errors': 'uniform physical operator norm for actual unitary controls; numerical residual is a separate diagnostic',
            'control_precision_bound': 'eta_op <= sum_p ||D_p|| |delta c_p| by unitary Duhamel',
            'records': 'success, every native coin/current or contrast, export failures, and reset timeouts; no conditioning denominator',
        },
        'source_sha256': {name: hashlib.sha256((RER / name).read_bytes()).hexdigest() for name in source_files},
        'implementation_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in (HERE / 'construction.py', HERE / 'finite_compiler.py', HERE / 'test_compiler.py')},
    }
