"""Finite reciprocal geometry--field heat bath with explicit energy records.

Compatibility is specified independently: the joint Gibbs state of two
declared finite field Hamiltonians, at a supplied inverse temperature.
The primary case uses vacuum_repair's SAME three-oscillator D2, lambda6
Hamiltonians at r=9/10,11/10. Native C/E events on an auxiliary copy of the
two-state geometry orbit propose a geometry and a field-energy label.
The system geometry is changed only by an energy-balanced bath collision.
Rejected proposals and bath outputs remain recorded.

For each allowed cross-geometry pair i,j there is a fixed-gap thermal bath
qubit. Swapping |low,1> with |high,0> conserves joint system+bath energy.
Its upward/downward probabilities obey exact detailed balance with the
same joint Gibbs state. The resulting energy-basis measurement/feedback
channel has that full quantum state as a fixed point; its population law
is irreducible on the explicitly admitted sector. Source controls execute
one complete collision and its purified thermal input, without replacing
source arrays by the target matrices.

This replaces a driven uniform geometry/vacuum reset by reciprocal heat
exchange. It supplies no isolated autonomous cooling: beta, the Hamiltonian,
fixed-gap reservoir bank, energy readout, coherent controls and fresh thermal
bath inputs are premises. The regional bare-area stationarity obstruction
is not resolved by thermodynamic probability stationarity. No metric or
Einstein law is inferred. The optional fermion case is fixed-N=1 thermal
matter, not a vacuum or a grandcanonical attraction claim.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, eigh
import sympy as s

from tools import strict_json
from . import vacuum_repair
from ..observer import native
from ..observer.source_discrimination import circuit, experiment as native_source
from ..observer.source_scalar import compiler


ROOT = Path(__file__).resolve().parents[2]
GEOMETRIES = (F(9, 10), F(11, 10))
BETA = 1/s.sqrt(14)
SOURCE_ADMISSION = 2e-8
REPLAY_TOLERANCE = 2e-9


def _freeze(value):
    value = np.asarray(value, complex)
    return np.frombuffer(value.tobytes(), complex).reshape(value.shape)


def default_hamiltonians():
    """Exact genuine P_D H P of the same D2 lambda6 scalar control."""
    terms = tuple(s.Matrix(term.astype(int)) for term in vacuum_repair._terms(2))
    return tuple(s.ImmutableMatrix(sum((coefficient*term for coefficient, term in
        zip(vacuum_repair.coefficients(2, r), terms, strict=True)), s.zeros(8))) for r in GEOMETRIES)


@dataclass(frozen=True)
class FiniteModel:
    exact_hamiltonians: tuple
    beta: s.Expr
    field_dimension: int
    hamiltonian: np.ndarray
    frame: np.ndarray
    energies: np.ndarray
    charges: tuple
    active: tuple
    weights: np.ndarray
    probabilities: np.ndarray
    gibbs: np.ndarray
    transition: np.ndarray
    spectral_diagnostics: tuple


def finite_model(hamiltonians=None, *, beta=BETA, charge=None, sector=None):
    """One real finite Hamiltonian pair; optional EXACT conserved charge sector."""
    hs = default_hamiltonians() if hamiltonians is None else hamiltonians
    if type(hs) is not tuple or len(hs) != 2:
        raise ValueError('two exact finite Hamiltonians in one fixed basis required')
    hs = tuple(s.ImmutableMatrix(h) for h in hs)
    d = hs[0].rows
    if (d not in (2, 4, 8) or any(h.shape != (d, d) or h != h.T or h.free_symbols
            or h.atoms(s.Float) or any(x.is_real is not True for x in h) for h in hs)):
        raise ValueError('real exact symmetric field Hamiltonians of dimension2,4,8 required')
    if isinstance(beta, (bool, np.bool_)):
        raise ValueError('an exact positive inverse temperature required')
    b = s.sympify(beta)
    if b.free_symbols or b.atoms(s.Float) or b.is_positive is not True or b.is_finite is not True:
        raise ValueError('a declared exact positive inverse temperature required')
    if charge is None:
        if sector is not None:
            raise ValueError('a fixed sector requires an explicit conserved charge')
        labels = (0,)*d
    else:
        q = s.ImmutableMatrix(charge)
        if (q.shape != (d, d) or q != s.diag(*q.diagonal()) or any(x.is_Integer is not True for x in q.diagonal())
                or any(h*q != q*h for h in hs) or type(sector) is not int):
            raise ValueError('an exact diagonal conserved integer charge and original fixed sector required')
        labels = tuple(map(int, q.diagonal()))
        if sector not in labels:
            raise ValueError('the admitted charge sector is empty')
    values, frames, charges, spectral_diagnostics = [], [], [], []
    for h in hs:
        numeric = np.asarray(h.evalf(20), float)
        columns = []
        for q in sorted(set(labels)):
            indices = [i for i, value in enumerate(labels) if value == q]
            block = numeric[np.ix_(indices, indices)]
            exact_block = h.extract(indices, indices)
            if exact_block == s.diag(*exact_block.diagonal()):
                order = sorted(range(len(indices)), key=lambda i: (block[i, i], i))
                eigenvalues, eigenvectors = np.diag(block)[order], np.eye(len(indices))[:, order]
                gap = None
            else:
                eigenvalues, eigenvectors = eigh(block)
                gap = float(min(np.diff(eigenvalues))) if len(eigenvalues) > 1 else None
                if gap is not None and gap <= 1e-8*max(1., np.linalg.norm(block, 2)):
                    raise ValueError('degenerate/unresolved nondiagonal spectral sector; supply finer exact symmetry sectors')
            pivots = []
            for energy, vector in zip(eigenvalues, eigenvectors.T, strict=True):
                full = np.zeros(d)
                full[indices] = vector
                pivot = int(np.flatnonzero(abs(full) > 1e-8)[0])
                if abs(full[pivot]) < 1e-6:
                    raise ValueError('spectral phase pivot is numerically unresolved')
                if full[pivot] < 0:
                    full = -full
                pivots.append(pivot)
                columns.append((float(energy), q, full))
            spectral_diagnostics.append({'geometry_block': len(frames), 'charge_label': q,
                'sector_dimension': len(indices), 'minimum_resolved_sector_gap': gap,
                'phase_rule': 'first computational component above1e-8 is positive, with magnitude at least1e-6',
                'phase_pivots': pivots, 'diagonal_degeneracy_uses_fixed_computational_basis': gap is None})
        columns.sort(key=lambda row: (row[0], row[1]))
        frame = np.column_stack([row[2] for row in columns])
        energy = np.array([row[0] for row in columns])
        if (np.linalg.norm(numeric@frame-frame*energy, 2) > 1e-10
                or np.linalg.norm(frame.T@frame-np.eye(d), 2) > 1e-10):
            raise ArithmeticError('declared finite Hamiltonian spectral frame is unresolved')
        frames.append(frame); values.extend(energy); charges.extend(row[1] for row in columns)
    energies = np.asarray(values)
    active = tuple(i for i, q in enumerate(charges) if sector is None or q == sector)
    weights = np.zeros(2*d)
    if not math.isfinite(float(b)) or float(b)*(max(energies[list(active)])-min(energies[list(active)])) > 80:
        raise ValueError('declared finite thermal spectral window is unresolved numerically')
    weights[list(active)] = np.exp(-float(b)*(energies[list(active)]-min(energies[list(active)])))
    probabilities = weights/weights.sum()
    frame = block_diag(*frames)
    h = block_diag(*(np.asarray(value.evalf(20), float) for value in hs))
    transition = np.eye(2*d)
    for i in active:
        for j in active:
            if i//d == j//d:
                continue
            probability = weights[j]/(weights[i]+weights[j])/(2*d)
            transition[j, i] = probability
            transition[i, i] -= probability
    if (min(np.diag(transition)) <= 0 or np.linalg.norm(transition.sum(axis=0)-1, np.inf) > 1e-12
            or np.linalg.norm(transition@probabilities-probabilities, np.inf) > 1e-12):
        raise ArithmeticError('derived reciprocal heat-bath population law failed')
    gibbs = frame@np.diag(probabilities)@frame.T
    return FiniteModel(hs, b, d, _freeze(h), _freeze(frame), np.asarray(energies), tuple(charges),
        active, np.asarray(weights), np.asarray(probabilities), _freeze(gibbs), np.asarray(transition), tuple(spectral_diagnostics))


def exact_balance_identity():
    """Exact finite proof, independent of numerical eigenvalues/temperature."""
    wi, wj, proposal = s.symbols('wi wj proposal', positive=True)
    forward, reverse = proposal*wj/(wi+wj), proposal*wi/(wi+wj)
    if s.cancel(wi*forward-wj*reverse) != 0:
        raise ArithmeticError('finite Gibbs pair detailed balance failed')
    e0, e1 = s.symbols('E_low E_high', real=True)
    # Order: systemlow/high is low bit, bath0/1 is high bit.
    total = s.diag(e0, e1, e0+(e1-e0), e1+(e1-e0))
    swap = s.eye(4); swap[1, 1] = swap[2, 2] = 0; swap[1, 2] = swap[2, 1] = 1
    if swap*total != total*swap or swap.T*swap != s.eye(4):
        raise ArithmeticError('thermal collision does not conserve joint energy')
    return {'population_detailed_balance_exact': True,
        'identity': 'w_i*q*w_j/(w_i+w_j)=w_j*q*w_i/(w_i+w_j), w_i=exp(-beta*E_i)',
        'collision_energy_commutator_exact': True,
        'full_quantum_Gibbs_fixed_point_reason': '[U,H_system+H_bath]=0 implies U exp(-beta*(H_system+H_bath)) U.H is unchanged',
        'nearest_theorems': ['https://arxiv.org/abs/0911.3635', 'https://arxiv.org/abs/2311.09207'],
        'general_noncommutative_Gibbs_sampler_claimed': False}


def proposal_records(law):
    """Four ACTUAL native proposal coins; all16 histories have weight1/16.

    The proposal register copies the already measured classical geometry.
    It is auxiliary: an unbalanced primitive is never applied first to the
    system geometry. Its Hamiltonian/record service are declared degenerate.
    """
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native C or E proposal law required')
    states = []
    for bit in (0, 1):
        loads = [0]*12; loads[bit] = loads[5] = 1
        states.append(tuple(loads))
    records = np.zeros((2, 16, 4, 3), np.int64)
    destinations = np.zeros((2, 16, 2), np.int64)
    for initial in (0, 1):
        for address in range(16):
            state, bits = states[initial], []
            for event in range(4):
                coin = (address >> event) & 1
                after, (_, current) = native.native(state, (0, 1), coin)
                label = current if law == 'C' else state[0]-state[1]
                final = states.index(after)
                if final != 1-coin:
                    raise ArithmeticError('native geometry proposal is not its declared fair orbit')
                records[initial, address, event] = coin, label, final
                bits.append(final); state = after
            destinations[initial, address] = bits[0], sum(bits[k+1] << k for k in range(3))
        if {tuple(row) for row in destinations[initial]} != {(g, m) for g in (0, 1) for m in range(8)}:
            raise ArithmeticError('native proposal records omit a geometry/field label')
    return records, destinations


def records(model, law='C'):
    """Complete accepted/rejected finite histories and heat, without selection."""
    native_records, proposed = proposal_records(law)
    rows, weights, heat = [], [], []
    d = model.field_dimension
    if d != 8:
        raise ValueError('this literal native proposal grammar uses the fixed8-level field')
    for i in model.active:
        g = i//d
        for address, (h, m) in enumerate(proposed[g]):
            j = int(h*d+m)
            allowed = h != g and j in model.active
            if allowed:
                low, high = sorted((i, j), key=lambda state: (model.energies[state], state))
                p1 = model.weights[high]/(model.weights[low]+model.weights[high])
                gap = model.energies[high]-model.energies[low]
            else:
                low = high = i; p1 = .5; gap = 0.
            for before in (0, 1):
                accepted = allowed and (before == int(i == low))
                final = j if accepted else i
                after = 1-before if accepted else before
                probability = (p1 if before else 1-p1)/16
                q = float(model.energies[final]-model.energies[i])
                if abs(q+gap*(after-before)) > 2e-12:
                    raise ArithmeticError('retained transition has an unaccounted work term')
                rows.append((i, address, int(h), int(m), before, after, final, int(accepted)))
                weights.append(probability); heat.append(q)
    rows, weights, heat = np.asarray(rows, np.int64), np.asarray(weights), np.asarray(heat)
    reconstructed = np.zeros_like(model.transition)
    for row, weight in zip(rows, weights, strict=True):
        reconstructed[row[6], row[0]] += weight
    if np.linalg.norm(reconstructed[:, model.active]-model.transition[:, model.active], np.inf) > 1e-12:
        raise ArithmeticError('complete native/bath record law differs from derived heat bath')
    return rows, weights, heat, native_records


def apply_channel(model, rho):
    """Full quantum energy measurement plus reciprocal recorded population law."""
    if np.ma.isMaskedArray(rho):
        raise ValueError('masked quantum states are inadmissible')
    value = np.asarray(rho)
    if value.shape != model.gibbs.shape or value.dtype.kind not in 'fc' or not np.isfinite(value).all():
        raise ValueError('finite quantum joint state in the admitted fixed basis required')
    if np.linalg.norm(value-value.conj().T, 2) > 1e-10 or np.linalg.eigvalsh(value)[0] < -1e-10:
        raise ValueError('positive Hermitian quantum state required')
    if abs(np.trace(value)-1) > 1e-10:
        raise ValueError('normalized input state required')
    energy = model.frame.conj().T@value@model.frame
    projector = np.zeros(len(model.energies)); projector[list(model.active)] = 1
    if np.linalg.norm(energy-projector[:, None]*energy*projector[None, :], 2) > 1e-10:
        raise ValueError('input lies outside the explicitly admitted charge sector')
    population = np.diag(energy).real
    return model.frame@np.diag(model.transition@population)@model.frame.conj().T


def _entropy(rho):
    values = np.linalg.eigvalsh(rho)
    if values[0] < -1e-10:
        raise ValueError('nonpositive state in entropy comparison')
    return -math.fsum(float(x)*math.log(float(x)) for x in values if x > 0)


def thermodynamics(model):
    d, beta = model.field_dimension, float(model.beta)
    probabilities = model.probabilities
    geometry = [float(sum(probabilities[g*d:(g+1)*d])) for g in (0, 1)]
    free = [-math.log(sum(model.weights[g*d:(g+1)*d]))/beta for g in (0, 1)]
    # A shared constant energy shift in weights cancels every response.
    potential_derivative = -beta*geometry[0]*geometry[1]
    gradient = free[1]-free[0]+math.log(geometry[1]/geometry[0])/beta
    hessian = (1/geometry[0]+1/geometry[1])/beta
    currents = model.transition*probabilities[None, :]
    balance = float(np.linalg.norm(currents-currents.T, np.inf))
    stationary_heat = sum(probabilities[i]*model.transition[j, i]*(model.energies[j]-model.energies[i])
                          for i in model.active for j in model.active)
    naive = probabilities.copy()
    for g in (0, 1):
        naive[g*d:(g+1)*d] /= 2*geometry[g]
    drive_cost = sum(naive[i]*math.log(naive[i]/probabilities[i]) for i in model.active)/beta
    naive_defect = float(np.linalg.norm(model.transition@naive-naive, 1))
    entropy = -sum(probabilities[i]*math.log(probabilities[i]) for i in model.active)
    geometry_entropy = -sum(p*math.log(p) for p in geometry)
    conditional_entropy = sum(-probabilities[i]*math.log(probabilities[i]/geometry[i//d]) for i in model.active)
    active = np.asarray(model.active)
    roots = np.sqrt(probabilities[active])
    reversible = model.transition[np.ix_(active, active)]*roots[None, :]/roots[:, None]
    spectral_gap = 1-float(np.linalg.eigvalsh((reversible+reversible.T)/2)[-2])
    mean_energy = float(probabilities@model.energies)
    geometry_energy_covariance = float(sum(probabilities[i]*(i//d-geometry[1])*(model.energies[i]-mean_energy) for i in model.active))
    if balance > 1e-12 or abs(stationary_heat) > 1e-12 or abs(gradient) > 1e-12 or hessian <= 0:
        raise ArithmeticError('reciprocal stationary energy/free-energy identities failed')
    return {'geometry_probabilities': geometry, 'conditional_free_energies_up_to_one_common_constant': free,
        'geometry_free_energy_first_variation': gradient, 'geometry_free_energy_second_variation': hessian,
        'geometry_probability_response_to_added_branch1_energy': potential_derivative,
        'population_detailed_balance_residual': balance, 'stationary_mean_system_heat': float(stationary_heat),
        'stationary_external_drive_work_in_ideal_collision': 0,
        'joint_quantum_entropy': float(entropy), 'geometry_Shannon_entropy': float(geometry_entropy),
        'mean_conditional_field_entropy': float(conditional_entropy),
        'entropy_chain_rule_residual': float(entropy-geometry_entropy-conditional_entropy),
        'geometry_energy_covariance': geometry_energy_covariance,
        'geometry_probability_beta_response': -geometry_energy_covariance,
        'mean_field_energy_branch1_potential_response': -beta*geometry_energy_covariance,
        'population_gap_diagnostic_on_admitted_sector': spectral_gap,
        'uniform_geometry_conditional_Gibbs_state_stationarity_defect_L1': naive_defect,
        'isothermal_free_energy_cost_of_forcing_that_uniform_state': float(drive_cost),
        'bare_area_or_continuous_metric_entropy_stationarity_claimed': False}


def source_collision(model, edge=None):
    """Literal32-column energy-balanced collision plus retained bath purifier."""
    d = model.field_dimension
    edges = [(i, j) for i in model.active for j in model.active if i < d <= j]
    if edge is None:
        # Fixed lexical choice avoids optimizing a floating gap or breaking
        # ties between symmetry-related energies on different LAPACK hosts.
        candidates = [pair for pair in edges if abs(model.energies[pair[1]]-model.energies[pair[0]]) > 1e-6]
        if not candidates:
            raise ValueError('no resolved nonzero bath gap for the literal source control')
        edge = candidates[0]
        if abs(model.energies[edge[1]]-model.energies[edge[0]]) < 1e-5:
            raise ValueError('source edge gap is too close to its selection threshold')
    i, j = edge
    if (i, j) not in edges:
        raise ValueError('one allowed cross-geometry collision edge required')
    low, high = sorted((i, j), key=lambda state: (model.energies[state], state))
    size = 2*d
    if size != 16:
        raise ValueError('literal source collision fixes the16-level system')
    diagonalize = compiler.compile_real_orthogonal(np.asarray(model.frame.real))
    exchange = circuit.basis_swap_word(low+size, high, 5)
    word = compiler.inverse(diagonalize)+exchange+diagonalize
    actual = compiler.execute(word, 5)
    swap = np.eye(32); swap[:, [low+size, high]] = swap[:, [high, low+size]]
    frame = np.kron(np.eye(2), model.frame)
    target = frame@swap@frame.conj().T
    gap = float(model.energies[high]-model.energies[low])
    total_energy = np.kron(np.eye(2), model.hamiltonian)+np.kron(np.diag([0., gap]), np.eye(size))
    error = float(np.linalg.norm(actual-target, 2))
    commutator = float(np.linalg.norm(actual@total_energy-total_energy@actual, 2))
    endpoint_energy_defect = float(np.linalg.norm(actual.conj().T@total_energy@actual-total_energy, 2))
    if error > SOURCE_ADMISSION or max(commutator, endpoint_energy_defect) > 100*SOURCE_ADMISSION:
        raise ArithmeticError('source thermal collision violates all-input energy/target admission')
    p1 = float(model.weights[high]/(model.weights[low]+model.weights[high]))
    purification = np.array([math.sqrt(1-p1), 0., 0., math.sqrt(p1)])
    bath_word = compiler.remap(compiler.compile_real_state(purification), (4, 5))
    incoming = np.zeros((64, size), complex); incoming[:size] = np.eye(size)
    purified = compiler.execute(bath_word+word, 6, incoming)
    ideal_purified = np.zeros_like(purified)
    for reference, bath in ((0, 0), (1, 1)):
        weight = math.sqrt(1-p1) if bath == 0 else math.sqrt(p1)
        initial = np.zeros((32, size), complex)
        initial[bath*size:(bath+1)*size] = weight*np.eye(size)
        ideal_purified[reference*32:(reference+1)*32] = target@initial
    purified_error = float(np.linalg.norm(purified-ideal_purified, 2))
    if purified_error > SOURCE_ADMISSION:
        raise ArithmeticError('source purified bath/input-reference dilation failed')
    bath = np.diag([1-p1, p1])
    joint = np.kron(bath, model.gibbs)
    invariance = float(np.linalg.norm(actual@joint@actual.conj().T-joint, 2))
    if invariance > SOURCE_ADMISSION:
        raise ArithmeticError('full quantum Gibbs state is not fixed by source collision')
    return {'edge': [i, j], 'bath_gap': gap, 'thermal_excited_population': p1,
        'full_collision_input_columns': 32, 'purified_joint_system_input_columns': 16,
        'source_operator_error': error, 'energy_commutator_diagnostic': commutator,
        'literal_endpoint_total_energy_defect': endpoint_energy_defect,
        'full_quantum_Gibbs_stationarity_error': invariance,
        'purified_all_column_operator_error': purified_error,
        'native_pulses': compiler.resources(word)['native_pulses'],
        'bath_purification_native_pulses': compiler.resources(bath_word)['native_pulses'],
        'proper_code_service_events_required': 2*(compiler.resources(word)['native_pulses']+compiler.resources(bath_word)['native_pulses']),
        'thermal_bath_state_and_controls_are_supplied': True,
        'source_all_input_operator_admission': SOURCE_ADMISSION,
        'source_energy_commutator_admission': 100*SOURCE_ADMISSION,
        'physical_control_work_or_error_certified': False}, {'collision_unitary': actual, 'purified_collision_isometry': purified}


def source_energy_readout(model):
    """All16 input columns and all16 records of nondemolition energy readout."""
    word = compiler.compile_real_orthogonal(np.asarray(model.frame.real))
    read = compiler.inverse(word)+[event for bit in range(4) for event in compiler.cnot(bit, 4+bit)]+word
    incoming = np.zeros((256, 16), complex); incoming[:16] = np.eye(16)
    actual = compiler.execute(read, 8, incoming)
    target = np.concatenate([np.outer(model.frame[:, i], model.frame[:, i].conj()) for i in range(16)])
    error = float(np.linalg.norm(actual-target, 2))
    if error > SOURCE_ADMISSION:
        raise ArithmeticError('source energy readout loses an input column, record or reference')
    return {'full_joint_input_columns': 16, 'complete_energy_record_values': 16,
            'record_labels_refine_energy_by_declared_geometry_and_charge': True,
            'source_all_column_operator_error': error, 'native_pulses': compiler.resources(read)['native_pulses'],
            'new_committed_energy_record_bits': 4,
            'measurement_control_work_derived': False}, actual


def fermion_model():
    """Same law on the native-derived open CAR chain, in fixedN=1 only."""
    from ..common_source import repair_field_attachment as field
    hs = tuple(field.chain_hamiltonian(0, strength=r, cells=3) for r in GEOMETRIES)
    return finite_model(hs, beta=BETA, charge=field.occupation(3), sector=1)


def _source_proposal():
    states = []
    for bit in (0, 1):
        loads = [0]*12; loads[bit] = loads[5] = 1
        states.append(tuple(loads))
    columns = np.zeros((512, 2), complex)
    for bit, state in enumerate(states):
        columns[circuit.encode(state), bit] = 1
    arrays, errors = {}, {}
    for law in ('C', 'E'):
        branches = native_source.repair_branches(law, columns)
        outputs = np.stack([branch['amplitude'] for branch in branches])
        target = np.zeros_like(outputs)
        for branch, (coin, raw) in enumerate((coin, raw) for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law])):
            for initial, state in enumerate(states):
                final, (_, current) = native.native(state, (0, 1), coin)
                record = {0: 0, 1: 1, -1: 2}[current]
                if law == 'E':
                    record += 4*((state[0]+state[1]) % 2)
                if raw == record:
                    target[branch, circuit.encode(final), initial] = 1/math.sqrt(2)
        error = float(np.linalg.norm((outputs-target).reshape(-1, 2), 2))
        if error > SOURCE_ADMISSION:
            raise ArithmeticError('actual native proposal source differs from all-raw-record target')
        arrays['native_proposal_'+law] = outputs
        errors[law] = error
    return errors, arrays


def build():
    model = finite_model()
    theorem = exact_balance_identity()
    arrays = {}
    for law in ('C', 'E'):
        rows, weights, heat, native_records = records(model, law)
        arrays.update({law+'_record_labels': rows, law+'_record_weights': weights,
                       law+'_record_system_heat': heat, law+'_native_proposal_records': native_records})
    source, quantum = source_collision(model)
    arrays.update(quantum)
    energy_read, arrays['energy_readout_isometry'] = source_energy_readout(model)
    proposal_errors, proposals = _source_proposal(); arrays.update(proposals)
    channel_error = float(np.linalg.norm(apply_channel(model, model.gibbs)-model.gibbs, 2))
    if channel_error > 1e-12:
        raise ArithmeticError('complete quantum recorded heat bath is not stationary')
    commutator = model.exact_hamiltonians[0]*model.exact_hamiltonians[1]-model.exact_hamiltonians[1]*model.exact_hamiltonians[0]
    fermion = fermion_model()
    fermion_source, fermion_arrays = source_collision(fermion)
    arrays.update({'fermion_'+key: value for key, value in fermion_arrays.items()})
    for law in ('C', 'E'):
        rows, weights, heat, _ = records(fermion, law)
        arrays.update({'fermion_'+law+'_record_labels': rows, 'fermion_'+law+'_record_weights': weights,
                       'fermion_'+law+'_record_system_heat': heat})
    charge_energy_basis = np.diag(fermion.charges)
    charge_original = fermion.frame@charge_energy_basis@fermion.frame.conj().T
    U = fermion_arrays['collision_unitary']
    Q = np.kron(np.eye(2), charge_original)
    charge_error = float(np.linalg.norm(U@Q-Q@U, 2))
    if charge_error > SOURCE_ADMISSION:
        raise ArithmeticError('neutral heat bath changes the conserved encoded occupation')
    return {'schema': 'oph.reciprocal_vacuum.v1',
        'classification': 'finite reciprocal thermal geometry-field closure with supplied reservoirs',
        'independent_compatibility': 'joint Gibbs state exp(-beta*H_joint)/Z; fixed field Hamiltonians and geometry values before updates',
        'scalar_control': {'geometries': list(map(str, GEOMETRIES)), 'field_Fock_cutoff': 2,
            'geometry_only_potential': '0 on both branches, declared independently',
            'native_angular_band_eigenvalue': '6', 'field_oscillators': 3,
            'inverse_temperature_supplied': str(BETA), 'joint_quantum_dimension': 16,
            'spectral_frame_diagnostics': list(model.spectral_diagnostics),
            'conditional_field_Hamiltonians_do_not_commute_exactly': commutator != s.zeros(8),
            'full_recorded_quantum_channel_stationarity_error': channel_error,
            'thermodynamics': thermodynamics(model)},
        'exact_finite_theorem': theorem, 'native_proposal_source_errors': proposal_errors,
        'source_collision': source,
        'source_energy_readout': energy_read,
        'native_CAR_matter_control': {'field': 'repair_field_attachment.chain_hamiltonian(sector=0,strength=r,cells=3)',
            'declared_geometry_coupling': 'one shared r multiplies both native exchange bonds',
            'fixed_total_encoded_occupation': 1, 'admitted_joint_quantum_dimension': len(fermion.active),
            'same_supplied_inverse_temperature': str(BETA),
            'fixed_gap_bath_bank_edges': 9,
            'spectral_frame_diagnostics': list(fermion.spectral_diagnostics),
            'thermodynamics': thermodynamics(fermion), 'source_collision': fermion_source,
            'source_encoded_charge_commutator': charge_error,
            'bath_particle_exchange_or_chemical_potential_assumed': False,
            'unique_grandcanonical_state_from_arbitrary_charge_weights_claimed': False,
            'finite_temperature_one_particle_response_is_vacuum_gravity': False},
        'finite_bath_contract': {'one_fixed_gap_qubit_per_cross_geometry_energy_edge': 64,
            'bath_pair_orientation_including_degeneracy': 'fixed order by(energy,joint_basis_index), independent of transition direction',
            'additional_degenerate_noop_record_bit': 1,
            'all_bath_qubits_have_the_same_supplied_beta': True,
            'bath_Hamiltonians_retuned_during_update': False,
            'fresh_thermal_inputs_or_rethermalization_required_for_repeated_Markov_law': True,
            'all_rejections_bath_inputs_outputs_and_auxiliary_native_records_retained': True,
            'measurement_and_control_records_are_additional_memory_services': True,
            'bath_selection_and_rethermalization_are_not_free_or_derived': True,
            'compiled_endpoint_energy_conservation_is_not_per_pulse_or_autonomous_controller_energy_conservation': True,
            'proposal_geometry_is_auxiliary_until_energy_balanced_acceptance': True,
            'initial_energy_measurement_and_auxiliary_copy_of_classical_geometry_are_supplied': True,
            'stationary_external_drive_replaced_by_explicit_thermal_heat_exchange': True},
        'scope': {'isolated_autonomous_cooling_or_reservoir_preparation_derived': False,
            'energy_temperature_units': 'one declared model-coordinate generator for both controls; no physical clock, joule or kelvin identification',
            'stationarity_is_an_ensemble_density_while_record_histories_continue': True,
            'source_resource_count_scope': 'one energy-readout, nativeproposal primitive and one endpointcollision/bathpurification; full64edgecontroller andreservoirrefresh notsourceexecuted',
            'all_measurement_bath_selection_and_feedback_service_costs_included': False,
            'finite_temperature_control_is_an_exact_zero_temperature_vacuum': False,
            'scalar_D8_or_full36_oscillator_equilibrium_constructed_here': False,
            'bare_area_or_Einstein_stationarity_derived': False,
            'Hamiltonian_temperature_metric_and_physical_clock_dictionary_derived': False,
            'source_floating_diagnostics_are_outward_certificates': False}}, arrays


def source_bindings():
    from ..common_source import parameter_selection
    names = set(vacuum_repair.vacuum.SOURCE_PATHS)|set(compiler.source_summary()['source_sha256'])|{
        'model/native_equilibration/vacuum_repair.py', 'model/native_equilibration/reciprocal_vacuum.py',
        'model/native_equilibration/test_vacuum_repair.py',
        'model/native_equilibration/test_reciprocal_vacuum.py', 'model/observer/source_scalar/compiler.py',
        'model/observer/native.py', 'model/observer/source_discrimination/circuit.py',
        'model/observer/source_discrimination/experiment.py', 'tools/strict_json.py', 'requirements.txt',
        'model/common_source/repair_field_attachment.py', 'model/common_source/test_repair_field_attachment.py'}
    names.update(parameter_selection.source_bindings())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _compare(expected, actual, path='report'):
    if type(expected) is not type(actual):
        raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif type(expected) is float:
        exact = path.endswith(('source_all_input_operator_admission', 'source_energy_commutator_admission'))
        if not math.isfinite(actual) or (expected != actual if exact else not math.isclose(expected, actual, rel_tol=2e-9, abs_tol=2e-10)):
            raise ValueError(path+' numerical diagnostic differs')
    elif expected != actual:
        raise ValueError(path+' differs')


def check(path):
    path = Path(path); retained = strict_json.load(path)
    if type(retained) is not dict or retained.get('source_bindings') != source_bindings():
        raise ValueError('reciprocal finite source custody differs')
    artifact = path.with_suffix('.npz')
    if retained.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('reciprocal complete source artifact custody differs')
    fresh, arrays = build()
    _compare(fresh, {k: v for k, v in retained.items() if k not in ('source_bindings', 'artifact')})
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)) or set(archive.files) != set(arrays):
            raise ValueError('reciprocal complete artifact keys differ')
        for name, target in arrays.items():
            value = archive[name]
            if value.dtype != target.dtype or value.shape != target.shape or not np.isfinite(value).all():
                raise ValueError('reciprocal complete artifact schema differs: '+name)
            if target.dtype.kind in 'iu':
                if not np.array_equal(value, target):
                    raise ValueError('reciprocal native/bath records differ: '+name)
            elif 'record_weights' in name or 'record_system_heat' in name:
                if not np.allclose(value, target, rtol=2e-12, atol=2e-12):
                    raise ValueError('reciprocal classical heat probabilities differ: '+name)
            else:
                if 'native_proposal_' in name:
                    distance = float(np.linalg.norm((value-target).reshape(-1, 2), 2))
                    error = fresh['native_proposal_source_errors'][name[-1]]
                else:
                    distance = float(np.linalg.norm(value-target, 2))
                    case = fresh['native_CAR_matter_control']['source_collision'] if name.startswith('fermion_') else fresh['source_collision']
                    error = (fresh['source_energy_readout']['source_all_column_operator_error'] if name == 'energy_readout_isometry'
                             else case['purified_all_column_operator_error'] if name.endswith('purified_collision_isometry')
                             else case['source_operator_error'])
                if distance > REPLAY_TOLERANCE or distance+error > SOURCE_ADMISSION:
                    raise ValueError('retained source array exceeds unchanged all-input target admission: '+name)
        for prefix, model, contract in (('', finite_model(), fresh['source_collision']),
                ('fermion_', fermion_model(), fresh['native_CAR_matter_control']['source_collision'])):
            i, j = contract['edge']
            low, high = sorted((i, j), key=lambda state: (model.energies[state], state))
            gap = model.energies[high]-model.energies[low]
            total = np.kron(np.eye(2), model.hamiltonian)+np.kron(np.diag([0., gap]), np.eye(16))
            u = archive[prefix+'collision_unitary']
            if np.linalg.norm(u.conj().T@total@u-total, 2) > 100*SOURCE_ADMISSION:
                raise ValueError('retained collision fails unchanged endpoint energy admission')
            p = model.weights[high]/(model.weights[low]+model.weights[high])
            joint = np.kron(np.diag([1-p, p]), model.gibbs)
            if np.linalg.norm(u@joint@u.conj().T-joint, 2) > SOURCE_ADMISSION:
                raise ValueError('retained collision fails unchanged quantum Gibbs admission')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check); print('reciprocal finite geometry/field thermal collision replay: PASS'); return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('fresh receipt and complete source artifact required')
    bindings = source_bindings()
    report, arrays = build()
    if source_bindings() != bindings:
        raise RuntimeError('reciprocal executable sources changed during production')
    report['source_bindings'] = bindings
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')
    print('reciprocal joint geometry probabilities:', report['scalar_control']['thermodynamics']['geometry_probabilities'])


if __name__ == '__main__':
    main()
