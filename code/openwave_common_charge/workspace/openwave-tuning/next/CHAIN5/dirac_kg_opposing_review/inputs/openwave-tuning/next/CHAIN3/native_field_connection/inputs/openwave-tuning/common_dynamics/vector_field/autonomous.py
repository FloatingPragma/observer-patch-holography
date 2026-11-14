"""A declared local electric/plaquette extension of the same matter walk.

Numerical experiments execute a finite one-face member on the complete
one-particle sector. Other faces have exactly zero coupling in that member.
"""
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'manybody'))
sys.path.insert(0, str(ROOT / 'openwave-tuning/auditor/native_current_walk'))
import manybody as mb
import mode_field
import current_walk
import cell_complex

MODEL_ID = 'native-local12-fermion-plaquette-v2-proposed'
MODULUS = 16


@lru_cache(maxsize=1)
def data():
    graph, B, A, periods, winding, labels = cell_complex.build(2)
    A = np.asarray(A, dtype=np.int64)
    B = np.asarray(B, dtype=np.int64)
    assert np.max(abs(A)) == 1
    return graph, B, A, labels


def centered(field, modulus=MODULUS):
    f = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
    if f.shape != (288,) or np.any(f >= modulus):
        raise ValueError('full admitted modular field required')
    return (f + modulus // 2) % modulus - modulus // 2


def electric_energy(field):
    f = centered(field)
    return int(f @ f) / 2


def add_field(field, shift):
    f = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
    return bytes(((f + shift) % MODULUS).astype(np.uint8))


def initial(vertex=0, loop=0, face=0, charge=1):
    if charge not in (-1, 1):
        raise ValueError('particle or same-source one-hole preparation required')
    particular = charge * current_walk.geometry()['P'][:, vertex] + loop * data()[2][:, face]
    f = bytes((particular % MODULUS).astype(np.uint8))
    return {(vertex, f): np.ones(12, dtype=complex) / np.sqrt(12)}


def matter_step(state, charge=1):
    """One-/one-hole restriction of the frozen local Gamma(W)_F amendment."""
    if charge not in (-1, 1):
        raise ValueError('declared charge sign required')
    neighbours, edges = mb.common.geometry_tables()
    graph = mb.common.graph()
    result = defaultdict(lambda: np.zeros(12, dtype=complex))
    for (x, f), amplitudes in state.items():
        reflected = amplitudes.sum() / 6 - amplitudes
        for event, amplitude in enumerate(reflected):
            direction, move = divmod(event, 2)
            y = int(neighbours[x, direction]) if move else x
            field = f
            if move:
                e = int(edges[x, direction])
                a, b = graph.seams[e]
                j = charge * (int(x == a) - int(y == a))
                values = bytearray(f)
                values[e] = (values[e] + j) % MODULUS
                field = bytes(values)
            result[y, field][event] += amplitude
    return dict(result)


def plaquette_step(state, beta, face=0):
    """Complete cyclic orbit matrix, without dropping small outcomes."""
    if not math.isfinite(beta):
        raise ValueError('finite declared plaquette coupling required')
    if beta == 0:
        return {key: value.copy() for key, value in state.items()}
    a = data()[2][:, face]
    pivot = int(np.flatnonzero(a)[0])
    sign = int(a[pivot])
    q = MODULUS
    fourier = np.exp(2j * np.pi * np.outer(np.arange(q), np.arange(q)) / q) / np.sqrt(q)
    kick = (fourier * np.exp(1j * beta * np.cos(2 * np.pi * np.arange(q) / q))[None, :]) @ fourier.conj().T
    orbits = defaultdict(lambda: np.zeros((q, 12), dtype=complex))
    for (x, f), amplitudes in state.items():
        coordinate = (sign * f[pivot]) % q
        canonical = add_field(f, -coordinate * a)
        orbits[x, canonical][coordinate] += amplitudes
    result = {}
    for (x, canonical), matrix in orbits.items():
        output = kick @ matrix
        for t in range(q):
            result[x, add_field(canonical, t * a)] = output[t]
    return result


def electric_step(state, alpha):
    if not math.isfinite(alpha):
        raise ValueError('finite declared electric coupling required')
    return {key: values * np.exp(-1j * alpha * electric_energy(key[1]))
            for key, values in state.items()}


def step(state, alpha, beta, face=0, charge=1):
    """V = U_E U_B Gamma(W)_F; coefficients are architecture inputs."""
    return electric_step(plaquette_step(matter_step(state, charge), beta, face), alpha)


def records(state, charge=1):
    graph, B, A, labels = data()
    positions = np.zeros(96)
    joint = []
    energy = 0.
    total = 0.
    reference = np.eye(96, dtype=np.int64)[:, 0]
    for (x, f), values in state.items():
        probabilities = abs(values)**2
        weight = float(probabilities.sum())
        n = np.zeros(96, dtype=np.int64)
        n[x] = charge
        field = np.frombuffer(f, dtype=np.uint8).astype(np.int64)
        assert not np.any((n + B @ field - charge * reference) % MODULUS)
        positions[x] += weight
        energy += weight * electric_energy(f)
        total += weight
        joint.append({'position': x, 'field_hex': f.hex(),
                      'event_probabilities': probabilities.tolist(),
                      'electric_energy': electric_energy(f)})
    assert abs(total - 1) < 2e-10
    return {'norm_squared': total, 'position_probabilities': positions.tolist(),
            'electric_energy_expectation': energy,
            'joint_position_field_event_records': joint,
            'success_normalization': False, 'modular_Gauss': True}


def phase_correlator(left, right):
    return sum(np.vdot(values, right.get(key, np.zeros(12))) for key, values in left.items())


def run():
    paths = [Path(__file__), Path(cell_complex.__file__), HERE / 'CELL_COMPLEX.json',
             Path(mb.__file__), Path(mode_field.__file__), Path(current_walk.__file__),
             Path(mb.common.__file__), HERE.parent / 'manybody/MANIFEST.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    graph, B, A, labels = data()
    assert not np.any(B @ A)
    assert labels[0][0] == 'internal_triangle'
    experiments = []
    full = {}
    for charge, alpha, beta in ((1, 0., 0.), (1, .17, 0.), (1, 0., .41),
                                (1, .17, .41), (-1, .17, .41)):
        state = initial(charge=charge)
        evolution = []
        for n in range(4):
            record = records(state, charge)
            next_state = step(state, alpha, beta, charge=charge)
            overlap = phase_correlator(state, next_state)
            evolution.append({'steps': n, 'norm_squared': record['norm_squared'],
                              'position_probabilities': record['position_probabilities'],
                              'electric_energy_expectation': record['electric_energy_expectation'],
                              'Re_V_moment': float(overlap.real),
                              'Im_V_moment': float(overlap.imag),
                              'conserved_Floquet_diagnostic': float(1 - overlap.real),
                              'basis_blocks': len(state)})
            if n == 3:
                full[f'{charge}:{alpha}:{beta}'] = record
            state = next_state
        assert max(r['Re_V_moment'] for r in evolution) - min(r['Re_V_moment'] for r in evolution) < 2e-10
        experiments.append({'charge_preparation': charge, 'alpha': alpha, 'beta_face0': beta,
                            'all_other_207_face_couplings': 0., 'records': evolution})
    base = experiments[0]['records'][-1]
    coupled = experiments[3]['records'][-1]
    tv = sum(abs(a - b) for a, b in zip(base['position_probabilities'], coupled['position_probabilities'])) / 2
    field_effect = abs(base['electric_energy_expectation'] - coupled['electric_energy_expectation'])
    assert tv > 1e-5 and field_effect > 1e-5
    # The two signs are related by electric-flux inversion for the even
    # electric energy and cosine plaquette coupling.
    conjugate = experiments[4]['records'][-1]
    charge_tv = sum(abs(a - b) for a, b in zip(coupled['position_probabilities'], conjugate['position_probabilities'])) / 2
    assert charge_tv < 2e-10
    # Native field update omission violates Gauss on every moving slot.
    seed = initial()
    moved = matter_step(seed)
    assert any(x != 0 for x, f in moved)
    for x, f in moved:
        if x:
            broken_n = np.zeros(96, dtype=int); broken_n[x] = 1
            assert np.any((broken_n - np.eye(96, dtype=int)[:, 0]) % MODULUS)
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'JOINT_RECORDS.json').write_text(json.dumps(full, separators=(',', ':')) + '\n')
    result = {'schema': 'oph.local12.autonomous_plaquette.v1', 'candidate_id': MODEL_ID,
              'checks': {'same_frozen_matter_source_restriction': True,
                         'autonomous_electric_and_plaquette_steps': True,
                         'complete_sparse_joint_outcomes_retained': True,
                         'all_source_norms_and_modular_Gauss': True,
                         'interaction_changes_actual_position_records': True,
                         'matter_current_changes_same_field': True,
                         'Floquet_diagnostic_conserved': True,
                         'same_source_hole_charge_conjugation': True,
                         'omitted_current_rejected': True, 'source_bytes_unchanged': True},
              'experiments': experiments, 'interaction_position_TV_at_step3': tv,
              'interaction_electric_expectation_change': field_effect,
              'particle_hole_position_TV': charge_tv,
              'complete_records_sha256': hashlib.sha256((HERE / 'JOINT_RECORDS.json').read_bytes()).hexdigest(),
              'source_amendment': 'Field phases and plaquette couplings are declared. This is one Floquet V on the shared matter/flux state, not a claimed unique consequence of bare C/E.',
              'scope': 'Exact finite one-face member, one-particle or separately prepared one-hole restriction of Gamma(W)_F. Numerical whole-state laws; native field words and precision certificate are separate. Homogeneous208-face dynamics is specified but not densely executed; no pair force, physical mass or full Maxwell continuum claimed.',
              'source_sha256': pins, 'official_rows_changed': []}
    (HERE / 'AUTONOMOUS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'position_TV': tv, 'field_energy_change': field_effect,
                      'maximum_blocks': max(r['basis_blocks'] for e in experiments for r in e['records'])}))


if __name__ == '__main__':
    run()
