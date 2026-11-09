"""Exact native-memory diffusion on the coherent family harmonics.

This is the stationary unread h*G1 coefficient at the fixed native rates.
It is not a physical self-energy or a replacement for the complete archive.
"""
from functools import lru_cache
import json

import sympy as s

from model.common_source import native_family_record_memory as memory

POLICIES = ('always', 'C', 'E', 'CE')
SOURCE_INPUTS = (
    'model/common_source/native_family_record_memory.py',
    'model/common_source/slow_geometry.py',
    'model/particle_memory_response/certificate.py',
    'model/particle_memory_response/PROOF.md',
)


def _clean(matrix):
    return s.ImmutableMatrix(s.Matrix(matrix).applyfunc(
        lambda x: s.factor(s.radsimp(s.expand(x)))))


def _comm(a, b):
    return a*b-b*a


def _entries(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def _first_nonzero(matrix):
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            if matrix[i, j] != 0:
                return {'row_zero_based': i, 'column_zero_based': j,
                        'value': str(matrix[i, j])}
    return None


@lru_cache(maxsize=1)
def _source():
    source = memory.architecture()
    if (source['p'], source['q']) != (s.Rational(1, 12), s.Rational(1, 6)):
        raise ArithmeticError('the fixed native event rates changed')
    geometry, (gc, ge) = memory.geometry()
    frame = geometry['frame']; norm = (frame.T*frame)[0, 0]
    def restrict(g):
        return _clean(frame.T*memory.objects.permutation_matrix(g)*frame/norm)
    r = _clean(restrict(ge)*restrict(gc))
    if _clean(r**3-s.eye(3)) != s.zeros(3) or _clean(r.T*r-s.eye(3)) != s.zeros(3):
        raise ArithmeticError('the actual family cycle changed')
    transition = s.ImmutableMatrix(source['T'])
    pi = s.ImmutableMatrix(source['stationary'])
    stationary = pi*s.ones(1, 12)
    inverse = _clean((s.eye(12)-transition+stationary).inv()-stationary)
    if _clean((s.eye(12)-transition)*inverse-s.eye(12)+stationary) != s.zeros(12):
        raise ArithmeticError('the full source group inverse failed')
    if len(source['records']) != 128:
        raise ArithmeticError('the complete native record count changed')
    return source, transition, pi, inverse, r


@lru_cache(maxsize=4)
def _response(duty):
    if duty not in POLICIES:
        raise ValueError('one of the four specified duties is required')
    source, transition, pi, inverse, r = _source()
    enabled = s.zeros(12); raw_transition = s.zeros(12)
    count = 0
    for old in range(12):
        for cc, cr, ec, er, _, weight in source['records']:
            c, e = cc == 0 and cr == 1, ec == 0 and er in (4, 5)
            dest = memory._destination(old, c, e, 2, source['gC'], source['gE'])
            raw_transition[dest, old] += weight
            if duty == 'always' or duty == 'C' and c or duty == 'E' and e or duty == 'CE' and c and e:
                enabled[dest, old] += weight
            count += 1
    if raw_transition != transition or count != 1536:
        raise ArithmeticError('complete raw records do not reconstruct the source')
    weight = enabled*pi
    k = _clean(s.diag(*weight)/2-weight*weight.T/2+enabled*inverse*s.diag(*weight))
    diffusion = _clean(k+k.T)
    projectors = tuple(s.ImmutableMatrix(p) for p in source['projectors'])
    return {'T': transition, 'pi': pi, 'inverse': inverse,
            'enabled': s.ImmutableMatrix(enabled), 'weight': _clean(weight),
            'K': k, 'D': diffusion, 'projectors': projectors, 'R': r,
            'raw_records_over_all_old_inputs': count}


def exact_response(duty='always'):
    """Fresh map of immutable exact matrices; all twelve inputs are retained."""
    return dict(_response(duty))


def diffusion(matrix, duty='always'):
    """Heisenberg diffusion coefficient on I or X tensor a family matrix."""
    matrix = s.Matrix(matrix)
    if matrix.shape != (3, 3) or any(x.has(s.Float) for x in matrix):
        raise ValueError('an exact three-by-three family matrix is required')
    row = _response(duty); p, d = row['projectors'], row['D']
    value = -sum((d[i, j]*_comm(p[i], _comm(p[j], matrix))
                  for i in range(12) for j in range(12) if d[i, j] != 0), s.zeros(3))/2
    return _clean(value)


@lru_cache(maxsize=1)
def _certificate_json():
    rows = []
    for duty in POLICIES:
        response = _response(duty); r = response['R']
        h, j = _clean(r+r*r), _clean(s.I*(r-r*r))
        dh, dj, dr = (diffusion(b, duty) for b in (h, j, r))
        lh, lj = _clean(_comm(dh, r)), _clean(_comm(dj, r))
        if _first_nonzero(lh) is None or _first_nonzero(lj) is None:
            raise ArithmeticError('a native harmonic unexpectedly closes')
        if _clean(dh.conjugate()-dh) != s.zeros(3) or _clean(dj.conjugate()+dj) != s.zeros(3):
            raise ArithmeticError('the source real/imaginary decomposition changed')
        if any(s.trace(b) != 0 for b in (dh, dj, dr)) or diffusion(s.eye(3), duty) != s.zeros(3):
            raise ArithmeticError('the diffusion lost trace preservation or unitality')
        mean = _clean(sum((response['weight'][i]*response['projectors'][i]
                           for i in range(12)), s.zeros(3)))
        rows.append({'duty': duty,
                     'original_duty_weight': str(sum(response['weight'])),
                     'all_old_inputs_times_raw_records': response['raw_records_over_all_old_inputs'],
                     'diffusion_H': _entries(dh), 'diffusion_J': _entries(dj),
                     'diffusion_R': _entries(dr),
                     'commutator_diffusion_H_R': _entries(lh),
                     'commutator_diffusion_J_R': _entries(lj),
                     'H_leak_witness': _first_nonzero(lh),
                     'J_leak_witness': _first_nonzero(lj),
                     'mean_H_commutator': _entries(_clean(_comm(mean, h))),
                     'scalar_C3_damping_mode': False})
    result = {'fixed_native_rates': {'p': '1/12', 'q': '1/6', 'gain': '1'},
              'native_memory_dimension': 12, 'native_raw_records': 128,
              'native_R': _entries(_source()[-1]), 'duties': rows,
              'complete_original_records_reconstructed': True,
              'success_conditioning_used': False,
              'all_lags_summed_by_exact_group_inverse': True,
              'scope': 'stationary unread h*G1 coefficient in the declared pulse clock',
              'positive_chamber_literal_scalar_damping_excluded': True,
              'different_rates_or_controller_excluded': False,
              'isospectral_match_after_family_rotation_excluded': False,
              'isolated_time_coincidence_excluded': False,
              'chi_response_rate_derived': False,
              'physical_selfenergy_derived': False,
              'physical_mass_assigned': None}
    return json.dumps(result, sort_keys=True, allow_nan=False)


def certificate():
    """JSON-safe detached exact source result, without a parent producer replay."""
    return json.loads(_certificate_json())
