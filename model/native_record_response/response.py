"""Complete native record response on the fixed weighted H66 source.

Only the source's exported C/E records are read. The matrices describe the
actual discrete repair process, with its readout backaction and all outcomes.
No coherent source control, spatial coordinate map or field action is added.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from model.local_repair_family import local_bias as local, mixing
from model.native_program import interpreter
from tools import strict_json

DIMENSION = 66
PRIME = 65521
ROOT = Path(interpreter.__file__).resolve().parents[2]
HERE = ROOT/'model/native_record_response'
OUTPUT = HERE/'results.json'
SOURCE_FILES = (
    'model/native_record_response/__init__.py',
    'model/native_record_response/response.py',
    'model/native_record_response/PROOF.md',
    'model/native_record_response/test_response.py',
    'model/native_program/interpreter.py',
    'model/observer/source_discrimination/circuit.py',
    'model/observer/source_scalar/compiler.py',
    'model/requirements.txt',
    'tools/strict_json.py',
)


@dataclass(frozen=True)
class Branch:
    kind: str
    law: str | None
    seam: tuple[int, int] | None
    coin: int | None
    raw: int | None
    label: int | None
    current: int | None
    weight: F
    pairs: tuple[tuple[int, int], ...]

    def feature(self, seam):
        """Unconditional score; no branch is discarded or renormalized."""
        return (60*self.label if self.law == 'E' and self.seam == seam
                and self.label is not None and self.pairs else 0)


def source():
    value = interpreter.program()
    expected = {'port_weights': [25, 25]+[8]*10, 'seam_weights': [1]*30,
                'c_weight': '1/2', 'repair_probability': '1', 'exploration_probability': '1/2'}
    command = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1}
    if (value['carrier'] != {'kind': 'local', 'shell': 'H66'}
            or value['parameters'] != expected or value['commands'] != [command]
            or value['memory'] != {'capacity': 1, 'overflow': 'export_oldest'}):
        raise ArithmeticError('the authoritative native member changed')
    return value


def parameters():
    p = source()['parameters']
    return local.Parameters(tuple(p['port_weights']), tuple(p['seam_weights']),
                            F(p['c_weight']), F(p['repair_probability']))


def branches():
    """Read all 721 interpreter outcomes using a classical reference marker.

    The operator sum_i |i,i><i,i| tags every input, without postselection.
    The interpreter is linear, so normalization is immaterial for this audit.
    Its positive partial injections also specify the full quantum branch maps.
    """
    program, p = source(), parameters()
    marker = {(i, i, i, i): 1 for i in range(DIMENSION)}
    outcomes = interpreter._outcomes(program, program['commands'][0], marker)
    result = []
    for outcome in outcomes:
        r = outcome.record
        weight = F() if r.kind == 'hold' else F(1, 60)*local.coin_probability(p, r.seam, r.coin)
        pairs = []
        for (j, i, k, other), coefficient in outcome.density.items():
            if j != k or i != other or coefficient != weight:
                raise ArithmeticError('a native marker output changed its meaning')
            pairs.append((i, j))
        if len({i for i, _ in pairs}) != len(pairs) or len({j for _, j in pairs}) != len(pairs):
            raise ArithmeticError('a native branch ceased to be a partial injection')
        result.append(Branch(r.kind, r.law, r.seam, r.coin, r.raw, r.label,
                             r.current, weight, tuple(sorted(pairs))))
    if len(result) != 721:
        raise ArithmeticError('complete raw record alphabet changed')
    return tuple(result)


def data():
    p, rows, states = parameters(), branches(), local.basis()
    pi = s.ImmutableMatrix(local.reference(p))
    P = s.zeros(DIMENSION)
    B = [s.zeros(DIMENSION) for _ in local.SEAMS]
    for row in rows:
        for i, j in row.pairs:
            P[j, i] += row.weight
            if row.law == 'E':
                B[local.SEAMS.index(row.seam)][j, i] += row.weight*row.feature(row.seam)
    D = s.ImmutableMatrix([[x[a]-x[b] for x in states] for a, b in local.SEAMS])
    G = s.ImmutableMatrix([[F(p.port_weights[a]-p.port_weights[b],
                             p.port_weights[a]+p.port_weights[b])*int(x[a]+x[b] == 1)
                           for x in states] for a, b in local.SEAMS])
    P = s.ImmutableMatrix(P)
    B = tuple(map(s.ImmutableMatrix, B))
    return {'parameters': p, 'branches': rows, 'states': states, 'pi': pi,
            'P': P, 'B': B, 'D': D, 'G': G, 'mean': D*pi}


def validate(d):
    P, pi, D, G, B = (d[key] for key in ('P', 'pi', 'D', 'G', 'B'))
    one = s.ones(1, DIMENSION)
    expected = {('hold', None, None, None, None)} | {
        ('repair', law, seam, coin, raw) for law in ('C', 'E') for seam in local.SEAMS
        for coin in (0, 1) for raw in range(4 if law == 'C' else 8)}
    labels = [(row.kind, row.law, row.seam, row.coin, row.raw) for row in d['branches']]
    if len(labels) != len(expected) or set(labels) != expected:
        raise ArithmeticError('complete raw record slots were dropped or duplicated')
    for row in d['branches']:
        if row.kind == 'repair' and row.label != interpreter.circuit.record_values(row.law, row.coin).get(row.raw):
            raise ArithmeticError('raw native record decoding changed')
    if one*P != one or P*pi != pi:
        raise ArithmeticError('native probability normalization or reference failed')
    if P*s.diag(*pi) != s.diag(*pi)*P.T:
        raise ArithmeticError('native detailed balance failed')
    if P != s.Matrix(local.kernel(d['parameters'])).T:
        raise ArithmeticError('independent endpoint population law disagrees')
    for e in range(30):
        if one*B[e] != D[e, :] or B[e]*pi != s.diag(*pi)*G[e, :].T:
            raise ArithmeticError('native readout backaction identity failed')
    for row in d['branches']:
        for i, j in row.pairs:
            a, b = row.seam
            if row.law == 'E' and row.label != d['states'][i][a]-d['states'][i][b]:
                raise ArithmeticError('the E record is not its incoming contrast')
            if any(d['states'][j][k]-d['states'][i][k] != row.current*(int(k == b)-int(k == a))
                   for k in range(12)):
                raise ArithmeticError('complete native current continuity failed')
    return True


def apply_branch(row, density):
    """Exact sparse source/reference output of one retained native branch.

    This is a linear operator API; positivity and normalization are properties
    of the caller's input. Rational and complex-linear inputs use the same map.
    """
    mapping = dict(row.pairs)
    return {(mapping[i], a, mapping[j], b): row.weight*value
            for (i, a, j, b), value in density.items() if i in mapping and j in mapping and value}


def history_output(rows, density):
    """Unnormalized source/reference output; every raw path has an outcome."""
    output = dict(density)
    for row in rows:
        output = apply_branch(row, output)
    return output


def _count(value, name='attempt count'):
    if type(value) is not int or value < 0:
        raise ValueError('a nonnegative integer '+name+' is required')
    return value


def _edge(value):
    if type(value) is not int or not 0 <= value < 30:
        raise ValueError('an integer seam index from zero through 29 is required')
    return value


def covariance(d, lag):
    """Exact 30 by 30 matrix E[Y_e(t+lag)Y_f(t)] minus mean products."""
    _count(lag)
    if lag == 0:
        second = s.diag(*(60*(d['D'][e, :].applyfunc(lambda x: x*x)*d['pi'])[0]
                          for e in range(30)))
    else:
        second = d['D']*(d['P']**(lag-1))*s.diag(*d['pi'])*d['G'].T
    return s.ImmutableMatrix(second-d['mean']*d['mean'].T)


def ordered_moment(d, events):
    """Complete raw-record moment at strictly increasing integer attempts.

    Each event is (attempt, seam_index). Omitted attempts execute the same
    unread native channel; the measured attempts retain their actual B maps.
    """
    if type(events) is not tuple or not events:
        raise ValueError('a nonempty tuple of (attempt, seam index) is required')
    if any(type(row) is not tuple or len(row) != 2 for row in events):
        raise ValueError('integer attempt and seam-index pairs are required')
    for t, edge in events:
        _count(t); _edge(edge)
    if any(events[i][0] >= events[i+1][0] for i in range(len(events)-1)):
        raise ValueError('event attempts must be strictly increasing')
    vector, previous = d['pi'], events[0][0]-1
    for t, edge in events:
        vector = d['B'][edge]*(d['P']**(t-previous-1))*vector
        previous = t
    return sum(vector)


def resolvent(d, z):
    """Exact centered record response sum_{k>=1} z^-k C(k), for rational z>1."""
    if type(z) is not F or z <= 1:
        raise ValueError('an exact rational resolvent parameter greater than one is required')
    A = s.diag(*d['pi'])*d['G'].T-d['pi']*d['mean'].T
    return s.ImmutableMatrix(d['D']*(z*s.eye(DIMENSION)-d['P']).inv()*A)


def _modular(matrix):
    return np.array([[int(v.p)*pow(int(v.q), -1, PRIME) % PRIME for v in row]
                     for row in matrix.tolist()], dtype=np.int64)


def independent_rows_modular(matrix):
    """Exact Gaussian elimination in F_65521; no numerical rank tolerance."""
    matrix = np.asarray(matrix, dtype=np.int64)
    if matrix.ndim != 2 or matrix.shape[1] != DIMENSION:
        raise ValueError('a matrix with 66 columns is required')
    pivots, selected = {}, []
    for index, row in enumerate(matrix):
        row = row.copy() % PRIME
        for column, pivot in sorted(pivots.items()):
            row = (row-row[column]*pivot) % PRIME
        nonzero = np.flatnonzero(row)
        if len(nonzero):
            column = int(nonzero[0])
            row = row*pow(int(row[column]), -1, PRIME) % PRIME
            pivots[column] = row
            selected.append(index)
    return tuple(selected)


def response_ranks(d, depth=12):
    """Exact full state observability and centered reachability certificates.

    Nonzero minors modulo 65521 prove the corresponding rational lower bounds.
    Ambient dimension and exact zero column sums give the matching upper bounds.
    """
    _count(depth, 'response depth')
    P, D = _modular(d['P']), _modular(d['D'])
    A = s.diag(*d['pi'])*d['G'].T-d['pi']*d['mean'].T
    if s.ones(1, DIMENSION)*A != s.zeros(1, 30):
        raise ArithmeticError('centered input has a stationary component')
    O, R = np.ones((1, DIMENSION), dtype=np.int64), np.empty((0, DIMENSION), dtype=np.int64)
    left, right = D, _modular(A)
    observations, reachability = [], []
    for k in range(depth+1):
        O, R = np.concatenate((O, left)), np.concatenate((R, right.T))
        observations.append(len(independent_rows_modular(O)))
        reachability.append(len(independent_rows_modular(R)))
        left, right = left@P % PRIME, P@right % PRIME
    return {'prime': PRIME, 'depth': depth, 'observable_rank_by_depth': observations,
            'reachable_rank_by_depth': reachability,
            'observation_minor_rows': list(independent_rows_modular(O)),
            'reachability_minor_columns': list(independent_rows_modular(R)),
            'full_population_dimension': DIMENSION, 'centered_dimension': DIMENSION-1}


def hierarchy(d):
    """Actual biased occupation drift and a four-configuration nonclosure witness."""
    states, p = d['states'], d['parameters']
    N = s.Matrix([[x[k] for x in states] for k in range(12)])
    predicted = N.copy()
    for e, (a, b) in enumerate(local.SEAMS):
        drift = (d['D'][e, :]-d['G'][e, :])/60
        predicted[a, :] -= drift
        predicted[b, :] += drift
    if N*d['P'] != predicted:
        raise ArithmeticError('biased occupation hierarchy failed')
    # An affine single-occupation function has vanishing rectangle contrast.
    index = {tuple(i for i, value in enumerate(x) if value): n for n, x in enumerate(states)}
    for port in range(12):
        row = N*d['P']-N
        for a in range(12):
            for b in range(a+1, 12):
                for c in range(b+1, 12):
                    for e in range(c+1, 12):
                        pairs = ((a, b), (c, e), (a, c), (b, e))
                        contrast = sum(sign*row[port, index[pair]]
                                       for sign, pair in zip((1, 1, -1, -1), pairs, strict=True))
                        if contrast:
                            return {'observable': 'occupation of port '+str(port),
                                    'pairs': [list(pair) for pair in pairs],
                                    'signs': [1, 1, -1, -1], 'drift_contrast': str(contrast),
                                    'mean_dimension_before': N.rank(),
                                    'mean_dimension_after_one_attempt': N.col_join(N*d['P']).rank()}
    raise ArithmeticError('the weighted nonclosure witness vanished')


def spectral_certificate(d):
    """Exact polynomial of the finite response, with no eigenvalue fitting."""
    C = 1980*d['P']
    if any(value.q != 1 for value in C):
        raise ArithmeticError('the native transition denominator changed')
    polynomial = C.charpoly().as_poly()
    factors = s.factor_list(polynomial)[1]
    if s.gcd(polynomial, polynomial.diff()).degree() != 0 or polynomial.eval(1980) != 0:
        raise ArithmeticError('the exact relaxation spectrum ceased to be simple')
    if sorted((factor.degree(), exponent) for factor, exponent in factors) != [
            (1, 1), (13, 1), (15, 1), (15, 1), (22, 1)]:
        raise ArithmeticError('the native spectral factorization changed')
    return {'integer_scale': 1980,
            'integer_transition_characteristic_coefficients': [str(v) for v in polynomial.all_coeffs()],
            'irreducible_factors': [{'degree': factor.degree(), 'multiplicity': exponent,
                                    'coefficients': [str(v) for v in factor.all_coeffs()]}
                                   for factor, exponent in factors],
            'squarefree': True, 'stationary_root': 1980,
            'distinct_centered_record_response_poles': 65,
            'positive_real_poles_below_one': True,
            'pole_units': 'one native attempt',
            'physical_mass_interpretation': False,
            'residues': 'nonzero real rank-one ordered-response matrices; no positive quantum spectral measure asserted'}


def preparation_error(attempts=1600):
    _count(attempts)
    summary = mixing.mixing_bound_summary(parameters(), attempts)
    return summary


def certificate():
    from . import transfer, short_records
    d = data(); validate(d)
    ranks = response_ranks(d)
    if ranks['observable_rank_by_depth'][-1] != 66 or ranks['reachable_rank_by_depth'][-1] != 65:
        raise ArithmeticError('the complete finite response lost observability or reachability')
    equal = [e for e, (a, b) in enumerate(local.SEAMS)
             if d['parameters'].port_weights[a] == d['parameters'].port_weights[b]]
    if any(d['G'][e, :] != s.zeros(1, DIMENSION) or d['mean'][e] != 0 for e in equal):
        raise ArithmeticError('equal-weight readout erasure failed')
    stationary = Counter(str(v) for v in d['pi'])
    return {'schema': 'oph.native_record_response.v1',
            'source': source(),
            'full_raw_slots': len(d['branches']),
            'zero_slots_including_global_hold': sum(not row.pairs for row in d['branches']),
            'abort_slots': sum(row.kind == 'repair' and row.label is None for row in d['branches']),
            'reference': dict(sorted(stationary.items())),
            'record_features': 'Y_e = 60*1[law=E,seam=e]*incoming_contrast',
            'same_process_backaction': 'B_e*pi = diag(pi)*G_e; G_e=b_e*1[one endpoint occupied]',
            'zero_future_correlation_seams': [list(local.SEAMS[e]) for e in equal],
            'equal_weight_seam_count': len(equal),
            'contrast_response_dimensions': [30, 30],
            'minimum_native_hold_probability': str(min(d['P'][i, i] for i in range(66))),
            'rational_spectral_lower_bound': str(min(2*d['P'][i, i]-1 for i in range(66))),
            'pole_scope': 'real positive native-attempt relaxation factors; centered response excludes the stationary pole',
            'response_ranks': ranks,
            'complete_relaxation_spectrum': spectral_certificate(d),
            'positive_record_transfer': transfer.certificate(d, ranks),
            'short_record_population_realization': short_records.certificate(),
            'occupation_hierarchy': hierarchy(d),
            'preparation': preparation_error(),
            'same_time_variance_seam_0': str(covariance(d, 0)[0, 0]),
            'next_time_covariance_seam_0': str(covariance(d, 1)[0, 0]),
            'cross_response_lag_1_0_1': str(covariance(d, 1)[0, 1]),
            'physical_premises_discharged': [],
            'spatial_fourier_map_supplied': False,
            'coherent_hamiltonian_or_mass_supplied': False,
            'native_program_modified': False}


def source_bindings():
    from .transfer import SOURCE_FILES as transfer_sources
    from .short_records import SOURCE_FILES as short_record_sources
    paths = set(mixing.source_bindings()) | set(SOURCE_FILES) | set(transfer_sources) | set(short_record_sources)
    if any(path.startswith('docs/') for path in paths):
        raise ArithmeticError('research notes cannot be direct source evidence')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def report():
    bindings = source_bindings()
    value = certificate()
    if source_bindings() != bindings:
        raise ArithmeticError('native response sources changed during calculation')
    return {**value, 'source_sha256': bindings}


def check(path=OUTPUT):
    retained = strict_json.load(Path(path))
    actual = report()
    if json.dumps(retained, sort_keys=True, separators=(',', ':')) != json.dumps(actual, sort_keys=True, separators=(',', ':')):
        raise ValueError('native record response receipt or source bindings differ')
    return actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', nargs='?', const=str(OUTPUT))
    group.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Complete native record response: PASS')
    else:
        value = json.dumps(report(), indent=2, sort_keys=True)+'\n'
        if args.write:
            Path(args.write).write_text(value)
            print('Complete native record response receipt written')
        else:
            print(value, end='')


if __name__ == '__main__':
    main()
