"""Integer spectral exhaustion and independent composite energy sum rules."""
from fractions import Fraction as F

import numpy as np

from ..centered_particle_spectrum import certificate as inherited
from ..centered_particle_spectrum import packet as ground_packet
from . import source

sqrt_bounds = inherited.sqrt_bounds
decimal_bounds = inherited.decimal_bounds
VECTOR_BITS = inherited.VECTOR_BITS
ENERGY_BITS = inherited.ENERGY_BITS


def eigensystem(witness):
    if set(witness) != {'columns', 'energies', 'vectors'}:
        raise ValueError('complete composite witness arrays required')
    columns = np.asarray(witness['columns'])
    if columns.ndim != 1 or columns.dtype != np.dtype(np.int64):
        raise ValueError('exact int64 projector columns required')
    block = source.block(tuple(int(i) for i in columns))
    n = block['dimension']
    ei = inherited._integer_array(witness['energies'], (n,), 'energy numerators')
    vi = inherited._integer_array(witness['vectors'], (n, n), 'vector numerators')
    energies = tuple(F(int(e), 1 << ENERGY_BITS) for e in ei)
    if any(a >= b for a, b in zip(energies, energies[1:])):
        raise ValueError('strictly ordered spectral proposals required')
    y = block['basis'].astype(object) @ vi
    ty = np.zeros(y.shape, dtype=object)
    for i, row in enumerate(block['coefficient_source']):
        for j in np.flatnonzero(row):
            ty[i] += int(row[j]) * y[j]
    residual = ty * (1 << ENERGY_BITS) - 16 * y * ei[None, :]
    gram = y.T @ (block['sizes'][:, None] * y)
    norms = [F(int(gram[i, i]), 1 << (2*VECTOR_BITS)) for i in range(n)]
    if any(v <= 0 for v in norms):
        raise ArithmeticError('zero witness vector')
    defect = gram - (1 << (2*VECTOR_BITS)) * np.eye(n, dtype=object)
    delta = F(max(sum(abs(int(v)) for v in row) for row in defect),
              1 << (2*VECTOR_BITS))
    if not delta < F(1, 1000000):
        raise ArithmeticError('complete Gram defect outside admitted domain')
    squares = [F(sum(int(d)*int(v)**2 for d, v in
                     zip(block['sizes'], col, strict=True)),
                 1 << (2*(VECTOR_BITS+ENERGY_BITS+4))) for col in residual.T]
    radius = sqrt_bounds(sum(squares, F(0))/(1-delta))[1]
    if not radius < F(1, 10**9):
        raise ArithmeticError('complete spectral radius outside admitted domain')
    errors = []
    intervals = []
    for i, e in enumerate(energies):
        separation = min(abs(e-other) for j, other in enumerate(energies) if j != i) - radius
        r = sqrt_bounds(squares[i]/norms[i])[1]
        if separation <= 2*radius or r >= separation:
            raise ArithmeticError('composite spectral line not isolated')
        errors.append(2*r/separation)
        intervals.append((e-r, e+r))
    return dict(block=block, vectors=vi, norms=norms, energies=energies,
                vector_errors=errors, energy_intervals=intervals,
                gram_defect=delta, complete_radius=radius)


def sparse_integer_action(matrix, vector):
    """Independent exact source action without a reduced eigensystem."""
    matrix = matrix.tocsr()
    return np.array([sum(int(matrix.data[j])*int(vector[matrix.indices[j]])
                         for j in range(matrix.indptr[i], matrix.indptr[i+1]))
                     for i in range(matrix.shape[0])], dtype=object)


def source_moments(ground):
    """Enclose <A*A> and <A*[K,A]> directly on all source configurations."""
    g = ground['block']['full_basis'].astype(object) @ ground['vector_numerators'][:, 0]
    a = source.probe()
    q = sparse_integer_action(a, g)
    kg = source.parent.centered(4, 4)
    kc = source.parent.centered(5, 2)
    hq = sparse_integer_action(kc, q)
    akg = sparse_integer_action(a, sparse_integer_action(kg, g))
    norm = sum(int(v)**2 for v in g)
    weight = F(sum(int(v)**2 for v in q), norm)
    first = F(sum(int(v)*(int(h)-int(z)) for v, h, z in zip(q, hq, akg, strict=True)), 16*norm)
    # The certified normalized ground vector differs by epsilon.  For a
    # bounded operator B, its expectation changes by at most 2||B|| epsilon.
    epsilon = ground['vector_errors'][0]
    kg_bound = F(int(np.asarray(abs(kg).sum(axis=1)).max()), 16)
    kc_bound = F(int(np.asarray(abs(kc).sum(axis=1)).max()), 16)
    weight_error = 18*epsilon
    first_error = 18*(kg_bound+kc_bound)*epsilon
    return dict(weight=(max(F(0), weight-weight_error), weight+weight_error),
                first=(first-first_error, first+first_error),
                source_row_sum_norm_bounds=(kg_bound, kc_bound),
                source_rows=(len(g), len(q)),
                independent_of_composite_eigenvectors=True)


def intersect(left, right):
    result = max(left[0], right[0]), min(left[1], right[1])
    if result[0] > result[1]:
        raise ArithmeticError('spectral and direct source sum rules disagree')
    return result


def verify(witness):
    previous = ground_packet.load_witness()
    ground = inherited.eigensystem('ground', previous['ground_columns'],
                                  previous['ground_energies'], previous['ground_vectors'])
    local = eigensystem(witness)
    selection = source.selection_certificate(ground['block'], local['block'])
    transition = source.transition(ground['block'], local['block']).astype(object)
    amplitudes = local['vectors'].T @ (transition @ ground['vector_numerators'][:, 0])
    rows = []
    weights = []
    gaps = []
    for i, num in enumerate(amplitudes):
        rootlo, roothi = sqrt_bounds(ground['norms'][0]*local['norms'][i])
        if rootlo <= 0:
            raise ArithmeticError('normalization interval touches zero')
        raw = F(abs(int(num)), 1 << (2*VECTOR_BITS))
        error = 3*(ground['vector_errors'][0]+local['vector_errors'][i])
        lo, hi = max(F(0), raw/roothi-error)**2, (raw/rootlo+error)**2
        elo = local['energy_intervals'][i][0]-ground['energy_intervals'][0][1]
        ehi = local['energy_intervals'][i][1]-ground['energy_intervals'][0][0]
        if elo <= 0:
            raise ArithmeticError('positive charged excitation not certified')
        weights.append((lo, hi))
        gaps.append((elo, ehi))
        rows.append({'index': i, 'energy_interval': decimal_bounds(elo, ehi, 14),
                     'weight_interval': decimal_bounds(lo, hi, 28),
                     'strictly_positive_weight_proved': lo > 0,
                     'energy_diagnostic': float((elo+ehi)/2),
                     'weight_diagnostic': float(raw*raw/(ground['norms'][0]*local['norms'][i]))})
    spectral_weight = tuple(sum((w[j] for w in weights), F(0)) for j in (0, 1))
    spectral_first = tuple(sum((w[j]*e[j] for w, e in zip(weights, gaps, strict=True)), F(0)) for j in (0, 1))
    direct = source_moments(ground)
    weight = intersect(spectral_weight, direct['weight'])
    first = intersect(spectral_first, direct['first'])
    if weight[0] <= 0 or first[0] <= 0:
        raise ArithmeticError('nonzero composite response not certified')
    mean = first[0]/weight[1], first[1]/weight[0]
    return {'reference': 'unique centered N8 hidden-spin-zero eta-zero ground',
            'energy_unit': 'dimensionless source-clock energy; hbar_model=1',
            'source_selection': selection,
            'complete_block_lines': len(rows),
            'strictly_positive_lines': sum(r['strictly_positive_weight_proved'] for r in rows),
            'weight_positivity_unresolved_lines': sum(not r['strictly_positive_weight_proved'] for r in rows),
            'omitted_block_lines': 0, 'weight_outside_reducing_block': '0',
            'rows': rows,
            'ground_energy_interval': decimal_bounds(*ground['energy_intervals'][0], 14),
            'total_weight_interval': decimal_bounds(*weight, 14),
            'first_energy_moment_interval': decimal_bounds(*first, 14),
            'conditional_mean_energy_interval': decimal_bounds(*mean, 14),
            'scaled_probe_success_interval': decimal_bounds(weight[0]/9, weight[1]/9, 14),
            'independent_sum_rules': {
                'spectral_weight_interval': decimal_bounds(*spectral_weight, 14),
                'direct_AstarA_interval': decimal_bounds(*direct['weight'], 14),
                'spectral_first_energy_moment_interval': decimal_bounds(*spectral_first, 14),
                'direct_Astar_commutator_interval': decimal_bounds(*direct['first'], 14),
                'source_rows': list(direct['source_rows']),
                'source_norm_bounds': [str(v) for v in direct['source_row_sum_norm_bounds']],
                'independent_of_composite_eigenvectors': True},
            'exact_certificate': {
                'Gram_defect_upper': str(local['gram_defect']),
                'complete_spectrum_radius_upper': str(local['complete_radius']),
                'maximum_unit_eigenvector_error_upper': str(max(local['vector_errors']))},
            'all_dark_charged_states_exhausted': False,
            'finite_gaps_identified_as_relativistic_rest_masses': False,
            'physical_particle_species_identified': False,
            'native_detector_executed': False,
            'autonomous_localized_defect_proved': False,
            'continuum_or_gravity_conclusion': False,
            'eligible_empirical_comparisons_added': 0}
