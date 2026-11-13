"""Actual shifted Gram correction and complete-complement inertia bounds."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import json
import sympy as sy
from model.openwave_refined_haar_spectrum.certificate import data, offblock_bound
from model.openwave_refined_haar_source.gauss import electric
from model.openwave_coupled_haar_spectrum.certificate import source_digest, algebraic_bounds
from model.openwave_colour_string.fixed_witness import spectral_enclosure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


@lru_cache(None)
def blocks(which):
    source = data(which)
    retained = set(source['states'])
    rows = defaultdict(list)
    for j, column in enumerate(source['actions']):
        for target, value in column:
            if target not in retained:
                rows[target].append((j, value))
    grouped = {}
    for target, row in rows.items():
        energy = electric(target) if which == 'refined' else F(sum(k*(k+2) for k in target), 4)
        diagonal = source['alpha']*energy + source['faces']*source['beta']
        terms = grouped.setdefault(diagonal, {})
        for j, x in row:
            for k, y in row:
                terms[j, k] = terms.get((j, k), 0) + sy.Rational(source['beta']**2)*x*y
    grouped = {d: {ij: sy.simplify(v) for ij, v in entries.items() if v}
               for d, entries in grouped.items()}
    norm, count = offblock_bound(source['states'], source['actions'])
    assert count == len(rows)
    return source, grouped, norm*source['beta']**2, count


@lru_cache(None)
def weighted_remainder(which, z):
    """Use the actual squared-resolvent Gram, not only the smallest Q energy."""
    z = F(z)
    source, grouped, norm2, count = blocks(which)
    floor = source['floor']; magnetic = source['faces']*source['beta']
    assert z < floor
    entries = {}
    for diagonal, terms in grouped.items():
        for ij, value in terms.items():
            entries[ij] = entries.get(ij, 0) + value/sy.Rational((diagonal-z)**2)
    entries = {ij: sy.simplify(value) for ij, value in entries.items()}
    rows = [F(0) for _ in source['states']]
    for (i, j), value in entries.items():
        lo, hi = algebraic_bounds(value)
        rows[i] += max(abs(lo), abs(hi))
    gram = sy.SparseMatrix(source['H'].rows, source['H'].cols, entries)
    assert gram == gram.T
    upper = max(rows)
    ratio = magnetic/(floor+magnetic-z)
    weighted = magnetic*upper/(1-ratio)
    unweighted = norm2*magnetic/((floor+magnetic-z)*(floor-z))
    return min(weighted, unweighted), {'squared_resolvent_Gram_sha256': source_digest(gram),
            'squared_resolvent_Gram_norm_upper_exact': str(upper),
            'weighted_full_self_energy_error_upper_exact': str(weighted),
            'unweighted_full_self_energy_error_upper_exact': str(unweighted)}


def shifted(which, z):
    z = F(z)
    source, grouped, norm2, count = blocks(which)
    floor = source['floor']
    assert z < floor
    entries = {}
    for diagonal, terms in grouped.items():
        for ij, value in terms.items():
            entries[ij] = entries.get(ij, 0) + value/sy.Rational(diagonal-z)
    matrix = sy.SparseMatrix(source['H'].rows, source['H'].cols,
                             {ij: sy.simplify(v) for ij, v in entries.items()})
    assert matrix == matrix.T
    magnetic = source['faces']*source['beta']
    d0 = floor+magnetic
    error, _ = weighted_remainder(which, z)
    return source['H']-matrix, matrix, error


def tail(norm2, floor, magnetic, z, degree):
    """Remainder after words 0..degree; every omitted Haar state remains."""
    z = F(z)
    assert z < floor and type(degree) is int and degree >= 0
    return norm2*magnetic**(degree+1)/((floor+magnetic-z)**(degree+1)*(floor-z))


def patch(length):
    """Exact local incidence and a full-Q orthogonal two-loop lower bound."""
    vertices = tuple(product(range(length+1), repeat=3))
    faces = []
    for a in range(3):
        for b in range(a+1, 3):
            for x in vertices:
                if x[a] == length or x[b] == length:
                    continue
                corners = []
                for i, j in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    y = list(x); y[a] += i; y[b] += j; corners.append(tuple(y))
                faces.append(frozenset(corners))
    p = faces[0]
    far = [f for f in faces if not f & p]
    alpha, beta = F(32, length), F(length**3, 32)
    magnetic = beta*len(faces)
    # P=1_{E<5}, so C>=5alpha on every full-Haar representation.
    floor = 5*alpha
    return {'subdivisions': length, 'vertices': len(vertices), 'elementary_faces': len(faces),
            'spacing_exact': str(F(1, length)), 'declared_g_squared_exact': str(F(64, length**2)),
            'alpha_exact': str(alpha), 'beta_exact': str(beta),
            'same_clock_electric_subdivision_identity': length*alpha == 32,
            'plaquettes_disjoint_from_fixed_boundary_probe': len(far),
            'full_QHP_norm_squared_lower_exact': str(beta**2*len(far)/4),
            'global_zero_energy_Neumann_ratio_upper_exact': str(magnetic/(floor+magnetic)),
            'Q_floor_for_E_less_than_5_projection_exact': str(floor),
            'local_Lieb_Robinson_word_rate_exact': str(32*beta),
            'third_patch_is_incidence_and_Haar_lower_bound_not_spectrum': length == 3}


def endpoint(which, label, packet):
    z = F(packet['energy_exact'])
    effective, gram, error = shifted(which, z)
    assert source_digest(effective) == packet['exact_matrix_sha256']
    frame = spectral_enclosure(effective, packet['complete_frame'])
    intervals = [tuple(map(F, row)) for row in frame['all_eigenvalue_intervals']]
    negative = sum(hi-z+error < 0 for lo, hi in intervals)
    nonpositive_possible = sum(lo-z-error <= 0 for lo, hi in intervals)
    # Congruence with positive C-z transfers exact inertia, without a cutoff H.
    if label == 'lower':
        assert negative == nonpositive_possible == 1
    else:
        assert negative >= 2
    return {'energy_exact': str(z), 'energy': float(z),
            'full_self_energy_error_upper_exact': str(error),
            'weighted_complement_certificate': weighted_remainder(which, z)[1],
            'actual_shifted_Gram_sha256': source_digest(gram),
            'actual_shifted_Gram_nonzero_entries': len(gram.todok()),
            'negative_eigenvalues_lower': negative,
            'negative_eigenvalues_upper': nonpositive_possible,
            'finite_complete_frame_certificate': frame}


def run():
    witnesses = json.loads((HERE/'numerical_witnesses.json').read_text())
    inherited = json.loads((ROOT/'model/openwave_refined_haar_spectrum/results.json').read_text())
    records = []
    for which in ('coarse', 'refined'):
        source, grouped, norm2, count = blocks(which)
        lower = endpoint(which, 'lower', witnesses[which]['lower'])
        upper = endpoint(which, 'upper', witnesses[which]['upper'])
        interval = (F(lower['energy_exact']), F(upper['energy_exact']))
        prior = next(r for r in inherited['records'] if r['source'] == which)
        ground = tuple(map(F, prior['full_ground_energy_interval_exact']))
        oldgap = tuple(map(F, prior['full_gap_interval_exact']))
        gap = (interval[0]-ground[1], interval[1]-ground[0])
        # Intersect independently valid full-source enclosures, retaining their ancestry.
        gap = (max(gap[0], oldgap[0]), min(gap[1], oldgap[1]))
        assert 0 < gap[0] < gap[1]
        oldwidth, width = oldgap[1]-oldgap[0], gap[1]-gap[0]
        assert width < oldwidth
        magnetic, floor = source['faces']*source['beta'], source['floor']
        comparison_z = F(57) if which == 'refined' else F(96)
        old_error = tail(norm2, floor, magnetic, 0, 0) + norm2*comparison_z/(floor*(floor-comparison_z))
        new_error = tail(norm2, floor, magnetic, comparison_z, 0)
        assert new_error < old_error
        records.append({'source': which, 'complete_retained_dimension': source['H'].rows,
            'all_direct_full_Q_action_rows': count, 'distinct_omitted_electric_energies': len(grouped),
            'full_offblock_norm_squared_upper_exact': str(norm2),
            'full_Q_floor_exact': str(floor), 'bounded_magnetic_norm_upper_exact': str(magnetic),
            'lower_endpoint': lower, 'upper_endpoint': upper,
            'full_first_excited_energy_interval_exact': list(map(str, interval)),
            'inherited_full_ground_energy_interval_exact': list(map(str, ground)),
            'full_gap_interval_exact': list(map(str, gap)), 'full_gap_interval': list(map(float, gap)),
            'prior_gap_interval': list(map(float, oldgap)), 'gap_interval_width_improvement_factor': float(oldwidth/width),
            'comparison_energy_exact': str(comparison_z), 'old_unshifted_bound_exact': str(old_error),
            'shifted_bound_exact': str(new_error),
            'finite_word_tail_bounds_degrees_0_through_3_exact': [str(tail(norm2, floor, magnetic, comparison_z, degree)) for degree in range(4)],
            'words_beyond_degree_zero_executed': False,
            'all_omitted_representations_and_intertwiners_retained': True})
    patches = [patch(n) for n in (1, 2, 3)]
    assert all(r['same_clock_electric_subdivision_identity'] for r in patches)
    assert patches[2]['plaquettes_disjoint_from_fixed_boundary_probe'] > patches[1]['plaquettes_disjoint_from_fixed_boundary_probe']
    return {'schema': 'oph.energy_dependent_full_Haar_Schur.v1',
            'verdict': 'PASS_SHIFTED_FULL_COMPLEMENT_INERTIA_AND_REFINEMENT_OBSTRUCTION',
            'records': records, 'local_refinement_controls': patches,
            'resolvent_domain': 'real z strictly below the stated full-Q floor',
            'finite_word_extension': 'b2 M^(N+1)/((d+M-z)^(N+1)*(d-z)) for N>=0; only N=0 matrix executed',
            'global_volume_obstruction': '||QHP||^2 >= beta^2 N_far/4, N_far >= N_faces-48 on a cubic graph',
            'uniform_global_offblock_bound_claimed': False,
            'local_algebra_dynamics_shortcut': 'Bounded plaquette incidence and existing interaction-picture Lieb-Robinson theorem; no interacting-vacuum replacement',
            'nontrivial_continuum_or_Clay_result': False, 'physical_GeV_calibration': None,
            'new_empirical_fit_inputs': []}
