"""Exact spatial residues, visible multiplicities and source rate sum rules."""
from fractions import Fraction as F
from math import comb

import numpy as np

from ..charged_composite_search import certificate as prior_certificate
from ..charged_composite_search import packet as prior_packet
from ..centered_particle_spectrum import certificate as ground_certificate
from ..centered_particle_spectrum import packet as ground_packet
from . import source

decimal = prior_certificate.decimal_bounds


def interval(value):
    return tuple(map(F, value))


def multiply(a, b):
    products = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(products), max(products)


def add_intervals(values):
    values = tuple(values)
    return tuple(sum((value[j] for value in values), F(0)) for j in (0, 1))


def dot(a, b):
    return sum(int(x)*int(y) for x, y in zip(a, b, strict=True))


def spatial_moments(local):
    block = local['block']
    y = block['full_basis'].astype(object) @ local['vectors']
    moments = []
    for shift in range(8):
        rows, signs = source.parent.fock_action(5, 2, (0, 1, 2), shift)
        values = np.sum(y*(signs[:, None]*y[rows]), axis=0)
        moments.append([F(int(v), 1 << 120)/norm
                        for v, norm in zip(values, local['norms'], strict=True)])
    dimensions = []
    enclosures = []
    for i, error in enumerate(local['vector_errors']):
        square_sum = sum((row[i]**2 for row in moments), F(0))
        candidates = [d for d in range(1, 9)
                      if square_sum-32*error <= F(8, d) <= square_sum+32*error]
        if len(candidates) != 1:
            raise ArithmeticError('visible spatial rank is not uniquely certified')
        dimensions.append(candidates[0])
        enclosures.append([(max(F(-1), row[i]-2*error), min(F(1), row[i]+2*error))
                           for row in moments])
    return dimensions, enclosures


def direct_moments(ground):
    # These independent full-source identities attach the Clifford triangle
    # bounds to the sparse occupation matrices used for the moment readout.
    source.kinetic_source_anchor()
    norm_certificate = source.commutator_certificate()
    if any(F(row['kinetic_rate_norm_bound']) > F(21, 16) or
           F(row['centered_contact_rate_norm_bound']) > F(3, 2)
           for row in norm_certificate['site_checks']):
        raise ArithmeticError('source rate norm exceeds the certified moment budget')
    g = ground['block']['full_basis'].astype(object) @ ground['vector_numerators'][:, 0]
    normalization = dot(g, g)
    epsilon = ground['vector_errors'][0]
    maps = [source.direct_rate_maps(x) for x in range(8)]
    action = prior_certificate.sparse_integer_action
    q = [action(row['probe'], g) for row in maps]
    rate16 = [action(row['full_rate16'], g) for row in maps]
    matrices = [[] for _ in range(3)]
    errors = (18*epsilon, F(135, 8)*epsilon, F(2025, 128)*epsilon)
    for order in range(3):
        for x in range(8):
            row = []
            for y in range(8):
                if order == 0:
                    value = F(dot(q[x], q[y]), normalization)
                elif order == 1:
                    value = F(dot(q[x], rate16[y]), 16*normalization)
                else:
                    value = F(dot(rate16[x], rate16[y]), 256*normalization)
                row.append((value-errors[order], value+errors[order]))
            matrices[order].append(row)
    kinetic = action(maps[0]['kinetic_rate16'], g)
    force = action(maps[0]['centered_contact_rate16'], g)
    definitions = {
        'kinetic_first_moment': (F(dot(q[0], kinetic), 16*normalization), F(63, 8)),
        'centered_contact_first_moment': (F(dot(q[0], force), 16*normalization), F(9)),
        'kinetic_rate_norm_squared': (F(dot(kinetic, kinetic), 256*normalization), F(441, 128)),
        'centered_contact_rate_norm_squared': (F(dot(force, force), 256*normalization), F(9, 2)),
        'kinetic_contact_interference': (F(2*dot(kinetic, force), 256*normalization), F(63, 8)),
    }
    components = {name: (value-coefficient*epsilon, value+coefficient*epsilon)
                  for name, (value, coefficient) in definitions.items()}
    if components['centered_contact_rate_norm_squared'][0] <= 0:
        raise ArithmeticError('nonzero source contact contribution not certified')
    prior_certificate.intersect(add_intervals(components[name] for name in
                                ('kinetic_first_moment', 'centered_contact_first_moment')),
                                matrices[1][0][0])
    prior_certificate.intersect(add_intervals(components[name] for name in
                                ('kinetic_rate_norm_squared', 'centered_contact_rate_norm_squared',
                                 'kinetic_contact_interference')), matrices[2][0][0])
    prior_certificate.intersect(components['centered_contact_rate_norm_squared'],
                                multiply((F(1, 4), F(1, 4)), matrices[0][0][0]))
    return matrices, components


def verify():
    witness = prior_packet.load_witness()
    inherited = prior_certificate.verify(witness)
    local = prior_certificate.eigensystem(witness)
    previous = ground_packet.load_witness()
    ground = ground_certificate.eigensystem('ground', previous['ground_columns'],
                                           previous['ground_energies'], previous['ground_vectors'])
    # The union count consumes the actual prior linear-CAR certificate,
    # including every positive weight and its independently certified rank.
    linear = ground_certificate.verify(previous)
    dimensions, moments = spatial_moments(local)
    source.spatial_commutant_certificate()
    adjacency = source.parent.ground.band_data()[0]
    sqrtlo, sqrthi = prior_certificate.sqrt_bounds(F(3))
    signs = source.parent.shift_gram_signs()
    rows = []
    residue_intervals = []
    energies = []
    for line, dimension, overlaps in zip(inherited['rows'], dimensions, moments, strict=True):
        weight = interval(line['weight_interval'])
        energy = interval(line['energy_interval'])
        if overlaps[1][0] > 0:
            band = 1
        elif overlaps[1][1] < 0:
            band = -1
        else:
            raise ArithmeticError('actual kinetic band sign is not certified')
        if adjacency[0,1] != 1:
            raise ArithmeticError('source band sign anchor changed')
        residue = []
        for x in range(8):
            entries = []
            for y in range(8):
                translated = multiply(weight, multiply((F(int(signs[x,y])), F(int(signs[x,y]))),
                                                       overlaps[x ^ y]))
                coefficient = ((F(1), F(1)) if x == y else
                               multiply((F(band*int(adjacency[x,y])), F(band*int(adjacency[x,y]))),
                                        (1/sqrthi, 1/sqrtlo)))
                projected = multiply(weight, coefficient)
                entries.append(prior_certificate.intersect(translated, projected))
            residue.append(entries)
        rate = [[multiply((energy[0]**2, energy[1]**2), item) for item in row] for row in residue]
        residue_intervals.append(residue)
        energies.append(energy)
        rows.append({'index': line['index'], 'energy_interval': line['energy_interval'],
                     'origin_weight_interval': line['weight_interval'],
                     'visible_spatial_rank': dimension,
                     'actual_source_kinetic_band_sign': band,
                     'exact_residue_formula': 'w_i (I8 + sign_i A/sqrt(3))',
                     'visible_states_per_charge_with_spin_ladder': 4*dimension,
                     'residue_matrix_intervals': [[decimal(*entry, 28) for entry in row] for row in residue],
                     'rate_residue_matrix_intervals': [[decimal(*entry, 28) for entry in row] for row in rate],
                     'normalized_signed_translation_overlaps': [decimal(*v, 20) for v in overlaps]})
    direct, components = direct_moments(ground)
    moments_report = []
    for order in range(3):
        spectral = []
        intersection = []
        for x in range(8):
            spectral_row = []
            intersected_row = []
            for y in range(8):
                values = (multiply((energy[0]**order, energy[1]**order), residue[x][y])
                          for energy, residue in zip(energies, residue_intervals, strict=True))
                value = add_intervals(values)
                spectral_row.append(decimal(*value, 18))
                intersected_row.append(decimal(*prior_certificate.intersect(value, direct[order][x][y]), 18))
            spectral.append(spectral_row)
            intersection.append(intersected_row)
        moments_report.append({'energy_power': order, 'spectral_matrix_intervals': spectral,
                               'independent_source_matrix_intervals': [[decimal(*entry, 18) for entry in row]
                                                                       for row in direct[order]],
                               'intersection_matrix_intervals': intersection,
                               'all_64_entries_checked': True})
    spatial = sum(dimensions)
    composite = 4*spatial
    old = linear['reachable_N7_dimension_all_linear_CAR_probes']
    full = comb(16, 7)
    if dimensions != [4]*138 or composite+old > full:
        raise ArithmeticError('unexpected complete composite multiplicity')
    return {'state': inherited['reference'], 'energy_unit': inherited['energy_unit'],
            'rows': rows, 'positive_lines': len(rows),
            'actual_source_kinetic_band_counts': {str(sign): sum(row['actual_source_kinetic_band_sign'] == sign
                                                               for row in rows) for sign in (-1, 1)},
            'visible_spatial_dimension_fixed_highest_spin': spatial,
            'visible_N7_dimension_with_spin_ladder': composite,
            'visible_N9_dimension_with_spin_and_eta_ladders': composite,
            'linear_CAR_N7_dimension': old, 'orthogonal_union_N7_dimension': composite+old,
            'full_N7_dimension': full, 'unseen_N7_dimension_after_union': full-composite-old,
            'unseen_N9_dimension_after_ladder_union': full-composite-old,
            'source_current_moment_matrices': moments_report,
            'origin_kinetic_contact_decomposition': {key: decimal(*value, 18) for key, value in components.items()},
            'all_138_poles_and_all_64_spatial_entries_retained': True,
            'rate_adds_no_cyclic_states_to_same_probe_family': True,
            'full_energy_eigenspace_multiplicities_claimed': False,
            'all_dark_charged_states_exhausted': False,
            'spatial_labels_are_continuum_momentum': False,
            'hidden_spin_is_Lorentz_spin': False,
            'source_clock_calibrated_in_SI': False,
            'physical_particle_species_identified': False,
            'native_detector_execution_claimed': False,
            'continuum_interaction_or_gravity_proved': False,
            'eligible_empirical_comparisons_added': 0}
