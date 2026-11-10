"""Exact bare neutral transfer moments and a full-H remainder bound."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ..model import clean
from ..periodic_family.model import Configuration, apply_state, norm2, kinetic_terms, local_yukawa
from ..periodic_family.checks import norm_bound
from ..neutral_mobility.model import preparation, arrival_projector
from ..neutral_mobility.checks import local_preparation_checks
from ..oriented_links.rotations import axial_rotation_column
from ...physical_attachment.chiral import representation
from .clebsch import native_charges
from .model import column, hopping_column, gauss_column


def inner(left, right):
    return s.expand(sum(s.conjugate(value) * right.get(key, 0) for key, value in left.items()))


def spectator_key(key, sites, links, scalars):
    """Exact unchanged coordinates of an actual Hamiltonian summand."""
    mask, flux, scalar = key
    for site in sites:
        mask &= ~(((1 << 30) - 1) << (30 * site))
    return (mask, tuple(v for j, v in enumerate(flux) if j not in links),
            tuple(v for j, v in enumerate(scalar) if j not in scalars))


def local_column(key, site, config):
    mask, links, scalars = key
    offset = 30 * site; local = (mask >> offset) & ((1 << 30) - 1)
    out = {}
    for (target, scalar), value in local_yukawa(local, scalars[site],
            (config.y_up, config.y_down, config.y_lepton)).items():
        changed = list(scalars); changed[site] = scalar
        out[mask ^ ((local ^ target) << offset), links, tuple(changed)] = value
    return out


def supported_bilinear(left, right, config=Configuration()):
    """Compute <left,H right>, excluding only exact spectator mismatches.

    Every diagonal, onsite and kinetic summand is enumerated. No source-state
    coordinate is deleted from either input vector.
    """
    diagonal = s.Integer(0)
    intersection = left.keys() & right.keys()
    for key in intersection:
        _, links, scalars = key
        energy = (config.link_penalty * sum(v != 0 for v in links)
                  + config.scalar_penalty * sum(v != 0 for v in scalars))
        diagonal += s.conjugate(left[key]) * right[key] * energy
    rows = [{'kind': 'diagonal', 'candidate_columns': len(intersection),
             'contribution': str(s.expand(diagonal))}]
    total = diagonal
    summands = [('Yukawa', site, (site,), (), (site,),
                 lambda key, site=site: local_column(key, site, config))
                for site in range(config.side ** 3)]
    summands.extend(('kinetic', j, (term.left_cell, term.right_cell), term.path, (),
                     lambda key, term=term: hopping_column(key, (term,)))
                    for j, term in enumerate(kinetic_terms(config)))
    for kind, index, sites, links, scalars, action in summands:
        left_spectators = {spectator_key(key, sites, links, scalars) for key in left}
        candidates = [key for key in right
                      if spectator_key(key, sites, links, scalars) in left_spectators]
        contribution = s.Integer(0); nonzero = 0
        fields = {q: s.Integer(0) for q in native_charges()}
        for key in candidates:
            value = s.Integer(0)
            for target, amplitude in action(key).items():
                if target not in left:
                    continue
                term = s.conjugate(left[target]) * amplitude
                value += term
                if kind == 'kinetic':
                    changed = key[0] ^ target[0]
                    mode = (changed & -changed).bit_length() - 1
                    charge = int(6 * representation()[2][-1][mode % 15, mode % 15])
                    fields[charge] += right[key] * term
            value = s.expand(value)
            if value:
                nonzero += 1
                contribution += right[key] * value
        contribution = s.expand(contribution)
        row = {'kind': kind, 'index': index, 'sites': list(sites),
            'path': list(links), 'candidate_columns': len(candidates),
            'nonzero_contributing_columns': nonzero,
            'contribution': str(contribution)}
        if kind == 'kinetic':
            assert s.expand(sum(fields.values())) == contribution
            row['middle_field_contributions_six_Y'] = {str(q): str(s.expand(v)) for q, v in fields.items()}
        rows.append(row)
        total += contribution
    return s.expand(total), rows


@lru_cache(maxsize=1)
def transfer_checks():
    config = Configuration()
    psi = preparation(config=config)
    hs = apply_state(psi, lambda key: column(key, config))
    hhs = apply_state(hs, lambda key: column(key, config))
    assert len(hs) == 433 and norm2(hs) == s.Rational(1651, 8)
    assert len(hhs) == 62799 and norm2(hhs) == s.Rational(15494801, 128)
    physical = 0
    for state in (psi, hs):
        for site in range(8):
            for current in range(12):
                assert not apply_state(state, lambda key: gauss_column(key, site, current, config))
                physical += 1
    rows = []; target_data = []
    for spin in (0, 1):
        target = preparation(cell=4, spin=spin, config=config)
        ht = apply_state(target, lambda key: column(key, config))
        hht = apply_state(ht, lambda key: column(key, config))
        assert len(hht) == len(hhs) and norm2(hht) == norm2(hhs)
        moments = [inner(target, psi), inner(target, hs), inner(ht, hs),
                   inner(ht, hhs), inner(hht, hhs)]
        assert moments == [0] * 5
        fifth, support = supported_bilinear(hht, hhs, config)
        assert fifth == (0 if spin == 0 else -151 * s.I / 32)
        target_data.append(target)
        rows.append({'target_spin': spin, 'moments_orders_zero_through_five':
                     [str(value) for value in (*moments, fifth)],
                     'complete_second_image_terms': len(hht),
                     'complete_middle_summands': support})
    # Rotation about z fixes both N=2 sites geometrically. The AP seam at
    # x=1 reverses the target's spin-up eigenvalue, hence this amplitude is
    # exactly zero at every time, not just through the tested powers.
    rotate = lambda key: axial_rotation_column(key, 2, config)
    assert apply_state(psi, rotate) == {key: -s.I * value for key, value in psi.items()}
    assert apply_state(target_data[0], rotate) == {key: s.I * value for key, value in target_data[0].items()}
    assert apply_state(target_data[1], rotate) == {key: -s.I * value for key, value in target_data[1].items()}
    assert not arrival_projector(psi, 4, config)
    for target in target_data:
        assert arrival_projector(target, 4, config) == target
    bound = norm_bound(config)
    assert bound == 1896
    # The inherited word norm estimate uses only ||U_ab||<=1. This holds
    # here on all119 states by the complete transporter Gram calculation.
    sixth_moment_bound = bound ** 2 * norm2(hhs)
    moment = s.Rational(151, 32)
    time = s.cancel(3 * moment / sixth_moment_bound)
    leading = s.cancel(moment * time ** 5 / s.factorial(5))
    remainder = s.cancel(sixth_moment_bound * time ** 6 / s.factorial(6))
    assert remainder == leading / 2
    lower, upper = s.cancel(leading ** 2 / 4), s.cancel(9 * leading ** 2 / 4)
    assert 0 < lower < upper < 1
    return {'source_configuration': {'side': 2, 'spacing': 'pi',
            'all_three_Yukawa_coefficients': '1', 'scalar_and_link_penalties': '1'},
        'input': 'D_origin,up^dagger Omega, normalized, with the spin-singlet fermionic reference retained',
        'origin': [0, 0, 0], 'destination': [1, 0, 0],
        'full_source_first_image_terms': len(hs),
        'full_source_first_image_norm_squared': str(norm2(hs)),
        'full_source_second_image_terms': len(hhs),
        'full_source_second_image_norm_squared': str(norm2(hhs)),
        'all_cell_Gauss_checks_on_preparation_and_first_image': physical,
        'target_spin_moments': rows,
        'same_spin_transfer_all_times': '0 by exact pi-z rotation eigenvalues and the AP seam',
        'bare_transfer_amplitude_leading_term': '-151*t^5/3840',
        'complete_source_operator_norm_bound': str(bound),
        'mixed_sixth_derivative_absolute_bound': str(sixth_moment_bound),
        'amplitude_remainder_bound': 't^6 * ||H^3 target|| * ||H^3 source|| / 720',
        'controlled_time': str(time), 'leading_amplitude_absolute_value': str(leading),
        'outward_amplitude_remainder': str(remainder),
        'arrival_projector': 'sum_spin |D_destination,spin^dagger Omega><D_destination,spin^dagger Omega|',
        'detector_support': 'complete finite source: every other matter, scalar and link memory is blank',
        'complete_arrival_probability_interval': [str(lower), str(upper)],
        'complete_no_arrival_probability_interval': [str(1 - upper), str(1 - lower)],
        'local_packet_and_star_projector_contains_exclusive_arrival': True,
        'complete_local_arrival_probability_interval': [str(lower), '1'],
        'complete_local_no_arrival_probability_interval': ['0', str(1 - lower)],
        'local_detector_support': 'destination site and its incident links',
        'retained_outcomes': ['bare neutral arrival', 'all other records'],
        'reference_retained_probability': '1',
        'both_input_spins_same_bounds': 'exact pi-x rotation exchanges spin preparations',
        'extra_neutral_hopping_term': False,
        'discarded_source_states_or_dynamical_terms': False,
        'Omega_is_stationary_or_ground': False,
        'physical_mass_or_relativistic_propagation_inferred': False}


def preparation_checks():
    return local_preparation_checks()
