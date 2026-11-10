"""Native irreps, physical static-reference records and full-H tails."""
from collections import Counter
from functools import lru_cache
from math import isqrt

import sympy as s

from model.gauge_defect_transport.transport import complete_partial_transport
from ...chiral_source.model import local_encoding_checks
from ..oriented_links.checks import zero
from ..peter_weyl_links.clebsch import transporter
from .model import (representation, source, irreps, native_charges, Configuration,
    apply_state, norm2, empty_key, preparation, preparation_pulse_column, column,
    gauss_column, spin_column, directed_arrivals, complete_readout, arrival_projector,
    local_arrival_projector, reference_current, one_body_column, sequence)


@lru_cache(maxsize=1)
def dictionary_checks():
    currents = representation()[2]
    fundamental = representation()[1]
    color_metric = s.ImmutableMatrix(8, 8, lambda a, b: 2 * s.trace(fundamental[a] * fundamental[b]))
    inverse_color_metric = color_metric.inv()
    names = {-4: 'u^c', 1: 'Q', 6: 'e^c', -3: 'L', 2: 'd^c'}
    rows = []
    for charge in native_charges():
        block = irreps()[charge]['native_indices']; d = len(block)
        small = tuple(t.extract(block, block) for t in currents)
        color = sum((inverse_color_metric[a, b] * small[a] * small[b]
                     for a in range(8) for b in range(8)), s.zeros(d))
        weak = sum((t * t for t in small[8:11]), s.zeros(d))
        color_value, weak_value = color[0, 0], weak[0, 0]
        assert color == color_value * s.eye(d) and weak == weak_value * s.eye(d)
        assert small[11] == s.Rational(charge, 6) * s.eye(d)
        weights = Counter(tuple(str(t[a, a]) for t in (small[6], small[7], small[10], small[11]))
                          for a in range(d))
        # Reuse the complete-instrument lemma with the actual119-state matrix.
        # Its operation angle is not used as the source Hamiltonian clock.
        maps = complete_partial_transport(transporter(charge))
        size = maps['forward'].rows
        identity = s.ImmutableSparseMatrix(size, size, {(j, j): 1 for j in range(size)})
        for moved, held in (('forward', 'forward_hold'), ('reverse', 'reverse_hold')):
            assert zero(maps[moved].H * maps[moved] + maps[held].H * maps[held] - identity)
        assert zero(maps['reversible'].H * maps['reversible'] - identity)
        rows.append({'conventional_native_label': names[charge], 'six_Y': charge,
            'native_indices': list(block), 'dimension': d,
            'color_quadratic_Casimir': str(color_value),
            'weak_quadratic_Casimir': str(weak_value),
            'actual_joint_Cartan_weights': [{'weights': list(w), 'multiplicity': n}
                                            for w, n in sorted(weights.items())],
            'Q_equals_T3_plus_Y_eigenvalues': sorted({str(small[10][a, a] + small[11][a, a]) for a in range(d)}),
            'dual_reference_gauge_current': '-T^transpose',
            'actual119_hardware_complete_transport_instrument': True})
    return {'native_families': 1, 'native_internal_dimensions_sum': sum(q['dimension'] for q in rows),
        'native_irrep_rows': rows, 'actual_current_count': 12,
        'color_current_metric_twice_fundamental_trace': [[str(v) for v in row] for row in color_metric.tolist()],
        'complete_reference_CAR_modes': 30,
        'common_source_native_encoding': local_encoding_checks(),
        'representation_names_are_measured_particle_identifications': False}


@lru_cache(maxsize=1)
def preparation_and_detector_checks():
    config = Configuration(); vacuum = {empty_key(): s.Integer(1)}
    vacuum_first = apply_state(vacuum, column)
    rows = []
    for charge in native_charges():
        psi = preparation(charge); h = apply_state(psi, column)
        assert norm2(psi) == 1
        assert apply_state(vacuum, lambda key: preparation_pulse_column(key, charge)) == psi
        assert apply_state(psi, lambda key: preparation_pulse_column(key, charge)) == vacuum
        assert all((key[0][0].bit_count() + key[1].bit_count()) % 2 == 0 for key in psi)
        for axis in range(3):
            assert not apply_state(psi, lambda key: spin_column(key, axis))
        current_checks = 0
        for state in (psi, h):
            for site in range(8):
                for current in range(12):
                    assert not apply_state(state, lambda key: gauss_column(key, site, current))
                    current_checks += 1
        targets = directed_arrivals(charge)
        for _, target in targets:
            assert norm2(target) == 1
            for site in (0, 4):
                for current in range(12):
                    assert not apply_state(target, lambda key: gauss_column(key, site, current))
        arrived, other = complete_readout(h, charge)
        assert norm2(arrived) == s.Rational(1, 8)
        assert norm2(h) == norm2(arrived) + norm2(other)
        assert arrival_projector(arrived, charge) == arrived
        assert not arrival_projector(other, charge)
        assert local_arrival_projector(h, charge) == arrived
        assert not local_arrival_projector(psi, charge)
        assert not local_arrival_projector(vacuum, charge)
        assert not local_arrival_projector(vacuum_first, charge)
        for _, target in targets:
            assert local_arrival_projector(target, charge) == target
        # Matrix generators preserve every selector used by the local detector.
        block = irreps()[charge]['native_indices']
        matter_selector = s.diag(*[int(a in block) for a in range(15)])
        scalar_selector = s.diag(1, 0, 0, 0, 0)
        from ...chiral_higgs.model import higgs_memory
        for current, t in enumerate(representation()[2]):
            assert zero(t * matter_selector - matter_selector * t)
            scalar = higgs_memory()[1][current]
            assert zero(scalar * scalar_selector - scalar_selector * scalar)
            for sign in (1, -1):
                projector = s.ImmutableSparseMatrix(119, 119, {(j, j): 1
                    for j, label in enumerate(source.labels())
                    if label[0] == sign and label[1] in block and label[2] in block})
                for end in (0, 1):
                    generator = source.endpoint_current(current, end)
                    assert zero(generator * projector - projector * generator)
        # The opposite reference orientation gives a nonzero physical residual.
        bad = apply_state(psi, lambda key: {(key[0], ref): -2 * value
                    for ref, value in one_body_column(key[1], reference_current(11)).items()})
        assert norm2(bad) == s.Rational(charge * charge, 9)
        # An odd intervening spectator fixes the graded preparation sign.
        state, ref = empty_key(); mask, links, scalars = state
        spectator_key = ((1 << 30, links, scalars), ref)
        pulse = preparation_pulse_column(spectator_key, charge)
        d = len(block)
        for key, value in pulse.items():
            occupied = (key[0][0] ^ (1 << 30)).bit_length() - 1
            reference = key[1].bit_length() - 1
            step = sequence(1 << 30, ((occupied, True), (240 + reference, True)))
            assert step is not None
            spin = occupied // 15
            assert value == (-1) ** spin * step[1] / s.sqrt(2 * d)
        assert apply_state(pulse, lambda key: preparation_pulse_column(key, charge)) == {spectator_key: s.Integer(1)}
        rows.append({'six_Y': charge, 'preparation_terms': len(psi),
            'complete_first_image_terms': len(h), 'complete_first_image_norm_squared': str(norm2(h)),
            'joint_Gauss_checks_on_preparation_and_first_image': current_checks,
            'exclusive_arrival_basis_rank': 2, 'exclusive_arrival_first_image_terms': len(arrived),
            'each_directed_string_derivative_norm_squared': '1/16',
            'complete_local_and_exclusive_arrival_derivative_norm_squared': '1/8',
            'vacuum_local_arrival_derivative_norm_squared': '0',
            'wrong_reference_orientation_residual_norm_squared': str(norm2(bad)),
            'odd_intervening_spectator_preparation_sign_checked': True})
    return {'reference_preparation': 'even spin-singlet with a static dual-R15 CAR reference at origin',
        'pulse': 'exp(-i*pi*K/2), K=|Phi_r><vac|+|vac><Phi_r| on the active local pair sector',
        'orthogonal_complement_of_preparation_orbit': 'identity under the unitary pulse',
        'complete_global_readout': ['one of two compensated one-link arrival strings', 'all other records'],
        'complete_local_readout': ['one native-r fermion, blank scalar and matching joining flux', 'all other records'],
        'local_detector_support': 'destination site and its six incident stored links',
        'local_detector_reference_support': 'identity on the complete static reference',
        'all_detector_selectors_commute_with_actual_gauge_currents': True,
        'static_reference_discarded': False, 'rows': rows}


@lru_cache(maxsize=1)
def full_powers():
    vacuum_key = source.empty_key(); vacuum_first = source.column(vacuum_key)
    vacuum_second = apply_state(vacuum_first, source.column)
    vacuum_norm = norm2(vacuum_second)
    assert vacuum_norm == 141832 and len(vacuum_second) == 33841
    rows = []
    for mode in range(30):
        key = (1 << mode, vacuum_key[1], vacuum_key[2])
        first = source.column(key)
        second = apply_state(first, source.column)
        rows.append({'mode': mode, 'complete_first_image_terms': len(first),
            'complete_first_image_norm_squared': str(norm2(first)),
            'complete_second_image_terms': len(second),
            'complete_second_image_norm_squared': str(norm2(second))})
    families = []
    expected = {-4: 17893585, 1: 17894609, 6: 17896657, -3: 18057745, 2: 17895633}
    for charge in native_charges():
        modes = [15 * spin + a for spin in range(2) for a in irreps()[charge]['native_indices']]
        values = [s.Rational(rows[mode]['complete_second_image_norm_squared']) for mode in modes]
        assert all(v == s.Rational(expected[charge], 128) for v in values)
        families.append({'six_Y': charge,
            'complete_reference_resolved_second_image_terms': sum(rows[mode]['complete_second_image_terms'] for mode in modes),
            'complete_prepared_second_image_norm_squared': str(sum(values) / len(values))})
    maximum = max([vacuum_norm] + [s.Rational(row['complete_prepared_second_image_norm_squared']) for row in families])
    bound = isqrt(int(maximum)) + 1
    assert bound == 377 and maximum < bound ** 2
    time = s.Rational(1, 4 * bound)
    lower, upper = 25 * time ** 2 / 576, 25 * time ** 2 / 64
    background = time ** 2 / 64
    contrast = time ** 2 / 36
    assert lower - background == contrast > 0
    return {'streaming_rule': 'H acts identically on orthogonal fixed reference labels; every one of30 source components is applied twice',
        'source_mode_powers': rows, 'five_prepared_reference_resolved_powers': families,
        'vacuum_first_image_terms': len(vacuum_first),
        'vacuum_first_image_norm_squared': str(norm2(vacuum_first)),
        'vacuum_second_image_terms': len(vacuum_second),
        'vacuum_second_image_norm_squared': str(vacuum_norm),
        'common_outward_second_image_norm_bound': bound, 'common_source_time': str(time),
        'common_full_state_first_order_remainder': str(time / 8),
        'prepared_exclusive_arrival_probability_interval': [str(lower), str(upper)],
        'prepared_exclusive_complement_probability_interval': [str(1 - upper), str(1 - lower)],
        'prepared_local_arrival_probability_interval': [str(lower), '1'],
        'prepared_local_complement_probability_interval': ['0', str(1 - lower)],
        'preparation_off_local_arrival_probability_interval': ['0', str(background)],
        'preparation_off_local_complement_probability_interval': [str(1 - background), '1'],
        'strict_local_record_contrast_lower_bound': str(contrast),
        'contrast_applies_to_all_five_native_blocks_same_clock': True,
        'base_Yukawa_or_quantum_link_terms_removed': False,
        'source_or_reference_states_discarded': False,
        'prepared_states_identified_as_ground_or_free_particles': False}
