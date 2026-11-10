"""Full matrix, representation, path-word and occupied-source controls."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from model.particle_gauge_band.quantum_link_completion import wedge_sign
from ...physical_attachment.chiral import representation
from ..model import clean
from ..periodic_family.model import Configuration, native_blocks, apply_state, norm2
from ..oriented_links.checks import zero, rotation_word_checks
from ..oriented_links.rotations import axial_rotation_column
from .clebsch import (irreps, native_charges, intertwiner, all_channels,
                      coordinate, link, transporter, carrier_index)
from .model import (labels, label_index, endpoint_current, orientation,
                    path_transitions, empty_key, column, gauss_column,
                    translation_column)


@lru_cache(maxsize=1)
def clebsch_checks():
    channels, census = all_channels()
    assert len(census) == 55 and len(channels) == 37
    for q in irreps():
        commutant, info = intertwiner(0, q, q)
        assert len(commutant) == 1 and info['multiplicity'] == 1
    subsets = representation()[0]
    masks = tuple(sum(1 << j for j in subset) for subset in subsets)
    wedge_channels = []
    for (r, source, target), embeddings in channels.items():
        if source not in native_charges() or target not in native_charges():
            continue
        ri, si, ti = (irreps()[q]['native_indices'] for q in (r, source, target))
        raw = s.zeros(len(ri) * len(si), len(ti))
        for a, left in enumerate(ri):
            for b, right in enumerate(si):
                for c, output in enumerate(ti):
                    if masks[left] | masks[right] == masks[output]:
                        raw[len(si) * a + b, c] = wedge_sign(masks[left], masks[right])
        gram = raw.H * raw
        assert gram[0, 0] > 0 and gram == gram[0, 0] * s.eye(len(ti))
        wedge = raw / s.sqrt(gram[0, 0])
        assert len(embeddings) == 1
        assert wedge * wedge.H == embeddings[0] * embeddings[0].H
        wedge_channels.append([r, source, target])
    assert len(wedge_channels) == 5
    rows = []
    for (r, source, target), embeddings in channels.items():
        rows.append({'field': r, 'source': source, 'target': target,
            'dimensions': [irreps()[q]['dimension'] for q in (r, source, target)],
            'multiplicity': len(embeddings),
            'orthonormal_embeddings': [
                {'rows': c.rows, 'cols': c.cols,
                 'entries': [[int(i), int(j), str(value)]
                             for (i, j), value in sorted(c.todok().items())]}
                for c in embeddings]})
    return {'actual_native_currents': 12, 'retained_irreps': 11,
        'exact_scalar_commutants': 11,
        'field_source_pairs_exhausted': 55,
        'nonzero_retained_channels': 37,
        'nonvacuum_to_nonvacuum_channels': 27,
        'all_channels_multiplicity_one': True,
        'full_current_intertwiner_identities': 37 * 12,
        'independent_native_exterior_wedge_projectors': wedge_channels,
        'census': list(census), 'channels': rows}


@lru_cache(maxsize=1)
def matrix_checks():
    j = orientation(); identity = s.eye(119)
    assert j.H * j == identity and j * j == identity
    electric = s.diag(0, *([1] * 118))
    parity = identity - 2 * electric
    ward = reversal = 0
    for current, t in enumerate(representation()[2]):
        left, right = (endpoint_current(current, end) for end in (0, 1))
        assert left.H == left and right.H == right and j * left * j == right
        assert zero(left * electric - electric * left)
        assert zero(right * electric - electric * right)
        for block in native_blocks():
            for a in block:
                for b in block:
                    u = link(a, b)
                    expected_left = sum((-t[a, c] * link(c, b) for c in block), s.zeros(119))
                    expected_right = sum((t[c, b] * link(a, c) for c in block), s.zeros(119))
                    assert zero(left * u - u * left - expected_left)
                    assert zero(right * u - u * right - expected_right)
                    ward += 2
    for block in native_blocks():
        for a in block:
            for b in block:
                assert j * link(a, b) * j == link(b, a).H
                reversal += 1
    rows = []
    for r in native_charges():
        d = irreps()[r]['dimension']; v = transporter(r)
        p, q = v.H * v, v * v.H
        assert p.H == p and q.H == q
        assert zero(p * p - p) and zero(q * q - q)
        rank = sum(irreps()[source]['dimension'] * irreps()[target]['dimension'] * len(cs)
                   for (field, source, target), cs in all_channels()[0].items() if field == r)
        assert s.trace(p) == s.trace(q) == rank
        blank = s.ImmutableSparseMatrix(119 * d, d, {(119 * a, a): 1 for a in range(d)})
        assert blank.H * p * blank == s.eye(d)
        rows.append({'six_Y': r, 'native_dimension': d,
                     'full_internal_link_dimension': 119 * d,
                     'transporter_nonzero_entries': len(v.todok()),
                     'exact_Gram_projection_rank': int(rank), 'operator_norm': '1'})
    ec = coordinate(6, 0, 0)
    source = carrier_index(-3, 0, 0)
    targets = [(row, value) for (row, col), value in ec.todok().items() if col == source]
    assert len(targets) == 1 and targets[0][0] != 0 and abs(targets[0][1]) == 1
    assert not zero(parity * ec + ec * parity)
    wrong = s.ImmutableSparseMatrix(119, 119, {(0, 0): 1, **{
        (label_index()[-sign, a, b], col): 1
        for col, (sign, a, b) in enumerate(labels()[1:], 1)}})
    a, b = native_blocks()[0][:2]
    assert wrong * link(a, b) * wrong != link(b, a).H
    scaled = s.sqrt(6) * transporter(1)
    assert s.trace(scaled.H * scaled) == 6 * s.trace(transporter(1).H * transporter(1))
    return {'full_link_basis': 119, 'all_matrix_columns_retained': True,
        'complete_endpoint_Ward_identities': ward,
        'orientation_coordinate_identities': reversal,
        'unitary_orientation_and_endpoint_exchange': True,
        'blocks': rows,
        'excited_to_excited_witness': {'field_six_Y': 6, 'source_label': list(labels()[source]),
            'target_label': list(labels()[targets[0][0]]), 'amplitude': str(targets[0][1])},
        'link_excitation_parity_anticommutation': False,
        'wrong_orientation_transpose_rejected': True,
        'omitted_Haar_normalization_norm_squared_mutation': '6'}


def path_checks():
    checked = 0
    for adjoint in (False, True):
        for value in range(119):
            actual = {(a, b, targets[0]): coefficient
                      for a, b, targets, coefficient in path_transitions((value,), adjoint)}
            expected = {}
            for block in native_blocks():
                for a in block:
                    for b in block:
                        matrix = link(a, b).H if adjoint else link(a, b)
                        for (row, col), coefficient in matrix.todok().items():
                            if col == value:
                                expected[a, b, row] = coefficient
            assert actual == expected
            checked += 1
    rows = []
    for length in (1, 2, 3):
        sums = defaultdict(lambda: s.Integer(0))
        for a, b, targets, value in path_transitions((0,) * length):
            sums[b] += value ** 2
        assert set(sums) == set(range(15)) and all(s.simplify(v) == 1 for v in sums.values())
        rows.append({'length': length, 'all_fifteen_input_norms_squared': '1'})
    return {'exhaustive_one_link_columns_both_adjoint_branches': checked,
        'blank_path_executions': rows,
        'arbitrary_path_operator_norm': '<=1 by composition of complete-carrier contractions',
        'distinct_edge_blank_path_operator_norm': '1'}


def actual_hamiltonian_checks():
    config = Configuration(y_up=s.Rational(2, 3), y_down=s.Rational(3, 5),
                           y_lepton=s.Rational(5, 7), scalar_penalty=s.Rational(7, 11),
                           link_penalty=s.Rational(11, 13))
    base = empty_key(config); cases = []
    for block in native_blocks():
        a, b = block[0], block[-1]
        for sign in (1, -1):
            links, scalars = list(base[1]), list(base[2])
            links[0] = label_index()[sign, a, b]
            links[4] = label_index()[-sign, b, a]
            scalars[0], scalars[3] = 1, 4
            mask = (1 << a) | (1 << (30 + 15 + b)) | (1 << (120 + b))
            cases.append((mask, tuple(links), tuple(scalars)))
    cases.append(base)
    checked = 0
    for key in cases:
        h = column(key, config)
        for axis in range(3):
            rotate = lambda state: axial_rotation_column(state, axis, config)
            assert clean(apply_state(rotate(key), lambda state: column(state, config))) == clean(apply_state(h, rotate))
            checked += 1
    key = cases[2]
    for axis in range(3):
        translate = lambda state: translation_column(state, axis, config)
        assert clean(apply_state(translate(key), lambda state: column(state, config))) == clean(apply_state(column(key, config), translate))
    for current in (0, 9, 11):
        gauss = lambda state: gauss_column(state, 0, current, config)
        assert clean(apply_state(gauss(key), lambda state: column(state, config))) == clean(apply_state(column(key, config), gauss))
    # An excited-to-excited e^c hop changes both endpoint stars. The basis
    # witness complements the full-coordinate anticommutator check.
    from ..neutral_mobility.model import star_signature
    links = list(base[1]); links[0] = carrier_index(-3, 0, 0)
    witness = 1 << (30 * 4 + 15 + 9), tuple(links), base[2]
    before = star_signature(witness, config)[0]
    commutator = {target: value * (before - star_signature(target, config)[0])
                  for target, value in column(witness, config).items()
                  if before != star_signature(target, config)[0]}
    assert commutator and norm2(commutator) > 0
    return {'complete_occupied_basis_rotation_controls': checked,
        'complete_basis_translation_controls': 3, 'actual_H_G_controls': 3,
        'star_commutator_witness_terms': len(commutator),
        'star_commutator_witness_norm_squared': str(norm2(commutator)),
        'star_witness_is_physical_preparation': False,
        'all_native_blocks_and_both_link_orientations_tested': True}
