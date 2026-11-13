import json
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
from .transport import (constant_identity, differential_identity, fock_operators,
                        pair_dictionary, transform, full_flux_records)


def test_all_modes_are_retained_and_inverse_roundtrip():
    dictionary, pairs = pair_dictionary()
    assert len(dictionary) == 48 and len(pairs) == 24
    for i, (j, sign) in dictionary.items():
        for creation in (False, True):
            poly = {((i, creation),): sy.Integer(1)}
            assert transform(transform(poly), inverse=True) == poly
    i, j, _ = pairs[0]
    density = {((i, True), (i, False)): sy.Integer(1)}
    image = transform(density)
    clipped = {w: c for w, c in image.items() if len(w) != 2 or w[0][1] != w[1][1]}
    assert transform(clipped, inverse=True) != density


def test_complete_reference_and_nonzero_work():
    result = constant_identity()
    assert result['work_preserved_exact'] == '-1/2'
    assert result['entangled_reference_full_joint_output_identity']
    assert sum(map(F, result['complete_local_records_exact'])) == 1


def test_differential_pair_and_square_control():
    result = differential_identity()
    assert F(result['field_connection_pair_creation_norm_squared_exact']) == F(1296, 625)
    a, b = fock_operators()
    charge = a.H*a-b.H*b
    assert charge*(a.H*b.H)-(a.H*b.H)*charge == sy.zeros(4)


def test_full_integer_flux_and_original_input_outcomes():
    for row in full_flux_records():
        ground = tuple(map(F, row['full_four_state_Fock_ground_interval_exact']))
        assert ground[1] < F(row['all_omitted_flux_lower_exact'])
        assert F(row['connection_omitted_n0_offdiagonal_squared_exact']) > 0
        for result in row['complete_original_input_records']:
            assert result['original_input_probabilities_diagnostic'][0] == .5
            assert abs(sum(result['original_input_probabilities_diagnostic'])-1) < 1e-12


def test_receipt_covers_actual_full_polynomial_without_source_substitution():
    receipt = json.loads(Path(__file__).with_name('results.json').read_text())
    native = receipt['native_transport']
    assert len(native['all_source_factor_records']) == 114
    assert native['all_24_actual_weak_currents_retain_pair_terms']
    assert not native['kinetic_or_gauge_terms_added_to_native_H']
    assert not receipt['source_equivalence_obstruction']['onto_unitary_between_these_complete_carriers']
    assert receipt['additional_common_source_poles_identified'] == []
