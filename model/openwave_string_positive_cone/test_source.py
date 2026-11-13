from fractions import Fraction as F
from math import comb
from .source import (local_colour_census, uniform_bound, shared_length_attachment,
                     incidence_certificates, three_quark_control)
from .certificate import complete_span_records, six_quark_negative_colour_control


def test_complete_local_carrier_and_sharp_bound():
    d = local_colour_census()
    assert d['number_sector_dimensions'] == [comb(24, n) for n in range(25)]
    assert len(d['irreducible_rows']) == 547
    assert d['sharp_nonvacuum_equality_rows'][0]['highest_weight'] == [4, 4]
    assert all(F(r['Casimir']) <= 3 * r['number'] for r in d['irreducible_rows'])


def test_physical_attachment_inside_full_positive_interval():
    d = shared_length_attachment()
    assert d['entire_required_interval_number_floor'] > F('15.510795')
    assert d['quadrupole_fm2_interval'][0] <= F('0.285699') <= d['quadrupole_fm2_interval'][1]
    assert uniform_bound(d['kappa'])['adjusted_local_colour_coefficient'] < 0
    assert uniform_bound(1200)['number_floor'] < 0
    assert d['a_fm'] * d['tension_GeV_per_fm'] == F(4, 3) * d['kappa'] * d['electron_unit_GeV']


def test_Gauss_bound_keeps_native_six_quark_negative_control():
    assert len(incidence_certificates()) == 16
    assert three_quark_control()['six_quark_counterexample']['electric_minus_half_local_Casimir'] == '-1/3'
    d = six_quark_negative_colour_control()
    assert d['colour_remainder_eigenvalue'] == '-1/3'
    assert d['all_eight_native_Gauss_currents_annihilate_vector_exactly']


def test_complete_span_records_include_preparation_failure():
    d = shared_length_attachment()
    for L in (2, 3, 5, 16):
        r = complete_span_records(L, d['kappa'], d['a_fm'], d['electron_unit_GeV'])
        assert sum(F(x['original_probability']) for x in r['records']) + F(r['retained_coloured_preparation_failure']) == 1
        assert len(r['records']) == L
        assert all(F(x['string_energy_GeV']) == d['tension_GeV_per_fm'] * F(x['separation_fm']) for x in r['records'])
