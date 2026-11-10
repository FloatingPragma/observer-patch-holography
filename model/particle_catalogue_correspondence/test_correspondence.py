from collections import Counter
import sympy as s
from .model import exterior_rows, centre_information, adjoint_witness, pair_witness, conjugation_witness


def test_full_native_weights_and_multiplicities():
    out = centre_information()
    assert [r['exterior_dimension'] for r in out['centre_classes']] == [4,6,6,4,6,6]
    match = out['actual_fifteen_current_match']
    assert sum(len(r['source_masks']) for r in match.values()) == 15
    assert match['charged_lepton_right_conjugate']['q6Y'] == '6'
    assert match['neutrino_left']['Q'] == '0'


def test_centre_loses_hypercharge_and_electric_charge():
    zero = [r for r in exterior_rows() if r['centre_class']==0]
    assert {r['q6Y'] for r in zero} == {-6,0,6}
    three = [r for r in exterior_rows() if r['centre_class']==3]
    assert {r['Q'] for r in three} == {'-1','0','1'}
    assert {r['q6Y'] for r in three} == {-3,3}


def test_centre_neutral_does_not_pass_full_gauss():
    out = adjoint_witness()
    assert out['singlet_rank']==1 and out['non_singlet_rank']==8
    assert out['Casimir_octet']=='3'


def test_equal_currents_do_not_override_physical_parity():
    out = pair_witness()
    assert out['total_q6Y']==2 and out['pair_CAR_parity']==0
    assert out['single_d_conjugate_CAR_parity']==1
    assert out['pair_baryon_number'] != out['single_d_conjugate_baryon_number']


def test_self_inverse_centre_is_not_self_conjugate_representation():
    out = conjugation_witness()
    assert out['complement_pairs']==16 and out['class_three_nonzero_q6Y']==[-3,3]
    assert out['centre_self_inverse_implies_Majorana_or_self_conjugate'] is False
