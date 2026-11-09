from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_momentum import momentum as m


@pytest.fixture(scope='module')
def member():
    return contract.member(2)


def test_exact_fourier_instrument_and_original_gram_law():
    c=m.exact_control()
    assert [r['p'] for r in c['rows']] == [0,s.Rational(3,16),0,s.Rational(1,16),s.Rational(3,4),0,0,0]
    assert sum(r['p'] for r in c['rows'])==1
    assert c['instrument']['complete_effect']==s.eye(2)
    for k in range(4):
        value=sum(s.I**(k*(j-l))*c['gram'][j,l] for j in range(4) for l in range(4))/16
        assert s.simplify(value-c['rows'][k]['p'])==0


def test_coherent_address_differs_from_classical_mixing_and_PVM():
    c=m.exact_control()
    assert c['rows'][1]['p'] != c['classical_address_mixture_blank_each_probability']
    effect=apparatus.clean(c['rows'][1]['K'].H*c['rows'][1]['K'])
    assert set(effect.eigenvals()) == {0,s.Rational(1,4)}
    assert effect*effect != effect
    assert apparatus.clean(c['H']*c['rho']-c['rho']*c['H'])==s.zeros(2)


def test_entangled_reference_is_retained_and_unchanged_on_complete_sum():
    c=m.exact_control()
    assert c['rho'].det()>0
    assert c['reference_after']==c['rho'].T
    assert c['instrument']['outcomes']==8
    assert c['instrument']['source_translation_operator_used'] is False


def test_nonzero_aperture_edge_identity_is_not_periodic_covariance():
    c=m.exact_control(); probes=c['probes']; fifth=s.zeros(2)
    shifted=m.fourier_instrument((*probes[1:],fifth))
    for k in range(4):
        phase=s.I**k
        defect=apparatus.clean(shifted['rows'][k]['K']-phase*c['rows'][k]['K'])
        assert defect==phase*(fifth-probes[0])/4
        assert defect != s.zeros(2)
    assert m.aperture_bound()['operator_phase_covariance_error_upper']==F(1,2)


def test_fixed_physical_aperture_does_not_remove_generator_edge_error():
    for n in range(2,9):
        row=m.aperture_bound(2**n,1,F(1,2**n))
        assert row['aperture_period_length']==1
        assert row['one_step_error_divided_by_spacing_upper']==2
    assert m.aperture_bound(4,5)['operator_phase_covariance_error_upper']==2


def test_zero_and_identity_probes_keep_every_label():
    for b,branch in [(s.zeros(2),1),(s.eye(2),0)]:
        out=m.fourier_instrument([b]*4)
        assert len(out['rows'])==8
        for row in out['rows']:
            assert row['K']==(s.eye(2) if row['branch']==branch and row['wavevector_index']==0 else s.zeros(2))


def test_invalid_probe_and_large_exact_control_rejected():
    for menu in [[],[2*s.eye(2)],[s.Matrix([[0,1],[0,0]])],[s.eye(2),s.eye(3)],[s.eye(2)]*17]:
        with pytest.raises((ValueError,ArithmeticError)):
            m.fourier_instrument(menu)


def test_noninvariant_faithful_state_has_actual_spatial_implementation():
    c=m.implementation_control(); U=c['source_rotation']; V=c['standard_GNS_unitary']
    assert apparatus.clean(V.H*V)==s.eye(4)
    for A in [s.Matrix([[0,1],[0,0]]),s.diag(1,0),s.Matrix([[0,-s.I],[s.I,0]])]:
        assert apparatus.clean(V*s.kronecker_product(s.eye(2),A)*V.H-s.kronecker_product(s.eye(2),U*A*U.H))==s.zeros(4)
    assert c['prepared_vector_fixed'] is False
    assert c['mean_ground_excess']==0
    assert c['momentum_variance']==2-s.sqrt(3)


def test_phase_shift_ambiguity_preserves_observable_covariance():
    c=m.implementation_control(); P=c['standard_GNS_generator']
    A=s.kronecker_product(s.eye(2),s.diag(1,-1))
    assert (P+7*s.eye(4))*A-A*(P+7*s.eye(4)) == P*A-A*P
    assert s.simplify((c['prepared_vector'].H*(P+7*s.eye(4))*c['prepared_vector'])[0])==7


@pytest.mark.parametrize('side',[3,8,32])
def test_actual_open_incidence_is_not_cyclic_translation(side):
    out=m.open_box_geometry(side)
    assert out['native_bond_before'][1][0]-out['native_bond_before'][0][0]==1
    assert out['cyclic_image_bond'][1][0]-out['cyclic_image_bond'][0][0]==1-side
    assert out['cyclic_image_is_native_bond'] is False
    assert out['incoming_links_at_origin'] != out['incoming_links_at_shifted_origin']


def test_authenticated_primary_aperture_with_complete_reference_error(member):
    row=m.common_payload(member)
    assert row['fine_outcomes']==8
    assert all(p['active'] and p['commutes_all_local_Gauss_constraints'] for p in row['aperture_probes'])
    assert F(row['complete_error_to_ideal_filtered_aperture_upper'])<F(1,10)
    assert row['actual_source_momentum_PVM_proved'] is False
    assert row['selected_preparation_momentum_concentration_proved'] is False
    assert row['source_translation_unitary_applied'] is False
    witness=row['separate_unfiltered_original_blank_weight_each_label_lower']
    assert F(witness['coefficient'])==F(1,360)
    assert witness['negative_exponent']>0
    assert row['unfiltered_positive_bound_transferred_to_fixed_time_filter'] is False
    json.dumps(row,allow_nan=False)


def test_source_and_false_momentum_mutations_fail(member):
    receipt=m.bind_member(member);m.verify_bound(member,receipt)
    fake=deepcopy(receipt);fake['result']['actual_source_momentum_PVM_proved']=True
    fake['result_sha256']=contract.digest(fake['result'])
    with pytest.raises(ValueError,match='differs'):m.verify_bound(member,fake)
    bad=deepcopy(member);bad['source']['spacing']='1/8'
    with pytest.raises(ValueError):m.common_payload(bad)
    with pytest.raises(ValueError,match='distinct active'):m.common_payload(contract.member(0),spacing='1')
    with pytest.raises(ValueError,match='integer lattice step'):m.common_payload(member,spacing='1/5')


def test_all_source_paths_exist():
    for path in m.SOURCE_INPUTS: assert Path(path).is_file(),path


def test_distinct_shared_edge_plaquettes_have_original_Haar_witness():
    anchors=[(j,0,0) for j in range(4)]
    row=m.blank_loop_geometry(anchors)
    assert len(row['pairwise_unique_edge_witnesses'])==6
    assert set(map(tuple,row['edge_sets'][0])) & set(map(tuple,row['edge_sets'][1]))
    # The inherited ten orthogonal character coefficients each have amplitude1/30.
    assert 10*F(1,30)**2==F(row['normalized_blank_Gram_diagonal'])
    for witness in row['pairwise_unique_edge_witnesses']:
        i,j=witness['pair']; edge=witness['unique_edge']
        assert edge in row['edge_sets'][j] and edge not in row['edge_sets'][i]
    with pytest.raises(ValueError):m.blank_loop_geometry([(0,0,0),(0,0,0)])
