"""Full local CAR, actual cutoff-star locality and normalized spatial events."""
from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s
from scipy.sparse import csr_matrix, eye, kron

from . import model as m


def test_all_six_square_zero_local_fields_have_exact_signed_complete_effects():
    control=m.local_anticommutators();assert control['complete_components']==6
    for family in m.spin.local_controls()['families'].values():
        for A in family['operators']:
            assert A*A==s.zeros(64)
            assert max((A.H*A+A*A.H).eigenvals())==1
    assert all(sum(row['anticommutator_spectrum'].values())==64 for row in control['rows'])
    assert control['occupation_lower_constant']==6 and control['occupation_lower_slope']==3


def test_actual_open_cube_has_four_disjoint_odd_residual_star_colors():
    graph=m.cube_stars();assert len(graph['cells'])==8 and len(graph['edges'])==12
    assert len(graph['colors'])==4
    flattened=[x for color in graph['colors'] for x in color];assert set(flattened)==set(graph['cells'])
    for x,y in graph['colors']:
        assert sum(a!=b for a,b in zip(x,y))==3
        assert len(graph['stars'][x]['matter'])==4 and len(graph['stars'][x]['links'])==3
        assert not graph['stars'][x]['matter']&graph['stars'][y]['matter']
        assert not graph['stars'][x]['links']&graph['stars'][y]['links']


def test_actual_odd_fields_have_graded_disjoint_site_anticommutators():
    local=m.spin.local_controls()['families'];A=csr_matrix(np.array(local['G1']['operators'][0],complex))
    B=csr_matrix(np.array(local['H']['operators'][1],complex))
    parity=csr_matrix(np.diag([(-1)**i.bit_count() for i in range(64)]))
    left=kron(A,eye(64),format='csr');right=kron(parity,B,format='csr')
    for second in (right,right.getH()):
        difference=left@second+second@left;difference.eliminate_zeros()
        assert not difference.nnz


def test_arbitrary_global_compression_would_not_preserve_this_locality_argument():
    create0=m.spin.odd.creation_word((0,),2);create1=m.spin.odd.creation_word((1,),2)
    assert create0*create1.H+create1.H*create0==s.zeros(4)
    P=s.diag(1,1,1,0);A=P*create0*P;B=P*create1.H*P
    assert A*B+B*A!=s.zeros(4)  # actual product-link cutoff plus reducing Gauss is essential


def test_full_twelve_mode_adjoint_residuals_have_the_same_summed_caps():
    for label,cap in [('G1',F(19,3)),('H',F(12))]:
        value=m.adjoint_residual_control(label)
        assert value['adjoint_Gram_is_exact_particle_hole_conjugate']
        assert value['one_sign_squared_residual_cap']==cap
        assert value['both_sign_squared_residual_cap']==2*cap


def test_complete_Walsh_rotation_preserves_successes_and_requires_row_normalization():
    control=m.spatial_menu_control();assert control['total_records_before_clocks_and_meters']==104
    assert control['signed_success_records']==96 and control['retained_address_failure_records']==8
    cells=m.cube_stars()['cells'];W=np.array([[(-1)**sum(a*b for a,b in zip(q,x)) for x in cells] for q in cells])/np.sqrt(8)
    rng=np.random.default_rng(9103);K=rng.normal(size=(8,6,3,3))+1j*rng.normal(size=(8,6,3,3))
    transformed=np.einsum('qx,xjab->qjab',W,K)/np.sqrt(48)
    selected=np.einsum('qx,xjab->qjab',W,K)/np.sqrt(6)
    assert np.linalg.norm(transformed-selected/np.sqrt(8))<1e-13
    before=sum(A.conj().T@A for A in K.reshape(-1,3,3))/48
    after=sum(A.conj().T@A for A in transformed.reshape(-1,3,3))
    assert np.linalg.norm(before-after)<1e-13
    # Complex Fourier phases require conjugating the negative pointer rotation.
    U=np.exp(2j*np.pi*np.outer(np.arange(8),np.arange(8))/8)/np.sqrt(8)
    chosen=np.einsum('x,xab->ab',U[1],K[:,0])
    correct=np.einsum('x,xab->ab',U[1].conj(),K[:,0].conj().transpose(0,2,1))
    wrong=np.einsum('x,xab->ab',U[1],K[:,0].conj().transpose(0,2,1))
    assert np.linalg.norm(correct-chosen.conj().T)<1e-13
    assert np.linalg.norm(wrong-chosen.conj().T)>1


@pytest.fixture(scope='module')
def actual():return m.bind_member(m.contract.member(0))


def test_original_fixed_spatial_spin_events_have_controlled_positive_energy(actual):
    r=actual['result'];assert F(r['complete_signed_centered_second_moment_upper'])==F(110,3)
    assert F(r['all_negative_transfer_at_least_half_weight_upper'])==F(1,8)
    assert F(r['finite_reported_original_creation_weight_lower'])>F(1,5)
    assert list(map(F,r['finite_reported_creation_interval']))==[F(1,4),F(71,4)]
    menu=r['complete_Walsh_menu']
    assert F(menu['each_reported_positive_spatial_row_weight_lower'])>F(1,40)
    assert F(menu['each_reported_positive_spatial_row_coarse_global_error_lower'])>F(1,40)
    assert F(menu['sum_of_eight_simple_disjoint_event_lowers'])==F(1,5)
    for channel,lower in zip(r['spin_resolved_channels'],(F(1,16),F(9,64))):
        assert F(channel['reported_original_positive_creation_weight_lower'])>lower
        assert F(channel['complete_Walsh_each_positive_spatial_row_lower'])>lower/8
        assert F(channel['negative_transfer_weight_upper'])==F(channel['dimension'],48)
    kernel=r['source_spatial_kernel_reduction']
    assert kernel['source_cube_orbits']['ordered_pair_orbit_sizes']==[8,24,24,8]
    assert kernel['complete_menu_per_row_eigenmeasure_divisor']==48
    assert not kernel['actual_primary_off_diagonal_shell_values_computed']
    assert not kernel['time_transforms_claimed_real']
    assert r['source_input_is_original_selected_density']
    assert r['all_signed_failure_memory_reference_and_clock_outputs_retained']
    assert not r['source_conjugation_implies_equal_retained_reference_output_states']
    assert not r['fixed_spatial_analyzer_is_conserved_momentum'] and r['physical_mass_table_GeV']==[]


def test_binding_rejects_refinement_or_unearned_component_promotion(actual):
    m.verify_bound(m.contract.member(0),actual)
    changed=deepcopy(actual);changed['result']['each_individual_internal_component_has_positive_fixed_spatial_weight_proved']=True
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),changed)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))


def test_certificate_is_canonical_and_pins_the_literal_source():
    certificate=m.certificate();assert json.loads(json.dumps(certificate,sort_keys=True))==certificate
    assert len(certificate['members'])==1
    assert 'model/particle_primary_selected_spin_response/PROOF.md' in certificate['source_sha256']
