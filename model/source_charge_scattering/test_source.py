"""Independent controls for source pair algebra and complete finite response."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from . import algebra,response,packet


@pytest.fixture(scope='module')
def data():return algebra.data()


@pytest.fixture(scope='module')
def receipt():return packet.build()


def test_exact_sylvester_accepts_and_rejects():
    assert algebra.positive_leading_minors([[2,1],[1,2]])==(2,3)
    with pytest.raises(ArithmeticError):algebra.positive_leading_minors([[1,2],[2,1]])
    with pytest.raises(ValueError):algebra.positive_leading_minors([[1,2],[0,1]])


def test_actual_H_blocks_reduce_whole_source(data):
    op,g,norm,error,checked,creation,sectors=data
    for q,s in sectors.items():
        assert len(s['labels'])==28
        assert np.array_equal(op['H'].astype(object)@s['W'],s['W']@s['block'])
        assert np.array_equal(op['twice_charge'].astype(object)@s['W'],4*q*s['W'])
        assert s['true_gram_floor']>F(9,10**10)


def test_all_pair_inputs_have_source_charge_and_complete_instruments(data):
    op,g,norm,error,checked,creation,sectors=data
    identity=np.eye(256,dtype=np.int64)
    for x,y in algebra.PAIRS:
        complete=np.zeros_like(identity)
        for q in(-1,1):
            for r in(-1,1):
                k=creation[x,q]@creation[y,r]
                assert np.array_equal(k,creation[y,r]@creation[x,q])
                assert np.array_equal(op['twice_charge']@k-k@op['twice_charge'],2*(q+r)*k)
                complete+=k.T@k
        assert np.array_equal(complete,16*identity)


def test_connected_source_walk_and_edge_unit_extraction(data):
    sectors=data[-1]
    for s in sectors.values():
        walk=s['walk'];h=s['block']
        assert set(walk)==set(range(28))
        for i,j in zip(walk,walk[1:]):
            pi=np.zeros((28,28),dtype=object);pj=pi.copy();pi[i,i]=1;pj[j,j]=1
            compressed=pi@h@pj
            expected=np.zeros((28,28),dtype=object);expected[i,j]=h[i,j]
            assert h[i,j]!=0 and np.array_equal(compressed,expected)
    bad=np.zeros((3,3),dtype=int)
    with pytest.raises(ArithmeticError):algebra.spanning_walk(bad)


def test_Gram_floor_agrees_with_independent_singular_values(data):
    for s in data[-1].values():
        coords=np.array(s['coordinates'],dtype=float)/np.sqrt(float(s['denominator']))
        singular=np.linalg.svd(coords,compute_uv=False)
        assert min(singular)**2>float(s['true_gram_floor'])
        assert min(singular)>7e-5


def test_complete_pair_detector_has_no_dark_sector_vector(data):
    for s in data[-1].values():
        c=np.array(s['coordinates'],dtype=float)/np.sqrt(float(s['denominator']))
        effect=c@c.T/28
        levels=np.linalg.eigvalsh(effect)
        assert levels[0]>float(s['true_gram_floor']/28)
        assert levels[-1]<1


@pytest.mark.parametrize('q',[-1,1])
def test_all_pair_clicks_against_full_source_exponential(q,data):
    op,g,norm,error,checked,creation,sectors=data
    gn=np.array(g,dtype=float)/np.sqrt(float(norm))
    b=np.array([creation[x,q]@creation[y,q]@gn/4 for x,y in algebra.PAIRS]).T
    eig,v=np.linalg.eigh(op['H']/op['scale'])
    incoming=algebra.PAIRS.index((0,7))
    for t in(F(),F(1,4)):
        full=v@(np.exp(-1j*float(t)*eig)*(v.T@b[:,incoming]))
        kernel=b.T@full;probs=abs(kernel)**2/28
        report=response.pair_readout(q,t)
        for index,row in enumerate(report['all28_detector_outcomes']):
            lo,hi=map(float,row['original_input_herald_and_detector_probability'])
            assert lo-1e-14<=probs[index]<=hi+1e-14
        dark=np.dot(b[:,incoming],b[:,incoming])-sum(probs)
        lo,hi=map(float,report['same_charge_input_and_dark_probability'])
        assert lo-1e-14<=dark<=hi+1e-14


@pytest.mark.parametrize('y',[1,3,7])
def test_connected_response_against_independent_H_evolution(y,data):
    op,g,norm,error,checked,creation,sectors=data
    gn=np.array(g,dtype=float)/np.sqrt(float(norm));indices=np.arange(256)
    zx=1-2*(indices&1);zy=1-2*((indices>>y)&1)
    e,v=np.linalg.eigh(op['H']/op['scale']);u=v@(np.exp(-1j*e/8)*(v.T@gn))
    w=v@(np.exp(-1j*e/8)*(v.T@(zy*gn)))
    chi=2*np.vdot(u,zx*w).imag
    report=response.connected_response(y)
    lo,hi=map(float,report['retarded_Z_response_interval'])
    assert lo-1e-14<=chi<=hi+1e-14
    assert (lo>0) or (hi<0)
    p=(1+chi/2)/2
    plo,phi=map(float,report['complete_interference_Y_plus_probability'])
    assert plo-1e-14<=p<=phi+1e-14


def test_actual_bond_mixed_energy_and_nonedge_exact_cancellation():
    for y in(1,2,4):
        row=response.connected_response(y)
        assert row['actual_source_edge']
        assert F(row['mixed_Z_kick_energy_interval'][1])<-1
    for y in(3,5,6,7):
        row=response.connected_response(y)
        assert not row['actual_source_edge']
        assert list(map(F,row['mixed_Z_kick_energy_interval']))==[0,0]


def test_finite_closure_not_promoted_to_large_source_or_particles(receipt):
    rows=receipt['actual_source_algebra']['finite_dimensional_scaling_boundary']
    assert rows[0]['two_local_pair_vectors']==rows[0]['total_charge2_sector_dimension']==28
    assert rows[2]['two_local_pair_vectors']==66 and rows[2]['total_charge2_sector_dimension']==495
    scope=receipt['complete_response']['scope']
    assert not scope['asymptotic_scattering_on_this_selected_source_proved']
    assert not scope['other_Weyl_contact_source_S_matrix_transferred']
    assert scope['physical_species_mass_or_cross_section'] is None
    assert not scope['pair_detector_and_connected_response_native_tapes_executed']


@pytest.mark.parametrize('target',['gram','walk','outcome','scope','binding'])
def test_full_producer_rejects_changed_claims_and_outcomes(target,receipt):
    bad=deepcopy(receipt)
    if target=='gram':bad['actual_source_algebra']['same_charge_pair_sectors'][0]['true_selected_source_Gram_lower']='0'
    elif target=='walk':bad['actual_source_algebra']['same_charge_pair_sectors'][0]['source_nonzero_transition_walk'].pop()
    elif target=='outcome':bad['complete_response']['all_pair_frame_readouts'][0]['all28_detector_outcomes'].pop()
    elif target=='scope':bad['complete_response']['scope']['asymptotic_scattering_on_this_selected_source_proved']=True
    else:bad['source_bindings']['Lean/QFT/SourceOperatorGeneration.lean']='0'*64
    with pytest.raises(ValueError):packet.validate(bad)


def test_requested_path_is_actually_checked(tmp_path,receipt):
    path=tmp_path/'requested.json';path.write_text(packet.canonical(receipt))
    command=[sys.executable,'-m','model.source_charge_scattering','--check',str(path)]
    process=subprocess.run(command,cwd=packet.ROOT,capture_output=True,text=True)
    assert process.returncode==0,process.stdout+process.stderr
    path.write_text('{"schema":"forged"}')
    process=subprocess.run(command,cwd=packet.ROOT,capture_output=True,text=True)
    assert process.returncode!=0


def test_source_bindings_cover_full_fresh_runtime():
    code='''
import json,sys
from pathlib import Path
from model.source_charge_scattering import packet
receipt=packet.build();root=packet.ROOT.resolve();imports=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        path=Path(name).resolve()
        if path.suffix=='.py' and path.is_relative_to(root):imports.add(path.relative_to(root).as_posix())
assert not imports-set(receipt['source_bindings']),sorted(imports-set(receipt['source_bindings']))
assert 'model/protected_defect_catalog/charge.py' in imports
print(json.dumps({'count':len(imports)}))
'''
    process=subprocess.run([sys.executable,'-c',code],cwd=packet.ROOT,capture_output=True,text=True)
    assert process.returncode==0,process.stdout+process.stderr
    assert json.loads(process.stdout)['count']>10


def test_continuum_pair_keeps_all_nine_records_and_original_weights():
    from . import continuum
    row=continuum.pair_bounds()
    k=F(row['same_source_single_charge_original_weight_lower'])
    assert k==F(1,1757600)
    assert F(row['spatial_box_mean_pair_probability_lower'])==k*k
    assert F(row['arbitrarily_far_rational_shift_pair_probability_lower'])==k*k/4
    assert len(row['all_preparation_records'])==9
    assert set(map(tuple,row['all_preparation_records']))=={(x,y) for x in(-1,0,1) for y in(-1,0,1)}
    assert not row['clustering_or_Fock_factorization_assumed']
    assert not row['bound_valid_at_every_large_separation']


def test_pair_energy_window_keeps_zero_energy_and_unnormalized_tail():
    from . import continuum
    row=continuum.pair_bounds(F(1,8))
    original=F(row['arbitrarily_far_rational_shift_pair_probability_lower'])
    energy=F(row['each_pair_unnormalized_first_energy_moment_upper'])
    cutoff=F(row['same_source_pair_upper_energy_cut'])
    kept=F(row['original_pair_weight_below_cut_lower'])
    assert energy==F(27,32)
    assert kept==original-energy/cutoff==original/2
    assert row['cut_includes_zero_energy']


def test_mean_lemma_cannot_assert_pointwise_clustering_or_all_shifts():
    from . import continuum
    control=continuum.exact_mean_controls()
    values=list(map(F,control['periodic_four_shift_correlations']))
    assert sum(values)/4==F(control['periodic_mean_equals_squared_first_moment'])>0
    assert values.count(0)==3
    assert F(control['invariant_component_mean_correlation'])>F(control['squared_mean_in_invariant_component_control'])
    assert not control['controls_claimed_as_actual_source_states']


def test_pair_bounds_follow_the_actual_replayed_order_certificate(receipt):
    row=receipt['separated_continuum_charge_pairs']
    source=row['replayed_selected_source_order_certificate']
    assert source['selected_magnetization_strict_lower']=='1/10'
    assert source['controls'][1]['existing_exact_infrared_producer']['divisions']==256
    assert row['positive_original_weight_at_arbitrarily_large_separation']
    assert not row['zero_energy_pair_weight_excluded']
    assert not row['asymptotic_particle_states_or_scattering_matrix_constructed']


def test_exact_positive_cube_ceiling_controls_rounding():
    from . import continuum
    for value in(F(),F(1),F(8),F(9),F(27,8),F(10**20,7)):
        k=continuum.positive_cube_ceiling(value)
        assert k>=1 and k**3>=value
        assert k==1 or (k-1)**3<value
    with pytest.raises(ValueError):continuum.positive_cube_ceiling(F(-1))
    with pytest.raises(ValueError):continuum.positive_cube_ceiling(.5)


@pytest.mark.parametrize('n',[1,2,3,4,8])
def test_complete_n_creation_search_box_and_energy_budget(n):
    from . import continuum
    row=continuum.configuration_bounds(n,F(3,2),F(1,8))
    p=F(1,1757600);mean=p**n;k=row['least_positive_integer_cube_multiplier']
    assert F(row['explicit_cube_half_width'])==F(3,2)*k
    assert F(row['anchored_relative_configuration_half_width'])==3*k
    assert row['first_center_can_be_translated_to_zero']
    assert F(row['near_tuple_original_weight_upper'])<=mean/4
    assert F(row['separated_tuple_integral_lower'])>=3*mean/4
    assert F(row['some_rational_separated_tuple_original_weight_lower'])==mean/4
    assert row['unit_sign_record_count']==2**n
    assert row['all_complete_preparation_record_count']==3**n
    assert row['nonzero_unit_sign_total_charges']==list(range(-n,n+1,2))
    energy=F(row['unnormalized_first_energy_moment_upper'])
    cutoff=F(row['same_source_upper_energy_cut'])
    assert energy==F(3,2)*(1-F(1,2**n))**2
    assert mean/4-energy/cutoff==F(row['original_weight_below_energy_cut_lower'])==mean/8
    assert row['zero_energy_included']
    assert not row['specified_good_tuple_computed']


def test_general_two_creation_bound_matches_retained_pair_convention():
    from . import continuum
    pair=continuum.pair_bounds(F(1,8));general=continuum.configuration_bounds(2,F(1),F(1,8))
    assert general['some_rational_separated_tuple_original_weight_lower']==pair['arbitrarily_far_rational_shift_pair_probability_lower']
    assert general['unnormalized_first_energy_moment_upper']==pair['each_pair_unnormalized_first_energy_moment_upper']
    assert general['same_source_upper_energy_cut']==pair['same_source_pair_upper_energy_cut']
    assert general['original_weight_below_energy_cut_lower']==pair['original_pair_weight_below_cut_lower']


def test_configuration_count_is_not_a_particle_count(receipt):
    result=receipt['separated_continuum_charge_pairs']
    assert result['every_integer_source_charge_sector_nonzero']
    assert result['each_integer_charge_energy_infimum_zero']
    assert not result['zero_energy_eigenstate_or_positive_frequency_particle_inferred']
    assert all(not row['charge_count_identified_with_particle_number'] for row in result['finite_configuration_bounds'])
