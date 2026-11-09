from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_joint_spectroscopy import spectroscopy as q
from model.particle_primary_variational_readout import variational
from model.particle_primary_endpoint_bounds import endpoint


@pytest.fixture(scope='module')
def member():return contract.member(2)


def test_exact_complete_joint_records_have_original_up_down_and_failures():
    c=q.exact_control();rows=c['rows']
    assert len(rows)==64
    observed={(r['initial_energy_bin'],r['branch'],r['wavevector_index'],r['final_energy_bin']):r['p']
              for r in rows if r['p']!=0}
    assert observed=={(0,0,1,1):s.Rational(3,16),(1,0,1,0):s.Rational(1,16),
                      (0,1,0,0):s.Rational(9,16),(1,1,0,1):s.Rational(3,16)}
    assert c['complete_original_probability_sum']==1
    assert sum((r['K'].H*r['K'] for r in rows),s.zeros(2))==s.eye(2)


def test_positive_two_kernel_law_matches_every_direct_matrix_probability():
    c=q.exact_control();Frows=c['energy_rows'];rho=c['rho']
    spatial={(r['branch'],r['wavevector_index']):r['K'] for r in c['spatial_rows']}
    for row in c['rows']:
        K=spatial[row['branch'],row['wavevector_index']]
        a=row['initial_energy_bin'];b=row['final_energy_bin']
        expected=sum(rho[i,i]*abs(K[j,i])**2*abs(Frows[a][i,i])**2*abs(Frows[b][j,j])**2
                     for i in range(2) for j in range(2))
        assert s.simplify(expected-row['p'])==0


def test_source_stationarity_does_not_preserve_purification_after_initial_record():
    c=q.exact_control();rho=c['rho']
    assert sum((f*rho*f.H for f in c['energy_rows']),s.zeros(2))==rho
    before=c['initial_reference_joint_before'];after=c['initial_reference_joint_after']
    assert before!=after
    assert s.trace(before*before)==1
    assert s.trace(after*after)==s.Rational(5,8)
    # The untouched reference marginal is unchanged, although correlations changed.
    marginal=s.Matrix([[sum(after[j*2+i,k*2+i] for i in range(2)) for k in range(2)] for j in range(2)])
    assert marginal==rho.T


def test_finite_clock_rejects_nonunitary_and_incomplete_inputs():
    with pytest.raises(ValueError):q.energy_instrument(s.diag(1,2),4)
    with pytest.raises(ValueError):q.energy_instrument(s.eye(2),32)
    c=q.exact_control()
    with pytest.raises(ValueError):q.joint_instrument(c['energy_rows'][:1],c['spatial_rows'])
    with pytest.raises(ValueError):q.joint_instrument(c['energy_rows'],c['spatial_rows'][:1])


def test_actual_alias_free_clock_and_error_schedule(member):
    plan=q.clock_plan(member)
    assert F(plan['clock_step'])*F(plan['source_norm_upper'])==F(1,4)
    assert plan['alias_probability_upper']==0
    assert F(plan['clock_plan']['complete_clock_tail_upper'])<=F(1,1024)
    assert F(plan['one_complete_clock_error_upper'])<=F(1,1024)
    J=plan['clock_plan']['clock_size'];r=plan['Cayley_substeps_per_base_tick']
    assert F(plan['one_complete_clock_error_upper'])==F(J-1,768*r*r)
    assert not plan['source_energy_zero_changed']
    with pytest.raises(ValueError):q.clock_plan(member,resolution=str(8*F(plan['source_norm_upper'])))


def test_clock_tail_formula_by_direct_noncentral_phase_kernel():
    # Exact algebraic phase between four grid centers, independently of the bound.
    U=s.diag(s.exp(-s.I*s.pi/4))
    rows=q.energy_instrument(U,4)
    probabilities=[s.simplify((f.H*f)[0]) for f in rows]
    assert s.simplify(sum(probabilities)-1)==0
    # A turn collar1/4 excludes both nearest centers, at distance1/8.
    distant=s.simplify(probabilities[2]+probabilities[3])
    bound=q.clock_bounds.tail_bound(4,F(1,4))
    assert float(distant)<=float(bound)
    assert distant>0


def test_original_record_adapter_keeps_clock_cooling_and_event_errors(member):
    result=q.observed_transfer_to_excitation(member,('1','2'),'1/2',additional_event_error='1/100')
    assert list(map(F,result['final_excitation_interval']))==[F(31,32),F(69,32)]
    allowances=sum(F(result[k]) for k in ['two_energy_tail_allowance','complete_clock_allowance',
                                         'additional_event_allowance','initial_excitation_tail_allowance'])
    assert F(result['original_final_excitation_weight_lower'])==max(F(0),F(1,2)-allowances)
    assert F(result['initial_excitation_tail_allowance'])==F(1,8)
    assert result['numerical_ground_energy_required'] is False
    assert result['supplied_event_probability_executed_here'] is False
    assert result['rest_frame_or_particle_pole_proved'] is False


def test_empty_or_inconsistent_record_bounds_are_not_normalized(member):
    row=q.observed_transfer_to_excitation(member,('0','1'),'0')
    assert row['original_final_excitation_weight_lower']=='0'
    with pytest.raises(ValueError):q.observed_transfer_to_excitation(member,('2','1'),'1/2')
    with pytest.raises(ValueError):q.observed_transfer_to_excitation(member,('-10','-9'),'1/2')
    with pytest.raises(ValueError):q.observed_transfer_to_excitation(member,('0','1'),'2')


def test_actual_broad_band_and_huge_clock_cost_are_factored(member):
    row=q.actual_zero_wavevector_band(member)
    w=row['original_positive_band_weight_lower'];N=w['negative_exponent'];cw=F(w['coefficient'])
    h=(N+1)//2;a=F(member['source']['spacing']);tau=F(row['clock_step'])
    j=row['clock_bins']['exponent'];r=row['Cayley_substeps_per_tick']['exponent']
    assert j>=N+h+apparatus._ceil_log2(16/(cw*tau/(28*a)))
    assert 2*r>=j+N+apparatus._ceil_log2(1/(48*cw))
    observed=row['proved_reported_original_band_weight_lower']
    assert F(observed['coefficient'])==3*cw/4 and observed['negative_exponent']==N
    assert row['wavevector_zero_identified_as_rest_frame'] is False
    assert row['broad_support_is_resolved_ridge_or_particle_branch'] is False
    assert row['huge_register_or_pulse_tape_allocated'] is False
    assert len(json.dumps(row))<10000


def test_complete_primary_plan_and_forged_scope_rejection(member):
    receipt=q.bind_member(member);q.verify_bound(member,receipt)
    row=receipt['result'];J=row['energy_plan']['clock_plan']['clock_size']
    assert row['complete_fine_outcomes']==8*J*J
    assert row['two_energy_register_dimensions']==[J,J]
    assert row['compact_two_bit_Ramsey_carrier_bound_applies'] is False
    assert row['initial_source_marginal_preserved'] is True
    assert row['initial_energy_record_preserves_joint_purification'] is False
    fake=deepcopy(receipt);fake['result']['spatial_label_is_conserved_momentum']=True
    fake['result_sha256']=contract.digest(fake['result'])
    with pytest.raises(ValueError):q.verify_bound(member,fake)
    assert row['continuum_limit_required_for_this_finite_instrument'] is False


def test_final_energy_adapter_regenerates_physical_variational_upper():
    member=contract.member(0);upper=variational.bind_member(member)
    bracket=q.ground_energy_bracket(member,upper_receipt=upper)
    lo,hi=map(F,bracket['ground_energy_interval'])
    assert lo==F(member['preparation_plan']['source']['spectral_lower'])
    assert 0<=hi-F(upper['result']['ground_energy_upper'])<F(1,1<<32)
    assert hi < -90*member['source']['side']**3/F(member['source']['spacing'])
    row=q.observed_final_energy_to_excitation(member,(hi+1,hi+2),'1/2',upper_receipt=upper)
    assert list(map(F,row['final_excitation_interval']))==[F(63,64),hi+F(129,64)-lo]
    error=sum(F(row[k]) for k in ['final_energy_tail_allowance','complete_two_clock_allowance',
                                 'additional_event_allowance'])
    assert F(row['original_final_excitation_weight_lower'])==F(1,2)-error
    assert not row['initial_cooling_tail_required_for_this_inference']
    assert not row['supplied_event_probability_executed_here']
    forged=deepcopy(upper);forged['result']['ground_energy_upper']=str(lo-1)
    forged['result_sha256']=contract.digest(forged['result'])
    with pytest.raises(ValueError):q.ground_energy_bracket(member,upper_receipt=forged)
    with pytest.raises(ValueError):q.observed_final_energy_to_excitation(member,(lo-2,lo-1),'1/2',upper_receipt=upper)


def test_endpoint_and_ritz_bounds_require_the_identical_member():
    member=contract.member(0);lower=endpoint.bind_member(member);upper=variational.bind_member(member)
    bracket=q.ground_energy_bracket(member,lower,upper)
    floor,ceiling=map(F,bracket['ground_energy_interval'])
    assert floor==F(lower['result']['physical_ground_lower'])
    assert floor>F(member['preparation_plan']['source']['spectral_lower'])
    assert floor<ceiling<-720
    forged=deepcopy(lower);forged['result']['physical_ground_lower']=str(ceiling)
    forged['result_sha256']=contract.digest(forged['result'])
    with pytest.raises(ValueError):q.ground_energy_bracket(member,forged,upper)
    other=contract.member(1)
    with pytest.raises(ValueError):q.ground_energy_bracket(other,lower,upper)
    row=q.observed_final_energy_to_excitation(member,(ceiling+1,ceiling+2),'1/4',lower,upper)
    assert F(row['final_excitation_interval'][1])==ceiling+F(129,64)-floor
    assert row['verified_ground_bracket']['trial_replaces_selected_preparation'] is False


def test_all_proof_and_source_paths_exist():
    for name in q.SOURCE_INPUTS:assert Path(name).is_file(),name
