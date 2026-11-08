import copy
from fractions import Fraction as F
import json
import subprocess
import sys

import numpy as np
import pytest

from . import algebra,bounds,packet,rotations


def test_whole_projection_domain_has_same_clock_incompatibility():
    for edges in (1,2,12,37):
        for ratio in (4*edges,8*edges,100*edges):
            for p in (F(1),F(1,16),F(1,256)):
                j=F(3,7);c=j*ratio;u=j*ratio*ratio;q=1/(j*p)
                row=bounds.common_clock(c,u,edges,p,q)
                assert row['ballistic_spin_clock_rate']==1
                assert row['rescaled_charge_gap_lower']==ratio*(F(ratio,2)-edges)/p
                assert row['rescaled_charge_gap_lower']>=4*edges*edges/p
                assert row['complete_charged_weight_upper']<=p/F(4*edges*edges)


def test_chemical_shift_cannot_close_both_conjugate_branches():
    for shift in range(-10,11):
        row=bounds.chemical_pair(F(3),F(4),F(shift),F(2))
        assert row['conjugate_sum']==14
        assert row['same_reference_addition_energy']+row['same_reference_removal_energy']==14
        assert row['same_reference_stable_against_these_two_sectors']==(-4<=shift<=3)
        assert not row['negative_excitation_reinterpreted_as_ground_gap']


def test_complete_energy_weight_bound_on_independent_full_dimer_state():
    creators=[np.asarray(algebra.creation(i),dtype=float)for i in range(4)]
    number=sum(a@a.T for a in creators)
    h=np.zeros((16,16));u=1.;c=.125
    for a,b in ((0,2),(1,3)):
        h+=c/2*(creators[a]@creators[b].T+creators[b]@creators[a].T)
    h+=u*((creators[0]@creators[0].T)@(creators[1]@creators[1].T)+(creators[2]@creators[2].T)@(creators[3]@creators[3].T))-u*number/2
    eig,vec=np.linalg.eigh(h);ground=vec[:,0];relative=h-eig[0]*np.eye(16)
    charged=np.diag([int(mask.bit_count()!=2)for mask in range(16)])
    gap_lower=.5-c
    assert np.linalg.eigvalsh(relative-gap_lower*charged)[0]>-1e-14
    for amplitude in (.05,.2,.5):
        created=creators[0]@ground
        weight=np.vdot(created,created).real
        assert abs(weight-.5)<1e-14
        # Original preparation probability and energy stay unnormalized.
        v=amplitude*created
        probability=np.vdot(v,v).real
        energy=np.vdot(v,relative@v).real
        assert probability<=energy/gap_lower


def test_internal_register_and_source_spin_grading_are_distinct():
    row=algebra.dimer()
    assert sum(row['complete_source_eigenvalues'].values())==16
    assert row['source_odd_creator_compression_to_single_occupancy_zero']
    assert all(x['addition_weight']==x['removal_weight']=='1/2'for x in row['complete_addition_removal_weights'])
    assert row['finite_odd_and_spin_excitations_both_nonzero']


def test_actual_spatial_cube_action_and_spin_four_countercontrol():
    row=rotations.certificate()
    assert row['proper_spatial_rotations']==24
    assert row['source_vertices']==64 and row['all_source_bonds']==192
    assert row['H_J_bonds_and_staggered_selection_preserved']
    assert row['local_source_charge_creators_are_even_and_disjointly_commuting']
    assert row['continuous_rotation_derivative_of_countercontrol']!='0'
    assert not row['cubic_singlet_alone_proves_spin_zero']
    assert not row['relativistic_spin_statistics_theorem_instantiated']


@pytest.mark.parametrize('args',[(1,3,1,F(1,2),1),(True,4,1,F(1,2),1),(1,4,0,F(1,2),1),(1,4,1,0,1),(1,4,1,F(1,2),0),(1,4,1,F(2),1)])
def test_unproved_or_invalid_domain_is_refused(args):
    with pytest.raises(ValueError):bounds.common_clock(*args)


@pytest.fixture(scope='module')
def result():return packet.build()


def test_actual_finite_gauge_countercontrol_and_internal_currents(result):
    row=result['actual_internal_current_and_complete_gauge_control']
    assert len(row['internal_source_components'])==15
    assert all(x['occupation_mask'].bit_count()%2==0 for x in row['internal_source_components'])
    control=row['same_singlet_gauge_transport_countercontrol']
    assert control['reference_dimension']==3 and control['endpoint_singlet_checks']==48
    assert control['same_input_and_detector_rotation_probabilities']==['16/25','144/169']
    assert F(144,169)-F(16,25)==F(896,4225)


def test_retained_receipt_is_current(result):
    assert result==json.loads((packet.HERE/'results.json').read_text())


def test_full_ground_local_records_and_aligned_reference_are_separate(result,monkeypatch):
    local=result['local_interacting_matter_history']
    charged=result['centered_original_charged_bands']
    aligned=result['aligned_hole_position_dispersion_and_continuum']
    assert local['same_centered_thermodynamic_ground_subsequence']
    assert not local['global_raw_archive_equality']
    assert charged['cubic_readable_parameter_band']['each_original_preparation_weight']=='1/2'
    assert charged['actual_dimer']['actual_composite_clock']=='B=-3h'
    assert not aligned['full_interacting_ground_reference_selected']
    assert not aligned['physical_particle_labels'] and not aligned['physical_masses']
    composite=result['source_generated_composite_CAR_response']
    assert composite['independent_CAR_flavours_added']==0
    assert not composite['all_volume_residual_discarded']
    assert composite['finite_full_ground_controls'][0]['residual_squared_trace']=='0'
    assert F(composite['finite_full_ground_controls'][1]['residual_squared_trace'])>0
    for proof in ('LOCAL_MATTER_READOUT.md','CHARGED_SPECTRAL_BANDS.md','ALIGNED_HOLE_DISPERSION.md','COMPOSITE_CAR.md'):
        assert 'model/source_family_compatibility/'+proof in result['source_bindings']
    monkeypatch.setattr(packet,'build',lambda:result)
    for name,key in (('local_interacting_matter_history','global_raw_archive_equality'),
                     ('centered_original_charged_bands','charged_pole_mass_or_Lorentz_covariance_derived'),
                     ('source_generated_composite_CAR_response','all_volume_residual_discarded'),
                     ('aligned_hole_position_dispersion_and_continuum','full_interacting_ground_reference_selected')):
        changed=copy.deepcopy(result);changed[name][key]=True
        with pytest.raises(ValueError):packet.validate(changed)


def test_full_interacting_budget_retains_leakage_and_no_physical_assignment(result,monkeypatch):
    row=result['full_interacting_defect_budget']
    assert row['finite_source']['actual_interacting_time_control']['defect_expectation']=='9/25'
    assert row['full_Hubbard_propagation_in_same_repair_domain'] is True
    assert row['interacting_relativistic_particle_limit_proved'] is False
    assert 'model/source_family_compatibility/DEFECT_BUDGET.md' in result['source_bindings']
    monkeypatch.setattr(packet,'build',lambda:result)
    for key in ('capacity_action_or_physical_clock_selected','interacting_relativistic_particle_limit_proved'):
        value=copy.deepcopy(result);value['full_interacting_defect_budget'][key]=True
        with pytest.raises(ValueError):packet.validate(value)


@pytest.mark.parametrize('key',[k for k,v in packet.claims().items()if not v])
def test_scope_cannot_be_promoted(result,monkeypatch,key):
    value=copy.deepcopy(result);value['claims'][key]=True
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(value)


@pytest.mark.parametrize('kind',['gap','weight','chemical','gauge','type','hash'])
def test_scientific_or_custody_mutation_is_rejected(result,monkeypatch,kind):
    value=copy.deepcopy(result)
    if kind=='gap':value['same_clock_controls'][0]['rescaled_charge_gap_lower']='0'
    if kind=='weight':value['same_clock_controls'][0]['complete_charged_weight_upper']='1'
    if kind=='chemical':value['chemical_reference_controls'][0]['conjugate_sum']='0'
    if kind=='gauge':value['actual_internal_current_and_complete_gauge_control']['same_singlet_gauge_transport_countercontrol']['same_input_and_detector_rotation_probabilities'][0]='0'
    if kind=='type':value['all_pass']=1
    if kind=='hash':value['source_bindings'][next(iter(value['source_bindings']))]='0'*64
    monkeypatch.setattr(packet,'build',lambda:result)
    with pytest.raises(ValueError):packet.validate(value)


def test_actual_check_PATH_is_consumed(tmp_path):
    missing=tmp_path/'missing.json'
    run=subprocess.run([sys.executable,'-m','model.source_family_compatibility','--check',str(missing)],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode and str(missing)in run.stderr


def test_clean_runtime_source_imports_are_bound():
    script='''
from pathlib import Path
import sys
from model.source_family_compatibility import packet
packet.build();bindings=packet.bindings();missing=[]
for module in list(sys.modules.values()):
    source=getattr(module,'__file__',None)
    if not source:continue
    p=Path(source).resolve()
    if p.suffix=='.py'and p.is_relative_to(packet.ROOT)and p.relative_to(packet.ROOT).as_posix()not in bindings:missing.append(str(p))
assert not missing,missing
'''
    run=subprocess.run([sys.executable,'-W','error','-c',script],cwd=packet.ROOT,capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
