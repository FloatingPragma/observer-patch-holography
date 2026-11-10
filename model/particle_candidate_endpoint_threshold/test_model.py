from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import model as m


def test_complete_geodesic_geometry():
    faces=m.plaquettes()
    assert len(faces)==6 and all(len(face)==4 for face in faces)
    for length,count in ((1,24),(2,48),(3,48)):
        geo=m.geometry(length)
        assert geo['path_count']==count
        for path in m.paths(length):
            assert len(set(path))==length+1
            assert sum(a!=b for a,b in zip(path[0],path[-1]))==length
            assert all(sum(x!=y for x,y in zip(a,b))==1 for a,b in zip(path,path[1:]))
            incident={face for face in faces if any(m.edge(a,b) in face for a,b in zip(path,path[1:]))}
            assert len(incident)==length+1
    for bad in (0,4,-1,True,1.0):
        with pytest.raises(ValueError):m.paths(bad)


@pytest.mark.parametrize('length',[1,2,3])
def test_exact_full_chain_charge_and_flux(length):
    row=m.gauss_chain_control(length)
    assert row['control_dimension']==4*3**length
    assert row['Gauss_generators']==length+1
    assert row['onsite_electric_commutator_center']==6+F(63*length,4)


def test_literal_endpoint_columns_and_full_fock_bound():
    row=m.endpoint_columns_control()
    assert row['literal_endpoint_column_cases']==240
    assert row['each_endpoint_full_Fock_norm_upper']==1
    assert row['summed_off_path_norm_upper']==2


@pytest.mark.parametrize('length',[1,2,3])
def test_path_hopping_blocks_have_orthogonal_inputs_and_outputs(length):
    # Independent finite link model tests the support argument, including a
    # hopping matrix with constant-to-character transitions and complex entries.
    I=s.eye(3);P=s.diag(1,0,0);Q=s.diag(0,1,0)
    L=s.zeros(3);L[1,0]=1
    T=s.Matrix([[0,1,s.I],[1,0,2],[-s.I,2,3]])
    A=s.kronecker_product(*([L]*length)) if length>1 else L
    Ppath=s.kronecker_product(*([P]*length)) if length>1 else P
    Qpath=s.kronecker_product(*([Q]*length)) if length>1 else Q
    left=[];right=[]
    for i in range(length):
        factors=[I]*length;factors[i]=T
        Ti=s.kronecker_product(*factors) if length>1 else T
        left.append(Ti*A);right.append(A*Ti)
    for i in range(length):
        assert Qpath*left[i]==s.zeros(3**length)
        assert right[i]*Ppath==s.zeros(3**length)
        for j in range(length):
            assert left[i].H*right[j]==s.zeros(3**length)
            assert left[i]*right[j].H==s.zeros(3**length)
            if i!=j:
                assert left[i].H*left[j]==s.zeros(3**length)
                assert right[i]*right[j].H==s.zeros(3**length)


def test_original_dressed_complete_record_bounds():
    member,bare,odd,w,epsilon=m.original_inputs()
    assert F(bare['result']['original_success_weight_lower'])==w
    assert epsilon==F(odd['result']['weak_plus_e_defect_count_upper'])/48
    result=m.bind_member(member)['result']
    assert result['complete_unfiltered_probe_outcomes_with_setting']==483
    assert result['physical_mass_table_GeV']==[]
    intervals=(('55/32','1337/32'),('239/32','2161/32'),('263/32','3145/32'))
    for row,interval,margin in zip(result['separations'],intervals,(F(4,25),F(2,25),F(7,250))):
        assert tuple(row['reported_transfer_interval'])==interval
        assert F(row['original_filtered_recorded_weight_strict_lower'])>margin
        assert F(row['uniform_three_setting_joint_simple_strict_lower'])==margin/3
        assert not row['weight_is_conditioned_on_success']
        assert row['failure_and_reference_records_retained']
        assert F(row['complete_two_clock_event_loss_upper'])<=F(1,8192)
        assert F(row['filter_compiler']['distinguished_original_event_probability_error_upper'])<=F(6,65536)


def test_primary_binding_rejects_capacity_clock_and_mass_substitution():
    member=m.contract.member(0)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    changed=deepcopy(member);changed['ids']['clock_id']='substituted'
    with pytest.raises(ValueError):m.bind_member(changed)
    receipt=m.bind_member(member);changed=deepcopy(receipt)
    changed['result']['physical_mass_table_GeV']=[{'species':'electron','mass':0.000511}]
    changed['result_sha256']=m.contract.digest(changed['result'])
    with pytest.raises(ValueError):m.verify_bound(member,changed)


def test_pair_data_do_not_determine_endpoint_mass_without_potential_convention():
    mass=F(7,5);shift=F(3,7);potential=[F(3,10),F(2,5),F(9,20)]
    energies=[2*mass+v for v in potential]
    assert energies==[2*(mass+shift)+(v-2*shift) for v in potential]
    assert mass!=mass+shift


def test_sqrt_rounds_down():
    for x in (F(1,3),F(7,11),F(2),F(0)):
        lo=m.rational_sqrt_lower(x)
        assert lo*lo<=x<(lo+F(1,10**6))**2


def test_increment_meter_requires_inside_temporal_placement():
    control=m.controlled_increment_meter_control()
    assert control['temporal_interference_preserved']
    assert not control['moving_meter_outside_temporal_conjugation_is_equivalent']
    assert control['all_source_charge_meter_records']==390625
