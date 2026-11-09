"""Exact local path and full-source controls for bulk collision calibration."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import pytest

from . import bulk_action_identification as b
from model.strong_coupling_refinement import algebra


def test_all_parities_and_directions_are_exhausted_with_actual_pi_flux():
    row=b.local_path_certificate()
    assert row['parity_and_direction_cases']==24
    for case in row['cases']:
        assert len(case['length_three_paths'])==15
        assert case['positive_paths']==11 and case['negative_paths']==4
        assert case['contact_subtraction_coefficient']==25
        assert case['collar_sites']==12


def test_removing_pi_flux_changes_the_observable_fourth_coefficient():
    # All15 walks would have relative sign+ in an unsigned lattice.
    wrong=8*2*(F(15,16)+3*F(3,2)*F(1,4))
    assert wrong==33 and wrong !=25


def test_complete_raw_collision_detector_including_filled_source_and_frame():
    from model.defect_field_bridge.full_domain import branch_action
    row=b.collision_effect()
    assert row['full_local_Fock_dimension']==16
    assert row['uniform_law_coin_raw_records_per_seam']==64
    assert row['all_three_actual_frame_effect_checks']
    assert not row['full_fine_instrument_covariance_asserted']
    for law,coin in product(('C','E'),(0,1)):
        assert branch_action(15,law,coin,sites=2,flavors=2,seam=(0,1))[:3]==(8,15,1)


@pytest.mark.parametrize('c,u',[(F(1,8),F(1)),(F(3,7),F(5,11))])
def test_bulk_taylor_coefficients_from_exact_two_particle_source(c,u):
    x,y=(3,3,3),(4,3,3)
    sites={x,y,*b.neighbors(x),*b.neighbors(y)}
    first=b.pair_hop({(x,y):F(1)},sites,c,u)
    second=b.pair_hop(first,sites,c,u)
    third=b.pair_hop(second,sites,c,u)
    quadratic=sum(first.get((z,z),F())**2 for z in (x,y))
    quartic=sum(second.get((z,z),F())**2/4-first.get((z,z),F())*third.get((z,z),F())/3
                 for z in (x,y))
    assert quadratic==c*c/2
    assert quartic==-c*c*(u*u+25*c*c)/24


def test_sparse_coordinate_sector_matches_actual_car_hop_signs():
    x=(3,3,3);y=(4,3,3)
    sites=sorted({x,y,*b.neighbors(x),*b.neighbors(y)})
    index={p:i for i,p in enumerate(sites)}
    def mask_state(pair):
        a,d=pair;up=2*index[a];down=2*index[d]+1
        return (1<<up)|(1<<down),1 if up<down else -1
    pair={(x,y):F(1)}
    for _ in range(4):
        raw={mask_state(k)[0]:mask_state(k)[1]*v for k,v in pair.items()}
        next_raw=defaultdict(F)
        for mask,value in raw.items():
            for i,p in enumerate(sites):
                if (mask>>(2*i))&3==3:next_raw[mask]+=value
                for q in b.neighbors(p):
                    if q not in index:continue
                    j=index[q]
                    for spin in (0,1):
                        hop=algebra.hop(mask,2*i+spin,2*j+spin)
                        if hop:
                            target,sign=hop
                            next_raw[target]+=value*b.sign(p,q)*sign/16
        pair=b.pair_hop(pair,sites,F(1,8),F(1))
        expected={mask_state(k)[0]:mask_state(k)[1]*v for k,v in pair.items()}
        assert {k:v for k,v in next_raw.items() if v}==expected


def test_actual_eight_cube_boundary_is_not_the_bulk_coefficient():
    sites,edges=algebra.geometry((2,2,2))
    masks,twice,_=algebra.operators(sites,edges,2)
    # Same source edge in the first coordinate direction.
    coordinates=list(product(range(2),repeat=3))
    x=coordinates.index((0,0,0));y=coordinates.index((1,0,0))
    initial=(1<<(2*x))+(1<<(2*y+1));j=masks.index(initial)
    one=twice[:,j].toarray().ravel()
    three=(twice@twice@twice[:,j]).toarray().ravel()
    targets=[masks.index(3<<(2*z)) for z in (x,y)]
    inner=F(sum(int(one[k])*int(three[k]) for k in targets),16)
    assert inner==F(3,2) and 8*inner==12


def test_true_bulk_rational_records_recover_default_ratio():
    widths=[]
    for n in (10,12,14):
        t=F(1,2**n)
        p1=b.probability_enclosure(2,F(1,8),F(1),t,6)
        p2=b.probability_enclosure(2,F(1,8),F(1),2*t,6)
        out=b.identify(t,p1,p2,F(1,8),F(1))
        assert out['coupling_squared'].contains(F(1,64))
        ratio=out['contact_over_coupling_squared']
        assert ratio.contains(F(64))
        widths.append(ratio.high-ratio.low)
    assert widths[2]<widths[1]<widths[0]


def test_finite_collar_planner_controls_both_times_without_global_vacuum():
    t=F(1,2**12)
    plan=b.collar_plan(t,F(1,8),t**6)
    assert plan['both_time_settings_local_probability_error']<=t**6
    assert plan['arbitrary_exterior_state_allowed']
    assert not plan['global_vacuum_preparation_required']
    assert plan['pure_two_CAR_collar_input_required']
    assert not plan['reset_pulse_tape_executed']


def test_common_clock_scaling_leaves_exact_taylor_records_unchanged():
    c,u,t,scale=F(1,8),F(1),F(1,1024),F(3,2)
    assert b.probability_enclosure(2,c,u,t,6)==b.probability_enclosure(2,c*scale,u*scale,t/scale,6)


def test_uniform_global_attempts_preserve_denominator_and_sampling_cost():
    one=b.uniform_attempt_budget(1,F(1,100),F(1,100))
    large=b.uniform_attempt_budget(300,F(1,100),F(1,100))
    assert large['trials_per_time']==300**2*one['trials_per_time']
    assert large['score_multiplier']==300
    assert not large['other_seam_records_discarded_or_postselected']
    assert large['repair_attempts_between_preparation_and_readout']==0


def test_joined_certificate_has_only_finite_calibration_claims():
    row=b.certificate()
    assert row['known_local_calibration_preparation_consumed']
    assert not row['interacting_ground_preparation_inferred']
    assert not row['capacity_selects_action_or_clock']
    assert not row['large_native_experiment_executed']
    assert row['physical_particles']==row['physical_masses']==[]
    assert all(a['ratio_squared_enclosure'][0]<=64<=a['ratio_squared_enclosure'][1]
               for a in row['finite_time_ratio_witnesses'])
    row['configured_action']['c']._numerator=999
    assert b.certificate()['configured_action']['c']==F(1,8)


@pytest.mark.parametrize('origin,axis',[((0,0,0),0),((2,2,2),3),((2.,2,2),0)])
def test_invalid_path_domains_rejected(origin,axis):
    with pytest.raises(ValueError):b.paths(origin,axis)


def test_underresolved_collision_data_cannot_identify_contact():
    p=b.Interval(F(),F(1))
    with pytest.raises(ValueError):b.identify(F(1,10),p,p,F(1),F(8))
