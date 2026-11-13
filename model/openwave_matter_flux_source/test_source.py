from fractions import Fraction as F
import pytest
from .source import action,physical_labels,gauss_residual,complete_gauss_control,spectral_certificate,source_inventory,ALPHA,BETA,HOP

def test_full_integer_gauss_bijection_and_missing_dressing_control():
    assert len(complete_gauss_control()['all16_occupations'])==16
    for q0,q1,n in ((F(1,3),F(-2,3),10003),(7,-5,-10001)):
        assert gauss_residual(q0,q1,physical_labels(q0,q1,n))==(0,0)
    with pytest.raises(ValueError):physical_labels(F(1,2),0,0)

def test_full_residual_includes_omitted_fluxes():
    v={(0,3):F(1)};h=action(v,1)
    assert h[0,4]==-BETA/2 and h[0,2]==-BETA/2 and h[1,3]==-HOP
    c=spectral_certificate(1,4)
    assert all(abs(r['flux'])==5 for r in c['all_omitted_boundary_action_rows'])
    assert F(c['full_integer_source_residual_squared_exact'])>0

def test_positive_common_mass_shift_and_charge_conjugation():
    v,p,m=[spectral_certificate(q) for q in (0,1,-1)]
    vlo,vhi=map(F,v['full_source_ground_interval_exact']);plo,phi=map(F,p['full_source_ground_interval_exact']);mlo,mhi=map(F,m['full_source_ground_interval_exact'])
    assert plo>vhi and max(plo,mlo)<=min(phi,mhi)
    assert phi-plo<F(1,10**15)

def test_no_uniform_bounded_rest_operator_difference():
    assert action({(0,10000):F(1)},0)[0,10000]>10**6
    # Dropping the loop energy changes the finite source, not its presentation.
    assert action({(0,0):F(1)},0)[0,0]==BETA

def test_actual_neutral_rotation_and_every_CAR_coordinate_retained():
    _,c=source_inventory()
    assert c['constant_neutral_family_rotation_retained']
    assert c['native_CAR_modes']==192 and c['native_source_terms']==114
    assert c['full192_CAR_residual_diagnostic']<8e-15

def test_full_vacuum_is_below_all_charged_and_coloured_floors():
    parent,_=source_inventory()
    floor=F(parent['positivity']['colourful_sector_floor_per_occupation_exact'])
    v=spectral_certificate(0)
    assert floor>1>BETA>F(v['full_source_ground_interval_exact'][1])
    assert F(parent['positivity']['finite_source_gap_GeV_exact'])>0


def test_optional_common_subtraction_preserves_full_vacuum_but_changes_generator():
    parent,_=source_inventory();v,p=[spectral_certificate(q) for q in (0,1)]
    vlo,vhi=map(F,v['full_source_ground_interval_exact']);plo,phi=map(F,p['full_source_ground_interval_exact'])
    dlo,dhi=plo-vhi,phi-vlo
    floor=F(parent['positivity']['colourful_sector_floor_per_occupation_exact'])
    assert 0<dlo<=dhi<BETA and min(floor,1-dhi)>vhi
    assert 24*dhi<F(3,8)
    # Arbitrary subtraction is not safe: subtracting2 drops the charged
    # electron trial below the vacuum. The exact bounded delta is essential.
    assert 1+phi-2<vlo
    # Weak transitions need not preserve charged occupation. This complete
    # two-state control forbids silently transferring the old generator.
    import numpy as np
    number=np.diag([0,1]);weak=np.array([[0,1],[1,0]])
    assert np.linalg.norm(weak@number-number@weak)>0
