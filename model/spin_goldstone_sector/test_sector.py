from copy import deepcopy
from fractions import Fraction as F
import math
import pytest
from . import algebra,bounds,experiment,infrared,packet


@pytest.fixture(scope='module')
def value():return packet._produce()


def test_exact_integral_certifies_order_and_retains_singularity():
    r=infrared.certificate()
    assert F(r['integral_upper'])<F(383,1000)<F(2,5)
    assert F(r['origin_cell_upper'])==F(33,64**2)
    assert r['nonorigin_cells']==64**3-1
    assert F(r['regular_upper'])+F(r['origin_cell_upper'])==F(r['integral_upper'])


def test_coarse_mesh_does_not_falsely_certify_order():
    assert infrared.certificate(32)['integral_below_two_fifths'] is False


def test_exact_trigonometric_enclosures():
    lo,hi=infrared.pi_interval()
    assert float(lo)==pytest.approx(math.pi) and float(hi)==pytest.approx(math.pi)
    assert lo<hi and hi-lo<F(1,10**40)
    for i in range(1,64):
        assert infrared.cosine_upper(lo*i/64,1<<40)/(1<<40)>=math.cos(math.pi*i/64)-1e-15


def test_shell_count_and_origin_ball_control():
    for r in range(1,100):
        assert (2*r+1)**3-(2*r-1)**3==24*r*r+2<=26*r*r
    assert 6*F(22,7)*F(7,4)<=33 and F(7,4)**2>3


def test_exact_current_ward_and_free_spin_countercontrol():
    r=algebra.certificate()
    assert len(r['exact_spin_current_and_moment_axes'])==3
    assert r['free_spin_counterexample']['gap_above_ground_space']=='1'


def test_complete_finite_cube_normalization_control():
    r=experiment.check()
    assert r['full_spin_dimension']==256 and r['nonsingular_susceptibilities_checked']==7


def test_window_scaling_is_not_an_inferred_pole():
    a=bounds.window(F(1,64),F(1,1000),F(3))
    b=bounds.window(F(1,128),F(1,1000),F(3))
    for key in ('filter_energy_lower','filter_energy_upper','filtered_weight_lower','conditional_mean_energy_upper'):
        assert F(a[key])==2*F(b[key])
    assert F(a['first_moment_upper'])==4*F(b['first_moment_upper'])
    assert F(a['first_moment_upper'])/F(a['filtered_weight_lower'])==F(a['conditional_mean_energy_upper'])


@pytest.mark.parametrize('p,m,j',[(0,F(1,10),1),(True,F(1,10),1),(F(1,2),F(3,4),1),(2,F(1,10),1),(F(1,2),F(1,10),float('nan'))])
def test_invalid_window_inputs_rejected(p,m,j):
    with pytest.raises(ValueError):bounds.window(p,m,j)


def test_selection_defect_vanishes_on_diagonal():
    old=None
    for e in (12,96,768):
        r=bounds.selection_defect(F(1,e),3,e)
        x=F(r['unit_norm_local_ground_defect_upper'])
        assert x==F(9,e)
        if old is not None:assert x<old
        old=x


def test_receipt_has_honest_inherited_custody(value):
    packet.validate(value)
    assert value['new_checks_pass'] is True
    assert value['all_pass'] is False
    assert value['prior_custody'][packet.PRIORS[0]]['historical_all_pass'] is False


@pytest.mark.parametrize('key,wrong',[('physical_particle_species_identified',True),('sharp_one_particle_pole_proved',True),('field_selected_ground_itself_proved_pure',True),('internal_spin_identified_with_spatial_spin',0),('ordered_ground_sector_constructed',1)])
def test_scope_or_type_laundering_rejected(value,key,wrong):
    v=deepcopy(value);v['theorem'][key]=wrong
    with pytest.raises(ValueError):packet.validate(v)


def test_omitting_singular_cell_rejected(value):
    v=deepcopy(value);v['direct_checks']['infrared']['origin_cell_upper']='0'
    with pytest.raises(ValueError):packet.validate(v)


def test_inherited_custody_cannot_be_promoted(value):
    v=deepcopy(value);v['all_pass']=v['dependency_custody_pass']=True
    with pytest.raises(ValueError):packet.validate(v)
