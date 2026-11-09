"""Outward source normalization and sharp necessary-matching controls."""
from fractions import Fraction as F
from copy import deepcopy
import json
import mpmath as mp
import pytest

from . import measured_mass_join as m


@pytest.mark.parametrize('q',[0,1,8,F(1,27),F(3,7),F(10**50+1,10**70)])
def test_integer_cube_enclosures(q):
    for bits in (8,32,100):
        lo,hi=m.cube_bounds(q,bits)
        assert lo**3<=q<=hi**3 and hi-lo<=F(1,1<<bits)


@pytest.mark.parametrize('bad',[True,1.2,float('nan')])
def test_exact_input_requirement(bad):
    with pytest.raises(ValueError):m.cube_bounds(bad)


def test_shared_denominator_requires_cross_ratio_distortion():
    r=[[1,1],[2,2],[4,4]];ref=[[1,1],[F(9,5),F(9,5)],[F(22,5),F(22,5)]]
    d=m.necessary_budgets(r,ref)
    assert d['minimal_component_multiplier_spread_interval']==['11/9','11/9']
    assert d['necessary_relative_matching_error_lower']=='1/10'
    assert d['minimal_additive_E_over_source_light_mass_interval']==['8/31','8/31']
    w=m.additive_witness([1,2,4],ref)
    x=list(map(F,w['shifted_singular_values']));eta=F(w['eta'])
    assert x[1]/x[0]==F(9,5) and x[2]/x[0]==F(22,5)
    assert max(abs(a-b) for a,b in zip(x,(1,2,4)))==eta


def test_common_scale_cannot_change_an_excluded_ratio():
    r=[[1,1],[2,2],[4,4]];ref=[[1,1],[3,3],[6,6]]
    d=m.necessary_budgets(r,ref)
    assert F(d['minimal_component_multiplier_spread_interval'][0])>1
    assert d['actual_source_correction_upper_bound_available'] is False


def test_exact_compatible_windows_need_no_correction():
    r=[[1,1],[2,2],[4,4]];ref=[[1,1],[F(19,10),F(21,10)],[F(39,10),F(41,10)]]
    d=m.necessary_budgets(r,ref)
    assert d['minimal_component_multiplier_spread_interval']==['1','1']
    assert d['minimal_additive_E_over_source_light_mass_interval']==['0','0']


def test_interval_lower_bound_is_safe_on_all_corners():
    r=[[1,1],[F(19,10),F(21,10)],[F(39,10),F(41,10)]]
    ref=[[1,1],[F(22,10),F(23,10)],[F(42,10),F(43,10)]]
    envelope=m.necessary_budgets(r,ref)
    for a in r[1]:
        for b in r[2]:
            point=m.necessary_budgets([[1,1],[a,a],[b,b]],ref)
            for key in ('minimal_component_multiplier_spread_interval','minimal_additive_E_over_source_light_mass_interval'):
                lo,hi=map(F,envelope[key]);p=F(point[key][0])
                assert lo<=p<=hi


@pytest.mark.parametrize('r',[
    [[2,2],[3,3],[4,4]],[[1,1],[0,0],[4,4]],[[1,1],[3,2],[4,4]],
    [[1,1],[2,3],[3,4]],[[1,1],[2,2]],[[1,1],[True,True],[4,4]]])
def test_ratios_require_order_and_common_denominator(r):
    with pytest.raises(ValueError):m.necessary_budgets(r,[[1,1],[2,2],[4,4]])


@pytest.fixture(scope='module')
def report():
    return m.certificate()


def test_all_four_source_branches_and_no_invented_physical_bound(report):
    assert [r['branch'] for r in report['branches']]==list(m.BRANCHES)
    json.dumps(report,allow_nan=False)
    assert report['computed_interacting_correction'] is False
    assert report['proved_physical_correction_upper_bounds']==report['physical_particle_map']==[]
    for row in report['branches']:
        assert row['physical_mass_predictions']==[]
        assert row['complete_memory_inverse_branches']==2
        assert row['matching_budgets']['actual_source_correction_upper_bound_available'] is False
        assert row['additional_exact_balance_hypothesis']['assumed_for_reference_comparison'] is False


def test_observed_window_needs_differ_by_branch(report):
    values={r['branch']:r['matching_budgets'] for r in report['branches']}
    assert values['measured_endpoint']['a_nonzero_correction_is_necessary_in_this_constant_matching_class'] is False
    for name in ('forward','mssm_beta4_proxy','sm_color_beta7'):
        assert values[name]['a_nonzero_correction_is_necessary_in_this_constant_matching_class'] is True


def test_higgs_units_agree_with_existing_point_producer(report):
    from model.mcpr_native_memory_bridge import construction
    for row in report['branches'][:2]:
        target=construction.target(row['branch'])['masses_over_v']
        with mp.workdps(90):
            for pair,x in zip(row['masses_over_declared_Higgs_amplitude'],target):
                lo,hi=map(F,pair)
                assert mp.mpf(lo.numerator)/lo.denominator<=x<=mp.mpf(hi.numerator)/hi.denominator


def test_higgs_units_overlap_owning_closure_intervals(report):
    retained=json.loads((m.ROOT/'model/pixel_mass_closure/results.json').read_text())
    for row,old in zip(report['branches'][2:],retained['closures']):
        assert row['branch']==old['design']
        for pair,box in zip(row['masses_over_declared_Higgs_amplitude'],old['implicit_and_mass_enclosures']['masses_over_v']):
            lo,hi=map(F,pair)
            assert max(lo,F(box['lo']))<=min(hi,F(box['hi']))


def test_common_gain_cancels_in_ratios_but_not_higgs_units(report):
    for row in report['branches']:
        masses=[tuple(map(F,r)) for r in row['masses_over_declared_Higgs_amplitude']]
        for pair,ratio in zip(masses[1:],row['original_action_ratios'][1:]):
            lo,hi=map(F,ratio)
            assert max(pair[0]/masses[0][1],lo)<=min(pair[1]/masses[0][0],hi)


def test_calibration_cannot_promote_source_scope():
    capacity=m.source.certificate();ref=json.loads((m.ROOT/m.source.REFERENCE).read_text())
    bad=deepcopy(capacity);bad['source_parameters_fitted_by_this_certificate']=True
    with pytest.raises(ValueError):m.summarize(bad,ref)
    bad=deepcopy(capacity);bad['branches'].pop()
    with pytest.raises(ValueError):m.summarize(bad,ref)
