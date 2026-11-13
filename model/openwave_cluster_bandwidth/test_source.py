from fractions import Fraction as F
from itertools import product
import math
from .source import local_completion,angle,run
from .lambda_source import source
from .meson import native_two_cell_check, degree_ancestry, data


def test_native_sources_are_actual_parent_restrictions():
    assert source(2)['D']==64
    for strange in(False,True):
        assert native_two_cell_check(strange)['intertwiner_residual_diagnostic']<1e-9
    from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
    from model.openwave_hadron_refinement_flow.positivity import inputs
    p,_,_=inputs();dU,dn=degree_coefficients();ancestry=degree_ancestry()
    c=F(3,4)*(p['Jw']+p['eta'])
    base=2*p['mbar']+2*p['t']+c
    edge=base+2*p['A']/9+8*p['U']/3+F(3,2)*p['Js']+3*p['spin_gradient']
    a,b,_,_,_=data(False)
    assert a-base==2*dn and b-edge==2*dn+8*dU/3
    assert F(ancestry['bulk_colocated_pair_shift_source_exact'])==2*dn


def test_local_completion_is_bounded_and_cluster_additive():
    rows=list(product(range(3),repeat=2))
    for counts in product(rows,repeat=3):
        value=local_completion(counts)
        assert 0<=value<=sum(sum(row)for row in counts)
    for a,b in product(rows,repeat=2):
        assert local_completion((a,(0,0),(0,0),b))==min(a)+min(b)
    assert local_completion(((1,0),(0,1),(1,0)))==2


def test_energy_shell_and_original_bandwidth_control():
    r=run();s=r['energy_shell_bracket']
    import json
    from model.openwave_pure_su2_source.custody import equivalent
    assert equivalent(json.loads(json.dumps(r,allow_nan=False)),r)
    assert F(s['below_margin_GeV_exact'])>0 and F(s['above_margin_GeV_exact'])>0
    assert r['locality_only_control']['insufficient_bandwidth']
    assert not r['native_weak_scope']['complete_outgoing_weak_decay_or_lifetime_proved']
    assert not r['all_parent_two_cell_mass_fits_are_bulk_poles']


def test_arctangent_interval_bounds_independent_small_series():
    for t in(F(0),F(1,6),F(1,5),F(1,2)):
        lo,hi=angle(t);assert lo<=hi
        assert abs(float((lo+hi)/2)-2*math.atan(float(t)))<1e-15
        assert hi-lo<F(1,10**20)
