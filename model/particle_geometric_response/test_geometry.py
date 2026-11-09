"""Independent finite source, covariance and full-record response controls."""
import numpy as np
import pytest
import sympy as s

from . import geometry as g


def test_flat_native_source_and_all_six_covariant_strains():
    report=g.native_control()
    assert report['complete_one_body_columns']==60
    assert report['all_six_strain_labels']==['00','11','22','01','02','12']
    assert len(set(report['all_six_strain_HS_squared']))>1
    assert all(x>0 for x in report['all_six_strain_HS_squared'])
    assert report['flat_Wilson_source_equality'] and report['flat_magnetic_source_equality']
    assert report['native_gauge_covariance'] and report['constant_frame_rotation_covariance']
    assert not report['coordinate_symbol_is_a_normalizable_quantum_preparation']


def test_six_metric_derivatives_against_independent_symbolic_inverse():
    t=s.symbols('t',real=True)
    for _,S in g.strain_basis():
        V=s.eye(3)+t*S; G=V.T*V; v=V.det()
        assert (V.inv().T.diff(t).subs(t,0)+S)==s.zeros(3)
        assert (G.inv().diff(t).subs(t,0)+2*S)==s.zeros(3)
        assert s.simplify(v.diff(t).subs(t,0)-s.trace(S))==0
        assert ((G/v).diff(t).subs(t,0)-2*S+s.trace(S)*s.eye(3))==s.zeros(3)


def test_dilation_metric_weights_and_native_regulator_term():
    v=g.coframe_data(2*s.eye(3))
    assert v['volume']==8
    assert v['inverse_frame']==s.eye(3)/2
    assert v['inverse_metric']==s.eye(3)/4
    assert v['electric_magnetic_coefficient']==s.eye(3)/2
    link=g.link_data(g.native_links()); gamma=g.density.gamma_data(); beta=gamma[3]
    W=sum((s.kronecker_product(beta,link['L'][i][i])/2 for i in range(3)),s.zeros(60))
    flat=g.native_symbol(s.eye(3))
    dilation=g.native_strain(s.eye(3))
    assert g.clean(dilation['total']+flat['total']+W)==s.zeros(60)
    assert W!=s.zeros(60)


@pytest.mark.parametrize('bad',[s.zeros(3),s.diag(-1,1,1),s.eye(2),s.eye(3)*s.I,
                                s.diag(s.Float(1.0),1,1)])
def test_coframe_requires_exact_real_positive_orientation(bad):
    with pytest.raises(ValueError): g.coframe_data(bad)


def test_full_native_rotor_work_records_keep_every_occupation_and_flux():
    r=g.rotor_record_control(); d=g.rotor_strains()
    assert r['complete_dimension']==208 and r['all_CAR_occupations']==16
    assert r['all_electric_fluxes']==list(range(-6,7))
    assert r['target_state_trace_error_upper']<s.Rational(1,10**15)
    for row in r['stress_instruments']:
        assert sum(row['original_probabilities'])==1
        assert all(p>0 for p in row['original_probabilities'])
        assert sum(a*p for a,p in zip(row['record_labels'],row['original_probabilities']))==row['exact_prepared_work_slope']
        assert row['operator_row_sum_bound']<=28 and row['Gauss_commutator_zero']
    assert d['longitudinal']-d['isotropic']==2*d['source']['HE']
    assert r['pressure_over_energy_density_prepared']!=s.Rational(1,3)
    assert r['pressure_over_energy_density_prepared']!=s.Rational(1,2)


def test_ideal_gibbs_free_energy_and_entropy_derivatives_use_same_H():
    from model.particle_field_state_response.state_response import exact_control
    d=exact_control(); b=1/64
    he,hw,hd=(np.asarray(d[k].tolist(),dtype=np.complex128) for k in ('HE','HW','HD'))
    H=he+hw+hd; K=-he-hd-2*hw
    vals,vecs=np.linalg.eigh(H)
    weights=np.exp(-b*vals);weights/=weights.sum()
    rho=(vecs*weights)@vecs.conj().T
    mean=float(np.trace(rho@K).real)
    cov=float((np.trace(rho@H@K)-np.trace(rho@H)*mean).real)
    def quantities(t):
        w=np.linalg.eigvalsh(np.exp(-t)*(he+hd)+np.exp(-2*t)*hw)
        e=np.exp(-b*w);Z=e.sum();p=e/Z
        return -np.log(Z)/b, -np.sum(p*np.log(p))
    h=1e-4; plus=quantities(h);minus=quantities(-h)
    assert abs((plus[0]-minus[0])/(2*h)-mean)<2e-7
    assert abs((plus[1]-minus[1])/(2*h)+b*b*cov)<2e-8
    report=g.rotor_record_control()
    assert abs(float(report['pressure_prepared_at_declared_volume1'])+mean/3)<2e-13


def test_flat_records_do_not_fix_metric_counterterm():
    t,xi=s.symbols('t xi',real=True)
    potential=xi*s.log(s.exp(3*t))
    assert potential.subs(t,0)==0
    assert s.simplify(potential.diff(t).subs(t,0))==3*xi
    assert s.simplify(-potential.diff(t).subs(t,0)/3)==-xi


def test_coframe_causal_rate_recovers_flat_and_bounds_rational_shear():
    from model.particle_physical_causality.cone import family
    flat=g.coframe_causal_bound(s.eye(3))
    assert flat['weighted_rate_times_spacing_upper']==family('wilson')['weighted_rate_times_spacing_upper']
    V=s.Matrix([[1,s.Rational(1,3),0],[0,1,s.Rational(1,5)],[0,0,1]])
    compact=g.compact_coframe_family()
    member=g.coframe_family_member(V,compact)
    assert member['single_coframe_rate']['weighted_rate_times_spacing_upper']==s.Rational(19183024,5625)
    assert compact['weighted_rate_times_spacing_upper']==s.Rational(275136,5)
    assert compact['physical_speed_rational_upper']==s.Rational(2201088,5)
    assert member['all_fourteen_principal_minors_nonnegative']
    assert len(member['lower_metric_principal_minors'])==len(member['upper_metric_principal_minors'])==7
    exported=g.coframe_causality_certificate()
    assert exported['native_coordinate_control_magnetic_coefficient']==s.Rational(1,7)
    assert exported['compact_family']['magnetic_scale_upper']==s.Rational(1,7)
    assert exported['flat_source_rate']['weighted_rate_times_spacing_upper']==s.Rational(5232,7)


def test_coframe_causal_mixed_star_and_collapsed_bound():
    u=s.Matrix([1,-1,0]);v=s.Matrix([1,0,-1]);L=(u*v.T+v*u.T)/2
    assert L.eigenvals()=={s.Rational(3,2):1,s.Rational(-1,2):1,s.Integer(0):1}
    assert 2*15*sum(abs(e)*n for e,n in L.eigenvals().items())==60
    # A one-cell periodic example with both links -I gives L=4I.
    assert 2*15*4==120
    pair_support={(0,0,0),(1,0,0),(0,1,0)}
    magnetic_support=pair_support|{(0,0,1)}
    distance=lambda x,y:sum(abs(a-b) for a,b in zip(x,y))
    assert len(pair_support)**2*120==1080
    assert len(magnetic_support)**2*4==64
    assert max(distance(x,y) for x in magnetic_support for y in magnetic_support)==2


def test_coframe_causal_compact_membership_and_periodic_distance_scope():
    compact=g.compact_coframe_family()
    for V in (s.eye(3)/2,2*s.eye(3)):
        assert g.coframe_family_member(V,compact)['all_fourteen_principal_minors_nonnegative']
    for V in (s.diag(s.Rational(1,3),1,1),s.diag(3,1,1)):
        with pytest.raises(ValueError):g.coframe_family_member(V,compact)
    bad=dict(compact);bad['weighted_rate_times_spacing_upper']/=2
    with pytest.raises(ValueError):g.coframe_family_member(s.eye(3),bad)
    # On a length-eight ring, sites0 and7 are one edge apart, not seven.
    assert min(abs(7+8*k) for k in (-1,0,1))==1
    assert 'quotient' in compact['physical_distance']
    assert not compact['proper_geometry_Lorentz_covariance_or_limiting_state_proved']


@pytest.mark.parametrize('kwargs',[{'min_singular':0},{'min_singular':2,'max_singular':1},
                                  {'magnetic_scale':-1},{'weight':1},{'weight':2.0}])
def test_coframe_causal_family_rejects_invalid_inputs(kwargs):
    with pytest.raises(ValueError):g.compact_coframe_family(**kwargs)
