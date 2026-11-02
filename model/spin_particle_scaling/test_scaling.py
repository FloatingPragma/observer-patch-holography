from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import algebra as a, bounds as b, experiment as e, packet


def test_full_affine_and_onsite_checks():
    c=a.certificate()
    assert c['affine_cases_exact']==5184
    assert c['full_onsite_CAR_frame_cases']==192
    assert c['ordered_spatial_stabilizer_mod2_count']==32
    assert c['transverse_residual_zone_centers_in_pi_units']==[[0,0,1],[0,1,0]]


def test_translation_order_countercontrol():
    n=np.array([1,0,0])
    assert np.array_equal(a.order_translation((1,0,0))@n,-n)
    assert np.array_equal(a.order_translation((0,1,0))@n,n)
    assert not np.array_equal(a.frame((1,0,0)),np.eye(3))


@pytest.mark.parametrize('x',[(True,0,0),(0.,0,0),(0,0),[0,0,0]])
def test_inexact_coordinate_rejected(x):
    with pytest.raises(ValueError):a.frame(x)


def test_full_finite_spectral_matrix_and_spatial_characters():
    c=e.check()
    assert c['spectral_matrix_scalar_at_each_cluster']
    assert c['circular_channels_orthogonal']
    assert c['every_finite_energy_eigenspace_checked']
    assert c['significant_positive_transverse_response_checked']
    assert c['one_spin_U1_alone_does_not_force_equal_weights']


def test_corrupted_finite_eigenspace_rejected(monkeypatch):
    original=e.np.linalg.eigh
    def corrupted(matrix):
        values,vectors=original(matrix)
        vectors[:,1]*=2
        return values,vectors
    monkeypatch.setattr(e.np.linalg,'eigh',corrupted)
    with pytest.raises(ArithmeticError,match='finite eigensolver controls failed'):e.check()


def test_lattice_integral_and_diagonal_constants():
    c=b.certificate()
    assert F(c['lattice_sum_constant_exact'])==F(6561,1120)<6
    for n in (2,3,16,128):
        row=b.diagonal(n)
        assert row['boundary_power_of_two']<=-19*n**3
        E=row['edges']
        assert E==6*128*n**3*(256*n**3+1)**2
        assert F(row['source_trace_error_bound'])==F(24,E**3)+F(24*n*n,E**2)+F(18*n*n,E)


def test_normalized_error_requires_stronger_diagonal_than_unscaled_error():
    # The naive radius n² gives E~n6 and absolute error n^-4, which
    # fails after the independently derived n4 normalization cost.
    naive=[F(n**4*n*n,n**6) for n in (4,8,16)]
    assert naive==[1,1,1]
    rows=[b.diagonal(n) for n in (4,8,16)]
    scaled=[F(v['source_normalized_error_coefficient_K2_over_c']) for v in rows]
    assert scaled[1]<scaled[0]/7 and scaled[2]<scaled[1]/7


@pytest.mark.parametrize('n',[True,0,-1,2.,'2'])
def test_inexact_scale_rejected(n):
    with pytest.raises(ValueError):b.diagonal(n)


def test_normalization_fails_closed():
    assert b.normalized_ratio_error(F(1,100),F(1,100),F(1,10))==F(2,9)
    with pytest.raises(ValueError):b.normalized_ratio_error(F(1,100),F(1,10),F(1,10))
    with pytest.raises(ValueError):b.normalized_ratio_error(True,F(1,100),F(1,10))


def test_source_bindings_are_narrow_and_explicit():
    pins=packet.source_bindings()
    assert 15<=len(pins)<=25
    assert all(not p.startswith(('claims/','tracking/','model/completion/')) for p in pins)
    assert not set(packet.PRIORS).intersection(pins)
    assert all(path in pins for path,_ in packet.SNAPSHOTS.values())
    assert all(row['checkpoint']=='ba312dbf0' for row in packet.inherited().values())


def test_historical_context_needs_no_git_history_and_rejects_changed_snapshot(monkeypatch,tmp_path):
    def no_git(*args,**kwargs):raise AssertionError('runtime git history is forbidden')
    monkeypatch.setattr(packet.subprocess,'check_output',no_git)
    assert packet.inherited()
    bad=tmp_path/'bad.json';bad.write_text('{}')
    snapshots=dict(packet.SNAPSHOTS)
    snapshots[packet.PRIORS[0]]=(str(bad),snapshots[packet.PRIORS[0]][1])
    monkeypatch.setattr(packet,'SNAPSHOTS',snapshots)
    with pytest.raises(ValueError):packet.inherited()


@pytest.mark.parametrize('mutation',['scope','numeric_bool','energy','custody'])
def test_receipt_rejects_scope_or_type_laundering(mutation):
    value=packet.produce();bad=deepcopy(value)
    if mutation=='scope':bad['theorem']['sharp_particle_pole_or_Lorentz_helicity_proved']=True
    if mutation=='numeric_bool':bad['theorem']['two_orthogonal_equal_weight_transverse_channels']=1
    if mutation=='energy':bad['direct_checks']['scaling']['examples'][0]['spectral_weight_lower']='1'
    if mutation=='custody':bad['dependency_custody_pass']=True
    with pytest.raises(ValueError):packet.validate(bad)
