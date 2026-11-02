from fractions import Fraction as F
import copy
import numpy as np
import pytest
from . import algebra, bounds, controls, packet


def test_exact_noncommuting_full_CAR_and_Sherman_Morrison():
    result=algebra.certificate()
    assert all(result['checks'].values())
    assert result['checks']['noncommuting_covariance_perturbation']


def test_full_Fock_density_Hessian_and_finite_amplitude():
    result=controls.finite_checks()
    assert result['full_Fock_Hessian_and_entropy_checked']
    assert result['opposite_resolvent_sign_mutation_rejected']


def test_natural_endpoint_and_unnatural_negative_control():
    assert controls.endpoint_checks()['arbitrary_fixed_norm_probe_extension_rejected']


@pytest.mark.parametrize('e,b,B',[(F(1,2),F(1,2),6),(1,0,1),(F(1,100),F(9,10),3)])
def test_exact_symbolic_modulus_budget(e,b,B):
    result=bounds.modulus(e,b,B)
    assert F(result['hessian_error_upper'])<=F(e)/2
    assert result['mesh_float_materialized'] is False
    assert result['mesh_strictly_positive'] is True
    assert result['mesh_power_exponent']==result['series_terms']+1


def test_positive_actual_shell_constants_and_finite_spacing():
    result=bounds.shell_bounds(F(1,2**18));m=F(result['mass_lower'])
    assert F(result['continuum_hessian_lower'])==4*m*m>0
    assert F(result['finite_hessian_lower_with_extra_mesh'])==m*m
    assert F(result['positive_mesh_numerator_over_B_squared'])>0


def test_zero_direction_and_nonzero_amplitude():
    assert all(x==0 for x in bounds.entropy_bounds(0,F(1,4)).values())
    result=bounds.entropy_bounds(F(1,8),F(-1,4))
    assert result['relative_entropy_lower']>0
    assert result==bounds.entropy_bounds(F(1,8),F(1,4))


@pytest.mark.parametrize('call',[lambda:bounds.modulus(True),lambda:bounds.modulus(0),
    lambda:bounds.modulus(F(1,2),1),lambda:bounds.modulus(F(1,2),F(1,2),0),
    lambda:bounds.entropy_bounds(F(1,3),0),lambda:bounds.entropy_bounds(F(1,8),F(3,4)),
    lambda:bounds.shell_bounds(F(1,2**17)),lambda:bounds.regularity(False),
    lambda:bounds.shell_bounds(0.1)])
def test_invalid_inputs_rejected(call):
    with pytest.raises(ValueError):call()


@pytest.mark.parametrize('matrix',[np.array([[np.nan]]),np.array([[0.0]]),np.array([[1,1],[0,1]])])
def test_invalid_numerical_log_inputs_rejected(matrix):
    with pytest.raises(ValueError):controls.log_matrix(matrix)


@pytest.fixture(scope='module')
def certificate():return packet.build()


def test_current_certificate(certificate):packet.validate(certificate)


@pytest.mark.parametrize('path,value',[
    (('theorem','source_U'),0.0),
    (('theorem','source_valleys_retained'),1),
    (('theorem','actual_Araki_relative_entropy_limit'),1),
    (('theorem','area_dynamics_or_stationarity_derived'),True),
    (('theorem','uniform_covariance_gap_assumed'),True),
    (('bounds','modulus','mesh_strictly_positive'),False),
    (('algebra','checks','noncommuting_covariance_perturbation'),False),
    (('custody','predecessor_full_repository_receipts_replayed'),True)])
def test_false_claim_or_type_alias_rejected(certificate,path,value):
    result=copy.deepcopy(certificate);cursor=result
    for key in path[:-1]:cursor=cursor[key]
    cursor[path[-1]]=value
    with pytest.raises(ArithmeticError):packet.validate(result)


def test_unknown_root_claim_rejected(certificate):
    result=copy.deepcopy(certificate);result['Einstein_equation']=True
    with pytest.raises(ArithmeticError,match='root fields'):packet.validate(result)


def test_bad_source_digest_rejected(certificate):
    result=copy.deepcopy(certificate);result['source_bindings'][next(iter(result['source_bindings']))]='0'*64
    with pytest.raises(ArithmeticError):packet.validate(result)


def test_source_race_rejected(monkeypatch):
    original=packet.source_bindings;calls=0
    def changed():
        nonlocal calls
        calls+=1;result=original()
        if calls>1:result[next(iter(result))]='0'*64
        return result
    monkeypatch.setattr(packet,'source_bindings',changed)
    with pytest.raises(ArithmeticError,match='source changed'):packet.build()


def test_nonfinite_quadrature_rejected(monkeypatch):
    monkeypatch.setattr(controls,'quad',lambda *a,**k:(float('nan'),0.0))
    with pytest.raises(ArithmeticError,match='nonfinite'):
        controls.finite_checks()


def test_nonfinite_eigendecomposition_rejected(monkeypatch):
    monkeypatch.setattr(controls.np.linalg,'eigh',lambda a:(np.array([float('nan')]),np.eye(1)))
    with pytest.raises(ArithmeticError,match='nonfinite'):
        controls.log_matrix(np.eye(1))
