from copy import deepcopy
from fractions import Fraction as F
import pytest
from . import algebra,bounds,controls,packet


def test_exact_spinor_fourier_and_two_hessian_calculations():
    c=algebra.certificate()
    assert all(c['checks'].values())
    assert c['generic_natural_vector_bound_is_insufficient'] is True


def test_independent_spectral_integrals():
    d=controls.diagnostics()
    assert d['hessian_relative_error']<1e-6
    assert d['Stieltjes_max_error']<1e-8


def test_reject_nonfinite_quadrature(monkeypatch):
    monkeypatch.setattr(controls,'quad',lambda *a,**k:(float('nan'),0))
    with pytest.raises(ArithmeticError,match='nonfinite'):controls.diagnostics()


def test_exact_geometric_and_error_constants():
    assert bounds.domain()['exponential_norm_coefficient']==2**25
    assert bounds.domain()['maximum_rho']==F(1,8)
    small=bounds.estimates(F(1,64),F(1,4))
    half=bounds.estimates(F(1,128),F(1,4))
    assert half['norm_squared_upper']*8==small['norm_squared_upper']
    assert half['relative_entropy_upper']*64==small['relative_entropy_upper']
    assert bounds.neumann_budget(F(1,2**36))<F(1,4)


def test_symbolic_joint_mesh_never_materializes_extreme_power():
    c=bounds.joint_mesh(F(1,8),F(1,2),True)
    assert c['epsilon']==str(F(1,8)**7)
    assert c['mesh_float_materialized'] is False
    assert c['linear']['mesh_float_value_supplied'] is False
    assert c['hessian']['mesh_float_materialized'] is False
    assert c['leading_coefficient_error_after_rho6_division']=='1/8'


@pytest.mark.parametrize('value',[True,0.1,float('nan'),'1/8'])
def test_reject_inexact_or_boolean_parameters(value):
    with pytest.raises(ValueError):bounds.estimates(value,F(1,4))


@pytest.mark.parametrize('args',[(0,F(1,4)),(F(1,4),F(1,4)),(F(1,64),1)])
def test_reject_outside_proved_domain(args):
    with pytest.raises(ValueError):bounds.estimates(*args)


def test_generic_natural_vectors_do_not_gain_rho6_without_new_estimate():
    r=F(1,64);lam=r**3;mass=lam*(1-lam)
    chi=mass**2*(1/lam+1/(1-lam))
    assert chi==mass
    assert chi/r**6>2**17


@pytest.fixture(scope='module')
def receipt():return packet.build()


def test_current_complete_receipt(receipt):
    packet.validate(receipt)
    assert receipt['theorem']['source_valleys_retained']==8
    assert receipt['custody']['full_repository_replay_performed'] is False


@pytest.mark.parametrize('flag',[
    'uniform_covariance_gap_assumed','arbitrary_L_R_schedule_o_R4_asserted',
    'sharp_boosted_lattice_algebra_constructed','area_Newton_or_stationarity_derived',
    'interacting_Einstein_endpoint_closed','nonzero_U_entropy_theorem',
    'native_entropy_readout_executed','measured_particle_species_identified',
])
def test_reject_missing_premise_promotions(receipt,flag):
    bad=deepcopy(receipt);bad['theorem'][flag]=True
    with pytest.raises(ArithmeticError):packet.validate(bad)


@pytest.mark.parametrize('where',['theorem','bounds','algebra','root','source'])
def test_reject_forged_exact_or_scientific_payload(receipt,where):
    bad=deepcopy(receipt)
    if where=='theorem':bad['theorem']['source_valleys_retained']=8.0
    elif where=='bounds':bad['bounds']['constant_spinor_hessian_pi4_coefficient']='32/35'
    elif where=='algebra':bad['algebra']['checks']['Stieltjes_second_derivative']=1
    elif where=='root':bad['physical_mass_identified']=True
    else:bad['source_bindings'][next(iter(bad['source_bindings']))]='0'*64
    with pytest.raises(ArithmeticError):packet.validate(bad)


def test_reject_straddled_source_snapshot(monkeypatch):
    original=packet.source_bindings();calls=0
    def changing():
        nonlocal calls
        calls+=1;copy=dict(original)
        if calls>1:copy[next(iter(copy))]='0'*64
        return copy
    monkeypatch.setattr(packet,'source_bindings',changing)
    with pytest.raises(ArithmeticError,match='source changed'):packet.build()
