import numpy as np
from scipy.integrate import solve_ivp
from . import source as s


def test_constraint_propagation_and_spatial_trace_are_exact_identities():
    assert s.constraint_identity()['symbolic_residuals_exact']==['0','0']


def test_actual_initial_data_satisfy_both_constraints():
    for kind in ('capacity','lambda'):
        for k in (.1,10.,100.):
            initial=s.initial(k,kind,1e-7)
            values=s.components(np.log(1e-7),initial,k,kind)
            assert abs(values[10])<2e-14
            assert abs(initial[1]+initial[6]-1.5)<2e-14
            assert values[9]<0  # The growing potential already evolves in the mixed era.


def test_metric_is_dynamical_and_scalar_control_has_distinct_transfer():
    result,_=s.integrate(10.,rtol=2e-11)
    other,_=s.integrate(10.,kind='lambda',rtol=2e-11)
    assert result['maximum_sampled_normalized_constraint_residual']<2e-6
    assert abs(result['records'][-1]['metric_potential_N_transfer'])>1e-3
    assert abs(result['records'][-1]['metric_potential_transfer']-other['records'][-1]['metric_potential_transfer'])>1e-3


def test_einstein_de_sitter_exact_growing_solution(monkeypatch):
    # An independent analytic cosmology checks the transformed coordinates.
    # In pure dust Phi=1, U_m=2/3, Delta_m=-2*k^2*a/3 exactly.
    def dust_background(kind='capacity'):
        return lambda N: (1.,0.,0.,0.,0.,0.,-1.5,np.exp(-3*N))
    monkeypatch.setattr(s,'background',dust_background)
    k=3.; N0=-4.
    state=np.array([-2*k*k*np.exp(N0)/3,2/3,0.,0.,0.,0.,1.])
    solution=solve_ivp(lambda N,v:s.rhs(N,v,k,'capacity'),(N0,0.),state,
                       method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    assert solution.success
    for N in np.linspace(N0,0.,31):
        v=solution.sol(N)
        assert abs(v[6]-1.)<2e-11
        assert abs(v[1]-2/3)<2e-11
        assert abs(v[0]+2*k*k*np.exp(N)/3)<2e-10
        assert abs(s.rhs(N,v,k,'capacity')[0]/v[0]-1)<2e-10
        assert abs(s.components(N,v,k,'capacity')[10])<3e-10
