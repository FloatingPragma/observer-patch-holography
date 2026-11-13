from fractions import Fraction as F
import sympy as sy
from . import source as s

def test_general_curvature_and_gauge_identities():
    assert all(s.symbolic_identities().values())

def test_full_halo_is_required_and_stencil_error_retained():
    f={(0,0,0):F(1),(1,1,0):F(2,3),(-2,0,1):F(-3,5)}
    for gamma in (F(-1),F(0),F(1),F(7,5)):
        G=s.einstein_static(f,gamma)
        assert all(not s.divergence_static(G,j) for j in range(4))
    G=s.einstein_static(f,F(0));G[1,2]={}
    assert any(s.divergence_static(G,j) for j in range(4))
    assert s.central_laplacian(f)!=s.nearest_laplacian(f)

def test_actual_lapse_shape_and_parameter_condition():
    result=s.run();rows=result['finite_geometry']['controls']
    assert [r['gamma'] for r in rows if r['pressureless_static_tensor_condition']]==['1']
    interval=list(map(F,result['native_clock_response']['far_field_response_over_small_q_approximation_interval']))
    assert F(1)<interval[0]<interval[1]<F(101,100)
    for f in (F(0),F(1,10),F(1)):
        lo,hi=s.clock_slope(f)
        assert lo<=hi<=0
    assert s.clock_slope(F(0))==(0,0)

def test_clock_derivative_of_actual_two_level_ground_energy():
    # Independently differentiate the exact parent ground eigenvalue with
    # respect to probe energy, then source energy, retaining the vacuum lapse.
    E, P, q, f, omega = sy.symbols('E P q f omega', positive=True)
    ground = omega/2 + E*sy.cosh(q) + P*sy.cosh(q*f) - sy.sqrt(
        omega**2/4 + (E*sy.sinh(q) + P*sy.sinh(q*f))**2)
    lapse = sy.diff(ground, P).subs(P, 0)
    assert sy.simplify(lapse.subs(E, 0)-sy.cosh(q*f)) == 0
    slope = sy.diff(lapse, E).subs(E, 0)/lapse.subs(E, 0)
    assert sy.simplify(slope + 2*sy.sinh(q)*sy.tanh(q*f)/omega) == 0
