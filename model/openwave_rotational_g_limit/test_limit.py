import numpy as np
from .exact import certificate
from .limit import bounds,stable_kernel,oscillator,operators
from model.openwave_angular_torque.source import kernel
def test_exact_all_mode_polynomials():assert certificate()['rows']==40
def test_stable_kernel_retains_both_signs():
    data=oscillator();m=7.;b=np.array([.02,-.01,.03]);ext=np.array([.01,.04,-.03]);ops=operators(data,m,1.)
    D=ops[0]+sum((v*x for v,x in zip(b,ops[1])),np.zeros_like(ops[0]))+sum((v*x for v,x in zip(ext,ops[2])),np.zeros_like(ops[0]))
    K=stable_kernel(data,m,b,ext)[0]
    assert np.linalg.norm(K+m*np.eye(len(K))-kernel(D)[0]-ops[3])<2e-12
def test_uniform_rational_route_improves():
    values=[bounds(m)for m in (64,256,1024,4096,16384)]
    assert all(x['self_field_contraction']<1 and x['axial_stationary_branch_contraction']<1 for x in values)
    assert all(a['quasistatic_precession_g_error_from_two']>b['quasistatic_precession_g_error_from_two']for a,b in zip(values,values[1:]))
    assert values[-1]['quasistatic_precession_g_error_from_two']<.004863
def test_wrong_spin_coefficient_fails():
    data=oscillator();K,D,E,A=stable_kernel(data,10.,np.zeros(3),np.zeros(3))
    e=operators(data,10.,1.)[2][2];vac=np.zeros(len(A));vac[0]=1
    assert abs(np.vdot(vac,(A@e+e@A)@vac)-.5)>.49
