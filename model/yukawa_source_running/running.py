"""Conditional minimal-Standard-Model one-loop matrix running.

All scales and gauge boundary values are dimensionless declared controls.
The source does not select this continuum, matching map or running interval.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

LOOP = 16*np.pi*np.pi
GAUGE_B = np.array([41/10, -19/6, -7.])
GAUGE_C = np.array([[17/20,9/4,8.],[1/4,9/4,8.],[9/4,9/4,0.]])


def pack(gauge, yukawas):
    return np.concatenate([np.asarray(gauge,complex),np.asarray(yukawas,complex).reshape(-1)])


def unpack(values):
    return values[:3], values[3:].reshape(3,3,3)


def beta(_t, values):
    g, yy = unpack(values)
    grams = yy@yy.conj().transpose(0,2,1)
    trace = np.trace(3*grams[0]+3*grams[1]+grams[2]).real
    alpha = trace-GAUGE_C@(g*g)
    left = np.asarray([1.5*(grams[0]-grams[1]),1.5*(grams[1]-grams[0]),1.5*grams[2]])
    left += alpha[:,None,None]*np.eye(3)
    return pack(GAUGE_B*g*g*g,left@yy)/LOOP


def evolve(yukawas, *, end=-40., gauge=(.5,.5,.5), rtol=2e-11, atol=2e-13):
    initial = pack(gauge,yukawas)
    solution = solve_ivp(beta,(0.,end),initial,method='DOP853',rtol=rtol,atol=atol,
                         t_eval=np.linspace(0.,end,81))
    if not solution.success or not np.isfinite(solution.y).all():
        raise ArithmeticError('finite one-loop flow failed')
    gs = solution.y[:3].T
    yy = solution.y[3:].T.reshape(-1,3,3,3)
    sv = np.linalg.svd(yy,compute_uv=False)[...,::-1]
    gauge_exact = 1/np.sqrt(1/np.asarray(gauge)**2-2*GAUGE_B*end/LOOP)
    if max(abs(gs[-1]-gauge_exact)) > 1e-9:
        raise ArithmeticError('analytic gauge running control')
    ratio = sv[:,2,1:]/sv[:,2,:-1]
    if np.max(np.diff(ratio,axis=0)) > 1e-7:
        raise ArithmeticError('charged-lepton IR contraction violated')
    reverse = solve_ivp(beta,(end,0.),solution.y[:,-1],method='DOP853',rtol=rtol,atol=atol)
    error = float(np.max(abs(reverse.y[:,-1]-initial)))
    if error > 2e-9:
        raise ArithmeticError('finite-flow injectivity control')
    return {'log_scale_ratio':end,'gauge_initial':list(gauge),
            'gauge_final':gs[-1].real.tolist(),'gauge_exact_error':float(max(abs(gs[-1]-gauge_exact))),
            'initial_singular_values':sv[0].tolist(),'final_singular_values':sv[-1].tolist(),
            'initial_charged_lepton_adjacent_ratios':ratio[0].tolist(),
            'final_charged_lepton_adjacent_ratios':ratio[-1].tolist(),
            'inverse_flow_max_error':error,'function_evaluations':solution.nfev,
            'source_selected_gauge_or_scale':False}, {'times':solution.t,'gauge':gs,'yukawas':yy,'singular_values':sv}


def algebra_controls():
    rng = np.random.default_rng(1701)
    unitary = []
    for _ in range(5):
        q,_ = np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        unitary.append(q)
    yy = (rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3)))/7
    left = (unitary[0],unitary[0],unitary[3])
    right = (unitary[1],unitary[2],unitary[4])
    changed = np.asarray([l@y@r.conj().T for l,y,r in zip(left,yy,right)])
    _,b = unpack(beta(0,pack((.5,.5,.5),yy)))
    _,c = unpack(beta(0,pack((.5,.5,.5),changed)))
    covariance = max(np.max(abs(cj-l@bj@r.conj().T)) for cj,l,bj,r in zip(c,left,b,right))
    zero = np.zeros((3,3,3),complex)
    if np.any(unpack(beta(0,pack((.5,.5,.5),zero)))[1]):
        raise ArithmeticError('radiative zero source')
    rank_one = zero.copy();rank_one[:,0,0] = [.5,.25,.125]
    _,rb = unpack(beta(0,pack((.5,.5,.5),rank_one)))
    rb[:,0,0] = 0
    if np.any(rb) or covariance > 1e-14:
        raise ArithmeticError('rank or flavour covariance')
    # In finite IR time, even the scalar attractive fixed-point model has
    # strictly positive sensitivity to every positive boundary coefficient.
    a,b,time = 4.5,1.,40.
    fixed = b/a;decay=np.exp(-2*b*time/LOOP)
    fixed_rows=[]
    for initial in (0.,.125,.5,1.,2.):
        squared = initial*initial
        final = 0. if not squared else fixed/(1+(fixed/squared-1)*decay)
        sensitivity = decay*fixed*fixed/(squared*(1-decay)+fixed*decay)**2
        fixed_rows.append({'initial':initial,'final':float(np.sqrt(final)),
                           'squared_boundary_sensitivity':float(sensitivity)})
    return {'flavour_covariance_max':float(covariance),'zero_Yukawa_beta_exact':True,
            'rank_one_diagonal_support_preserved_exact':True,
            'scalar_fixed_point_control':{'a':a,'b':b,'IR_time':time,'rows':fixed_rows}}
