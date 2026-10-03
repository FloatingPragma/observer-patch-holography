"""Conditional top criticality with arbitrary positive two-anchor capacities.

No mass targets or old source-root/RSCC closures are read. Only independent
alpha calibration and the declared energy-unit/top transport specification
are parsed. The one-loop band is an analytic monotonicity consequence;
two-loop rows are explicit controls, not a certified continuum enclosure.
"""
import argparse
from decimal import Decimal
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, root

REPO = Path(__file__).resolve().parents[3]
CAL = REPO / 'code/particles/calibration'
SPEC = CAL / 'top_positive_capacity_band_spec.json'
ALPHA = REPO / 'code/P_derivation/codata_2022_alpha_fixture.json'
OUT = REPO / 'code/particles/runs/calibration/top_positive_capacity_band.json'
sys.path.insert(0, str(CAL))
sys.path.insert(0, str(REPO / 'code/P_derivation'))
import sm_two_loop_rge_engine as rge
from paper_math import PaperMathContext

REVIEWED_SPEC_SHA256 = 'fda35fc4b193454fbaf7738f747b02ba29972d13feefd3a69becd5701192977d'
PIN_PATHS = [SPEC, ALPHA, CAL/'sm_two_loop_rge_engine.py',
             REPO/'code/P_derivation/paper_math.py', Path(__file__).resolve()]

pi = math.pi
k = 1 / (16*pi*pi)
B = np.array([41/6, -19/6, -7.])
C = np.array([17/12, 9/4, 8.])


def require(ok, message):
    if not ok:
        raise ValueError(message)


def reviewed_spec():
    require(hashlib.sha256(SPEC.read_bytes()).hexdigest() == REVIEWED_SPEC_SHA256,
            'unreviewed specification, input or promotion change')
    return json.loads(SPEC.read_text())


class GaugeOnlyContext(PaperMathContext):
    def _derive_stage5_integer_vectors(self):
        return {}

    def forbidden(self, *args, **kwargs):
        raise ValueError('Flavor-dependent closure or mass template is forbidden')

    diagonal_quark_masses = forbidden
    charged_lepton_masses = forbidden
    alpha_external_from_d10 = forbidden
    structured_thomson_running = forbidden
    structured_thomson_running_asymptotic = forbidden
    solve_closure = forbidden


def source_packet(spec):
    fixture = json.loads(ALPHA.read_text())['inverse_fine_structure_constant']
    require(fixture == {'value': spec['inverse_alpha'],
                       'standard_uncertainty': spec['inverse_alpha_standard_uncertainty']},
            'independent alpha fixture differs')
    d = spec['D10']
    ctx = GaugeOnlyContext(precision=d['precision'], su2_cutoff=d['su2_cutoff'],
                           su3_cutoff=d['su3_cutoff'])
    pp = ctx.p_from_inverse_alpha(Decimal(fixture['value']))
    packet = ctx.build_d10_from_p(pp)
    P, E = float(pp), float(spec['E_star_display_GeV'])
    return dict(P=P, E_star=E, alpha_U=float(packet.alpha_u),
                mu_U=math.exp(-2*pi)*P**(1/6)*E, E_cell=E/math.sqrt(P),
                v=float(packet.v)*E, mz=float(packet.mz_run)*E,
                low_g_squared=[4*pi*float(packet.alpha_y_mz),
                               4*pi*float(packet.alpha2_mz),4*pi*float(packet.alpha3_mz)],
                alpha_inverse=fixture['value'])


def qcd_proxy(m,alpha):
    a = alpha/pi
    return m*(1+4*a/3+(13.4434-1.0414*5)*a*a+(190.595-26.655*5+.653*25)*a**3)


def build_payload():
    spec = reviewed_spec()
    packet = source_packet(spec)
    t0 = math.log(packet['mz'])
    low = np.array(packet['low_g_squared'])
    def den(t): return 1/low - 2*k*B*(t-t0)
    def gauge(t): return 1/den(t)
    def crit(t):
        a,b,_ = gauge(t)
        return ((a*a+2*a*b+3*b*b)/16)**.25
    def evolve_factor(t,s): return math.exp(float(np.sum((-C/B)*np.log(den(t)/den(s)))))
    def yukawa(t,tb,eps=2e-12):
        val,err = quad(lambda s:evolve_factor(t,s),t,tb,epsabs=eps,epsrel=eps)
        z = evolve_factor(t,tb)/crit(tb)**2+9*k*val
        assert z > 0
        return z**(-.5), 9*k*err
    def mass_1loop(tb):
        def fn(m):return yukawa(math.log(m),tb)[0]*packet['v']/math.sqrt(2)-m
        m=brentq(fn,50.,packet['v'],xtol=1e-11)
        y,err=yukawa(math.log(m),tb)
        # An independent differential-equation replay checks the quadrature.
        def beta(t,y):return k*y*(4.5*y*y-float(C@gauge(t)))
        sol=solve_ivp(beta,(tb,math.log(50)),[crit(tb)],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        assert sol.success
        mm=brentq(lambda mu:sol.sol(math.log(mu))[0]*packet['v']/math.sqrt(2)-mu,50.,packet['v'],xtol=1e-11)
        a,b,c=gauge(tb);ap,bp,cp=2*k*B*gauge(tb)**2
        X=a*a+2*a*b+3*b*b
        Fprime_over_F=(2*a*ap+2*ap*b+2*a*bp+6*b*bp)/(4*X)
        sensitivity_seed=crit(tb)*(Fprime_over_F-k*(4.5*crit(tb)**2-float(C@gauge(tb))))
        return dict(mu_boundary=math.exp(tb),y_boundary=crit(tb),mass_coordinate=m,
                    QCD_pole_proxy=qcd_proxy(m,gauge(math.log(m))[2]/(4*pi)),
                    self_scale_residual=fn(m),quadrature_error_estimate_in_inverse_y_squared=err,
                    independent_ODE_mass=mm,ODE_difference=mm-m,
                    boundary_sensitivity_seed=sensitivity_seed)

    tU,tE=math.log(packet['mu_U']),math.log(packet['E_cell'])
    capacity_rows=[]
    # Fixed fractions show continuity and numerical direction. They do not
    # prove the continuum band; that result follows from the analytic bounds.
    for u in spec['anchors']['control_fractions']:
        row=mass_1loop((1-u)*tU+u*tE)
        row['weight_E_over_sum']=u
        capacity_rows.append(row)

    domain_gauge=[gauge(t) for t in [math.log(50),tU,tE]]
    proof_checks={
        'positive_inverse_gauge_minimum':float(min(np.min(den(t)) for t in [math.log(50),tE])),
        'c_over_b_on_anchor_endpoints':[float(g[2]/g[1]) for g in domain_gauge[1:]],
        'c_over_a_on_domain_endpoints':[float(g[2]/g[0]) for g in [domain_gauge[0],domain_gauge[-1]]],
        'positive_sensitivity_coefficient':4-(11/12+9*math.sqrt(3)/8),
        'negative_Cprime_coefficient':697/36-28,
        'uniform_lower_y_at_50':(3*gauge(tE)[1]**2/16)**.25,
        'required_y_at_50':50*math.sqrt(2)/packet['v'],
        'uniform_upper_y_at_mu_U':math.sqrt(2*float(C@gauge(tU))/9),
        'required_y_at_mu_U':packet['mu_U']*math.sqrt(2)/packet['v'],
    }
    assert proof_checks['positive_inverse_gauge_minimum']>0
    assert min(proof_checks['c_over_b_on_anchor_endpoints'])>.5
    assert min(proof_checks['c_over_a_on_domain_endpoints'])>.5
    assert proof_checks['uniform_lower_y_at_50']>proof_checks['required_y_at_50']
    assert proof_checks['uniform_upper_y_at_mu_U']<proof_checks['required_y_at_mu_U']

    # Fully coupled two-loop control with low gauge anchors held fixed.
    gl=np.sqrt(low*np.array([5/3,1.,1.]))
    def weak_yt(g):
        f=lambda y:rge.beta_2loop(*g,y,0.)[4]
        y=brentq(f,.05,2.,xtol=1e-14)
        a,b,c=np.array(g)**2
        coeff=[30*k*k,-6*k+k*k*(-32*c-8*a/5),
               k*k*(-9*b*b/4-171*a*a/100+63*a*b/10),
               k*3*(2*b*b+(b+3*a/5)**2)/8+k*k*(305*b**3/16-3411*a**3/2000-289*b*b*a/80-1677*b*a*a/400)]
        rr=np.roots(coeff)
        assert max(abs(r.imag) for r in rr)<1e-9
        positive=sorted(r.real for r in rr if r.real>0)
        assert abs(y*y-positive[0])<1e-10
        return y
    coupled=[]
    for u in spec['anchors']['control_fractions']:
        tb=(1-u)*tU+u*tE
        def evolve(g,tol):
            sol=solve_ivp(lambda t,y:rge.beta_2loop(*y),(tb,math.log(50)),[*g,weak_yt(g),0.],method='DOP853',rtol=tol,atol=tol/100,dense_output=True)
            assert sol.success
            return sol
        guess=np.sqrt(gauge(tb)*np.array([5/3,1.,1.]))
        shots=[]
        for tol in spec['transport']['two_loop_rtol']:
            def residual(g):return evolve(g,tol).sol(t0)[:3]-gl
            shot=root(residual,guess,tol=1e-10)
            assert shot.success and np.max(np.abs(residual(shot.x)))<1e-9
            guess=shot.x
            sol=evolve(shot.x,tol)
            fn=lambda m:sol.sol(math.log(m))[3]*packet['v']/math.sqrt(2)-m
            m=brentq(fn,50.,packet['v'],xtol=1e-11)
            _,_,g3,y,lam=sol.sol(math.log(m))
            shots.append(dict(rtol=tol,mass_coordinate=m,QCD_pole_proxy=qcd_proxy(m,g3*g3/(4*pi)),
                              Higgs_tree_proxy=packet['v']*math.sqrt(2*lam),
                              gauge_anchor_max_residual=float(np.max(abs(residual(shot.x)))),
                              criticality_residual=rge.beta_2loop(*shot.x,weak_yt(shot.x),0.)[4],
                              self_scale_residual=fn(m)))
        coupled.append(dict(weight_E_over_sum=u,mu_boundary=math.exp(tb),rows=shots,
                            tolerance_mass_difference=abs(shots[0]['mass_coordinate']-shots[1]['mass_coordinate'])))
    return dict(schema='oph.top_positive_capacity_band.v1', specification=spec,
                scope=spec['scope'],packet=packet,
                proof_checks=proof_checks,one_loop_monotone_band_endpoints=[capacity_rows[0]['mass_coordinate'],capacity_rows[-1]['mass_coordinate']],
                one_loop_controls=capacity_rows,coupled_two_loop_controls=coupled,
                numerical_interpretation=spec['uncertainty'],
                source_pins={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in PIN_PATHS})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    result=build_payload()
    out=args.output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'output':str(out),'one_loop_band':result['one_loop_monotone_band_endpoints'],
                      'coupled_controls':[row['rows'][-1]['mass_coordinate'] for row in result['coupled_two_loop_controls']]}))
