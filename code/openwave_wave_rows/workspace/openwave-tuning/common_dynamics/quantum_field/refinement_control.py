"""Phase-sensitive finite source versus the analytic quantum Gaussian."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import quantum_state as qs

HERE=Path(__file__).resolve().parent


def gaussian_coefficients(data):
    a=1/(2*data['var_p'])+1j*data['cov_qp']/data['var_p']
    b=a*data['mean_p']-1j*data['mean_q']
    C=(2*math.pi*data['var_p'])**(-.25)*np.exp(-a*data['mean_p']**2/2+1j*data['mean_q']*data['mean_p'])
    return a,b,C


def gaussian_step(a,b,C,delta,beta):
    denominator=1+1j*delta*a
    C=C/np.sqrt(denominator)*np.exp(.5j*delta*b*b/denominator)*np.exp(1j*beta)
    b=b/denominator
    a=a/denominator+1j*delta*qs.A
    return a,b,C


def run():
    paths=[Path(__file__),Path(qs.__file__),HERE/'GAUSSIAN_BOUNDS.json']
    pins={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    certificates=json.loads((HERE/'GAUSSIAN_BOUNDS.json').read_text())['cases']
    rows=[];records={}
    for case in certificates:
        q=case['modulus'];eps2=F(case['epsilon_squared']);N=case['substeps'];cap=case['initial_cap']
        state,values,phases,data=qs.gaussian_cycle(q,eps2,N,cap=cap)
        a,b,C=gaussian_coefficients(data)
        p=data['epsilon']*qs.centered(q)
        initial_gaussian=np.sqrt(data['epsilon'])*C*np.exp(-a*p*p/2+b*p)
        initial_norm=float(np.linalg.norm(initial_gaussian))
        target=initial_gaussian/initial_norm
        initial_error=float(np.linalg.norm(state-target))
        for _ in range(N):
            state=qs.source_step(state,data['alpha'],data['beta'])
            a,b,C=gaussian_step(a,b,C,float(data['delta']),float(data['beta']))
        target=np.sqrt(data['epsilon'])*C*np.exp(-a*p*p/2+b*p)
        target/=np.linalg.norm(target)
        error=float(np.linalg.norm(state-target))
        bound=case['complete_finite_state_to_periodized_Gaussian_norm_upper']
        assert error<=bound+2e-11
        # No phase alignment is permitted: both exponentials carry every
        # exp(i beta) and the metaplectic square-root phase.
        rows.append({'modulus':q,'epsilon_squared':str(eps2),'substeps':N,
                     'initial_state_error':initial_error,'phase_sensitive_final_norm_error':error,
                     'analytic_norm_bound':bound,'phase_alignment_performed':False,
                     'electric_TV':float(np.sum(abs(abs(state)**2-abs(target)**2))/2),
                     'angle_TV':float(np.sum(abs(abs(np.fft.fft(state,norm='ortho'))**2-abs(np.fft.fft(target,norm='ortho'))**2))/2)})
        records[f'{q}:{N}']={'actual_source':qs.complete_records(state),'Gaussian_reference':qs.complete_records(target)}
    assert pins=={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'REFINEMENT_RECORDS.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
    result={'schema':'oph.local12.quantum_gaussian_refinement_control.v1','checks':{
      'actual_finite_cosine_source_executed':True,'metaplectic_Gaussian_evolution_executed':True,
      'every_scalar_phase_retained':True,'no_output_phase_alignment':True,
      'analytic_state_bound_controls_all_cases':True,'complete_final_field_records':True,
      'source_bytes_unchanged':True},'cases':rows,
      'scope':'Floating independent controls of the separately outward-certified theorem. The Gaussian reference is an optical cycle oscillator, not an interacting photon vacuum. Native q16 apparatus is supplied separately; larger alphabets are ideal source executions with no inherited native pulse precision.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'REFINEMENT_CONTROL.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'cases':rows}),flush=True)


if __name__=='__main__':run()
