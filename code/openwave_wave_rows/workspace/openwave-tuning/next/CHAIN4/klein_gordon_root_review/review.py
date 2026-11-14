"""Independent source, native record and refinement review of the KG candidate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PRODUCER = ROOT/'openwave-tuning/next/CHAIN3/klein_gordon_completion'
sys.path.insert(0, str(PRODUCER))
import native_spectroscopy as ns


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def custody():
    path = PRODUCER/'MANIFEST.json'
    assert sha(path) == 'de7ae79937061bdd1686c779d7660158416d73b5505c819dc3c2ef9d4371a05e'
    d = json.loads(path.read_text())
    pins = {str(path.relative_to(ROOT)): sha(path)}
    for key in ('files_sha256', 'source_sha256'):
        for name, expected in d[key].items():
            p = PRODUCER/name if key == 'files_sha256' else ROOT/name
            assert sha(p) == expected, name
            pins[str(p.relative_to(ROOT))] = expected
    return pins


def independent_algebra():
    perm = np.arange(144)
    for direction, pairs in enumerate(ns.a.common.INTERNAL_MATCHINGS):
        for v, w in pairs:
            x, y = 12*v+2*direction+1, 12*w+2*direction+1
            perm[x], perm[y] = y, x
    signs = np.zeros((3,144), dtype=np.int64)
    for axis, pairs in enumerate(ns.a.spatial.links()):
        for v, w in pairs:
            x, y = 12*v+11, 12*w+11
            perm[x], perm[y] = y, x
            signs[axis,x], signs[axis,y] = 1, -1
    I = np.eye(144, dtype=np.int64)
    R6 = np.kron(np.eye(12,dtype=np.int64), np.ones((12,12),dtype=np.int64)-6*np.eye(12,dtype=np.int64))
    W6 = R6[perm]
    C12 = W6 + W6.T
    P = I.copy()
    for root in (4,6,8,-12):
        P = P @ (C12-root*I)
    denominator = 4608
    assert np.array_equal(P@P,denominator*P)
    assert np.trace(P) == 26*denominator
    assert np.array_equal(R6@P,denominator//12*np.ones((144,144),dtype=np.int64)-6*P)
    # Project the literal displacement derivative, retaining exact numerators.
    H = [P @ np.diag(s) @ P for s in signs]
    Hden = denominator**2
    one = np.ones(144,dtype=np.int64)
    g = [h@one for h in H]
    for i in range(3):
        assert np.array_equal(144*H[i],g[i][:,None]+g[i][None,:])
        assert not np.any(R6@H[i]+H[i]@R6)
        for j in range(3):
            dot = sum(int(x)*int(y) for x,y in zip(g[i],g[j]))
            assert dot == (3*Hden**2 if i == j else 0)
    return perm, signs, [x.astype(float)/(12*Hden) for x in g]


def apply_walk(state,perm,phase,mu):
    # Direct local sum/reflection and actual seam permutation, no producer W.
    reflection = (state.reshape(12,12).sum(axis=1)[:,None]/6-state.reshape(12,12)).reshape(-1)
    return phase*(math.cos(mu)*reflection-1j*math.sin(mu)*state)[perm]


def main():
    pins = custody()
    perm, signs, gradients = independent_algebra()
    e = np.ones(144)/12
    result = json.loads((PRODUCER/'RESULT.json').read_text())
    worst = 0.
    cases = []
    for mode in ((1,0,0),(0,1,1),(2,1,1)):
        p = 2*np.pi*np.array(mode)
        K = F(7 if max(mode)==1 else 13)
        for mass in (F(0),F(3,4),F(1)):
            energy = math.sqrt(float(mass)**2+float(p@p)/48)
            gradient = sum(p[i]*gradients[i] for i in range(3))
            limit = math.cos(energy)*e+math.sin(energy)/energy*(1j*gradient-1j*float(mass)*e)
            for size in (64,256,1024):
                h = F(1,size)
                phase = np.exp(1j*(float(h)*p@signs))
                state = e.astype(complex)
                for _ in range(size): state=apply_walk(state,perm,phase,float(h*mass))
                amplitude = np.vdot(e,state)
                target = math.cos(energy)-1j*float(mass)*math.sin(energy)/energy
                error = float(np.linalg.norm(state-limit))
                C = F(5,4)*K*(K+abs(mass)+K/4+abs(mass))+F(13,16)*(K+abs(mass))**2+(K/4+abs(mass))**2/2
                bound = F(5,2)*h*K+size*h*h*C
                assert error <= float(bound)+1e-11
                assert abs(np.vdot(e,limit)-target)<1e-12
                assert abs(np.vdot(state,state)-1)<2e-12
                worst=max(worst,error)
                cases.append({'mode':mode,'mass':str(mass),'size':size,'full_state_error':error,
                              'exact_upper':str(bound),'return':[amplitude.real,amplitude.imag]})
    # All recorded terminal laws retain both apparatus families and failures.
    settings = 0
    for row in json.loads((PRODUCER/'COMPLETE_RECORDS.json').read_text()):
        size, steps, mass = row['size'], row['steps'], float(F(row['mass_exact']))
        phase=np.exp(1j*(2*np.pi*np.array(row['mode'])/size@signs))
        state=e.astype(complex)
        for _ in range(steps):state=apply_walk(state,perm,phase,mass/size)
        c=np.vdot(e,state)
        for setting,law in row['laws'].items():
            assert all(x['probability']>=-1e-14 for x in law['matter']+law['clock']+law['joint'])
            assert abs(sum(x['probability'] for x in law['joint'])+sum(law['preparation_export_failures'])-1)<3e-11
            observed=sum((1-2*x['record']['bit'])*x['probability'] for x in law['matter'] if x['outcome']=='native_pointer')
            assert abs(observed-(c.real if setting[0]=='X' else c.imag))<3e-11
            clock=sum((1-2*x['record']['bit'])*x['probability'] for x in law['clock'] if x['outcome']=='native_pointer')
            assert abs(clock-math.cos(steps/(8*size)))<3e-11
            for x in law['joint']:
                assert abs(x['probability']-law['matter'][x['matter']]['probability']*law['clock'][x['clock']]['probability'])<1e-15
            settings+=1
    exact_bounds=[]
    for row in result['formal_complete_native_family']:
        h=F(row['h_exact']);n=row['steps'];K=F(7);m=F(row['mass_exact'])
        C=F(5,4)*K*(K+abs(m)+K/4+abs(m))+F(13,16)*(K+abs(m))**2+(K/4+abs(m))**2/2
        limit=F(5,2)*h*K+n*h*h*C
        reset=row['reset'];timeout=(reset['blank_carriers']+reset['ready_carriers'])*F(reset['base_failure_exact'])**reset['blocks_per_register']
        native=F(row['delta_preparation_exact'])+n*F(row['delta_controlled_step_exact'])+F(row['delta_clock_and_readout_exact'])+timeout
        assert timeout==F(reset['timeout_bound_exact']) and timeout<=h*h
        assert native+limit==F(row['native_vs_continuum_uncapped_exact'])
        assert row['source_time_in_tau0_exact']=='1' and row['large_word_executed'] is False
        exact_bounds.append(native+limit)
    assert all(x>y for x,y in zip(exact_bounds,exact_bounds[1:])) and exact_bounds[-1]<F(1,40)
    # Actual native scalar control is essential against the inactive reference.
    word=ns.controlled_event('I',(0,),ns.a.angle(0,0,F(2,7)),1)
    native=ns.b.complete_run(word,2,np.eye(4),96);ns.b.require_complete(native)
    expected=np.diag([1,1,np.exp(2j/7),np.exp(2j/7)])
    scalar_error=float(np.max(abs(native['success']-expected)))
    assert scalar_error<3e-12 and np.max(abs(expected-np.eye(4)))>.2
    completeness=native['success'].conj().T@native['success']+sum(native['failure_grams'])
    assert np.max(abs(completeness-np.eye(4)))<3e-12
    assert custody()==pins
    output={'schema':'oph.kg.root.full_criterion_review.v1','reviewer':'Codex root KG reviewer',
            'producer_independent':True,'candidate_id':'native-local12-mass-coin-rational-spectroscopy-v1',
            'criterion_sha256':'85edfe6f1eb9c189c6c83ec00cfb45c5bb63b863ecd5c4fe6529bdbe3a77cea5',
            'full_criterion_met':True,'manifest_sha256':sha(PRODUCER/'MANIFEST.json'),
            'independent_projector_rank':26,'independent_gradient_Gram':'delta_ab/48',
            'independent_full_state_cases':cases,'complete_terminal_settings_checked':settings,
            'exact_family_bounds':[str(x) for x in exact_bounds],
            'native_scalar_control_complete_columns':4,'native_scalar_control_error':scalar_error,
            'source_sha256':pins,'all_provider_bytes_unchanged':True,
            'scope':'Declared coherent scalar LOCAL12 source, bandlimited scalar preparation, common macro clock, complete native record convergence; supplied mass and coherent services.',
            'excluded':['all flat bands obeying KG','interacting mass pole','observed particle masses','scalar-to-Dirac source equality','literal enormous pulse execution','official platform award']}
    (HERE/'REVIEW.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'full_criterion_met':True,'full_state_cases':len(cases),'terminal_settings':settings,'native_scalar_error':scalar_error}),flush=True)


if __name__=='__main__':main()
