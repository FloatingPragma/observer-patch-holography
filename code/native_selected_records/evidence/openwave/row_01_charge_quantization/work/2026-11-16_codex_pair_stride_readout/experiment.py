"""Registered new pair-record analysis; source evidence retained unchanged."""
from pathlib import Path
from itertools import product
import gzip,hashlib,json,math,time
import numpy as np
from reader import read

HERE=Path(__file__).resolve().parent
FAMILY=HERE.parent/'2026-11-16_codex_common_refinement_family'
ETA=.01;DELTA=1e-5;STRIDE=8;SCALE=511*math.sqrt(5)/2
OLD_PAIR=1024.4508487023415
def unpack(x):return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)
def probabilities(z):return {'x_plus':float((1+z.real)/2),'x_minus':float((1-z.real)/2),'y_plus':float((1+z.imag)/2),'y_minus':float((1-z.imag)/2)}
def unique(result):
    assert result['unique_component'],'multiple or empty compatible branches'
    return result['frequency_components'][0]
def contains(interval,x):return interval[0]-1e-10<=x<=interval[1]+1e-10
def altered(record,sx,sy):
    q=dict(record)
    for axis,s in (('x',sx),('y',sy)):
        q[axis+'_plus']+=s*DELTA;q[axis+'_minus']-=s*DELTA
    return q

def run():
    begin=time.monotonic();bindings=json.loads((HERE/'INPUTS.json').read_text())
    for name,digest in bindings['source_files'].items():
        assert hashlib.sha256((FAMILY/name).read_bytes()).hexdigest()==digest,name
    rows=[];controls=[];all_columns=[]
    gates={'single_source_columns':True,'pair_off_band_certificate':True,'all_aliases_enumerated_unique':True,
      'nominal_and_error_intervals_contain_source_frequency':True,'all_corner_perturbations_contain_source_frequency':True,
      'k6_excludes_old_pair_energy':True,'k6_corner_envelope_excludes_old_pair_energy':True,
      'k6_shifted_clock_distinguished':True,'higher_state_distinguished':True}
    for k in range(7):
        evidence=json.loads(gzip.decompress((FAMILY/f'k{k}/EVIDENCE.json.gz').read_bytes()));base=evidence['base_execution']
        eps=2.**-k;R=np.eye(33);columns=[]
        for p in range(3):
            for e,(a,b) in enumerate(((0,1),(0,1),(0,2),(0,2),(1,2))):
                for c in (0,1):
                    q=[int(x==p) for x in range(3)];d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2
                    q[a]-=j;q[b]+=j;y=q.index(1);cc=c^bool(j)
                    R[:,11*p+2*e+c]=0;R[11*y+2*e+cc,11*p+2*e+c]=1
                    columns.append([p,e,c,j,y,int(cc)])
        a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps**2)]
        A=np.kron(np.eye(3),a[:,None]);P=A.T@R@A;W=R@(2*A@A.T-np.eye(33))
        gates['single_source_columns'] &= np.max(abs(W-unpack(base['clock']['W'])))<1e-12
        # Native word powers, not a target Hamiltonian exponential.
        word=np.eye(33)
        for _ in range(2**k):word=W@word
        later=np.eye(33)
        for _ in range(STRIDE):later=word@later
        joint_words={1:np.kron(word,word),STRIDE:np.kron(later,later)}
        f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6)
        def lifts(x):
            lam=float(x@P@x);aux=(R@A@x-lam*A@x)/math.sqrt(1-lam*lam)
            return (A@x-1j*aux)/math.sqrt(2),(A@x+1j*aux)/math.sqrt(2)
        plus,minus=lifts(f);source_frequency=2*math.acos(float(f@P@f))/eps
        all_columns.append({'k':k,'columns':columns,'R_events_per_first_record':2*2**k,'R_events_per_stride8_record':16*2**k})
        for case in base['cases']:
            rho=unpack(case['clock_density']);joint=np.kron(rho,rho)
            z={m:complex(np.einsum('ab,ba->',joint,U)) for m,U in joint_words.items()}
            records={m:probabilities(value) for m,value in z.items()}
            pp=float(np.vdot(plus,rho@plus).real);pm=float(np.vdot(minus,rho@minus).real)
            eta=1-pp**2-pm**2;gates['pair_off_band_certificate'] &= -1e-10<=eta<=ETA
            nominal=read(records[1],records[STRIDE],probability_error=0)
            bounded=read(records[1],records[STRIDE],probability_error=DELTA)
            envelope=read(records[1],records[STRIDE],probability_error=2*DELTA)
            gates['all_aliases_enumerated_unique'] &= all(q['unique_component'] for q in (nominal,bounded,envelope))
            intervals=[unique(q) for q in (nominal,bounded,envelope)]
            gates['nominal_and_error_intervals_contain_source_frequency'] &= all(contains(q,source_frequency) for q in intervals)
            corners=[]
            for signs in product((-1,1),repeat=4):
                first=altered(records[1],*signs[:2]);last=altered(records[STRIDE],*signs[2:]);result=read(first,last)
                interval=unique(result);gates['all_corner_perturbations_contain_source_frequency'] &= contains(interval,source_frequency)
                corners.append({'probability_error_signs':signs,'first_record':first,'stride8_record':last,'reader':result})
            if k==6:
                gates['k6_excludes_old_pair_energy'] &= SCALE*intervals[1][1]<OLD_PAIR
                gates['k6_corner_envelope_excludes_old_pair_energy'] &= SCALE*intervals[2][1]<OLD_PAIR
                shift=(OLD_PAIR-1022)/SCALE
                shifted={m:probabilities(z[m]*np.exp(1j*m*shift)) for m in (1,STRIDE)}
                changed=read(shifted[1],shifted[STRIDE]);changed_interval=unique(changed)
                # Only the +/+ ray is the desired eigenphase of this supplied shift.
                shifted_eta=1-pp**2;assert shifted_eta<=ETA
                gates['k6_shifted_clock_distinguished'] &= contains(changed_interval,source_frequency+shift) and not contains(changed_interval,source_frequency)
                controls.append({'case':{x:case[x] for x in ('M','turns','direction','law')},'label':'supplied common phase shift of actual clock',
                   'integrity_sensitivity_only':True,'shift_radians':shift,'records':shifted,'reader':changed,'separate_off_band_certificate':shifted_eta})
            rows.append({'k':k,'case':{x:case[x] for x in ('M','turns','direction','law')},'pair_records':records,
              'nominal_reader':nominal,'bounded_reader':bounded,'worst_corner_envelope_reader':envelope,
              'corner_perturbations':corners,'numerical_source_certificate':{'eta':eta,'fixed_upper_bound':ETA,'positive_pair_weight':pp**2,'negative_pair_weight':pm**2},
              'post_readout_comparison':{'source_frequency':source_frequency,'source_prediction_keV':SCALE*source_frequency,
                 'bounded_interval_keV':[SCALE*x for x in intervals[1]],'worst_corner_interval_keV':[SCALE*x for x in intervals[2]],'old_pair_keV':OLD_PAIR if k==6 else None},
              'original_four_pole_diagnostic':case['reader']})
        if k==6:
            gg=lifts(g)[0];joint=np.outer(np.kron(gg,gg),np.kron(gg,gg).conj())
            rec={m:probabilities(complex(np.einsum('ab,ba->',joint,U))) for m,U in joint_words.items()}
            result=read(rec[1],rec[STRIDE]);interval=unique(result)
            gates['higher_state_distinguished'] &= not contains(interval,source_frequency)
            controls.append({'label':'supplied higher-mode pair','integrity_sensitivity_only':True,'records':rec,'reader':result})
    return {'candidate':'A-native-common-pair-stride-readout-v1','fixed_stride':STRIDE,'fixed_eta':ETA,'fixed_probability_error':DELTA,
      'one_limiting_calibration_keV_per_frequency':SCALE,'rows':rows,'source_columns':all_columns,'controls':controls,
      'gates':{k:bool(v) for k,v in gates.items()},'passed':all(gates.values()),'seconds':time.monotonic()-begin,
      'scope':'New registered analysis of known source evidence; exact outcome laws, numerical certificates and deterministic perturbations, not a sampled rare-herald experiment.',
      'finite_shot_confidence_claim':False,'external_G6':False,'physical_mass_award':False}

if __name__=='__main__':
    result=run();print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
