"""Verification-only source certificate and bounded reader; no producer imports."""
from pathlib import Path
import gzip,hashlib,json,math,time
import numpy as np
from reader import bounded_frequency

HERE=Path(__file__).resolve().parent
COMMON=HERE.parent/'2026-11-15_codex_common_selected_source'
ALT=HERE.parents[2]/'row_02_electron_rest_energy/work/2026-11-15_codex_native_work_vacuum/refinement/ALTERNATIVE_OUTPUT.json'
ETA=.01
def unpack(x):return np.array(x['real'])+1j*np.array(x['imag'])
def outcomes(z):return {'x_plus':float((1+z.real)/2),'x_minus':float((1-z.real)/2),'y_plus':float((1+z.imag)/2),'y_minus':float((1-z.imag)/2)}
def inside(interval,x):return interval[0]-1e-11<=x<=interval[1]+1e-11
def certificate(arc,projector):
    eta=float(1-np.trace(arc@projector).real)
    assert -1e-10<=eta<=ETA,'numerical source certificate fails'
    return max(0.,eta)
def verify_source_moment(record,arc,word):
    z=np.trace(arc@word)
    assert max(abs(record[k]-outcomes(z)[k]) for k in record)<1e-10,'source record mismatch'
    return z

def run():
    start=time.monotonic()
    manifest=COMMON/'MANIFEST.sha256'
    assert hashlib.sha256(manifest.read_bytes()).hexdigest()=='89bd1d876a6c428e20f1f0477f3b539c84b1e27292bb32f05e209f97f85887fe'
    for line in manifest.read_text().splitlines():
        digest,name=line.split(None,1);assert hashlib.sha256((COMMON/name.lstrip('*')).read_bytes()).hexdigest()==digest
    v=json.loads(gzip.decompress((COMMON/'v2_run/EVIDENCE.json.gz').read_bytes()));r=v['base_execution']
    edges=((0,1),(0,1),(0,2),(0,2),(1,2));R=np.eye(33);columns=[]
    for p in range(3):
        for e,(a,b) in enumerate(edges):
            for c in (0,1):
                q=[int(i==p) for i in range(3)];d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2
                q[a]-=j;q[b]+=j;y=q.index(1);cc=c^bool(j)
                R[:,11*p+2*e+c]=0;R[11*y+2*e+cc,11*p+2*e+c]=1
                columns.append([p,e,c,j,y,int(cc)])
    eps=.25;a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps**2)];A=np.kron(np.eye(3),a[:,None])
    P=A.T@R@A;W=R@(2*A@A.T-np.eye(33));word=np.linalg.matrix_power(W,4)
    assert np.max(abs(W-unpack(r['clock']['W']) if isinstance(r['clock']['W'],dict) else W-np.array(r['clock']['W'])))<1e-12
    f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6)
    def lifts(x):
        lam=float(x@P@x);J=(R@A@x-lam*A@x)/math.sqrt(1-lam**2)
        return (A@x-1j*J)/math.sqrt(2),(A@x+1j*J)/math.sqrt(2)
    plus,minus=lifts(f);projector=np.outer(plus,plus.conj())+np.outer(minus,minus.conj())
    expected=4*math.acos(float(f@P@f));rows=[]
    for case in r['cases']:
        arc=unpack(case['clock_density']);raw=case['ramsey_records'][1]
        record={'x_plus':raw['x_plus'],'x_minus':1-raw['x_plus'],'y_plus':raw['y_plus'],'y_minus':1-raw['y_plus']}
        verify_source_moment(record,arc,word);eta=certificate(arc,projector);read=bounded_frequency(record,ETA)
        assert inside(read['frequency_interval'],expected)
        assert read['frequency_interval'][1]-read['frequency_interval'][0]<.06
        positive=float(np.vdot(plus,arc@plus).real);negative=float(np.vdot(minus,arc@minus).real)
        noise=[]
        for sign in (-1,1):
            perturbed=dict(record);perturbed['x_plus']+=sign*1e-6;perturbed['x_minus']-=sign*1e-6
            bound=bounded_frequency(perturbed,ETA,1e-6);assert inside(bound['frequency_interval'],expected);noise.append(bound)
        rows.append({'case':{k:case[k] for k in ('M','turns','direction','law')},'record':record,'reader':read,
            'separate_numerical_source_audit':{'actual_off_band_population':eta,'margin_to_fixed_certificate':ETA-eta,'positive_target_weight':positive,'negative_target_weight':negative},
            'deterministic_perturbation_checks':noise,'original_four_pole_diagnostic':case['reader']})
    controls=[]
    for label,vector in [('same-band-negative',minus),('different-positive-band',lifts(g)[0])]:
        record=outcomes(np.vdot(vector,word@vector));read=bounded_frequency(record,ETA)
        omega=math.acos(float((f if label.startswith('same') else g)@P@(f if label.startswith('same') else g)))/eps
        assert inside(read['frequency_interval'],omega)
        if label.startswith('different'):assert not inside(read['frequency_interval'],expected)
        else:assert read['imaginary_moment_separate_sign_record']<0
        controls.append({'label':label,'supplied_spectral_state_sensitivity_only':True,'record':record,'reader':read,'source_frequency_audit':omega})
    alternative=json.loads(ALT.read_text());alt_eta=1-alternative['compact_fidelity'];assert alt_eta>ETA
    controls.append({'label':'retained-always-swap-source','source_file_sha256':hashlib.sha256(ALT.read_bytes()).hexdigest(),
       'off_band_population':alt_eta,'certificate_accepted':False,'reason':'The retained alternative source does not satisfy this fixed narrow-band certificate; no licensed narrow interval is reported.',
       'original_rows':[{'k':x['k'],'frequency':x['inferred_lowest_positive_frequency'],'record':x['records'][1]} for x in alternative['rows']]})
    detectors=[]
    def reject(label,action):
        try:action()
        except (AssertionError,ValueError):detectors.append({'label':label,'rejected':True})
        else:raise AssertionError('detector accepted '+label)
    base=rows[0]['record'];bad=dict(base);bad['x_minus']+=.1
    reject('binary normalization',lambda:bounded_frequency(bad,ETA))
    reject('invalid certificate',lambda:bounded_frequency(base,-.1))
    bad=dict(base);bad['x_plus']+=.1;bad['x_minus']-=.1
    reject('normalized but false source moment',lambda:verify_source_moment(bad,unpack(r['cases'][0]['clock_density']),word))
    reject('false source band certificate',lambda:certificate(np.outer(lifts(g)[0],lifts(g)[0].conj()),projector))
    return {'candidate':'A-native-common-bounded-cosine-reader-v1','scope':'Verification-only bounded readout of sealed finite outcome laws; no new source experiment or finite-shot precision.',
      'fixed_eta':ETA,'native_source_columns':columns,'rows':rows,'controls':controls,'detectors':detectors,
      'maximum_original_fit_condition':max(x['original_four_pole_diagnostic']['fit_condition'] for x in rows),
      'maximum_interval_width':max(x['reader']['frequency_interval'][1]-x['reader']['frequency_interval'][0] for x in rows),
      'passed':True,'seconds':time.monotonic()-start,'external_G6':False,'physical_mass_award':False}

if __name__=='__main__':print(json.dumps(run(),indent=2))
