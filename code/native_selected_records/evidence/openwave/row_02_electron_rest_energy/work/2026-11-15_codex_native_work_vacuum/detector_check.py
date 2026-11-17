"""Verification-only scientific detector controls, separate from integrity mutations."""
from pathlib import Path
import json,math
import numpy as np
p=Path(__file__).resolve().parent/'refinement'
main=json.loads((p/'OUTPUT.json').read_text());alt=json.loads((p/'ALTERNATIVE_OUTPUT.json').read_text())

def recover(z):
    # Overdetermined recurrence differs from the producer's first-four-equation fit.
    matrix=np.array([[z[n+j] for j in range(4)] for n in range(9)])
    c=np.linalg.lstsq(matrix,-z[4:13],rcond=None)[0]
    roots=np.roots(np.r_[1.,c[::-1]])
    return sorted(float(np.angle(v)) for v in roots)

results=[]
for r,a in zip(main['refinements'],alt['rows']):
    k=r['k'];eps=2.**(-k)
    z=np.array([(2*t['x_plus']-1)+1j*(2*t['y_plus']-1) for t in r['ramsey_records']])
    actual=recover(z);planted=recover(z*np.exp(.1j*np.arange(len(z))))
    assert max(abs(np.array(planted)-np.array(actual)-.1))<1e-6
    competing=np.array([(2*t['x_plus']-1)+1j*(2*t['y_plus']-1) for t in a['records']])
    poles=recover(competing);expected=sorted([sgn*2*math.asin(eps*math.sqrt(t))/eps for t in (.4,.6) for sgn in (-1,1)])
    assert max(abs(np.array(poles)-expected))<1e-7
    for n,v in enumerate(competing):
        predicted=(4/7)*math.cos(n*2*math.asin(eps*math.sqrt(.4))/eps)+(3/7)*math.cos(n*2*math.asin(eps*math.sqrt(.6))/eps)
        assert abs(v-predicted)<1e-9
    Q=len(a['all_counter_probabilities']);grid=2*math.pi*np.arange(Q)/Q;prob=np.zeros(Q)
    for weight,t in ((4/7,.4),(3/7,.6)):
        theta=2*math.asin(eps*math.sqrt(t))
        for sign in (-1,1):
            angle=(sign*theta-grid)/2;den=np.sin(angle);ratio=np.ones(Q);good=abs(den)>1e-12
            ratio[good]=np.sin(Q*angle[good])/(Q*den[good]);prob+=weight*ratio**2/2
    assert max(abs(prob-np.array(a['all_counter_probabilities'])))<1e-10
    results.append({'k':k,'main_signed_frequencies':actual,'planted_shift_detected':.1,'competing_signed_frequencies':poles,'all_changed_law_counter_bins_checked':Q})
print(json.dumps({'verdict':'PASS','scope':'Relevant always-swap source control and separately labelled planted phase sensitivity. No emergence credit for the plant.','results':results},indent=2))
