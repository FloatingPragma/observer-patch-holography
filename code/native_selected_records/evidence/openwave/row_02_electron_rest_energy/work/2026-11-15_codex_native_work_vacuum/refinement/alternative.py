"""Same finite reader and geometry under a reversible always-swap source."""
import json,math
import numpy as np
EDGES=[(0,1),(0,1),(0,2),(0,2),(1,2)]
R=np.eye(33);records=[]
for p in range(3):
 for e,(a,b) in enumerate(EDGES):
  for c in (0,1):
   q=[int(x==p) for x in range(3)];out=q.copy();out[a],out[b]=out[b],out[a]
   j=q[a]-out[a];y=out.index(1);i=11*p+2*e+c;o=11*y+2*e+c
   R[:,i]=0;R[o,i]=1;records.append([p,e,c,j,y,c])
a=np.r_[np.ones(10)/math.sqrt(10),0.];A=np.kron(np.eye(3),a[:,None]);P=A.T@R@A
s=(np.eye(3)-P)@np.array([0.,1.,0.])/2
for n in range(8):s=P@s
survival=float(s@s);s/=math.sqrt(survival);f=np.array([0.,1.,-1.])/math.sqrt(2)
rows=[]
for k in range(7):
 eps=2.**(-k);stride=2**k;a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps*eps)]
 A=np.kron(np.eye(3),a[:,None]);C=np.kron(np.eye(3),2*np.outer(a,a)-np.eye(11));W=R@C
 clock=np.eye(33,dtype=complex)
 for j in range(stride):clock=W@clock
 u=A@s;live=u.astype(complex);moments=[];data=[]
 for n in range(17):
  z=np.vdot(u,live);data.append({'time':n,'words':n*stride,'x_plus':float((1+z.real)/2),'y_plus':float((1+z.imag)/2)})
  moments.append((2*data[-1]['x_plus']-1)+1j*(2*data[-1]['y_plus']-1));live=clock@live
 z=np.array(moments);coeff=np.linalg.solve(np.array([[z[n+j] for j in range(4)] for n in range(4)]),-z[4:8])
 roots=np.roots(np.r_[1.,coeff[::-1]]);roots=roots[np.argsort(np.angle(roots))]
 weights=np.linalg.solve(np.array([[p**n for p in roots] for n in range(4)]),z[:4])
 fit=np.array([sum(w*p**n for p,w in zip(roots,weights)) for n in range(17)])
 positive=[float(np.angle(p)) for p,w in zip(roots,weights) if np.angle(p)>0 and w.real>.01]
 frequency=min(positive);baseline=2*math.asin(eps/math.sqrt(5))/eps
 Q=256*stride;live=u.astype(complex);waves=[]
 for n in range(Q):waves.append(live.copy());live=W@live
 amps=np.fft.fft(np.array(waves),axis=0)/Q;prob=np.sum(abs(amps)**2,axis=1)
 rows.append({'k':k,'records':data,'signed_poles':[float(np.angle(p)) for p in roots],
  'weights':[[float(w.real),float(w.imag)] for w in weights],
  'inferred_lowest_positive_frequency':frequency,'baseline_frequency':baseline,
  'source_frequency_difference':frequency-baseline,'heldout_error':float(max(abs(fit[8:]-z[8:]))),
  'all_counter_probabilities':prob.tolist(),'complete_counter_mass':float(sum(prob)),
  'baseline_uniform_frequency_bound_rejected':frequency>1.2})
output={'candidate':'A-native-triangle-refined-work-v1-always-swap-control',
 'scope':'Relevant changed-law experiment; no literal unchanged-native calls are claimed.',
 'source_columns':records,'cooling_success':survival,'compact_fidelity':float((f@s)**2),'rows':rows}
print(json.dumps(output,indent=2,sort_keys=True))
assert abs(float((f@s)**2)-4/7)<1e-12
assert all(r['source_frequency_difference']>.3 and r['heldout_error']<1e-8 and abs(r['complete_counter_mass']-1)<1e-10 for r in rows)
