"""Independent algebraic-source interval residual and complete-current audit."""
from pathlib import Path
import json,math,time
import numpy as np
import mpmath as mp

ROOT=Path(__file__).parent;mp.mp.dps=65;mp.iv.dps=55;iv=mp.iv
started=time.time();r=json.loads((ROOT/'RESULT.json').read_text())
f=json.loads((ROOT/'FULL_TAKAGI.json').read_text())
def interval_matrix(d):
 return [[iv.mpc(iv.mpf(a),iv.mpf(b))for a,b in zip(ar,br)]for ar,br in zip(d['real'],d['imag'])]
def complex_array(d):return np.array(d['real'],float)+1j*np.array(d['imag'],float)
def upper(x):return mp.mpf(x.b)
def lower(x):return mp.mpf(x.a)
def conjugate(x):return iv.mpc(x.real,-x.imag)

# Reconstruct the source independently from literal port vectors. No
# producer module, stored source matrix or numerical eigensolver is used.
phi=(1+iv.sqrt(5))/2;rt2=iv.sqrt(2);p=iv.mpf(1)/12
ports=[[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],
 [0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]]
orbit=[1,5,11,10,7]
P=[]
for port in orbit:
 v=ports[port];norm=sum(x*x for x in v)
 P.append([[v[i]*v[j]/norm for j in range(3)]for i in range(3)])
Y=[[sum(iv.mpf(12)**(4-k)*P[k][i][j]/22621 for k in range(5))for j in range(3)]for i in range(3)]
M=[[iv.mpc(0)for j in range(63)]for i in range(63)]
def put(i,j,v):M[i][j]+=v
for stage in range(5):
 for branch in range(2):
  rr=3+12*stage+3*branch;ll=rr+6
  for i in range(3):put(rr+i,ll+i,1);put(ll+i,rr+i,1)
  for i in range(3):
   for j in range(3):
    V=iv.mpf(int(i==j))-(1+iv.j)*P[stage][i][j]if branch else iv.mpf(int(i==j))
    for parent in range(1 if stage==0 else 2):
     aa=0 if stage==0 else 9+12*(stage-1)+3*parent
     value=p*V/(rt2 if stage==0 else 2)
     put(aa+i,rr+j,value);put(rr+j,aa+i,value)
for a in range(2):
 for b in range(2):
  for i in range(3):
   for j in range(3):put(57+3*a+i,57+3*b+j,Y[i][j]/2)
U=interval_matrix(f['takagi']);masses=[iv.mpf(x)for x in r['full_positive_masses']]
assert len(masses)==63 and len(r['full_signed_doubled_spectrum'])==126
max_res=mp.mpf(0);max_orth=mp.mpf(0)
for i in range(63):
 support=[k for k in range(63)if upper(abs(M[i][k]))>0]
 for j in range(63):
  res=sum((M[i][k]*U[k][j]for k in support),iv.mpc(0))-conjugate(U[i][j])*masses[j]
  max_res=max(max_res,upper(abs(res)))
  gram=sum((conjugate(U[k][i])*U[k][j]for k in range(63)),iv.mpc(0))-int(i==j)
  max_orth=max(max_orth,upper(abs(gram)))
 if i%16==0:print('interval column',i,'seconds',time.time()-started,flush=True)

# The exact source has ||M||<=1+||Y||+2p<=13/6. The polar factor Q
# of the decimal trial U is exactly unitary. The complete Weyl bound
# controls every singular value, including all three tiny light poles.
eta=iv.mpf(str(max_orth))*63
assert upper(eta)<1
polar_error=eta/(1+iv.sqrt(1-eta))
spectral_error=iv.mpf(str(max_res))*63+iv.mpf(13)/3*polar_error
eps=upper(spectral_error)
assert eps<mp.mpf('1e-48')
intervals=[[mp.nstr(lower(m)-eps,62),mp.nstr(upper(m)+eps,62)]for m in masses]
assert all(mp.mpf(intervals[j][1])<mp.mpf(intervals[j+1][0])for j in range(2))
assert mp.mpf(intervals[0][0])>0 and mp.mpf(intervals[3][0])>mp.mpf(13)/30

# Independently recover the charged row frame from the full current and
# full Takagi frame, then check its exact C3 character assignment.
u=complex_array(f['takagi']);B=complex_array(r['full_current']);active=u[:3,:]
rowframe=B@active.conj().T
phi0=(1+math.sqrt(5))/2
vv=np.array([[-1,phi0,0],[1,phi0,0],[-1,-phi0,0],[1,-phi0,0],[0,-1,phi0],[0,1,phi0],
 [0,-1,-phi0],[0,1,-phi0],[phi0,0,-1],[phi0,0,1],[-phi0,0,-1],[-phi0,0,1]])
rot=vv[[5,11,4]].T@np.linalg.inv(vv[[4,5,11]].T)
characters=rowframe@rot@rowframe.conj().T
assert np.max(abs(characters-np.diag(np.diag(characters))))<1e-12
assert np.max(abs(np.diag(characters)-np.array([np.exp(2j*np.pi/3),np.exp(4j*np.pi/3),1])))<1e-12
assert np.max(abs(rowframe@rowframe.conj().T-np.eye(3)))<1e-12
assert np.max(abs(B-rowframe@active))<1e-12
completeness=np.max(abs(B@B.conj().T-np.eye(3)))
assert completeness<1e-12
assert np.max(abs(np.sum(abs(B[:,3:])**2,axis=1)-np.array(r['all_heavy_current_weights'],float)))<1e-12

assert len(r['complete_records'])==18
max_time=max(mp.mpf(q['source_time'])for q in r['complete_records'])
for q in r['complete_records']:
 probs=[mp.mpf(x)for x in q['complete_probabilities']]
 assert len(probs)==5 and min(probs)>-mp.mpf('1e-45')
 assert abs(sum(probs)-1)<mp.mpf('1e-45')
phase_error=max_time*eps
assert phase_error<mp.mpf('1e-18')
# Individual eigenvector labels need not be stable inside a degenerate
# heavy band. Compare complete operators instead. For p>=1, the positive
# energy square-root Sylvester bound gives dE <=(13/3+eps)*eps/2.
# W=(I+p/E)/2 changes by at most dE/(2p), and Duhamel bounds the
# complete evolution by t*dE. The factor50 conservatively includes
# the decimal-to-polar current frame, all five probability outcomes,
# and both preparation and return terms without resolving degeneracies.
energy_operator_error=(mp.mpf(13)/3+eps)*eps/2
record_error=50*(upper(polar_error)+(max_time+1)*energy_operator_error)
assert record_error<mp.mpf('1e-17')

result={'producer_code_imported':False,'stored_source_matrix_used':False,
 'full_Weyl_modes':63,'full_signed_bands':126,'interval_digits':55,
 'maximum_Takagi_equation_entry_residual_upper':mp.nstr(max_res,40),
 'maximum_Gram_entry_residual_upper':mp.nstr(max_orth,40),
 'decimal_frame_to_exact_unitary_norm_bound':mp.nstr(upper(polar_error),40),
 'every_actual_pole_absolute_error_bound':mp.nstr(eps,40),
 'all63_actual_pole_intervals':intervals,
 'charged_character_assignment_reconstructed':True,'all_current_completeness_error':float(completeness),
 'complete_five_outcome_records':18,'maximum_phase_error_from_certified_poles':mp.nstr(phase_error,40),
 'complete_probability_operator_error_upper':mp.nstr(record_error,40),
 'physical_interacting_poles_certified':False,'elapsed_seconds':time.time()-started}
(ROOT/'REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items()if k!='all63_actual_pole_intervals'},indent=2))
