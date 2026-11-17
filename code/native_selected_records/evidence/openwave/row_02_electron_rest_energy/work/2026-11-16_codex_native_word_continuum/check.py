"""Verification of an existing native-family theorem; no files/network/processes."""
import json,math
import numpy as np
EDGES=[(0,1),(0,1),(0,2),(0,2),(1,2)]
R=np.eye(33);records=[]
for p in range(3):
 for e,(a,b) in enumerate(EDGES):
  for c in (0,1):
   q=[int(x==p) for x in range(3)];d=q[a]-q[b];j=(d-(2*c-1)*(abs(d)%2))//2
   q[a]-=j;q[b]+=j;y=q.index(1);nc=c^int(j!=0)
   R[:,11*p+2*e+c]=0;R[11*y+2*e+nc,11*p+2*e+c]=1;records.append([p,e,c,j,y,nc])
a=np.r_[np.ones(10)/math.sqrt(10),0.];A=np.kron(np.eye(3),a[:,None]);L=10*(np.eye(3)-A.T@R@A)
values,V=np.linalg.eigh(L);values[np.abs(values)<1e-12]=0
assert np.max(abs(values-np.array([0.,4.,6.])))<1e-12
energies=np.sqrt(values/5);rows=[]
for k in range(7):
 eps=2.**(-k);a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps*eps)]
 A=np.kron(np.eye(3),a[:,None]);W=R@np.kron(np.eye(3),2*np.outer(a,a)-np.eye(11))
 chart=[]
 for lam,v in zip(values,V.T):
  u=A@v
  if lam==0:chart.append(u)
  else:
   transfer=1-eps*eps*lam/10;w=(R@u-transfer*u)/math.sqrt(1-transfer*transfer)
   chart.append((u-1j*w)/math.sqrt(2))
 J=np.array(chart).T;live=J.copy();worst=0.;count=0;all_times=[]
 assert np.max(abs(J.conj().T@J-np.eye(3)))<1e-11
 for n in range(int(8/eps)+1):
  observed=J.conj().T@live
  for fraction in (0.,.25,.75):
   t=(n+fraction)*eps
   if t>8:continue
   target=np.diag(np.exp(1j*t*energies))
   error=float(np.linalg.norm(observed-target,ord=2));bound=.1*t*eps*eps+1.2*eps
   assert error<=bound+1e-11
   all_times.append([t,n,error,bound]);worst=max(worst,error);count+=1
  live=W@live
 rows.append({'k':k,'epsilon':eps,'checked_times':count,'maximum_operator_error':worst,'uniform_T8_bound':.8*eps*eps+1.2*eps,'records':all_times})
print(json.dumps({'kind':'Verification-only of the exact native-word limit theorem','source_endpoint_columns':records,'native_formula_evaluations':30,'generator_eigenvalues':energies.tolist(),'coherent_word_matrix_powers_not_literal_native_calls':True,'rows':rows,'verdict':'PASS'},indent=2))
