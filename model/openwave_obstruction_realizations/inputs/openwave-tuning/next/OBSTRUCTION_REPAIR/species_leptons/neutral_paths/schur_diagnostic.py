"""Zero-frequency approximation only; the full source is run.py."""
from pathlib import Path
import json,numpy as np
r=(1+np.sqrt(5))/2
v=np.array([[-1,r,0],[1,r,0],[-1,-r,0],[1,-r,0],[0,-1,r],[0,1,r],
 [0,-1,-r],[0,1,-r],[r,0,-1],[r,0,1],[-r,0,-1],[-r,0,1]])
P=[np.outer(v[i],v[i])/(v[i]@v[i])for i in [1,5,11,10,7]]
Y=sum(12**(4-k)*P[k]/22621 for k in range(5));K=np.eye(3,dtype=complex)
for q in P:K=K@(np.eye(3)+(-1-1j)*q/2)
S=K@Y@K.T;m2,U=np.linalg.eigh(S.conj().T@S)
R=v[[5,11,4]].T@np.linalg.inv(v[[4,5,11]].T)
Q=np.eye(3)+(np.exp(2j/9)*R+np.exp(-2j/9)*R.T)/np.sqrt(2)
_,Fe=np.linalg.eigh(Q@Q);B=Fe.conj().T@U;w=abs(B)**2
out={'diagnostic_only':'zero-frequency Schur singular vectors, no residue or full-pole substitution',
 'unscaled_Schur_singular_values':np.sqrt(m2).tolist(),
 'squared_splitting_ratio':float((m2[1]-m2[0])/(m2[2]-m2[0])),
 'orientation_sin2_12_13_23':[float(w[0,1]/(1-w[0,2])),float(w[0,2]),float(w[1,2]/(1-w[0,2]))]}
(Path(__file__).parent/'SCHUR_DIAGNOSTIC.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
