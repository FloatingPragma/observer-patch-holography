"""Independent rational check of the finite-wall residual, without source imports."""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

I=np.eye(2,dtype=int).astype(object)
X=np.array([[0,1],[1,0]],dtype=object)
Z=np.array([[1,0],[0,-1]],dtype=object)
J=np.array([[0,1],[-1,0]],dtype=object) # i sigma_y
beta=np.kron(Z,I);gamma5=np.kron(X,I)
Pp=(np.eye(4,dtype=object)+gamma5)*F(1,2)
Pm=(np.eye(4,dtype=object)-gamma5)*F(1,2)
T=-(np.kron(Z,I)+np.kron(J,I))*F(1,2)
rows=[]
for S in (2,4,8,16):
 for delta in (F(1,4),F(1,2)):
  H=np.kron(np.eye(S,dtype=object),delta*beta)
  for j in range(S-1):
   H[4*j:4*j+4,4*j+4:4*j+8]=T
   H[4*j+4:4*j+8,4*j:4*j+4]=T.T
  E=np.vstack([delta**j*Pp+delta**(S-1-j)*Pm for j in range(S)])
  C=sum(delta**(2*j)for j in range(S))
  HE=H@E
  assert np.all(E.T@E==C*np.eye(4,dtype=object))
  assert np.all(HE.T@HE==delta**(2*S)*np.eye(4,dtype=object))
  assert np.all(E.T@HE==delta**S*beta)
  rows.append({'S':S,'delta':str(delta),'normalization_squared_inverse':str(C),
   'rho_squared':str(delta**(2*S)/C),
   'effective_boundary_mass':str(delta**S/C),
   'all_three_residual_identities_exact':True})

root=Path(__file__).parent
r=json.loads((root/'RESULT.json').read_text())
bands=json.loads((root/'ALL_BANDS.json').read_text())
records=json.loads((root/'ROTOR_RECORDS.json').read_text())
assert len(bands)==4*64
assert all(len(b['eigenvalues'])==4*b['S']==len(b['boundary_weights'])for b in bands)
assert all(len(x['records'])==16*(2*x['K']+1)for x in records)
assert all(abs(sum(p['probability']for p in x['records'])-1)<1e-12 for x in records)
assert all(all(p['occupation'].bit_count()==2 or p['probability']<1e-25 for p in x['records'])for x in records)
assert r['euclidean_gw']['link_current_divergence_ward_error']<1e-12
assert not r['acceptance']['physical_muon_ratio_repaired']
result={'exact_residuals':rows,'all_band_and_record_counts_pass':True,
 'all_record_probabilities_normalized':True,'fixed_particle_number_pass':True,
 'source_module_imported':False,'physical_acceptance_overclaim':False}
(root/'EXACT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'exact_rational_checks':len(rows),'complete_bands':len(bands),
 'complete_rotor_records':sum(len(x['records'])for x in records)},indent=2))
