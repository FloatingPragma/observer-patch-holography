"""Every number sector, complete particle/hole spectra and flux-tail bounds."""
import json,math,time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
import source as s

ROOT=Path(__file__).parent
def save(name,x):(ROOT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def norm(a):return float(np.linalg.norm(a))
def root_low(a,d,c):return (a+d-math.sqrt((a-d)**2+4*c*c))/2

def diagonalize(H):
 D,Q=eigh(H.toarray(),driver='evd')
 n=len(D);u=np.finfo(float).eps/2
 R=H@Q-Q*D
 gram=Q.conj().T@Q-np.eye(n)
 # Conservative standard-IEEE product envelopes. Each complex dot is
 # bounded by16n real operations; no interval library is substituted.
 gamma=16*n*u/(1-16*n*u)
 qf=norm(Q);hf=norm(H.data)
 eta=norm(gram)+gamma*qf*qf+1e-12
 residual=norm(R)+gamma*(hf*qf+norm(Q*D))+1e-10
 similarity_error=residual/math.sqrt(1-eta)
 # Orthogonal-polar comparison is used for the boundary projector bound.
 polar_error=residual/math.sqrt(1-eta)+2*math.sqrt(1+eta)*max(abs(D))*(1/math.sqrt(1-eta)-1)
 assert eta<1e-3
 return D,Q,{'orthogonality_bound':eta,'residual_bound':residual,
  'finite_eigenvalue_error_bound':similarity_error,'polar_operator_error_bound':polar_error,
  'roundoff_model':'IEEE double; complex dot <=16n real operations; assembly margin1e-10'}

def tail_bound(D,Q,cert,S,N,K,m):
 amin,c,d=s.exterior_bound(S,N,K,m)
 # Outward numerical envelope for the elementary analytic exterior bound.
 d-=1e-9;c+=1e-9
 e=cert['finite_eigenvalue_error_bound'];pe=cert['polar_operator_error_bound']
 coarse=root_low(float(D[0]-e),d,c)
 # Declare an eigenvalue-certification tolerance; it does not select a
 # physical parameter. Test positivity at x=finite_estimate−1e−4.
 x=float(D[0]-1e-4)
 result={'A_lower_formula_value':amin,'coupling_norm_upper':c,'discarded_block_lower':d,
  'coarse_infinite_ground_lower':coarse,'finite_ground_upper':float(D[0]+e),
  'candidate_lower':x,'candidate_verified':False}
 if d<=x:return result
 delta=c*c/(d-x)
 # Hfull−x has Schur complement >=HK−x−delta P_edge.
 # The complete finite eigensystem supplies a low/complement estimate
 # for this boundary perturbation, including every omitted finite level.
 J=int(np.searchsorted(D,D[0]+2*delta+1));J=max(1,min(J,len(D)-1))
 block=math.comb(2*S,N)
 edge=Q[np.r_[0:block,len(D)-block:len(D)],:J]
 eta=cert['orthogonality_bound']
 polar_move=math.sqrt(1+eta)*(1/math.sqrt(1-eta)-1)
 p=min(1,(norm(edge)+math.sqrt(J)*polar_move)**2)
 low=float(D[0])-delta*p
 complement=float(D[J])-delta
 bound=root_low(low,complement,delta*math.sqrt(p))-pe
 result.update({'edge_schur_delta':delta,'low_space_dimension':J,
  'edge_low_space_norm_squared_upper':p,'perturbed_finite_lower':bound,
  'candidate_verified':bool(bound>=x),
  'certified_infinite_ground_lower':x if bound>=x else coarse})
 return result

results=[];spectra=[]
for S in [2,4]:
 for masslabel,m in [('zero',0.),('register_electron',6.0300795404423236e-5)]:
  for K in [4,8,16]:
   start=time.time();blocks={};DQs={};sectors=[]
   for N in range(S+1):
    H,A,C,E=s.block(S,N,K,m);D,Q,cert=diagonalize(H)
    bound=tail_bound(D,Q,cert,S,N,K,m)
    blocks[N]=H;DQs[N]=(D,Q,cert)
    sectors.append({'N':N,'dimension':len(D),'ground':float(D[0]),
      'first_gap':float(D[1]-D[0]) if len(D)>1 else None,
      'ground_residual':norm(H@Q[:,0]-D[0]*Q[:,0]),
      'finite_certificate':cert,'infinite_tail':bound})
   # The exact CAR particle-hole transform U=reflection_s tensor sigma_y
   # maps N to2S−N, leaves the loop flux unchanged, and preserves H.
   for N in range(S+1,2*S+1):
    mirror=2*S-N
    row=dict(sectors[mirror]);row['N']=N;row['particle_hole_of']=mirror
    sectors.append(row)
   vac=min(sectors,key=lambda r:r['ground']);N0=vac['N']
   assert N0==S
   D,Q,cert=DQs[S];state=Q[:,0];E0=float(D[0]);H=blocks[S]
   dimN=math.comb(2*S,S);wf=state.reshape(2*K+1,dimN)
   cos=float(np.vdot(wf[1:],wf[:-1]).real)
   sin=float(np.vdot(wf[1:],wf[:-1]).imag)
   qone=np.diag([1]*S+[-1]*S);Qwall=sparse.kron(sparse.eye(2*K+1),sparse.csr_matrix(s.lift(qone,S)))
   axial_vector=Qwall@state
   midpoint_norm=norm(H@axial_vector-E0*axial_vector)
   # Particle and hole sectors are diagonalized independently here to
   # retain actual boundary-current residues in both sectors.
   Dminus,Qminus,cminus=DQs[S-1]
   Hplus,Ap,Cp,E=s.block(S,S+1,K,m)
   Dplus,Qplus,cplus=diagonalize(Hplus)
   creations=[s.creation(S,S,E[:,a])for a in range(2)]
   annihilations=[s.creation(S,S-1,E[:,a]).conj().T for a in range(2)]
   vp=np.column_stack([sparse.kron(sparse.eye(2*K+1),a)@state for a in creations])
   vm=np.column_stack([sparse.kron(sparse.eye(2*K+1),a)@state for a in annihilations])
   ampP=Qplus.conj().T@vp;ampM=Qminus.conj().T@vm
   gapP=Dplus-E0;gapM=Dminus-E0
   assert min(gapP)>0 and min(gapM)>0
   # Retarded zero-frequency Green matrix: addition negative, removal positive.
   GP=-(ampP.conj().T/gapP)@ampP
   GM=((ampM.conj().T/gapM)@ampM).T
   G0=GP+GM;inverse=np.linalg.inv(G0)
   boundary_mass_matrix=-inverse
   coefficients={name:float(np.trace(P@boundary_mass_matrix).real/2)for name,P in [('I',np.eye(2)),('X',s.X),('Y',s.Y),('Z',s.Z)]}
   residuesP=np.einsum('ia,ib->iab',ampP.conj(),ampP)
   residuesM=np.einsum('ib,ia->iab',ampM.conj(),ampM)
   total=sum(residuesP)+sum(residuesM)
   def pack(ds,res):
    return [{'gap':float(e),'residue_real':R.real.tolist(),'residue_imag':R.imag.tolist()}for e,R in zip(ds,res)]
   spectra.append({'S':S,'mass':masslabel,'K':K,'particle':pack(gapP,residuesP),'hole':pack(gapM,residuesM)})
   # Both charged transitions and vacuum use certified infinite-sector
   # intervals when the edge Schur test succeeds.
   lows={r['N']:r['infinite_tail'].get('certified_infinite_ground_lower',r['infinite_tail']['coarse_infinite_ground_lower'])for r in sectors}
   highs={r['N']:r['infinite_tail']['finite_ground_upper']for r in sectors}
   vac_identified=all(lows[n]>highs[S]for n in range(2*S+1)if n!=S)
   item={'S':S,'mass_label':masslabel,'mass':m,'K':K,'all_number_sectors':sectors,
    'vacuum_number':S,'global_number_sector_certified':vac_identified,
    'vacuum_energy':E0,'holonomy_cos_mean':cos,'holonomy_sin_mean':sin,
    'vacuum_edge_flux_probability':float(norm(wf[0])**2+norm(wf[-1])**2),
    'axial_mean':float(np.vdot(state,axial_vector).real),
    'midpoint_commutator_on_vacuum_norm':midpoint_norm,
    'particle_lowest_gap':float(gapP[0]),'hole_lowest_gap':float(gapM[0]),
    'charged_transition_infinite_interval':[lows[S-1]-highs[S],highs[S-1]-lows[S]],
    'boundary_CAR_sum_error':float(np.max(abs(total-np.eye(2)))),
    'complete_particle_levels':len(gapP),'complete_hole_levels':len(gapM),
    'zero_frequency_Green_real':G0.real.tolist(),'zero_frequency_Green_imag':G0.imag.tolist(),
    'minus_G_inverse_Pauli_coefficients':coefficients,
    'electron_pole_identification':False,
    'elapsed_seconds':time.time()-start}
   results.append(item);save('RESULT.json',results);save('ALL_BOUNDARY_SPECTRA.json',spectra)
   print(json.dumps({'S':S,'mass':masslabel,'K':K,'E0':E0,'gap':float(gapP[0]),'cos':cos,
    'mass_matrix_Z':coefficients['Z'],'CAR_error':item['boundary_CAR_sum_error'],
    'number_certified':vac_identified,'gap_interval':item['charged_transition_infinite_interval'],
    'seconds':item['elapsed_seconds']}),flush=True)
