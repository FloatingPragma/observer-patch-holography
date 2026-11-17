"""Bounded full-rotor charged response of complete variational preparations."""
import json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh,cg,LinearOperator
import source as s
ROOT=Path(__file__).parent
def save(n,x):(ROOT/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def norm(x):return float(np.linalg.norm(x))
potential=json.loads((ROOT/'POTENTIAL_BOUNDS.json').read_text())
results=[];warm={}

for label,m in [('zero',0.),('register_electron',6.0300795404423236e-5)]:
 for K in [2,4]:
  start=time.time();neutral=s.Sector(8,K,m);charged=s.Sector(7,K,m)
  if label=='zero' and K==2 and (ROOT/'PILOT_K2.npz').exists():
   pilot=np.load(ROOT/'PILOT_K2.npz');D=pilot['eigenvalues'];Q=pilot['eigenvectors']
  else:
   if K in warm:v0=warm[K]
   elif K==4 and 2 in warm:
    v0=np.zeros((9,neutral.n),complex);v0[2:7]=warm[2].reshape(5,neutral.n);v0=v0.ravel()
   else:v0=np.cos(np.arange(neutral.dim)+1)+1j*np.sin((np.arange(neutral.dim)+1)/3)
   v0=v0/norm(v0)
   generic=np.cos(np.arange(neutral.dim)+1)+1j*np.sin((np.arange(neutral.dim)+1)/3)
   generic-=v0*np.vdot(v0,generic);generic/=norm(generic)
   v0=v0+generic/16
   D,Q=eigsh(neutral.operator,k=2,ncv=24,which='SA',tol=1e-10,maxiter=3000,v0=v0)
   order=np.argsort(D);D=D[order];Q=Q[:,order]
  psi=Q[:,0]/norm(Q[:,0]);warm[K]=psi.copy()
  hpsi=neutral.matvec(psi);energy=float(np.vdot(psi,hpsi).real)
  trial_upper=energy+1e-6
  exterior_state=neutral.exterior(psi)
  state_res=math.sqrt(norm(hpsi-energy*psi)**2+norm(exterior_state)**2)
  sector_bounds=[]
  for row in potential['all_number_sectors']:
   N=row['N'];shift=0 if N in [0,16]else 4*abs(m)
   lower=row.get('best_Barta_infinite_rotor_lower',row['infinite_rotor_sector_lower_m0'])-shift
   sector_bounds.append({'N':N,'infinite_rotor_lower':lower,'dimension_at_this_cutoff':math.comb(16,N)*(2*K+1)})
  assert sum(x['dimension_at_this_cutoff']for x in sector_bounds)==(1<<16)*(2*K+1)
  vacuum_sector_certified=all(x['infinite_rotor_lower']>trial_upper for x in sector_bounds if x['N']!=8)
  charged_gap=min(sector_bounds[7]['infinite_rotor_lower'],sector_bounds[9]['infinite_rotor_lower'])-trial_upper
  _,_,E=s.one_particle(m)
  vp=np.column_stack([s.between(s.create(8,E[:,a]),psi,neutral.n,neutral.f)for a in range(2)])
  vm=np.column_stack([s.between(s.create(7,E[:,a]).conj().T,psi,neutral.n,neutral.f)for a in range(2)])
  PH=s.particle_hole_fermion(7)
  # Map the full particle sector to the isospectral hole sector, including
  # flux reversal; no particle or hole current is identified by assumption.
  vpminus=np.column_stack([(PH.conj().T@vp[:,a].reshape(neutral.f,charged.n)[::-1].T).T.ravel()for a in range(2)])
  rhs=np.column_stack([vpminus,vm]);shifted=charged.shifted(energy)
  diag=(charged.diagonal-energy).ravel()
  assert min(diag)>0
  pre=LinearOperator((charged.dim,charged.dim),matvec=lambda x:x/diag,dtype=np.complex128)
  solved=[];checks=[]
  for a in range(4):
   calls=[0]
   def count(_):calls[0]+=1
   sol,info=cg(shifted,rhs[:,a],M=pre,rtol=1e-10,atol=1e-12,maxiter=2000,callback=count)
   rin=shifted@sol-rhs[:,a];rout=charged.exterior(sol)
   total=math.sqrt(norm(rin)**2+norm(rout)**2)+1e-8
   bound=total/charged_gap if charged_gap>0 else None
   checks.append({'channel':'particle'if a<2 else'hole','Dirac_probe':a%2,'cg_info':int(info),
    'iterations':calls[0],'finite_residual':norm(rin),'outside_flux_residual':norm(rout),
    'full_infinite_residual_upper':total,'full_infinite_solution_error_upper':bound})
   solved.append(sol)
  solved=np.column_stack(solved)
  GP=-vpminus.conj().T@solved[:,:2];GM=(vm.conj().T@solved[:,2:]).T;G=GP+GM
  CAR=vp.conj().T@vp+(vm.conj().T@vm).T
  Gerror=None;inverse_error=None
  if charged_gap>0:
   Gerror=norm(vpminus)*math.sqrt(sum(c['full_infinite_solution_error_upper']**2 for c in checks[:2]))+norm(vm)*math.sqrt(sum(c['full_infinite_solution_error_upper']**2 for c in checks[2:]))+1e-8
  inverse=-np.linalg.inv(G);invnorm=float(np.linalg.norm(inverse,2))
  if Gerror is not None and invnorm*Gerror<1:
   inverse_error=invnorm**2*Gerror/(1-invnorm*Gerror)
  wave=psi.reshape(neutral.f,neutral.n)
  cos=float(np.vdot(wave[1:],wave[:-1]).real)
  probability=abs(wave)**2
  energy_electric=float(np.sum(probability*neutral.electric))
  # The lower free momentum satisfies |k|<=pi/4 for every holonomy.
  # The interval probability for W_low<=1/4 is a fixed observer event.
  angle=min(np.pi,4*math.acos(3/4))
  differences=neutral.flux[None,:]-neutral.flux[:,None]
  kernel=np.full(differences.shape,angle/np.pi);mask=differences!=0
  kernel[mask]=np.sin(angle*differences[mask])/(np.pi*differences[mask])
  density=wave.conj()@wave.T
  small_wall_probability=float(np.sum(density*kernel).real)
  name=f'{label}_K{K}'
  np.savez_compressed(ROOT/(name+'_PREPARATION.npz'),amplitudes=psi,words=neutral.words,flux=neutral.flux,
   Gauss_prefix_charge=neutral.q,probabilities=probability,finite_low_eigenvalues=D)
  np.savez_compressed(ROOT/(name+'_RESOLVENTS.npz'),rhs=rhs,solutions=solved,words=charged.words,flux=charged.flux,
   Gauss_prefix_charge=charged.q)
  item={'case':name,'S':2,'L':4,'CAR_modes':16,'mass':m,'K':K,
   'every_sector':sector_bounds,'global_vacuum_number_certified':vacuum_sector_certified,
   'vacuum_number':8,'neutral_sector_dimension':neutral.dim,'charged_sector_dimension':charged.dim,
   'variational_energy':energy,'variational_energy_upper':trial_upper,
   'energy_excess_above_true_vacuum_upper':trial_upper-sector_bounds[8]['infinite_rotor_lower'],
   'finite_neutral_low_eigenvalues':D.tolist(),'finite_neutral_ground_order_certified':False,
   'seed_control':'generic component1/16 added to each warm seed; initial pilot has generic seed',
   'finite_preparation_eigen_residual':norm(hpsi-energy*psi),'full_preparation_spectral_residual':state_res,
   'exact_neutral_vacuum_state_certified':False,
   'charged_infinite_resolvent_gap_lower':charged_gap,
   'holonomy_cos_mean_in_preparation':cos,'electric_energy_in_preparation':energy_electric,
   'probability_lower_free_W_le_quarter':small_wall_probability,
   'free_low_W_uniform_upper':1-math.cos(np.pi/4),
   'edge_flux_probability':float(sum(probability[0])+sum(probability[-1])),
   'full_CAR_sum_error':float(np.max(abs(CAR-np.eye(2)))),
   'complete_resolvent_checks':checks,'Green_real':G.real.tolist(),'Green_imag':G.imag.tolist(),
   'infinite_prepared_Green_error_upper':Gerror,
   'inverse_response_Pauli':{name:float(np.trace(P@inverse).real/2)for name,P in [('I',np.eye(2)),('X',s.X),('Y',s.Y),('Z',s.Z)]},
   'infinite_prepared_inverse_error_upper':inverse_error,
   'true_vacuum_propagator_certified':False,'physical_lepton_pole_claim':False,
   'roundoff_admissions':{'Rayleigh_upper_margin':1e-6,'full_residual_margin':1e-8,'Green_product_margin':1e-8},
   'elapsed_seconds':time.time()-start}
  results.append(item);save('RESULT.json',results)
  print(json.dumps({'case':item['case'],'energy':energy,'global_N8':vacuum_sector_certified,
   'charged_gap_lower':charged_gap,'state_residual_full':state_res,'cos':cos,
   'Gerror':Gerror,'inverse_beta':item['inverse_response_Pauli']['Z'],
   'inverse_error':inverse_error,'seconds':item['elapsed_seconds']}),flush=True)
