"""Complete finite quantum geometry coupled to native excitation energy."""
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.linalg import eigh,expm
from model.openwave_energy_lapse.source import native
KAPPA=10000

def basis(cap):return tuple(product(range(cap+1),repeat=2))
@lru_cache(maxsize=4)
def geometry(cap):
 if cap not in(1,2,3,4):raise ValueError('complete product mediator caps1through4 required')
 states=basis(cap);where={v:i for i,v in enumerate(states)};n=len(states);a=[]
 for axis in range(2):
  m=np.zeros((n,n))
  for j,v in enumerate(states):
   if v[axis]:w=list(v);w[axis]-=1;m[where[tuple(w)],j]=np.sqrt(v[axis])
  a.append(m)
 q0=(a[0]+a[0].T)/(2*np.sqrt(KAPPA));q1=(a[1]+a[1].T)/(2*np.sqrt(KAPPA)*3**.25)
 Q=(q0+q1,q0-q1);hg=np.diag([i+np.sqrt(3)*j for i,j in states]);expQ=tuple(expm(q)for q in Q)
 return {'basis':states,'a':a,'Q':Q,'Hg':hg,'expQ':expQ,'vacuum':np.eye(n)[:,0]}
@lru_cache(maxsize=1)
def energies():return np.array([float((g.lo+g.hi)/2)for g in native()[2]])
def block(cap,i,j):
 if i not in range(8)or j not in range(8):raise ValueError('both complete native level labels required')
 g=geometry(cap);e=energies();return g['Hg']+e[i]*g['expQ'][0]+e[j]*g['expQ'][1]
def evolve(cap,i,j,t):
 H=block(cap,i,j);w,v=eigh(H);return v@(np.exp(-1j*t*w)*(v.T@geometry(cap)['vacuum']))
def partial_transpose(rho):return rho.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)
def clock_state(cap,t,level=1):
 branches=np.stack([evolve(cap,i,j,t)for i,j in product((0,level),repeat=2)])
 rho=branches@branches.conj().T/4
 # rho_ab=<psi_b|psi_a>/4 in the source energy basis.
 return branches,rho

def complete_records(cap,t,level=1):
 g=geometry(cap);branches,rho=clock_state(cap,t,level);dim=len(g['basis']);records=[]
 X=np.array([[1,1],[1,-1]],complex)/np.sqrt(2);Y=np.array([[1,-1j],[1,1j]])/np.sqrt(2);Z=np.eye(2)
 for nameA,A in(('X',X),('Y',Y),('Z',Z)):
  for nameB,B in(('X',X),('Y',Y),('Z',Z)):
   vectors=np.kron(A,B)@branches/2;p=abs(vectors)**2
   assert abs(p.sum()-1)<2e-13
   records.append({'axes':nameA+nameB,'active_four_clock_lines_times_all_geometry_occupations':p.tolist(),'other60_native_pair_lines_each_geometry_occupation_probability':'0','complete_outcome_count':64*dim})
 energy0=sum(np.vdot(g['vacuum'],block(cap,i,j)@g['vacuum']).real for i,j in product((0,level),repeat=2))/4
 energy1=sum(np.vdot(v,block(cap,i,j)@v).real for v,(i,j) in zip(branches,product((0,level),repeat=2)))/4
 geometry_energy=sum(np.vdot(v,g['Hg']@v).real for v in branches)/4
 qmeans=[float(sum(np.vdot(v,q@v).real for v in branches)/4)for q in g['Q']]
 # Complete native X/X plus geometry number dephasing need not conserve H.
 A=np.kron(X,X);measured=A@branches/2;recovered=[]
 after=0.
 for a,b in product(range(2),repeat=2):
  for n in range(dim):
   p=abs(measured[2*a+b,n])**2
   for i,j in product(range(2),repeat=2):
    weight=abs(A[2*a+b,2*i+j])**2
    after+=p*weight*block(cap,level*i,level*j)[n,n]
 return {'cap':cap,'source_time':t,'native_clock_level':level,'all9_complete_joint_instruments':records,
  'reduced_clock_density_real':rho.real.tolist(),'reduced_clock_density_imag':rho.imag.tolist(),
  'partial_transpose_eigenvalues':np.linalg.eigvalsh(partial_transpose(rho)).tolist(),
  'source_energy_before':float(energy0),'source_energy_after':float(energy1),'geometry_energy_after':float(geometry_energy),
  'mean_geometry_coordinates':qmeans,'complete_X_X_number_measurement_work':float(after-energy1),
  'native_level_populations_conserved':True,'autonomous_source_evolution':True}

def independent_geometry_control(cap,t):
 """Remove the spatial edge: two independent local native-field factors."""
 n=cap+1;a=np.diag(np.sqrt(np.arange(1,n)),1);frequency=np.sqrt(2)
 q=(a+a.T)/np.sqrt(2*KAPPA*frequency);hg=frequency*np.diag(np.arange(n));vac=np.eye(n)[:,0]
 local=[]
 for level in(0,1):
  e,U=eigh(hg+energies()[level]*expm(q));local.append(U@(np.exp(-1j*t*e)*(U.T@vac)))
 rho=np.array([[np.vdot(local[j],local[i])/2 for j in range(2)]for i in range(2)])
 product_rho=np.kron(rho,rho)
 return {'spatial_Laplacian':[[2,0],[0,2]],'cap':cap,'source_time':t,
  'partial_transpose_eigenvalues':np.linalg.eigvalsh(partial_transpose(product_rho)).tolist(),
  'exact_factorization_for_all_times_and_native_levels':True,
  'product_native_clock_preparation_cannot_entangle':True}
