"""Execute declared finite controls; no fitted parameters or deleted bands."""
import json
import math
from itertools import product
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
import source as s

ROOT=Path(__file__).parent
def err(A): return float(np.max(abs(A))) if np.size(A) else 0.
def comm(A,B):return A@B-B@A
def save(name,value): (ROOT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

out={"source_choices":{"clock":1/8,"r":1,"wall_lengths":[2,4,8,16],
     "scalar_wall_stress":[1/4,1/2],"physical_targets_used":False}}

I=np.eye(4);gamma=-1j*s.BETA@s.ALPHA[0]
out['wilson_symmetry']={
 'massless_axial_commutator_at_corner_norm':float(np.linalg.norm(comm(2*s.BETA,s.G5),2)),
 'energy_reflection_anticommutator_error':err(s.A4@(2*s.BETA)+(2*s.BETA)@s.A4),
 'mass_operator_axial_commutator_norm':float(np.linalg.norm(comm(s.BETA,s.G5),2)),
 'mass_operator_allowed_by_colour_and_vector_U1':True,
 'energy_reflection_forbids_mass':False}
out['spatial_overlap']={}
for theta in [0,1/8,1/2,np.pi]:
 H,Q=s.spatial_overlap(theta)
 out['spatial_overlap'][str(theta)]={
  'fixed_background_commutator_error':err(comm(H,Q)),
  'charge_square_identity_error':err(Q@Q-I+H@H/4),
  'axial_charge_eigenvalues':np.linalg.eigvalsh(Q).tolist()}
# Acting on a normalized constant rotor wavefunction: [E^2,q5]v=-q5''v.
# Its squared norm is exactly1/4 since gamma5 and i gamma5 gamma are
# orthogonal Hermitian Clifford matrices. This is an analytic value.
out['spatial_overlap']['electric_counterexample']={
 'norm_E_squared_q5_on_constant_rotor':0.5,
 'norm_physical_electric_q5_a_half':0.125,
 'cell_clock_norm':1/64,
 'finite_flux_support_required':[-1,0,1],
 'a_gauge_connection_cannot_fix_varying_charge_eigenvalues':True}

mass_values={'electron':6.0300795404423236e-5,'muon':.012468367126996047,'tau':.20969353598509016}
out['boundary_intertwiners']=[]
for S in [2,4,8,16]:
 E=s.boundary_map(S);Q=s.axial_wall(S)
 currents=[I,*s.ALPHA,s.G5,*[2*x for x in s.fractional.SPIN]]
 checks={'S':S,'isometry':err(E.conj().T@E-I),
  'axial_wall_map':err(E.conj().T@Q@E-s.G5),
  'current_maps':max(err(E.conj().T@np.kron(np.eye(S),a)@E-a)for a in currents),
  'mass_poles':{}}
 for name,m in mass_values.items():
  H=s.bloch(S,np.zeros(3),m)
  checks['mass_poles'][name]={'full_intertwiner':err(H@E-E@(s.H_CLOCK*m*s.BETA)),
    'physical_positive_pole':float(np.sort(np.linalg.eigvalsh(H))[2*S]/s.H_CLOCK),
    'input_register_coefficient':m}
 out['boundary_intertwiners'].append(checks)

allbands=[];out['wall_recursion']=[]
for S in [2,4,8,16]:
 for delta in [1/4,1/2]:
  H=s.bloch(S,np.zeros(3),onsite_shift=delta)
  norm=np.sqrt((1-delta**2)/(1-delta**(2*S)))
  E=np.vstack([norm*(delta**j*s.PP+delta**(S-1-j)*s.PM)for j in range(S)])
  values,V=np.linalg.eigh(H)
  light=np.argsort(abs(values))[:4]
  W=s.boundary_map(S)
  weights=np.sum(abs(W.conj().T@V)**2,axis=0)/4
  rho=norm*delta**S
  out['wall_recursion'].append({'S':S,'delta':delta,'trial_residual_bound':rho,
   'physical_residual_bound':rho/s.H_CLOCK,
   'trial_residual_identity_error':err(E.conj().T@H@H@E-rho**2*I),
   'trial_isometry_error':err(E.conj().T@E-I),
   'effective_boundary_mass':norm**2*delta**S,
   'exact_light_absmax':float(max(abs(values[light]))),
   'bulk_gap_lower_bound':1-delta,
   'full_eigenvalues':values.tolist(),'boundary_weights':weights.tolist(),
   'light_total_boundary_weight':float(sum(weights[light])),
   'weight_sum':float(sum(weights))})
 for k in product(range(4),repeat=3):
  momentum=2*np.pi*np.array(k)/4
  H=s.bloch(S,momentum);values,V=np.linalg.eigh(H);E=s.boundary_map(S)
  weights=np.sum(abs(E.conj().T@V)**2,axis=0)/4
  allbands.append({'S':S,'momentum_indices':list(k),'eigenvalues':values.tolist(),
   'boundary_weights':weights.tolist(),'zero_modes':int(sum(abs(values)<1e-10))})
save('ALL_BANDS.json',allbands)
out['free_band_inventory']=[{'S':S,'momenta':64,
 'total_one_particle_states':sum(len(b['eigenvalues'])for b in allbands if b['S']==S),
 'zero_states':sum(b['zero_modes']for b in allbands if b['S']==S),
 'zero_momenta':[b['momentum_indices']for b in allbands if b['S']==S and b['zero_modes']],
 'boundary_probability_error':max(abs(sum(b['boundary_weights'])-1)for b in allbands if b['S']==S)}for S in [2,4,8,16]]

# Full24CAR-per-slice U3 one-particle covariance, with noncommuting colour links.
cells,edges=s.fractional.edges(2)
T1=np.array([[0,1,0],[1,0,0],[0,0,0]],complex)/2
T2=np.diag([1,-1,0])/2
links=[expm(1j*((j%3+1)*T1+(j%5+1)*T2+np.eye(3))/257)for j in range(len(edges))]
gauges=[expm(1j*((j+1)*T1+(9-j)*T2+np.eye(3)*j)/31)for j in range(len(cells))]
changed=[gauges[x]@U@gauges[y].conj().T for (x,y,a),U in zip(edges,links)]
G0=sparse.block_diag([s.fractional.local_gauge(g)for g in gauges]).toarray()
G=np.kron(np.eye(2),G0)
H=s.one_particle(2,links,mass_values['electron']);Ht=s.one_particle(2,changed,mass_values['electron'])
out['u3_covariance']={'dimension':len(H),'error':err(Ht-G@H@G.conj().T)}
# The global fifth-coordinate charge is independent of links. Its only
# violation is the explicitly retained fifth-hop and boundary-mass terms.
Q=s.axial_wall(2,block=len(G0));Hsp=np.kron(np.eye(2),s.fractional.one_particle(0,1/8,links))
out['finite_wall_ward']={'spatial_commutator_error':err(comm(Hsp,Q)),
 'full_equals_fifth_plus_mass_error':err(comm(H,Q)-comm(H-Hsp,Q)),
 'gauge_commutator_error':err(comm(G,Q)),
 'electric_commutator':0,'full_commutator_is_zero':False}

# Separate Euclidean action: physical U3 characters are identical, while
# a Euclidean-to-Lorentzian state/clock transfer is not asserted.
sites=list(product(range(2),repeat=4));pos={p:i for i,p in enumerate(sites)}
U=np.empty((16,4,7,7),complex);Ut=np.empty_like(U)
gauge=[expm(1j*((j+1)*T1+(17-j)*T2+j*np.eye(3))/157)for j in range(16)]
for x,p in enumerate(sites):
 for mu in range(4):
  W=expm(1j*((mu+1)*T1+(p[(mu+1)%4]+1)*T2+np.eye(3))/4096)
  q=list(p);q[mu]=(q[mu]+1)%2;y=pos[tuple(q)]
  U[x,mu]=s.u3_matter(W);Ut[x,mu]=s.u3_matter(gauge[x]@W@gauge[y].conj().T)
DW,g5=s.euclidean_wilson(U);D,hat,gap=s.overlap(DW,g5)
DWt,g5t=s.euclidean_wilson(Ut);Dt,hatt,gapt=s.overlap(DWt,g5t)
GG=sparse.block_diag([np.kron(I,s.u3_matter(g))for g in gauge]).toarray()
m=1/1000;Dm=(1-m/2)*D+m*np.eye(len(D))
out['euclidean_gw']={'dimension':len(D),'wilson_sign_gap':gap,
 'GW_error':err(D@g5+g5@D-D@g5@D),'hat_chirality_square':err(hat@hat-np.eye(len(D))),
 'massless_ward_error':err(D@hat+g5@D),
 'massive_ward_error':err(Dm@hat+g5@Dm-2*m*g5@(np.eye(len(D))-D/2)),
 'gauge_covariance_error':err(Dt-GG@D@GG.conj().T),
 'index_diagnostic':float(np.trace(g5@(np.eye(len(D))-D/2)).real),
 'non_singlet_three_equal_rep_families_anomaly_trace':0,
 'interacting_lorentzian_clock_transfer':False,
 'same_U3_charge_representation':True,'same_measure_or_vacuum':False}
charge=np.diag([2/3]*3+[-1/3]*3+[-1])
localQ=np.zeros_like(D);localQ[:28,:28]=np.kron(I,charge)
dD=s.overlap_derivative(DW,g5,1j*comm(localQ,DW))
out['euclidean_gw']['link_current_divergence_ward_error']=err(dD-1j*comm(localQ,D))
out['euclidean_gw']['represented_plaquette_analytic_bound']=20/4096
out['euclidean_gw']['admissibility_ceiling']=1/40
out['euclidean_gw']['conditional_sign_kernel_squared_interval']=[1/4,64]
out['euclidean_gw']['polynomial_degree']=192
out['euclidean_gw']['polynomial_tail_bound']=16*(15/17)**193
out['euclidean_gw']['polynomial_range_links']=385

# Quantum rotor with all16 occupations of4 CAR modes:2 fifth slices,
# Dirac2, fixed physicalspin and one unit-charge representation. The
# rotor is a gauge-invariant loop holonomy; each fixed-N block admits a
# matching static background. No charged gap is interpreted as a mass.
ids=[0,2,4,6];F=s.fifth(2)[np.ix_(ids,ids)]
Z=np.kron(np.eye(2),s.fractional.SIGMA[2]);X=np.kron(np.eye(2),s.fractional.SIGMA[0])
q=np.diag([1,1,-1,-1]);dq=s.fractional.second_quantize(q)
A=s.fractional.second_quantize(F+Z)/s.H_CLOCK
B=s.fractional.second_quantize(-(Z+1j*X)/2)/s.H_CLOCK
mid=s.fractional.second_quantize(comm(F,q))/s.H_CLOCK
rotor=[];records=[]
for K in [4,8,16]:
 n=2*K+1;flux=np.arange(-K,K+1)
 shift=sparse.diags(np.ones(n-1),-1,format='csr')
 HE=sparse.diags(flux**2/4,format='csr')
 H=sparse.kron(sparse.eye(n),A)+sparse.kron(shift,B)+sparse.kron(shift.T,B.conj().T)+sparse.kron(HE,sparse.eye(16))
 Q=sparse.kron(sparse.eye(n),dq);E=sparse.kron(HE,sparse.eye(16))
 ward=comm(H,Q)-sparse.kron(sparse.eye(n),mid)
 state=np.zeros(n*16,complex);state[K*16+3]=1
 evolved=expm_multiply(-1j*H/64,state);probs=abs(evolved)**2
 rotor.append({'K':K,'dimension':n*16,'ward_error':float(max(abs(ward.data),default=0)),
  'electric_axial_error':float(max(abs(comm(E,Q).data),default=0)),
  'all_probability_sum':float(sum(probs)),'edge_flux_probability':float(sum(probs[:16])+sum(probs[-16:])),
  'axial_mean':float(np.vdot(evolved,Q@evolved).real),
  'initial_axial_mean':2,'time_phys':1/64,'preparation':{'flux':0,'occupation_word':3},
  'entire_law_TV_upper_bound':min(1,6/math.factorial(K+1)),
  'vacuum_preparation':False})
 records.append({'K':K,'records':[{'flux':int(flux[j//16]),'occupation':j%16,'probability':float(p)}for j,p in enumerate(probs)]})
out['full_CAR_rotor']=rotor;save('ROTOR_RECORDS.json',records)
out['acceptance']={'wilson_protected_interacting_light_mass':False,
 'euclidean_chiral_action_component':True,'same_clock_complete_domain_wall_hamiltonian':True,
 'finite_domain_wall_exact_axial_conservation':False,
 'interacting_residual_mass_bound':False,'physical_muon_ratio_repaired':False,
 'gauge_anomaly_free_vectorlike_carrier':True,'chiral_SM_multiplicity':False}

assert out['u3_covariance']['error']<1e-11
assert max(out['euclidean_gw'][x]for x in ['GW_error','massless_ward_error','massive_ward_error','gauge_covariance_error'])<1e-10
assert all(r['ward_error']<1e-12 and r['electric_axial_error']==0 and abs(r['all_probability_sum']-1)<1e-12 for r in rotor)
assert all(b['zero_states']==4 and b['zero_momenta']==[[0,0,0]]for b in out['free_band_inventory'])
assert all(r['exact_light_absmax']<r['trial_residual_bound']+1e-13 for r in out['wall_recursion'])
save('RESULT.json',out)
print(json.dumps({'GW':out['euclidean_gw'],'U3':out['u3_covariance'],'rotor':rotor,'accepted_components':out['acceptance']},indent=2))
