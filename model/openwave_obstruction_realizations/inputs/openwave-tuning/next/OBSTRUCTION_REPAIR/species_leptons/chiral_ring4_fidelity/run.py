"""All-Fock heat trace, exact eigen-count criterion and vacuum consumers."""
from pathlib import Path
from itertools import combinations
import importlib.util,json,math,time
import numpy as np
import mpmath as mp
ROOT=Path(__file__).parent;GRID=1024
BETAS=[.25,.5,.75,1.,1.25,1.5,2.,3.];OFFSETS=[1/16,1/8,1/4,1/2,1]
mp.iv.dps=40
spec=importlib.util.spec_from_file_location('frozen_source',ROOT/'inputs/source.py')
src=importlib.util.module_from_spec(spec);spec.loader.exec_module(src)
prepared=json.loads((ROOT/'inputs/prepared_results.json').read_text())
points=json.loads((ROOT/'inputs/pointwise_controls.json').read_text())
point={r['grid_index']:r for r in points if r['N']==8}
coulomb=json.loads((ROOT/'inputs/coulomb_results.json').read_text())
anchor=next(r['variational_energy_upper']for r in prepared if r['case']=='zero_K4')
occ=np.array(list(combinations(range(16),8)),dtype=int)
words=np.array([sum(1<<int(i)for i in row)for row in occ],dtype=np.int32)
assert len(words)==12870
A,C,_=src.one_particle(0);bands=[];traces=[];started=time.time()
for j in range(GRID):
 theta=2*np.pi*j/GRID
 H=A+C*np.exp(1j*theta)+C.conj().T*np.exp(-1j*theta)
 values=np.linalg.eigvalsh(H)
 energies=values[occ].sum(axis=1)
 # Omitting any one free state is at least as conservative as omitting
 # the exact free ground when its Boltzmann weight is replaced by the
 # lower bound on the interacting pointwise ground.
 removed=int(np.argmin(energies));excited=np.delete(energies,removed)
 ground_lower=point[j]['improved_pointwise_lower']-1e-8
 row=[]
 for beta in BETAS:
  excited_upper=float(np.sum(np.exp(beta*(anchor-excited+1e-8))))
  ground_upper=math.exp(beta*(anchor-ground_lower))
  row.append((ground_upper+excited_upper)*(1+1e-10))
 traces.append(row);bands.append(values)
 if j%256==0:print('grid',j,'seconds',time.time()-started,flush=True)
traces=np.array(traces);bands=np.array(bands)
np.savez_compressed(ROOT/'COMPLETE_TRACE_INPUTS.npz',one_particle_bands=bands,N8_occupation_indices=occ,
 N8_occupation_words=words,pointwise_scaled_trace_upper=traces,betas=np.array(BETAS))
trace_bounds=[]
for ib,beta in enumerate(BETAS):
 b=mp.iv.mpf(beta)
 rotor=1+2*sum((mp.iv.exp(-b*n*n)for n in range(1,17)),mp.iv.mpf(0))
 rotor+=mp.iv.exp(-b*16*16)/(b*16)
 rotor_upper=np.nextafter(float(rotor.b),np.inf)
 quadrature=mp.iv.exp(b*16*mp.iv.sqrt(2)*mp.iv.pi/GRID)
 quadrature_upper=np.nextafter(float(quadrature.b),np.inf)
 mean=float(np.mean(traces[:,ib]))*(1+1e-10)
 scaled=float(rotor_upper*quadrature_upper*mean*(1+1e-9))
 trace_bounds.append({'beta':beta,'rotor_heat_trace_upper':float(rotor_upper),
  'relative_quadrature_upper':float(quadrature_upper),'mean_pointwise_scaled_trace_upper':mean,
  'scaled_at_anchor_trace_upper':scaled})

alltests=[];consumers=[]
for case in prepared:
 m=case['mass'];mass_shift=4*abs(m)
 source_lower={r['N']:r['complete_source_lower']-(0 if r['N']in[0,16]else mass_shift)for r in coulomb['every_number_sector']}
 admitted=[]
 for offset in OFFSETS:
  cut=case['variational_energy_upper']+offset
  others_excluded=all(v>cut for n,v in source_lower.items()if n!=8)
  for base in trace_bounds:
   beta=base['beta'];count=float(base['scaled_at_anchor_trace_upper']*math.exp(beta*(cut-anchor+mass_shift))*(1+1e-10))
   # A normalized trial lies below this cut, so a valid trace upper
   # bound must admit at least one eigenvalue.
   assert count>1-1e-8
   passed=bool(count<2 and others_excluded)
   row={'case':case['case'],'beta':beta,'offset':offset,'cut':cut,
    'eigenvalue_count_upper':count,'all_other_number_sectors_excluded':others_excluded,
    'exactly_one_below_cut_certified':passed}
   alltests.append(row)
   if passed:admitted.append(row)
 consumer={'case':case['case'],'unique_actual_vacuum_certified':bool(admitted),
  'source_Hamiltonian_changed':False,'physical_lepton_pole_claim':False}
 if admitted:
  best=max(admitted,key=lambda r:(r['offset'],-r['eigenvalue_count_upper']))
  delta=best['offset'];sigma=case['full_preparation_spectral_residual']+1e-8
  excited_weight=min(1,(sigma/delta)**2)
  energy_error=2e-6+sigma*sigma/delta
  lower_E=case['variational_energy']-1e-6-sigma*sigma/delta
  actual_gap=source_lower[7]-case['variational_energy_upper']
  prepared_receipt=next(x for x in coulomb['prepared_response_consumers']if x['case']==case['case'])
  prep_error=prepared_receipt['prepared_Green_error_upper']
  state_error=2*math.sqrt(excited_weight)/actual_gap
  energy_response_error=energy_error/(actual_gap**2)
  total_error=prep_error+state_error+energy_response_error
  G=np.array(case['Green_real'])+1j*np.array(case['Green_imag'])
  invnorm=float(np.linalg.norm(np.linalg.inv(G),2));inverse_error=None
  if invnorm*total_error<1:inverse_error=invnorm**2*total_error/(1-invnorm*total_error)
  consumer.update({'count_witness':best,'neutral_spectral_gap_lower':delta,
   'actual_ground_energy_interval':[lower_E,case['variational_energy_upper']],
   'true_ground_fidelity_lower':1-excited_weight,'prepared_full_residual_upper':sigma,
   'charged_resolvent_gap_lower':actual_gap,'prepared_Green_error_upper':prep_error,
   'ground_state_response_error_upper':state_error,'ground_energy_response_error_upper':energy_response_error,
   'actual_vacuum_Green_error_upper':total_error,'actual_vacuum_inverse_error_upper':inverse_error,
   'finite_volume_true_vacuum_propagator_certified':True})
 consumers.append(consumer)

result={'unchanged_source':'9cbfad0966962e59ad6424a59925cae61ca53b7f4dbb553e62d37228d507d0ff',
 'coulomb_bound_source':'deed27e01f002d986253614d8dca0f1d5bc1e153de04f2b2b139440007dd9135',
 'neutral_occupations_retained':12870,'holonomy_grid':1024,'complete_free_fock_energies_evaluated':12870*1024,
 'all_one_particle_bands_per_grid':16,'trace_anchor_energy':anchor,
 'all_heat_trace_bounds':trace_bounds,'all_eigen_count_tests':alltests,'vacuum_consumers':consumers,
 'numerical_energy_envelope':1e-8,'source_parameters_fitted':False,'elapsed_seconds':time.time()-started}
(ROOT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'trace_bounds':trace_bounds,'vacuum_consumers':consumers,'elapsed':time.time()-started},indent=2))
