"""Distinct positive canonical trace and exact count-consumer review."""
from pathlib import Path
import json,math
import numpy as np
import mpmath as mp
ROOT=Path(__file__).parent;mp.iv.dps=50
r=json.loads((ROOT/'RESULT.json').read_text());data=np.load(ROOT/'COMPLETE_TRACE_INPUTS.npz')
bands=data['one_particle_bands'];occ=data['N8_occupation_indices'];words=data['N8_occupation_words'];traces=data['pointwise_scaled_trace_upper'];betas=data['betas']
assert bands.shape==(1024,16) and occ.shape==(12870,8)
assert len(set(map(int,words)))==12870 and all(int(w).bit_count()==8 for w in words)
assert all(sum(1<<int(i)for i in row)==int(w)for row,w in zip(occ,words))
points={x['grid_index']:x for x in json.loads((ROOT/'inputs/pointwise_controls.json').read_text())if x['N']==8}
anchor=r['trace_anchor_energy'];maximum=0
def elementary(weights):
 coeff=np.zeros(9);coeff[0]=1
 for w in weights:
  for degree in range(8,0,-1):coeff[degree]+=w*coeff[degree-1]
 return coeff
rebuilt=np.zeros_like(traces)
for j,e in enumerate(bands):
 assert np.max(abs(e+e[::-1]))<1e-11
 e0=float(sum(e[:8]))
 for ib,beta in enumerate(betas):
  holes=elementary(np.exp(beta*e[:8]));particles=elementary(np.exp(-beta*e[8:]))
  # Sum all r-hole/r-particle sectors r1..8. The binomial identity
  # sum_r choose(8,r)^2=12870 includes every neutral Fock occupation.
  excited=float(sum(holes[k]*particles[k]for k in range(1,9)))
  value=(math.exp(beta*(anchor-points[j]['improved_pointwise_lower']))+
         math.exp(beta*(anchor-e0))*excited)*math.exp(beta*1e-8)*(1+1e-10)
  rebuilt[j,ib]=value
  maximum=max(maximum,abs(value-traces[j,ib]))
assert maximum<1e-8
assert sum(math.comb(8,k)**2 for k in range(9))==12870

# Every actual Gauss mean shift is covered; the proof applies to all
# real shifts. This numerical check retains the complete shift set.
shifts=set()
for w in map(int,words):
 q=[(w&((1<<(4*(x+1)))-1)).bit_count()-2*(x+1)for x in range(4)]
 shifts.add(sum(q)/4)
maxratio=0
for beta in betas:
 center=sum(math.exp(-beta*l*l)for l in range(-64,65))
 for q in shifts:
  shifted=sum(math.exp(-beta*(l+q)**2)for l in range(-64,65))
  maxratio=max(maxratio,shifted/center)
assert maxratio<1+1e-12

for ib,base in enumerate(r['all_heat_trace_bounds']):
 beta=mp.iv.mpf(base['beta'])
 rotor=1+2*sum((mp.iv.exp(-beta*n*n)for n in range(1,17)),mp.iv.mpf(0))
 rotor+=mp.iv.exp(-beta*256)/(beta*16)
 assert base['rotor_heat_trace_upper']>=float(rotor.b)
 q=mp.iv.exp(beta*16*mp.iv.sqrt(2)*mp.iv.pi/1024)
 assert base['relative_quadrature_upper']>=float(q.b)
 recomputed=base['rotor_heat_trace_upper']*base['relative_quadrature_upper']*float(np.mean(rebuilt[:,ib]))*(1+1e-10)*(1+1e-9)
 assert abs(recomputed-base['scaled_at_anchor_trace_upper'])<1e-8

prepared=json.loads((ROOT/'inputs/prepared_results.json').read_text());pmap={p['case']:p for p in prepared}
coulomb=json.loads((ROOT/'inputs/coulomb_results.json').read_text())
cmap={p['case']:p for p in coulomb['prepared_response_consumers']}
charged_lower=next(x['complete_source_lower']for x in coulomb['every_number_sector']if x['N']==7)
for test in r['all_eigen_count_tests']:
 p=pmap[test['case']];base=next(x for x in r['all_heat_trace_bounds']if x['beta']==test['beta'])
 expected=base['scaled_at_anchor_trace_upper']*math.exp(test['beta']*(test['cut']-anchor+4*abs(p['mass'])))*(1+1e-10)
 assert abs(expected-test['eigenvalue_count_upper'])<1e-10
 assert test['exactly_one_below_cut_certified']==(expected<2 and test['all_other_number_sectors_excluded'])
for consumer in r['vacuum_consumers']:
 assert consumer['unique_actual_vacuum_certified']
 w=consumer['count_witness'];p=pmap[consumer['case']]
 assert 1<w['eigenvalue_count_upper']<2 and w['cut']>p['variational_energy_upper']
 sigma=p['full_preparation_spectral_residual']+1e-8;delta=w['offset']
 weight=min(1,(sigma/delta)**2)
 assert abs(consumer['true_ground_fidelity_lower']-(1-weight))<1e-14
 ground_lower=p['variational_energy']-1e-6-sigma*sigma/delta
 assert abs(consumer['actual_ground_energy_interval'][0]-ground_lower)<1e-12
 assert consumer['actual_ground_energy_interval'][1]==p['variational_energy_upper']
 gap=charged_lower-4*abs(p['mass'])-p['variational_energy_upper']
 state_error=2*math.sqrt(weight)/gap
 energy_error=(2e-6+sigma*sigma/delta)/(gap*gap)
 total_error=cmap[p['case']]['prepared_Green_error_upper']+state_error+energy_error
 assert abs(consumer['charged_resolvent_gap_lower']-gap)<1e-12
 assert abs(consumer['ground_state_response_error_upper']-state_error)<1e-12
 assert abs(consumer['ground_energy_response_error_upper']-energy_error)<1e-12
 assert abs(consumer['actual_vacuum_Green_error_upper']-total_error)<1e-12
 G=np.array(p['Green_real'])+1j*np.array(p['Green_imag'])
 invnorm=1/min(np.linalg.svd(G,compute_uv=False))
 if invnorm*total_error<1:
  inverse_error=invnorm*invnorm*total_error/(1-invnorm*total_error)
  assert abs(consumer['actual_vacuum_inverse_error_upper']-inverse_error)<1e-9
 else:assert consumer['actual_vacuum_inverse_error_upper']is None
 if consumer['case'].endswith('K4'):
  assert consumer['true_ground_fidelity_lower']>.9999977
  assert consumer['actual_vacuum_Green_error_upper']<.001046
 assert not consumer['physical_lepton_pole_claim']
review={'producer_source_imported':False,'all_positive_particle_hole_trace_controls':1024*8,
 'complete_Fock_occupations_per_control':12870,'trace_reconstruction_max_error':maximum,
 'all_Gauss_mean_shifts':sorted(shifts),'shifted_kernel_max_ratio':maxratio,
 'outward_rotor_and_holonomy_bounds_pass':True,'all160_count_tests_rechecked':len(r['all_eigen_count_tests'])==160,
 'all_vacuum_fidelity_and_current_consumers_pass':True,'physical_mass_not_promoted':True}
(ROOT/'REVIEW.json').write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps(review,indent=2))
