import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from model.openwave_six_quark.source import reduction,MODES
from model.openwave_six_quark.__main__ import pins as nuclear_pins
from model.openwave_baryon_ordering.__main__ import pins as baryon_pins
from model.particle_common_contract.provenance import import_closure
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def partitions(n,minimum=1):
 yield (n,)
 for i in range(minimum,n//2+1):
  for rest in partitions(n-i,i):yield (i,)+rest

def certificate():
 d=reduction();path=ROOT/'model/openwave_six_quark/results.json';r=json.loads(path.read_text());bpath=ROOT/'model/openwave_baryon_ordering/results.json';b=json.loads(bpath.read_text())
 assert r['source_sha256']==nuclear_pins() and b['source_sha256']==baryon_pins()
 t=r['positive_bound_state'];v=list(map(int,t['full_Cartan_sector_integer_coefficients']));norm=sum(x*x for x in v);assert norm==int(t['exact_trial_squared_norm'])
 moment=F(0);charge=F(0);N6=0;fixed_colour=0
 for index,amp in zip(d['keep'],v):
  mask=d['basis'][int(index)];occ=[m for i,m in enumerate(MODES)if mask>>i&1];assert len(occ)==6
  q=[F(2,3)if m[2]==0 else -F(1,3)for m in occ];X=F(sum(m[0]for m in occ),6)
  assert sum(q)==1 and sum(1 if m[1]==0 else -1 for m in occ)==2
  assert all(sum(m[3]==c for m in occ)==2 for c in range(3))
  p=F(amp*amp,norm);moment+=p*2*sum(qi*(F(m[0])-X)**2 for qi,m in zip(q,occ));charge+=p*sum(q);N6+=p*len(occ)
 assert charge==1 and N6==6
 assert moment==F(t['exact_trial_intrinsic_axial_quadrupole'])
 assert F(t['ground_quadrupole_interval'][0])>0
 assert F(t['all_ground_vectors_quadrupole_lower_from_binding'])>0
 trials=F(t['exact_Rayleigh_energy']);lowerB=F(b['exact_equal_mass_spectrum_certificate']['lowest_interval_rational'][0]);rows=[]
 for part in partitions(6):
  if len(part)==1:continue
  if part==(3,3):floor=2*lowerB
  else:
   nonzero=sum(n%3!=0 for n in part);assert nonzero>=2;floor=F(6)+F(4,3)*nonzero
  assert floor>trials;rows.append({'partition':list(part),'threshold_lower':str(floor),'trial_binding_margin_lower':str(floor-trials)})
 # A globally singlet pair of non-singlet three-quark fragments requires two
 # conjugate nontrivial colour representations. Each fragment costs >=8/3
 # beyond number energies: local3+0 costs >=3; local2+1 costs >=8/3.
 nonsinglet33=F(6)+F(16,3);assert nonsinglet33>trials
 assert r['additional_degenerate_partner'] and r['conditional_common_electron_unit']['physical_deuteron_mass_or_binding_identified']is False
 assert F(r['unbalanced_S1_I0_ground_above_two_baryons_lower'])>0
 success=t['complete_success_occupation_probabilities'];failure=t['complete_failure_occupation_probabilities'];assert len(success)==len(failure)==11700
 assert min(success)>=0 and min(failure)>=0 and abs(sum(success)+sum(failure)-1)<1e-12
 seed=[str(p.relative_to(ROOT))for p in HERE.glob('*.py')];sources=set(import_closure(seed));sources.update(r['source_sha256']);sources.update(b['source_sha256']);sources.update((str(path.relative_to(ROOT)),str(bpath.relative_to(ROOT))))
 pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(sources)}
 return {'verdict':'PASS_INDEPENDENT_COMPLETE_SIX_QUARK_RECORDS_AND_FRAGMENT_REVIEW','source_sha256':pins,
 'reviewed_receipt_sha256':{str(path.relative_to(ROOT)):pins[str(path.relative_to(ROOT))],str(bpath.relative_to(ROOT)):pins[str(bpath.relative_to(ROOT))]},
 'independent_native_occupation_charge_exact':str(charge),'independent_particle_number_exact':str(N6),'independent_direct_native_quadrupole_exact':str(moment),
 'every_nontrivial_partition':rows,'noncolour_singlet_three_plus_three_threshold_lower':str(nonsinglet33),
 'complete_success_and_failure_native_records':23400,'unbalanced_no_weak_cost_control_retained':True,'extra_spin_zero_isotriplet_retained':True,
 'review_scope':'Exact occupation-based charge and number-centre quadrupole, all fragment partitions including non-singlet3+3, complete records, unchanged scale and architecture ancestry; producer supplies exact27-state projector and integer spectral residual witnesses.',
 'physical_deuteron_identity_claimed':False}

def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=certificate();path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:assert json.loads(path.read_text())==r
 print(r['verdict'])
if __name__=='__main__':main()
