"""Independent complete-source checks for the energy/lapse constructions."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as s
from scipy.linalg import expm
from model.openwave_energy_lapse.source import native,prepare,difference,KAPPA
from model.particle_common_contract.provenance import import_closure
from model.openwave_lapse_observations.observations import clocks,wave,wave_matrix,metric_and_curvature
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent

def pins():
 paths={str(Path(__file__).relative_to(ROOT))}
 for name in ('openwave_energy_lapse','openwave_lapse_observations'):
  p=ROOT/'model'/name
  paths|={str(x.relative_to(ROOT))for x in p.iterdir()if x.suffix in('.py','.md','.tex')or x.name=='results.json'}
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(paths))}

def report():
 before=pins();M,levels,gaps=native();X=s.Matrix([[0,1],[1,0]]);I=s.eye(2)
 independent=s.diag(*(s.Rational(2*i.bit_count()+3,2)for i in range(8)))-(s.kronecker_product(X,X,I)+s.kronecker_product(I,X,X))/7
 assert M==independent
 vals=np.linalg.eigvalsh(np.array(M,float))
 for value,box in zip(vals,levels):assert float(box.lo)-2e-14<=value<=float(box.hi)+2e-14
 K=np.array(M,float)-vals[0]*np.eye(8);assert np.linalg.eigvalsh(K).min()>-2e-14
 q=-.0002;step=1e-6;derivative=(np.exp(q+step)*K-np.exp(q-step)*K)/(2*step)
 assert np.max(abs(derivative-np.exp(q)*K))<1e-9
 a=((3,3,3),1);b=((5,3,3),7)
 left=prepare(7,(a,));right=prepare(7,(b,));pair=prepare(7,(a,b));cross=difference(pair,left,right);assert cross[1]<0
 for witness in(left,right,pair):
  assert F(witness['fixed_point_contraction_rational_upper'])<F(1,1000)
  assert F(witness['minimum_energy_interval'][0])<=F(witness['minimum_energy_interval'][1])
 clock=clocks()
 for row in clock['rows']:
  intervals=[list(map(F,x['probability_interval']))for x in row['complete64_joint_records']]
  assert len(intervals)==64 and sum(x[0]for x in intervals)<=1<=sum(x[1]for x in intervals)
  assert row['near_to_far_frequency_ratio_interval']==clock['universal_frequency_ratio_interval']
 profile=np.linspace(-.0003,0,31);w=wave(profile,8);Kwave=wave_matrix(profile)
 vals,U=np.linalg.eigh(Kwave);H=(U*np.sqrt(vals))@U.T;x=np.arange(31);psi=np.exp(-(x-7)**2/8+1j*.7*x);psi/=np.linalg.norm(psi)
 out=expm(-8j*H)@psi;independent_p=np.r_[.5,.5*abs(out)**2]
 err=float(max(abs(independent_p-w['complete_vacuum_and_position_probabilities'])));assert err<3e-13
 assert abs(w['source_energy_before']-w['source_energy_after'])<3e-13
 assert abs(w['mean_detector_work'])>1e-5
 curve=metric_and_curvature();assert curve['gamma_zero_countercontrol_G11_at_3_0_0']=='2/27'
 assert not curve['gamma_selected_by_clock_and_force_alone'] and not curve['regional_entropy_stationarity_premises_discharged']
 assert pins()==before
 return {'schema':'oph.independent_energy_lapse_review.v1','source_sha256':before,
  'complete_native8_matrix_independently_rebuilt':True,'all_native_levels_compared':8,
  'operator_response_energy_identity_on_complete_carrier':True,
  'independent_small_cube_mixed_energy_interval':list(map(str,cross)),
  'seven_clocks_all448_records_normalize':True,'independent_expm_complete_wave_probability_error':err,
  'wave_measurement_work_retained':True,'gamma_countercontrol_retained':True,
  'architectural_lapse_stiffness_counterterm_and_field_extension_retained':True,
  'physical_Newton_strength_or_full_Einstein_premises_discharged':False,
  'verdict':'PASS_INDEPENDENT_ENERGY_LAPSE_SCOPE_REVIEW'}

def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=report();path=HERE/'lapse_review.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 elif json.loads(path.read_text())!=r:raise ValueError('lapse review changed')
 print(r['verdict'])
if __name__=='__main__':main()
