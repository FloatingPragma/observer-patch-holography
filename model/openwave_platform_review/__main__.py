"""Full sparse evolution, independently of the portable Fourier evaluator."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as s
from scipy.sparse import csr_matrix,eye,kron
from scipy.sparse.linalg import expm_multiply
from model.openwave_lepton_mass_shell.spectrum import native_frames
from model.lepton_mass_compiler.analyzer import rational_effect_intervals
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];TARGET=ROOT/'model/openwave_platform_kernel'

def pins():
 seeds=set(SOURCE_FILES)
 seeds.update(str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md'))
 seeds.update(str(p.relative_to(ROOT))for p in TARGET.iterdir()if p.suffix in('.py','.md','.tex')or p.name in('source_packet.json','results.json'))
 return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}

def run():
 before=pins();packet=json.loads((TARGET/'source_packet.json').read_text());receipt=json.loads((TARGET/'results.json').read_text())
 for p,h in receipt['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 assert receipt['full_native_execution']['source_packet_sha256']==hashlib.sha256((TARGET/'source_packet.json').read_bytes()).hexdigest()
 native=native_frames();intervals=rational_effect_intervals()['effect_intervals_k0_k1_k2'];masses=[F(a)+F(b)for a,b in intervals];masses=[m/2 for m in masses]
 substitution={symbol:s.Rational(m.numerator,m.denominator)for symbol,m in zip(native['mass_symbols'],masses)}
 H=native['H'].subs(substitution);origin={p:0 for p in native['momenta']}
 exact=[H.subs(origin)]+[H.diff(p)for p in native['momenta']];matrices=[]
 for index,matrix in enumerate(exact):
  numeric=np.asarray(matrix,complex);exported=np.zeros((180,180),complex)
  for i,j,re,im in packet['native_H0_and_three_kinetic_sparse_matrices'][index]:exported[i,j]=re+1j*im
  assert np.max(abs(numeric-exported))<2e-16
  matrices.append(numeric)
 for i in range(1,4):
  assert exact[i]*exact[0]+exact[0]*exact[i]==s.zeros(180)
  for j in range(1,4):assert exact[i]*exact[j]+exact[j]*exact[i]==(2*s.eye(180)if i==j else s.zeros(180))
 x=s.symbols('x',real=True);polynomial=1-s.I*x-x*x/2+s.I*x**3/6+x**4/24
 assert s.expand(polynomial*s.conjugate(polynomial))==1-x**6/72+x**8/576
 results=[]
 for record in receipt['full_native_execution']['runs']:
  side=record['side'];volume=side**3;shape=(side,)*3;action=kron(eye(volume,format='csr'),csr_matrix(matrices[0]),format='csr')
  for axis in range(3):
   rows=[];cols=[]
   for address in np.ndindex(shape):
    neighbor=list(address);neighbor[axis]=(neighbor[axis]+1)%side
    rows.append(np.ravel_multi_index(address,shape));cols.append(np.ravel_multi_index(tuple(neighbor),shape))
   shift=csr_matrix((np.ones(volume),(rows,cols)),shape=(volume,volume));difference=(shift-shift.T)/(2j)
   action+=kron(difference,csr_matrix(matrices[axis+1]),format='csr')
  assert (action-action.conj().T).nnz==0
  initial=np.zeros((volume*180,3),complex)
  for family,frame in enumerate(native['frames']):
   seed=np.asarray(frame[:,0],complex).ravel();momentum=np.zeros((*shape,180),complex)
   for address in np.ndindex(shape):
    p=np.sin(2*np.pi*np.array(address)/side);p[np.abs(p)<1e-14]=0;symbol=matrices[0]+sum(q*g for q,g in zip(p,matrices[1:]));energy=np.sqrt(p@p+float(masses[family])**2)
    momentum[address]=(seed+symbol@seed/energy)/(2*np.sqrt(volume))
   initial[:,family]=np.fft.ifftn(momentum,axes=(0,1,2),norm='ortho').ravel()
  evolved=expm_multiply(-1j*action,initial,traceA=0)
  assert np.max(abs(np.sum(abs(evolved)**2,axis=0)-.5))<2e-14, np.sum(abs(evolved)**2,axis=0)
  maximum=0.;tv=0.;minimum_total=2.;maximum_total=0.
  for case in record['cases']:
   family=case['family'];reference=abs(evolved[:,family])**2;saved=np.asarray(case['exact_source_coordinate_records']);approx=np.asarray(case['normalized_numerical_records'])
   assert np.isfinite(saved).all()and np.isfinite(approx).all()and min(saved)>=0 and min(approx)>=0
   difference=float(np.max(abs(saved-reference)));total=float(np.sum(abs(saved-reference))/2)
   assert difference<2e-13 and total<1e-12
   maximum=max(maximum,difference);tv=max(tv,total);minimum_total=min(minimum_total,float(saved.sum()+.5));maximum_total=max(maximum_total,float(saved.sum()+.5))
   bound=F(2**5,120*case['steps']**4);assert abs(float(bound)-case['analytic_RK4_operator_error_bound'])<1e-18
   for key in('vector_error','normalization_drift','complete_detector_TV'):
    assert np.isfinite(case[key])
   assert case['vector_error']<float(bound)+2e-12 and case['complete_detector_TV']<2*float(bound)+2e-12
   assert abs(case['normalization_drift'])<2*float(bound)+2e-12
  assert record['full_carrier_derivative_error']<3e-14 and record['kinetic_sign_mutation_difference']>.1
  results.append({'side':side,'independently_evolved_full_Nambu_rows':action.shape[0],
    'independent_sparse_action_nonzero_entries':action.nnz,'charged_preparations':3,
    'maximum_independent_every_coordinate_record_difference':maximum,
    'maximum_independent_complete_record_TV':tv,'complete_original_probability_total_range':[minimum_total,maximum_total]})
 assert pins()==before
 return {'schema':'oph.independent_full_native_platform_review.v1','verdict':'PASS_INDEPENDENT_FULL_SOURCE_PLATFORM_REVIEW',
  'source_sha256':before,'exact_native_symbol_and_Clifford_checks':True,'exact_RK4_norm_polynomial_verified':True,
  'independent_method':'Full periodic sparse180-row action and scipy expm_multiply, independently of exported charged-frame Fourier evolution',
  'runs':results,'complete_preparation_failure_and_full_complement_checked':True,
  'roundoff_scope':'Independent numerical records are diagnostics. The analytic integration bound is exact arithmetic, not a hardware rounding enclosure.',
  'scope':'Same finite free source with retained Nambu conjugates, sine valleys and massless coefficients; no interacting-spectrum or physical-clock transfer.'}

def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();v=run();path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
 else:assert equivalent(json.loads(path.read_text()),v),'independent platform review differs'
 print(v['verdict'])
if __name__=='__main__':main()
