from pathlib import Path
import hashlib,json
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh,expm
from model.openwave_dynamic_wave.source import shape,geometry,block,field,packet,KAPPA
from model.openwave_dynamic_stress.response import operators
ROOT=Path.cwd(); HERE=ROOT/'model/openwave_current_averaging/audits'; HERE.mkdir(exist_ok=True)
paths=['model/openwave_dynamic_wave/results.json','model/openwave_dynamic_stress/results.json']
receipts=[json.loads((ROOT/p).read_text())for p in paths];pins={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths}
pins[str(Path(__file__).resolve().relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
for r in receipts:
 for p,h in r['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
s=shape();g=geometry();f=field();o=operators();L=s['L_exact'];W=s['F_exact'];sites=s['source_sites']
assert W[sites,:].det()==1
assert W.T*L*W==s['S_exact'] and W.T*W==s['M_exact']
assert all((L*W)[i,j]==0 for i in range(49)if i not in sites for j in range(2))
# Independent derivative from central differences of the actual source roots.
def root(q):
 d=np.exp(2*q);A=d[:,None]*s['L']*d[None,:];e,V=eigh(A);return (V*np.sqrt(e))@V.T
errors=[]
for a in range(2):
 D=s['F'][:,a];h=2**-12
 for q,R in zip(f['qprofiles'],f['coordinate_blocks']):
  derivative=(root(q+h*D)-root(q-h*D))/(2*h)
  # Frechet derivative computed spectrally, independently of solve_sylvester.
  ev,V=eigh(R);K=R@R;rhs=2*(D[:,None]*K+K*D[None,:]);G=V@((V.T@rhs@V)/(ev[:,None]+ev[None,:]))@V.T
  errors.append(float(np.linalg.norm(derivative-G)))
assert max(errors)<2e-6
# Exact finite commutator boundary must be retained, including a nonzero control.
ward=[]
for a in range(2):
 identity=o['actual_momentum_rate'][a]-o['spring'][a]+o['native_force'][a]+o['field_force'][a]-o['finite_boundary_correction'][a]
 ward.append(float(np.linalg.norm(identity)));assert ward[-1]<1e-10
 assert np.linalg.norm(o['finite_boundary_correction'][a])>.01
# Fresh full-block diagonalization and complete sequential instruments.
means=[];arrival=[];energyerrors=[]
for level in (0,1):
 H=block(left=level)['wave'];e,V=eigh(H);psi=np.kron(g['vacuum'],packet());v=V@(np.exp(-4j*e)*(V.T@psi));mean=float(abs(v.reshape(4,49)) .__pow__(2).sum(axis=0)@np.tile(np.arange(7),7));means.append(mean)
 rho=v.reshape(4,49)@v.reshape(4,49).conj().T
 if level==1:assert np.trace(rho@rho).real<1
 energyerrors.append(float(abs(np.vdot(v,H@v)-np.vdot(psi,H@psi))))
 U=(V*np.exp(-.5j*e))@V.T;mask=np.tile(np.repeat(np.arange(7)>=5,7),4);v=psi.copy();hit=[]
 for step in range(16):v=U@v;hit.append(float(np.vdot(v[mask],v[mask]).real));v[mask]=0
 assert abs(sum(hit)+np.vdot(v,v).real-1)<1e-11
 arrival.append(sum((i+1)*.5*p for i,p in enumerate(hit))/sum(hit))
wave=receipts[0];delay=arrival[1]-arrival[0];bend=means[1]-means[0]
assert wave['same_source_operational_delay_interval'][0]<delay<wave['same_source_operational_delay_interval'][1]
a=wave['same_source_bending_intervals'][3]['source_minus_zero_native_control_shift_interval'];assert a[0]<bend<a[1]
# Every native ground block and seven ratios, without importing the interval producer.
stress=receipts[1]['dressed_clocks'];ground=np.zeros((8,8))
for c in stress['all64_native_vacuum_block_certificates']:
 i,j=c['native_levels'];ev=np.linalg.eigvalsh(block(left=i,right=j)['vacuum']);ground[i,j]=ev[0]
 for x,iv in zip(ev,c['all_eigenvalue_intervals']):assert float(F(iv[0]))-2e-14<=x<=float(F(iv[1]))+2e-14
ratios=[]
for r in stress['all7_same_source_clock_ratios']:
 j=r['native_probe_level'];ratio=(ground[1,j]-ground[1,0])/(ground[0,j]-ground[0,0]);ratios.append(ratio)
 lo,hi=map(F,r['same_source_clock_frequency_ratio_interval']);assert float(lo)-2e-14<ratio<float(hi)+2e-14
assert max(ratios)-min(ratios)>87e-6
assert F(stress['exact_test_clock_theorem']['coefficient_exact'])==F(17,250)**2/F(9,50)
result={'verdict':'PASS_INDEPENDENT_AUTONOMOUS_WAVE_STRESS_REVIEW','source_sha256':pins,'harmonic_projection_and_induced_mass_stiffness_exact':True,'all64_native_clock_blocks_checked':64,'full196_wave_blocks_recomputed':2,'finite_difference_force_error_max':max(errors),'boundary_corrected_Ehrenfest_residuals':ward,'complete_hit_timeout_records_retained':True,'independent_native_delay':delay,'independent_transverse_bending':bend,'energy_error_max':max(energyerrors),'clock_ratio_spread_ppm':(max(ratios)-min(ratios))*1e6,'scope':'Independent numerical controls of exact source identities and positive interval margins. Harmonic geometry, scalar metric, cap and preparation remain declared; no physical G, Einstein law, continuum, observed lensing or universal finite-strength redshift inferred.'}
out=HERE/'gravity_review.json';out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
