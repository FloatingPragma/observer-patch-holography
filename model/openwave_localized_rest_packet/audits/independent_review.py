"""Independent finite shifts, spectral inequalities and complete phase records."""
from pathlib import Path
import hashlib,json
from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from model.openwave_localized_rest_packet.packet import one_axis_moments,moment_upper_bounds,finite_execution,certify
ROOT=Path(__file__).resolve().parents[3]
def run():
 maximum=0.;counts=0
 for L in range(4,25):
  for N in(L+5,2*(L+5)):
   shift=np.roll(np.eye(N),1,axis=0);D=(shift-shift.T)/(2j)
   g=np.zeros(N);j=np.arange(1,L+1);g[j]=np.sqrt(8/(3*(L+1)))*np.sin(np.pi*j/(L+1))**2
   d2=float(np.vdot(D@g,D@g).real);d4=float(np.vdot(D@D@g,D@D@g).real)
   a,b=one_axis_moments(L);maximum=max(maximum,abs(d2-a),abs(d4-b));counts+=2
   upper=moment_upper_bounds(L);assert d2<=float(upper['one_axis_second'])and d4<=float(upper['one_axis_fourth'])
 assert maximum<1e-13
 # Reconstruct scalar Fourier spectral probabilities independently, then
 # directly evolve a pointer plus every retained eigenline.
 row=finite_execution(length=4,side=9,time=.7,mass=1.);N=9;L=4
 g=np.zeros(N);j=np.arange(1,L+1);g[j]=np.sqrt(8/(3*(L+1)))*np.sin(np.pi*j/(L+1))**2
 fourier=np.array([sum(g[x]*np.exp(-2j*np.pi*k*x/N)for x in range(N))/np.sqrt(N)for k in range(N)])
 p=np.sin(2*np.pi*np.arange(N)/N);E=np.sqrt(1+p[:,None,None]**2+p[None,:,None]**2+p[None,None,:]**2)
 weights=abs(fourier[:,None,None]*fourier[None,:,None]*fourier[None,None,:])**2*(1+1/E)/2
 w=weights.sum();C=np.sum(weights*np.exp(-1j*.7*E))/w
 assert abs(w-(1-row['original_preparation_failure']))<1e-13
 assert abs(C-complex(*row['conditional_complex_return']))<1e-13
 # Full pure branch and orthogonal-return effects, without the producer law.
 v=np.sqrt(weights.ravel()/w);psi=np.r_[1.,v*np.exp(-1j*.7*E.ravel())]/np.sqrt(2)
 expected=[]
 for phase in(1,-1j):
  for sign in(-1,1):
   arm=np.r_[1.,sign*phase*v]/np.sqrt(2);expected.append(w*abs(np.vdot(arm,psi))**2/2)
  expected.append(w*(1-abs(np.vdot(np.r_[1.,v]/np.sqrt(2),psi))**2-abs(np.vdot(np.r_[1.,-v]/np.sqrt(2),psi))**2)/2)
 expected.append(1-w)
 actual=[x['original_probability']for x in row['complete_clock_records']]
 error=float(np.max(abs(np.array(expected)-actual)));assert error<1e-13
 # Random independent positive spectral laws exercise all energy/phase
 # inequalities, including the denominator of the frequency estimate.
 rng=np.random.default_rng(337);trials=0
 for _ in range(100):
  k2=rng.uniform(0,.0001,17);m=rng.uniform(.5,2);z=rng.random(17);z/=z.sum();E=np.sqrt(m*m+k2)
  w=np.dot(z,(1+m/E)/2);law=z*(1+m/E)/(2*w);d=np.dot(z,k2)/(2*m*w)
  t=1.;C=np.dot(law,np.exp(-1j*t*E));deriv=np.dot(law,-1j*E*np.exp(-1j*t*E));freq=-(deriv/C).imag
  assert abs(C-np.exp(-1j*m*t))<=t*d+1e-14
  assert abs(freq-m)<=d/(1-t*d)+1e-14;trials+=1
 big=certify();assert float(F(big['mean_relative_rest_energy_excess_upper']))<3.014e-9
 assert float(F(big['conditional_outside_cube_uniform_time_upper']))<1.103e-7
 compiler=big['finite_full_Fock_compiler'];volume=big['periodic_side']**3
 # 45 positive/negative Weyl pairs per site, |pi|<=3; six pairing terms.
 norm_bound=2*(45*3+6)*volume
 assert norm_bound==int(compiler['physical_excitation_energy_operator_norm_upper'])
 transferred=norm_bound*F(compiler['full_source_unhalved_trace_bound'])
 assert transferred==F(compiler['conditional_on_initial_success_physical_energy_mean_transfer_error_upper'])
 assert transferred<=F(big['mean_source_energy_excess_upper'])/10**6
 assert transferred<F(big['mean_source_energy_excess_lower'])
 assert big['autonomous_clock']['instantaneous_frequency_derivative_transferred_to_controlled_pulses']is False
 refs=['model/openwave_localized_rest_packet/'+f for f in ['packet.py','PROOF.md','results.json']]
 return {'verdict':'PASS','independent_shift_moments':counts,'maximum_shift_error':maximum,
  'complete_pointer_probability_error':error,'independent_spectral_inequality_trials':trials,
  'physical_marginal_energy_error_budget_checked':True,
  'scope':'Analytic finite-carrier existence and complete small-box records; enormous volume not numerically simulated. The physical energy marginal retains every ancillary output. Autonomy and phase derivatives concern the declared quadratic target; the compiler uses continuing controls and transfers sampled records. Free persistence is not interacting stability; prepared phase is not passive onset.',
  'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in refs}}
if __name__=='__main__':
 v=run();Path(__file__).with_suffix('.json').write_text(json.dumps(v,indent=2,sort_keys=True)+'\n');print(v['verdict'])
