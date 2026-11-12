"""A nonzero infinite-line emission interval with complete source error."""
from fractions import Fraction as F
from math import isqrt
import numpy as np
import sympy as sy
from mpmath import iv
from .window import certified_universal,taylor_first,compare_native,construction,reduction,G,finite_volume_bound
from model.openwave_six_quark.certificate import spectral_enclosure,sqrt_upper_fraction
from model.openwave_quantum_rotational.enclosure import scalar
from model.native_flavour_transport.transport import source,word_unitary

def family_bound():
 c=spectral_enclosure(source()['y']);D=c['dyadic_denominator'];O=sy.Matrix([[sy.Rational(int(v),D)for v in row]for row in c['dyadic_eigenframe']]);z=(O[:,0].T*word_unitary((2,3,10))*O[:,1])[0]
 iv.dps=60;scalar.cache_clear();squares=F(0);components=[]
 for q in(sy.re(z),sy.im(z)):
  x=scalar(sy.expand(q));lo=F.from_float(float(np.nextafter(float(x.a),-np.inf)));hi=F.from_float(float(np.nextafter(float(x.b),np.inf)));a=min(abs(lo),abs(hi))if lo*hi>0 else F(0);squares+=a*a;components.append([str(lo),str(hi)])
 e=[F(x,D)for x in map(int,c['dyadic_eigenvalues'])];gap=min(e[i+1]-e[i]for i in range(2))-2*F(c['all_eigenvalue_error_upper']);error=16*F(c['all_eigenvalue_error_upper'])/gap+8*F(c['exact_Gram_error_upper']);lower=F(19,100)
 assert squares>(lower+error)**2
 return {'actual_native_word':[2,3,10],'absolute_C01_lower':str(lower),'complete_native_family_spectral_certificate':c,'dyadic_matrix_element_component_intervals':components,'eigenframe_amplitude_error_upper':str(error)}

def norm_lower(v):
 D=2**60;a=[int(round(float(x)*D))for x in v.real];b=[int(round(float(x)*D))for x in v.imag];square=sum(x*x+y*y for x,y in zip(a,b));n=isqrt(square);error=F(isqrt(2*len(v))+1,2*D)
 return F(n,D)-error

def finite_control(L,t):
 from scipy.sparse.linalg import expm_multiply
 from .window import first_light,universal
 d=construction(L);B,h,_=reduction(d);a=B.conj().T@d['initial'];state=expm_multiply(-1j*float(t)*h,a);u=universal(L)
 # First order in the block off-diagonal coupling, with full D and all families.
 b=np.zeros(len(state)-6,complex);n=u['dimension']
 for f in range(3):
  for ell in range(3):
   v=dict(u);v['D']=u['D']+(d['mass'][f]+d['mass'][ell]-2)*__import__('scipy').sparse.eye(n)
   b[(3*f+ell)*n:(3*f+ell+1)*n]=d['C'][f,1]*first_light(v,t)
 delta=(F(1,64)*sqrt_upper_fraction(F(24))*F(t))**3/6
 b*=np.exp(-1j*float(t)*u['parent_contact_energy']);error=float(np.linalg.norm(state[6:]-b));assert error<float(delta)*1.0001+1e-13
 return {'ring_sites':L,'source_time':str(t),'complete_prepared_dimension':h.shape[0],'exact_pair_probability_floating':float(np.linalg.norm(state[6:])**2),'full_nine_channel_first_emission_probability_floating':float(np.linalg.norm(b)**2),'actual_vector_difference_floating':error,'proved_cubic_remainder_upper':str(delta),'all_channels_in_comparison':True}

def certificate():
 family=family_bound();rows=[];lmb=F(1,64)*sqrt_upper_fraction(F(24));native_controls=[compare_native(L)for L in(2,4)]
 for L,t in((32,F(1,100)),(128,F(1,10))):
  u=certified_universal(L);v,num=taylor_first(u,t);lower=norm_lower(v)-F(num['total_numeric_vector_error_upper_exact']);cubic=(lmb*t)**3/6;vol=finite_volume_bound(L,t);tail=F(vol['prepared_domain_vector_error_upper_exact']);amplitude=F(19,100)*lower-cubic-tail
  assert amplitude>0;prob=amplitude*amplitude;success=construction(2)['success'];assert success>F(6,25)
  row={'ring_sites':L,'source_time':str(t),'complete_prepared_dimension':6+36*L*L,'universal_emitted_block_dimension':u['dimension'],'all_nine_emitted_families_retained_by_exact_scalar_shifts':True,'first_light_probability_floating':float(abs(construction(2)['C'][0,1])**2*np.linalg.norm(v)**2),'certified_universal_first_amplitude_norm_lower':str(lower),'full_source_cubic_amplitude_error_upper':str(cubic),'numeric_certificate':num,'coefficient_surrogate_error_upper':str(u['operator_error']),'injection_surrogate_error_upper':str(u['injection_error']),'finite_volume_bound':vol,'infinite_line_light_channel_probability_lower':str(prob),'infinite_line_light_channel_probability_lower_decimal':float(prob),'original_input_probability_lower':str(F(6,25)*prob),'original_input_probability_lower_decimal':float(F(6,25)*prob),'preparation_success_lower':'6/25','complete_coarse_instrument':['preparation_failure_quantum_branch','prepared_light_family_pair_yes','prepared_all_other_pair_channels','prepared_no_pair'],'all_orthogonal_source_complements_retained':True,'conditional_pair_readout_absolute_work_upper':str(2*lmb*lmb*t)}
  rows.append(row)
 return {'schema':'oph.native_meson_controlled_memory_window.v1','verdict':'PASS_NONZERO_FULL_SOURCE_INFINITE_LINE_EMISSION_WINDOW','source':'same native bound-meson current, filled neutral reference and original fixed resonance-selected scalar box as openwave_meson_ring; no kinetic or interaction retuning','coupling':'1/64','exact_memory_injection_Gram':'C^dagger C=24 I6','exact_emitted_family_factorization':'D_fell=D00+(m_f+m_ell-2)I; injected column C_f1 psi with ||psi||^2=8','exact_first_emission_remainder':'||b(t)-b1(t)|| <= (sqrt24 G t)^3/6','family_amplitude_bound':family,'complete_native_small_volume_controls':native_controls,'full_unfiltered_small_volume_time_controls':[finite_control(L,F(1,10))for L in(2,4)],'controlled_windows':rows,'full_source_mean_energy_conserved':True,'preparation_nonselective_mean_work':'0 for the same complete native seed/filter instrument as the ring source','readout_work':'Pair/no-pair dephasing removes only the source off-diagonal contact, W=-2G Re<a,C^dagger b>; |W|<=2 (sqrt24 G)^2 t. Every emitted family projector is a subprojector; its complete work is the exact full-source Lüders expectation.','physical_strange_decay_or_irreversible_rate_claimed':False,'failed_small_box_peak_retained':'The time720 L2 peak falls under volume refinement; it is not used as a width. Finite rings recur. The theorem supplies a finite interval and volume control, not irreversible exponential decay.','scope':'Complete same-source finite-time emission, with binary64 propagation error, local exact coefficient error, cubic memory feedback and ring-to-line error separated. The infinite-line transfer concerns complete spatial/local instruments with their outside complement.'}
