"""Exact reduced spectra, common-source binding and complete finite records."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from mpmath import iv
from model.openwave_baryon_ordering.source import balanced_baryon,spectral_certificate as baryon_certificate,frobenius_ceil
from model.openwave_quantum_rotational.enclosure import scalar
from model.common_source.chiral_candidate.action import SOURCE_FILES
from .source import reduction,MODES,native_Fierz,exact_product


def spectral_enclosure(exact):
 n=exact.rows;h=np.asarray(exact.evalf(18),float);D=2**40;iv.dps=55;scalar.cache_clear();H=np.rint(h*D).astype(np.int64);entry=F(0)
 for i in range(n):
  for j in range(n):
   if exact[i,j]==0:continue
   x=scalar(exact[i,j]);lo=F.from_float(float(np.nextafter(float(x.a),-np.inf)));hi=F.from_float(float(np.nextafter(float(x.b),np.inf)));mid=F(int(H[i,j]),D);entry=max(entry,abs(lo-mid),abs(hi-mid))
 delta=n*entry;e,Q=eigh(h);O=np.array([[int(round(float(x)*D))for x in row]for row in Q],dtype=object);E=np.array([int(round(float(x)*D))for x in e],dtype=object)
 eta=F(frobenius_ceil(O.T@O-D*D*np.eye(n,dtype=object)),D*D);assert eta<F(1,2)
 residual=np.asarray(H,dtype=object)@O-O*E;rho=F(frobenius_ceil(residual),D*D)+delta*(1+eta);K=(F(max(sum(abs(int(x))for x in row)for row in H),D)+delta).__ceil__()+1;error=2*rho+4*K*eta
 assert max(abs(int(x))for x in E)<K*D
 return {'complete_reduced_dimension':n,'dyadic_denominator':D,'source_operator_rounding_upper':str(delta),'exact_Gram_error_upper':str(eta),'exact_residual_upper':str(rho),'all_eigenvalue_error_upper':str(error),
  'all_eigenvalue_intervals':[[str(F(int(x),D)-error),str(F(int(x),D)+error)]for x in E],
  'dyadic_eigenframe':[[str(int(x))for x in row]for row in O],'dyadic_eigenvalues':[str(int(x))for x in E]}

@lru_cache(None)
def old_baryon_enclosure():
 d=balanced_baryon();B=d['Binteger'];A=exact_product(B.T,exact_product(d['Hnum']-d['weak_term_numerator'],B)).toarray();g=d['Gram'];n=len(g)
 exact=sy.Matrix(n,n,lambda i,j:sy.Integer(int(A[i,j]))/(48*sy.sqrt(int(g[i])*int(g[j]))))
 return spectral_enclosure(exact)

def exact_matvec(A,v):
 return [sum(int(A.data[k])*v[int(A.indices[k])]for k in range(A.indptr[i],A.indptr[i+1]))for i in range(A.shape[0])]

def sqrt_upper_fraction(value,bits=80):
 numerator=value.numerator*2**(2*bits);denominator=value.denominator;n=isqrt(numerator//denominator)
 if n*n*denominator<numerator:n+=1
 return F(n,2**bits)


def trial_records(data,spectrum):
 e,U=eigh(data['h']);coefficient=data['coefficient_transform']@U[:,0];integers=[int(round(float(x)*2**56))for x in coefficient];w=exact_matvec(data['Binteger'],integers);D=sum(x*x for x in w);Hw=exact_matvec(data['Hnum'],w);A=sum(x*y for x,y in zip(w,Hw));E=F(A,48*D)
 residual_squared=F(sum((D*y-A*x)**2 for x,y in zip(w,Hw)),48**2*D**3);residual=sqrt_upper_fraction(residual_squared);next_floor=F(spectrum['all_eigenvalue_intervals'][1][0]);gap=next_floor-E;assert gap>0;trace=2*residual/gap
 full=np.zeros(data['full_dimension']);prob=[];qnumerator=0;colocation=0
 for i,index in enumerate(data['keep']):
  mask=data['basis'][int(index)];p=F(w[i]*w[i],D);full[int(index)]=float(p);n0=sum(1 for j,v in enumerate(MODES)if v[0]==0 and(mask>>j)&1);n1=6-n0
  qnumerator+=w[i]*w[i]*n0*n1
  if n0 in(0,6):colocation+=w[i]*w[i]
 quadrupole=F(qnumerator,18*D);delta=F(29,4)-E;assert delta>0
 # H_Q>=6, ||P H Q||<=6t=3/4; p<=sqrt(p) gives this rational bound.
 outside_lower=(delta/F(11,4))**2;Qlower=F(5,18)*outside_lower
 return {'trial_basis_integer_coefficients':list(map(str,integers)),'full_Cartan_sector_integer_coefficients':list(map(str,w)),'exact_trial_squared_norm':str(D),
  'exact_Rayleigh_energy':str(E),'Rayleigh_energy_diagnostic':float(E),'exact_residual_squared':str(residual_squared),'ground_line_trace_distance_upper':str(trace),
  'exact_trial_intrinsic_axial_quadrupole':str(quadrupole),'intrinsic_axial_quadrupole_diagnostic':float(quadrupole),
  'ground_quadrupole_interval':[str(max(F(0),quadrupole-trace/2)),str(min(F(1,2),quadrupole+trace/2))],
  'all_ground_vectors_quadrupole_lower_from_binding':str(Qlower),'number_centre_operator':'X_N=(sum_i z_i N_i)/6; z_i=0,1',
  'I_zero_exact_charge_quadrupole_identity':'Qzz=2 sum q_i(z_i-X_N)^2 has expectation <N0*N1>/18',
  'variational_trial_co_location_probability':float(F(colocation,D)),
  'full_sector_occupation_masks':[str(x)for x in data['basis']],
  'complete_initial_uniform_sector_dimension':data['full_dimension'],
  'exact_preparation_success_probability':str(F(1,data['full_dimension'])),
  'complete_success_occupation_probabilities':list(map(float,full/data['full_dimension'])),
  'complete_failure_occupation_probabilities':list(map(float,(1-full)/data['full_dimension'])),
  'full_fixed_sector_complement_outcome_probability':0.,
  'complete_records_total_variation_to_exact_ground_filter_upper':str(trace/data['full_dimension']),
  'trial_complete_records_are_exact_stationary_ground_records':False,
  'conditional_trial_vs_stationary_ground_trace_bound_uniform_all_times':str(trace),
  'maximally_mixed_initial_energy':str(F(data['full_H_trace_numerator'],48*data['full_dimension'])),
  'nonselective_complete_filter_mean_work_exact':'0 on the maximally mixed admitted input',
  'retained_idle_reference_charge':'-1',
  'quadrupole_unit':'squared source-cell separation; no physical length calibration'}


def certificate():
 d=reduction();balanced=spectral_enclosure(d['H_exact']);unbalanced=spectral_enclosure(d['H0_exact']);baryon=baryon_certificate();old=old_baryon_enclosure();trial=trial_records(d,balanced)
 lowerB=F(baryon['lowest_interval_rational'][0]);upperOldB=F(old['all_eigenvalue_intervals'][0][1]);upperN=F(trial['exact_Rayleigh_energy']);gap=2*lowerB-upperN;assert gap>F(6,25)
 unbound=F(unbalanced['all_eigenvalue_intervals'][0][0])-2*upperOldB;assert unbound>F(9,10000)
 # Every partition except3+3 contains compensating nonzero triality fragments.
 other=F(6)+F(8,3);assert other>upperN
 return {'schema':'oph.complete_native_six_quark_binding.v1','source':'full60-CAR two-cell action mN+positive graph K+sum C_x^2+(1/4)sum S_x^2+(1/4)sum I_colourful,x^2; m1,t1/8',
  'exact_native_current_attachment':native_Fierz(),
  'architecture_selection':'Equal spin/weak-current costs impose a declared spin-weak exchange symmetry; the1/4 coefficient is inherited from the existing spin cost. This equality is an architectural input, not a derived physical coupling. The no-weak-cost member is retained as a failed S1,I0 binding control.',
  'same_native_weak_current_only_on_Q_doublet':True,'inserted_composite_projector_or_transition_action':False,
  'exact_dimensions':{'fixed_N6_Nu3_Nspinup4':d['full_dimension'],'zero_colour_Cartan':d['zero_Cartan_dimension'],'complete_colour0_S1_I0_Sz1_Iz0':27},
  'projector_denominator':d['Pden'],'projector_steps':d['projector_steps'],'independent_integer_basis_columns':d['basis_columns'],
  'exact_basis_Gram':[[str(x)for x in row]for row in d['Gram'].tolist()],'positive_LDL_diagonal':[str(x)for x in d['LDL_diagonal']],
  'complete_balanced_spectrum':balanced,'complete_unbalanced_control_spectrum':unbalanced,
  'same_balanced_baryon_spectrum_certificate':baryon,'same_unbalanced_baryon_spectrum_certificate':old,
  'strict_binding_below_two_baryons_lower':str(gap),'strict_binding_diagnostic':float(gap),
  'all_other_fragment_partition_threshold_lower':str(other),
  'unbalanced_S1_I0_ground_above_two_baryons_lower':str(unbound),'unbalanced_control_diagnostic':float(unbound),
  'positive_bound_state':trial,
  'conditional_common_electron_unit':{'electron_anchor_GeV':'0.00051099895069','baseline_equal_bare_mass_six_quark_energy_GeV':float(F('0.00051099895069')*upperN),'binding_gap_lower_GeV':str(F('0.00051099895069')*gap),'physical_deuteron_mass_or_binding_identified':False},
  'additional_degenerate_partner':'At equal spin and weak-current costs, swapping spin and weak labels maps S1,I0 to S0,I1. The isotriplet has charge0,1,2; it is retained, not identified with observed nuclei.',
  'unequal_bare_up_down_mass_shift_cancels_binding':'N6 has3u3d and a proton/neutron pair has3u3d, so their common number-energy shift cancels exactly. This does not match physical binding without the missing scale and dynamics.',
  'finite_perturbation_domain':'For perturbations preserving the admitted occupation sectors, if both one-source six-quark and each one-source baryon generator change by norm at most epsilon, binding survives3epsilon<the certified gap.',
  'scope':'An exact finite native-current six-quark composite below every nontrivial quark-fragment partition, with positive intrinsic graph charge quadrupole. Internal spin is not spatial Noether angular momentum; no physical deuteron identity, fm length or observed nuclear binding is inferred.',
  'verdict':'PASS_COMPLETE_NATIVE_SIX_QUARK_BINDING_AND_POSITIVE_GRAPH_QUADRUPOLE'}
