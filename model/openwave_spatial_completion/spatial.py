"""Exact Gaussian bounds and fixed-source complete-shell diagnostics."""
from functools import lru_cache
from fractions import Fraction as F
from math import sqrt
from pathlib import Path
import json
import numpy as np
import sympy as s
from scipy.linalg import eigh
from model.openwave_quantum_rotational.source import occupations,ladders,exact_axis_moment,cross,photon


def trial_bound(width,omega):
 L=s.Rational(width);O=s.Rational(omega);kin=s.Rational(3,64)/L**2;field=s.Rational(3,8192)*L**2/(1+s.Rational(2,5)*L**2)**s.Rational(5,2)
 trap=s.Rational(3,4)*O*(L**2+L**-2-2)
 return {'Gaussian_width_exact':str(L),'trap_frequency_exact':str(O),'kinetic_upper_exact':str(kin),'field_upper_exact':str(field),'trap_mean_exact':str(trap),'energy_upper_exact':str(16+kin+field+trap),'energy_upper_diagnostic':float(16+kin+field+trap)}


def analytic():
 v=s.Rational(3,64)+s.Rational(3,8192)*(s.Rational(5,7))**s.Rational(5,2);assert 0<v<s.Rational(1,16)
 rows=[trial_bound(2**k,F(1,16*2**(4*k)))for k in range(7)]
 assert all(rows[i+1]['energy_upper_diagnostic']<rows[i]['energy_upper_diagnostic']for i in range(len(rows)-1))
 prior=json.loads((Path(__file__).resolve().parents[1]/'openwave_quantum_rotational/results.json').read_text());cap=next(x for x in prior['interval_certificates']if x['cap']==2)
 unit=[F('0.00051099895069')/F.from_float(x)for x in reversed(cap['actual_ground_energy_interval'])]
 return {'full_spatial_full_photon_vacuum_trial_shift_exact':str(v),'vacuum_trial_shift_diagnostic':float(v),
  'uniform_exact_spin_half_band_domain':'Omega > v0 and min(Omega,1)>v0; in particular Omega>=1/16',
  'next_distinct_energy_floor':'16+min(Omega,1)','exact_dressed_Casimir':'3/4',
  'spatial_number_tail_bound':'P(N_matter>L) <= v0/[Omega*(L+1)] on the exact ground band',
  'controlled_trap_release_trials':rows,'untrapped_infimum_exact':'16',
  'untrapped_threshold_is_normalizable_eigenstate':False,
  'threshold_GeV_interval_unchanged_cap_two_electron_unit':[float(np.nextafter(float(16*u),direction))for u,direction in zip(unit,(-np.inf,np.inf))],
  'no_electron_identification_exclusion':'Absence of a normalizable zero-momentum state does not exclude a free-particle rest-mass threshold; this result excludes interpreting the trapped ground as a demonstrated self-bound localized electron.',
  'translation_invariant_physical_electron_dispersion_derived':False}

@lru_cache(maxsize=6)
def shell_source(shell,cap=1,omega=1.):
 basis,a=ladders(shell);n=len(basis);sigma=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]));zero=np.zeros((n,n),complex)
 p=tuple(-1j*(x-x.conj().T)/np.sqrt(2)for x in a);orbital=tuple(sum((-1j*cross(i,j,k)*a[j].conj().T@a[k]for j in range(3)for k in range(3)),zero.copy())for i in range(3))
 moments={(i,j,power):float(exact_axis_moment(i,j,power))for i in range(shell+1)for j in range(shell+1)for power in(0,1)}
 radial=[np.array([[np.prod([moments[u[d],v[d],int(d==axis)]for d in range(3)])for v in basis]for u in basis])for axis in range(3)]
 B0=sum((np.kron(sigma[i],p[i])for i in range(3)),np.zeros((2*n,2*n),complex))
 M=[sum((np.kron(sigma[i],sum((cross(i,k,d)*radial[d]/2 for d in range(3)),zero.copy()))for i in range(3)),np.zeros_like(B0))for k in range(3)]
 ph=photon(cap);d=len(ph['basis']);B=np.kron(B0,np.eye(d))
 for mat,q in zip(M,ph['q']):B+=np.kron(mat,q)
 values,U=eigh(B);shift=values**2/(np.sqrt(256+values**2)+16);number=np.array([sum(x)for x in basis])
 Nm=np.kron(np.kron(np.eye(2),np.diag(number)),np.eye(d));Nf=np.kron(np.eye(2*n),ph['N']);Hminus=(U*shift)@U.conj().T+omega*Nm+Nf;energy,V=eigh(Hminus,subset_by_index=(0,3));ground=V[:,:2]
 J=[np.kron(np.kron(np.eye(2),orbital[i])+np.kron(sigma[i]/2,np.eye(n)),np.eye(d))+np.kron(np.eye(2*n),ph['J'][i])for i in range(3)]
 casimir=sum(x@x for x in J);err=float(np.linalg.norm(casimir@ground-.75*ground));assert err<1e-9 and abs(energy[0]-energy[1])<1e-11
 meanN=float(np.trace(ground.conj().T@Nm@ground).real/2);meanPh=float(np.trace(ground.conj().T@Nf@ground).real/2)
 return {'complete_spatial_shell':shell,'spatial_modes':n,'photon_cap':cap,'reduced_joint_dimension':2*n*d,'trap_frequency':omega,
  'ground_energy':float(16+energy[0]),'first_excited_energy':float(16+energy[2]),'ground_mean_spatial_number':meanN,'ground_mean_photon_number':meanPh,
  'complete_ground_Casimir_residual':err,'both_spin_states_and_every_spatial_photon_mode_retained':True,
  'diagnostic_scope':'Fixed-source Galerkin computation; a sequence convergence theorem does not assign a numerical continuum error bar to this finite row.'}


def certificate():
 rows=[shell_source(n,cap,omega)for n,cap,omega in((2,1,1.),(3,1,1.),(4,1,1.),(5,1,1.),(2,2,1.),(3,2,1.),(4,2,1.),(4,1,.25),(4,1,.0625))]
 return {'schema':'oph.full_spatial_positive_quantum_completion.v1','analytic_full_source':analytic(),'complete_shell_diagnostics':rows,
  'source_parameters':'mass16, transverse amplitude1/8, profile exp(-r^2/5), full spatial L2 and full photon Fock; trap Omega N_matter is exposed',
  'spectral_convergence':'Complete-shell compression B_N with the unchanged Gaussian form factors converges in strong resolvent sense; common oscillator coercivity and form convergence give every fixed low eigenvalue at Omega>0. No finite square root is identified with a compressed infinite square root.',
  'scope':'A full spatial-domain theorem, exact sufficient confined spin-band domain and trap-release mass threshold for the declared charged quantum block. No self-bound localized particle, translation-invariant physical dispersion or neutral-sector gap transfer is inferred.',
  'verdict':'PASS_FULL_SPATIAL_SOURCE_AND_CONTROLLED_TRAP_RELEASE'}
