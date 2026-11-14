"""Prepared quantum states on an exact native field-cycle sector.

The cycle is an actual internal face of the frozen local complex. Its
oscillator is optical. No low-momentum photon or fixed-coupling continuum
claim is inferred from this reducing sector.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy.special import erfc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/'vector_field'))
sys.path.insert(0, str(ROOT/'openwave-tuning/auditor/field_cutoff'))
import autonomous as af
import native_field as nf
import cutoff

RHO = F(2,3)
FACE = 0
CYCLE_NORM2 = 3
A = float(F(CYCLE_NORM2)/RHO)
OMEGA = math.sqrt(A)


def parameters(epsilon_squared=F(1,3), substeps=8):
    eps2 = F(epsilon_squared)
    if eps2 <= 0 or type(substeps) is not int or substeps < 2:
        raise ValueError('positive coupling and at least two field substeps required')
    delta = F(1,substeps)
    return {'epsilon_squared': eps2, 'delta': delta,
            'alpha': delta*eps2/RHO, 'beta': delta/eps2,
            'speed_squared': RHO*(delta*eps2/RHO)*(delta/eps2)*substeps**2,
            'charge_coupling_squared': eps2}


def gaussian_data(epsilon_squared, substeps, displacement=1.):
    """Metaplectic Floquet vacuum and an explicitly displaced Gaussian.

    In Q=q/A**(1/4), P=A**(1/4)*p, the positive symplectic squeeze S
    satisfies S S^T = [[1,s/2],[s/2,1]]/kappa. Its displacement is S(d,0).
    The momentum wavefunction's global phase is fixed by positive value
    at its mean; the cosine/quadratic source phase is never dropped.
    """
    par = parameters(epsilon_squared, substeps)
    s = float(par['delta'])*OMEGA
    kappa = math.sqrt(1-s*s/4)
    u,v = math.sqrt((1+s/2)/kappa), math.sqrt((1-s/2)/kappa)
    S = np.array([[(u+v)/2,(u-v)/2],[(u-v)/2,(u+v)/2]])
    mu = S @ np.array([float(displacement),0.])
    covariance = np.array([[1,s/2],[s/2,1]])/(2*kappa)
    scale = A**.25
    return {**par,'epsilon':math.sqrt(float(epsilon_squared)), 'kappa':kappa,
            's':s,'omega':2*math.asin(s/2),'S':S,
            'mean_q':scale*mu[0],'mean_p':mu[1]/scale,
            'var_q':scale*scale*covariance[0,0],
            'var_p':covariance[1,1]/(scale*scale),
            'cov_qp':covariance[0,1], 'displacement':float(displacement)}


def centered(modulus):
    if type(modulus) is not int or modulus < 4 or modulus & (modulus-1):
        raise ValueError('power-of-two field alphabet required')
    return (np.arange(modulus)+modulus//2)%modulus-modulus//2


def gaussian_cycle(modulus=16, epsilon_squared=F(1,3), substeps=8,
                   displacement=1., cap=None, amplitude_bits=44):
    """Exact integer-amplitude state plus literal binary64 phase array.

    The Gaussian comparison is separate from this exact finite target.
    No data-dependent success normalization occurs in the experiment.
    """
    if type(amplitude_bits) is not int or not 12 <= amplitude_bits <= 48:
        raise ValueError('declared integer amplitude precision required')
    z=centered(modulus); data=gaussian_data(epsilon_squared,substeps,displacement)
    p=data['epsilon']*z
    offset=p-data['mean_p']
    gaussian=np.exp(-offset**2/(4*data['var_p']))
    values=[int(round(float(value)*(1<<amplitude_bits))) for value in gaussian]
    if cap is not None:
        if type(cap) is not int or not 0 <= cap < modulus//2:
            raise ValueError('strictly interior support cap required')
        values=[value if abs(int(e))<=cap else 0 for value,e in zip(values,z)]
    norm2=sum(v*v for v in values)
    real=np.asarray(values,dtype=float)/math.sqrt(norm2)
    phases=-data['cov_qp']*offset**2/(2*data['var_p'])-data['mean_q']*offset
    state=real*np.exp(1j*phases)
    return state,values,phases,data


def source_step(state, alpha, beta, quadratic=False):
    """Complete finite U_E U_B, retaining exp(i beta) in both variants."""
    modulus=len(state);z=centered(modulus)
    angles=2*np.pi*z/modulus
    angle_state=np.fft.fft(state,norm='ortho')
    potential=1-angles**2/2 if quadratic else np.cos(angles)
    angle_state*=np.exp(1j*float(beta)*potential)
    after=np.fft.ifft(angle_state,norm='ortho')
    return after*np.exp(-.5j*float(alpha)*CYCLE_NORM2*z*z)


def complete_records(state):
    field=abs(state)**2
    angle=abs(np.fft.fft(state,norm='ortho'))**2
    assert abs(float(field.sum())-1)<2e-11 and abs(float(angle.sum())-1)<2e-11
    return {'electric_probabilities':field.tolist(),'angle_probabilities':angle.tolist(),
            'unconditional_norm_squared':float(field.sum()),'renormalized':False}


def eight_moment(mean, variance):
    return mean**8+28*mean**6*variance+210*mean**4*variance**2+420*mean**2*variance**3+105*variance**4


def oscillator_bound(epsilon_squared, substeps, macros=1, displacement=1.):
    """Explicit analytic compact-rotor to periodized Gaussian bound.

    Values are numerical evaluations of formulas proved in PROOF.md;
    the closed symbolic formula, not these rounded values, is the bound.
    Native and finite-register errors are supplied separately.
    """
    data=gaussian_data(epsilon_squared,substeps,displacement)
    eps=data['epsilon'];L=2*math.pi/eps
    v=data['var_q'];eta=1/(4*v)
    M=np.array([[1-A/float(substeps**2),A/substeps],[-1/substeps,1.]])
    means=np.array([data['mean_q'],data['mean_p']])
    moments=[]
    steps=substeps*macros
    for _ in range(steps):
        mu=float(means[0]);D=1-2*eta*v
        moment=math.exp(eta*mu*mu/D)/math.sqrt(D)*eight_moment(mu/D,v/D)
        moments.append(moment);means=M@means
    c=data['var_p']*L*L/2
    tau=2*math.exp(-c)/(1-math.exp(-3*c))
    periodic_weight=1+2*math.exp(-eta*L*L/4)/(1-math.exp(-2*eta*L*L))
    normalization_lower=math.sqrt(max(0.,1-tau))
    nonlinear=float(data['delta'])*float(epsilon_squared)*math.sqrt(periodic_weight)*sum(math.sqrt(x) for x in moments)/24
    if normalization_lower:
        nonlinear/=normalization_lower
        normalization_drift=2*tau/((1+normalization_lower)*normalization_lower)
    else:
        nonlinear=normalization_drift=2.
    r=.5*math.atanh(data['s']/2)
    nbar=displacement**2/2
    phase_bound=steps*data['s']**3/(24*(1-data['s']**2/4)**1.5)
    splitting=r*math.sqrt(nbar*nbar+4*nbar+2)+phase_bound*math.sqrt(nbar*nbar+2*nbar+.25)
    return {'nonlinear_periodized_state_norm_bound':nonlinear,
            'periodization_normalization_drift_bound':normalization_drift,
            'compact_to_normalized_Floquet_Gaussian_bound':min(2.,nonlinear+normalization_drift),
            'squeeze_and_time_splitting_bound':splitting,
            'periodization_norm_squared_deviation':tau,
            'nonzero_reference_one_quantum_weight':nbar*math.exp(-nbar),
            'weak_coupling_parameter':float(epsilon_squared),
            'scalar_phase_per_substep':float(data['beta']),
            'oscillator_Floquet_phase':data['omega'],'oscillator_frequency':OMEGA,
            'cell_time':macros,'native_clock_phase':macros/8}


def cycle_checks():
    graph,B,faces,labels=af.data();face=faces[:,FACE]
    assert int(face@face)==CYCLE_NORM2 and not np.any(B@face)
    assert labels[FACE][0]=='internal_triangle'
    rows=[]
    for modulus in (16,32,64):
        for value in range(modulus):
            field=(value*face)%modulus
            e=(field+modulus//2)%modulus-modulus//2
            z=int(centered(modulus)[value])
            assert int(e@e)==CYCLE_NORM2*z*z
            assert not np.any((B@field)%modulus)
            shifted=(field+face)%modulus
            assert np.array_equal(shifted,((value+1)*face)%modulus)
        rows.append({'modulus':modulus,'all_cycle_columns':modulus})
    return face,rows


def run():
    paths=[Path(__file__),Path(af.__file__),Path(af.cell_complex.__file__),Path(nf.__file__),
           Path(cutoff.__file__),HERE.parent/'vector_field/MANIFEST.json',
           HERE.parent/'spin/MANIFEST.json',HERE.parent/'units/MANIFEST.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    face,columns=cycle_checks()
    experiments=[];all_records={};overlaps=[]
    for q,eps2 in ((16,F(1,3)),(128,F(1,12)),(512,F(1,48)),(2048,F(1,192))):
        cap=min(q//2-1,math.ceil(6/math.sqrt(float(eps2))))
        state,values,phases,data=gaussian_cycle(q,eps2,8,cap=cap)
        cos=state.copy();quad=state.copy();telescope=0.;records=[]
        angles=2*np.pi*centered(q)/q
        for n in range(9):
            records.append({'substeps':n,'cosine':complete_records(cos),'quadratic':complete_records(quad)})
            if n==8:break
            # This is the exact finite-state Taylor/remainder control;
            # its numerical evaluation is independently bounded by the
            # continuum Gaussian estimate, with explicit register tails.
            angle_probs=abs(np.fft.fft(quad,norm='ortho'))**2
            telescope+=float(data['beta'])*math.sqrt(float(angle_probs@(angles**8)))/24
            cos=source_step(cos,data['alpha'],data['beta'])
            quad=source_step(quad,data['alpha'],data['beta'],True)
        difference=float(np.linalg.norm(cos-quad))
        assert difference<=telescope+1e-11
        electric_tv=float(np.sum(abs(abs(cos)**2-abs(quad)**2))/2)
        angle_tv=float(np.sum(abs(abs(np.fft.fft(cos,norm='ortho'))**2-abs(np.fft.fft(quad,norm='ortho'))**2))/2)
        overlap=np.vdot(state,np.roll(state,1))
        shift_error=float(np.linalg.norm(np.roll(state,1)-state))
        assert abs(shift_error**2-(2-2*overlap.real))<1e-12
        overlaps.append({'modulus':q,'epsilon_squared':str(eps2),
                         'one_cycle_translation_overlap_real':float(overlap.real),
                         'one_cycle_translation_overlap_imag':float(overlap.imag),
                         'one_axis_code_intertwiner_norm_defect':shift_error,
                         'field_energy_expectation':CYCLE_NORM2*float(abs(state)**2@(centered(q)**2))/2,
                         'scaled_electric_kick_expectation':float(data['alpha'])*CYCLE_NORM2*float(abs(state)**2@(centered(q)**2))/2})
        register=cutoff.transfer_bound(q,8,data['beta'],cap,0)
        row={'modulus':q,'epsilon_squared':str(eps2),'substeps':8,
             'alpha_exact':str(data['alpha']),'beta_exact':str(data['beta']),
             'initial_support_cap':cap,'numeric_state_error_cosine_quadratic':difference,
             'finite_state_Taylor_telescope':telescope,'electric_TV':electric_tv,'angle_TV':angle_tv,
             'register_to_integer_rotor_norm_bound':float(register),
             'register_bound_exact':str(register), 'analytic_oscillator':oscillator_bound(eps2,8)}
        experiments.append(row);all_records[str(q)]=records
    assert all(experiments[i+1]['numeric_state_error_cosine_quadratic']<experiments[i]['numeric_state_error_cosine_quadratic'] for i in range(3))
    # Uniform phase has no small-angle guarantee; F0 is exactly uniform
    # under the angle instrument. Treating it as the Gaussian is rejected.
    theta=2*np.pi*centered(16)/16
    flat_angle=np.ones(16)/16
    narrow=complete_records(gaussian_cycle()[0])['angle_probabilities']
    flat_control=float(np.sum(abs(flat_angle-np.array(narrow)))/2)
    assert flat_control>.4
    # Removing exp(i beta) preserves separate field probabilities but
    # changes a controlled interferometer, so its phase is retained.
    beta=float(F(3,8));global_control=abs(1-np.exp(1j*beta))
    assert global_control>.3
    splitting=[{'substeps':n,**oscillator_bound(F(1,192),n)} for n in (8,16,32,64)]
    assert all(splitting[i+1]['squeeze_and_time_splitting_bound']<splitting[i]['squeeze_and_time_splitting_bound'] for i in range(3))
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'QUANTUM_CYCLE_RECORDS.json').write_text(json.dumps(all_records,separators=(',',':'))+'\n')
    result={'schema':'oph.local12.quantum_field_cycle.v1',
      'candidate':'pure-gauge reducing sector of native-local12-fermion-plaquette-v2-proposed',
      'checks':{'actual_face_cycle_and_Gauss':True,'all_modular_cycle_columns':True,
                'same_288edge_electric_and_plaquette_restriction':True,
                'complete_electric_and_angle_outcomes':True,'weak_angle_state_error_decreases':True,
                'scalar_cosine_phase_retained':True,'separate_substep_refinement':True,
                'finite_energy_gaussian_shift_defect_computed':True,
                'definite_F0_not_small_angle_gaussian':True,'source_bytes_unchanged':True},
      'face':FACE,'face_edges':[(int(e),int(face[e])) for e in np.flatnonzero(face)],
      'cycle_norm_squared':CYCLE_NORM2,'full_harmonic_rho':str(RHO),
      'common_parameter_family':'delta=1/N; alpha=delta*epsilon^2/rho; beta=delta/epsilon^2. Macro speed1 and charge-coupling-squared epsilon^2.',
      'cycle_oscillator':'p=epsilon*z,q=theta/epsilon; H=(A*p^2+q^2)/2,A=9/2; physical cell-scale frequency=(3/sqrt2)/(tau0*h)',
      'column_checks':columns,'experiments':experiments,'substep_controls':splitting,
      'shift_overlaps':overlaps,'F0_vs_Gaussian_angle_TV':flat_control,
      'dropped_scalar_phase_vector_defect':global_control,
      'records_sha256':hashlib.sha256((HERE/'QUANTUM_CYCLE_RECORDS.json').read_bytes()).hexdigest(),
      'scope':'Actual prepared quantum one-face loop with controlled weak-angle oscillator comparison. The oscillator is optical. The weak-angle sequence sends the charge coupling to zero. Full transverse photon preparation, homogeneous ultraviolet state control, interacting continuum, matter vacuum and calibrated forces are not established.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'QUANTUM_CYCLE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'experiments':[{k:v for k,v in row.items() if k!='register_bound_exact'} for row in experiments],
                      'overlaps':overlaps,'F0_angle_TV':flat_control}),flush=True)


if __name__=='__main__':run()
