"""Quantum covariance of every retained coexact field mode on the L2 graph."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import quantum_state as qs

HERE=Path(__file__).resolve().parent


def run():
    paths=[Path(__file__),Path(qs.__file__),Path(qs.af.cell_complex.__file__),
           HERE.parent/'vector_field/LINEAR_RESPONSE.json']
    pins={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    graph,B,A,labels=qs.af.data();K=A@A.T
    lam,Q=np.linalg.eigh(K.astype(float))
    keep=lam>1e-10;lam=lam[keep];Q=Q[:,keep]
    assert len(lam)==190 and np.max(abs(B@Q))<2e-12
    rows=[]
    for eps2 in (F(1,3),F(1,12),F(1,48)):
        par=qs.parameters(eps2,8);alpha=float(par['alpha']);beta=float(par['beta'])
        kappa=np.sqrt(1-alpha*beta*lam/4)
        assert min(kappa)>0
        theta=(Q*(.5*np.sqrt(alpha/(beta*lam))/kappa))@Q.T
        electric=(Q*(.5*np.sqrt(beta*lam/alpha)/kappa))@Q.T
        cross=(Q*(.25*np.sqrt(alpha*beta*lam)/kappa))@Q.T
        V=np.block([[theta,cross],[cross,electric]])
        M=np.block([[np.eye(288)-alpha*beta*K,alpha*np.eye(288)],[-beta*K,np.eye(288)]])
        invariant=float(np.max(abs(M@V@M.T-V)))
        assert invariant<1e-11
        # This is a Gaussian on the190-dimensional coexact physical phase
        # space; the95 gauge and3 harmonic coordinates are separate.
        theta_modes=Q.T@theta@Q;electric_modes=Q.T@electric@Q;cross_modes=Q.T@cross@Q
        uncertainty=float(np.max(abs(theta_modes@electric_modes-cross_modes@cross_modes-np.eye(190)/4)))
        assert uncertainty<2e-12
        face_variance=np.sum(A*(theta@A),axis=0)
        assert min(face_variance)>0
        nonlinear_expectation=beta*sum(np.exp(-face_variance/2)-1+face_variance/2)
        quadratic_norm_budget=beta*math.sqrt(105)*float(sum(face_variance**2))/24
        rows.append({'epsilon_squared':str(eps2),'alpha':str(par['alpha']),'beta':str(par['beta']),
                     'coexact_modes':190,'all208_face_variances':face_variance.tolist(),
                     'minimum_face_variance':float(min(face_variance)),
                     'maximum_face_variance':float(max(face_variance)),
                     'cosine_minus_quadratic_kick_generator_mean':float(nonlinear_expectation),
                     'one_step_Gaussian_Taylor_state_norm_upper':quadratic_norm_budget,
                     'unscaled_electric_energy':float(np.trace(electric)/2),
                     'scaled_electric_kick_mean':float(alpha*np.trace(electric)/2),
                     'covariance_invariance_residual':invariant,'pure_uncertainty_residual':uncertainty})
    assert abs(rows[1]['one_step_Gaussian_Taylor_state_norm_upper']/rows[0]['one_step_Gaussian_Taylor_state_norm_upper']-.25)<1e-12
    assert pins=={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.local12.quantum_field_uv_state.v1','checks':{
        'actual288edge_curvature':True,'all190coexact_modes_retained':True,
        'Gauss_transverse_support':True,'full_Floquet_Gaussian_covariance_invariant':True,
        'pure_state_uncertainty':True,'all208plaquette_angle_variances':True,
        'weak_coupling_and_fixed_coupling_domains_separate':True,'source_bytes_unchanged':True},
        'harmonic_zero_mode_scope':'95gauge coordinates and3global harmonic coordinates are not normalized oscillators. The calculated Gaussian is on the190-dimensional coexact phase space. No infinite uniform zero-mode vector is used.',
        'cases':rows,
        'uniform_refinement_obligation':'At fixed epsilon^2 the local ultraviolet plaquette variance is not set to zero by an assumed low-momentum excitation. A crude global Taylor norm sum scales with the number of active faces and observation steps. This bound therefore requires quantified vacuum/renormalization control or a joint weak-coupling scaling; it is not a theorem forbidding a fixed-coupling Maxwell limit.',
        'quadratic_successor_contract':{'state_space':'Real coexact canonical field plus explicitly treated harmonic/charge sectors on the same B,A complex',
            'magnetic_generator':'one local quadratic function of each actual A_f^T theta, with its full branch/domain declared',
            'Gauss':'Electric currents translate the same field; BA=0 and routed endpoint continuity retained',
            'matter_coupling':'same Dirac route Wilson phase/current translation, with a proved native register approximation',
            'required_outputs':'native preparation, complete readout, clock, controlled register/spatial refinement and nonzero transverse spectral weight',
            'changes':'A noncompact/quadratic generator is an explicit architecture amendment. Compact cosine receipts do not transfer automatically.'},
        'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'UV_STATE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'cases':[{'eps2':r['epsilon_squared'],
      'variance_range':[r['minimum_face_variance'],r['maximum_face_variance']],
      'nonlinear_mean':r['cosine_minus_quadratic_kick_generator_mean'],
      'norm_budget':r['one_step_Gaussian_Taylor_state_norm_upper']} for r in rows]}),flush=True)


if __name__=='__main__':run()
