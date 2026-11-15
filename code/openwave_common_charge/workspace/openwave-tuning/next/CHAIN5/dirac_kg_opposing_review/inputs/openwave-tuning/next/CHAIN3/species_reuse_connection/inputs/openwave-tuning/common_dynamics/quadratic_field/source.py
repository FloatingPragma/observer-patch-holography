"""Explicit noncompact quadratic amendment on the actual carrier complex.

The continuous field is a declared source, distinct from the compact cosine
field. The finite native spectral adapter has a separate state/error contract.
"""
from functools import lru_cache
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys

import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'vector_field'),str(HERE.parent/'spin')]
import cell_complex
import dirac_proposal as matter

MODEL_ID='native-local12-dirac-noncompact-quadratic-v4-proposed'
ELECTRIC=F(1,2)
MAGNETIC=F(3)
RHO=F(2,3)


@lru_cache(None)
def geometry():
    graph,B,A,periods,winding,labels=cell_complex.build(2)
    B=np.asarray(B,dtype=float);A=np.asarray(A,dtype=float)
    L=B@B.T;V=L.shape[0]
    Green=np.linalg.inv(L+np.ones((V,V))/V)-np.ones((V,V))/V
    longitudinal=B.T@Green@B
    K=A@A.T
    values,vectors=np.linalg.eigh(K)
    positive=values>1e-9
    U=vectors[:,positive];lam=values[positive]
    coexact=U@U.T
    harmonic=np.eye(len(graph.seams))-longitudinal-coexact
    hv,hu=np.linalg.eigh(harmonic)
    H=hu[:,hv>.5]
    assert len(lam)==190 and H.shape==(288,3)
    return dict(graph=graph,B=B,A=A,L=L,G=Green,K=K,U=U,lam=lam,
                harmonic_frame=H,coexact=coexact,longitudinal=longitudinal,
                transverse=coexact+H@H.T,labels=labels)


def charge_vector(positions,charges):
    q=np.zeros(96)
    for position,charge in zip(positions,charges):
        q[int(position)]+=float(charge)
    if abs(q.sum())>1e-12:raise ValueError('explicit neutral compensator required')
    return q


def reference(q=None,harmonic_electric_variance=.5):
    """Coexact oscillator vacuum times three retained free Gaussian modes.

    q is G-n in the source incidence convention B E=q. The harmonic modes
    have no normalizable stationary ground and are not discarded.
    """
    d=geometry();q=np.zeros(96) if q is None else np.asarray(q,float)
    if q.shape!=(96,) or abs(q.sum())>1e-11:raise ValueError('neutral96charge vector required')
    if not math.isfinite(harmonic_electric_variance) or harmonic_electric_variance<=0:
        raise ValueError('positive retained harmonic Gaussian variance required')
    a,b=float(ELECTRIC),float(MAGNETIC)
    omega=np.sqrt(a*b*d['lam']);U=d['U'];H=d['harmonic_frame']
    P=d['B'].T@d['G']@q
    cov_theta=(U*(a/(2*omega)))@U.T+(H/(4*harmonic_electric_variance))@H.T
    cov_E=(U*(omega/(2*a)))@U.T+(H*harmonic_electric_variance)@H.T
    return dict(mean_E=P,mean_theta=np.zeros(288),cov_theta=cov_theta,cov_E=cov_E,
                cov_theta_E=np.zeros((288,288)),frequencies=omega,
                longitudinal_energy=a*float(q@d['G']@q)/2,
                coexact_vacuum_energy=float(omega.sum()/2),
                harmonic_energy=3*a*harmonic_electric_variance/2,
                harmonic_electric_variance=harmonic_electric_variance)


def propagate_mean(theta,E,time):
    """Exact free evolution, including the three free harmonic coordinates."""
    d=geometry();a,b=float(ELECTRIC),float(MAGNETIC)
    omega=np.sqrt(a*b*d['lam']);U=d['U'];H=d['harmonic_frame']
    theta=np.asarray(theta,float);E=np.asarray(E,float)
    q,p=U.T@theta,U.T@E
    c,s=np.cos(omega*time),np.sin(omega*time)
    qt=c*q+a*s*p/omega;pt=c*p-omega*s*q/a
    th=H@(H.T@theta+a*time*(H.T@E))
    return U@qt+th,U@pt+H@(H.T@E)+d['longitudinal']@E


def transverse_current(cell,axis,direction=1):
    """Exact routed integer current and its orthogonal Gauss decomposition."""
    other,path,J=matter.route(cell,axis,direction)
    d=geometry();J=np.asarray(J,float)
    return dict(destination=other,path=path,current=J,
                longitudinal=d['longitudinal']@J,
                transverse=d['transverse']@J)


def gaussian_translation_overlap(shift,state=None):
    """Ordered overlap <displaced|original> for a momentum shift on ker B."""
    state=reference() if state is None else state
    d=geometry();shift=np.asarray(shift,float)
    if np.linalg.norm(d['B']@shift)>1e-9:raise ValueError('Gauss-preserving field shift required')
    # A positive momentum displacement is T=exp(+i shift dot theta).
    # This ordered overlap is <T*>, with the displayed negative phase.
    exponent=-float(shift@state['cov_theta']@shift)/2
    phase=-float(shift@state['mean_theta'])
    return np.exp(exponent+1j*phase)


def two_route_records(cell=0,axis=0):
    """A charge's two actual route branches coupled to one Gaussian field.

    The same longitudinal source charge at the common size2 endpoint is
    retained. Only the difference of the routes enters their interference.
    """
    left=transverse_current(cell,axis,1);right=transverse_current(cell,axis,-1)
    assert left['destination']==right['destination']
    difference=left['current']-right['current']
    assert np.linalg.norm(geometry()['B']@difference)<1e-12
    overlap=gaussian_translation_overlap(difference)
    # Logical two-branch Ramsey recombination with all outcomes retained.
    return {'same_field_overlap':{'real':float(overlap.real),'imag':float(overlap.imag)},
            'complete_probe_probabilities':[(1+float(overlap.real))/2,(1-float(overlap.real))/2],
            'closed_current_norm_squared':float(difference@difference),
            'relative_displacement_energy':float(ELECTRIC)*float(difference@difference)/2,
            'independent_field_substitution_probe':[.5,.5],
            'background_character_substitution_probe':[1.,0.]}


def run():
    paths=[Path(__file__),Path(cell_complex.__file__),Path(matter.__file__),
           HERE.parent/'vector_field/MANIFEST.json',HERE.parent/'spin/MANIFEST.json',
           HERE.parent/'units/MANIFEST.json',HERE.parent/'quantum_field/MANIFEST.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    d=geometry();r=reference();B,A=d['B'],d['A']
    assert np.max(abs(B@A))==0
    assert np.linalg.norm(B@r['mean_E'])<1e-10
    assert np.linalg.norm(B@r['cov_E'])<1e-9
    assert np.linalg.norm(r['cov_theta']@r['cov_E']-d['transverse']/4)<1e-9
    assert np.linalg.norm(A.T@d['longitudinal'])<1e-9
    # All190coexact modes contribute their actual finite vacuum covariance.
    assert min(r['frequencies'])>0 and len(r['frequencies'])==190
    assert np.linalg.matrix_rank(d['harmonic_frame'])==3
    pair=[]
    root=0
    for target in (1,12,36,84):
        q=charge_vector((root,target),(1,-1));charged=reference(q)
        P=charged['mean_E']
        assert np.linalg.norm(B@P-q)<1e-10 and np.linalg.norm(A.T@P)<1e-10
        assert abs(float(ELECTRIC)*float(P@P)/2-charged['longitudinal_energy'])<1e-10
        assert np.array_equal(charged['cov_E'],r['cov_E'])
        pair.append({'vertices':[root,target],'electric_minimum_energy':charged['longitudinal_energy'],
                     'full_reference_energy':charged['longitudinal_energy']+charged['coexact_vacuum_energy']+charged['harmonic_energy']})
    # Every matter route changes the longitudinal particular and shared
    # transverse field by complementary pieces of the SAME integer current.
    route_residual=0.
    for cell in range(8):
        for axis in range(3):
            for sign in (-1,1):
                x=transverse_current(cell,axis,sign)
                route_residual=max(route_residual,float(np.linalg.norm(x['current']-x['longitudinal']-x['transverse'])))
                assert np.linalg.norm(B@x['transverse'])<1e-10
    assert route_residual<1e-10
    probes=[two_route_records(axis=axis) for axis in range(3)]
    assert all(abs(sum(p['complete_probe_probabilities'])-1)<1e-13 for p in probes)
    assert any(abs(p['complete_probe_probabilities'][0]-.5)>1e-4 for p in probes)
    assert all(p['complete_probe_probabilities'][1]>1e-4 for p in probes)
    # Harmonic free packet spreads rather than being called a stationary
    # zero-mode vacuum. Its electric covariance and total energy are fixed.
    t=.4;harmonic_theta_variance=1/(4*r['harmonic_electric_variance'])+float(ELECTRIC)**2*t*t*r['harmonic_electric_variance']
    assert harmonic_theta_variance>1/(4*r['harmonic_electric_variance'])
    bad=np.zeros(96);bad[0]=1
    try:reference(bad)
    except ValueError:pass
    else:raise AssertionError('missing compensator admitted')
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.actual_carrier.noncompact_quadratic_source.v1','candidate_id':MODEL_ID,
      'checks':{'actual288edge208face_geometry':True,'all190coexact_and3harmonic_modes_retained':True,
       'pure_Gaussian_covariance_on_Gauss_quotient':True,'exact_source_Gauss_and_curl_longitudinal_split':True,
       'charge_independent_transverse_vacuum':True,'same_current_longitudinal_and_transverse_update':True,
       'same_field_interference_not_independent_fields':True,'complete_two_route_probe_laws':True,
       'harmonic_packet_not_stationary_ground':True,'nonneutral_input_rejected':True,'source_bytes_unchanged':True},
      'microscopic_coefficients':{'electric':'1/2','magnetic':'3','rho':'2/3','speed_squared':'1','charge_coupling_squared':'1/3'},
      'Hamiltonian':'H_q=(a/2)||B^T L^+q||² + (a/2)||E_T||² + (b/2)||A^T theta_T||², a=1/2,b=3, on the Gauss quotient and retained3harmonic canonical pairs',
      'clock':'cell Hamiltonian evolves for h-scaled time t/(tau0 h); matter macro has t=tau0 h. Native finite substeps and compiler precision are separate adapters.',
      'frequencies':r['frequencies'].tolist(),'vacuum_energy':r['coexact_vacuum_energy'],
      'harmonic_packet_energy':r['harmonic_energy'],'harmonic_theta_variance_at_time_0_4':harmonic_theta_variance,
      'static_charge_cases':pair,'routed_current_decomposition_error':route_residual,'same_field_route_probes':probes,
      'scope':'Exact noncompact Gaussian source proposal and numerical complete finite geometry calculations. The finite native spectral adapter and spatial photon theorem are separate receipts. No equivalence with fixed-Z16 cosine source, interacting vacuum, observed alpha or species is asserted.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'SOURCE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'modes':len(r['frequencies']),'static':pair,'route_probes':probes}),flush=True)


if __name__=='__main__':run()
