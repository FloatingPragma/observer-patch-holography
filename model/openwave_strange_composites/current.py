"""Complete native strange-current instruments, including every failure record."""
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
import json
from pathlib import Path
from .source import source,select,ORBITALS,INDEX,lifted_generator,parent


def mass_parameters():
    d=json.loads((Path(__file__).parent/'inputs.json').read_text())
    dm,pair=F(d['strange_number_coefficient']),F(d['strange_pair_coefficient'])
    assert dm>0 and pair>0
    return dm,pair


def weak(sine=F('0.22431'),axial=False):
    d=source('baryon');s=float(sine);c=np.sqrt(1-s*s)
    C=np.array([[c,s],[-s,c]])
    assert np.max(abs(C.T@C-np.eye(2)))<2e-16
    a=np.zeros((16,16));z=np.zeros(16)
    for i,(x,f,spin,w) in enumerate(ORBITALS):
        z[i]=(1-2*w)/2
        if w==0:
            for ff in range(2):a[i,INDEX[x,ff,spin,1]]=C[f,ff]*(1-2*spin if axial else 1)
    return lifted_generator(d,a),lifted_generator(d,np.diag(z))


def instrument(sine=F('0.22431')):
    d=source('baryon');dm,pair=mass_parameters()
    L=select('baryon',(1,1,0,1),(0,0,0,0),.5,0,float(dm),float(pair))
    proton=select('baryon',(2,1,0,0),(0,0,0,0),.5,.5,float(dm),float(pair))
    J,Z=weak(sine);Jz,_=weak(sine,True);M=J.conj().T@J;eye=sp.eye(len(d['labels']),format='csr')
    algebra=J@J.conj().T-M-2*Z
    assert np.max(abs(algebra.data),initial=0)<3e-14
    polynomial=M@(M-eye)@(M-3*eye)@(M-4*eye)
    assert np.max(abs(polynomial.data),initial=0)<3e-12
    failure=sp.csr_matrix(eye.shape,dtype=complex)
    for x in (0,1,3,4):
        P=eye.copy()
        for y in (0,1,3,4):
            if y!=x:P=P@(M-y*eye)/(x-y)
        failure+=np.sqrt(1-x/9)*P
    assert np.max(abs((failure.conj().T@failure+M/9-eye).data),initial=0)<3e-13
    psi=L['state'];success=J@psi/3;failed=failure@psi
    assert abs(np.vdot(success,success)+np.vdot(failed,failed)-1)<2e-13
    H=d['H0']+sp.diags(float(dm)*d['Nstrange']+float(pair)*d['strange_pair']+10000*d['Ncharm'])
    current_energy=float(np.vdot(success,H@success).real+np.vdot(failed,H@failed).real)
    measured_energy=float(np.dot(H.diagonal().real,abs(success)**2+abs(failed)**2))
    records={}
    for q in sorted(set(map(tuple,d['flavour_q']))):
        ids=np.flatnonzero(np.all(d['flavour_q']==q,axis=1))
        records[','.join(map(str,q))]={'success_probability':float(np.sum(abs(success[ids])**2)),
                                      'failure_probability':float(np.sum(abs(failed[ids])**2))}
    raw_success=d['B']@success;raw_failure=d['B']@failed
    scale=1/len(d['labels'])
    return {'Cabibbo_sine_input':str(sine),'mass_coefficients':[str(dm),str(pair)],
            'Lambda_to_proton_vector_amplitude':complex(np.vdot(proton['state'],J@psi)),
            'Lambda_to_proton_axial_amplitude':complex(np.vdot(proton['state'],Jz@psi)),
            'exclusive_proton_probability':float(abs(np.vdot(proton['state'],success))**2),
            'all_flavour_records':records,
            'complete_original_native_occupation_records':{'preparation_success':scale,'preparation_failure':1-scale,
                'success_and_current_success_occupation':(scale*abs(raw_success)**2).tolist(),
                'success_and_current_failure_occupation':(scale*abs(raw_failure)**2).tolist(),
                'outside_reducing_source':'identity failure; prepared probability0'},
            'work_source_units':{'current_conditional':current_energy-L['energy'],
                                 'occupation_readout_conditional':measured_energy-current_energy,
                                 'unconditional_total':scale*(measured_energy-L['energy'])},
            'reference':'An initially charge-balanced apparatus lowers its charge by one on current success and is unchanged on failure; its declared number energy is zero.',
            'scope':'A complete native charged-current intervention on every three-quark colour singlet. This is not an autonomous semileptonic lifetime or outgoing neutrino computation.'}


if __name__=='__main__':
    import json
    a=instrument();b=instrument(F(0))
    assert a['exclusive_proton_probability']>0 and b['exclusive_proton_probability']<1e-25
    print(json.dumps({'selected':a,'zero_Cabibbo_control_exclusive_proton':b['exclusive_proton_probability']},default=str),flush=True)
