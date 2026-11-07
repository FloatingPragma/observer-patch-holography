"""Complete recorded interferometry for an actual interacting addition response.

Observers infer response features from retained outcomes. A list of Hamiltonian
eigenvalues alone is never identified as particles. Measurements are configured
settings on the same full carrier and include their coherent back-action.
"""
from __future__ import annotations
from fractions import Fraction
import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse import issparse,block_diag,csr_matrix
from scipy.sparse.linalg import expm_multiply


def dense(a):
    return a.toarray() if issparse(a) else np.asarray(a)


def hermitian(a):
    a = csr_matrix(a, dtype=complex,copy=True)
    a.sum_duplicates();a.eliminate_zeros()
    if a.shape[0] != a.shape[1] or not a.shape[0] or not np.isfinite(a.data).all():
        raise ValueError('A finite nonempty square Hamiltonian is required')
    with np.errstate(over='ignore',invalid='ignore'):
        scale=np.linalg.norm(a.data)
        error=np.linalg.norm((a-a.getH()).data)
    if not np.isfinite(scale) or not np.isfinite(error):
        raise ValueError('Hamiltonian magnitude exceeds the numerical domain')
    if error > 1e-11*max(1.,scale):
        raise ValueError('Hamiltonian is not Hermitian')
    return a


def partial_isometry(b, n0, n1):
    b = csr_matrix(b, dtype=complex,copy=True)
    b.sum_duplicates();b.eliminate_zeros()
    if b.shape != (n1,n0) or not np.isfinite(b.data).all():
        raise ValueError('Finite creation probe with the sector dimensions required')
    gram = b.getH()@b
    with np.errstate(over='ignore',invalid='ignore'):
        scale=np.linalg.norm(gram.data)
        error=np.linalg.norm((gram@gram-gram).data)
    if not np.isfinite(scale) or not np.isfinite(error):
        raise ValueError('Creation probe magnitude exceeds the numerical domain')
    if error > 1e-10*max(1.,scale):
        raise ValueError('The quarter-turn probe requires a partial isometry')
    return b


def checked_spectrum(h, supplied):
    if supplied is None:
        return eigh(h)
    e,u = map(np.asarray,supplied)
    n = len(h)
    if (e.shape != (n,) or u.shape != (n,n) or not np.isrealobj(e)
            or not np.isfinite(e).all() or not np.isfinite(u).all()
            or np.any(np.diff(e)<0)):
        raise ValueError('A complete ordered finite eigensystem is required')
    with np.errstate(over='ignore',invalid='ignore'):
        orthogonality=np.linalg.norm(u.conj().T@u-np.eye(n),'fro')
        residual=np.linalg.norm(h@u-u*e,'fro')
    if (not np.isfinite(orthogonality) or not np.isfinite(residual)
            or orthogonality > 1e-9*max(1.,np.sqrt(n))
            or residual > 1e-9*max(1.,np.linalg.norm(h,'fro'))):
        raise ValueError('Supplied eigensystem does not reconstruct the Hamiltonian')
    return e,u


def spectral_matrix(h, vectors, reference_energy=0.0):
    """Matrix spectral measure of all supplied detector-created states."""
    if (np.ndim(reference_energy)!=0 or isinstance(reference_energy,(bool,np.bool_,complex,np.complexfloating))
            or not np.isfinite(float(reference_energy))):
        raise ValueError('Finite real reference energy required')
    e,u=eigh(dense(hermitian(h)))
    vectors=np.asarray(vectors,complex)
    if vectors.ndim != 2 or vectors.shape[0] != len(e) or not np.isfinite(vectors).all():
        raise ValueError('Finite detector columns on the full carrier required')
    amplitudes=u.conj().T@vectors
    weights=np.einsum('ni,nj->nij',amplitudes.conj(),amplitudes)
    gram=vectors.conj().T@vectors
    error=float(np.linalg.norm(weights.sum(axis=0)-gram,2))
    if error>1e-10:raise ArithmeticError('Complete detector addition sum rule')
    return {'energies':e,'frequencies':e-reference_energy,'basis':u,
            'amplitudes':amplitudes,'matrix_weights':weights,
            'Gram':gram,'sum_rule_error':error}


def record_law(h0,h1,creation,times, *, spectral0=None,spectral1=None):
    """Two-path, forward-only complete Ramsey records for every detector pair.

    creation[j]: H0-sector -> H1-sector, with an explicit compensating detector
    reference. Each creation operator must be a partial isometry.
    W_s = I-G_j^2-i*s*G_j is its quarter-turn complete unitary on H0+H1.
    A sign difference of complete detector records extracts the addition
    correlator; the blocked component is retained and cancels algebraically.
    """
    h0,h1=dense(hermitian(h0)),dense(hermitian(h1))
    e0,u0=checked_spectrum(h0,spectral0)
    e1,u1=checked_spectrum(h1,spectral1)
    if len(e0)>1 and e0[1]-e0[0]<1e-10:raise ValueError('This run needs a specified unique prepared ground or a density extension')
    ground=u0[:,0];energy=float(e0[0])
    creation=[dense(partial_isometry(x,len(h0),len(h1))) for x in creation]
    if not creation:raise ValueError('At least one detector probe is required')
    excited=np.column_stack([b@ground for b in creation])
    blocked=np.column_stack([ground-b.conj().T@(b@ground) for b in creation])
    aa=u1.conj().T@excited;bb=u0.conj().T@blocked
    raw_times=np.asarray(times)
    if (raw_times.ndim!=1 or np.iscomplexobj(raw_times)
            or any(isinstance(t,(bool,np.bool_,complex,np.complexfloating)) for t in times)):
        raise ValueError('Real observation times required')
    tau=np.asarray(times,dtype=float)
    if tau.ndim != 1 or not len(tau) or not np.isfinite(tau).all() or np.any(tau<0):
        raise ValueError('A nonempty finite forward observation schedule is required')
    addition=np.einsum('ni,tn,nj->tij',aa.conj(),np.exp(-1j*tau[:,None]*(e1-energy)),aa,optimize=True)
    background=np.einsum('ni,tn,nj->tij',bb.conj(),np.exp(-1j*tau[:,None]*(e0-energy)),bb,optimize=True)
    signals=np.stack((background+addition,background-addition),axis=1)
    # axes: time, relative sign (+,-), detector i, detector j, quadrature(X,Y), outcome(+,-)
    means=np.stack((signals.real,signals.imag),axis=-1)
    probabilities=np.stack(((1+means)/2,(1-means)/2),axis=-1)
    if probabilities.min() < -2e-10 or probabilities.max()>1+2e-10:
        raise ArithmeticError('A complete Ramsey probability left [0,1]')
    raw_probabilities=probabilities.copy()
    probabilities=np.clip(probabilities,0,1)
    clipping=float(np.max(abs(probabilities-raw_probabilities)))
    reconstructed=(signals[:,0]-signals[:,1])/2
    gram=excited.conj().T@excited
    partial_errors=[float(np.linalg.norm((b.conj().T@b)@(b.conj().T@b)-b.conj().T@b,2)) for b in creation]
    if max(partial_errors)>2e-9:raise ArithmeticError('Literal detector maps are not partial isometries')
    return {'times':tau,'addition_correlator':addition,'blocked_correlator':background,
            'probabilities':probabilities,'addition_vectors':excited,'blocked_vectors':blocked,
            'ground':ground,'ground_energy':energy,'energies0':e0,'basis0':u0,
            'energies1':e1,'basis1':u1,'amplitudes':aa,'Gram':gram,
            'checks':{'complete_probability_sum_error':float(np.max(abs(probabilities.sum(axis=-1)-1))),
                      'interference_sign_difference_error':float(np.max(abs(reconstructed-addition))),
                      'maximum_probability_clipping_correction':clipping,
                      'partial_isometry_errors':partial_errors,
                      'addition_weight_plus_blocked_weight_error':float(np.max(abs(np.diag(gram+blocked.conj().T@blocked)-1))),
                      'all_complementary_detector_outcomes_retained':True,
                      'postselection':False,'backward_time_evolution_required':False}}


def ramsey_kraus_actions(h0,h1,creation_left,creation_right,state,time,
                        relative_sign=1,quadrature="X"):
    """Apply both Kraus maps to any vector (or columns of a purification).

    Returns unnormalized postmeasurement vectors. Their squared norms are
    the probabilities when the input is normalized; both outputs are kept.
    M^X_±=(W_j U ± U W_sk)/2 and M^Y_±=(W_j U ∓ i U W_sk)/2.
    The caller may normalize an individual output only after recording it.
    """
    if type(relative_sign) is not int or relative_sign not in (-1,1) or quadrature not in ("X","Y"):
        raise ValueError("relative_sign must be +/-1 and quadrature X or Y")
    h0,h1=hermitian(h0),hermitian(h1)
    n0,n1=h0.shape[0],h1.shape[0]
    bj=partial_isometry(creation_left,n0,n1)
    bk=partial_isometry(creation_right,n0,n1)
    state=np.asarray(state,complex)
    if state.ndim not in (1,2) or state.shape[0]!=n0+n1 or not np.isfinite(state).all():
        raise ValueError("Finite full direct-sum input required")
    if (np.ndim(time)!=0 or isinstance(time,(bool,np.bool_,complex,np.complexfloating))
            or not np.isfinite(float(time)) or float(time)<0):
        raise ValueError('A finite nonnegative time is required')
    def quarter(b,z,sgn):
        u,v=z[:n0],z[n0:]
        return np.concatenate((u-b.getH()@(b@u)-1j*sgn*(b.getH()@v),
                               v-b@(b.getH()@v)-1j*sgn*(b@u)),axis=0)
    generator=-1j*float(time)*block_diag((h0,h1),format="csr")
    path0=quarter(bj,expm_multiply(generator,state),1)
    path1=expm_multiply(generator,quarter(bk,state,relative_sign))
    if quadrature=="X":return ((path0+path1)/2,(path0-path1)/2)
    return ((path0-1j*path1)/2,(path0+1j*path1)/2)


def sample_diagonal_records(law,shots=4096,seed=20261108):
    """Retain every binary result in packed form, plus sufficient counts."""
    if type(shots) is not int or shots<=0:raise ValueError('Positive integer shot count required')
    p=law['probabilities'];nt,ns,nf,_,nq,no=p.shape
    probabilities=np.empty((nt,ns,nf,nq))
    for j in range(nf):probabilities[:,:,j,:]=p[:,:,j,j,:,0]
    rng=np.random.default_rng(seed)
    records=rng.random((*probabilities.shape,shots))<probabilities[...,None]
    counts=records.sum(axis=-1)
    means=2*counts/shots-1
    estimated=(means[:,0]-means[:,1])/2
    truth=np.stack((np.diagonal(law['addition_correlator'],axis1=1,axis2=2).real,
                    np.diagonal(law['addition_correlator'],axis1=1,axis2=2).imag),axis=-1)
    # Hoeffding and a union bound cover every recorded +/- mean. Taking half
    # their difference preserves the same uniform radius for each quadrature.
    alpha=.01
    radius=float(np.sqrt(2*np.log(2*probabilities.size/alpha)/shots))
    return {'packed_binary_records':np.packbits(records,axis=-1),
            'plus_counts':counts,'plus_probabilities':probabilities,
            'estimated_addition_quadratures':estimated,
            'exact_addition_quadratures':truth,
            'receipt':{'seed':seed,'shots_per_setting':shots,
                       'binary_records_retained':int(records.size),
                       'settings':int(probabilities.size),
                       'simulated_records_not_experimental_data':True,
                       'joint_Hoeffding_failure_probability_upper':alpha,
                       'uniform_quadrature_error_radius':radius,
                       'observed_maximum_quadrature_error':float(np.max(abs(estimated-truth))),
                       'all_outcomes_recorded':True,
                       'sample_error_is_not_a_physical_intrinsic_linewidth':True}}


def spectral_shift_tv_squared(times,shots,settings_per_time,shift):
    """Exact upper bound on squared record TV for a commuting sector shift.

    Independent freshly prepared trials, the fixed Ramsey instruments, and
    H1' = H1 - shift*P with P a spectral projector are required. The result
    bounds any subsequent processing of the retained records.
    """
    if (type(shots) is not int or shots<=0 or type(settings_per_time) is not int
            or settings_per_time<=0 or type(shift) is not Fraction):
        raise ValueError('Exact shift and positive integer trial counts required')
    times=tuple(times)
    if not times or any(type(t) is not Fraction for t in times):
        raise ValueError('A nonempty rational observation schedule is required')
    return min(Fraction(1),shift**2*shots*settings_per_time*sum(t*t for t in times)/4)


def independent_full_control():
    """Check phase/sign convention with full ancilla unitaries on a small CAR example."""
    h0=np.array([[.2,.13j],[-.13j,.9]])
    h1=np.array([[1.1,.2],[.2,1.8]])
    b=np.array([[1.,0.],[0.,0.]])
    t=.37
    law=record_law(h0,h1,[b],[t])
    g=np.block([[np.zeros_like(h0),b.conj().T],[b,np.zeros_like(h1)]])
    whole=np.block([[h0,np.zeros_like(h0)],[np.zeros_like(h1),h1]])
    psi=np.r_[law['ground'],np.zeros(2,complex)]
    u=expm(-1j*t*whole)
    x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    errors=[]
    for sidx,sgn in enumerate((1,-1)):
        w0=expm(-1j*np.pi*g/2)
        w1=expm(-1j*sgn*np.pi*g/2)
        branch0=w0@u@psi
        branch1=u@w1@psi
        final=np.r_[branch0,branch1]/np.sqrt(2)
        for qidx,pauli in enumerate((x,y)):
            mean=float(np.vdot(final,np.kron(pauli,np.eye(4))@final).real)
            expected=float(law['probabilities'][0,sidx,0,0,qidx,0]*2-1)
            errors.append(abs(mean-expected))
    if max(errors)>1e-12:raise ArithmeticError('Recorded X/Y sign or forward-order mismatch')
    return {'full_joint_dimension':8,'maximum_probability_mean_error':max(errors),
            'complete_forward_only_ancilla_control_checked':True}
