"""Finite diagnostics supplementing, but not replacing, the analytic proof."""
import numpy as np
from .algebra import finite_example


def hermitian_log(matrix):
    values, vectors = np.linalg.eigh(matrix)
    if not np.isfinite(matrix).all() or np.min(values) <= 0:
        raise ValueError('positive definite finite matrix required')
    return (vectors*np.log(values))@vectors.conj().T


def entropy(matrix):
    return float(-np.trace(matrix@hermitian_log(matrix)).real)


def relative_entropy(left, right):
    return float(np.trace(left@(hermitian_log(left)-hermitian_log(right))).real)


def quasifree_density(covariance):
    if covariance.shape != (2, 2):
        raise ValueError('two regional modes required')
    d=np.linalg.det(covariance)
    out=np.zeros((4,4),complex)
    out[0,0]=np.linalg.det(np.eye(2)-covariance)
    out[3,3]=d
    out[1:3,1:3]=covariance-d*np.eye(2)
    return out


def finite_entropy_checks():
    e=finite_example();t=e['theta']
    A=np.array(e['A'],dtype=float);F=np.array(e['F'],dtype=float).reshape(2)
    base=np.array(e['rho'].subs(t,0),dtype=float)
    k=hermitian_log(np.eye(2)-A)-hermitian_log(A)
    derivative=-float(F@k@F)
    max_wick=max_gibbs=max_difference=0.0
    for theta in (-0.5,-0.125,0.125,0.5):
        rho=np.array(e['rho'].subs(t,theta),dtype=float)
        covariance=A-theta*np.outer(F,F)
        max_wick=max(max_wick,float(np.linalg.norm(rho-quasifree_density(covariance))))
        d=relative_entropy(rho,base)
        max_gibbs=max(max_gibbs,abs(entropy(rho)-entropy(base)-theta*derivative+d))
        if d < -1e-12 or d > theta*theta+1e-12:
            raise ArithmeticError('finite binary-channel entropy bound failed')
    for theta in (0.5,0.125,0.03125):
        plus=np.array(e['rho'].subs(t,theta),dtype=float)
        minus=np.array(e['rho'].subs(t,-theta),dtype=float)
        difference=(entropy(plus)-entropy(minus))/(2*theta)
        max_difference=max(max_difference,abs(difference-derivative)/theta)
    if max_wick>1e-12 or max_gibbs>1e-12 or max_difference>1:
        raise ArithmeticError('finite entropy identity failed')
    # A genuinely distinct raw cell norm. It must divide the raw covariance
    # update, not be replaced by the occupation probability of that mode.
    raw_norm=9/16
    raw_F=np.sqrt(raw_norm)*F
    normalized=np.outer(raw_F,raw_F)/raw_norm
    wrong=np.outer(raw_F,raw_F)/float(e['occupation'])
    if np.linalg.norm(normalized-np.outer(F,F))>1e-14 or np.linalg.norm(wrong-normalized)<1e-3:
        raise ArithmeticError('normalization negative control is ineffective')
    return {'full_16_dimensional_Fock_reduction_checked':True,
            'rank_one_quasifree_density_checked':True,
            'finite_Gibbs_entropy_identity_checked':True,
            'two_sided_entropy_difference_checked':True,
            'occupation_in_place_of_cell_norm_rejected':True,
            'numerical_identity_absolute_tolerance':'1e-12',
            'numerical_values_used_as_proof':False}


def endpoint_checks():
    count=0
    for x in (1e-12,1e-8,1e-4,0.1,0.5,0.9,1-1e-8):
        for delta in (1e-10,1e-6,0.001,0.25):
            k=np.log1p(-x)-np.log(x)
            kd=np.log(1-x+delta)-np.log(x+delta)
            if np.sqrt(x*(1-x))*abs(k-kd)>2*np.sqrt(delta)+1e-14:
                raise ArithmeticError('weighted endpoint inequality failed')
            if np.sqrt(x*(1-x))*abs(k)>2:
                raise ArithmeticError('natural logarithm domain bound failed')
            count+=1
    # The same assertion for arbitrary unweighted probes is false.
    x=1e-12;delta=1e-6
    if abs(np.log1p(-x)-np.log(x)-np.log(1-x+delta)+np.log(x+delta))<=2*np.sqrt(delta):
        raise ArithmeticError('unweighted endpoint negative control failed')
    return {'weighted_endpoint_pairs_checked':count,
            'arbitrary_unweighted_probe_claim_rejected':True,
            'uniform_covariance_gap_used':False}


def projector_checks():
    sigma=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    count=0
    for bits in range(8):
        signs=np.array([(-1)**((bits>>j)&1) for j in range(3)])
        for raw in ((1,0,0),(1,2,3),(-2,1,4),(1,-1,1)):
            direction=np.array(raw,dtype=float);direction/=np.linalg.norm(direction)
            for radius in (0.01,0.2,0.5,1.0):
                p=radius*direction
                d=np.sin(p)
                exact=(np.eye(2)-np.einsum('j,jkl->kl',signs*p/np.linalg.norm(p),sigma))/2
                lattice=(np.eye(2)-np.einsum('j,jkl->kl',signs*d/np.linalg.norm(d),sigma))/2
                if np.linalg.norm(exact-lattice,2)>radius*radius/6+1e-14:
                    raise ArithmeticError('all-valley low-frequency projector bound failed')
                count+=1
    return {'all_eight_valleys_checked':True,'finite_low_frequency_pairs':count,
            'high_frequency_aliases_removed_from_source':False}


def certificate():
    return {'entropy':finite_entropy_checks(),'endpoints':endpoint_checks(),
            'projectors':projector_checks()}
