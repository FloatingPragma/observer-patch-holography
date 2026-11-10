"""Exact spin and finite routing-current identities on the owned Clifford coin."""
from itertools import product
from functools import lru_cache
import sympy as s

from ..spatial.dirac_route import clifford, symbol


def clean(matrix):
    return s.ImmutableMatrix(matrix.applyfunc(s.simplify))


@lru_cache(maxsize=1)
def spin_data():
    alpha,beta=clifford()
    spin=tuple(clean(-s.I*alpha[(i+1)%3]*alpha[(i+2)%3]/2) for i in range(3))
    chirality=clean(-s.I*alpha[0]*alpha[1]*alpha[2])
    return alpha,beta,spin,chirality


def rotation(axis,cos_half,sin_half):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('a coordinate rotation axis is required')
    c,d=map(s.sympify,(cos_half,sin_half))
    if (any(x.has(s.Float) or x.free_symbols or x.is_real is not True
            or x.is_finite is not True for x in (c,d)) or s.simplify(c*c+d*d-1)!=0):
        raise ValueError('exact real unit-circle half-angle coordinates required')
    _,_,spin,_=spin_data()
    unitary=clean(c*s.eye(4)-2*s.I*d*spin[axis])
    cross=s.Matrix(3,3,lambda i,j:s.LeviCivita(i,axis,j))
    coordinate=s.eye(3)+(2*c*d)*cross+2*d*d*(cross*cross)
    return unitary,clean(coordinate)


def one_axis_shift(state,axis,length):
    """Sparse position/coin execution of one actual directional route gate."""
    if type(axis) is not int or axis not in range(3) or type(length) is not int or length<3:
        raise ValueError('one axis and a periodic length at least three required')
    alpha=spin_data()[0][axis]
    out={}
    for (site,coin),value in state.items():
        if type(site) is not int or not 0<=site<length or type(coin) is not int or not 0<=coin<4:
            raise ValueError('a site and coin basis value are required')
        value=s.sympify(value)
        if (value.has(s.Float) or value.free_symbols or value.is_finite is not True):
            raise ValueError('finite exact amplitudes required')
        for direction in (1,-1):
            projector=(s.eye(4)+direction*alpha)/2
            for row in range(4):
                key=((site+direction)%length,row)
                out[key]=s.expand(out.get(key,0)+projector[row,coin]*value)
    return {key:value for key,value in out.items() if value!=0}


def density(state,site,matrix=None):
    v=s.Matrix([state.get((site,j),0) for j in range(4)])
    return s.simplify((v.H*(s.eye(4) if matrix is None else matrix)*v)[0])


def continuity_certificate():
    # A coherent exact input with support on three sites and nontrivial phases.
    initial={(0,0):s.Integer(1),(1,1):s.I,(2,2):s.Rational(1,2),(2,3):-s.I/2}
    normalization=s.Rational(5,2);rows=[]
    for axis in range(3):
        alpha=spin_data()[0][axis]
        forward=(s.eye(4)+alpha)/2;backward=(s.eye(4)-alpha)/2
        out=one_axis_shift(initial,axis,4)
        flux=[s.simplify((density(initial,x,forward)-density(initial,(x+1)%4,backward))/normalization)
              for x in range(4)]
        before=[density(initial,x)/normalization for x in range(4)]
        after=[density(out,x)/normalization for x in range(4)]
        assert sum(before)==sum(after)==1
        assert all(s.simplify(after[x]-before[x]-flux[(x-1)%4]+flux[x])==0 for x in range(4))
        rows.append({'axis':axis,'position_before':list(map(str,before)),
                     'position_after':list(map(str,after)),
                     'directed_edge_flux':list(map(str,flux))})
    return rows


def certificate():
    alpha,beta,spin,chi=spin_data();checks=0
    for i,j in product(range(3),repeat=2):
        assert clean(spin[i]*spin[j]-spin[j]*spin[i])==sum(
            (s.I*s.LeviCivita(i,j,k)*spin[k] for k in range(3)),s.zeros(4))
        assert clean(spin[i]*alpha[j]-alpha[j]*spin[i])==sum(
            (s.I*s.LeviCivita(i,j,k)*alpha[k] for k in range(3)),s.zeros(4))
        checks+=2
    assert sum((j*j for j in spin),s.zeros(4))==s.Rational(3,4)*s.eye(4)
    assert chi==chi.H and chi*chi==s.eye(4)
    assert all(chi*a==a*chi for a in alpha)
    assert chi*beta==-beta*chi
    assert all(j*beta==beta*j for j in spin)
    rotation_rows=[]
    for axis in range(3):
        u,r=rotation(axis,s.Rational(3,5),s.Rational(4,5))
        assert clean(u.H*u)==s.eye(4) and r.T*r==s.eye(3) and r.det()==1
        for j in range(3):
            assert clean(u*alpha[j]*u.H-sum((r[i,j]*alpha[i] for i in range(3)),s.zeros(4)))==s.zeros(4)
        assert clean(u*beta*u.H)==beta
        assert rotation(axis,-1,0)[0]==-s.eye(4)
        rotation_rows.append({'axis':axis,'half_angle_cos':'3/5','half_angle_sin':'4/5',
                              'coordinate_rotation':[[str(v) for v in r.row(i)] for i in range(3)]})
    momentum=(s.pi/2,s.pi/3,0)
    turn,coordinate=rotation(2,s.sqrt(2)/2,s.sqrt(2)/2)
    rotated=tuple(coordinate*s.Matrix(momentum))
    residual=clean(turn*symbol(momentum)*turn.H-symbol(rotated))
    residual_norm=s.simplify(s.trace(residual.H*residual))
    assert residual_norm>0
    # Pointwise boosts preserve bar(psi) psi, but do not preserve density.
    boost_checks=0
    currents=(s.eye(4),*alpha)
    for axis in range(3):
        boost=s.Rational(5,3)*s.eye(4)+s.Rational(4,3)*alpha[axis]
        assert clean(boost.H*beta*boost)==beta
        assert clean(boost.H*boost)!=s.eye(4)
        lorentz=s.eye(4)
        lorentz[0,0]=lorentz[axis+1,axis+1]=s.Rational(41,9)
        lorentz[0,axis+1]=lorentz[axis+1,0]=s.Rational(40,9)
        metric=s.diag(1,-1,-1,-1)
        assert lorentz.T*metric*lorentz==metric
        for mu in range(4):
            target=sum((lorentz[mu,nu]*currents[nu] for nu in range(4)),s.zeros(4))
            assert clean(boost.H*currents[mu]*boost-target)==s.zeros(4)
            boost_checks+=1
    return {'spin_commutator_identities':checks,'spin_squared':'3/4 I4',
            'two_pi_rotation':'-I4','exact_rotation_controls':rotation_rows,
            'chirality_eigenvalue_multiplicities':{'-1':2,'1':2},
            'massless_route_preserves_chirality':True,
            'beta_mass_interchanges_chiralities':True,
            'finite_ordered_route_quarter_turn_residual_HS_squared':str(residual_norm),
            'finite_axis_order_exact_continuous_rotation_symmetry':False,
            'free_rotation_covariance_error_bound':'t*a*(C(p)+C(R*p)); C(p)=sum_(r<s)|h_r*h_s|, h=(p1,p2,p3,m)',
            'finite_complete_position_and_flux_records':continuity_certificate(),
            'boost_preserves_beta_bilinear':True,
            'boost_preserves_density_bilinear':False,
            'exact_Lorentz_vector_current_identities':boost_checks,
            'pointwise_boost_is_equal_time_CAR_unitary':False,
            'native_occupation_coupling_is_automatically_Lorentz_scalar':False,
            'parent_wave_internal_j_identified_as_physical_spin':False}
