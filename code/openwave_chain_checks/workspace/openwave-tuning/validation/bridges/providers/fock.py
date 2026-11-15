"""Complete CAR occupations and exact finite particle/hole readouts."""
from functools import lru_cache
from itertools import product
import sympy as s


def act(mask,index,create,size):
    if type(size) is not int or size<1 or type(mask) is not int or not 0<=mask<1<<size:
        raise ValueError('finite full occupation mask required')
    if type(index) is not int or index not in range(size) or type(create) is not bool:
        raise ValueError('mode index and creation flag required')
    occupied=bool(mask&(1<<index))
    if occupied==create:return None
    sign=(-1)**((mask&((1<<index)-1)).bit_count())
    return mask^(1<<index),sign


@lru_cache(maxsize=32)
def operator(size,index,create):
    matrix=s.zeros(1<<size)
    for mask in range(1<<size):
        result=act(mask,index,create,size)
        if result is not None:
            target,sign=result;matrix[target,mask]=sign
    return matrix


def second_quantize(matrix):
    n=matrix.rows
    if matrix.cols!=n or n>6:raise ValueError('square exact control block of at most six modes required')
    return sum((matrix[i,j]*operator(n,i,True)*operator(n,j,False)
                for i,j in product(range(n),repeat=2) if matrix[i,j]!=0),s.zeros(1<<n))


def exterior_matrix(unitary,number):
    """Exact degree-N exterior lift, retaining every occupation of that degree."""
    from itertools import combinations
    n=unitary.rows
    if unitary.cols!=n or type(number) is not int or number not in range(n+1):
        raise ValueError('square matrix and admitted occupation number required')
    states=tuple(combinations(range(n),number))
    return s.Matrix([[unitary.extract(rows,columns).det() for columns in states] for rows in states])


def spectral_projector(number_operator,value,largest):
    result=s.eye(number_operator.rows)
    for other in range(largest+1):
        if other!=value:result=result*(number_operator-other*s.eye(number_operator.rows))/(value-other)
    return result.applyfunc(s.simplify)


def rest_instrument():
    from .spectrum import symbol,energy,projectors
    a=s.Rational(1,2);p=(0,0,0);h=symbol(a,p);e=energy(a,p)
    plus,minus=projectors(a,p)
    total=second_quantize(s.eye(4));np=second_quantize(plus)
    nh=2*s.eye(16)-second_quantize(minus)
    raw=second_quantize(h);positive=raw+2*e*s.eye(16)
    assert positive==e*(np+nh)
    ground=positive.nullspace()
    assert len(ground)==1
    ground=ground[0];norm=(ground.H*ground)[0]
    assert total*ground==2*ground
    creator=operator(4,0,True);annihilator=creator.H;n0=creator*annihilator
    # A complete, number-conserving fermionic SWAP with a reservoir mode.
    # Its recorded initial reservoir setting has probabilities 1/2,1/2.
    kraus={'add_from_reservoir':creator/s.sqrt(2),
           'remove_to_reservoir':annihilator/s.sqrt(2),
           'empty_no_change':(s.eye(16)-n0)/s.sqrt(2),
           'occupied_no_change':-n0/s.sqrt(2)}
    assert sum((k.H*k for k in kraus.values()),s.zeros(16))==s.eye(16)
    projectors_p=[spectral_projector(np,j,2) for j in range(3)]
    projectors_h=[spectral_projector(nh,j,2) for j in range(3)]
    pvms={(j,k):projectors_p[j]*projectors_h[k] for j,k in product(range(3),repeat=2)}
    assert sum(pvms.values(),s.zeros(16))==s.eye(16)
    rows=[]
    for setting,kraus_matrix in kraus.items():
        state=kraus_matrix*ground
        for (particles,holes),projector in pvms.items():
            probability=s.simplify((state.H*projector*state)[0]/norm)
            rows.append({'preparation_record':setting,'particles':particles,'holes':holes,
                         'relative_pair_charge':particles-holes,
                         'excitation_energy':str(e*(particles+holes)),
                         'probability':str(probability)})
    assert sum(s.Rational(row['probability']) for row in rows)==1
    populated=[row for row in rows if row['probability']!='0']
    assert len(populated)==6
    return {'spin_modes':4,'full_Fock_dimension':16,
            'sea_occupation':2,'sea_energy':str(-2*e),
            'positive_excitation_gap':str(e),'fixed_charge_neutral_gap':str(2*e),
            'unique_filled_sea_ground':True,
            'complete_preparation_and_particle_hole_records':rows,
            'retained_zero_probability_records':len(rows)-len(populated),
            'one_add_and_one_remove_record_probability':'1/4 each',
            'source_reservoir_number_conserved':True}
