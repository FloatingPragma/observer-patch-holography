"""The complete uniform two-cell lepton/link/two-Higgs Hamiltonian.

Every column belongs to the full 614400-state carrier.  No vacuum-cyclic
restriction or scalar replacement enters this operator.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import sympy as s
from ..model import active_terms, sequence, clean

DIMENSION = 2**12 * 6 * 5**2
PARAMETERS = {'left_y_e':'1','right_y_e':'1','omega':'1','kappa':'1',
              'y_u':'0','y_d':'0','spectral_z_edge_coefficient':'I*sigma3/2'}


def validate_key(key):
    if len(key)!=4 or not all(isinstance(n,int) for n in key):
        raise ValueError('key must contain four exact integer coordinates')
    mask,link,hl,hr=key
    if not (0<=mask<4096 and 0<=link<6 and 0<=hl<5 and 0<=hr<5):
        raise ValueError('key outside full two-cell carrier')
    return key


def index(key):
    mask,link,hl,hr=validate_key(key)
    return ((mask*6+link)*5+hl)*5+hr


def from_index(i):
    if not isinstance(i,int) or not 0<=i<DIMENSION:
        raise ValueError('index outside full two-cell carrier')
    t,hr=divmod(i,5);t,hl=divmod(t,5);mask,link=divmod(t,6)
    return mask,link,hl,hr


def conserved_labels(key):
    mask,_,_,_=validate_key(key)
    lepton=sum((-1 if j%3==0 else 1) for j in range(12) if mask>>j&1)
    twice_spin=sum((1 if j%6<3 else -1) for j in range(12) if mask>>j&1)
    return lepton,twice_spin


@lru_cache(maxsize=320)
def local_yukawa(mask,scalar):
    """Both actual charged-lepton tensors on every six-CAR/scalar column."""
    _,cs,fields,_=active_terms();out=defaultdict(lambda:s.Integer(0))
    for k,c in enumerate(cs):
        phi=fields[k].H
        for i,j in product(range(3),repeat=2):
            if not c[i,j]:continue
            operators=((i,False),(3+j,False))
            for matrix,ops,coef in ((phi,operators,c[i,j]),
                (phi.H,((3+j,True),(i,True)),s.conjugate(c[i,j]))):
                step=sequence(mask,ops)
                if step is None:continue
                for destination in range(5):
                    if matrix[destination,scalar]:
                        out[step[0],destination]+=coef*step[1]*matrix[destination,scalar]
    return clean(out)


@lru_cache(maxsize=4096*6)
def hopping(mask,link):
    _,_,_,labels=active_terms();out=defaultdict(lambda:s.Integer(0))
    for label,(a,b) in enumerate(labels[1:],1):
        for spin in range(2):
            coefficient=s.I*s.Integer(-1)**spin/2
            lmode=3*spin+a;rmode=6+3*spin+b
            if link==0:
                step=sequence(mask,((lmode,True),(rmode,False)))
                if step:out[step[0],label]+=coefficient*step[1]
            if link==label:
                step=sequence(mask,((rmode,True),(lmode,False)))
                if step:out[step[0],0]+=s.conjugate(coefficient)*step[1]
    return clean(out)


@lru_cache(maxsize=65536)
def column(key):
    mask,link,hl,hr=validate_key(key);out=defaultdict(lambda:s.Integer(0))
    out[key]+=int(link!=0)+int(hl!=0)+int(hr!=0)
    for (target,newlink),coefficient in hopping(mask,link).items():
        out[target,newlink,hl,hr]+=coefficient
    lm=mask&63;rm=mask>>6
    for (target,scalar),coefficient in local_yukawa(lm,hl).items():
        out[target+(rm<<6),link,scalar,hr]+=coefficient
    # The local Yukawa operator has even CAR parity, so its two offset
    # Jordan-Wigner factors cancel exactly at the right cell.
    for (target,scalar),coefficient in local_yukawa(rm,hr).items():
        out[lm+(target<<6),link,hl,scalar]+=coefficient
    result=clean(out)
    assert all(conserved_labels(target)==conserved_labels(key) for target in result)
    return result


def apply_state(state):
    out=defaultdict(lambda:s.Integer(0))
    for key,coefficient in state.items():
        for target,value in column(key).items():out[target]+=coefficient*value
    return clean(out)


def matrix_on(keys):
    """Exact restriction, rejecting any omitted outgoing coordinate."""
    keys=tuple(keys);positions={key:j for j,key in enumerate(keys)}
    if len(positions)!=len(keys):raise ValueError('duplicate coordinate')
    entries={}
    for j,key in enumerate(keys):
        for target,value in column(key).items():
            if target not in positions:raise ValueError(('noninvariant coordinate menu',key,target))
            entries[positions[target],j]=value
    h=s.ImmutableSparseMatrix(len(keys),len(keys),entries)
    assert h==h.H
    return h


@lru_cache(maxsize=1)
def structural_checks():
    """Complete elementary domains, including all occupation alternatives."""
    local_entries=0
    for mask in range(64):
        for scalar in range(5):
            for (target,dest),value in local_yukawa(mask,scalar).items():
                assert local_yukawa(target,dest).get((mask,scalar),0)==s.conjugate(value)
                assert conserved_labels((mask,0,scalar,0))==conserved_labels((target,0,dest,0))
                local_entries+=1
    hopping_entries=0
    for mask in range(4096):
        for link in range(6):
            for (target,dest),value in hopping(mask,link).items():
                assert hopping(target,dest).get((mask,link),0)==s.conjugate(value)
                assert conserved_labels((mask,link,0,0))==conserved_labels((target,dest,0,0))
                hopping_entries+=1
    from ...chiral_higgs.model import quantum_lepton_block
    from ..model import path_ward_checks
    local_h,_,_,_,local_ward=quantum_lepton_block()
    local_entries_matrix={}
    for mask in range(64):
        for scalar in range(5):
            j=5*mask+scalar
            if scalar:local_entries_matrix[j,j]=s.Integer(1)
            for (target,dest),value in local_yukawa(mask,scalar).items():
                local_entries_matrix[5*target+dest,j]=value
    assert s.ImmutableSparseMatrix(320,320,local_entries_matrix)==local_h
    path_ward=path_ward_checks()
    assert DIMENSION==614400
    return {'full_carrier_dimension':DIMENSION,'local_CAR_scalar_columns_checked':320,
            'hopping_CAR_link_columns_checked':24576,'exact_local_adjoint_entries':local_entries,
            'exact_hopping_adjoint_entries':hopping_entries,
            'both_local_quantum_Yukawa_terms_present':True,
            'full_H_is_Hermitian':True,'L_and_twoSz_conserved':True,
            'identical_to_primary_full_320_state_local_Higgs_operator':True,
            'exact_local_quantum_Gauss_commutators':local_ward,
            'native_quantum_path_Ward_checks_reused':path_ward['all_endpoint_and_internal_vertex_Ward_checks'],
            'full_operator_norm_upper_bound':'29',
            'parameters':PARAMETERS,'quark_occupations':'retained empty reducing spectators'}
