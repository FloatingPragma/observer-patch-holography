"""Exact finite scalar and closed two-link operators, with complete failures."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..model import active_terms,clean
from ..two_cell.basis import current_data

SETTINGS=('scalar_number','singlet_loop','doublet_loop')


@lru_cache(maxsize=None)
def operators(setting):
    if setting=='scalar_number':
        a=s.diag(0,1,1,1,1);return a,a,s.eye(5)-a
    if setting not in SETTINGS:raise ValueError('specified scalar or link-loop setting required')
    labels=active_terms()[3];block=(0,) if setting=='singlet_loop' else (1,2)
    r=s.zeros(36)
    for a in block:
        for b in block:r[6*labels.index((a,b))+labels.index((b,a)),0]=1
    a=(r+r.H)/len(block);a=s.ImmutableSparseMatrix(a)
    assert a.H==a and a**3==a and s.trace(a*a)==2
    success=a/2;failure=s.eye(36)+(s.sqrt(3)/2-1)*a*a
    assert (success.H*success+failure.H*failure).applyfunc(s.expand)==s.eye(36)
    return a,success,s.ImmutableSparseMatrix(failure)


def operator_certificate():
    _,scalar,edge=current_data();rows=[]
    for name in SETTINGS:
        a,k,f=operators(name);dimension=a.rows
        assert (k.H*k+f.H*f).applyfunc(s.expand)==s.eye(dimension)
        gauge=0
        if name=='scalar_number':
            for q in scalar:assert q*a==a*q;gauge+=1
            multiplicities={'0':1,'1':4};translation='local scalar moved between the two cells'
        else:
            for site in (0,1):
                for current in range(12):
                    q=s.kronecker_product(edge[site][current],s.eye(6))+s.kronecker_product(s.eye(6),edge[1-site][current])
                    assert q*a==a*q;gauge+=1
            swap=s.ImmutableSparseMatrix(36,36,{(6*b+a,6*a+b):1 for a in range(6) for b in range(6)})
            assert swap*a==a*swap
            multiplicities={'-1':1,'0':34,'1':1};translation='commutes with exact periodic translation'
        unitary=f.row_join(-k.H).col_join(k.row_join(f))
        assert (unitary.H*unitary).applyfunc(s.expand)==s.eye(2*dimension)
        rows.append({'setting':name,'complete_operator_dimension':dimension,
                     'Hermitian_operator_eigenvalue_multiplicities':multiplicities,
                     'Hermitian_operator_norm':'1','success_operator':'A' if name=='scalar_number' else 'A/2',
                     'failure_operator':'I-A' if name=='scalar_number' else 'I+(sqrt(3)/2-1) A^2',
                     'exact_native_current_commutators':gauge,
                     'translation_action':translation,'explicit_unitary_dilation_dimension':2*dimension,
                     'all_link_and_scalar_alternatives_retained':True})
    return {'settings':rows,'setting_probabilities':['1/3','1/3','1/3'],
            'loop_definition':'A_B=(sum_ab U01_ab U10_ba+h.c.)/d_B',
            'native_retained_block_dimensions':{'singlet_loop':1,'doublet_loop':2},
            'original_quantum_link_matrix_units_used':True,
            'colored_link_channels_included':False,'source_Hamiltonian_changed':False,
            'complete_instrument':'each selected setting retains success and failure',
            'apparatus_parity_and_gauge_charge':'even and gauge neutral'}


@lru_cache(maxsize=None)
def operator_column(key,setting,outcome='success'):
    a,k,f=operators(setting);op={'raw':a,'success':k,'failure':f}[outcome]
    out=defaultdict(lambda:s.Integer(0))
    coordinate=key[4] if setting=='scalar_number' else 6*key[1]+key[2]
    for (row,col),value in s.ImmutableSparseMatrix(op).todok().items():
        if col!=coordinate:continue
        target=list(key)
        if setting=='scalar_number':target[4]=row
        else:target[1],target[2]=divmod(row,6)
        out[tuple(target)]+=value
    return clean(out)
