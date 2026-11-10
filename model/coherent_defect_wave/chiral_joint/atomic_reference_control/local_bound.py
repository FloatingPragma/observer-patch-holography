"""Unrestricted collective weight blocks, retaining all charged local sectors."""
from collections import defaultdict,Counter
from functools import lru_cache
import sympy as s
from ..local_reference.reduction import quark_module,lepton_operators
from ..local_reference.physical import labels,columns
from ..local_reference.characters import inventory

SQUARED_BOUND=s.Rational(595,4)


@lru_cache(maxsize=1)
def blank_blocks():
    qweights,lweights,_=labels();blocks=defaultdict(list)
    for q,(qz,qy) in enumerate(qweights):
        for mask,(n,lz,ly,spin) in enumerate(lweights):
            blocks[n,spin,qz+lz,qy+ly].append((q,mask))
    assert sum(map(len,blocks.values()))==140*64
    assert len(blocks)==1686 and max(map(len,blocks.values()))==30
    return {label:tuple(keys) for label,keys in sorted(blocks.items())}


@lru_cache(maxsize=1)
def transition_columns():
    q=quark_module();l=lepton_operators()
    assert all(q['metric']*q['transitions'][(r+2)%4]==q['transitions'][r].T*q['metric'] for r in range(4))
    assert all(l['transitions'][(r+2)%4]==l['transitions'][r].T for r in range(4))
    return tuple(columns(m) for m in q['transitions']),tuple(columns(m) for m in l['transitions'])


def transition(key,r):
    qc,lc=transition_columns();q,mask=key;out=defaultdict(lambda:s.Integer(0))
    for target,value in qc[r].get(q,()):out[target,mask]+=value
    for target,value in lc[r].get(mask,()):out[q,target]+=value
    return {key:s.expand(value) for key,value in out.items() if value}


def squared_column(key):
    out=defaultdict(lambda:s.Integer(0))
    for r in range(4):
        for middle,value in transition(key,r).items():
            for target,weight in transition(middle,(r+2)%4).items():out[target]+=value*weight
    return {key:s.expand(value) for key,value in out.items() if value}


@lru_cache(maxsize=None)
def exact_block(label):
    keys=blank_blocks()[label];positions={key:i for i,key in enumerate(keys)};entries={};gram={}
    qg=columns(quark_module()['metric'])
    for j,key in enumerate(keys):
        for target,value in squared_column(key).items():
            assert target in positions,('omitted unrestricted squared-transfer output',label,key,target)
            entries[positions[target],j]=value
        for q,value in qg.get(key[0],()):
            target=(q,key[1]);assert target in positions
            gram[positions[target],j]=value
    transfer=s.ImmutableSparseMatrix(len(keys),len(keys),entries)
    metric=s.ImmutableSparseMatrix(len(keys),len(keys),gram)
    form=s.ImmutableSparseMatrix(metric*transfer)
    assert form==form.T
    return keys,metric,form,transfer


def positive_pivots(matrix):
    """Exact rational LDL pivots; rejects zero or negative directions."""
    if matrix!=matrix.T:raise ValueError('a symmetric rational form is required')
    a=s.MutableDenseMatrix(matrix);n=a.rows;pivots=[]
    for k in range(n):
        pivot=a[k,k]
        if pivot.is_Rational is not True or pivot<=0:
            raise ValueError(('nonpositive exact LDL pivot',k,str(pivot)))
        pivots.append(pivot)
        for i in range(k+1,n):
            if not a[i,k]:continue
            factor=a[i,k]/pivot
            for j in range(i,n):
                a[j,i]-=factor*a[j,k]
                if j!=i:a[i,j]=a[j,i]
    return tuple(pivots)


@lru_cache(maxsize=1)
def maximal_module_bound():
    rows=[]
    for label in blank_blocks():
        keys,g,k,h=exact_block(label)
        gp=positive_pivots(g)
        positive=positive_pivots(SQUARED_BOUND*g-k)
        rows.append({'transformed_lepton_number':label[0],'two_spin':label[1],
            'twice_weak_weight':label[2],'twice_hypercharge_weight':label[3],
            'complete_scalar_blank_dimension':len(keys),
            'Gram_LDL_pivots':[str(p) for p in gp],
            'upper_bound_LDL_pivots':[str(p) for p in positive]})
    assert sum(row['complete_scalar_blank_dimension'] for row in rows)==8960
    return {'quark_module_dimension':140,'quark_Casimir':'54','lepton_Fock_dimension':64,
            'full_scalar_states':5,'complete_module_dimension':44800,
            'complete_scalar_blank_dimension':8960,'complete_scalar_excited_dimension':35840,
            'unrestricted_weight_blocks':rows,'largest_block':30,
            'all_nonzero_weak_and_hypercharge_weights_retained':True,
            'Gauss_kernel_or_highest_weight_projection_applied':False,
            'metric_self_adjoint_squared_transfer':True,
            'strict_squared_transfer_upper_bound':str(SQUARED_BOUND),
            'exact_all_block_positive_definiteness_verified':True}


def other_modules_bound():
    rows=inventory();allmodules=[(row['transformed_quark_number'],item) for row in rows for item in row['modules']]
    assert sum(row['full_dimension'] for row in rows)==2**24
    maximal=[(n,item) for n,item in allmodules if s.Rational(item['Casimir'])==54]
    assert len(maximal)==1 and maximal[0][0]==12 and maximal[0][1]['dimension']==140 and maximal[0][1]['multiplicity']==1
    others=[s.Rational(item['Casimir']) for _,item in allmodules if s.Rational(item['Casimir'])!=54]
    assert max(others)==s.Rational(95,2)
    # (sqrt95+sqrt6)^2=101+2sqrt570 <595/4.
    assert SQUARED_BOUND>101 and ((SQUARED_BOUND-101)/2)**2>570
    return {'full_quark_Fock_dimension':2**24,'full_local_CAR_dimension':2**30,
            'complete_character_inventory_reused':True,'color_singlet_filter_applied':False,
            'largest_nonmaximal_Casimir':'95/2','all_nonmaximal_quark_squared_transfer_bound':'95',
            'complete_lepton_squared_transfer_bound':'6',
            'joint_squared_transfer_bound':'101+2 sqrt(570)',
            'strict_rational_squared_transfer_upper_bound':str(SQUARED_BOUND),
            'charged_and_nonzero_baryon_sectors_retained':True}


def certificate():
    maximal=maximal_module_bound();other=other_modules_bound()
    lower=s.Rational(-1171,100)
    assert lower<0 and lower*(lower-1)>SQUARED_BOUND
    return {'maximal_quark_module':maximal,'every_other_quark_module':other,
            'full_unconstrained_local_dimension':5*2**30,
            'exact_scalar_star_Hamiltonian':'[[0,B^dagger],[B,I]]',
            'all_unconstrained_squared_transfer_strict_upper_bound':str(SQUARED_BOUND),
            'algebraic_unconstrained_energy_lower_bound':'(1-sqrt(596))/2',
            'strict_rational_unconstrained_energy_lower_bound':str(lower),
            'all_scalar_dark_states_retained':True,'all_local_Gauss_charge_sectors_covered':True,
            'numerical_eigensolver_used':False,'exact_unconstrained_ground_energy_identified':False}
