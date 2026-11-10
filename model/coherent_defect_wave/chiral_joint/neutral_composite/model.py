"""Gauge-neutral composite, full-domain contraction and retained failure operator."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..higgs_composites.model import (neutral_composite,adjoint,multiply,add,scale,
    monomial,commutator,gauge_current,clean,column)
from ..model import sequence


@lru_cache(maxsize=None)
def local_instrument(kind='remove',spin=0):
    if kind not in ('add','remove') or spin not in (0,1):raise ValueError('specified neutral probe required')
    j=neutral_composite(spin)
    if kind=='add':j=adjoint(j)
    a=multiply(adjoint(j),j);identity=monomial()
    # Exact complete 2-CAR/scalar effect spectrum is {0,1,2}.
    p0=scale(multiply(add(a,scale(identity,-1)),add(a,scale(identity,-2))),s.Rational(1,2))
    p1=scale(multiply(a,add(a,scale(identity,-2))),-1)
    p2=scale(multiply(a,add(a,scale(identity,-1))),s.Rational(1,2))
    projectors=(p0,p1,p2)
    assert add(*projectors)==identity
    for k,p in enumerate(projectors):
        assert p==adjoint(p) and multiply(p,p)==p
        assert multiply(a,p)==scale(p,k)
        for other in projectors:
            if other is not p:assert not multiply(p,other)
    failure=add(p0,scale(p1,s.sqrt(3)/2),scale(p2,s.sqrt(2)/2))
    success=scale(j,s.Rational(1,2))
    assert add(multiply(adjoint(success),success),multiply(adjoint(failure),failure))==identity
    assert failure==adjoint(failure)
    return j,success,failure,projectors


def instrument_certificate():
    rows=[]
    for kind in ('remove','add'):
        j,success,failure,projectors=local_instrument(kind)
        keys=[(sum(1<<i for bit,i in enumerate((10,11)) if m>>bit&1),h) for m in range(4) for h in range(5)]
        positions={key:i for i,key in enumerate(keys)}
        matrices=[]
        for op in (j,success,failure,*projectors):
            entries={}
            for col,key in enumerate(keys):
                for target,value in column(op,*key).items():entries[positions[target],col]=value
            matrices.append(s.ImmutableSparseMatrix(20,20,entries))
        jm,sm,fm,*pm=matrices
        assert (sm.H*sm+fm.H*fm).applyfunc(s.expand)==s.eye(20)
        multiplicities=[int(s.trace(p)) for p in pm];assert multiplicities==[14,4,2]
        opposite=local_instrument('add' if kind=='remove' else 'remove')[2]
        reverse_entries={}
        for col,key in enumerate(keys):
            for target,value in column(opposite,*key).items():reverse_entries[positions[target],col]=value
        reverse=s.ImmutableSparseMatrix(20,20,reverse_entries)
        unitary=fm.row_join(-sm.H).col_join(sm.row_join(reverse))
        assert (unitary.H*unitary).applyfunc(s.expand)==s.eye(40)
        assert (unitary*unitary.H).applyfunc(s.expand)==s.eye(40)
        parity=s.diag(*[(-1)**mask.bit_count() for mask,scalar in keys])
        number=s.diag(*[mask.bit_count() for mask,scalar in keys])
        reservoir_charge=1 if kind=='remove' else -1
        total_number=s.diag(number,number+reservoir_charge*s.eye(20))
        total_parity=s.diag(parity,-parity)
        assert (unitary*total_number-total_number*unitary).applyfunc(s.expand)==s.zeros(40)
        assert (unitary*total_parity-total_parity*unitary).applyfunc(s.expand)==s.zeros(40)
        for a in range(12):
            assert not commutator(gauge_current(a),j)
            assert not commutator(gauge_current(a),failure)
        rows.append({'kind':kind,'local_complete_dimension':20,'effect_eigenvalues':['0','1','2'],
                     'effect_multiplicities':multiplicities,'exact_complete_local_columns':20,
                     'full_F30_CAR_identity':True,'failure_operator_positive':True,
                     'success_operator':'J/2','failure_operator':'P0+sqrt(3)/2 P1+sqrt(2)/2 P2',
                     'success_effect_norm_upper_bound':'1/2','all_twelve_gauge_currents_preserved':True,
                     'failure_computed_before_translation_projection':True,
                     'explicit_unitary_dilation_dimension':40,'two_sided_unitarity_checked':True,
                     'total_reservoir_source_parity_and_lepton_number_conserved':True,
                     'reservoir_success_lepton_charge':reservoir_charge,
                     'reservoir_success_spin_z':str(s.Rational(reservoir_charge,2))})
    return {'settings_probability':'1/2 each, recorded','branches':rows,
            'complete_total_Kraus_identity':'sum_kind(Ksuccess^dagger Ksuccess+Kfailure^dagger Kfailure)/2=I',
            'success_changes_fermion_parity':True,'reservoir_retains_complementary_parity_lepton_number_and_spin':True,
            'reservoir_gauge_charge':'0','controller_is_declared_finite_apparatus':True,
            'controller_is_idle_Hamiltonian':False,'conditioned_success_normalization_applied':False}


@lru_cache(maxsize=None)
def operator_column(key,kind='remove',site=1,spin=0,outcome='raw'):
    j,success,failure,_=local_instrument(kind,spin)
    op={'raw':j,'success':success,'failure':failure}[outcome]
    mapping={9:0,10:1,11:2,24:3,25:4,26:5}
    out=defaultdict(lambda:s.Integer(0))
    for (word,row,col),value in op.items():
        if key[3+site]!=col:continue
        operators=tuple((6*site+mapping[mode],create) for mode,create in word)
        step=sequence(key[0],operators)
        if step is not None:
            target=list(key);target[0]=step[0];target[3+site]=row
            out[tuple(target)]+=value*step[1]
    return clean(out)
