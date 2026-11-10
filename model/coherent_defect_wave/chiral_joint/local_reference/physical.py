"""Every weak/hypercharge singlet and spin multiplet of the maximal module."""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from .reduction import quark_module,lepton_operators,coordinate_kernel,restrict
from ...chiral_higgs.model import higgs_memory


def columns(matrix):
    result=defaultdict(list)
    for (i,j),value in matrix.todok().items():result[j].append((i,value))
    return result


@lru_cache(maxsize=1)
def labels():
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    lepton=tuple((mask.bit_count(),int(2*l['currents'][2][mask,mask]),
                  int(2*l['currents'][3][mask,mask]),int(2*l['spin'][2][mask,mask])) for mask in range(64))
    scalar=tuple((int(2*scalar[10][i,i]),int(2*scalar[11][i,i])) for i in range(5))
    return q['weights'],lepton,scalar


@lru_cache(maxsize=None)
def weight_basis(number,two_spin,twice_weak=0):
    qw,lw,sw=labels()
    return tuple((i,mask,h) for i,(qz,qy) in enumerate(qw)
                 for mask,(n,lz,ly,spin) in enumerate(lw) if (n,spin)==(number,two_spin)
                 for h,(sz,sy) in enumerate(sw) if qz+lz+sz==twice_weak and qy+ly+sy==0)


def tensor_map(source,target,quark=None,lepton=None,scalar=None):
    positions={key:i for i,key in enumerate(target)};entries=defaultdict(lambda:s.Integer(0))
    factors=tuple(None if matrix is None else columns(matrix) for matrix in (quark,lepton,scalar))
    for j,key in enumerate(source):
        for axis,factor in enumerate(factors):
            if factor is None:continue
            for destination,value in factor.get(key[axis],()):
                changed=list(key);changed[axis]=destination;changed=tuple(changed)
                assert changed in positions,('omitted current image',key,changed)
                entries[positions[changed],j]+=value
    return s.ImmutableSparseMatrix(len(target),len(source),entries)


def metric_on(keys):
    qg=columns(quark_module()['metric']);positions={key:i for i,key in enumerate(keys)};entries={}
    for j,(q,mask,h) in enumerate(keys):
        for destination,value in qg.get(q,()):
            target=(destination,mask,h)
            assert target in positions
            entries[positions[target],j]=value
    return s.ImmutableSparseMatrix(len(keys),len(keys),entries)


def hamiltonian_on(keys):
    q=quark_module();l=lepton_operators()
    qc=tuple(columns(m) for m in q['transitions']);lc=tuple(columns(m) for m in l['transitions'])
    positions={key:i for i,key in enumerate(keys)};entries=defaultdict(lambda:s.Integer(0))
    for j,(a,mask,h) in enumerate(keys):
        if h:entries[j,j]+=1
        branches=tuple((r,r+1) for r in range(4)) if h==0 else (((h+1)%4,0),)
        for r,target_scalar in branches:
            for destination,value in qc[r].get(a,()):
                target=(destination,mask,target_scalar)
                assert target in positions,('omitted quark/scalar Hamiltonian image',target)
                entries[positions[target],j]+=value
            for destination,value in lc[r].get(mask,()):
                target=(a,destination,target_scalar)
                assert target in positions,('omitted lepton/scalar Hamiltonian image',target)
                entries[positions[target],j]+=value
    return s.ImmutableSparseMatrix(len(keys),len(keys),entries)


@lru_cache(maxsize=None)
def physical_block(number,two_spin):
    keys=weight_basis(number,two_spin)
    if not keys:return keys,s.zeros(0),s.zeros(0),s.zeros(0),s.zeros(0)
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    weak=tensor_map(keys,weight_basis(number,two_spin,2),
        quark=q['currents'][0]+s.I*q['currents'][1],
        lepton=l['currents'][0]+s.I*l['currents'][1],scalar=scalar[8]+s.I*scalar[9])
    spin=tensor_map(keys,weight_basis(number,two_spin+2),lepton=l['spin'][0]+s.I*l['spin'][1])
    basis,coordinates=coordinate_kernel(weak.col_join(spin))
    expected=(len(keys)-len(weight_basis(number,two_spin,2))
              -len(weight_basis(number,two_spin+2))+len(weight_basis(number,two_spin+2,2)))
    assert basis.cols==expected>=0
    metric=metric_on(keys);h=hamiltonian_on(keys)
    assert metric*h==h.H*metric
    generator=restrict(h,basis,coordinates)
    gram=s.ImmutableSparseMatrix(basis.H*metric*basis)
    energy=s.ImmutableSparseMatrix(gram*generator)
    assert energy==energy.H
    return keys,basis,gram,energy,generator


def scalar_split(keys,basis,generator):
    """Check the full scalar decomposition before reducing the squared transfer."""
    blank=[];excited=[]
    for j in range(basis.cols):
        support={int(keys[i][2]!=0) for i in range(basis.rows) if basis[i,j]}
        if len(support)!=1:raise ValueError('physical basis mixes scalar excitation sectors')
        (excited if support.pop() else blank).append(j)
    if generator.extract(blank,blank)!=s.zeros(len(blank)):
        raise ValueError('nonzero scalar-blank diagonal block')
    if generator.extract(excited,excited)!=s.eye(len(excited)):
        raise ValueError('scalar-excited diagonal block is not the unit penalty')
    squared=s.ImmutableSparseMatrix(generator.extract(blank,excited)*generator.extract(excited,blank))
    if squared.det()==0 and blank:raise ValueError('squared transfer has a dark scalar-blank state')
    if len(excited)<len(blank):raise ValueError('invalid scalar transfer dimensions')
    energy=s.Symbol('E');x=s.Symbol('x')
    polynomial=s.Poly(squared.charpoly(x).as_expr(),x)
    expected=(energy-1)**(len(excited)-len(blank))*polynomial.as_expr().subs(x,energy*(energy-1))
    if s.expand(generator.charpoly(energy).as_expr()-expected)!=0:
        raise ValueError('complete scalar characteristic polynomial mismatch')
    return blank,excited,squared,polynomial


@lru_cache(maxsize=1)
def certificate():
    records=[];energy=s.Symbol('E');total=0
    for number in range(7):
        for twice_spin in range(4):
            keys,basis,gram,form,generator=physical_block(number,twice_spin)
            if not basis.cols:continue
            polynomial=s.Poly(generator.charpoly(energy).as_expr(),energy)
            blank,excited,squared,squared_polynomial=scalar_split(keys,basis,generator)
            total+=(twice_spin+1)*basis.cols
            records.append({'transformed_lepton_number':number,'raw_L':number-2,
                'twice_spin':twice_spin,'spin_multiplicity':twice_spin+1,
                'complete_highest_spin_dimension':basis.cols,'weight_zero_dimension':len(keys),
                'exact_complete_H_and_Gauss_intertwiner':True,
                'scalar_blank_dimension':len(blank),'scalar_excited_dimension':len(excited),
                'dark_energy_one_dimension':len(excited)-len(blank),
                'squared_transfer_polynomial':str(squared_polynomial.as_expr()),
                'squared_transfer_factorization':str(s.factor(squared_polynomial.as_expr())),
                'characteristic_polynomial':str(polynomial.as_expr()),
                'characteristic_factorization':str(s.factor(polynomial.as_expr())),
                'metric':[[str(value) for value in row] for row in gram.tolist()],
                'Hamiltonian_generator':[[str(value) for value in row] for row in generator.tolist()]})
    direct=sum(len(weight_basis(n,spin))-len(weight_basis(n,spin,2))
               for n in range(7) for spin in range(-3,4))
    assert total==direct
    return {'maximal_quark_dimension':140,'retained_lepton_Fock_dimension':64,
            'scalar_dimension':5,'full_maximal_module_dimension':140*64*5,
            'full_maximal_module_physical_dimension':total,'complete_spin_blocks':records,
            'all_lepton_occupations_and_all_scalar_states_retained':True,
            'exact_independent_weak_and_spin_character_count':True}
