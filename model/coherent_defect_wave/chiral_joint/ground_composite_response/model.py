"""Actual J/J-dagger maps on complete physical target modules."""
from collections import defaultdict
from functools import lru_cache
import sympy as s

from ..local_reference import physical
from ..local_reference.algebra import LEPTONS,particle_hole
from ..local_reference.reduction import coordinate_kernel,lepton_operators,quark_module
from ..higgs_composites import model as car
from ..neutral_composite.model import local_instrument
from ..model import sequence
from ...chiral_higgs.model import higgs_memory


def reduced(keys,basis):
    coordinates=tuple(basis.T.rref()[1])
    left=basis.extract(coordinates,range(basis.cols)).inv()
    metric=physical.metric_on(keys)
    gram=s.ImmutableSparseMatrix(basis.H*metric*basis)
    image=physical.hamiltonian_on(keys)*basis
    h=s.ImmutableSparseMatrix(left*image.extract(coordinates,range(basis.cols)))
    assert basis*h==image and gram*h==h.H*gram
    return {'keys':keys,'basis':basis,'coordinates':coordinates,'left':left,'gram':gram,'H':h}


@lru_cache(maxsize=None)
def source(two_sz):
    if two_sz not in (-1,1):raise ValueError('ground doublet magnetic component required')
    keys,basis,*_=physical.physical_block(3,1)
    if two_sz==-1:
        target=physical.weight_basis(3,-1);spin=lepton_operators()['spin']
        lowering=physical.tensor_map(keys,target,lepton=spin[0]-s.I*spin[1])
        basis=lowering*basis;keys=target
    result=reduced(keys,basis)
    assert result['basis'].cols==16
    return result


@lru_cache(maxsize=None)
def target(number,two_sz):
    if number not in (2,4) or two_sz not in (-2,0,2):
        raise ValueError('complete even target charge and magnetic component required')
    keys=physical.weight_basis(number,two_sz)
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    weak=physical.tensor_map(keys,physical.weight_basis(number,two_sz,2),
        quark=q['currents'][0]+s.I*q['currents'][1],
        lepton=l['currents'][0]+s.I*l['currents'][1],scalar=scalar[8]+s.I*scalar[9])
    basis,_=coordinate_kernel(weak)
    assert basis.cols==len(keys)-len(physical.weight_basis(number,two_sz,2))
    assert basis.cols==(16 if two_sz==0 else 5)
    return reduced(keys,basis)


def polynomial_map(source_keys,target_keys,operator):
    """Quarks precede leptons in the tensor convention; Nq'=12 has even parity."""
    operator=particle_hole(operator)
    mapping={mode:index for index,mode in enumerate(LEPTONS)}
    positions={key:i for i,key in enumerate(target_keys)}
    entries=defaultdict(lambda:s.Integer(0))
    for column,(q,mask,h) in enumerate(source_keys):
        for (word,row,col),coefficient in operator.items():
            if col!=h:continue
            if any(mode not in mapping for mode,_ in word):
                raise ValueError('operator leaves the retained quark module')
            step=sequence(mask,tuple((mapping[mode],create) for mode,create in word))
            if step is None:continue
            destination=(q,step[0],row)
            if destination not in positions:raise ValueError('omitted full composite image')
            entries[positions[destination],column]+=coefficient*step[1]
    return s.ImmutableSparseMatrix(len(target_keys),len(source_keys),entries)


@lru_cache(maxsize=None)
def complete_map(kind,spin):
    if kind not in ('remove','add') or spin not in (0,1):raise ValueError('specified native probe required')
    j,success,failure,projectors=local_instrument(kind,spin)
    assert success==car.scale(j,s.Rational(1,2))
    number=2 if kind=='remove' else 4
    sources=[source(m) for m in (-1,1)]
    targets=[target(number,m) for m in (-2,0,2)]
    h=s.ImmutableSparseMatrix(s.diag(*(block['H'] for block in sources)))
    g=s.ImmutableSparseMatrix(s.diag(*(block['gram'] for block in sources)))
    ht=s.ImmutableSparseMatrix(s.diag(*(block['H'] for block in targets)))
    gt=s.ImmutableSparseMatrix(s.diag(*(block['gram'] for block in targets)))
    matrix=s.zeros(26,32);row_offsets={-2:0,0:5,2:21}
    for source_index,m in enumerate((-1,1)):
        shift=(-1 if spin==0 else 1)*(1 if kind=='remove' else -1)
        target_m=m+shift;sb=sources[source_index];tb=targets[(target_m+2)//2]
        raw=polynomial_map(sb['keys'],tb['keys'],j)
        image=raw*sb['basis']
        projected=tb['left']*image.extract(tb['coordinates'],range(sb['basis'].cols))
        assert tb['basis']*projected==image
        start=row_offsets[target_m]
        matrix[start:start+tb['basis'].cols,16*source_index:16*(source_index+1)]=projected
        # The adjoint is evaluated from the actual full CAR/scalar polynomial.
        reverse=polynomial_map(tb['keys'],sb['keys'],car.adjoint(j))
        reverse_image=reverse*tb['basis']
        # This image can contain the excluded spin-3/2 source; only its
        # orthogonal spin-1/2 projection enters this adjoint relation.
        left_form=sb['basis'].H*physical.metric_on(sb['keys'])*reverse_image
        assert left_form==projected.H*tb['gram']
    matrix=s.ImmutableSparseMatrix(matrix)
    adjoint=s.ImmutableSparseMatrix(g.inv()*matrix.H*gt)
    assert g*adjoint==matrix.H*gt
    assert g*h==h.H*g and gt*ht==ht.H*gt
    return {'source_H':h,'source_gram':g,'target_H':ht,'target_gram':gt,
            'J':matrix,'J_adjoint':adjoint,'target_raw_L':number-2,
            'source_dimension':32,'target_dimension':26,
            'all_composite_images_and_metric_adjoint_checked':True}


@lru_cache(maxsize=None)
def moments(kind,spin):
    data=complete_map(kind,spin);x=s.Symbol('x')
    p=s.Poly(data['source_H'].charpoly(x).as_expr(),x).sqf_part()
    q=s.Poly(data['target_H'].charpoly(x).as_expr(),x).sqf_part()
    assert p.degree()==q.degree()==9
    powers=[s.eye(32)]
    for _ in range(p.degree()-1):powers.append(data['source_H']*powers[-1])
    image=data['J'];rows=[]
    for _ in range(q.degree()):
        pullback=data['J_adjoint']*image
        rows.append([s.trace(power*pullback) for power in powers])
        image=data['target_H']*image
    # The squarefree polynomials annihilate the entire respective blocks.
    for polynomial,matrix in ((p,data['source_H']),(q,data['target_H'])):
        validate_annihilator(matrix,polynomial)
    return data,p,q,tuple(tuple(row) for row in rows)


def validate_annihilator(matrix,polynomial):
    residual=s.zeros(matrix.rows)
    for coefficient in polynomial.all_coeffs():residual=matrix*residual+coefficient*s.eye(matrix.rows)
    if residual!=s.zeros(matrix.rows):raise ValueError('polynomial omits part of the full physical spectrum')
