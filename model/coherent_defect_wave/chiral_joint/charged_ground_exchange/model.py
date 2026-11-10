"""Complete charged response blocks and separate-spin source identities."""
from functools import lru_cache
import sympy as s
from ..local_reference import physical,algebra
from ..local_reference.reduction import coordinate_kernel,quark_module,lepton_operators
from ..ground_composite_response.model import reduced,source,polynomial_map,validate_annihilator
from ..higgs_composites import model as car
from ...chiral_higgs.model import higgs_memory
from ...chiral_source.model import PAULI


@lru_cache(maxsize=1)
def separate_spin():
    parts,couplings=car.quantum_yukawa()
    h=car.clean({key:value.subs({c:1 for c in couplings.values()})
                 for key,value in car.add(*parts.values()).items()})
    h=car.add(h,car.monomial(matrix=higgs_memory()[2]))
    checks=0
    for sector in (algebra.LEPTONS,algebra.QUARKS):
        internal=s.diag(*(int(a in sector) for a in range(15)))
        for sigma in PAULI:
            current=car.one_body(s.kronecker_product(sigma/2,internal))
            assert car.commutator(current,h)=={}
            for mode in set(range(30))-set(sector):
                assert car.commutator(current,car.monomial(((mode,False),)))=={}
                checks+=1
    return {'independent_quark_and_lepton_spin_H_commutators':6,
            'opposite_sector_elementary_field_commutators':checks,
            'all_three_quark_spin_response_functions':'0 for every semigroup time',
            'source_identity_uses_all_30_CAR_modes_and_quantum_scalar_states':True}


@lru_cache(maxsize=None)
def weight_basis(number,two_spin,twice_weak,twice_hyper):
    qw,lw,sw=physical.labels()
    return tuple((i,mask,h) for i,(qz,qy) in enumerate(qw)
                 for mask,(n,lz,ly,spin) in enumerate(lw) if (n,spin)==(number,two_spin)
                 for h,(sz,sy) in enumerate(sw)
                 if qz+lz+sz==twice_weak and qy+ly+sy==twice_hyper)


@lru_cache(maxsize=None)
def target(family,two_spin):
    if family not in ('ec','L') or two_spin not in (0,2):
        raise ValueError('specified elementary removal response required')
    number,weak,hyper=(4,0,-2) if family=='ec' else (2,1,1)
    keys=weight_basis(number,two_spin,weak,hyper)
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    raised=weight_basis(number,two_spin,weak+2,hyper)
    weak_map=physical.tensor_map(keys,raised,
        quark=q['currents'][0]+s.I*q['currents'][1],
        lepton=l['currents'][0]+s.I*l['currents'][1],scalar=scalar[8]+s.I*scalar[9])
    basis,_=coordinate_kernel(weak_map)
    assert basis.cols==len(keys)-len(raised)
    return reduced(keys,basis)


@lru_cache(maxsize=None)
def transfer(family,spin):
    if spin not in (0,1):raise ValueError('native spin component required')
    sb=source(1);tb=target(family,2*spin)
    # The native lower weak component has a highest-weight conjugate image.
    mode=(9 if family=='ec' else 11)+15*spin
    raw=polynomial_map(sb['keys'],tb['keys'],car.monomial(((mode,False),)))
    image=raw*sb['basis']
    matrix=s.ImmutableSparseMatrix(tb['left']*image.extract(tb['coordinates'],range(sb['basis'].cols)))
    assert tb['basis']*matrix==image
    adjoint=s.ImmutableSparseMatrix(sb['gram'].inv()*matrix.H*tb['gram'])
    return {'H':tb['H'],'gram':tb['gram'],'map':matrix,'adjoint':adjoint,
            'full_weight_dimension':len(tb['keys']),'complete_charged_irrep_dimension':tb['basis'].cols,
            'native_removed_mode':mode}


@lru_cache(maxsize=None)
def moment_data(family):
    sb=source(1);x=s.Symbol('x')
    p=s.Poly(sb['H'].charpoly(x).as_expr(),x).sqf_part()
    blocks=[transfer(family,spin) for spin in range(2)]
    q=s.Poly(1,x)
    for block in blocks:
        q=s.lcm(q,s.Poly(block['H'].charpoly(x).as_expr(),x).sqf_part())
        validate_annihilator(block['H'],q)
    source_powers=[s.eye(sb['H'].rows)]
    for _ in range(p.degree()-1):source_powers.append(sb['H']*source_powers[-1])
    images=[block['map'] for block in blocks];moments=[]
    for _ in range(q.degree()):
        pullback=blocks[0]['adjoint']*images[0]-blocks[1]['adjoint']*images[1]
        moments.append(tuple(s.trace(power*pullback) for power in source_powers))
        images=[block['H']*image for block,image in zip(blocks,images,strict=True)]
    return p,q,tuple(moments),blocks
