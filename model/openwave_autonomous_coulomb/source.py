"""Literal PW119 coordinate and complete invariant pair/loop sector."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import sympy as s
import numpy as np
from model.openwave_current_coulomb.source import native_density,coulomb_diagonal,cochains
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms,local_yukawa
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import coordinate
from model.coherent_defect_wave.chiral_joint.model import sequence
from model.coherent_defect_wave.physical_attachment.chiral import representation

SIDE=4
REFERENCE_MASK=(1<<30)-1-(1<<9)
CONFIG=Configuration(side=SIDE,link_penalty=0,y_up=0,y_down=0,y_lepton=0)


@lru_cache(maxsize=1)
def ingredients():
    U=coordinate(6,0,0);labels=(0,46,105);flux=pw.endpoint_current(11,0)
    for r,c in U.todok():assert (r in labels)==(c in labels)
    assert U.extract(labels,labels)==s.Matrix([[0,0,1],[1,0,0],[0,0,0]])
    for current in range(11):
        assert pw.endpoint_current(current,0)[:,labels]==s.zeros(119,3)
        assert pw.endpoint_current(current,1)[:,labels]==s.zeros(119,3)
    terms=tuple(t for t in kinetic_terms(CONFIG)if t.axis==2 and t.displacement==1 and t.left_cell in range(4))
    assert len(terms)==4
    weights=representation()[2][11].diagonal()
    reference_y=sum(weights[j%15]for j in range(30)if REFERENCE_MASK>>j&1)
    assert reference_y==-1 and weights[9]==1
    return U,labels,tuple(flux[i,i]for i in labels),terms


def column(key,kappa=s.Integer(1),phase_edge=None,theta=0,reference_mask=REFERENCE_MASK):
    mask,links,scalars=pw.validate_key(key,CONFIG);U,labels,flux,terms=ingredients()
    out=defaultdict(lambda:s.Integer(0))
    out[key]=2*mask.bit_count()+4*reference_mask.bit_count()+sum(v!=0 for v in scalars)+coulomb_diagonal(key,reference_mask,4,kappa)
    # This is the literal native ec matrix coordinate on every PW119 state.
    # Other species, axial paths and Yukawa vertices have declared zero weights.
    for term in terms:
        edge=term.path[0];phase=s.exp(s.I*theta)if edge==phase_edge else 1
        for adjoint in (False,True):
            matrix=U.H if adjoint else U
            for (row,col),factor in matrix.todok().items():
                if col!=links[edge]:continue
                for spin in range(2):
                    p=30*term.left_cell+15*spin+9;q=30*term.right_cell+15*spin+9
                    coefficient=term.coefficient*(1 if spin==0 else -1)*phase
                    if adjoint:p,q=q,p;coefficient=s.conjugate(coefficient)
                    target=sequence(mask,((p,True),(q,False)))
                    if target is not None:
                        other=list(links);other[edge]=row
                        out[target[0],tuple(other),scalars]+=coefficient*factor*target[1]
    return {key:s.simplify(v)for key,v in out.items()if v!=0}


@lru_cache(maxsize=1)
def basis():
    U,labels,flux,terms=ingredients();out=[]
    for spin,y,indices in product(range(2),range(4),product(range(3),repeat=4)):
        f=[flux[j]for j in indices]
        rho=[int(x==y)-int(x==0)for x in range(4)]
        if any(rho[x]+f[x]-f[(x-1)%4]!=0 for x in range(4)):continue
        links=[0]*192
        for x,j in enumerate(indices):links[3*x+2]=labels[j]
        key=(1<<(30*y+15*spin+9),tuple(links),(0,)*64)
        out.append((key,spin,y,tuple(f)))
    assert len(out)==18
    return tuple(out)


@lru_cache(maxsize=4)
def matrices(kappa=s.Integer(1)):
    rows=basis();where={key:i for i,(key,_,_,_)in enumerate(rows)};h=s.zeros(len(rows))
    for j,(key,spin,y,flux)in enumerate(rows):
        for target,value in column(key,kappa).items():
            assert target in where,'invariant physical source sector violated'
            h[where[target],j]+=value
    assert h==h.H
    v=s.diag(*[coulomb_diagonal(key,REFERENCE_MASK,4,kappa)for key,_,_,_ in rows])
    return h,v


@lru_cache(maxsize=1)
def native_certificate():
    rows=basis();checks=0
    for key,spin,y,flux in rows:
        rho=native_density(key,REFERENCE_MASK,4)
        assert sum(rho)==0 and rho[y]==(0 if y==0 else 1)
        for cell,current in product(range(4),range(12)):
            action=pw.gauss_column(key,cell,current,CONFIG)
            if cell==0 and current==11:
                action[key]=action.get(key,0)-1
                action={k:v for k,v in action.items()if v!=0}
            assert not action;checks+=1
    # The unchanged interacting onsite source does not preserve this sector.
    leakage=local_yukawa(1<<9,0)
    assert leakage and any(mask!=(1<<9)or scalar!=0 for mask,scalar in leakage)
    U,labels,flux,terms=ingredients()
    h,v=matrices();nn=len(rows)
    pin=s.diag(*[int(y==1)for _,_,y,_ in rows]);comm=h*pin-pin*h
    assert comm!=s.zeros(nn)
    return {'full_carrier':'2^(30*64) *5^64 *119^192 *2^30_reference',
        'native_coordinate_dimension':119,'literal_ec_coordinate_invariant_labels':list(labels),
        'actual_endpoint_hypercharges':list(map(str,flux)),'active_ec_particle_hypercharge':1,
        'actual_29_particle_reference_hypercharge':-1,'reference_number_energy':116,
        'complete_invariant_pair_loop_spin_rows':nn,'all_active_site_native_Gauss_checks':checks,
        'other_site_Gauss_constraints':'all incident states and matter/scalar factors are native singlets',
        'complete_loop_patterns':len({f for _,_,_,f in rows}),
        'native_nearest_axial_coefficients':[str(t.coefficient)for t in terms],
        'position_pin_commutator_nonzero':True,'unchanged_unit_Yukawa_leakage_terms':len(leakage),
        'source_changes':'only native ec nearest z-line kinetic terms are enabled; all other native hopping coefficients and Yukawa vertices are zero; positive number mass2 and quadratic native-current Green energy replace the interacting local/link energy',
        'transfer':'native PW119 coordinates, currents, CAR signs, scalar/reference factors and full Gauss domain are literal; interacting ground, clocks and Yukawa dynamics are not transferred'}
