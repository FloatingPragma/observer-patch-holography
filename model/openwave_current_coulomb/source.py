"""Exact native current, Gauss projection and an explicitly changed Hamiltonian."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.periodic_family.model import (
    Configuration,cells,cell_index,edges,edge_endpoints,kinetic_terms,local_yukawa)
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current
from model.openwave_reference_binding import strings


@lru_cache(maxsize=4)
def cochains(side=2):
    Configuration(side=side)
    volume=side**3
    d=s.zeros(3*volume,volume)
    for edge in range(3*volume):
        left,right=edge_endpoints(side,edge)
        d[edge,cell_index(side,left)]=-1
        d[edge,cell_index(side,right)]=1
    d=s.ImmutableSparseMatrix(d)
    lap=d.T*d;zero=s.ones(volume)/volume
    if side in (2,4):
        # Exact finite Fourier characters remove a dense symbolic inverse.
        # The independent incidence identity below verifies the result.
        cosines=(1,-1) if side==2 else (1,0,-1,0)
        kernel={}
        for displacement in cells(side):
            value=s.Integer(0)
            for momentum in cells(side):
                eigenvalue=2*sum(1-cosines[k] for k in momentum)
                if eigenvalue:
                    phase=sum(k*x for k,x in zip(momentum,displacement,strict=True))%side
                    value+=s.Rational(cosines[phase],volume*eigenvalue)
            kernel[displacement]=value
        green=s.Matrix(volume,volume,lambda i,j:kernel[tuple((x-y)%side
            for x,y in zip(cells(side)[i],cells(side)[j],strict=True))])
    else:
        green=(lap+zero).inv()-zero
    assert lap*green==s.eye(volume)-zero and green*s.ones(volume,1)==s.zeros(volume,1)
    return s.ImmutableSparseMatrix(d),s.ImmutableSparseMatrix(green)


def native_density(key,reference_mask=0,side=2):
    config=Configuration(side=side)
    mask,links,scalars=pw.validate_key(key,config)
    if type(reference_mask) is not int or not 0<=reference_mask<1<<30:
        raise ValueError('complete thirty-mode reference occupation required')
    weights=representation()[2][11].diagonal()
    scalar=higgs_memory()[1][11].diagonal()
    density=[sum(weights[j%15] for j in range(30) if mask>>(30*x+j)&1)+scalar[scalars[x]]
             for x in range(side**3)]
    density[0]+=sum(weights[j%15] for j in range(30) if reference_mask>>j&1)
    return s.Matrix(density)


def native_flux(link_labels,side=2):
    if len(link_labels)!=3*side**3 or any(type(x) is not int or x not in range(119) for x in link_labels):
        raise ValueError('every full PW119 link label is required')
    values=endpoint_current(11,0).diagonal()
    return s.Matrix([values[i] for i in link_labels])


def gauss_minimum(density,side=2):
    d,g=cochains(side)
    density=s.Matrix(density)
    if density.shape!=(side**3,1) or any(x.is_Rational is not True for x in density):
        raise ValueError('one exact native current value at every cell required')
    if sum(density)!=0:raise ValueError('periodic Gauss readout requires zero total source plus reference current')
    field=d*g*density
    energy=(density.T*g*density)[0]/2
    assert d.T*field==density
    assert (field.T*field)[0]/2==energy
    return field,s.cancel(energy)


def coulomb_diagonal(key,reference_mask=0,side=2,kappa=s.Integer(1)):
    kappa=s.sympify(kappa)
    if kappa.is_Rational is not True or kappa<=0:raise ValueError('positive exact field energy coefficient required')
    rho=native_density(key,reference_mask,side)
    # This defines the full-carrier operator. The zero mode is projected by G;
    # the physical domain has sum rho=0, and its field readout rejects violations.
    return s.cancel(kappa*(rho.T*cochains(side)[1]*rho)[0]/2)


def changed_column(key,reference_mask=0,side=2,kappa=s.Integer(1),coupling=s.Rational(1,2**34)):
    """H'=H_onsite+lambda T+4N_ref+kappa rho^T G rho/2, no link penalty."""
    coupling=s.sympify(coupling)
    if coupling.is_Rational is not True or coupling<0:raise ValueError('nonnegative exact hopping coefficient required')
    config=Configuration(side=side,link_penalty=0)
    result=defaultdict(lambda:s.Integer(0),pw.column(key,config))
    for target,coefficient in pw.hopping_column(key,kinetic_terms(config)).items():
        result[target]+=(coupling-1)*coefficient
    result[key]+=4*reference_mask.bit_count()+coulomb_diagonal(key,reference_mask,side,kappa)
    return {target:s.expand(coefficient) for target,coefficient in result.items() if coefficient!=0}


def path(side,vertices):
    route=[]
    for left,right in zip(vertices,vertices[1:]):
        match=[]
        for edge in range(3*side**3):
            x,y=edge_endpoints(side,edge)
            if (x,y)==(left,right):match.append((edge,1))
            elif (x,y)==(right,left):match.append((edge,0))
        if not match:raise ValueError('nearest-neighbor native route required')
        edge,forward=min(match)
        route.append((edge,cell_index(side,left),cell_index(side,right),forward))
    if len({step[0] for step in route})!=len(route):raise ValueError('simple stored-edge path required')
    return tuple(route)


def physical_path(side,vertices):
    route=path(side,vertices)
    state=strings.string_state('ec',route)
    for node in {0,*[x[1] for x in route],*[x[2] for x in route]}:
        for current in range(12):assert not strings.gauss_image('ec',route,state,node,current)
    assert len(state)==1
    (_,active,_),amplitude=next(iter(state.items()));assert amplitude==1
    links=[0]*(3*side**3)
    for row,label in zip(route,active,strict=True):links[row[0]]=label
    rho=s.zeros(side**3,1);rho[0]+=1;rho[cell_index(side,vertices[-1])]-=1
    b=native_flux(links,side);d,g=cochains(side)
    assert d.T*b==rho
    minimum,energy=gauss_minimum(rho,side)
    transverse=b-minimum
    assert d.T*transverse==s.zeros(side**3,1)
    assert (minimum.T*transverse)[0]==0
    assert (b.T*b)[0]==(minimum.T*minimum)[0]+(transverse.T*transverse)[0]
    return {'path':[list(x) for x in route],'full_native_link_labels':links,
        'length':len(route),'density':[str(x) for x in rho],
        'native_oriented_flux':[str(x) for x in b],
        'minimum_Gauss_field':[str(x) for x in minimum],
        'unit_kappa_Coulomb_energy':str(energy),
        'native_flux_squared_norm':str((b.T*b)[0]),
        'Gauss_minimum_squared_norm':str((minimum.T*minimum)[0]),
        'transverse_native_flux_squared_norm':str((transverse.T*transverse)[0]),
        'all_twelve_actual_endpoint_Gauss_constraints':True,
        'all_link_states_retained_in_Hamiltonian':True}
