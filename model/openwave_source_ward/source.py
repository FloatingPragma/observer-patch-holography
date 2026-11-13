"""Exact difference geometry, retaining the actual native lapse response."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sympy as sy
from model.openwave_universal_energy_geometry import source as native


def add(*fields):
    keys=set().union(*(set(v) for v in fields))
    return {p:z for p in keys if (z:=sum(v.get(p,F(0)) for v in fields))}


def scale(c,field):
    return {p:c*v for p,v in field.items() if c*v}


def shift(field,axis,amount):
    result={}
    for p,v in field.items():
        q=list(p);q[axis]+=amount;result[tuple(q)]=v
    return result


def derivative(field,axis):
    # E_i f(x)=f(x+e_i); dictionary shifts therefore have the opposite sign.
    return scale(F(1,2),add(shift(field,axis,-1),scale(-1,shift(field,axis,1))))


def second(field,axis):
    return add(shift(field,axis,1),shift(field,axis,-1),scale(-2,field))


def central_laplacian(field):
    return add(*(derivative(derivative(field,i),i) for i in range(3)))


def nearest_laplacian(field):
    return add(*(second(field,i) for i in range(3)))


def einstein_static(phi,gamma):
    """First metric variation h00=-2phi, hij=-2gamma phi delta_ij."""
    lap=central_laplacian(phi)
    out={(0,0):scale(2*gamma,lap)}
    for i,j in product(range(3),repeat=2):
        hess=derivative(derivative(phi,i),j)
        out[i+1,j+1]=scale(gamma-1,add(hess,scale(-int(i==j),lap)))
    for j in range(1,4):out[0,j]=out[j,0]={}
    return out


def divergence_static(tensor,j):
    return add(*(derivative(tensor[i+1,j],i) for i in range(3)))


def maximum(field):
    return max((abs(x) for x in field.values()),default=F(0))


def profile():
    packet=json.loads(Path(native.__file__).with_name('profile_fixture.json').read_text())
    n=packet['side'];den=packet['denominator'];coords=list(product(range(n),repeat=3))
    assert n==7 and den==2**40 and len(coords)==len(packet['integer_coordinates'])
    values={p:F(v,den) for p,v in zip(coords,packet['integer_coordinates'],strict=True)}
    assert all(0<=v<=1 for v in values.values()), 'native profile outside proved tanh error domain'
    center=(n//2,)*3;assert values[center]==1
    lap=scale(-1,nearest_laplacian(values))
    stiffness=sum(v*lap.get(p,F(0)) for p,v in values.items())
    assert stiffness==native.profile()['stiffness_exact']
    return values,center,stiffness


@lru_cache(None)
def symbolic_identities():
    """Independent full four-coordinate Ricci calculation, not a static ansatz."""
    p=sy.symbols('p0:4');eta=(-1,1,1,1)
    h=sy.zeros(4)
    for i in range(4):
        for j in range(i,4):h[i,j]=h[j,i]=sy.Symbol(f'h{i}{j}')
    box=sum(eta[i]*p[i]**2 for i in range(4))
    def curvature(matrix):
        trace=sum(eta[i]*matrix[i,i] for i in range(4))
        ricci=sy.Matrix(4,4,lambda a,b:sy.expand((sum(eta[c]*p[c]*(p[a]*matrix[c,b]+p[b]*matrix[c,a]) for c in range(4))-box*matrix[a,b]-p[a]*p[b]*trace)/2))
        scalar=sum(eta[i]*ricci[i,i] for i in range(4))
        return sy.Matrix(4,4,lambda a,b:sy.expand(ricci[a,b]-F(1,2)*int(a==b)*eta[a]*scalar))
    G=curvature(h)
    assert all(sy.expand(sum(eta[a]*p[a]*G[a,b] for a in range(4)))==0 for b in range(4))
    xi=sy.symbols('xi0:4');pure=sy.Matrix(4,4,lambda a,b:p[a]*xi[b]+p[b]*xi[a])
    assert curvature(pure)==sy.zeros(4)
    phi,gamma=sy.symbols('phi gamma')
    static=sy.diag(-2*phi,*([-2*gamma*phi]*3))
    sG=curvature(static).subs(p[0],0);lap=sum(x*x for x in p[1:])*phi
    assert sy.expand(sG[0,0]-2*gamma*lap)==0
    for i,j in product(range(1,4),repeat=2):
        assert sy.expand(sG[i,j]-(gamma-1)*(p[i]*p[j]*phi-int(i==j)*lap))==0
    return {'four_coordinate_bianchi_polynomial_zero':True,'pure_gauge_curvature_zero':True,
            'static_formula_from_full_ricci_variation':True,'commuting_derivatives_required':True}


def dyadic_interval(interval,bits=64):
    D=2**bits
    return F((interval.lo*D).__floor__(),D),F((interval.hi*D).__ceil__(),D)


@lru_cache(None)
def clock_slope(f):
    """d/dE log(<exp(FQ)>_E / <exp(FQ)>_0), at E=0."""
    sinhq,_=native.sinhcosh(native.QSIZE)
    sf,cf=native.sinhcosh(native.QSIZE*f)
    return dyadic_interval(-2*sinhq*sf/cf/native.OMEGA)


def run():
    values,center,S=profile();phi=scale(-1,values)
    # Zero extension is kept. No boundary/halo outcome is silently removed.
    nearest=nearest_laplacian(phi);central=central_laplacian(phi)
    fourth=scale(F(1,4),add(*(second(second(phi,i),i) for i in range(3))))
    assert add(central,scale(-1,nearest))==fourth
    controls=[]
    for gamma in (F(-1),F(0),F(1),F(2)):
        tensor=einstein_static(phi,gamma)
        assert all(not divergence_static(tensor,j) for j in range(4))
        anisotropic=max(maximum(tensor[i,j]) for i,j in product(range(1,4),repeat=2))
        controls.append({'gamma':str(gamma),'all_four_divergences_exact_zero':True,
                         'spatial_stress_residual_per_unit_potential':str(anisotropic),
                         'pressureless_static_tensor_condition':anisotropic==0})
    # An actual nonzero off-diagonal Hessian makes gamma=1 necessary here.
    mixed=derivative(derivative(phi,0),1);witness=max(mixed,key=lambda p:abs(mixed[p]))
    assert mixed[witness]!=0
    lap=scale(-1,nearest_laplacian(values))
    poisson={p:lap.get(p,F(0))/S-int(p==center) for p in values}
    sinhq,_=native.sinhcosh(native.QSIZE)
    strength=dyadic_interval(sinhq/native.QSIZE)
    clock=[{'site':list(p),'source_F':str(f),'log_lapse_source_derivative_interval':list(map(str,clock_slope(f)))} for p,f in sorted(values.items())]
    bound=dyadic_interval(2*sinhq*(native.QSIZE**3)/3/native.OMEGA)
    return {'schema':'oph.native_lapse_linearized_ward.v1','verdict':'PASS_NATIVE_LAPSE_AND_EXACT_LINEARIZED_BIANCHI',
      'symbolic':symbolic_identities(),'source_sites':len(values),'boundary_convention':'entire zero extension, including every difference halo',
      'source_stiffness_exact':str(S),'source_q_exact':str(native.QSIZE),
      'native_clock_response':{'formula':'-2 sinh(q) tanh(q F)/omega','type':'source-energy derivative at zero source strength of the actual weak-probe log lapse',
         'complete_source_reference':'all eight native source levels and both geometry outcomes from the parent source; derivative is linear in the native energy operator',
         'records':clock,'far_field_response_over_small_q_approximation_interval':list(map(str,strength)),
         'uniform_error_to_linear_F_profile_interval':list(map(str,bound)),
         'relative_profile_error_upper_exact':str(native.QSIZE**2/3)},
      'finite_geometry':{'metric_assignment':'h00=-2 phi; hij=-2 gamma phi delta_ij, as a declared first metric variation',
         'controls':controls,'gamma_one_is_forced_if_spatial_linearized_stress_vanishes':True,
         'nonzero_mixed_hessian_witness':{'site':list(witness),'value':str(mixed[witness])},
         'nearest_vs_central_stencil_identity':'Delta_c phi - Delta_nn phi = (a^2/4) sum_i (Delta_nn,i)^2 phi',
         'exact_native_profile_poisson_residual':str(maximum(poisson)),
         'all_halo_stencil_correction_max':str(maximum(fourth)),
         'all_halo_spatial_points':len(set().union(*(set(v) for v in einstein_static(phi,F(0)).values())))},
      'scope':{'native_clock_derivative_derived':True,'finite_linearized_bianchi_derived':True,
         'physical_metric_assignment_derived':False,'nonlinear_einstein_equation_derived':False,
         'physical_Newton_strength_matched':False,'no_new_empirical_fit':True}}
