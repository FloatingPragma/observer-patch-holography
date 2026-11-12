"""Native energy, nonlinear geometry minima and complete finite certificates."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.sparse import diags,kron,eye
from scipy.sparse.linalg import spsolve
from model.openwave_native_force.force import source_matrices,laplacian,green_solve
from model.finite_area_response.response import Interval,exponential,spectral_enclosures

KAPPA=F(10000);WINDOW=F(1,1000);BITS=40

@lru_cache(maxsize=1)
def native():
    _,_,matrix,_=source_matrices()
    levels=spectral_enclosures(matrix)
    gaps=(Interval.point(0),)+tuple(level-levels[0] for level in levels[1:])
    return matrix,levels,gaps


def exp_interval(value):
    if isinstance(value,Interval):return Interval(exponential(value.lo).lo,exponential(value.hi).hi)
    return exponential(F(value))


def energy_interval(level):
    if type(level)is not int or not 1<=level<=7:raise ValueError('an excited native level from1to7 is required')
    return native()[2][level]


def prepare(n=31,sources=None,kappa=KAPPA,bits=BITS,window=WINDOW):
    """Prepare a dyadic variational witness; no geometry-update implementation."""
    if sources is None:sources=(((n//2,)*3,1),)
    kappa=F(kappa);window=F(window)
    if type(n)is not int or n<5 or n%2!=1 or kappa<=0 or not 0<window<=F(1,1000) or not 20<=bits<=40:
        raise ValueError('odd cube, positive stiffness and supported finite geometry required')
    if len({tuple(site)for site,_ in sources})!=len(sources):raise ValueError('source sites must be distinct')
    rho=np.zeros((n,n,n));intervals={}
    for site,level in sources:
        site=tuple(site)
        if len(site)!=3 or any(type(x)is not int or not 0<=x<n for x in site):raise ValueError('source outside the finite cube')
        interval=energy_interval(level);intervals[site]=interval;rho[site]=float((interval.lo+interval.hi)/2)
    # On the invariant negative window, exp is1-Lipschitz. The positive
    # Green kernel bounds this fixed-point contraction by max(G rho)/kappa.
    linear=green_solve(rho,float(kappa));contraction=float(np.max(abs(linear)))
    if contraction>=float(window):raise ValueError('source exceeds the proved geometry window')
    denominator=2**bits
    linear_integer=np.rint(linear*denominator).astype(np.int64)
    linear_lap=laplacian(linear_integer);empty_linear=linear_lap.copy()
    for site in intervals:empty_linear[site]=0
    linear_residual=max(abs(kappa*F(int(empty_linear.min()),denominator)),abs(kappa*F(int(empty_linear.max()),denominator)))
    for site,e in intervals.items():
        r=Interval.point(kappa*F(int(linear_lap[site]),denominator))+e
        linear_residual=max(linear_residual,abs(r.lo),abs(r.hi))
    green_infinity=F((n+1)**2,8)
    contraction_upper=F(int(np.max(abs(linear_integer))),denominator)+green_infinity*linear_residual/kappa
    assert contraction_upper<window
    q=linear.copy()
    for _ in range(8):q=green_solve(rho*np.exp(q),float(kappa))
    denominator=2**bits;integer=np.rint(q*denominator).astype(np.int64)
    if np.max(abs(integer))>=float(window)*denominator:raise ValueError('dyadic witness outside the finite carrier')
    q=integer/denominator;lq=laplacian(integer)
    qlq=sum(int(x)*int(y)for x,y in zip(integer.ravel(),lq.ravel(),strict=True))
    energy=Interval.point(kappa*F(qlq,2*denominator**2))
    # The raw Laplacian at occupied sites includes the source cancellation;
    # mask it before bounding every source-free residual with exact integers.
    empty=lq.copy()
    for site in intervals:empty[site]=0
    maxres=max(abs(kappa*F(int(empty.min()),denominator)),abs(kappa*F(int(empty.max()),denominator)))
    source_rows=[]
    for site,e in intervals.items():
        value=F(int(integer[site]),denominator);factor=exp_interval(value)
        energy+=e*(factor-1)
        residual=Interval.point(kappa*F(int(lq[site]),denominator))+e*factor
        maxres=max(maxres,abs(residual.lo),abs(residual.hi))
        source_rows.append({'site':list(site),'native_gap_interval':e.encoded(),'geometry':str(value),'lapse_interval':factor.encoded(),'residual_interval':residual.encoded()})
    spectral_lower=F(12,(n+1)**2);strong=kappa*spectral_lower
    norm2=F(n**3)*maxres**2;error=norm2/(2*strong)
    # Strong convexity bounds the unrestricted minimum below; the dyadic
    # point is an admissible finite-grid upper witness.
    interval=[energy.lo-error,energy.hi]
    return {'n':n,'kappa':str(kappa),'dyadic_bits':bits,'window':str(window),
       'integer_geometry':integer,'source_rows':source_rows,
       'minimum_energy_interval':list(map(str,interval)),
       'witness_energy_interval':energy.encoded(),'residual_infinity_upper':str(maxres),
       'residual_squared_norm_upper':str(norm2),'finite_minimum_bracket_width':str(interval[1]-interval[0]),
       'continuum_relaxation_error_upper':str(error),
       'strong_convexity_lower':str(strong),'geometry_max_absolute':str(F(int(np.max(abs(integer))),denominator)),
       'fixed_point_contraction_diagnostic':contraction,'fixed_point_contraction_rational_upper':str(contraction_upper),
       'minimum_exactly_enumerated':False}


def serializable(preparation):
    import hashlib
    row=dict(preparation);integer=row.pop('integer_geometry')
    row['integer_geometry_sha256']=hashlib.sha256(integer.astype('<i8').tobytes()).hexdigest()
    return row


def difference(pair,left,right):
    p,a,b=(list(map(F,x['minimum_energy_interval']))for x in(pair,left,right))
    return [p[0]-a[1]-b[1],p[1]-a[0]-b[0]]


def run():
    matrix,levels,gaps=native();n=31;center=(15,15,15);one=prepare(n,((center,1),))
    rows=[];cross={}
    for distance in(4,5,6,7,8):
        site=(15+distance,15,15);other=prepare(n,((site,1),));pair=prepare(n,((center,1),(site,1)))
        value=difference(pair,one,other);assert value[1]<0
        cross[distance]=value
        rows.append({'separation':distance,'pair':serializable(pair),'second_source':serializable(other),'cross_energy_interval':list(map(str,value))})
    force=[(cross[7][0]-cross[5][1])/2,(cross[7][1]-cross[5][0])/2];assert force[0]>0
    d=float((gaps[1].lo+gaps[1].hi)/2)
    expected=d*d/(4*np.pi*float(KAPPA)*6**2)
    force_ratio=[float(x)/expected for x in force]
    # Small-grid Hessian demonstrates the exact inverse-positivity theorem.
    small=prepare(5,(((2,2,2),1),((3,2,2),7)))
    q=small['integer_geometry'].ravel()/2**BITS;axis=diags([-np.ones(4),2*np.ones(5),-np.ones(4)],[-1,0,1]);i=eye(5)
    l=kron(kron(axis,i),i)+kron(kron(i,axis),i)+kron(kron(i,i),axis)
    weights=np.zeros(125)
    for site,level in(((2,2,2),1),((3,2,2),7)):weights[np.ravel_multi_index(site,(5,5,5))]=float((gaps[level].lo+gaps[level].hi)/2)
    h=float(KAPPA)*l+diags(weights*np.exp(q));inverse=np.linalg.inv(h.toarray());assert inverse.min()>0
    a=np.ravel_multi_index((2,2,2),(5,5,5));b=np.ravel_multi_index((3,2,2),(5,5,5))
    mixed=-np.exp(q[a]+q[b])*inverse[a,b];assert mixed<0
    scale_rows=[]
    for stiffness in(KAPPA,2*KAPPA):
        p=prepare(15,(((7,7,7),1),),stiffness)
        scale_rows.append({'kappa':str(stiffness),'central_geometry':str(F(int(p['integer_geometry'][7,7,7]),2**BITS))})
    return {'native_complete_dimension':8,'native_reference':'actual lowest eigenvalue of M0; K=M0-E0 I',
       'native_M0_rational':[[str(x)for x in matrix.row(i)]for i in range(8)],
       'native_spectrum_intervals':[v.encoded()for v in levels],'native_ground_relative_gap_intervals':[v.encoded()for v in gaps],
       'declared_action':'kappa/2 q^T L q + sum_x (exp(q_x)-1) E_x',
       'controlled_matter_Hamiltonian':'sum_x exp(q_x) K_x on all eight native levels',
       'universal_operator_identity':'dH_x/dq_x = H_x = exp(q_x)K_x >= 0',
       'geometry_carrier':'finite dyadic coordinates in [-1/1000,1/1000], zero Dirichlet exterior',
       'one_source':serializable(one),'pair_rows':rows,
       'attractive_central_difference_force_at_six':list(map(str,force)),
       'force_over_native_E1_squared_over_4pi_kappa_r_squared':force_ratio,
       'Hessian_minimum_inverse_entry_diagnostic':float(inverse.min()),'mixed_minimum_energy_derivative_diagnostic':float(mixed),
       'stiffness_controls':scale_rows,'source_units_Newton_coefficient':'1/(4*pi*kappa)',
       'physical_Newton_constant_selected':False,'unchanged_reciprocal_coupling_equivalence':False,
       'geometry_minimum_preparation_is_declared':True,'geometry_register_has_autonomous_relaxation':False,
       'full_covariant_stress_and_Einstein_dynamics_derived':False,
       'verdict':'PASS_NATIVE_ENERGY_LAPSE_ATTRACTION'}
