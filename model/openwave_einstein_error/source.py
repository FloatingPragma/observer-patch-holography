from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import sympy as s
from model.openwave_source_ward import source as ward
from model.openwave_su2_spectral.spectrum import sqrt_upper

ETA=s.diag(-1,1,1,1)


def tracefree_matrix(x):
    return s.Matrix([[x[0],x[1],x[2],x[3]],
                     [x[1],x[4],x[5],x[6]],
                     [x[2],x[5],x[7],x[8]],
                     [x[3],x[6],x[8],x[0]-x[4]-x[7]]])


@lru_cache(None)
def design():
    t=s.sqrt(3)/3
    vectors=[s.Matrix(v) for v in ((1,1,0,0),(1,-1,0,0),(1,0,1,0),
        (1,0,-1,0),(1,0,0,1),(1,0,0,-1),(1,t,t,t),(1,t,t,-t),(1,t,-t,t))]
    assert all((v.T*ETA*v)[0]==0 for v in vectors)
    A=s.Matrix(9,9,lambda i,j:s.simplify((vectors[i].T*tracefree_matrix([int(k==j) for k in range(9)])*vectors[i])[0]))
    R=s.simplify(A.inv());assert R*A==s.eye(9)
    return vectors,A,R


def charges(B):
    return s.Matrix([s.simplify((v.T*B*v)[0]) for v in design()[0]])


def metric_coefficient(B):
    return s.trace(ETA*B)/4


def reconstruct(B):
    if B!=B.T:raise ValueError('symmetric tensor required')
    x=s.simplify(design()[2]*charges(B))
    R=tracefree_matrix(x)
    assert s.simplify(R-(B-metric_coefficient(B)*ETA))==s.zeros(4)
    return R


def interval_encode(lo,hi):
    return [str(lo),str(hi)]


def run():
    vectors,A,R=design()
    norms=[s.simplify(sum(abs(x) for x in R.row(i))) for i in range(9)]
    assert norms==[s.Rational(3,4),s.Rational(1,2),s.Rational(1,2),s.Rational(1,2),
                   s.Rational(5,4),2+s.sqrt(3),2+s.sqrt(3),s.Rational(5,4),(5+s.sqrt(3))/2]
    # The remaining spatial diagonal is the linear combination x0-x4-x7.
    last=R.row(0)-R.row(4)-R.row(7)
    assert sum(abs(v) for v in last)==s.Rational(5,4)
    C=2+s.sqrt(3)
    sign=s.Matrix([s.sign(x) for x in R.row(5)])
    assert s.simplify((R*sign)[5]-C)==0
    root_hi=sqrt_upper(F(3));root_lo=root_hi-F(1,2**70)
    assert root_lo**2<=3<=root_hi**2
    norm_bound=(F(2)+root_lo,F(2)+root_hi)
    epsilon=F(1,10**6)
    entry=(norm_bound[0]*epsilon,norm_bound[1]*epsilon)
    # Ward's actual finite-profile field; the declared metric assignment is
    # retained rather than calling these a measured physical curvature.
    profile,center,_=ward.profile();G=ward.einstein_static(ward.scale(-1,profile),F(1))
    sites=[center,(center[0]+1,center[1],center[2]),(0,0,0),(-1,0,0)]
    examples=[]
    for p in sites:
        B=s.Matrix(4,4,lambda i,j:s.Rational(G[i,j].get(p,F(0))))
        rem=reconstruct(B);q=charges(B)
        examples.append({'site':list(p),'null_charges_exact':[str(v) for v in q],
            'metric_coefficient_exact':str(metric_coefficient(B)),
            'tracefree_tensor_exact':[[str(x) for x in rem.row(i)] for i in range(4)]})
    # Null data alone leave the metric term entirely undetermined.
    assert charges(7*ETA)==s.zeros(9,1)
    return {'schema':'oph.quantitative_nine_null_tensor.v1','verdict':'PASS_SHARP_FINITE_TENSOR_ERROR_TRANSFER',
       'source_theorem':'Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean:tomographyDecoder_charge',
       'design_exact':[[str(x) for x in A.row(i)] for i in range(9)],
       'decoder_exact':[[str(x) for x in R.row(i)] for i in range(9)],
       'coordinate_error_constants_exact':[str(x) for x in norms],
       'sharp_maximum_entry_error_constant':'2 + sqrt(3)',
       'sharpness_error_sign_vector':[str(x) for x in sign],
       'operator_norm_error_upper_constant':'23/4 + 2*sqrt(3)',
       'example_uniform_charge_error_exact':str(epsilon),
       'maximum_tensor_entry_error_bound_interval':interval_encode(*entry),
       'native_profile_examples':examples,
       'ward_transport':{
           'forward_step_bound':'|lambda(y)-lambda(x)| <= a*w + 8*(2+sqrt(3))*epsilon',
           'central_two_step_bound':'|lambda(x+e)-lambda(x-e)| <= 2*a*w + 8*(2+sqrt(3))*epsilon',
           'path_bound':'sum the stated edge budgets along a path in the actual difference graph',
           'central_checkerboard_control':'lambda(x)=(-1)^x has zero central derivative and zero null data but is not constant on nearest-neighbour edges',
           'connectivity_required':'the actual derivative graph; central two-step components cannot be replaced by nearest-neighbour connectivity'},
       'scope':{'finite_error_premise_closed':True,'source_metric_assignment_proved':False,
          'physical_stress_identification_proved':False,'cosmological_constant_fixed_by_null_records':False,
          'continuum_or_nonlinear_einstein_theorem':False}}
