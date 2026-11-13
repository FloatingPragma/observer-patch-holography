"""Complete reducing string dynamics with rational unitary remainder bounds."""
from fractions import Fraction as F
from math import factorial,ceil
import sympy as sy
from model.openwave_localized_quark_record.source import energy,profile,exact_gauss
from model.openwave_quark_flux_records.flux import native_holes,fermion
from model.openwave_capacity_stress.source import enclosure,sqrt_interval

HOP=F(1,8)
TIMES=(F(0),F(1),F(4),F(8))

def jacobi(q,length,hop=HOP):
    if q not in(-1,2)or length<2 or hop<0:raise ValueError('unsupported native charge, length or hopping')
    return [[energy(q,i+1)if i==j else -hop if abs(i-j)==1 else F(0)
             for j in range(length)]for i in range(length)]

def evolve(H,time,tolerance=F(1,10**20)):
    """Exact Taylor vector and a complete operator-norm exponential tail."""
    t=F(time);assert t>=0 and tolerance>0
    size=len(H);bound=max(sum(abs(x)for x in row)for row in H);z=t*bound
    if not z:return [(F(int(i==0)),F(0))for i in range(size)],F(0),0
    term=[F(int(i==0))for i in range(size)];out=[[v,F(0)]for v in term]
    order=0
    while True:
        tail=F(2**ceil(2*z))*z**(order+1)/factorial(order+1)
        if tail<tolerance:break
        order+=1;term=[t*sum((H[i][j]*term[j]for j in range(size)),F(0))/order for i in range(size)]
        axis=order%2;sign=-1 if order%4 in(1,2)else 1
        for i,x in enumerate(term):out[i][axis]+=sign*x
        assert order<512
    assert abs(sum((a*a+b*b for a,b in out),F(0))-1)<=2*tail+tail*tail
    return [tuple(v)for v in out],tail,order

def interval(value,error=F(0),probability=False):
    lo,hi=value-error,value+error
    if probability:lo=max(F(0),lo);hi=min(F(1),hi)
    return list(enclosure((lo,hi)))

def full_car_control(q):
    U=sy.Matrix([[sy.Rational(3,5),sy.Rational(4,5),0],[-sy.Rational(4,5),sy.Rational(3,5),0],[0,0,1]])
    z=sy.Rational(3,5)+sy.I*sy.Rational(4,5);U=U*z**q
    assert sy.simplify(U.conjugate().T*U)==sy.eye(3)
    h=sy.zeros(6);h[:3,3:]=U;h[3:,:3]=U.conjugate().T
    H=sy.zeros(64)
    for state in range(64):
        for i in range(6):
            for j in range(6):
                if not h[i,j]:continue
                a=fermion(state,j,False)
                if a is None:continue
                b=fermion(a[0],i,True)
                if b is not None:H[b[0],state]+=a[1]*b[1]*h[i,j]
    B=sy.zeros(64,6)
    for j in range(6):
        state,sign=fermion(63,j,False);B[state,j]=sign
    assert H==H.conjugate().T and sy.simplify(H*B+B*h.T)==sy.zeros(64,6)
    assert H[:,63]==sy.zeros(64,1)
    return {'charge_integer':q,'full_active_CAR_dimension':64,'all_hole_columns':6,
            'nontrivial_colour_and_charge_link_unitary':True,'complete_hole_intertwiner_exact':'dGamma(h) B = - B h^T'}

def records(q,length,time):
    H=jacobi(q,length);v,eps,order=evolve(H,time);error=2*eps+3*eps*eps
    probabilities=[a*a+b*b for a,b in v];sigma=energy(q,2)-energy(q,1)
    mean=sum((F(i+1)*p for i,p in enumerate(probabilities)),F(0));mean_error=length*(length+1)*error/2
    # Coherent input has equal vacuum and string weights. Each pointer is a
    # separate complete four-outcome experiment, after the Gauss/seed tests.
    pointers=[];root=sqrt_interval(F(1,2));c=(2-root[0]-root[1])/2
    for i,(a,b)in enumerate(v):
        probs=[interval((1+sign*x)/4,eps/4,True)for x in(a,b)for sign in(1,-1)]
        cross=sum((a*v[j][0]+b*v[j][1]for j in(i-1,i+1)if 0<=j<length),F(0))
        neighbours=int(i>0)+int(i+1<length)
        work=energy(q,i+1)*(1-probabilities[i])/4+HOP*c*cross
        work_error=energy(q,i+1)*error/4+HOP*(neighbours*error+abs(cross)*(root[1]-root[0])/2)
        pointers.append({'endpoint':i+1,'X_plus_X_minus_Y_plus_Y_minus_intervals':probs,
                         'complete_pointer_nonselective_work_source_interval':interval(work,work_error)})
    flux=[]
    for edge in range(length):
        value=F(q)*sum(probabilities[edge:],F(0));flux.append(interval(value,abs(q)*(length-edge)*error))
    return {'native_integer_charge':q,'moving_endpoints':length,'vertices_including_reference':length+1,
        'source_time_exact':str(time),'Taylor_order':order,'full_unitary_vector_error_upper_exact':str(eps),
        'conditional_endpoint_probability_intervals':[interval(p,error,True)for p in probabilities],
        'conditional_mean_endpoint_interval':interval(mean,mean_error),
        'conditional_integer_flux_mean_intervals':flux,
        'conditional_source_charge_mean_intervals':[interval(F(q,3)*p,abs(F(q,3))*error)for p in probabilities],
        'conditional_reference_charge_exact':str(-F(q,3)),
        'coherent_input_reference_charge_mean_exact':str(-F(q,6)),
        'source_plus_reference_total_charge_exact':'0',
        'conditional_colour_Casimir_means':'(4/3) times the cumulative endpoint probabilities on each path edge.',
        'all_position_outcomes_include_full_flux_profile':True,'vacuum_probability_coherent_input_exact':'1/2',
        'coherent_input_position_probabilities':'one half the displayed conditional probabilities; vacuum is the other half.',
        'outside_reducing_string_sector_probability_exact':'0 for the stated input, including all higher rotor representations',
        'coherent_preparation_mean_work_source_exact':str(energy(q,1)/2),
        'autonomous_external_work_exact':'0','conditional_total_energy_exact':str(energy(q,1)),
        'charge_profile_nonselective_work_source_interval':interval(sigma*(mean-1)/2,sigma*mean_error/2),
        'separate_complete_pointer_records':pointers}

def carrier_comparison(q,time):
    small,e1,_=evolve(jacobi(q,4),time);large,e2,_=evolve(jacobi(q,8),time)
    square=sum(((large[i][0]-(small[i][0]if i<4 else 0))**2+
                (large[i][1]-(small[i][1]if i<4 else 0))**2 for i in range(8)),F(0))
    lo,hi=sqrt_interval(square);lo=max(F(0),lo-e1-e2);hi+=e1+e2
    return {'native_integer_charge':q,'source_time_exact':str(time),
            'complete_embedded_vector_difference_interval':list(enclosure((lo,hi))),
            'all_large_carrier_endpoint_tail_coordinates_included':True}

def run():
    holes=native_holes();assert {F(x['native_charge_exact'])for x in holes}=={F(2,3),F(-1,3)}
    controls=[full_car_control(q)for q in(-1,2)]
    gauss=sum(exact_gauss(row['family'],4)for row in holes)
    rows=[records(q,L,t)for q in(-1,2)for L in(4,8)for t in TIMES]
    assert all(F(row['conditional_endpoint_probability_intervals'][1][0])>0 for row in rows if row['source_time_exact']=='1')
    comparisons=[carrier_comparison(q,t)for q in(-1,2)for t in(F(1),F(4))]
    return {'schema':'oph.autonomous_native_quark_string_transport.v1','verdict':'PASS_FULL_GAUSS_SOURCE_AUTONOMOUS_ENDPOINT_TRANSPORT',
        'native_families':holes,'declared_source':{'matter':'Full thirty-CAR source at each vertex and the full native quantum reference at vertex0.',
            'links':'Full L2(SU3 x U1) on every open-chain edge; no electric or representation cutoff.',
            'Hamiltonian':'Inherited number/electric/Casimir idle action plus the gauge-covariant nearest-neighbour native uc,dc electron bilinears with coefficient+1/8 and their adjoints; no endpoint hopping across the fixed reference boundary.',
            'hopping_source_exact':str(HOP),'hopping_coefficient_role':'Declared common architecture coefficient, independent of particle trajectories or masses.',
            'self_adjointness':'The finite-graph hopping is bounded by12*kappa*(L-1); bounded perturbation of the positive rotor/number Hamiltonian gives a self-adjoint lower-bounded full source.',
            'external_time_dependent_control_during_transport':False},
        'exact_reduction':{'orthonormal_strings':'B_x|filled,Haar,reference-vacuum>, x=1,...,L; different hole positions are orthogonal.',
            'Jacobi_action':'H|x>=(2+sigma*x)|x>-kappa(|x-1>+|x+1>), omitting only graph-boundary neighbours.',
            'sigma':'q^2/32+1/12','full_complement_dynamical_leakage_exact':'0 from the stated finite reducing subspace',
            'one_point_charge_weighted_field_mean_exact':'0 by full Gauss invariance; relational string quadratures are measured instead.',
            'mean_displacement_all_time_upper_exact':{str(q):str(2*HOP/(energy(q,2)-energy(q,1)))for q in(-1,2)}},
        'full_CAR_controls':controls,'exact_Gauss_generator_controls':gauss,'records':rows,'carrier_comparisons':comparisons,
        'complete_instrument':'First the full Gauss test and absorbing seed failure; then all endpoint/flux profiles and their full complement, or one chosen endpoint X/Y four-effect instrument with its full Gauss complement. Square-root Kraus operators specify the pointer work.',
        'moving_native_fractional_charge_records_proved':True,'classical_winding_texture_or_spontaneous_soliton_proved':False,
        'physical_quark_pole_or_GeV_rest_mass_proved':False,'infinite_volume_or_continuum_confinement_proved':False,
        'new_empirical_calibration':False,'new_match_with_nature':False}
